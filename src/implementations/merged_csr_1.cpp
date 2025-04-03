#include "graph.hpp"
#include <cstdint>
#include <cstdio>
#include <fstream>
#include <limits>
#include <omp.h>
#include <sys/wait.h>

#define DEGREE(vertex) merged_csr[vertex]
#define DISTANCE(vertex) merged_csr[vertex + 1]
#define THREAD(vertex) merged_csr[vertex + 2]
#define METADATA 3

MergedCSR_1::MergedCSR_1(Graph *graph) : BFS_Impl(graph) {
  create_merged_csr();
}

MergedCSR_1::~MergedCSR_1() { delete[] merged_csr; }

// Create merged CSR from CSR
void MergedCSR_1::create_merged_csr() {
  merged_csr = new eidType[graph->M + METADATA * graph->N];
  merged_rowptr = new eidType[graph->N];
  eidType merged_index = 0;

  for (vidType i = 0; i < graph->N; i++) {
    eidType start = graph->rowptr[i];
    // Add degree to start of neighbor list
    merged_csr[merged_index++] = graph->rowptr[i + 1] - graph->rowptr[i];
    // Initialize distance
    merged_csr[merged_index++] = std::numeric_limits<weight_type>::max();
    // Initialize thread id
    merged_csr[merged_index++] = 0;
    // Copy neighbors
    for (eidType j = start; j < graph->rowptr[i + 1]; j++) {
      merged_csr[merged_index++] =
          graph->rowptr[graph->col[j]] + METADATA * graph->col[j];
    }
  }
  // Fix rowptr indices caused by adding the degree to the start of each
  // neighbor list
  for (vidType i = 0; i <= graph->N; i++) {
    merged_rowptr[i] = graph->rowptr[i] + METADATA * i;
  }
}

// Extract distances from merged CSR
void MergedCSR_1::compute_distances(weight_type *distances,
                                    vidType source) const {
#pragma omp parallel for simd schedule(static)
  for (vidType i = 0; i < graph->N; i++) {
    distances[i] = DISTANCE(merged_rowptr[i]);
    // Reset distance for next BFS
    DISTANCE(merged_rowptr[i]) = std::numeric_limits<weight_type>::max();
    THREAD(merged_rowptr[i]) = 0;
  }
  distances[source] = 0;
}

#pragma omp declare reduction(vec_add                                          \
:frontier : omp_out.insert(omp_out.end(), omp_in.begin(), omp_in.end()))

void MergedCSR_1::top_down_step(const frontier &this_frontier,
                                frontier &next_frontier,
                                const weight_type &distance) {
  for (const auto &v : this_frontier) {
    eidType end = v + METADATA + DEGREE(v);
// Iterate over neighbors
#pragma omp simd
    for (eidType i = v + METADATA; i < end; i++) {
      eidType neighbor = merged_csr[i];
      // If neighbor is not visited, add to frontier
      if (DISTANCE(neighbor) == std::numeric_limits<weight_type>::max()) {
        if (DEGREE(neighbor) != 1) {
          next_frontier.push_back(neighbor);
        }
        DISTANCE(neighbor) = distance;
      }
    }
  }
}

void MergedCSR_1::check_vertex(eidType v, frontier &next_frontier,
                               const weight_type &local_distance) {
  if (v != std::numeric_limits<vidType>::max()) {
    eidType end = v + METADATA + DEGREE(v);
    for (eidType j = v + METADATA; j < end; j++) {
      eidType neighbor = merged_csr[j];
      if (DISTANCE(neighbor) > local_distance) {
        next_frontier.push_back(neighbor);
      }
    }
  }
}

std::vector<std::vector<uint32_t>> frontier_sizes;

void MergedCSR_1::balance_threads(frontier &next_frontier,
                                  weight_type &local_distance) {
  if (next_frontier.size() > 50) {
    // Guarantee that only one thread will enter this block and that the reader
    // has finished reading
    if (waiting_threads > 0) {
      if (comm_lock->try_acquire()) {
        if (waiting_threads > 0) {
          waiting_threads--;
          size_t half_size = next_frontier.size() / 2;
          comm_frontier.clear();
          comm_frontier.insert(comm_frontier.end(),
                               next_frontier.begin() + half_size,
                               next_frontier.end());
          next_frontier.erase(next_frontier.begin() + half_size,
                              next_frontier.end());
          comm_distance = local_distance;
          readers_lock->release();
        } else {
          comm_lock->release();
        }
      }
    }
  } else if (next_frontier.size() == 0 || thread_errors[omp_get_thread_num()] > 100) {
    waiting_threads++;
    if (waiting_threads == omp_get_num_threads()) {
      readers_lock->release(omp_get_num_threads());
    }
    readers_lock->acquire();
    if (waiting_threads < omp_get_num_threads()) {
      next_frontier = std::move(comm_frontier);
      local_distance = comm_distance++;
      thread_errors[omp_get_thread_num()] = 0;
      comm_lock->release();
    }
  }
}

double waiting_time = 0;

void MergedCSR_1::critical_writeback(frontier &next_frontier,
                                     const weight_type &local_distance) {
  double start = omp_get_wtime();
#pragma omp critical
  {
    waiting_time += omp_get_wtime() - start;
    uint32_t thread_id = omp_get_thread_num();
    frontier_sizes[thread_id].push_back(next_frontier.size());
    for (auto &v : next_frontier) {
      if (DISTANCE(v) > local_distance) {
        if (DISTANCE(v) != std::numeric_limits<weight_type>::max()) {
          overwrites++;
          thread_errors[THREAD(v)]++;
        }
        THREAD(v) = thread_id;
        DISTANCE(v) = local_distance;
        if (DEGREE(v) == 1) {
          v = std::numeric_limits<eidType>::max();
        }
      } else {
        reduntants++;
        v = std::numeric_limits<eidType>::max();
      }
    }
  }
}

void MergedCSR_1::parallel_frontiers(const frontier &this_frontier,
                                     weight_type &distance) {
#pragma omp parallel
  {
    frontier this_frontier_loc, next_frontier_loc;
    weight_type local_distance = distance;
    vidType chunk_size = this_frontier.size() / omp_get_num_threads();
    // Iterate over vertices assigned to this thread
    int start = omp_get_thread_num() * chunk_size;
    int end = (omp_get_thread_num() == omp_get_num_threads() - 1)
                  ? this_frontier.size()
                  : start + chunk_size;
    for (vidType i = start; i < end; i++) {
      check_vertex(this_frontier[i], next_frontier_loc, distance);
    }
    critical_writeback(next_frontier_loc, local_distance);
    local_distance++;
    this_frontier_loc = std::move(next_frontier_loc);
    // Expand frontier
    while (!this_frontier_loc.empty()) {
      next_frontier_loc.clear();
      for (const auto &v : this_frontier_loc) {
        check_vertex(v, next_frontier_loc, local_distance);
      }
      critical_writeback(next_frontier_loc, local_distance);
      balance_threads(next_frontier_loc, local_distance);
      local_distance++;
      this_frontier_loc = std::move(next_frontier_loc);
    }
  }
}

void MergedCSR_1::BFS(vidType source, weight_type *distances) {
  frontier_sizes.resize(omp_get_max_threads());
  // Initialize synchronization variables
  waiting_threads = 0;
  readers_lock = new semaphore(0);
  comm_lock = new semaphore(1);

  overwrites = 0;
  reduntants = 0;
  waiting_time = 0;
  thread_errors.resize(omp_get_max_threads(), 0);

  frontier this_frontier;
  eidType start = merged_rowptr[source];
  this_frontier.push_back(start);
  DISTANCE(start) = 0;
  weight_type distance = 1;
  while (this_frontier.size() < 50 && !this_frontier.empty()) {
    frontier next_frontier;
    next_frontier.reserve(this_frontier.size());
    top_down_step(this_frontier, next_frontier, distance);
    distance++;
    this_frontier = std::move(next_frontier);
  }
  if (!this_frontier.empty()) {
    parallel_frontiers(this_frontier, distance);
  }
  compute_distances(distances, source);
  printf("Overwrites: %llu\n", overwrites);
  printf("Redundants: %llu\n", reduntants);
  printf("Waiting time: %f\n", waiting_time/omp_get_max_threads());
  // Store frontier sizes into a file
  // std::ofstream outfile("scripts/frontier_sizes.txt");
  // if (outfile.is_open()) {
  //   for (size_t i = 0; i < frontier_sizes.size(); i++) {
  //     outfile << "Thread " << i << ": ";
  //     for (const auto &size : frontier_sizes[i]) {
  //       outfile << size << " ";
  //     }
  //     outfile << "\n";
  //   }
  //   outfile.close();
  // }
  // frontier_sizes.clear();
}

bool MergedCSR_1::check_result(vidType source, weight_type *distances) {
  return BFS_Impl::check_distances(source, distances);
}