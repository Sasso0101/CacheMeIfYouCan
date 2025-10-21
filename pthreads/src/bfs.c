#define _GNU_SOURCE
#include "cli_parser.h"
#include "config.h"
#include "debug_utils.h"
#include "frontier.h"
#include "merged_csr.h"
#include "mmio_c_wrapper.h"
#include "mt19937-64.h"
#include "thread_pool.h"
#include <assert.h>
#include <pthread.h>
#include <stdatomic.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#ifdef USE_PAPI
#include "papi.h"
#endif

MergedCSR *merged_csr;
Frontier *f1, *f2;
uint32_t *distances;

atomic_int active_threads;

volatile uint32_t exploration_done;
volatile int distance;

int max_chunks;

thread_pool_t tp;

void top_down_chunk(MergedCSR *merged_csr, Frontier *next, Chunk *c,
                    Chunk **dest, int distance, int thread_id) {
  assert(c != NULL && "Chunk passed to top_down_chunk is NULL!");
  mer_t v = VERT_MAX;
  while ((v = chunk_pop_vertex(c)) != VERT_MAX) {
    mer_t end = v + DEGREE(merged_csr, v) + METADATA_SIZE;
    for (mer_t i = v + METADATA_SIZE; i < end; i++) {
      mer_t neighbor = merged_csr->merged[i];
      if (DISTANCE(merged_csr, neighbor) == UINT32_MAX) {
        DISTANCE(merged_csr, neighbor) = distance;
        if (DEGREE(merged_csr, neighbor) != 1) {
          if (*dest == NULL || (*dest)->next_free_index >= CHUNK_SIZE) {
            *dest = frontier_create_chunk(next, thread_id);
          }
          chunk_push_vertex(*dest, neighbor);
        }
      }
    }
  }
}

void top_down(MergedCSR *merged_csr, Frontier *current_frontier,
              Frontier *next_frontier, int distance, int thread_id) {
  Chunk *c = NULL;
  Chunk *next_chunk = NULL;
  Chunk **dest = &next_chunk;
  // Run top-down step for all chunks belonging to the thread
  while ((c = frontier_remove_chunk(current_frontier, thread_id)) != NULL) {
    top_down_chunk(merged_csr, next_frontier, c, dest, distance, thread_id);
  }
  // Work stealing from other threads when finished processing chunks of this
  // thread
  bool work_to_do = true;
  while (work_to_do) {
    work_to_do = false;
    for (int i = 0; i < MAX_THREADS; i++) {
      // Skip threads that have only one chunk left
      // This is a heuristic to avoid that threads that start with no chunks
      // steal from threads that have one chunk This situation is common when
      // there are more threads than chunks This is not a perfect solution, but
      // it works okay for now
      if (current_frontier->thread_chunks[i]->top_chunk > 1) {
        work_to_do = 1;
        if ((c = frontier_remove_chunk(current_frontier, i)) != NULL) {
          top_down_chunk(merged_csr, next_frontier, c, dest, distance,
                         thread_id);
        }
        i--;
      }
    }
  }
}

void finalize_distances(MergedCSR *merged_csr, int thread_id) {
  // Write distances from mergedCSR to distances array
  mer_t chunk_size = merged_csr->num_vertices / MAX_THREADS;
  mer_t start = thread_id * chunk_size;
  mer_t end = (thread_id == MAX_THREADS - 1) ? merged_csr->num_vertices
                                             : (thread_id + 1) * chunk_size;
  for (mer_t i = start; i < end; i++) {
    distances[i] = DISTANCE(merged_csr, merged_csr->row_ptr[i]);
    DISTANCE(merged_csr, merged_csr->row_ptr[i]) = UINT32_MAX;
  }
}

void *thread_main(void *arg) {
  int thread_id = *(int *)arg;

  while (!exploration_done) {
    int old = distance;
    top_down(merged_csr, f1, f2, distance, thread_id);
    if (atomic_fetch_sub(&active_threads, 1) == 1) {
      // Swap frontiers
      active_threads = MAX_THREADS;
      Frontier *temp = f2;
      f2 = f1;
      f1 = temp;
      int chunks = frontier_get_total_chunks(f1);
      if (chunks == 0)
        exploration_done = 1;
      if (chunks > max_chunks)
        max_chunks = chunks;
      // printf("%u \n", distance);
      // print_chunk_counts(f1);
      atomic_thread_fence(memory_order_seq_cst);
      distance++;
    }
    while (distance == old)
      ;
  }
  finalize_distances(merged_csr, thread_id);

  if (atomic_fetch_sub(&active_threads, 1) == 1) {
    // printf("Max distance: %u\n", distance);
    thread_pool_notify_parent(&tp);
  }
  return NULL;
}

void initialize_bfs(const mmio_csr_u32_f32_t *graph) {
  merged_csr = to_merged_csr(graph);
  f1 = frontier_create();
  f2 = frontier_create();
  init_thread_pool(&tp, thread_main);
  thread_pool_create(&tp);
}

void bfs(uint32_t source) {
  // Convert source vertex to mergedCSR index
  source = merged_csr->row_ptr[source];
  DISTANCE(merged_csr, source) = 0;
  Chunk *c = frontier_create_chunk(f1, 0);
  chunk_push_vertex(c, source);
  exploration_done = 0;
  active_threads = MAX_THREADS;
  distance = 1;
  max_chunks = 0;
  atomic_thread_fence(memory_order_seq_cst);
  thread_pool_start_wait(&tp);
}

uint32_t *generate_sources(const mmio_csr_u32_f32_t *graph, int runs,
                           uint32_t num_vertices, uint32_t source) {
  uint32_t *sources = (uint32_t *)malloc(runs * sizeof(uint32_t));
  if (source != UINT32_MAX) {
    for (int i = 0; i < runs; i++) {
      sources[i] = source;
    }
  } else {
    init_genrand64(SEED);
    for (int i = 0; i < runs; i++) {
      do {
        uint64_t gen = genrand64_int64();
        sources[i] = (uint32_t)gen % num_vertices;
      } while (graph->row_ptr[sources[i] + 1] - graph->row_ptr[sources[i]] ==
               0);
    }
  }
  return sources;
}

typedef struct {
  char *filename; // Will be allocated by the parser
  int runs;
  int source_id;
  bool check;
  bool output;
} AppArgs;

int main(int argc, char **argv) {
  AppArgs args = {.filename = NULL,
                  .runs = 1,
                  .source_id = -1,
                  .check = false,
                  .output = true};
  const CliOption options[] = {
      {'f', "file", "Load graph from file", ARG_TYPE_STRING, &args.filename,
       true},
      {'n', "runs", "Number of runs", ARG_TYPE_INT, &args.runs, false},
      {'s', "source", "ID of source vertex", ARG_TYPE_INT, &args.source_id,
       false},
      {'c', "check", "Checks BFS correctness", ARG_TYPE_BOOL, &args.check,
       false}};
  int num_options = sizeof(options) / sizeof(options[0]);
  const char *app_description =
      "An optimized BFS algorithm for large-diameter graphs.";
  int parse_result =
      cli_parse(argc, argv, options, num_options, app_description);

  // Check the result: 0 is success, 1 means help was printed, -1 is an error.
  if (parse_result != 0) {
    if (args.filename)
      free(args.filename);
    return (parse_result == 1) ? 0 : 1;
  }

  mmio_csr_u32_f32_t *graph = mmio_read_csr_u32_f32(args.filename, false);
  if (graph == NULL) {
    printf("Failed to import graph from file [%s]\n", args.filename);
    return -1;
  }

  uint32_t *sources =
    generate_sources(graph, args.runs, graph->nrows, args.source_id);

  distances = (uint32_t *)malloc(graph->nrows * sizeof(uint32_t));
  memset(distances, UINT32_MAX, graph->nrows * sizeof(uint32_t));
  initialize_bfs(graph);

  struct timespec start, end;
  double elapsed;
  for (int i = 0; i < args.runs; i++) {
    #ifndef USE_PAPI
    clock_gettime(CLOCK_MONOTONIC, &start);
    #endif
    bfs(sources[i]);
    #ifndef USE_PAPI
    clock_gettime(CLOCK_MONOTONIC, &end);
    long seconds = end.tv_sec - start.tv_sec;
    long nanoseconds = end.tv_nsec - start.tv_nsec;
    elapsed = seconds + nanoseconds * 1e-9;

    printf(
        "run_id=%d,diameter=%d,threads=%d,chunk_size=%d,max_chunks=%d,source=%d,%.4f\n",
        i, distance, MAX_THREADS, CHUNK_SIZE, max_chunks, sources[i], elapsed);

    if (args.check) {
      check_bfs_correctness(graph, distances, sources[i]);
    }

    memset(distances, UINT32_MAX, merged_csr->num_vertices * sizeof(uint32_t));
    #endif
  }
  // Terminate threads
  thread_pool_terminate(&tp);

  free(sources);
  free(graph->row_ptr);
  free(graph->col_idx);
  free(graph);
  frontier_destroy(f1);
  frontier_destroy(f2);
  destroy_thread_pool(&tp);
  destroy_merged_csr(merged_csr);
  free(distances);

  return 0;
}
