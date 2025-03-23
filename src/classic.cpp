#include "graph.hpp"

Classic::Classic(eidType *rowptr, vidType *col, uint64_t N, uint64_t M)
    : BaseGraph(rowptr, col, N, M) {
  visited = new bool[N];
}

inline void Classic::set_distance(vidType i, weight_type distance,
                                  weight_type *distances) {
  distances[i] = distance;
  visited[i] = true;
}

inline void Classic::add_to_frontier(frontier &frontier, vidType v,
                                     vidType &edges_frontier) {
  frontier.push_back(v);
  edges_frontier += rowptr[v + 1] - rowptr[v];
}

#pragma omp declare reduction(                                                 \
        vec_add : std::vector<vidType>,                                        \
            std::vector<std::pair<vidType, bool>> : omp_out.insert(            \
                    omp_out.end(), omp_in.begin(), omp_in.end()))

void Classic::bottom_up_step(frontier this_frontier, frontier &next_frontier,
                             weight_type distance, weight_type *distances,
                             vidType &edges_frontier) {
#pragma omp parallel for reduction(vec_add : next_frontier)                    \
    reduction(+ : edges_frontier) schedule(static)
  for (vidType i = 0; i < N; i++) {
    if (!visited[i]) {
      for (vidType j = rowptr[i]; j < rowptr[i + 1]; j++) {
        if (visited[col[j]] && distances[col[j]] == distance - 1) {
          // If neighbor is in frontier, add this vertex to next frontier
          if (rowptr[i + 1] - rowptr[i] > 1) {
            add_to_frontier(next_frontier, i, edges_frontier);
          }
          set_distance(i, distance, distances);
          break;
        }
      }
    }
  }
}

void Classic::top_down_step(frontier this_frontier, frontier &next_frontier,
                            weight_type &distance, weight_type *distances,
                            vidType &edges_frontier,
                            vidType edges_frontier_old) {
#pragma omp parallel for reduction(vec_add : next_frontier)                    \
    reduction(+ : edges_frontier)                                              \
    schedule(static) if (edges_frontier_old > 150)
  for (const auto &v : this_frontier) {
    for (vidType i = rowptr[v]; i < rowptr[v + 1]; i++) {
      vidType neighbor = col[i];
      if (!visited[neighbor]) {
        if (rowptr[neighbor + 1] - rowptr[neighbor] > 1) {
          add_to_frontier(next_frontier, neighbor, edges_frontier);
        }
        set_distance(neighbor, distance, distances);
      }
    }
  }
}

void Classic::BFS(vidType source, weight_type *distances) {
  eidType unexplored_edges = M;
  vidType edges_frontier_old = 0;
  frontier this_frontier;
  Direction dir = Direction::TOP_DOWN;
  vidType edges_frontier = 0;
  add_to_frontier(this_frontier, source, edges_frontier);
  set_distance(source, 0, distances);
  weight_type distance = 1;
  while (!this_frontier.empty()) {
    frontier next_frontier;
    next_frontier.reserve(this_frontier.size());
    if (dir == Direction::BOTTOM_UP && this_frontier.size() < N / BETA) {
      dir = Direction::TOP_DOWN;
    } else if (dir == Direction::TOP_DOWN &&
               edges_frontier > unexplored_edges / ALPHA) {
      dir = Direction::BOTTOM_UP;
    }
    unexplored_edges -= edges_frontier;
    edges_frontier_old = edges_frontier;
    edges_frontier = 0;
    if (dir == Direction::TOP_DOWN) {
      top_down_step(this_frontier, next_frontier, distance, distances,
                    edges_frontier, edges_frontier_old);
    } else {
      bottom_up_step(this_frontier, next_frontier, distance, distances,
                     edges_frontier);
    }
    distance++;
    this_frontier = std::move(next_frontier);
  }
}