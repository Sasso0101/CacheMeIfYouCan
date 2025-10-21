#ifndef FRONTIER_H
#define FRONTIER_H

/**
 * @brief Thread-local vertex chunk pool for parallel graph traversal.
 *
 * The Frontier data structure manages per-thread pools of reusable vertex
 * chunks for use in parallel graph algorithms. Each thread maintains its own
 * chunk pool to minimize contention and improve cache locality.
 *
 * ## Design Overview
 * - The Frontier is composed of MAX_THREADS thread-local chunk pools.
 * - Each pool (`ChunkPool`) contains a dynamically resizable array of
 *   pointers to `VertexChunk` blocks.
 * - A `VertexChunk` holds a fixed-size array of vertices (CHUNK_SIZE) and acts
 *   as a LIFO (stack-style) buffer.
 * - Threads acquire and release chunks from their own pool in a lock-protected
 * manner.
 * - Chunks are not freed when released — they are reused, minimizing dynamic
 * allocations.
 *
 * ## Concurrency Model
 * - Chunk acquisition and release are thread-safe for individual threads (each
 * thread locks only its own pool).
 * - No synchronization is required for accessing the internal content of a
 * `VertexChunk` if it is used exclusively by one thread.
 * - Cross-thread stealing is supported via inspection functions.
 */

#include "config.h"
#include <pthread.h>
#include <stdbool.h>

typedef mer_t ver_t;

typedef struct {
  ver_t vertices[CHUNK_SIZE];
  int next_free_index;
} Chunk;

typedef struct {
  Chunk **chunks;        // Array of pointers to chunks
  int chunks_size;       // Current size of the chunks array
  int top_chunk;         // Index of the next chunk to be allocated
  pthread_mutex_t lock;  // Mutex for thread-safe access
} ThreadChunks;

typedef struct {
  ThreadChunks **thread_chunks;
  int *thread_chunk_counts;
} Frontier;

/**
 * Creates and initializes a new Frontier structure. Allocates memory for
 * per-thread vertex chunk pools and sets up mutexes for thread-safe chunk
 * acquisition and release.
 */
Frontier *frontier_create();

/**
 * Destroys a Frontier structure and deallocates all associated resources.
 */
void frontier_destroy(Frontier *f);

/**
 * Creates a new chunk for the specified thread in the Frontier. Allocates a new
 * Chunk and returns a pointer to it.
 */
Chunk *frontier_create_chunk(Frontier *f, int thread_id);

/**
 * Removes a chunk from the Frontier for the specified thread. Returns a pointer
 * to the removed Chunk, or NULL if no chunks are available.
 */
Chunk *frontier_remove_chunk(Frontier *f, int thread_id);

/**
 * Pushes a vertex onto the specified chunk. The chunk must not be full.
 * Asserts if the chunk is NULL or full.
 */
void chunk_push_vertex(Chunk *c, mer_t v);

/**
 * Pops a vertex from the specified chunk. Returns VERT_MAX if the chunk is
 * empty.
 */
mer_t chunk_pop_vertex(Chunk *c);

/**
 * Gets the total number of chunks currently allocated in the Frontier.
 */
int frontier_get_total_chunks(Frontier *f);

#endif // FRONTIER_H
