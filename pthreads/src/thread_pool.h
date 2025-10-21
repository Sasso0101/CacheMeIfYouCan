#ifndef THREAD_POOL_H
#define THREAD_POOL_H

/**
 * @file thread_pool.h
 * @brief Implementation of a simple thread pool.
 *
 * This file contains the declaration of functions to initialize, manage,
 * and destroy a thread pool. The pool uses condition variables for synchronization
 * between the main thread and worker threads. It supports dispatching work cycles
 * and waiting for their completion, as well as graceful termination.
 * Includes Linux-specific thread pinning for potential performance optimization.
 */

#include "config.h"
#include <pthread.h>
#include <stdatomic.h>
#include <stdbool.h>
#include <sys/types.h>

typedef struct {
  pthread_cond_t cond_children;
  pthread_mutex_t mutex_children;
  pthread_t threads[MAX_THREADS];
  int thread_ids[MAX_THREADS];
  atomic_uint run_id; // Counter for work cycles
  atomic_bool stop_threads; // Flag to signal threads to terminate

  atomic_bool children_done; // Flag indicating if workers completed the current cycle
  pthread_cond_t cond_parent;
  pthread_mutex_t mutex_parent;

  void *(*routine)(void *); // Store the worker function pointer
} thread_pool_t;

extern thread_pool_t tp;

/**
 * @brief Initializes a thread pool structure.
 *
 * Sets up the necessary mutexes, condition variables, and initial state
 * flags and counters for the thread pool operation.
 *
 * @param tp Pointer to the thread_pool_t structure to initialize.
 * @param routine The function pointer for the task that worker threads will execute.
 */
void init_thread_pool(thread_pool_t *tp, void *(*routine)(void *));

/**
 * @brief Waits for a signal indicating new work or termination.
 *
 * Called by worker threads. Blocks the calling thread using a condition variable
 * until the main thread signals a new work cycle (by incrementing `tp->run_id`)
 * or signals termination (`tp->stop_threads`).
 *
 * @param tp Pointer to the thread pool structure.
 * @param run_id Pointer to the worker thread's current run ID. This value is
 *               compared against the pool's global run ID and incremented by
 *               this function when the wait condition is met.
 * @return Always returns 0 in this implementation if the thread should proceed
 *         (to either work or exit cleanly based on `stop_threads`).
 *         Calls pthread_exit directly if termination is signaled.
 */
int wait_for_work(thread_pool_t *tp, uint32_t *run_id);

/**
 * @brief Creates and launches the worker threads in the pool.
 *
 * Spawns `MAX_THREADS` worker threads. Each thread is configured to execute
 * the `thread_main_wrapper` function. The thread's index (0 to MAX_THREADS-1)
 * is passed as the argument to `thread_main_wrapper`, which is then forwarded
 * to the user's routine.
 * On Linux systems, each thread is pinned to a specific logical CPU core
 * sequentially (thread 0 to core 0, thread 1 to core 1, etc.) using
 * `pin_thread_to_cpu` for potential performance benefits.
 * Exits the program with an error message if thread creation fails.
 *
 * @param tp Pointer to the initialized thread_pool_t structure. Thread handles
 *           will be stored in `tp->threads`.
 */
void thread_pool_create(thread_pool_t *tp);

/**
 * @brief Signals worker threads to start a new work cycle and waits for completion notification.
 *
 * @param tp Pointer to the thread_pool_t structure.
 */
void thread_pool_start_wait(thread_pool_t *tp);

/**
 * @brief Signals all worker threads to stop and waits for their termination.
 *
 * @param tp Pointer to the thread_pool_t structure.
 */
void thread_pool_terminate(thread_pool_t *tp);

/**
 * @brief Notifies the main (parent) thread that the current work cycle is complete.
 *
 * This function is intended to be called by a worker thread (or potentially
 * coordinated among workers) after they have finished their tasks for the
 * current `run_id` cycle. It signals the parent thread waiting in
 * `thread_pool_start_wait`.
 *
 * @note The exact logic for *when* this is called depends on the application using
 *       the pool. It might be called by the last worker to finish, or after a
 *       specific synchronization point among workers. This implementation provides
 *       the basic notification mechanism.
 *
 * @param tp Pointer to the thread_pool_t structure.
 */
 void thread_pool_notify_parent(thread_pool_t *tp);

/**
 * @brief Destroys the synchronization primitives used by the thread pool.
 *
 * Releases the resources associated with the mutexes and condition variables
 * (both children's and parent's). This function should be called only *after*
 * all worker threads have been joined (e.g., after `thread_pool_terminate`
 * has completed successfully). Destroying these primitives while threads might
 * still be using them results in undefined behavior.
 *
 * @param tp Pointer to the thread_pool_t structure whose resources are to be freed.
 */
void destroy_thread_pool(thread_pool_t *tp);

#endif