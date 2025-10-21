// Define _GNU_SOURCE to enable non-standard GNU extensions,
// specifically needed for pthread_setaffinity_np used in thread pinning.
#define _GNU_SOURCE
#include "thread_pool.h"
#include <pthread.h>
#include <sched.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#ifdef USE_PAPI
#include <papi.h>
#endif

void init_thread_pool(thread_pool_t *tp, void *(*routine)(void *)) {
  // Initialize synchronization primitives for worker threads
  pthread_mutex_init(&tp->mutex_children, NULL);
  pthread_cond_init(&tp->cond_children, NULL);

  // Initialize synchronization primitives for the main (parent) thread
  pthread_mutex_init(&tp->mutex_parent, NULL);
  pthread_cond_init(&tp->cond_parent, NULL);

  tp->run_id = 0;
  tp->stop_threads = false;
  tp->children_done = false;
  tp->routine = routine;
}

int wait_for_work(thread_pool_t *tp, uint *run_id) {
  // Lock the mutex protecting shared state accessed by worker threads
  pthread_mutex_lock(&tp->mutex_children);

  // Wait while the worker's run_id is greater than the global run_id
  // This ensures the worker doesn't proceed until the main thread signals a new cycle
  while (tp->run_id < *run_id) {
    pthread_cond_wait(&tp->cond_children, &tp->mutex_children);
  }

  *run_id += 1;
  pthread_mutex_unlock(&tp->mutex_children);
  if (tp->stop_threads) {
    #ifdef USE_PAPI
    int retval = PAPI_hl_region_end("computation");
    if (retval != PAPI_OK) {
      printf("PAPI error %d: %s\n", retval, PAPI_strerror(retval));
      exit(1);
    }
    #endif
    pthread_exit(0);
  } else {
    #ifdef USE_PAPI
    if (*run_id == 3) { // skip first two iterations
      int retval = PAPI_hl_region_begin("computation");
      if ( retval != PAPI_OK ) {
        printf("PAPI error %d: %s\n", retval, PAPI_strerror(retval));
        exit(1);
      }
    }
    #endif
    return 0;
  }
}

/**
 * @brief The main function executed by each worker thread.
 *
 * This function serves as the entry point for threads created by `pthread_create`.
 * It runs a loop that repeatedly calls `wait_for_work` to wait for tasks.
 * If `wait_for_work` returns 0 (indicating work is available and termination is not
 * requested *yet*), it executes the user-provided `tp->routine` function, passing
 * the `arg` parameter received during thread creation.
 * The loop continues as long as `wait_for_work` returns 0. If `wait_for_work`
 * detects the `stop_threads` flag, it calls `pthread_exit`, terminating this loop.
 *
 * @param arg The argument passed to the thread during creation (e.g., thread ID).
 *            This is forwarded to the `tp->routine` function.
 * @return NULL (The function typically runs indefinitely until terminated via `pthread_exit`).
 * @internal
 */
void *thread_main_wrapper(void *arg) {
  // Initialize the worker's local run ID. Starts at 1 to wait for the first cycle (global run id is initalized at 0).
  uint run_id = 1;

  // Main worker loop: continue as long as wait_for_work indicates readiness
  // wait_for_work handles the blocking and the termination check (via pthread_exit)
  while (wait_for_work(&tp, &run_id) == 0) {
    tp.routine(arg);
  }
  return NULL;
}

/**
 * @brief Pins a thread to a specific logical CPU core.
 *
 * Platform-specific optimization for Linux to improve cache locality
 * and reduce contention by binding threads to specific cores.
 *
 * @param thread The pthread handle.
 * @param core The logical CPU core index to bind the thread to.
 *
 * @note Only effective on Linux systems where `pthread_setaffinity_np` is available.
 *       It's a no-op on other systems.
 * @internal
 * @internal
 */
static void pin_thread_to_cpu(pthread_t thread, int core) {
#ifdef __linux__
  cpu_set_t cpuset;
  CPU_ZERO(&cpuset);
  CPU_SET(core, &cpuset);
  pthread_setaffinity_np(thread, sizeof(cpu_set_t), &cpuset);
#endif
}

void thread_pool_create(thread_pool_t *tp) {
  // Spawn threads
  for (int i = 0; i < MAX_THREADS; i++) {
    tp->thread_ids[i] = i;
    if (pthread_create(&tp->threads[i], NULL, thread_main_wrapper,
                       &tp->thread_ids[i]) != 0) {
      perror("Failed to create thread");
      exit(1);
    }
    // In this context CPUs in cpuset refer to logical cores, therefore the
    // optimal allocation is having one thread per logical CPU
    pin_thread_to_cpu(tp->threads[i], i);
  }
}

/**
 * @brief Waits for all worker threads in the pool to terminate.
 *
 * Calls `pthread_join` for each of the `MAX_THREADS` worker threads, blocking
 * the calling thread (usually the main thread) until they have all exited.
 * This should typically be called after signaling the threads to stop via
 * `thread_pool_terminate`.
 * Exits the program with an error message if joining any thread fails.
 *
 * @param tp Pointer to the thread_pool_t structure containing the thread handles.
 * @internal
 */
void join_threads(thread_pool_t *tp) {
  for (int i = 0; i < MAX_THREADS; ++i) {
    if (pthread_join(tp->threads[i], NULL) != 0) {
      perror("Failed to join thread");
      exit(1);
    }
  }
}

void thread_pool_terminate(thread_pool_t *tp) {
  pthread_mutex_lock(&tp->mutex_children);
  tp->stop_threads = true;
  tp->run_id++;
  pthread_cond_broadcast(&tp->cond_children);
  pthread_mutex_unlock(&tp->mutex_children);
  join_threads(tp);
}

void thread_pool_start_wait(thread_pool_t *tp) {
  pthread_mutex_lock(&tp->mutex_children);
  pthread_mutex_lock(&tp->mutex_parent);
  tp->run_id++;
  tp->children_done = false;

  pthread_cond_broadcast(&tp->cond_children);
  pthread_mutex_unlock(&tp->mutex_children);
  while (!tp->children_done) {
    pthread_cond_wait(&tp->cond_parent, &tp->mutex_parent);
  }
  pthread_mutex_unlock(&tp->mutex_parent);
}

void thread_pool_notify_parent(thread_pool_t *tp) {
  pthread_mutex_lock(&tp->mutex_parent);
  tp->children_done = true;
  pthread_cond_signal(&tp->cond_parent);
  pthread_mutex_unlock(&tp->mutex_parent);
}

void destroy_thread_pool(thread_pool_t *tp) {
  pthread_mutex_destroy(&tp->mutex_children);
  pthread_cond_destroy(&tp->cond_children);
  pthread_mutex_destroy(&tp->mutex_parent);
  pthread_cond_destroy(&tp->cond_parent);
}
