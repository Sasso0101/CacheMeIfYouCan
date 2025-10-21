#ifndef CONFIG_H
#define CONFIG_H

#include <stdint.h>

#ifndef MAX_THREADS
#define MAX_THREADS 8
#endif

#ifndef CHUNK_SIZE
#define CHUNK_SIZE 64
#endif
#define INITIAL_CHUNKS_PER_THREAD 128

// Seed used for picking source vertices
// Using same seed as in GAP benchmark for reproducible experiments
// https://github.com/sbeamer/gapbs/blob/b5e3e19c2845f22fb338f4a4bc4b1ccee861d026/src/util.h#L22
#define SEED 27491095

// Selects weather to use 32-bit or 64-bit integers for the vertices in the
// merged CSR
typedef uint32_t mer_t;
#define VERT_MAX UINT32_MAX

#endif // CONFIG_H
