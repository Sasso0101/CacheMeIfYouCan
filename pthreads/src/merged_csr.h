#ifndef MERGEDCSR_H
#define MERGEDCSR_H

#include "config.h"
#include "mmio_c_wrapper.h"

typedef struct {
  uint32_t num_vertices;
  uint32_t num_edges;
  mer_t *row_ptr;
  mer_t *merged;
} MergedCSR;

#define METADATA_SIZE 3
#define DATA_SIZE 1
#define DEGREE(mer, i) mer->merged[i]
#define DISTANCE(mer, i) mer->merged[i + 1]
#define ID(mer, i) mer->merged[i+2]

/**
 * Converts the CSR graph into a modified merged CSR format with embedded
 * metadata. This layout allows efficient BFS traversal where each vertex's 
 * distance and neighbor information are stored contiguously, improving 
 * cache performance.
 */
MergedCSR *to_merged_csr(const mmio_csr_u32_f32_t *graph); 

void destroy_merged_csr(MergedCSR *merged_csr); 

#endif // MERGEDCSR_H
