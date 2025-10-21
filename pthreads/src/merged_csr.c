#include "merged_csr.h"
#include <stdlib.h>

MergedCSR *to_merged_csr(const mmio_csr_u32_f32_t *graph) {
  MergedCSR *merged_csr = (MergedCSR *)malloc(sizeof(MergedCSR));

  merged_csr->num_edges = graph->nnz;
  merged_csr->num_vertices = graph->nrows;
  merged_csr->row_ptr =
      (mer_t *)malloc((merged_csr->num_vertices + 1) * sizeof(mer_t));
  merged_csr->merged = (mer_t *)malloc(
      ((merged_csr->num_edges) * DATA_SIZE + (merged_csr->num_vertices) * METADATA_SIZE) *
      sizeof(mer_t));

  for (mer_t i = 0; i < merged_csr->num_vertices; i++) {
    mer_t merged_pos = graph->row_ptr[i] * DATA_SIZE + i * METADATA_SIZE;
    uint32_t degree = graph->row_ptr[i + 1] - graph->row_ptr[i];
    DEGREE(merged_csr, merged_pos) = degree;
    DISTANCE(merged_csr, merged_pos) = UINT32_MAX;
    ID(merged_csr, merged_pos) = i;
    merged_pos += METADATA_SIZE;
    for (mer_t j = graph->row_ptr[i]; j < graph->row_ptr[i + 1]; j++, merged_pos++) {
      // merged_csr->merged[merged_pos] = graph->col_idx[j]; // Original vertex ID
      merged_csr->merged[merged_pos] =
          graph->row_ptr[graph->col_idx[j]] * DATA_SIZE +
          graph->col_idx[j] * METADATA_SIZE; // Vertex position in Merged CSR
    }
  }
  // Create new row_ptr accounting for the shift added by the metadata in the
  // merged CSR
  for (mer_t i = 0; i < merged_csr->num_vertices + 1; i++) {
    merged_csr->row_ptr[i] = graph->row_ptr[i] * DATA_SIZE + i * METADATA_SIZE;
  }
  return merged_csr;
}

void destroy_merged_csr(MergedCSR *merged_csr) {
  free(merged_csr->merged);
  free(merged_csr->row_ptr);
  free(merged_csr);
}
