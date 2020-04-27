/*******************************************************************************
 * Utilities and helper functions for SMIOL
 *******************************************************************************/
#ifndef SMIOL_UTILS_H
#define SMIOL_UTILS_H

#include "smiol_types.h"

#define SMIOL_COMP_TO_IO 1
#define SMIOL_IO_TO_COMP 2


/*
 * Searching and sorting
 */
void sort_triplet_array(size_t n_arr, SMIOL_Offset *arr, int sort_entry);
SMIOL_Offset *search_triplet_array(SMIOL_Offset key,
                                   size_t n_arr, SMIOL_Offset *arr,
                                   int search_entry);

/*
 * Communication
 */
int transfer_field(const struct SMIOL_decomp *decomp, int dir,
                   size_t element_size, const void *in_field, void *out_field);

/*
 * Field decomposition
 */
int get_io_elements(int comm_rank, int num_io_tasks, int io_stride,
                    size_t n_io_elements, size_t *io_start, size_t *io_count);

/*
 * Debugging
 */
void print_lists(int comm_rank, SMIOL_Offset *comp_list, SMIOL_Offset *io_list);

#endif
