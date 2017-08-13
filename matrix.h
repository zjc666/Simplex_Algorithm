#ifndef __MATRIX__
#define __MATRIX__

#include "common.h"

typedef struct {
	int		row,col;
	float**	mem;
}MATRIX;

u8 create_matrix(MATRIX* mat, int row, int col);
void delete_matrix(MATRIX* mat);
void set_matrix(MATRIX* mat, float* val);
void copy_matrix(MATRIX* des_mat, MATRIX src_mat);
void matrix_expand_row(MATRIX* mat, u32 exp_num, float* val);
void matrix_expand_column(MATRIX* mat, u32 exp_num, float* val);
void matrix_remove_column(MATRIX* mat, int col_index);
void print_matrix(MATRIX mat);

#endif