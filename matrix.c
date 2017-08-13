#include "matrix.h"

u8 create_matrix(MATRIX* mat, int row, int col)
{
	int i;

	mat->mem = (float**)malloc(row * sizeof(float*));
	for(i = 0 ; i < row ; i++){
		mat->mem[i] = (float*)malloc(col * sizeof(float));	
	}

	if(mat->mem == NULL){
		printf("malloc fail\n");
		return 1;
	}
	mat->row = row;
	mat->col = col;

	return 0;
}

void delete_matrix(MATRIX* mat)
{
	int i;

	for(i = 0 ; i<mat->row ; i++)
		free(mat->mem[i]);
	free(mat->mem);
}

void set_matrix(MATRIX* mat, float* val)
{
	int row,col;

	for(row = 0 ; row < mat->row ; row++){
		for(col = 0 ; col < mat->col ; col++){
			mat->mem[row][col] = val[col + row * mat->col];
		}
	}
}

void copy_matrix(MATRIX* des_mat, MATRIX src_mat)
{
	int row, col;

	if( des_mat->row != src_mat.row || des_mat->col != src_mat.col ){
		printf("err, unmatched matrix copy\n");
		return;
	}

	for(row = 0 ; row < des_mat->row ; row++){
		for(col = 0 ; col < des_mat->col ; col++){
			des_mat->mem[row][col] = src_mat.mem[row][col];
		}
	}
}

void matrix_expand_column(MATRIX* mat, u32 exp_num, float* val)
{
	MATRIX temp;
	int row, col;
	u32 cnt = 0;

	create_matrix(&temp, mat->row, mat->col);
	copy_matrix(&temp, *mat);
	/* recreate matrix with expanded column */
	delete_matrix(mat);
	create_matrix(mat, temp.row, temp.col+exp_num);

	for(row = 0 ; row < temp.row ; row++){
		for(col = 0 ; col < temp.col ; col++){
			mat->mem[row][col] = temp.mem[row][col];
		}
		for( ; col < mat->col ; col++){
			mat->mem[row][col] = val[cnt++];
		}
	}

	delete_matrix(&temp);
}

void matrix_expand_row(MATRIX* mat, u32 exp_num, float* val)
{
	MATRIX temp;
	int row, col;
	u32 cnt = 0;

	create_matrix(&temp, mat->row, mat->col);
	copy_matrix(&temp, *mat);
	/* recreate matrix with expanded row */
	delete_matrix(mat);
	create_matrix(mat, temp.row+exp_num, temp.col);

	for(row = 0 ; row < temp.row ; row++){
		for(col = 0 ; col < temp.col ; col++){
			mat->mem[row][col] = temp.mem[row][col];
		}
	}
	for(row = temp.row ; row < mat->row ; row++){
		for(col = 0 ; col < temp.col ; col++){
			mat->mem[row][col] = val[cnt++];
		}
	}

	delete_matrix(&temp);
}

//index start from 0
void matrix_remove_column(MATRIX* mat, int col_index)
{
	MATRIX temp;
	int row, col;
	u32 cnt;

	if(col_index > mat->col-1){
		printf("err, matrix_remove_column\n");
		return;
	}

	create_matrix(&temp, mat->row, mat->col);
	copy_matrix(&temp, *mat);
	/* recreate matrix*/
	delete_matrix(mat);
	create_matrix(mat, temp.row, temp.col-1);

	for(row = 0 ; row < temp.row ; row++){
		cnt = 0;
		for(col = 0 ; col < temp.col ; col++){
			if(col == col_index)
				continue;
			mat->mem[row][cnt++] = temp.mem[row][col];
		}
	}

	delete_matrix(&temp);
}

void print_matrix(MATRIX mat)
{
	int row,col;

	printf("\n");
	for(row = 0 ; row < mat.row ; row++){
		for(col = 0 ; col < mat.col ; col++){
			printf("%.2f\t", mat.mem[row][col]);
		}
		printf("\n");
	}
}
