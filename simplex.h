#ifndef __SIMPLEX__
#define __SIMPLEX__

#include "common.h"
#include "matrix.h"

typedef struct {
	int bvar;	//basic variable: yi
	int nbvar;	//non-basic variable: xi
}PIVOT;

typedef enum {
	LP = 0,
	ILP
}STRATEGY;

u8 initSimplexModel(MATRIX mat, MATRIX g_mat, char* var_name[]);
void deleteSimplexModel(void);
u8 SimplexRun(STRATEGY strategy);
int getMaxIntValue(void);
void getIntSolution(int* solution);

#endif