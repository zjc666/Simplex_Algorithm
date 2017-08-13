#include "common.h"
#include "matrix.h"
#include "simplex.h"

MATRIX mat;
MATRIX goal_mat;

void Test1(STRATEGY strategy)
{
	char* var_name[] = {"","x1","x2", "x3", "x4"};
	//	 x1 +     + 2x3       <= 18
	//      + 2x2       - 7x4 <= 0
	//   x1 +  x2 +  x3 +  x4 <= 9			[ x1+x2=b <=> (x1+x2>=b && x1+x2<=b) ]
	// - x1 -  x2 -  x3 -  x4 <= -9
	//        -x2 +  x3 - 2x4 <= -1
	float val[] = 
	{
		18,  1,  0,  2,  0,
		 0,  0,  2,  0, -7,
		 9,  1,  1,  1,  1,
		-9, -1, -1, -1, -1,
		-1,  0,- 1,  1, -2,
	};
	//	maximize U = x1 + x2 + 3x3 - x4
	float goal[] = {0, 1, 1, 3, -1};
	int max;
	int solution[4];
	int i;

	create_matrix(&mat, 5, 5);
	set_matrix(&mat, val);

	create_matrix(&goal_mat, 1, 5);
	set_matrix(&goal_mat, goal);

	initSimplexModel(mat, goal_mat, var_name);
	SimplexRun(strategy);
	max = getMaxIntValue();
	getIntSolution(solution);
	deleteSimplexModel();

	delete_matrix(&mat);
 	delete_matrix(&goal_mat);
}

void Test2(STRATEGY strategy)
{
	u8 i;
	char* var_name[] = { "", "x1", "x2"};
	//  9x1  +  7x2  <=   56
	// 7x1  +  20x2  <=  70
	float val[] =
	{
	    56,   9,  7,
	   70,  7, 20,
	};
	//	maximize U = 40x1 + 90x2
	float goal[] = { 0, 40, 90 };
	int max;
	int solution[2];

	create_matrix(&mat, 2, 3);
	set_matrix(&mat, val);

	create_matrix(&goal_mat, 1, 3);
	set_matrix(&goal_mat, goal);

	initSimplexModel(mat, goal_mat, var_name);
	SimplexRun(strategy);
	max = getMaxIntValue();
	getIntSolution(solution);
	deleteSimplexModel();

	delete_matrix(&mat);
 	delete_matrix(&goal_mat);
}

void initSolver(void);
void solverTest(void);

void main(void)
{
	Test1(LP);
	Test2(ILP);

	printf("\nPress any key to exit");
	getchar();
}