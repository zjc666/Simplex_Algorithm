#include "simplex.h"
#include <string.h>
#include <stdlib.h>
#include <math.h>

#define MAX_VAL_NUM		6144		/* The maximul variable number */
#define MAX_VAL_NAME	8		/* The maximul variable name */

static MATRIX temp_mat;
static MATRIX goal_mat;
static char VarName[MAX_VAL_NUM][MAX_VAL_NAME];
static u32 VarName_Map[MAX_VAL_NUM] = {0};
static int BVarOff;
static float* coefbuff;
static u8 expand_flag;
static u8 using_name;

static int temp_z_col;

float accuracy = 0.01f;
#define	eq(a,b)		( (a-b)<accuracy && (a-b)>-accuracy )
#define	lt(a, b)	((b-a)>accuracy)
#define	le(a, b)	((b-a)>accuracy || ((b-a)>-accuracy && (b-a)<accuracy))
#define isInt(a)	( fabs((fabs(a)-(int)(fabs(a)+0.5))) < accuracy )
#define toInt(a)	(a>0 ? (int)(a+0.5f) : (int)(a-0.5f))

#define RESERVE_SPACE		3
#define SHOW_RESULT			1

void showResult(MATRIX sf_mat, u8 goal)
{
	int row,col;

	if(using_name){
		for(row = 0 ; row < sf_mat.row ; row++){
			if(goal)
				printf("\nU = ");
			else
				printf("\n%s = ", VarName[VarName_Map[row + BVarOff]]);
			for(col = 0 ; col < sf_mat.col ; col++){
				printf("%.2f%s	", sf_mat.mem[row][col], VarName[VarName_Map[col]]);
			}
		}
		printf("\n");
	}else{
		for(row = 0 ; row < sf_mat.row ; row++){
			if(goal)
				printf("\nU = ");
			for(col = 0 ; col < sf_mat.col ; col++){
				printf("%.2f	", sf_mat.mem[row][col]);
			}
		}
		printf("\n");
	}
}

/************************************************************************
	forall b_i>=0, return 0,   otherwise: return 1                                                              
************************************************************************/
u8 should_expand(MATRIX mat)
{
	int row;

	for(row = 0 ; row < mat.row ; row++){
		/* the first column is b_i */
		if(mat.mem[row][0] < 0)
			return 1;
	}

	return 0;
}

u8 addNonBasicVarName(MATRIX sf_mat, char* name)
{
	if(sf_mat.col >= BVarOff){
		printf("err, non-basic var reserve space is not enough\n");
		return 1;
	}

	if(using_name)
		strcpy(VarName[sf_mat.col], name);
	VarName_Map[sf_mat.col] = sf_mat.col;

	return 0;
}

void resetSimplexModel(MATRIX* sf_mat)
{
	int i;

	//expand_flag = 0;

	if(should_expand(*sf_mat)){
		float* val;

		if(!expand_flag){
			//not expand before, we need allocate name space for z
			addNonBasicVarName(*sf_mat, "z");
		}else{
			//VarName_Map[temp_z_col] = sf_mat->col-1;

			for(i = temp_z_col ; i < sf_mat->col ; i++){
				VarName_Map[i] = VarName_Map[i+1]; 
			}
			VarName_Map[i] = sf_mat->col;
		}

		val = (float*)malloc(sf_mat->row * sizeof(float));
		for(i = 0 ; i < sf_mat->row ; i++)
			val[i] = 1.0f;	/* each row add a fresh variable z, for z>=0 */

		matrix_expand_column(sf_mat, 1, val);
		free(val);

		expand_flag = 1;
	}else{
		expand_flag = 0;
	}
	
	free(coefbuff);
	coefbuff = (float*)malloc(sf_mat->row*sizeof(float));
}

/************************************************************************
Set goal form, calculate the maximum value of goal form.
The last column of mat is constant
************************************************************************/
void setGoalForm(MATRIX mat)
{
	copy_matrix(&goal_mat, mat);
	if(expand_flag){
		/* add -z for goal form */
		float val = -1.0f;
		//float val = 0.0f;
		matrix_expand_column(&goal_mat, 1, &val);
	}
}

/************************************************************************
mat: A * x_i <= b_i (x_i >=0)
the first column of matrix stores constant b_i

Slack Form: b_i - A * x_i >= 0
************************************************************************/
void createSlackForm(MATRIX* mat)
{
	int row, col;

	for(row = 0 ; row < mat->row ; row++)
		for(col = 1 ; col < mat->col ; col++)
			mat->mem[row][col] = -mat->mem[row][col];
}

/************************************************************************
Find a non-basic variable, which has a positive coefficient 
return: -1:doesn't have positive coefficient  otherwise: index of coefficient                                                         
************************************************************************/
int findNextNonBasicVar(void)
{
	int i;

	/* the first column is constant */
	for(i = 1 ; i < goal_mat.col ; i++){
		/* ignore z */
		if( expand_flag && VarName_Map[i] == goal_mat.col-1 )
			continue;
		if( goal_mat.mem[0][i] > accuracy )
  			return i;
	}

	return -1;
}

int findNextNonBasicVarWithZ(MATRIX sf_mat, int row)
{
	int col;

	/* we calculate the min(z) here, so find negative coefficient */
	for(col = 1 ; col < sf_mat.col ; col++){
		//TODO: need goal_mat.mem[0][col]>0?
		//if( lt(sf_mat.mem[row][col], 0.0f) && lt(0.0f, goal_mat.mem[0][col]) ){
		if( lt(sf_mat.mem[row][col], 0.0f) ){
			return col;
		}
	}

	return -1;
}

/************************************************************************
     find most negative b_i to be pivot with z                                                                
************************************************************************/
PIVOT findMostNegativePivot(MATRIX sf_mat)
{
	int row;
	float neg_b;
	PIVOT pivot;

	pivot.bvar = -1;	
	pivot.nbvar = sf_mat.col-1;	 /* swap z and basic variable */
	
	neg_b = -accuracy;
	for(row = 0 ; row < sf_mat.row ; row++){
		if(sf_mat.mem[row][0] < neg_b){
			neg_b = sf_mat.mem[row][0];
			pivot.bvar = row;
		}
	}

	return pivot;
}

PIVOT findPivot(MATRIX sf_mat, int nbvar)
{
	int row;
	float tight_constrain = 100000000.0f;
	float temp_constrain;
	PIVOT pivot;

	pivot.bvar = -1;	/* -1 means doesn't find a valid pivot */
	pivot.nbvar = nbvar;
	
	for(row = 0 ; row < sf_mat.row ; row++){
		if( eq(0.0f, sf_mat.mem[row][nbvar]) || sf_mat.mem[row][nbvar] > 0.0f )
			continue;
		temp_constrain = sf_mat.mem[row][0] / -sf_mat.mem[row][nbvar];
		if( temp_constrain < tight_constrain ){
			tight_constrain = temp_constrain;
			pivot.bvar = row;
		}
	}

	return pivot;
}

void swapNameSpace(PIVOT pivot)
{
	int map_index;

	map_index = VarName_Map[pivot.nbvar];
	VarName_Map[pivot.nbvar] = VarName_Map[pivot.bvar+BVarOff];
	VarName_Map[pivot.bvar+BVarOff] = map_index;
}

/************************************************************************
        swap basic var and non-basic var                                                              
************************************************************************/
void swapPivot(MATRIX* sf_mat, PIVOT pivot)
{
	int col;
	float temp;

	temp = sf_mat->mem[pivot.bvar][pivot.nbvar];
	for(col = 0 ; col < sf_mat->col ; col++){
		if(col == pivot.nbvar)
			sf_mat->mem[pivot.bvar][col] = 1.0f / temp;
		else
			sf_mat->mem[pivot.bvar][col] /= -temp;
	}
}

void replaceNewBasicVar(MATRIX* sf_mat, PIVOT pivot)
{
	int row, col;

	for(row = 0 ; row < sf_mat->row ; row++){
		if(row == pivot.bvar)
			continue;

		coefbuff[row] = sf_mat->mem[row][pivot.nbvar];
		
		for(col = 0 ; col < sf_mat->col ; col++){
			if(col == pivot.nbvar)
				sf_mat->mem[row][col] = coefbuff[row] * sf_mat->mem[pivot.bvar][col];
			else
				sf_mat->mem[row][col] += coefbuff[row] * sf_mat->mem[pivot.bvar][col];
		}
	}
}

void updateGoalForm(MATRIX* gf_mat, MATRIX* sf_mat, PIVOT pivot)
{
	int col;
	float temp;
	
	temp = gf_mat->mem[0][pivot.nbvar];

	for(col = 0 ; col < gf_mat->col ; col++){
		if(col == pivot.nbvar)
			gf_mat->mem[0][col] = temp * sf_mat->mem[pivot.bvar][col];
		else
			gf_mat->mem[0][col] += temp * sf_mat->mem[pivot.bvar][col];
	}
}

/************************************************************************
swap pivot and reformulate slack form                                                               
************************************************************************/
void reformulation(MATRIX* sf_mat, PIVOT pivot)
{
	/* step1: swap basic var and non-basic var */
	swapNameSpace(pivot);
	swapPivot(sf_mat, pivot);
	/* step2: replace new basic var with non-basic var */
	replaceNewBasicVar(sf_mat, pivot);
	/* step3: update goal form */
	updateGoalForm(&goal_mat, sf_mat, pivot);
}

void dispSolution(MATRIX sf_mat, u8 show_res)
{
	int row, col;

	if(show_res){
		showResult(goal_mat, 1);
		showResult(temp_mat, 0);
		printf("\nSolution:\n");
	}
	for(col = 1 ; col < sf_mat.col ; col++){
		float solution = 0.0f;

		for(row = 0 ; row < sf_mat.row ; row++){
			if(VarName_Map[row + BVarOff] == col){
				solution = sf_mat.mem[row][0];
			}
		}
		if(show_res){
			if(using_name)
				printf("%s=%.2f, ", VarName[col], solution);
			else
				printf("%.2f, ", solution);
		}
	}
	if(show_res)
		printf("U=%.2f\n", goal_mat.mem[0][0]);
}

void dispIntSolution(MATRIX sf_mat, u8 show_res)
{
	int row, col;

	if(show_res){
		showResult(goal_mat, 1);
		showResult(temp_mat, 0);
		printf("\nInterger Solution:\n");
	}

	for(col = 1 ; col < temp_mat.col-expand_flag ; col++){
		int sol = 0;

		for(row = 0 ; row < temp_mat.row ; row++){
			if(VarName_Map[row + BVarOff] == col){	
				sol = toInt(temp_mat.mem[row][0]);
			}
		}
		if(show_res){
			if(using_name)
				printf("%s=%d, ", VarName[col], sol);
			else
				printf("%d, ", sol);
		}
	}

	if(show_res)
		printf("U=%.2f\n", goal_mat.mem[0][0]);
}

int getNonIntIndex(MATRIX sf_mat, int row_index)
{
	int row, col;

	if(row_index >= sf_mat.row){
		return -1;
	}

	for(row = (int)row_index ; row < sf_mat.row ; row++){
		if( !isInt(sf_mat.mem[row][0]) ){
			for(col = 1 ; col < sf_mat.col ; col++){
				if( !isInt(sf_mat.mem[row][col]) )
					return row;	/* find a basic var which has non-interger bi and xi */
			}
		}
	}

	return -100;
}

/* add gomory constraint to calculate interger solution */
void addGomoryConstraint(MATRIX* s_mat,  int row_index)
{
	static int cons_cnt = 1;
	float* gomory_cons;
	int col;

	gomory_cons = (float*)malloc(s_mat->col * sizeof(float));

	/* extract basic var fractional part, with b>=0 */
	gomory_cons[0] = -(s_mat->mem[row_index][0] - (int)s_mat->mem[row_index][0]);

	/* extract non-basic var fractional part */
	for(col = 1 ; col < s_mat->col ; col++){
		if( eq(s_mat->mem[row_index][col], 0.0f) || isInt(s_mat->mem[row_index][col]) ){
			/* no fraction part */
			gomory_cons[col] = 0.0f;
			continue;
		}
		gomory_cons[col] = -(s_mat->mem[row_index][col] - (float)ceil(s_mat->mem[row_index][col]));
	}

	/* add gomory constraint */
	matrix_expand_row(s_mat, 1, gomory_cons);

	if(using_name)
		sprintf(VarName[s_mat->row-1+BVarOff], "g%d", cons_cnt++);
	VarName_Map[s_mat->row-1+BVarOff] = s_mat->row-1+BVarOff;

	free(gomory_cons);
}

/* remove z column */
int removeZ(MATRIX* s_mat)
{
	int z_col = -1;
	if(!expand_flag)
		return -1;

	for(z_col = 0 ; z_col < s_mat->col ; z_col++){
		if(VarName_Map[z_col] == s_mat->col-1){
			break;
		}
	}

	if(z_col >= 0){
		matrix_remove_column(s_mat, (int)z_col);
		matrix_remove_column(&goal_mat, (int)z_col);
	}

	return z_col;
}

u8 solveSimplex(void)
{
	u8 res = 0;
	int next_nbvar;
	PIVOT pivot;

	/* we should first calculate the min(z), if min(z)=0, then change z to non-basic var */
	if(expand_flag){
		int row_z;

		pivot = findMostNegativePivot(temp_mat);
		row_z = (int)pivot.bvar;
			
		if(pivot.bvar == -1){
			//printf("NA, findMostNegativePivot\n");
			return 0;
		}
		reformulation(&temp_mat, pivot);

		while(1){
			if( eq(temp_mat.mem[row_z][0], 0.0f) ){
				/* min(z)=0 */
				int col;
				u8 res = 0;
				for(col = 1 ; col < temp_mat.col ; col++ ){
					if( !eq(temp_mat.mem[row_z][col], 0.0f)){
						pivot.bvar = row_z;
						pivot.nbvar = col;
						reformulation(&temp_mat, pivot);
						res = 1;
						break;
					}
				}
				if(res)
					break;
				else{
					printf("some err here\n");
					while(1);
				}
				
			}

			next_nbvar = findNextNonBasicVarWithZ(temp_mat, row_z);
			if(next_nbvar == -1){
				//printf("NA, findNextNonBasicVarForZ\n");
				//showResult(goal_mat, 1);
				//showResult(temp_mat, 0);
				return 0;
			}
			pivot = findPivot(temp_mat, (int)next_nbvar);
			reformulation(&temp_mat, pivot);
			if(pivot.bvar == row_z){
				//TODO: is it correct?
				/* z becomes non-basic variable now, we can continue */
				break;
			}
		}
	}

	while(1){
		next_nbvar = findNextNonBasicVar();
		if( next_nbvar < 0 ){
			/* non positive coefficient, we get the maximum value */
			res = 1;
			break;
		}

		pivot = findPivot(temp_mat, (int)next_nbvar);
		if(pivot.bvar == -1){
			/* goal formula can be infinite, so doesn't have maximum value */
			//printf("NA, inf result\n");
			res = 0;
			break;
		}else{
			next_nbvar = 1;
		}
		reformulation(&temp_mat, pivot);
	}
	
	return res;
}

u8 incise(void)
{
	int nonIntIndex = -1;
	u8 res;

	while(1){
		//nonIntIndex = getNonIntIndex(temp_mat, nonIntIndex+1);
		nonIntIndex = getNonIntIndex(temp_mat, 0);
		/* Handle ILP */
		if(nonIntIndex >= 0){	
			temp_z_col = removeZ(&temp_mat);
			addGomoryConstraint(&temp_mat, (int)nonIntIndex);
			resetSimplexModel(&temp_mat);
			setGoalForm(goal_mat);

			res = solveSimplex();
			if(!res)	//NA
				return res;
		}else if(nonIntIndex == -1){
			//TODO: is it possible for this situation?
			/* used up all gomory constraints and can't find interger solution */
			//printf("NA, can't find interger solution\n");
			return 0;
		}else{
			/* Get interger soluiton */
			return 1;
		}
	}
}

int getMaxIntValue(void)
{
	return toInt(goal_mat.mem[0][0]);
}

void getIntSolution(int* solution)
{
	int row, col;

	for(col = 1 ; col < temp_mat.col-expand_flag ; col++){
		float sol = 0.0f;

		for(row = 0 ; row < temp_mat.row ; row++){
			if(VarName_Map[row + BVarOff] == col){	
				sol = temp_mat.mem[row][0];
			}
		}
		solution[col-1] = toInt(sol);
	}
}

u8 initSimplexModel(MATRIX mat, MATRIX g_mat, char* var_name[])
{
	u8 res = 0;
	int i;
	
	if(var_name == NULL)
		using_name = 0;
	else
		using_name = 1;
	expand_flag = 0;

	res |= create_matrix(&goal_mat, 1, mat.col);
	/* we will do calculation on temp_mat later */
	res |= create_matrix(&temp_mat, mat.row, mat.col);
	copy_matrix(&temp_mat, mat);

	if(should_expand(temp_mat)){
		float* val;
		int i;
		expand_flag = 1;

		val = (float*)malloc(temp_mat.row * sizeof(float));
		for(i = 0 ; i < temp_mat.row ; i++)
			val[i] = -1.0f;	/* each row add a fresh variable z, for z>=0 */

		matrix_expand_column(&temp_mat, 1, val);
		free(val);
	}

	coefbuff = (float*)malloc(mat.row*sizeof(float));

	for(i = 0 ; i < mat.col ; i++){
		if(using_name)
			strcpy(VarName[i], var_name[i]);
		VarName_Map[i] = i;
	}
	if(expand_flag){
		if(using_name)
			strcpy(VarName[mat.col], "z");
		VarName_Map[mat.col] = mat.col;
		BVarOff = mat.col + RESERVE_SPACE;
	}else{
		BVarOff = mat.col  + RESERVE_SPACE;	/* basic value offset of VarName array */
	}
	for(i = 0 ; i < mat.row ; i++){
		if(using_name)
			sprintf(VarName[i+BVarOff], "y%d", i+1);
		VarName_Map[i+BVarOff] = i+BVarOff;
	}

	setGoalForm(g_mat);
	createSlackForm(&temp_mat);

	return res;
}

void deleteSimplexModel(void)
{
	delete_matrix(&goal_mat);
	delete_matrix(&temp_mat);
	free(coefbuff);
}

/***************************************************************************
Input:  
	c_mat: constraint matrix   
	g_mat: goal matrix   
	var_name: variable name  
	strategy:
				LP:	Linear Program
			   ILP:	Interger Linear Program
Return:
	0: unsatisfy	otherwise: satisfy
***************************************************************************/
u8 SimplexRun(STRATEGY strategy)
{
	u8 res = 0;

	res = solveSimplex();

	if(res){
		if(strategy == LP){
			dispSolution(temp_mat, SHOW_RESULT);
		}else{
			res = incise();
			dispIntSolution(temp_mat, SHOW_RESULT);
		}
	}

	return res;
}
