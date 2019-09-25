
#include "matrix.h"
#include <stdio.h>
#include <iostream>
#include <fstream>

#pragma region vector_operations
inline void vectorAdd(std::vector<double>& a, std::vector<double>& b, std::vector<double>& result);
inline void vectorSubtract(std::vector<double>& a, std::vector<double>& b, std::vector<double>& result);
inline void vectorAssigment(std::vector<double>& a, std::vector<double>& b);
inline void vectorMultOnConst(std::vector<double>& a, double constant, std::vector<double>& result);
inline double scalarProduct(std::vector<double>& a, std::vector<double>& b);
inline double calcNorm(std::vector<double>& vector);
#pragma endregion

#pragma region  someUtilits
SparseMatrix::SparseMatrix() = default;
int SparseMatrix::getDim() { return n; }
int SparseMatrix::getAmountElems() { return ggl.size(); }
int SparseMatrix::getMethod() { return method; }

void SparseMatrix::copy(SparseMatrix source) {
	ig = source.ig;				// index elems - start of row
	jg = source.jg;				// column index of elem
	di = source.di;				// diag elems
	ggu = source.ggu;			// elem's of upper trinagle
	ggl = source.ggl;			// elem's of lower trinagle

	n = source.n;						// dimension
	max_iter = source.max_iter;			// max amount of iter
	eps = source.eps;					// min relative discrepancy
	method = source.method;
}

void SparseMatrix::load() {
	FILE *in;
	fopen_s(&in, "kuslau.txt", "r");
	fscanf_s(in, "%d %d %le", &n, &max_iter, &eps);
	fclose(in);

	fopen_s(&in, "method.txt", "r");
	fscanf_s(in, "%d", &method);
	fclose(in);

	di.resize(n);
	ig.resize(n + 1);

	fopen_s(&in, "ig.txt", "r");
	for (int i = 0; i < n + 1; i++) fscanf_s(in, "%d", &ig[i]);
	fclose(in);

	fopen_s(&in, "di.txt", "r");
	for (int i = 0; i < n; i++) fscanf_s(in, "%le", &di[i]);
	fclose(in);

	jg.resize(ig[n]);
	ggu.resize(ig[n]);
	ggl.resize(ig[n]);

	fopen_s(&in, "jg.txt", "r");
	for (int i = 0; i < jg.size(); i++) {
		fscanf_s(in, "%d", &jg[i]);
	}
	fclose(in);

	fopen_s(&in, "ggu.txt", "r");
	for (int i = 0; i < ggu.size(); i++) fscanf_s(in, "%le", &ggu[i]);
	fclose(in);

	fopen_s(&in, "ggl.txt", "r");
	for (int i = 0; i < ggl.size(); i++) fscanf_s(in, "%le", &ggl[i]);
	fclose(in);

	if (ig[0] == 1) {	// transform FORTRAN format --> C Format.
		for (int i = 0; i < ig.size(); i++) ig[i]--;
		for (int i = 0; i < jg.size(); i++) jg[i]--;
	}
}
#pragma endregion

void SparseMatrix::resetMatrix() {
	for (auto &el : ggl) el = 0;
	for (auto &el : ggu) el = 0;
	for (auto &el : di) el = 0;
}

void SparseMatrix::init(const int elemNum, const int _max_iter, const double _eps) {
	// 6 elems of first finite elem
	// 5 elements for each next finite elem
	int amountElemsInMatrix = elemNum * 5 + 1; // 6 + (num - 1)*5 = num*5 + 1:

	// 4 dim for first finite elem
	// 2 dim for each next finite elem
	n = elemNum * 2 + 2;	// 4 + (num - 1)*2 = num*2 + 2:
	max_iter = _max_iter;
	eps = _eps;

	// memory allocation
	di = std::vector<double>(n);
	ig = std::vector<int>(n + 1);
	jg = std::vector<int>(amountElemsInMatrix);
	ggu = std::vector<double>(amountElemsInMatrix);
	ggl = std::vector<double>(amountElemsInMatrix);

	// forming ig array
	ig[0] = 0;
	ig[1] = 0;
	ig[2] = 1;
	for (int pos = 3; pos < n; pos += 2) {
		ig[pos] = ig[pos - 1] + 2;
		ig[pos + 1] = ig[pos] + 3;
	}

	// formin jg array
	jg[0] = 0;
	for (int pos = 1, i = 0, column = 0; pos < amountElemsInMatrix; pos++) {
		jg[pos] = column;				// columns:
		jg[++pos] = column + 1;			// j j+1
		jg[++pos] = column;				// j j+1 j+2
		jg[++pos] = ++column;
		jg[++pos] = ++column;
	}
}

// need test it
void SparseMatrix::init(const int dim, std::vector<int>& _ig, std::vector<int>& _jg) {


	max_iter = 40000;
	eps = 1e-13;

	n = dim;

	di = std::vector<double>(n);
	ig = _ig;
	jg = _jg;

	ggu = std::vector<double>(jg.size());
	ggl = std::vector<double>(jg.size());
}

void SparseMatrix::addDiagElem(int n, double value) {
	di[n] += value;
}

void SparseMatrix::addElem(int row, int column, double value) {
	int index = ig[row];
	int maxIndex = ig[row + 1];

	for (; index < maxIndex; index++) {
		if (jg[index] == column) {
			ggl[index] += value;
			return;
		}
	}

	throw "Not found column";
}

void SparseMatrix::fillGGU() {
	ggu = ggl;
}

void SparseMatrix::addLocalMatrix(int elemNum, LocalMatrix M) {
	// consider diag elems
	const int diagPos = elemNum * 2;
	di[diagPos] += M[0][0];
	di[diagPos + 1] += M[1][1];
	di[diagPos + 2] += M[2][2];
	di[diagPos + 3] += M[3][3];

	const int elemPos = elemNum * 5;
	ggl[elemPos] += M[1][0];
	ggu[elemPos] += M[0][1];

	ggl[elemPos + 1] += M[2][0];
	ggu[elemPos + 1] += M[0][2];

	ggl[elemPos + 2] += M[2][1];
	ggu[elemPos + 2] += M[1][2];

	ggl[elemPos + 3] += M[3][0];
	ggu[elemPos + 3] += M[0][3];

	ggl[elemPos + 4] += M[3][1];
	ggu[elemPos + 4] += M[1][3];

	ggl[elemPos + 5] += M[3][2];
	ggu[elemPos + 5] += M[2][3];
}

#pragma region multiplications

void SparseMatrix::multiplicateWithVector(std::vector<double>& b,
	std::vector<double>& result) {
	for (int i = 0; i < result.size(); i++) result[i] = di[i] * b[i];	// init result vect

	for (int i = 0; i < ig.size() - 1; i++) {								// multiplicate:
		for (int j = ig[i]; j < ig[i + 1]; j++) {
			result[i] += ggl[j] * b[jg[j]];
			result[jg[j]] += ggu[j] * b[i];
		}
	}
}

void SparseMatrix::multiplicateTransWithVector(std::vector<double>& b,
	std::vector<double>& result) {
	for (int i = 0; i < result.size(); i++) result[i] = di[i] * b[i];		// init result vect

	for (int i = 0; i < ig.size() - 1; i++) {								// multiplicate:
		for (int j = ig[i]; j < ig[i + 1]; j++) {
			result[i] += ggu[j] * b[jg[j]];
			result[jg[j]] += ggl[j] * b[i];
		}
	}
}

void SparseMatrix::multiplicateDiagWithVector(
	std::vector<double>& b, std::vector<double>& result) {
	for (int i = 0; i < di.size(); i++) result[i] = di[i] * b[i];
}

void  SparseMatrix::
multiplicateUfactor(std::vector<double>& u, std::vector<double>& diag,
	std::vector<double>& b, std::vector<double>& result) {

	for (int i = 0; i < result.size(); i++) result[i] = diag[i] * b[i];	// diag elems

	for (int i = 0; i < ig.size() - 1; i++) {		// multiplicate:		
		for (int j = ig[i]; j < ig[i + 1]; j++) {	// U elems
			result[jg[j]] += u[j] * b[i];
		}
	}
}

#pragma endregion

#pragma region solve_methods

void SparseMatrix::directSolveLU(std::vector<double>& l, std::vector<double>& ud,
	std::vector<double>& u, std::vector<double>& temp, std::vector<double>& q, std::vector<double>& f) {
	// note:
	// this method work only if LU partial decompose equalent to exact LU decompose!
	LUpartialDecomposite(l, ud, u);
	forwardGauss(l, temp, f);
	backwardGauss(u, ud, q, temp);
}


void SparseMatrix::forwardGauss(std::vector<double>& matrix,
	std::vector<double>& x, std::vector<double>& f) {
	int i, j;
	double sum;

	//x[i] = f[i] - sum[0;i-1](Lki * Xi);

	for (i = 0; i < n; i++) {
		sum = 0;
		for (j = ig[i]; j < ig[i + 1]; j++) {
			sum += matrix[j] * x[jg[j]];	// 
		}
		x[i] = f[i] - sum;
	}
}

void SparseMatrix::forwardGauss(std::vector<double>& matrix,
	std::vector<double>& diag, std::vector<double>& x, std::vector<double>& f) {
	int i, j;
	double sum;

	//x[i] = f[i] - sum[0;i-1](Lki * Xi);

	for (i = 0; i < n; i++) {
		sum = 0;
		for (j = ig[i]; j < ig[i + 1]; j++) {
			sum += matrix[j] * x[jg[j]];	// 
		}
		x[i] = (f[i] - sum) / diag[i];
	}
}

void SparseMatrix::backwardGauss(std::vector<double>& matrix,
	std::vector<double>& x, std::vector<double>& f) {

	// method will change vector f
	for (int i = n - 1; i >= 0; i--) {
		x[i] = f[i];

		for (int j = ig[i]; j < ig[i + 1]; j++) {
			f[jg[j]] -= matrix[j] * x[i];
		}
	}
}

void SparseMatrix::backwardGauss(std::vector<double>& matrix,
	std::vector<double>& diag, std::vector<double>& x, std::vector<double>& f) {
	// method will change vector f
	for (int i = n - 1; i >= 0; i--) {
		x[i] = f[i] / diag[i];

		for (int j = ig[i]; j < ig[i + 1]; j++) {
			f[jg[j]] -= matrix[j] * x[i];
		}
	}
}

void SparseMatrix::diagGauss(std::vector<double>& diag, std::vector<double>& x,
	std::vector<double>& f) {
	for (int i = 0; i < diag.size(); i++) x[i] = f[i] / diag[i];
}

void SparseMatrix::diagDecomposite(std::vector<double>& diag) {
	for (int i = 0; i < diag.size(); i++) {
		if (diag[i] >= 0)
			diag[i] = sqrt(di[i]);
		else {
			printf("error, diag elem < 0");
			return;
		}
	}
}

void SparseMatrix::LUpartialDecomposite(std::vector<double>& l,
	std::vector<double>& ud, std::vector<double>& u) {

	int column, row;
	int to_column, to_row;
	int first_column, first_row;
	int amount_column, amount_row;

	double diag_sum, l_sum, u_sum;
	int i, j, k, t;

	for (i = 0; i < n; i++) {
		diag_sum = 0;

		first_column = ig[i];
		//to_column = 0;

		for (j = first_column, amount_column = 0; j < ig[i + 1]; j++, amount_column++) {
			l_sum = 0;
			u_sum = 0;

			column = jg[j];
			amount_row = ig[column + 1] - ig[column];
			first_row = ig[column];

			// find elems to multiplicate (row == column)
			to_column = to_row = 0;

			// while don't yet reviewed all the columns and rows and
			//	current row < max column
			while (to_column < amount_column && to_row < amount_row &&
				jg[first_row + to_row] < jg[j]) {	// j = first_column + amount_column

				if (jg[first_row + to_row] < jg[first_column + to_column]) to_row++;
				if (jg[first_row + to_row] > jg[first_column + to_column]) to_column++;
				if (jg[first_row + to_row] == jg[first_column + to_column]) {
					l_sum += l[first_column + to_column] * u[first_row + to_row];			// accumulate sum for calc L elems
					u_sum += u[first_column + to_column] * l[first_row + to_row];			// accumulate sum for calc U elems

					to_column++;
					to_row++;
				}
			}

			//printf("calc l elem: column: %d\n", column);
			l[j] = (ggl[j] - l_sum) / ud[column];		// L elem
			u[j] = ggu[j] - u_sum;						// U elem

			diag_sum += l[j] * u[j];
		}

		// Chech diag elem:

		if (abs(ud[i] - di[i] + diag_sum) < ud[i] * 1e-15) {
			printf_s("error in LDU decomposition. zero on diag.\n");
			return;
		}
		else {
			ud[i] = di[i] - diag_sum;						// Udiag elem
		}
	}
}

#include <time.h>
void SparseMatrix::LOS(std::vector<double>& x,
	std::vector<double>& f, std::vector<double>& r, std::vector<double>& z,
	std::vector<double>& p, std::vector<double>& temp) {

	std::ofstream solver("solver.txt");

	double a, b, scalar_pp, scalar_rr, discrepancy = 1e40, norm_f;
	int iter = 0;

	// r0 = f - Ax0:
	multiplicateWithVector(x, temp);		// Ax0
	vectorSubtract(f, temp, r);			// r0 = f - Ax0

	z = r;									// z0 = r0
	multiplicateWithVector(z, p);			// p0 = Az0	

	norm_f = calcNorm(f);
	scalar_rr = scalarProduct(r, r);

	double CalcTimeStart = clock() / (double)CLOCKS_PER_SEC;
	// iterations:
	do {
		if (iter > max_iter) {
			solver << "max iter exit: discrepance == " << discrepancy << std::endl;
			break;
		}
		//scalar_rr = scalar_product(r, r);
		scalar_pp = scalarProduct(p, p);		// scalar_pp = (pk-1,pk-1)
		a = scalarProduct(p, r) / scalar_pp; // ak = (pk-1, rk-1)/(pk-1,pk-1)

		vectorMultOnConst(z, a, temp);		// temp = a*zk-1
		vectorAdd(x, temp, x);				// x0 = x0 + a*zk-1

		// step out:
		if (calcNorm(temp) / calcNorm(x) < step_eps) {
#ifdef CALCULATE
			printf("\nstep out\n");
			solver << "step out on " << iter << " with  discrepance == " << discrepancy << std::endl;
#endif
			break;
		}

		vectorMultOnConst(p, a, temp);		// temp = a*pk-1
		vectorSubtract(r, temp, r);			//  rk = rk-1 - a*pk-1

		multiplicateWithVector(r, temp);		// temp = Ark
		b = -1.0 * scalarProduct(p, temp) / scalar_pp;	// bk = -(pk-1, Ark)/(pk-1,pk-1)

		vectorMultOnConst(z, b, temp);		// temp = b*zk-1
		vectorAdd(r, temp, z);				// zk  = zk-1  + b*zk-1

		// pk = Ark + Bkpk-1 <==> 1. pk = Bkpk-1; 2. pk += Ark
		vectorMultOnConst(p, b, p);
		multiplicateWithVector(r, temp);		// temp = Ark
		vectorAdd(p, temp, p);


		// check discrepancy
		scalar_rr -= a * a * scalar_pp;
		if (scalar_rr < 0) scalar_rr = scalarProduct(r, r);

		discrepancy = sqrt(scalar_rr) / norm_f;

		iter++;

	} while (discrepancy > eps);
	
	solver << "good solve, discrepance == " << discrepancy << std::endl;


	double CalcTimeStop = clock() / (double)CLOCKS_PER_SEC;
	//printf("\nEx Time: %f sec\n", (CalcTimeStop - CalcTimeStart));
}

void SparseMatrix::
LOS_LU(std::vector<double>& l, std::vector<double>& ud, std::vector<double>& u,
	std::vector<double>& x, std::vector<double>& f,
	std::vector<double>& r, std::vector<double>& z,
	std::vector<double>& p, std::vector<double>& temp,
	std::vector<double>& temp2) {

	int iter = 0;
	double scalar_pp, scalar_rr, discrepancy, a, b, norm_f;

	// partial LU decomposition	
	LUpartialDecomposite(l, ud, u);

	//r = ..
	multiplicateWithVector(x, temp);
	vectorSubtract(f, temp, temp2);
	forwardGauss(l, r, temp2);

	// z = ..
	vectorAssigment(temp, r);
	backwardGauss(u, ud, z, temp);

	// p = ...
	multiplicateWithVector(z, temp);
	forwardGauss(l, p, temp);

	norm_f = calcNorm(f);

	scalar_rr = scalarProduct(r, r);
	do {
		if (iter > max_iter) break;

		scalar_pp = scalarProduct(p, p);		// scalar_pp = (pk-1,pk-1)

		a = scalarProduct(p, r) / scalar_pp;

		// x = ...
		vectorMultOnConst(z, a, temp);
		vectorAdd(x, temp, x);

		// step out:
		if (calcNorm(temp) / calcNorm(x) < step_eps) {
#ifdef CALCULATE
			printf("\nstep out\n");
#endif
			break;
		}

		//r = ...
		vectorMultOnConst(p, a, temp);
		vectorSubtract(r, temp, r);

		// b = ...
		vectorAssigment(temp, r);
		backwardGauss(u, ud, temp2, temp);
		multiplicateWithVector(temp2, temp);
		forwardGauss(l, temp2, temp);

		b = scalarProduct(p, temp2) / scalar_pp;

		// p = ...
		vectorMultOnConst(p, b, temp);
		vectorAdd(temp2, temp, p);

		// z = ...
		vectorAssigment(temp, r);
		backwardGauss(u, ud, temp2, temp);
		vectorMultOnConst(z, b, temp);
		vectorAdd(temp, temp2, z);

		// check discrepancy
		scalar_rr -= a * a * scalar_pp;
		//scalar_rr = scalar_product(r, r);

		if (scalar_rr < 0) scalar_rr = scalarProduct(r, r);
		discrepancy = sqrt(scalar_rr) / norm_f;
		//discrepancy = scalar_rr / norm_f;

		//for (int i = 0; i < x0.size(); i++) printf("%.15lf\n", x0[i]);

		iter++;
#ifdef CALCULATE
		printf("\riter: %d\t discr: %.15e", iter, discrepancy);
#endif
	} while (discrepancy > eps);
}


void SparseMatrix::
LOS_DIAG(std::vector<double>& diag, std::vector<double>& x, std::vector<double>& f,
	std::vector<double>& r, std::vector<double>& z,
	std::vector<double>& p, std::vector<double>& temp,
	std::vector<double>& temp2) {

	int iter = 0;
	double scalar_pp, scalar_rr, discrepancy, a, b, norm_f;

	// L = U = sqrt(di)
	diagDecomposite(diag);

	//r = ..
	multiplicateWithVector(x, temp);
	vectorSubtract(f, temp, temp2);

	diagGauss(diag, r, temp2);

	// z = ..
	vectorAssigment(temp, r);
	diagGauss(diag, z, temp);

	// p = ...
	multiplicateWithVector(z, temp);
	diagGauss(diag, p, temp);

	norm_f = calcNorm(f);
	scalar_rr = scalarProduct(r, r);

	double CalcTimeStart = clock() / (double)CLOCKS_PER_SEC;
	do {
		if (iter > max_iter) break;

		scalar_pp = scalarProduct(p, p);		// scalar_pp = (pk-1,pk-1)
		a = scalarProduct(p, r) / scalar_pp;

		// x = ...
		vectorMultOnConst(z, a, temp);
		vectorAdd(x, temp, x);

		// step out:
		if (calcNorm(temp) / calcNorm(x) < step_eps) {
#ifdef CALCULATE
			printf("\nstep out\n");
#endif
			break;
		}


		//r = ...
		vectorMultOnConst(p, a, temp);
		vectorSubtract(r, temp, r);

		// b = ...
		vectorAssigment(temp, r);
		diagGauss(diag, temp2, temp);
		multiplicateWithVector(temp2, temp);
		diagGauss(diag, temp2, temp);

		b = scalarProduct(p, temp2) / scalar_pp;

		// p = ...
		vectorMultOnConst(p, b, temp);
		vectorAdd(temp2, temp, p);

		// z = ...
		vectorAssigment(temp, r);
		diagGauss(diag, temp2, temp);
		vectorMultOnConst(z, b, temp);
		vectorAdd(temp, temp2, z);

		// check discrepancy
		scalar_rr -= a * a * scalar_pp;

		if (scalar_rr < 0) printf("!"), scalar_rr = scalarProduct(r, r);
		discrepancy = sqrt(scalar_rr) / norm_f;

		//discrepancy = scalar_rr / norm_f;


		iter++;
#ifdef CALCULATE
		printf("\riter: %d\t discr: %.15e", iter, discrepancy);
#endif
	} while (discrepancy > eps);

	double CalcTimeStop = clock() / (double)CLOCKS_PER_SEC;
	printf("\nEx Time: %f sec\n", (CalcTimeStop - CalcTimeStart));
}


int SparseMatrix::MCG(std::vector<double>& x, std::vector<double>& x_min,
	std::vector<double>& f, std::vector<double>& r, std::vector<double>& z,
	std::vector<double>& temp, std::vector<double> temp2) {

	std::ofstream solver("solver.txt");

	double a, b, scalar_rr, scalar_rr_prev, norm_f, discrepancy, min_discrepancy;
	int iter = 0;
	// system symmetrization:
//vector_assigment(temp, f);				// temp = f	
//multiplicate_trans_with_vector(temp, f);	// f = At * f

	multiplicateWithVector(x, temp);			// r0 = At * (f -  A * x0)
	vectorSubtract(f, temp, temp);
	multiplicateTransWithVector(temp, r);

	vectorAssigment(z, r);						// z0 = r0

	min_discrepancy = std::numeric_limits<double>::infinity();

	norm_f = calcNorm(f);
	scalar_rr_prev = scalarProduct(r, r);		// (rk-1, rk-1)
	do {
		if (iter > max_iter) {
			solver << "max iter exit: discrepance == " << discrepancy << std::endl;
			break;
		}

		multiplicateWithVector(z, temp);			// temp2 = At * A * z
		multiplicateTransWithVector(temp, temp2);

		a = scalar_rr_prev / scalarProduct(temp2, z);	// a = ...

		vectorMultOnConst(z, a, temp);			// xk = xk-1 + ak * zk-1
		vectorAdd(x, temp, x);

		// step out:
		if (calcNorm(temp) / calcNorm(x) < step_eps) {
			//std::cout << "step OUT!\n"; 
			solver << "step out on " << iter << " with  discrepance == " << discrepancy << std::endl;
			break;
		}

		vectorMultOnConst(temp2, a, temp);		// rk = rk-1 - At * A * zk-1
		vectorSubtract(r, temp, r);


		scalar_rr = scalarProduct(r, r);			// b = (rk,rk)/(rk-1, rk-1)
		b = scalar_rr / scalar_rr_prev;
		scalar_rr_prev = scalar_rr;

		vectorMultOnConst(z, b, temp);			// zk = rk + b * zk-1
		vectorAdd(r, temp, z);

		discrepancy = calcNorm(r) / norm_f;		// discrepancy
		iter++;
		//#ifdef CALCULATE
		//		printf("\riter: %d\t discr: %.6e", iter, discrepancy);
		//#endif
		if (discrepancy < min_discrepancy) {
			min_discrepancy = discrepancy;
			vectorAssigment(x_min, x);
		}

#ifdef _DEBUG
		multiplicateWithVector(x, temp);
		vectorSubtract(temp, f, temp);
		discrepancy = calcNorm(temp);


		std::cout << "rr: " << discrepancy << " Discrepancy " << calcNorm(temp) << std::endl;

#endif


	} while (discrepancy > eps);

	solver << "good solve, discrepance == " << discrepancy << std::endl;
	vectorAssigment(x, x_min);
	return iter;
}

void SparseMatrix::
MCG_LU(std::vector<double>& l, std::vector<double>& ud, std::vector<double>& u,
	std::vector<double>& x, std::vector<double>& x_min, std::vector<double>& f,
	std::vector<double>& r, std::vector<double>& z,
	std::vector<double>& temp, std::vector<double>& temp2)
{
	double discrepancy, min_discrepancy, norm_f, scalar_rr, scalar_rr_prev;
	double a, b;
	int iter = 0;

	// partial LU decomposition	
	LUpartialDecomposite(l, ud, u);

	// By = g
	// B = U-t * At * L-t * L-1 * A * U-1
	// y = U * x
	// g = U-t * At * L-t * L-1 * f

	// r = U-t * At * L-t * L-1 * ( f - Ax )
	multiplicateWithVector(x, temp);
	vectorSubtract(f, temp, temp);
	forwardGauss(l, temp2, temp);
	backwardGauss(l, temp, temp2);
	multiplicateTransWithVector(temp, temp2);
	forwardGauss(u, ud, r, temp2);		// r = ...

	vectorAssigment(z, r);				//  z = r

	vectorAssigment(temp, x);
	multiplicateUfactor(u, ud, temp, x);			// x = Ux

	min_discrepancy = std::numeric_limits<double>::infinity();
	norm_f = calcNorm(f);
	scalar_rr_prev = scalarProduct(r, r);

	do {
		if (iter > max_iter) {
			vectorAssigment(x, x_min);
#ifdef CALCULATE
			printf("\riter > max_iter = %d", max_iter);
#endif
			break;
		}

		// ak = ...
		vectorAssigment(temp, z);
		backwardGauss(u, ud, temp2, temp);
		multiplicateWithVector(temp2, temp);
		forwardGauss(l, temp2, temp);
		backwardGauss(l, temp, temp2);
		multiplicateTransWithVector(temp, temp2);
		forwardGauss(u, ud, temp, temp2);

		a = scalar_rr_prev / scalarProduct(temp, z);

		// rk = ...
		vectorMultOnConst(temp, a, temp);
		vectorSubtract(r, temp, r);

		//x'k = x'k-1 - a * z
		vectorMultOnConst(z, a, temp);
		vectorAdd(x, temp, x);

		// step out:
		if (calcNorm(temp) / calcNorm(x) < step_eps) {
#ifdef CALCULATE
			printf("\nstep out\n");
#endif
			break;
		}

		// b = ...
		scalar_rr = scalarProduct(r, r);
		b = scalar_rr / scalar_rr_prev;
		scalar_rr_prev = scalar_rr;

		// z = ...
		vectorMultOnConst(z, b, temp);
		vectorAdd(r, temp, z);

		discrepancy = calcNorm(r) / norm_f;		// discrepancy
		iter++;
#ifdef CALCULATE
		printf("\riter: %d\t discr: %.6e", iter, discrepancy);
#endif
		if (discrepancy < min_discrepancy) {
			min_discrepancy = discrepancy;
			vectorAssigment(x_min, x);
		}
	} while (discrepancy > eps);

	///// x = U-1 * x'
	vectorAssigment(temp, x);
	backwardGauss(u, ud, x, temp);
}



void SparseMatrix::
MCG_DIAG(std::vector<double>& diag, std::vector<double>& x, std::vector<double>& x_min, std::vector<double>& f,
	std::vector<double>& r, std::vector<double>& z,
	std::vector<double>& temp, std::vector<double>& temp2) {

	double discrepancy, min_discrepancy, norm_f, scalar_rr, scalar_rr_prev;
	double a, b;
	int iter = 0;

	// L = U = sqrt(di)
	diagDecomposite(diag);

	// r = U-t * At * L-t * L-1 * ( f - Ax )
	multiplicateWithVector(x, temp);
	vectorSubtract(f, temp, temp);
	diagGauss(diag, temp2, temp);
	diagGauss(diag, temp, temp2);
	multiplicateTransWithVector(temp, temp2);
	diagGauss(diag, r, temp2);		// r = ...

	vectorAssigment(z, r);				//  z = r

	vectorAssigment(temp, x);
	//multiplicate_ufactor(u, ud, temp, x);			// x = Ux
	multiplicateDiagWithVector(temp, x);

	min_discrepancy = std::numeric_limits<double>::infinity();
	norm_f = calcNorm(f);
	scalar_rr_prev = scalarProduct(r, r);

	do {
		if (iter > max_iter) {
			vectorAssigment(x, x_min);
#ifdef CALCULATE
			printf("\niter > max_iter = %d\n", max_iter);
#endif
			break;
		}

		// ak = ...
		vectorAssigment(temp, z);
		diagGauss(diag, temp2, temp);
		multiplicateWithVector(temp2, temp);
		diagGauss(diag, temp2, temp);
		diagGauss(diag, temp, temp2);
		multiplicateTransWithVector(temp, temp2);
		diagGauss(diag, temp, temp2);

		a = scalar_rr_prev / scalarProduct(temp, z);

		// rk = ...
		vectorMultOnConst(temp, a, temp);
		vectorSubtract(r, temp, r);

		//x'k = x'k-1 - a * z
		vectorMultOnConst(z, a, temp);
		vectorAdd(x, temp, x);

		// step out:
		if (calcNorm(temp) / calcNorm(x) < step_eps) {
#ifdef CALCULATE
			printf("\nstep out\n");
#endif
			break;
		}
		// b = ...
		scalar_rr = scalarProduct(r, r);
		b = scalar_rr / scalar_rr_prev;
		scalar_rr_prev = scalar_rr;

		// z = ...
		vectorMultOnConst(z, b, temp);
		vectorAdd(r, temp, z);

		discrepancy = calcNorm(r) / norm_f;		// discrepancy
		iter++;
#ifdef CALCULATE
		printf("\riter: %d\t discr: %.6e", iter, discrepancy);
#endif
		if (discrepancy < min_discrepancy) {
			min_discrepancy = discrepancy;
			vectorAssigment(x_min, x);
		}
	} while (discrepancy > eps);

	// x = U-1 * x'
	vectorAssigment(temp, x);
	diagGauss(diag, x, temp);
}


int SparseMatrix::BCGSTAB(std::vector<double> &l, std::vector<double> & ud, std::vector<double> & u,
	std::vector<double> & x, std::vector<double> & f, std::vector<double> & r, std::vector<double> &r_0,
	std::vector<double> & z, std::vector<double> & temp, std::vector<double> & temp2, std::vector<double> &v, std::vector<double> &p,
	std::vector<double> & y, std::vector<double> &h, std::vector<double> &s, std::vector<double> &t) {

	//eps = 1e-15;

	double discrepancy;
	int iter = 0;

	// partial LU decomposition	
	LUpartialDecomposite(l, ud, u);

	// r0 = Ax - b
	multiplicateWithVector(x, temp);
	vectorSubtract(f, temp, r_0);

	vectorAssigment(r, r_0);
	// z0 = r0
	z = r_0;

	double ro = 1, roi;
	double a = 1;
	double w = 1;
	double b = 0;

	std::vector<double> xPrev = x;



	// iteration process:
	do {
		iter++;

		// pi = 
		roi = scalarProduct(r_0, r);

		// b = 
		b = (roi / ro) * (a / w);
		ro = roi;

		// p = r(i-1) + b(p(i-1) - w(i-1)*v(i-1))
		vectorMultOnConst(v, w, temp);
		vectorSubtract(p, temp, p);
		vectorMultOnConst(p, b, p);
		vectorAdd(p, r, p);

		// y = U^-1  * L^-1 * p <--> LUy = p
		forwardGauss(l, temp, p);
		backwardGauss(u, ud, y, temp);

		// v = Ay
		multiplicateWithVector(y, v);

		// 
		a = roi / scalarProduct(r_0, v);

		// 
		vectorMultOnConst(y, a, temp);
		vectorAdd(x, temp, h);

		/* // check if "h' is accurate enough then xi = h and QUIT
		multiplicate_with_vector(h, temp);
		vector_subtract(temp, f, temp);
		discrepancy = calc_norm(temp);


		if (discrepancy < eps) {
			vector_assigment(x, h);
			break;
		}
		*/

		// 
		vectorMultOnConst(v, a, temp);
		vectorSubtract(r, temp, s);

		//
		forwardGauss(l, temp, s);
		backwardGauss(u, ud, z, temp);

		// 
		multiplicateWithVector(z, t);

		// 
		forwardGauss(l, temp, t);
		forwardGauss(l, temp2, s);

		w = scalarProduct(temp, temp2) / scalarProduct(temp, temp);

		// 
		vectorMultOnConst(z, w, temp);
		vectorAdd(h, temp, x);

		// check step exit:
		vectorSubtract(xPrev, x, temp);
		double stepValue = calcNorm(temp) / calcNorm(x);

		if (stepValue < step_eps) {
			break;
		}


		vectorAssigment(xPrev, x);

		//
		vectorMultOnConst(t, w, temp);
		vectorSubtract(s, temp, r);

		discrepancy = calcNorm(r);
		std::cout << "iter:" << iter << "   dscr: " << discrepancy  << std::endl;

	
		if (discrepancy < eps) {
			break;
		}

	} while (iter < max_iter); // exit by max iter
	
	return iter;
}

#pragma endregion

void SparseMatrix::setFirstBoundaryCondition(int ind) {
	//set diag elem = 1
	di[ind] = 1;

	// zeroing off-diagonal elements before diag
	for (int i = ig[ind]; i < ig[ind + 1]; i++) {
		ggl[i] = 0;
	}

	// zeroing off-diagonal elements after diag
	for (int i = ig[ind+1]; i < ig[n]; i++) {
		if (jg[i] == ind) {
			ggu[i] = 0;
			//ggl[i] = 0;
		}
	}
}
