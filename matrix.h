#pragma once
#include <vector>
#include <functional>

//class Matrix {
//	using LocalMatrix = std::vector<std::vector<double>>;
//
//public:
//	void init(int size);
//	void reset();
//	void addLocalMatrix(const int elemNum, const LocalMatrix& M, const LocalMatrix& R, const double dt);
//	
//	void setFirstBoundaryConditionsLeft();
//	void setFirstBoundaryConditionsRight();
//
//	int get_dim();
//
//	double calc_sum(int, std::vector<double>&);// sum of multiplicate elems of row matrix with corresponding vector's elems
//	double calc_relative_discrepancy(std::vector<double> &, std::vector<double>&, std::vector<double>&);
//	std::vector<double> multiplicate_with_vector(std::vector<double>&, std::vector<double>&);
//
//	void LUdecompose();
//

//void solveSystem(const std::vector<double>& f, std::vector<double>& q, std::vector<double>& temp);
//
//private: 
//	static const int diags = 5;
//	
//	int n, block_size, max_iter;
//	double accuracy;
//};
#include <vector>

#define CALCULATE

class SparseMatrix {
	using LocalMatrix = std::vector<std::vector<double>>;
public:
	SparseMatrix();
	void init(const int elemNum, const int max_iter, const double eps);
	void init(const int dim, std::vector<int>& _ig, std::vector<int>& _jg);

	void load();
	void copy(SparseMatrix source);
	void addLocalMatrix(int elemNum, LocalMatrix M);
	int  getDim();
	int  getAmountElems();
	int  getMethod();

	void addElem(int row, int column, double value);
	void addDiagElem(int n, double value);
	void fillGGU();
	void setFirstBoundaryCondition(int n);
	void resetMatrix();

	void directSolveLU(std::vector<double>& l, std::vector<double>& ud,
		std::vector<double>& u, std::vector<double>& temp, std::vector<double>& q, std::vector<double>& f);

	void LUpartialDecomposite(std::vector<double>& l,
		std::vector<double>& ud, std::vector<double>& u);

	void diagDecomposite(std::vector<double>& diag);

	void forwardGauss(std::vector<double>& matrix, std::vector<double>& x,
		std::vector<double>& f);
	void forwardGauss(std::vector<double>& matrix, std::vector<double>& diag,
		std::vector<double>& x, std::vector<double>& f);

	void backwardGauss(std::vector<double>& matrix, std::vector<double>& x,
		std::vector<double>& f);
	void backwardGauss(std::vector<double>& matrix, std::vector<double>& diag,
		std::vector<double>& x, std::vector<double>& f);

	void diagGauss(std::vector<double>& diag, std::vector<double>& x, std::vector<double>& f);

	void multiplicateUfactor(std::vector<double>& u, std::vector<double>& diag,
		std::vector<double>& b, std::vector<double>& result);

	void LOS(std::vector<double>& x0, std::vector<double>& f,
		std::vector<double>& r, std::vector<double>& z,
		std::vector<double>& p, std::vector<double>& temp);

	void LOS_LU(std::vector<double>& l, std::vector<double>& ud, std::vector<double>& u,
		std::vector<double>& x0, std::vector<double>& f,
		std::vector<double>& r, std::vector<double>& z,
		std::vector<double>& p, std::vector<double>& temp,
		std::vector<double>& temp2);

	void LOS_DIAG(std::vector<double>& diag, std::vector<double>& x, std::vector<double>& f,
		std::vector<double>& r, std::vector<double>& z,
		std::vector<double>& p, std::vector<double>& temp,
		std::vector<double>& temp2);

	int MCG(std::vector<double>& x0, std::vector<double>& x_min,
		std::vector<double>& f, std::vector<double>& r, std::vector<double>& z,
		std::vector<double>& temp, std::vector<double> temp2);

	void MCG_LU(std::vector<double>& l, std::vector<double>& ud, std::vector<double>& u,
		std::vector<double>& x, std::vector<double>& x_min, std::vector<double>& f,
		std::vector<double>& r, std::vector<double>& z,
		std::vector<double>& temp, std::vector<double>& temp2);

	void MCG_DIAG(std::vector<double>& diag, std::vector<double>& x, std::vector<double>& x_min,
		std::vector<double>& f, std::vector<double>& r, std::vector<double>& z,
		std::vector<double>& temp, std::vector<double>& temp2);


	int BCGSTAB(std::vector<double> &l, std::vector<double> & ud, std::vector<double> & u,
		std::vector<double> & x, std::vector<double> & f, std::vector<double> & r, std::vector<double> &r_0,
		std::vector<double> & z, std::vector<double> & temp, std::vector<double> & temp2, std::vector<double> &v, std::vector<double> &p,
		std::vector<double> & y, std::vector<double> &h, std::vector<double> &s, std::vector<double> &t);

	void multiplicateWithVector(std::vector<double>& b, std::vector<double>& result);
	void multiplicateTransWithVector(std::vector<double>& b, std::vector<double>& result);
	void multiplicateDiagWithVector(std::vector<double>& b, std::vector<double>& result);

private:
	std::vector<int> ig;		// index elems - start of row
	std::vector<int> jg;		// column index of elem
	std::vector<double> di;		// diag elems
	std::vector<double> ggu;	// elem's of upper trinagle
	std::vector<double> ggl;	// elem's of lower trinagle

	int n;						// dimension
	int max_iter;				// max amount of iter
	double eps;					// min relative discrepancy
	const double step_eps = 1e-25;
	int method;
};