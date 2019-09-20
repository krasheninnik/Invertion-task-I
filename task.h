#pragma once
#include <vector>
#include "matrix.h"
#include <fstream>	

struct Node {
	double r;
	double z;
};

struct FiniteElem {
	int lb; // left bottom
	int rb; // right bottom
	int lt; // left top
	int rt; // right top
};


class Task {
	using LocalMatrix = std::vector<std::vector<double>>;

	using func = std::function<double(const double&)>;
	using func2 = std::function<double(const double&, const double&)>;
	using func3d1fun = std::function<double(const double&, const double&, const double&, const func2&)>;

public:
	void init();
	void setParams();
	void solve();
	void saveResult();

private:
	void formatingGlobalMatrixPortrait();

	void calculateGlobalMatrixAndRightPart();
	void setFirstBoundaryConditions();

	void calculateLocalMatrixOfRigid(uint32_t elemNum);
	void calculateLocalRightPart(uint32_t elemNum);

	void addLocalMatrixToGlobal(uint32_t elem);
	void addLocalRigtPartToGlobal(uint32_t elemNum);

	// Simple Iterations Method:
	bool SimpleIterationDiscrepOut();

private:
	// result output stream
	std::ofstream fout;
	
	// time grid:
	std::vector<double> times;			// delete it

	// space grid:
	std::vector<int> subareas;	// type of subareas of elems in which it located.
	std::vector<FiniteElem> elems;
	std::vector<Node> nodes;	// vector of Nodes
	int raxisSize;
	int zaxisSize;

	// solutions:
	std::vector<double> q;
	std::vector<double> qExact;

	// local matrix of mass and rigidity, and local vector;
	LocalMatrix massLocalMatrix;
	LocalMatrix rigidLocalMatrix;
	std::vector<double> fLocal;

	// global matrix
	std::vector<double> f;
	SparseMatrix globalMatrix;

	// parameters of equation in different subareas:
	int amountSubareas;
	func2 uExact;
	func2 fFunc;
	func2 fStart;

	std::vector<double> lambda;
	std::vector<func2> sigma;

	std::vector<int> firstBoundaryConditionsNodesIndex;

	double epsDiscrep;
};