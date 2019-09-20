#include "pch.h"
#include "task.h"
#include <fstream>
#include "matrix.h"
#include "assert.h"
#include <algorithm>
#include <string>
#include <iostream>

void Task::setParams() {
	lambda = std::vector<double>(amountSubareas);
	sigma = std::vector<func2>(amountSubareas);
	
	lambda[0] = 1;

	uExact = [](const double& r, const double& z) {return r*r; };
	fFunc = [](const double& r, const double& z) {return -4; };
}

void Task::formatingGlobalMatrixPortrait() {
	std::vector<std::vector<int>> temp(nodes.size());

	std::function<bool(std::vector<int>&, double)> hasValue = [] (std::vector<int>& v, double value) {
		if (std::find(v.begin(), v.end(), value) == v.end()) return false;
		else return true;
	};

	// considered links
	for (auto &el : elems) {
		if (!hasValue(temp[el.rb], el.lb)) temp[el.rb].push_back(el.lb);

		if (!hasValue(temp[el.lt], el.lb)) temp[el.lt].push_back(el.lb);
		if (!hasValue(temp[el.lt], el.rb)) temp[el.lt].push_back(el.rb);

		if (!hasValue(temp[el.rt], el.lb)) temp[el.rt].push_back(el.lb);
		if (!hasValue(temp[el.rt], el.rb)) temp[el.rt].push_back(el.rb);
		if (!hasValue(temp[el.rt], el.lt)) temp[el.rt].push_back(el.lt);
	}

	// sort columns
	for (auto& columns : temp) {
		std::sort(columns.begin(), columns.end());
	}

	// compute jj vector:
	std::vector<int> jg;
	for (auto &vec : temp) {
		std::copy(vec.begin(), vec.end(), back_inserter(jg));
	}

	// compute ig vector
	std::vector<int> ig;
	ig.push_back(0);
	for (auto &vec : temp) {
		ig.push_back(ig.back() + vec.size());
	}

	// init global matrix
	globalMatrix.init(nodes.size(), ig, jg);
}

void Task::init() {
#pragma region outputResultFile
	fout.open(R"(C:\Users\krasheninnik.2016\source\repos\2D_mfe_stationary\output/result.txt)");
	fout.precision(17);
#pragma endregion

#pragma region InitSpaceGrid 
{
	nodes = std::vector<Node>();
	elems = std::vector<FiniteElem>();
	subareas = std::vector<int>();

	std::ifstream fin(R"(C:\Users\krasheninnik.2016\source\repos\2D_mfe_stationary\input\grid_space.txt)");

	int numSpaceGridDividion = 1; // auto grid dividion 
	int numOfAreasX = 1; // amount areas with equalent coefficien of unevenness for x
	int numOfAreasY = 1; // amount areas with equlent coefficien of unevenness for y
	double xStart, yStart, numOfElems, step, coef;

	// auto grid divideion coefficient
	fin >> numSpaceGridDividion;
	int k = pow(2, numSpaceGridDividion - 1);
		
	// read axis
	fin >> numOfAreasX >> numOfAreasY;

	// read X axis
	fin >> xStart >> numOfElems >> step >> coef;

	// calculate grid parameters for unevenness:
	numOfElems *= k;
	coef = pow(coef, 1.0 / k);

	// calculate first x-step
	double stepsCoef = 0;
	for (int i = 0; i < k; i++) stepsCoef += pow(coef, i);
	step /= stepsCoef;

	// remember x-axis:
	std::vector<double> xaxis = std::vector<double>();

	double x = xStart;	
	xaxis.push_back(x);		// add x0 in nodes

	for (int j = 0; j < numOfAreasX; j++) {
		for (int i = 0; i < numOfElems; i++) {
			x += step;
			xaxis.push_back(x);
			step *= coef;					// change step
		}

		if (j != numOfAreasX - 1) {
			// read info about new area
			fin >> numOfElems >> coef;
			// consider grid dividion coef
			numOfElems *= k;
			coef = pow(coef, 1.0 / k);
		}
	}


	// input y axis params
	fin >> yStart >> numOfElems >> step >> coef;

	// calculate grid parameters for unevenness:
	numOfElems *= k;
	coef = pow(coef, 1.0 / k);

	// calculate first y-step
	/*double*/ stepsCoef = 0;
	for (int i = 0; i < k; i++) stepsCoef += pow(coef, i);
	step /= stepsCoef;

	// remember y-axis:
	std::vector<double> yaxis = std::vector<double>();

	double y = yStart;
	yaxis.push_back(y);		// add y0 in nodes

	for (int j = 0; j < numOfAreasY; j++) {
		for (int i = 0; i < numOfElems; i++) {
			y += step;
			yaxis.push_back(y);
			step *= coef;					// change step
		}

		if (j != numOfAreasX - 1) {
			// read info about new area
			fin >> numOfElems >> coef;
			// consider grid dividion coef
			numOfElems *= k;
			coef = pow(coef, 1.0 / k);
		}
	}

	// forming grid:
	for (int i = 0; i < yaxis.size(); i++) {
		for (int j = 0; j < xaxis.size(); j++) {
			nodes.push_back({ xaxis[j], yaxis[i] });
		}
	}

	// forming array of elems (relation with local and global numbers)
	elems.reserve((yaxis.size() - 1) * (xaxis.size() - 1));

	const int yshift = xaxis.size();
	int shift = 0;

	for (int i = 0; i < yaxis.size() - 1; i++) {

		for (int j = 0; j < xaxis.size() - 1; j++) {
			elems.push_back({ j + shift, j + shift + 1, j + shift + yshift, j + shift + yshift + 1});
		}
		
		shift += yshift;
	}

	raxisSize = xaxis.size();
	zaxisSize = yaxis.size();

	// calculating elems with first boundary conditions

	for (int i = 1; i < yaxis.size(); i++) {
		firstBoundaryConditionsNodesIndex.push_back(i * xaxis.size() - 1);
	}

	for (int i = xaxis.size() * (yaxis.size() - 1); i < nodes.size(); i++) {
		firstBoundaryConditionsNodesIndex.push_back(i);
	}
	// calculating elems with first boundary conditions
	///////////////////////////////////////////////////////////
	/////////////////////////FOR TEST//////////////////////////
	///////////////////////////////////////////////////////////
	for (int i = 0; i < yaxis.size(); i++) {
		firstBoundaryConditionsNodesIndex.push_back(i * xaxis.size());
	}

	for (int i = 0; i < xaxis.size(); i++) {
		firstBoundaryConditionsNodesIndex.push_back(i);
	}

	///////////////////////////////// this, that down need be tested ///////////////////////////////
	// fill subareas: (type of parameters in area)
	int numOfFiniteElems = 0;
	int sum = 0;

	fin >> amountSubareas;

	subareas = std::vector<int>(elems.size()); // by default type = 0

	for (int type = 1; type < amountSubareas; type++) {
		int bottom = 0, top = 0, left = 0, right = 0;
		fin >> bottom, top, left, right;

		// consider coeff of auto grid dividion:
		bottom*= k; top*= k; left*= k; right *= k;

		// change type of subarea
		int shift = yaxis.size() * bottom;

		for (int ytemp = bottom; ytemp <= left; ytemp++) {
			for (int xtemp = left; xtemp <= right; xtemp++) {
				subareas[xtemp + shift] = type;
			}
			shift += yaxis.size();
		}
	}	
	
	// init vector of params of equals in subareas:
	lambda = std::vector<double>(amountSubareas);
	sigma = std::vector<func2>(amountSubareas); ////////////// think about this !!!

	fin.close();
}
#pragma endregion

#pragma region MatrixInit
/////////////////////////////////////////////////////////////////////////////////////////////////////////////
	const int matrixDim = 3 + (elems.size() - 1) * 2;
	formatingGlobalMatrixPortrait();
#pragma endregion

#pragma region MemoryAllocation

	q = std::vector<double>(nodes.size());
	//qPrev = std::vector<double>(nodes.size());
	//qPrevTime = std::vector<double>(nodes.size());
	qExact = std::vector<double>(nodes.size());;


	// there should be vector<vvector<doub>> results -> for save q on each step!
	/*temp = std::vector<double>(nodes.size());*/

	// local matrix of mass and rigidity, and local vector;
	const int pointsToFiniteElem = 4;
	massLocalMatrix = LocalMatrix(pointsToFiniteElem);
	for (auto& el : massLocalMatrix) el = std::vector<double>(pointsToFiniteElem);
	
	rigidLocalMatrix = LocalMatrix(pointsToFiniteElem);
	for (auto& el : rigidLocalMatrix) el = std::vector<double>(pointsToFiniteElem);

	fLocal = std::vector<double>(pointsToFiniteElem);

	// global matrix!
	f = std::vector<double>(nodes.size());;
	
	//TapeMatrix globalMatirx;
#pragma endregion
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void Task::calculateGlobalMatrixAndRightPart() {

	// calculate global matrix and vector:
	for (int elemNum = 0; elemNum < elems.size(); elemNum++) {

		calculateLocalMatrixOfRigid(elemNum);
		calculateLocalRightPart(elemNum);

		
		addLocalMatrixToGlobal(elemNum);
		addLocalRigtPartToGlobal(elemNum);
	}

	globalMatrix.fillGGU();
	// set boundary conditions 
	// secondary boundary conditions dont influence, === 0!
	setFirstBoundaryConditions();
	//globalMatrix.fillGGU();
}



void Task::solve() {
	int matrixDim = nodes.size();
	int amountElems = globalMatrix.getAmountElems();
	std::vector<double> r = std::vector<double>(matrixDim);
	std::vector<double> r_0 = std::vector<double>(matrixDim);

	std::vector<double> z = std::vector<double>(matrixDim);
	std::vector<double> p = std::vector<double>(matrixDim);
	std::vector<double> temp = std::vector<double>(matrixDim);
	std::vector<double> temp2 = std::vector<double>(matrixDim);

	std::vector<double> v = std::vector<double>(matrixDim);
	std::vector<double> y = std::vector<double>(matrixDim);
	std::vector<double> h = std::vector<double>(matrixDim);
	std::vector<double> s = std::vector<double>(matrixDim);
	std::vector<double> t = std::vector<double>(matrixDim);

	std::vector<double> ud = std::vector<double>(matrixDim);
	std::vector<double> l = std::vector<double>(amountElems);
	std::vector<double> u = std::vector<double>(amountElems);

	// need input things:
	epsDiscrep = 1e-15;
	calculateGlobalMatrixAndRightPart();
	//int iter = globalMatrix.BCGSTAB(l, ud, u, q, f, r, r_0, z, temp, temp2, v, p, y, h, s, t);
	globalMatrix.LOS(q, f, r, z, p, temp);
	//globalMatrix.MCG(q, p, f, r, z, temp, temp2);

	//LOS_LU(std::vector<double>& l, std::vector<double>& ud, std::vector<double>& u,
	//	std::vector<double>& x, std::vector<double>& f,
	//	std::vector<double>& r, std::vector<double>& z,
	//	std::vector<double>& p, std::vector<double>& temp,
	//	std::vector<double>& temp2) {

	//globalMatrix.LOS_LU(l, ud, u, q, f, r, z, p, temp, temp2);


}

void Task::calculateLocalMatrixOfRigid(uint32_t elemNum) {	

	const auto& elem = elems[elemNum];
	const double coef = 1;

	const double R1 = nodes[elem.lb].r;
	const double R2 = nodes[elem.rb].r;
	const double Z1 = nodes[elem.lb].z;
	const double Z2 = nodes[elem.lt].z;
	const double hr = R2 - R1;
	const double hz = Z2 - Z1;

	/////////////////////////////////////////////////////////////
	// integrals:

	const double dR1dR1r = -1.0 * (R1 + R2) / (2.0 * (R1 - R2));
	const double dR1dR2r = (R1 + R2) / (2.0 * (R1 - R2));
	const double dR2dR2r = -1.0  * (R1 + R2) / (2.0 * (R1 - R2));
	const double R1R1r = -1.0 * (R1 - R2) * (3.0 * R1 + R2) / 12.0;
	const double R1R2r = (R2 * R2 - R1 * R1) / 12.0;
	const double R2R2r = -1.0 * (R1 - R2) * (R1 + 3.0 * R2) / 12.0;

	const double dZ1dZ1 = -1.0 / (Z1 - Z2);
	const double dZ1dZ2 = 1.0 / (Z1 - Z2);
	const double dZ2dZ2 = -1.0 / (Z1 - Z2);
	const double Z1Z1 = (Z2 - Z1 ) / 3.0;
	const double Z1Z2 = (Z2 - Z1) / 6.0;
	const double Z2Z2 = (Z2 - Z1) / 3.0;

	////////////////////////////////////////////////////
	// fill local matrix of rigid

	rigidLocalMatrix[0][0] = coef *  (dR1dR1r * Z1Z1 + R1R1r * dZ1dZ1);
	rigidLocalMatrix[0][1] = coef * (dR1dR2r * Z1Z1 + R1R2r * dZ1dZ1);
	rigidLocalMatrix[0][2] = coef * (dR1dR1r * Z1Z2 + R1R1r * dZ1dZ2);
	rigidLocalMatrix[0][3] = coef * (dR1dR2r * Z1Z2 + R1R2r * dZ1dZ2);

	rigidLocalMatrix[1][1] = coef * (dR2dR2r * Z1Z1 + R2R2r * dZ1dZ1);
	rigidLocalMatrix[1][2] = coef * (dR1dR2r * Z1Z2 + R1R2r * dZ1dZ2);
	rigidLocalMatrix[1][3] = coef * (dR2dR2r * Z1Z2 + R2R2r * dZ1dZ2);

	rigidLocalMatrix[2][2] = coef * (dR1dR1r * Z2Z2 + R1R1r * dZ2dZ2);
	rigidLocalMatrix[2][3] = coef * (dR1dR2r * Z2Z2 + R1R2r * dZ2dZ2);

	rigidLocalMatrix[3][3] = coef * (dR2dR2r * Z2Z2 + R2R2r * dZ2dZ2);

	///////////////////////////////////////////////////
}


void multiplicateFullSquareMatrixOnVector(std::vector<std::vector<double>> &M, std::vector<double> &a, std::vector<double> &res) {
	for (int i = 0; i < res.size(); i++) {
		res[i] = 0;

		for (int j = 0; j < a.size(); j++) {
			res[i] += M[i][j] * a[j];
		}
	}
}

void Task::calculateLocalRightPart(uint32_t num) {
	const auto& elem = elems[num];

	// calculate local matrix of mass:

	const double R1 = nodes[elem.lb].r;
	const double R2 = nodes[elem.rb].r;
	const double Z1 = nodes[elem.lb].z;
	const double Z2 = nodes[elem.lt].z;
	const double hr = R2 - R1;
	const double hz = Z2 - Z1;

	/////////////////////////////////////////////////////////////
	// integrals:

	const double R1R1r = -1.0 * (R1 - R2) * (3.0 * R1 + R2) / 12.0;
	const double R1R2r = (R2 * R2 - R1 * R1) / 12.0;
	const double R2R2r = -1.0 * (R1 - R2) * (R1 + 3.0 * R2) / 12.0;
		
	const double Z1Z1 = (Z2 - Z1) / 3.0;
	const double Z1Z2 = (Z2 - Z1) / 6.0;
	const double Z2Z2 = (Z2 - Z1) / 3.0;

	// fill local matrix of mass

	massLocalMatrix[0][0] = R1R1r * Z1Z1;
	massLocalMatrix[0][1] = R1R2r * Z1Z1;
	massLocalMatrix[0][2] = R1R1r * Z1Z2;
	massLocalMatrix[0][3] = R1R2r * Z1Z2;

	massLocalMatrix[1][1] = R2R2r * Z1Z1;
	massLocalMatrix[1][2] = R1R2r * Z1Z2;
	massLocalMatrix[1][3] = R2R2r * Z1Z2;

	massLocalMatrix[2][2] = R1R1r * Z2Z2;
	massLocalMatrix[2][3] = R1R2r * Z2Z2;

	massLocalMatrix[3][3] = R2R2r * Z2Z2;

	// consider symmetrical
	massLocalMatrix[1][0] = massLocalMatrix[0][1];

	massLocalMatrix[2][0] = massLocalMatrix[0][2];
	massLocalMatrix[2][1] = massLocalMatrix[1][2];

	massLocalMatrix[3][0] = massLocalMatrix[0][3];
	massLocalMatrix[3][1] = massLocalMatrix[1][3];
	massLocalMatrix[3][2] = massLocalMatrix[2][3];

	// vector f:
	static std::vector<double> fValues(4);
	fValues[0] = fFunc(R1, Z1);
	fValues[1] = fFunc(R2, Z1);
	fValues[2] = fFunc(R1, Z2);
	fValues[3] = fFunc(R2, Z2);

	// calculate local vector of right part
	multiplicateFullSquareMatrixOnVector(massLocalMatrix, fValues, fLocal);
}

void Task::addLocalMatrixToGlobal(uint32_t num) {
	const auto& elem = elems[num];

	globalMatrix.addDiagElem(elem.lb, rigidLocalMatrix[0][0]);
	globalMatrix.addDiagElem(elem.rb, rigidLocalMatrix[1][1]);
	globalMatrix.addDiagElem(elem.lt, rigidLocalMatrix[2][2]);
	globalMatrix.addDiagElem(elem.rt, rigidLocalMatrix[3][3]);

	globalMatrix.addElem(elem.rb, elem.lb, rigidLocalMatrix[0][1]);
	globalMatrix.addElem(elem.lt, elem.lb, rigidLocalMatrix[0][2]);
	globalMatrix.addElem(elem.lt, elem.rb, rigidLocalMatrix[1][2]);
	globalMatrix.addElem(elem.rt, elem.lb, rigidLocalMatrix[0][3]);
	globalMatrix.addElem(elem.rt, elem.rb, rigidLocalMatrix[1][3]);
	globalMatrix.addElem(elem.rt, elem.lt, rigidLocalMatrix[2][3]);
}

void Task::addLocalRigtPartToGlobal(uint32_t num) {
	const auto& elem = elems[num];

	f[elem.lb] += fLocal[0];
	f[elem.rb] += fLocal[1];
	f[elem.lt] += fLocal[2];
	f[elem.rt] += fLocal[3];


}

void Task::setFirstBoundaryConditions() {
	for(auto &elIndex: firstBoundaryConditionsNodesIndex) {
		// set first bounday conditions in left side:
		globalMatrix.setFirstBoundaryCondition(elIndex);
		f[elIndex] = uExact(nodes[elIndex].r, nodes[elIndex].z);
	}

	//////////////////
	//std::vector<int> leftSide{ 0, 5, 10, 1, 2, 3, 4};
	//for (auto &elIndex : leftSide) {
	//	// set first bounday conditions in left side:
	//	globalMatrix.setFirstBoundaryCondition(elIndex);
	//	f[elIndex] = uExact(nodes[elIndex].r, nodes[elIndex].z);
	//}

}

void vectorSubtraction(std::vector<double>& result, const std::vector<double>& a){
	for (int i = 0; i < result.size(); i++) result[i] -= a[i];
}

double calcNorm(const std::vector<double> &x) {
	double norm = 0;
	for (int i = 0; i < x.size(); i++) {
		norm += x[i] * x[i];
	}
	norm = sqrt(norm);
	return norm;
}


void Task::saveResult() {
	double sum = 0;

	for (int i = 0; i < qExact.size(); i++) {
		qExact[i] = uExact(nodes[i].r, nodes[i].z);
		sum += (qExact[i] - q[i])*(qExact[i] - q[i]);
	}
	sum = sqrt(sum);

	fout << "\tNormOfError: " << sum << std::endl;

	for (int i = 0; i < qExact.size(); i++) {
		fout << qExact[i] << " "; 
		if ((i + 1) % raxisSize == 0)  fout << std::endl;
	}
	fout << std::endl;

	for (int i = 0; i < q.size(); i++) {
		fout << q[i] << " ";
		if ((i + 1) % (raxisSize) == 0)  fout << std::endl;
	}
	fout << std::endl;
	fout << std::endl;
}


/*
bool Task::SimpleIterationDiscrepOut() {
	// || A(qi) * qi - b(qi) || / || b(qi) || < eps => out:

	temp = globalMatrix.multiplicate_with_vector(q, temp);
	vectorSubtraction(temp, f);

	double resultNorm = calcNorm(temp) / calcNorm(f);
	//std::cout << resultNorm << std::endl;

	return resultNorm < epsDiscrep;
}

*/