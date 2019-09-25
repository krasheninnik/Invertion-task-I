#define _USE_MATH_DEFINES
#include "inversion.h"
#include <fstream>
#include "math.h"
#include <iostream>

void InversionTask::init() {
	// set initianl params
	eps = 1e-12;
	maxIter = 1000;

	Iexact = 1;
	Istart = 0.1;
	sigma = 0.1;
	receivers = std::vector<XCoords>(3);

	// input params of sources and receivers
	std::ifstream fin(R"(C:\Users\krasheninnik.2016\source\repos\2D_mfe_stationary\input\inversion.txt)");
	fin >> Iexact >> Istart;

	fin >> source.x0 >> source.x1;
	fin >> receivers[0].x0 >> receivers[0].x1;
	fin >> receivers[1].x0 >> receivers[1].x1;
	fin >> receivers[2].x0 >> receivers[2].x1;

	// open file to output results
	std::ofstream fout("inversionResult.txt");
}


double f(double I, double sigma) {
	return 1;
}

// this function can calculate derivates? when I == 1, v = derivates vector!
void InversionTask::calculateValues(double I, std::vector<double>& v) {
	for (int i = 0; i < v.size(); i++) {
		v[i] = I / (2 * M_PI * sigma) *
			((1.0 / abs(source.x1 - receivers[i].x0) - 1.0 / abs(source.x0 - receivers[i].x0)) -
			(1.0 / abs(source.x1 - receivers[i].x1) - 1.0 / abs(source.x0 - receivers[i].x1)));
	}
}

double InversionTask::calculateFunctional(std::vector<double>& exact, std::vector<double>& curr) {
	double fine = 0;

	for (int i = 0; i < exact.size(); i++) {		
		fine += 1.0 / (exact[i] * exact[i]) * (exact[i] - curr[i]) *  (exact[i] - curr[i]);
	}

	return fine;
}

double InversionTask::formingMatrix(std::vector<double>& exact, std::vector<double>& der) {
	double res = 0;
	for (int i = 0; i < der.size(); i++) {
		res += (der[i] * der[i])  / (exact[i] * exact[i]) ;
	}
	return res;
}

double InversionTask::formingVector(std::vector<double>& exact, std::vector<double>& curr,
		std::vector<double>& der) {
	double res = 0;
	for (int i = 0; i < der.size(); i++) {
		res += -1.0 * ((curr[i]-exact[i]) * der[i]) / (exact[i] * exact[i]);
	}
	return res;
}

void InversionTask::analiticalSolve(std::vector<double>& currentValues, std::vector<double>& exactValues,
	std::vector<double>& derivates) {
	int iter = 0;
	double dI = 0;
	double Ilast = 0;
	double functional = 0;

	calculateValues(Iexact, exactValues);   // calculate exact values

	while (iter < maxIter) {
		calculateValues(Icurr, currentValues); // calculate current values
		functional = calculateFunctional(exactValues, currentValues);

		fout << "iter: " << iter++ << " I: " << Icurr << " functional: " << functional << std::endl;
		if (abs(functional) < eps) break;

		// calculate derivates
		calculateValues(1.0, derivates); // calculate exact values

		// caluclate system of equation (there n = 1 )
		double A = formingMatrix(exactValues, derivates);
		double b = formingVector(exactValues, currentValues, derivates);
		dI = b / A;
		Ilast = Icurr;
		Icurr += dI;
	}
}

void InversionTask::initTask() {
	Task T;
	T.init();
}

// this function can calculate derivates? when I == 1, v = derivates vector!
void InversionTask::calculateValuesNum(double I, std::vector<double>& v) {
	const int amountNumValues = 8;
	static std::vector<double> numValues(amountNumValues);

	T.setParams(I);
	T.getValues(numValues);

	const double A = numValues[0];
	const double B = numValues[1];
	const double M1 = numValues[2];
	const double N1 = numValues[3];
	const double M2 = numValues[4];
	const double N2 = numValues[5];
	const double M3 = numValues[6];
	const double N3 = numValues[7];

	/*v[0] = 
	v[1] =
	v[2] = */







	//for (int i = 0; i < v.size(); i++) {
	//	v[i] = I / (2 * M_PI * sigma) *
	//		((1.0 / abs(source.x1 - receivers[i].x0) - 1.0 / abs(source.x0 - receivers[i].x0)) -
	//		(1.0 / abs(source.x1 - receivers[i].x1) - 1.0 / abs(source.x0 - receivers[i].x1)));
	//}
}

void InversionTask::numericalSolve(std::vector<double>& currentValues, std::vector<double>& exactValues,
	std::vector<double>& derivates) {

	int iter = 0;
	double dI = 0;
	double Ilast = 0;
	double functional = 0;

	calculateValuesNum(Iexact, exactValues);   // calculate exact values

	while (iter < maxIter) {
		calculateValuesNum(Icurr, currentValues); // calculate current values
		functional = calculateFunctional(exactValues, currentValues);

		fout << "iter: " << iter++ << " I: " << Icurr << " functional: " << functional << std::endl;
		if (abs(functional) < eps) break;

		// calculate derivates
		calculateValuesNum(1.0, derivates); // calculate exact values

		// caluclate system of equation (there n = 1 )
		double A = formingMatrix(exactValues, derivates);
		double b = formingVector(exactValues, currentValues, derivates);
		dI = b / A;
		Ilast = Icurr;
		Icurr += dI;
	}

}

void InversionTask::solve() {
	const int amountReceivers = receivers.size();
	const int amountParameters = 1;

	bool analitical = true;
	double fine;

	std::vector<double> exactValues(amountReceivers);
	std::vector<double> currentValues(amountReceivers);
	std::vector<double> derivates(amountReceivers);

	double Ilast = Istart;
	double functional;
	
	double dI;

	double iter = 0;
	if (analitical) analiticalSolve(currentValues, exactValues, derivates);
	else {
		initTask();
		numericalSolve(currentValues, exactValues, derivates);
	}

	fout.close();
}



