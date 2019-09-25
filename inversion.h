#pragma once
#include "task.h"
#include <vector>

struct XCoords {
	double x0;
	double x1;
};

class InversionTask {

public:
	void init();
	void solve();

private:
	void calculateValues(double p, std::vector<double>& v);
	void calculateValuesNum(double I, std::vector<double>& v);

	double calculateFunctional(std::vector<double>& exact, std::vector<double>& curr);
	double formingMatrix(std::vector<double>& exact, std::vector<double>& der);
	double formingVector(std::vector<double>& exact, std::vector<double>& curr, std::vector<double>& der);

	void analiticalSolve(std::vector<double>& currentValues, std::vector<double>& exactValues,
		std::vector<double>& derivates);
	void numericalSolve(std::vector<double>& currentValues, std::vector<double>& exactValues,
		std::vector<double>& derivates);
	

	void initTask();

	std::ofstream fout; // file for output result

	int maxIter;
	double eps;

	double Iexact;
	double Istart;
	double Icurr;	// current
	double sigma;

	Task T;
	std::vector<XCoords> receivers;
	XCoords source;
};