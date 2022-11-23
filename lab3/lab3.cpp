﻿#include <fstream>
#include "NM.h"
#include "SolutionOfEquations.h"


double phi(double x){
    return 0;
}

double trueManagement(double t){
    return 100;
}


int main(){

    double L = 1, T = 2;
    double a = 2;
    double beta = 4;
    int numberOfPointsL = 100, numberOfPointsT = 100;

    //double h = L / (numberOfPointsL - 1);
    //double tau = T / (numberOfPointsT - 1);

    SolutionOfEquations SolutionForY(L, T, a, beta, numberOfPointsL, numberOfPointsT, [](double x) { return 100. * x; }, [](double t) { return 1000. * t * t *(1 - t); });

    vector<vector<double>> TrueU = SolutionForY.straightTask();

    //NM::printArray(A);
    //NM::printArray(b);





    SolutionOfEquations CurrentSolution(L, T, a, beta, numberOfPointsL, numberOfPointsT, [](double x) { return 100. * x; }, [](double t) { return 50 * t * t; });

    vector<vector<double>> u = CurrentSolution.straightTask();

    CurrentSolution.printFile(u, "USolution.txt");


    vector<vector<double>> psi = CurrentSolution.conjugateTask(u, TrueU[TrueU.size() - 1]);

    CurrentSolution.printFile(psi, "PsiSolution.txt");




}
