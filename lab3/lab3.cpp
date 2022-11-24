#include <fstream>
#include "NM.h"
#include "SolutionOfEquations.h"

int main(){

    double L = 1, T = 1;
    double a = 2;
    double beta = 4;
    int numberOfPointsL = 100, numberOfPointsT = 100;

    function<double(double)> phi = [](double x) { return 100. * x; };

    SolutionOfEquations SolutionForY(L, T, a, beta, numberOfPointsL, numberOfPointsT, phi, [](double t) { return 1000. * t * t * (1 - t); });

    vector<vector<double>> trueU = SolutionForY.straightTask();

    //NM::printArray(A);
    //NM::printArray(b);





    SolutionOfEquations CurrentSolution(L, T, a, beta, numberOfPointsL, numberOfPointsT, phi, [](double t) { return 1000 * sin(10 * t); });


    vector<double> approxY = CurrentSolution.calculation(trueU[trueU.size() - 1]);

    CurrentSolution.printFile(trueU[trueU.size() - 1], "trueY.txt");
    CurrentSolution.printFile(approxY, "approxY.txt");




}
