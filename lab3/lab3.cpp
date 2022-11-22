#include <fstream>
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

    double h = L / (numberOfPointsL - 1);
    double tau = T / (numberOfPointsT - 1);

    SolutionOfEquations solutionForY(L, T, a, beta, numberOfPointsL, numberOfPointsT, [](double x) { return 0.; }, [](double t) { return 100.; });

    vector<vector<double>> u = solutionForY.straightTask();

    //NM::printArray(A);
    //NM::printArray(b);

    fstream trueSolution("trueSolution.txt", ios::out);

    for(int i = 0; i < u.size(); i++) {
        for(int j = 0; j  < u[i].size(); j++) {
            trueSolution << i * tau << " " << j * h << " " << u[i][j] << endl;
        }
    }

    trueSolution.close();



    vector<double> y = u[u.size() - 1];



    fstream ySolution("ySolution.txt", ios::out);

    for (int i = 0; i < y.size(); i++) {
        ySolution << i * h << " " << y[i] << endl;
    }

    ySolution.close();




}
