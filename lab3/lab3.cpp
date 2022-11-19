#include "NM.h"


double phi(double x){
    return x;
}

double trueManagement(double t){
    return 1000 * t * t * (1 - t);
}

double derivative_U_X_end(double t, double lastUInL){
    double beta = 4;
    return beta * (trueManagement(t) - lastUInL);
}

int main(){

    double L = 1, T = 1;
    double a = 2;
    int numberOfPointsL = 50, numberOfPointsT = 2;






    double  a2 = a * a;
    double h = L / (numberOfPointsL - 1);
    double tau = L / (numberOfPointsT - 1);
    vector<vector<double>> u(numberOfPointsT, vector<double>(numberOfPointsL));

    for(int i = 0; i < u[0].size(); i++){
        u[0][i] = phi(h * i);
    }

    vector<vector<double>> A(numberOfPointsL - 2, vector<double>(numberOfPointsL - 2));
    vector<double> b(numberOfPointsL - 2);

    double ku_1 = a2 / (h * h);
    double ku_2 = -((2. * a2) / (h * h) + 1. / tau);
    double ku_3 = ku_1;

    double ku_4 = -1. / tau;


    for(int i = 1; i < numberOfPointsT; i++){

        double p1 = 0;

        A[0][0] = ku_2 + ku_1;
        A[0][1] = ku_3;

        b[0] = ku_4 + ku_1 * p1 * h;


            for(int j = 1; j < A.size() - 1; j++){
        
                A[j][j - 1] = ku_1;
                A[j][j] = ku_2;
                A[j][j + 1] = ku_3;

                b[j] = ku_4 *  u[i - 1][j];


            }


        double p2 = derivative_U_X_end(i * tau, u[i - 1][u[i - 1].size() - 1]);

        A[A.size() - 1][A.size() - 1] = ku_2 + ku_3;
        A[A.size() - 1][A.size() - 1 - 1] = ku_2 + ku_1;

        b[b.size() - 1] = ku_4 - ku_3 * p2 * h;


    }

    NM::printArray(A);
    NM::printArray(b);

}
