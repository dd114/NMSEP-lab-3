#pragma once
#include <iostream>
#include <functional>
#include <vector>
#include "NM.h"

using namespace std;
class SolutionOfEquations {
    double L, T;
    double a;
    double beta;
    int numberOfPointsL, numberOfPointsT;

    function<double (double)> phi;
    function<double (double)> management;


public:
	SolutionOfEquations(double L, double T, double a, double beta, int numberOfPointsL, int numberOfPointsT, double (*phi)(double), double (*management)(double)) {
        this->L = L; 
        this->T = T;
        this->a = a;
        this->beta = beta;
        this->numberOfPointsL = numberOfPointsL;
        this->numberOfPointsT = numberOfPointsT;

        this->phi = phi;
        this->management = management;
	}

	vector<vector<double>> straightTask() {

        double a2 = a * a;
        double h = L / (numberOfPointsL - 1);
        double tau = T / (numberOfPointsT - 1);
        vector<vector<double>> u(numberOfPointsT, vector<double>(numberOfPointsL));

        for (int i = 0; i < u[0].size(); i++) {
            u[0][i] = phi(h * i);
        }

        vector<vector<double>> A(numberOfPointsL - 2, vector<double>(numberOfPointsL - 2));
        vector<double> b(numberOfPointsL - 2);

        double ku_1 = a2 / (h * h);
        double ku_2 = -((2. * a2) / (h * h) + 1. / tau);
        double ku_3 = ku_1;

        double ku_4 = -1. / tau;


        for (int i = 1; i < numberOfPointsT; i++) {

            double p1 = 0;

            A[0][0] = ku_2 + ku_1;
            A[0][1] = ku_3;

            b[0] = ku_4 * u[i - 1][1] + ku_1 * p1 * h;


            for (int j = 1; j < A.size() - 1; j++) {

                A[j][j - 1] = ku_1;
                A[j][j] = ku_2;
                A[j][j + 1] = ku_3;

                b[j] = ku_4 * u[i - 1][j + 1]; // j + 1, тк нужное значение на временном слое i - 1, смещенно на 1 в u из-за граничной точки


            }


            //double p2 = derivative_U_X_end(i * tau, u[i - 1][u[i - 1].size() - 1]);

            A[A.size() - 1][A.size() - 1] = ku_2 + ku_3 * (1. / (1. + beta * h));
            A[A.size() - 1][A.size() - 1 - 1] = ku_1;

            b[b.size() - 1] = ku_4 * u[i - 1][u[i - 1].size() - 1 - 1] - ku_3 * ((1. / (1. + beta * h)) * beta * h * management(i * tau));

            vector<double> tempU = NM::tridiagonalSolution(A, b);




            for (int j = 0; j < tempU.size(); j++) {
                u[i][j + 1] = tempU[j];
            }

            u[i][0] = u[i][1] - p1 * h;
            u[i][u[i].size() - 1] = u[i][u[i].size() - 1 - 1] * (1. / (1. + beta * h)) + beta * h * management(i * tau);



            //NM::printArray(tempU);
            //NM::printArray(A);
        }


        return u;
	}


};

