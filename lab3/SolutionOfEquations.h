#pragma once
#include <iostream>
#include <functional>
#include <vector>
#include <fstream>
#include "NM.h"

using namespace std;
class SolutionOfEquations {
    double L, T;
    double a, a2;
    double beta;
    int numberOfPointsL, numberOfPointsT;

    double h;
    double tau;

    function<double (double)> phi;
    function<double (double)> managementFunction;

    //vector<vector<double>> u;
    //vector<vector<double>> psi;

    vector<double> management;

    //vector<double> y;

    bool newU = false;


public:
	SolutionOfEquations(double L, double T, double a, double beta, int numberOfPointsL, int numberOfPointsT, double (*phiFunction)(double), double (*managementFunction)(double)) {
        this->L = L; 
        this->T = T;
        this->a = a;
        this->beta = beta;
        this->numberOfPointsL = numberOfPointsL;
        this->numberOfPointsT = numberOfPointsT;

        this->phi = phiFunction;
        this->managementFunction = managementFunction;

        a2 = a * a;
        h = L / (numberOfPointsL - 1);
        tau = T / (numberOfPointsT - 1);


        management = vector<double>(numberOfPointsT);

        for (int i = 0; i < numberOfPointsT; i++) {
            management[i] = managementFunction(i * tau);
        }

	}

	vector<vector<double>> straightTask() {

        vector<vector<double>> u = vector<vector<double>>(numberOfPointsT, vector<double>(numberOfPointsL));

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

            b[b.size() - 1] = ku_4 * u[i - 1][u[i - 1].size() - 1 - 1] - ku_3 * ((1. / (1. + beta * h)) * beta * h * management[i]);

            vector<double> tempU = NM::tridiagonalSolution(A, b);




            for (int j = 0; j < tempU.size(); j++) {
                u[i][j + 1] = tempU[j];
            }

            u[i][0] = u[i][1] - p1 * h;
            u[i][u[i].size() - 1] = u[i][u[i].size() - 1 - 1] * (1. / (1. + beta * h)) + beta * h * management[i];



            //NM::printArray(tempU);
            //NM::printArray(A);
        }

        newU = true;

        return u;
	}


    vector<vector<double>> conjugateTask(const  vector<vector<double>>& u, const vector<double>& y) {

        if (!newU)
            throw("U isn't ready for calculations");

        vector<vector<double>> psi = vector<vector<double>>(numberOfPointsT, vector<double>(numberOfPointsL));


        //this->y = y;

        for (int j = 0; j < psi[psi.size() - 1].size(); j++) {
            psi[psi.size() - 1][j] = 2 * (u[u.size() - 1][j] - y[j]);
        }

        vector<vector<double>> A(numberOfPointsL - 2, vector<double>(numberOfPointsL - 2));
        vector<double> b(numberOfPointsL - 2);

        double ku_1 = - a2 / (h * h);
        double ku_2 = -( - (2. * a2) / (h * h) + 1. / tau);
        double ku_3 = ku_1;

        double ku_4 = -1. / tau;


        for (int i = (numberOfPointsT - 1) - 1; i >= 0; i--) {

            double p1 = 0;

            A[0][0] = ku_2 + ku_1;
            A[0][1] = ku_3;

            b[0] = ku_4 * psi[i + 1][1] + ku_1 * p1 * h;


            for (int j = 1; j < A.size() - 1; j++) {

                A[j][j - 1] = ku_1;
                A[j][j] = ku_2;
                A[j][j + 1] = ku_3;

                b[j] = ku_4 * psi[i + 1][j + 1]; // j + 1, тк нужное значение на временном слое i - 1, смещенно на 1 в u из-за граничной точки


            }


            //double p2 = derivative_U_X_end(i * tau, u[i - 1][u[i - 1].size() - 1]);

            A[A.size() - 1][A.size() - 1] = ku_2 + ku_3 * (1. / (1. + beta * h));
            A[A.size() - 1][A.size() - 1 - 1] = ku_1;

            b[b.size() - 1] = ku_4 * psi[i + 1][psi[i + 1].size() - 1 - 1];

            vector<double> tempPsi = NM::tridiagonalSolution(A, b);




            for (int j = 0; j < tempPsi.size(); j++) {
                psi[i][j + 1] = tempPsi[j];
            }

            psi[i][0] = psi[i][1] - p1 * h;
            psi[i][psi[i].size() - 1] = psi[i][psi[i].size() - 1 - 1] * (1. / (1. + beta * h));



            //NM::printArray(tempU);
            //NM::printArray(A);
        }

        newU = false;

        return psi;

    }



    void calculation(const vector<double>& y) {


    }






    //void printFileLastLayerU(string nameOfFile = "ySolution.txt") {
    //    fstream solution(nameOfFile, ios::out);

    //    for (int i = 0; i < u[u.size() - 1].size(); i++) {
    //        solution << i * h << " " << u[u.size() - 1][i] << endl;
    //    }

    //    solution.close();
    //}

    //void printFileU(string nameOfFile = "trueSolution.txt") {
    //    fstream solution(nameOfFile, ios::out);

    //    for (int i = 0; i < u.size(); i++) {
    //        for (int j = 0; j < u[i].size(); j++) {
    //            solution << i * tau << " " << j * h << " " << u[i][j] << endl;
    //        }
    //    }

    //    solution.close();
    //}


    //void printFilePsi(string nameOfFile = "PsiSolution.txt") {
    //    fstream solution(nameOfFile, ios::out);

    //    for (int i = 0; i < psi.size(); i++) {
    //        for (int j = 0; j < psi[i].size(); j++) {
    //            solution << i * tau << " " << j * h << " " << psi[i][j] << endl;
    //        }
    //    }

    //    solution.close();
    //}

    template <typename T>
    void printFile(const vector<T>& u, string nameOfFile = "ySolution.txt") {
        fstream solution(nameOfFile, ios::out);

        for (int i = 0; i < u[u.size() - 1].size(); i++) {
            solution << i * h << " " << u[u.size() - 1][i] << endl;
        }

        solution.close();
    }

    template <typename T>
    void printFile(const vector<vector<T>>& u, string nameOfFile = "trueSolution.txt") {
        fstream solution(nameOfFile, ios::out);

        for (int i = 0; i < u.size(); i++) {
            for (int j = 0; j < u[i].size(); j++) {
                solution << i * tau << " " << j * h << " " << u[i][j] << endl;
            }
        }

        solution.close();
    }

};

