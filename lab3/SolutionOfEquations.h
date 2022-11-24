#pragma once
#include <iostream>
#include <functional>
#include <vector>
#include <fstream>
#include <cmath>
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

    double p_min = 0, p_max = 4000. / 27.;

    //vector<vector<double>> u;
    //vector<vector<double>> psi;

    vector<double> initialManagement;

    //vector<double> y;

    bool newU = false;


public:
	SolutionOfEquations(double L, double T, double a, double beta, int numberOfPointsL, int numberOfPointsT, function<double(double)> phiFunction, function<double(double)> managementFunction) {
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


        initialManagement = vector<double>(numberOfPointsT);

        for (int i = 0; i < numberOfPointsT; i++) {
            initialManagement[i] = managementFunction(i * tau);
        }

	}

    vector<vector<double>> straightTask(vector<double> currentManagement = {}) {

        if (currentManagement.empty())
            currentManagement = initialManagement;
        

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

            A[A.size() - 1][A.size() - 1 - 1] = ku_1;
            A[A.size() - 1][A.size() - 1] = ku_2 + ku_3 * (1. / (1. + beta * h));

            b[b.size() - 1] = ku_4 * u[i - 1][u[i - 1].size() - 1 - 1] - ku_3 * ((1. / (1. + beta * h)) * beta * h * currentManagement[i]);

            vector<double> tempU = NM::tridiagonalSolution(A, b);




            for (int j = 0; j < tempU.size(); j++) {
                u[i][j + 1] = tempU[j];
            }

            u[i][0] = u[i][1] - p1 * h;
            u[i][u[i].size() - 1] = u[i][u[i].size() - 1 - 1] * (1. / (1. + beta * h)) + beta * h * currentManagement[i];



            //NM::printArray(tempU);
            //NM::printArray(A);
        }


        return u;
	}


    vector<vector<double>> conjugateTask(const  vector<vector<double>>& u, const vector<double>& y) {

        vector<vector<double>> psi = vector<vector<double>>(numberOfPointsT, vector<double>(numberOfPointsL));


        //this->y = y;

        for (int j = 0; j < psi[psi.size() - 1].size(); j++) {
            psi[psi.size() - 1][j] = 2 * (u[u.size() - 1][j] - y[j]);
        }

        vector<vector<double>> A(numberOfPointsL - 2, vector<double>(numberOfPointsL - 2));
        vector<double> b(numberOfPointsL - 2);

        double ku_1 = -a2 / (h * h);
        double ku_2 = -((-2. * a2) / (h * h) - 1. / tau);
        double ku_3 = ku_1;

        double ku_4 = 1. / tau;


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

            A[A.size() - 1][A.size() - 1 - 1] = ku_1;
            A[A.size() - 1][A.size() - 1] = ku_2 + ku_3 * (1. / (1. + beta * h));

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

        return psi;

    }



    vector<double> calculation(const vector<double>& y, double epsilon = 0.01) {

        vector<double> currentManagement = initialManagement;

        double nowL2Norm = 1e+10;
        int numberOfIteration = 0;

        while(nowL2Norm > epsilon) {

            numberOfIteration++;

            vector<vector<double>> u = straightTask(currentManagement);

            nowL2Norm = l2Norm(u[u.size() - 1], y, h);

            cout << "i = " << numberOfIteration << " nowL2Norm = " << nowL2Norm << endl;


            vector<vector<double>> psi = conjugateTask(u, y);

            vector<double> tempManagement(psi.size());

            for (int j = 0; j < tempManagement.size(); j++) {
                if (psi[j][psi[j].size() - 1] >= 0)
                    tempManagement[j] = p_min;
                else
                    tempManagement[j] = p_max;
            }


            vector<double> integrandOfNumerator(psi.size());

            for (int j = 0; j < integrandOfNumerator.size(); j++) {
                integrandOfNumerator[j] = a2 * beta * psi[j][psi[j].size() - 1] * (tempManagement[j] - currentManagement[j]);
            }

            double numerator = rectangleSquare(integrandOfNumerator, tau);




            vector<double> integrandOfDenominator(u[u.size() - 1].size());

            vector<vector<double>> tempU = straightTask(tempManagement);

            for (int j = 0; j < integrandOfDenominator.size(); j++) {
                integrandOfDenominator[j] = (tempU[tempU.size() - 1][j] - u[u.size() - 1][j]) * (tempU[tempU.size() - 1][j] - u[u.size() - 1][j]);
            }

            double denominator = rectangleSquare(integrandOfDenominator, h);


            double currentAlfa = min(-0.5 * (numerator / denominator), 1.);

            //currentAlfa = 0.001;

            //if (nowL2Norm < 600) {
            //    currentAlfa = 0.00001;
            //}

            for (int j = 0; j < currentManagement.size(); j++) {
                currentManagement[j] = currentManagement[j] + currentAlfa * (tempManagement[j] - currentManagement[j]);
            }


        }

        vector<vector<double>> answerU = straightTask(currentManagement);

        return answerU[answerU.size() - 1];

    }

    double rectangleSquare(const vector<double>& vec, const double step) &{
        double sum = 0;
        for (int i = 0; i < vec.size(); i++) {
            sum += vec[i] * step;
        }

        return sum;
    }

    double l2Norm(const vector<double>& vec, const vector<double>& y, const double step) {
        vector<double> tempVec(vec.size());
        for (int i = 0; i < tempVec.size(); i++) {
            tempVec[i] = (vec[i] - y[i]) * (vec[i] - y[i]);
        }

        return rectangleSquare(tempVec, step);
    }

    template <typename T>
    void printFile(const vector<T>& u, string nameOfFile = "ySolution.txt") {
        fstream solution(nameOfFile, ios::out);

        for (int i = 0; i < u.size(); i++) {
            solution << i * h << " " << u[i] << endl;
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

