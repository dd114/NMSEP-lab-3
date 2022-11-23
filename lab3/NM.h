#pragma once
#include <iostream>
#include <vector>
#include <cassert>

using namespace std;
class NM{

public:
	NM(){

	}

	static vector<double> tridiagonalSolution(const vector<vector<double>>& matrix1, const vector<double>& matrix2) {
    assert(matrix1.size() == matrix2.size() && "Sizes match");

    double y;
    int N = matrix2.size();
    int N1 = N - 1;
    vector<double> a(N), B(N), matRes(N);

    y = matrix1[0][0];
    a[0] = -matrix1[0][1] / y;
    B[0] = matrix2[0] / y;
    for (int i = 1; i < N1; i++) {
        y = matrix1[i][i] + matrix1[i][i - 1] * a[i - 1];
        a[i] = -matrix1[i][i + 1] / y;
        B[i] = (matrix2[i] - matrix1[i][i - 1] * B[i - 1]) / y;
    }

    matRes[N1] = (matrix2[N1] - matrix1[N1][N1 - 1] * B[N1 - 1]) / (matrix1[N1][N1] + matrix1[N1][N1 - 1] * a[N1 - 1]);
    for (int i = N1 - 1; i >= 0; i--) {
        matRes[i] = a[i] * matRes[i + 1] + B[i];
    }

    return matRes;

    }

    template <typename T>
    static void printArray(const vector<T>& matrix1){
        cout << endl;
        cout << "SIZE = " << matrix1.size() << endl;

        for(int i = 0; i < matrix1.size(); i++){
            cout << matrix1[i] << endl;
        }

    }

    template <typename T>
    static void printArray(const vector<vector<T>>& matrix1){
        cout << endl;
        cout << "SIZE 1 = " << matrix1.size() << endl;
        cout << "SIZE 2 = " << matrix1[0].size() << endl;

        for(int i = 0; i < matrix1.size(); i++){
            for(int j = 0; j < matrix1[i].size(); j++){
                cout << matrix1[i][j] << " ";
            }
        cout << endl;
        }
    }




};

