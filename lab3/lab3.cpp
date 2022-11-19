#include <iostream>
#include "NM.h"

float phi(float x){
    return 0;
}

float trueManagement(float t){
    return 1000 * t * t * (1 - t);
}

float derivative_U_X_end(float t, float lastUInL){
    float beta = 4;
    return beta * (sin(t / 1000) - lastUInL);
}

int main(){

    float L = 1, T = 1;
    float a = 2;
    int numberOfPointsL = 100, numberOfPointsT = 100;






    float  a2 = a * a;
    float h = L / (numberOfPointsL - 1);
    float tau = L / (numberOfPointsT - 1);
    vector<vector<float>> u(numberOfPointsT, vector<float>(numberOfPointsL));

    for(int i = 0; i < u[0].size(); i++){
        u[0][i] = phi(h * i);
    }

    vector<vector<float>> A(numberOfPointsL - 1, vector<float>(numberOfPointsL - 1));
    vector<float> b(numberOfPointsL - 1);

    float ku_1 = a2 / (h * h);
    float ku_2 = -((2. * a2) / (h * h) + 1. / tau);
    float ku_3 = ku_1;

    float ku_4 = -1. / tau;


    for(int i = 1; i < numberOfPointsT; i++){

        float p1 = 0;

        A[0][0] = ku_2 + ku_1;
        b[0] = ku_4 + ku_1 * p1 * h;



        float p2 = derivative_U_X_end(i * tau, u[i - 1][u[i - 1].size() - 1]);

        A[A.size() - 1][A.size() - 1] = ku_2 + ku_3;
        b[b.size() - 1] = ku_4 - ku_3 * p2 * h;



            for(int j = 1; j < (numberOfPointsL - 1); j++){
        
                A[]

            

                float b = u[i - 1][j] / tau;



            }




    }

}
