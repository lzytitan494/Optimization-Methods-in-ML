#include <bits/stdc++.h>
#include "matrix.h"
using namespace std;

/*
    Newton's Method:
    It finds the minimum of a quadratic function in a single step. (assumption: h_X is not singular and ∇f(x) != 0)
    --> X_k+1 = X_k - H(x_k)⁻¹ * ∇f(x_k)
*/

int main() {
    // A general 2 variable quadratic function - f(x) = a*x1*x1 + b*x2*x2 + c*x1*x2 + d*x1 + e*x2 + f
    // Function to be minimized - f(x) = 2*x1*x1 + x2*x2 + 2*x1*x2 + 1*x1 + -1*x2 + 0
    double a, b, c, d, e, f;
    a = 2;
    b = 1;
    c = 2;
    d = 1;
    e = -1;
    f = 0;

    double x1, x2;
    // Starting value
    vector<vector<double>> x = {{0, 0}};
    x1 = x[0][0];
    x2 = x[0][1];

    // Objective function
    double f_x = a*x1*x1 + b*x2*x2 + c*x1*x2 + d*x1 + e*x2 + f;

    // Partial derivative x1 - f'(x) = 2*a*x1 + c*x2 + d
    // Partial derivative x2 - f'(x) = 2*b*x2 + c*x1 + e
    double df_x1 = 2*a*x1 + c*x2 + d;
    double df_x2 = 2*b*x2 + c*x1 + e;

    // hessian matrix = [[2*a c][c 2*b]]
    vector<vector<double>> h_x = {{2*a, c}, {c, 2*b}};
    
    // Gradient matrix
    vector<vector<double>> df_x = {{df_x1, df_x2}};


    cout << "Newton's Method: \n";
    cout << "\n------------------------------\n";
    x = mx_addition(x,mx_constant_multiplication(-1, mx_multiplication(df_x, mx_inverse_2x2(h_x))));
    cout <<"x: ";
    print_matrix(x);
    cout<<"\n";
    x1 = x[0][0];
    x2 = x[0][1]; 

    df_x1 = 2*a*x1 + c*x2 + d;
    df_x2 = 2*b*x2 + c*x1 + e;
    df_x = {{df_x1, df_x2}};
    cout << "∇f(x): ";
    print_matrix(df_x);
    cout <<"------------------------------\n";

    f_x = a*x1*x1 + b*x2*x2 + c*x1*x2 + d*x1 + e*x2 + f;
    cout << "\nMin value is: " << f_x << "\n";

    return 0;
}