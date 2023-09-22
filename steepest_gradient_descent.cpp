#include <bits/stdc++.h>
#include "matrix.h"
using namespace std;

/*
    Steepest Gradient Descent:
    1. Start at a point X_k.
    2. Compute the function's gradient at x_k, ∇f(x_k).
    3. Take a step in the opposite direction of the gradient, 
    x_k+1 = x_k + alpha * d_k -- alpha: step-size, d_k = - ∇f(x_k)
    4. Repeat steps 2 and 3 until convergence.
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

    // Iteration count
    int itr = 10;

    cout << "Steepest Gradient Descent: \n";

    for(int i=1; i<=itr; i++) {
        cout << "\nK: " << i << "\n";
        cout << "------------------------------\n";
        vector<vector<double>> d_k = mx_constant_multiplication(-1, df_x);
        double alpha = mx_to_constant(mx_multiplication(d_k, mx_transpose(d_k))) /
        mx_to_constant(mx_multiplication(d_k, mx_multiplication(h_x, mx_transpose(d_k))));
        x = mx_addition(x, mx_constant_multiplication(alpha, d_k));
        
        // displaying result
        cout <<"x"<<i<<": ";
        print_matrix(x);
        cout <<"alpha: "<<alpha<<"\n";
        cout <<"d"<<i<<": ";
        print_matrix(d_k);
        cout << "\n";
        
        x1 = x[0][0];
        x2 = x[0][1]; 

        df_x1 = 2*a*x1 + c*x2 + d;
        df_x2 = 2*b*x2 + c*x1 + e;
        df_x = {{df_x1, df_x2}};
        cout << "∇f(x"<<i<<"): ";
        print_matrix(df_x);
        cout <<"------------------------------\n";
        // optimality
        if((df_x1 == 0 && df_x2 == 0) || (abs(df_x1) + abs(df_x2) <= 0.01)) {
            f_x = a*x1*x1 + b*x2*x2 + c*x1*x2 + d*x1 + e*x2 + f;
            cout << "\nOptimality Reached\n"; 
            cout << "Min value is: " << f_x << "\n";
            break;
        }
    }

    return 0;
}