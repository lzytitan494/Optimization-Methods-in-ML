#include <bits/stdc++.h>
#include "matrix.h"
using namespace std;

/*
    Conjugate Gradient Descent: (moves towards conjugate direction && solves exactly in n (n -> no. of variables) iterations)
    --> 1. calculate ∇f(x) and H(x), set d_0 = 0, starting point: X_0
    --> 2. Find the search direction d_k as: (d_k * H(x) * d_k-1 = 0)
        --> d_k  = -∇f(x_k)  + B_K * d_k-1, where B_k = |∇f(x_k)|² / |∇f(x_k-1)|²
    --> 3. Determine the optimal steplength:
        --> a_k =  ∇f(x_k) * ∇f(x_k) / d_k * H(x) * d_k
    --> 4. x_k+1 = x_k + a_k * d_k
    --> 5. Test the optimality for the newpoint x_k+1:
        -->  ∇f(x_k+1) = 0
        --> If not reached, repeat step 2 to 5.
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

    cout << "Conjugate Gradient Descent: \n";

    // Initialization
    vector<vector<double>> d_k = {{0, 0}};
    vector<vector<double>> df_x_prev = {{0, 0}}; // To store the previous df_x
    double B_k = 0;

    for(int i=1; i<=itr; i++) {
        cout << "\nK: " << i << "\n";
        cout << "------------------------------\n";
        if(i>1) B_k = vector_magnitude_1xn_squared(df_x) / vector_magnitude_1xn_squared(df_x_prev);

        vector<vector<double>> b_d_k = mx_constant_multiplication(B_k, d_k);
        d_k = mx_addition(mx_constant_multiplication(-1, df_x), b_d_k);

        double alpha = mx_to_constant(mx_multiplication(df_x, mx_transpose(df_x))) /
        mx_to_constant(mx_multiplication(d_k, mx_multiplication(h_x, mx_transpose(d_k))));

        x = mx_addition(x, mx_constant_multiplication(alpha, d_k));
        
        // displaying result
        cout <<"x"<<i<<": ";
        print_matrix(x);
        cout <<"alpha"<<i<<": "<<alpha<<"\n";
        cout <<"d"<<i<<": ";
        print_matrix(d_k);
        cout <<"B_k"<<i<<": "<<B_k<<"\n";
        cout << "\n";
        
        x1 = x[0][0];
        x2 = x[0][1]; 

        df_x_prev = df_x; // storing the previous
        df_x1 = 2*a*x1 + c*x2 + d;
        df_x2 = 2*b*x2 + c*x1 + e;
        df_x = {{df_x1, df_x2}};
        cout << "∇f(x"<<i<<"): ";
        print_matrix(df_x);
        if(i>1){
            cout << "∇f(x"<<i-1<<"): "; print_matrix(df_x_prev);
        }
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