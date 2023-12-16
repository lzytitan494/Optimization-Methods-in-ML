#include <bits/stdc++.h>
#include 'matrix.h'
using namespace std;

/* 
    Drawback's of newtons'method:
    --> 1. Needs computation of H⁻¹(x).
    --> 2. H(x) needs to be positive definite.
    --> 3. Prone to local convergence.
*/

/*
    Quasi Newton Method:
    --> 1. calculate g_k (g_k = ∇f(x)) and H(x), set d_0 = 0, starting point: X_0, set B = In (Identity matrix), set tolerance
    --> 2. while ||gk|| > tolerance: 
        --> d_k = - B_k * g_k (here, g_k ==> (1x2)) -- reason for transpose
        --> a_k = g_k * g_k.T / d_k * H(x_k) * d_k.T
    --> 3. Calculate next value:
        --> x_k+1 = x_k + a_k * d_k
        --> s_k = x_k+1 - x_k
        --> y_k = g_k+1 - g_k
        --> DFP (Davidon - Fletcher - Powell) udpdate:
            --> B_k+1 = B_k + s_k.T * s_k / s_k * y_k.T - (B_k * y_k).T * (B_k * y_k) / y_k * B_k * y_k.T   
    --> 4. Until tolerance met or iteration constraint reached repeat the process
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

    cout << "Quasi Newton Method: \n";

    // Initialization
    vector<vector<double>> B_k = {{1, 0}, {0, 1}}; // Identity matrix of 2x2
    vector<vector<double>> d_k = {{0, 0}};
    vector<vector<double>> g_k = df_x; // only for books references - choose whatever you like
    vector<vector<double>> g_k_prev = {{0, 0}}; // To store the previous g_k
    vector<vector<double>> s_k = {{0, 0}};
    vector<vector<double>> y_k = {{0, 0}};

    for(int i=1; i<=itr; i++) {
        cout << "\nK: " << i << "\n";
        cout << "------------------------------\n";
        d_k = mx_constant_multiplication(-1, mx_multiplication(g_k, B_k));

        double alpha = mx_to_constant(mx_multiplication(g_k, mx_transpose(g_k))) /
        mx_to_constant(mx_multiplication(d_k, mx_multiplication(h_x, mx_transpose(d_k))));

        x = mx_addition(x, mx_constant_multiplication(alpha, d_k));
        
        // displaying result
        cout <<"x"<<i<<": ";
        print_matrix(x);
        cout <<"alpha"<<i<<": "<<alpha<<"\n";
        cout <<"d"<<i<<": ";
        print_matrix(d_k);
        cout << "\n";
        
        // update values
        x1 = x[0][0];
        x2 = x[0][1]; 

        g_k_prev = g_k;

        df_x1 = 2*a*x1 + c*x2 + d;
        df_x2 = 2*b*x2 + c*x1 + e;
        g_k = {{df_x1, df_x2}};
        cout << "g(k"<<i-1<<"): ";
        print_matrix(g_k_prev);
        cout << "g(k"<<i<<"): ";
        print_matrix(g_k);
        cout << "\n";


        s_k = mx_constant_multiplication(alpha, d_k); // x_k+1 - x_k = a_K * d_k
        y_k = mx_addition(g_k, mx_constant_multiplication(-1, g_k_prev));

        B_k = mx_addition(B_k, mx_addition(mx_constant_multiplication(1 / mx_to_constant(mx_multiplication(s_k, mx_transpose(y_k))), mx_multiplication(mx_transpose(s_k), s_k)),
                           mx_constant_multiplication( -1 / mx_to_constant(mx_multiplication(y_k, mx_multiplication(B_k, mx_transpose(y_k)))), mx_multiplication(mx_transpose(mx_multiplication(y_k, B_k)), mx_multiplication(y_k, B_k)))));

        cout << "s(k"<<i<<"): ";
        print_matrix(s_k);
        cout << "y(k"<<i<<"): ";
        print_matrix(y_k);
        cout << "B(k"<<i<<"): \n";
        print_matrix(B_k);
        cout <<"------------------------------\n";

        // optimality
        if(vector_magnitude_1xn_squared(g_k) < 0.001 ) {
            f_x = a*x1*x1 + b*x2*x2 + c*x1*x2 + d*x1 + e*x2 + f;
            cout << "Optimality Reached\n"; 
            cout << "Min value is: " << f_x << "\n";
            break;
        }
    }

    return 0;
}