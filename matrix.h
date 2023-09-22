#include <iostream>
#include <vector>
using namespace std;

// Transpose of a matrix
vector<vector<double>> mx_transpose(vector<vector<double>> m) {
    int rows = m.size();
    int cols = m[0].size();

    vector<vector<double>> result(cols, vector<double> (rows));
    
    for(int i=0; i<rows; i++) {
        for(int j=0; j<cols; j++) {
            result[j][i] = m[i][j];
        }
    }

    return result;
}

// Matrix multiplication
vector<vector<double>> mx_multiplication(vector<vector<double>> m1, vector<vector<double>> m2) {
    int rows1 = m1.size();
    int cols1 = m1[0].size();

    int rows2 = m2.size();
    int cols2 = m2[0].size();

    vector<vector<double>> result(rows1, vector<double> (cols2));

    for(int i=0; i<rows1; i++) {
        for(int j=0; j<cols2; j++) {
            for(int k=0; k<cols1; k++) {
                result[i][j] += m1[i][k] * m2[k][j];
            }
        }
    }

    return result;
}

// Multiply a constant with a matrix
vector<vector<double>> mx_constant_multiplication(double k, vector<vector<double>> m) {
    int rows = m.size();
    int cols = m[0].size();

    vector<vector<double>> result(rows, vector<double> (cols));
    
    for(int i=0; i<rows; i++) {
        for(int j=0; j<cols; j++) {
            result[i][j] = k * m[i][j];
        }
    }

    return result;
}

// Addition of two matrix
vector<vector<double>> mx_addition(vector<vector<double>> m1, vector<vector<double>> m2) {
    int rows1 = m1.size();
    int cols1 = m1[0].size();

    int rows2 = m2.size();
    int cols2 = m2[0].size();

    if(rows1 == rows2 && cols1 == cols2) {
        vector<vector<double>> result(rows1, vector<double> (cols1));

        for(int i=0; i<rows1; i++) {
            for(int j=0; j<cols2; j++) {
                result[i][j] = m1[i][j] + m2[i][j];
            }   
        }

        return result;
    } else {
        return {{0, 0}};
    }   
}

// Determinant of 2x2 matrix
double mx_determinant_2x2(vector<vector<double>> m) {
    return m[0][0] * m[1][1] - m[0][1] * m[1][0];
}

// Inverse of 2x2 matrix
vector<vector<double>> mx_inverse_2x2(vector<vector<double>> m) {
    if(mx_determinant_2x2(m) == 0) {
        return {{0, 0}};
    } else {
        int rows = m.size();
        int cols = m[0].size();
        vector<vector<double>> result(rows, vector<double> (cols));
        result[0][0] = m[1][1];
        result[0][1] = -1 * m[0][1];
        result[1][0] = -1 * m[1][0];
        result[1][1] = m[0][0];

        return mx_constant_multiplication(1 / mx_determinant_2x2(m), result);
    }
}

// 1x1 matrix to constant
double mx_to_constant(vector<vector<double>> m) {
    if(m.size() == 1 && m[0].size() == 1) {
        return m[0][0];
    } else {
        return 0;
    }
}

// magnitude of a vector
double vector_magnitude_1xn_squared(vector<vector<double>> m) {
    int n = m[0].size();

    double result = 0;

    for(int i=0; i<n; i++) {
        result += m[0][i] * m[0][i];
    }

    return result;
}

// print matrix
void print_matrix(vector<vector<double>> m) {
    int rows = m.size();
    int cols = m[0].size();
    
    for(int i=0; i<rows; i++) {
        for(int j=0; j<cols; j++) {
            cout << m[i][j] << " ";
        }
        cout << "\n";
    }
}