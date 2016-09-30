#include <iostream>
#include <armadillo>

using namespace std;
using namespace arma;

void jacobi_rotation();
void find_largest();
double unit_symmetry();
double unit_orthogonality();
double unit_known();

int main(int argc, char *argv[])
{
    // declare variables
    // fetch cmd-line arguments
    // allocate initial matrix of basis and initial matrix of eigenvalue(starts equal to basis)
    // iterate through rotations
    // fetch eigenvalues
    // calculate eigenvectors
    // sort eigenvalues with appropriate eigenvalues
    return 0;
} //END: main

void solve_jacobi(){
} //END: solve_jacobi

void jacobi_rotation(mat& A, int k, int l){
    /* Take a matrix A and turn it to matrix B by jacobi_rotation. */
    double tau, t, s, c, a_kk, a_ll;
    int n;

    if (A.shape()[0] == A.shape()[1]){
        n  = A.shape()[0]
    } else {
        exit("Matrix dimensions don't match")
    }

    tau         = (A[l,l] - A[k,k])/(2.0*A[k,l]); // find cos(2 \theta)
    if (tau > 0){
        t       = -tau + sqrt(1 + tau*tau); // lowest root of tan(\theta)
    } else {
        t       = -tau - sqrt(1 + tau*tau); // lowest root of tan(\theta)
    }
    c           = 1.0/(1+t*t); // cos(\theta)
    s           = c*t; // sin(\theta)

    for (int i=0; i<n; i++){
        if (i != k or i != l){
            A[i,k]  = A[i,k]*c - A[i,l]*s;
            A[k,i]  = A[i,k]; //symmetric value
            A[i,l]  = A[i,l]*c + A[i,k]*s;
            A[l,i]  = A[i,l]; // symmetric value
        }
    }

    a_kk = A[k,k]; a_ll = A[l,l];
    A[k,k]  = a_kk*c*c- 2.0*A[k,l]*c*s + a_ll*s*s; //diagonal values
    A[l,l]  = a_ll*c*c- 2.0*A[k,l]*c*s + a_kk*s*s; //diagonal values
    A[k,l]  = 0; //set to zero
    A[l,k]  = 0; //set to zero
} //END: jacobi_rotation

double find_largest(mat& A, int& max_i, int& max_j, double& eps){
    /* Take matrix A and calculate the difference squared between the upper triangle
     * and the lower triangle, return difference.
     */
    if (A.shape()[0] == A.shape()[1]){
        int n  = A.shape()[0]
    } else {
        exit("Matrix dimensions don't match")
    }

    double A_val, diff_max
    diff_max = A[0,1];
    max_i = 0;
    max_j = 1
    eps = 0;

    for (int i=0; i<n; i++){
        for (int j=i+1; j<=n; j++){
            A_val = abs(A[i,j]);
            if (A_val >= diff_max){
                max_i = i;
                max_j = j;
                diff_max = A_val;
            }
            eps += 2*A_val*A_val;
        }
    }
    return diff_max;
} //END: find_largest
double unit_symmetry(){} //END: unit_symmetry
double unit_orthogonality(){} //END: unit_orthogonality
double unit_known(){} //END:unit_known
