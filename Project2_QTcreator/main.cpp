#include <iostream>
#include <cmath>
#include <armadillo>

using namespace std;
using namespace arma;

void singular_electron_matrix(mat& A, vec&rho, double h, int n);
void double_electron_matrix(mat& A, vec&rho, double h, int n);
void solve_jacobi( mat& A, mat& eigvec, vec& eigval, int n );
void jacobi_rotation(mat& A, mat& R, int k, int l, int n);
double find_largest(mat& A, int& max_i, int& max_j, int n);
int kronicker(int i,int j);
bool unit_symmetry(mat& A, int n, double tol);
double unit_orthogonality(mat* matrix, int shape, double tol);

int main(int argc, char *argv[])
{
    // fetch cmd-line arguments

    // declare variables
    double h;
    int n = 10; //size of array
    double rho_0 = 0.0; //starting-point of length-array
    double rho_n = 10.0; //end-point of length-array
    h = (rho_n-rho_0)/(n-1); //steplength
    vec rho = linspace<vec> (rho_0, rho_n, n); //length-array
    vec lambda = zeros<vec> (n); //eigenvalue-array

    // allocate initial matrix and initial matrix of eigenvalues
    mat A(n,n); //matrix
    mat R(n,n); //eigenvalues

    // start solving for single electron
    singular_electron_matrix(A,rho,h,n);

    //solve eigenvectors and eigenvalues using jacobi's method
    solve_jacobi(A, R, lambda, n);

    // start solving for two electrons
    //double_electron_matrix(A,rho,n);

    return 0;
} //END: main

void singular_electron_matrix(mat& A, vec& rho, double h, int n){
    /* remake the matrix A to match a single electron with
     * length array rho of length n*/
    double di, ei;
    double h2_inv = 1.0/(h*h); //precalculate 1/h^2

    //set all matrix elements to zero
    for (int i=0; i<n; i++){
        for (int j=0; j<n; j++){
            A(i,j) = 0.0;
        }
    }

    //calulate diagonal d = rho^2 + 2/h^2
    for (int i=0; i<n; i++){
        di = rho(i)*rho(i) + 2.0*h2_inv;
        A(i,i) = di;
    }

    //calculate tridiagonal e = -1/h^2
    for (int i=0; i<n-1; i++){
        ei = -1.0*h2_inv;
        A(i,i+1) = ei;
        A(i+1,i) = ei;
    }
} //END: singular_electron_matrix

void double_electron_matrix(mat& A, vec* rho, double h, int n){
    /* remake the matrix A to match a single electron with
     * length array rho of length n*/
    double di, ei;
    double h2_inv = 1.0/(h*h); //precalculate 1/h^2

    //set all matrix elements to zero
    for (int i=0; i<n; i++){
        for (int j=0; j<n; j++){
            A(i,j) = 0.0;
        }
    }

    //calulate diagonal d = rho^2 + 2/h^2
    for (int i=0; i<n; i++){
        di = rho(i)*rho(i) + 2.0*h2_inv;
        A(i,i) = di;
    }

    //calculate tridiagonal e = -1/h^2
    for (int i=0; i<n-1; i++){
        ei = -1.0*h2_inv;
        A(i,i+1) = ei;
        A(i+1,i) = ei;
    }
    cout << "ERROR: function double_electron_matrix not finished" << endl;
    exit(1);
} //END: double_electron_matrix

void solve_jacobi( mat& A, mat& eigvec, vec& eigval, int n){
    //declaeigvece local variables
    int k, l;
    bool sym, orth;
    double eps_iter;
    double eps_tol = 1.0e-8;
    double max_iter = n*n*n;
    int iterations = 0;
    // Setting up the eigenvector matrix as identy-matrix
    eigvec.eye();

    //check symmetry and orthogonality
    sym = unit_symmetry(A, n, 1e-10);
    orth = unit_orthogonality(eigvec, n, 1e-10);

    //calculate epsilon of off-diagonal elements (also k,l)
    eps_iter = find_largest(A, k, l, n);

    //perform jacobi's method until off-diagonal elements are close enough to zero.
    while(iterations<max_iter && sym && orth && eps_iter>eps_tol){
        //perform one jacobi rotation B=S^TAS
        jacobi_rotation(A, eigvec, k, l, n);
        //advance iterations, symmetry, orthogonality, epsilon, max-indeces(k,l)
        sym = unit_symmetry(A, n, 1e-10);
        orth = unit_orthogonality(eigvec, n, 1e-10);
        eps_iter = find_largest(A, k, l, n);
        iterations++;
    }
    if (eps_iter>eps_tol){
        for (int i=0; i<n; i++){
            eigval(i) = A(i,i);
        }
    }
    cout << "Number of iterations: " << iterations << "\n";
} //END: solve_jacobi

void jacobi_rotation(mat& A, mat& R, int k, int l, int n){
    /* Take a matrix A and turn it to matrix B by jacobi_rotation. */
    double tau, t, s, c, a_kk, a_ll, r_ik, r_il;

    tau         = (A(l,l) - A(k,k))/(2.0*A(k,l)); // find cos(2 \theta)
    if (tau >= 0){
        t       = -tau + sqrt(1 + tau*tau); // lowest root of tan(\theta)
    } else {
        t       = -tau - sqrt(1 + tau*tau); // lowest root of tan(\theta)
    }
    c           = 1.0/sqrt(1+t*t); // cos(\theta)
    s           = c*t; // sin(\theta)

    a_kk = A(k,k); a_ll = A(l,l);
    A(k,k)  = a_kk*c*c- 2.0*A(k,l)*c*s + a_ll*s*s; //diagonal values
    A(l,l)  = a_ll*c*c- 2.0*A(k,l)*c*s + a_kk*s*s; //diagonal values
    A(k,l)  = 0; //set to zero
    A(l,k)  = 0; //set to zero

    for (int i=0; i<n; i++){
        if (i != k && i != l){
            A(i,k)  = A(i,k)*c - A(i,l)*s;
            A(k,i)  = A(i,k); //symmetric value
            A(i,l)  = A(i,l)*c + A(i,k)*s;
            A(l,i)  = A(i,l); // symmetric value
        }
        // compute the new eigenvectors in R-matrix
        r_ik = R(i,k);
        r_il = R(i,l);
        R(i,k) = c*r_ik - s*r_il;
        R(i,l) = c*r_il + s*r_ik;
    }
} //END: jacobi_rotation

double find_largest(mat& A, int& max_i, int& max_j, int n){
    /* Take matrix A and calculate the difference squared between the upper triangle
     * and the lower triangle, return difference of non-diagonal values.
     */
    //declare local variables
    double A_val, max, eps;
    eps = 0;

    // start with largest value being first value of upper tridiagonal
    max = A[0,1];
    max_i = 0;
    max_j = 1;

    //loop through upper tridiagonal
    for (int i=0; i<n-1; i++){
        for (int j=i+1; j<=n; j++){
            A_val = abs(A(i,j));
            eps += 2.0*A_val*A_val
            if (A_val >= max){
                max_i = i;
                max_j = j;
            }
        }
    }
    return eps;
} //END: find_largest

int kronicker(int i,int j){
    if (i == j){
        return 1;
    } else {
        return 0;
    }
}

bool unit_symmetry(mat& A, int n, double tol){
    // return relative error in symmetry per matrix-element
    double diff;
    double eps = 0;
    //loop over upper triangular
    for (int i=0; i<n-1; i++){
        for (int j=i+1; j<n; j++){
            diff = A(i,j) - A(j,i);
            eps += 4.0*diff*diff;
        }
    }
    eps /= (n*n);
    if (eps <= tol){
        return true;
    } else {
        cout << "ERROR in matrix A: " << endl
             << "Matrix is not symmetric, deviation larger then tolerance." << endl
             << "(sum(a_ij-aji)^2)/(n^2)="
             << eps << endl;
        return false;
    }
} //END: unit_symmetry

bool unit_orthogonality(mat* matrix, int shape, double tol){
    /*
        Takes matrix of all eigenvectors, check dot-product of all vectors.
        Create new matrix of V_i^T*V_j that should be delta_ij.
        Return difference squared of each element between these two matrices.
        */
    double delta, diff;
    double* V_i, V_j;
    double eps = 0;
    for (int i=0; i<shape; i++){
        for (int j=0; j<shape; j++){
            //V_i = matrix.colptr(i);
            //V_j = matrix.colptr(j);
            //delta = dot(V_i,V_j);
            //diff = delta - kronicker(i,j);
            //eps += diff*diff
        }
    }
    if (eps <= tol){
        return true;
    } else {
        cout << "ERROR in eigenvectors: " << endl
             << "vectors are not orthogonal, deviation larger then tolerance." << endl
             << "sum(V_j^T * V_j - kronicker_delta(i,j)^2="
             << eps << endl;
        return false;
    }
} //END: unit_orthogonality

double unit_known(){} //END:unit_known
