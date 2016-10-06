#include <iostream>
#include <armadillo>
#include <cmath>
#include <string>

using namespace std;
using namespace arma;

void solve_jacobi( mat& A, mat& eigvec, vec& eigval, int n );
void jacobi_rotation(mat& A, mat& R, int k, int l, int n);

void singular_electron_matrix(mat& A, vec&rho, double h, int n);
void double_electron_matrix(mat& A, vec&rho, double h, int n);
double find_largest(mat& A, int& max_i, int& max_j, int n);
int kronicker(int i,int j);

bool unit_symmetry(mat& A, int n, double tol);
bool unit_orthogonality(mat& matrix, int shape, double tol);
int unit_known();
int unit_known_arma();

int main(int argc, char *argv[])
{
    //unit tests
    //unit_known();
    unit_known_arma();
    exit(0);
    //declare variables
    double h; //step-length,
    double rho_0 = 0.0; //starting point of array
    double rho_n = 10.0; //ending point of array
    double omega = 0.0; //HO-angular frequency
    int n = 5; //size of arrays
    bool interact = false; //interacting or non-interacting electrons

    // fetch cmd-line arguments
    if (argc == 2){
        n = atoi(argv[1]);
    } else if (argc == 3) {
        n = atoi(argv[1]);
        rho_n = atof(argv[2]);
    } else if (argc == 4) {
        n = atoi(argv[1]);
        rho_0 = atof(argv[2]);
        rho_n = atof(argv[3]);
    } else if (argc == 5) {
        n = atoi(argv[1]);
        rho_0 = atof(argv[2]);
        rho_n = atof(argv[3]);
        omega = atof(argv[4]);
    }else if (argc == 6) {
        n = atoi(argv[1]);
        rho_0 = atof(argv[2]);
        rho_n = atof(argv[3]);
        omega = atof(argv[4]);
        interact = true;
    }

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

    //find data-filename
    string filename = "../data/project2";
    //filename += "_test.dat";
    if (interact) {
        filename += "_interacting";
    } else {
        filename += "_noninteracting";
    }
    filename += "_rho0=" + to_string( (int) rho_0 );
    filename += "_rho_N=" + to_string( (int) rho_n);
    filename += "_N=" + to_string( (int)n);
    if (omega != 0) {
        filename += "_omega=" + to_string( (int) (1000.0*omega));
    }
    filename += ".dat";

    ofstream outfile;
    outfile.open(filename, std::ofstream::out);

    //sort eigenvalues
    uvec indeces_sorted = sort_index(lambda); // array of indeces of lambda when sorted
    for (int i=0; i<3; i++){
        int sort_i = indeces_sorted(i);
        outfile << lambda(sort_i) << endl;
        for (int j=0; j<n; j++){
            outfile << R(j,sort_i) << ", ";
        }
        outfile << endl;
    }
    outfile.close();

    return 0;
} //END: main

void solve_jacobi( mat &A, mat &eigvec, vec &eigval, int n){
    //declare local variables
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
        //sym = unit_symmetry(A, n, 1e-10);
        //orth = unit_orthogonality(eigvec, n, 1e-10);
        eps_iter = find_largest(A, k, l, n);
        iterations++;
    }
    if (eps_iter<eps_tol){
        for (int i=0; i<n; i++){
            eigval(i) = A(i,i);
        }
    } else {
        eigval.zeros();
        cout << "Epsilon did not reach low enough levels, eigenvalues not calculated." << endl;
    }
    cout << "Number of iterations: " << iterations << "\n";
} //END: solve_jacobi

void jacobi_rotation(mat &A, mat &R, int k, int l, int n){
    /* Take a matrix A and turn it to matrix B by jacobi_rotation. */
    double tau, t, s, c, a_kk, a_ll, a_ik, a_il, r_ik, r_il;
    if (A(k,l) == 0){
        c = 1;
        s = 0;
    } else {
        tau         = (A(l,l) - A(k,k))/(2.0*A(k,l)); // find cos(2 \theta)
        if (tau >= 0){
            t       = -tau + sqrt(1 + tau*tau); // lowest root of tan(\theta)
        } else {
            t       = -tau - sqrt(1 + tau*tau); // lowest root of tan(\theta)
        }
        c           = 1.0/sqrt(1+t*t); // cos(\theta)
        s           = c*t; // sin(\theta)
    }
    a_kk = A(k,k); a_ll = A(l,l);
    A(k,k)  = a_kk*c*c - 2.0*A(k,l)*c*s + a_ll*s*s; //diagonal values
    A(l,l)  = a_ll*c*c + 2.0*A(k,l)*c*s + a_kk*s*s; //diagonal values
    A(k,l)  = 0; //set to zero
    A(l,k)  = 0; //set to zero
    for (int i=0; i<n; i++){
        if (i != k && i != l){
            a_ik = A(i,k); a_il = A(i,l);
            A(i,k)  = a_ik*c - a_il*s;
            A(k,i)  = A(i,k); //symmetric value
            A(i,l)  = a_il*c + a_ik*s;
            A(l,i)  = A(i,l); // symmetric value
        }
        // compute the new eigenvectors in R-matrix
        r_ik = R(i,k);
        r_il = R(i,l);
        R(i,k) = c*r_ik - s*r_il;
        R(i,l) = c*r_il + s*r_ik;
    }
} //END: jacobi_rotation

void singular_electron_matrix(mat& A, vec& rho, double h, int n){
    /* remake the matrix A to match a single electron with
     * length array rho of length n*/
    double di, ei, V;
    double h2_inv = 1.0/(h*h); //precalculate 1/h^2

    //set all matrix elements to zero
    for (int i=0; i<n; i++){
        for (int j=0; j<n; j++){
            A(i,j) = 0.0;
        }
    }

    //calulate diagonal d = rho^2 + 2/h^2
    for (int i=0; i<n; i++){
        V = rho(i)*rho(i);
        di = V + 2.0*h2_inv;
        A(i,i) = di;
    }

    //calculate tridiagonal e = -1/h^2
    for (int i=0; i<n-1; i++){
        ei = -1.0*h2_inv;
        A(i,i+1) = ei;
        A(i+1,i) = ei;
    }
} //END: singular_electron_matrix

void double_electron_matrix(mat &A, vec &rho, double omega_r, double h, int n){
    /* remake the matrix A to match a single electron with
     * length array rho of length n*/
    double di, ei, V;
    double h2_inv = 1.0/(h*h); //precalculate 1/h^2

    //set all matrix elements to zero
    for (int i=0; i<n; i++){
        for (int j=0; j<n; j++){
            A(i,j) = 0.0;
        }
    }

    //calulate diagonal d = rho^2 + 2/h^2
    for (int i=0; i<n; i++){
        V = omega_r*rho(i)*rho(i) - 1.0/rho(i);
        di = V + 2.0*h2_inv;
        A(i,i) = di;
    }

    //calculate tridiagonal e = -1/h^2
    for (int i=0; i<n-1; i++){
        ei = -1.0*h2_inv;
        A(i,i+1) = ei;
        A(i+1,i) = ei;
    }
} //END: double_electron_matrix


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
        for (int j=i+1; j<n; j++){
            A_val = abs(A(i,j));
            eps += 2.0*A_val;
            if (A_val >= max){
                max_i = i;
                max_j = j;
                max = A_val;
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

bool unit_orthogonality(mat& matrix, int shape, double tol){
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

int unit_known(){
    /* This function will test the imbedded jacobi solver for A=[[2,1],[1,2]].
     * The result should be:
     * lamba_1 = 1      v_1 = [1,-1]
     * lambda_2 = 3     v_2 = [1,1]
     */
    mat A = { {2, 1}, {1, 2} };
    mat R_correct = { {1, 1}, {-1, 1} };
    mat R(2,2);
    vec lambda_correct = { 1, 3 };
    vec lambda(2);
    cout << "Testing Jacobi's method on a known solution:" << endl;
    A.print("Matrix A:");
    lambda_correct.print("known solution for eigenvalues:");
    R_correct.print("known solution for eigenvectors:");

    solve_jacobi(A,R,lambda, 2);

    lambda.print("calculated solution for eigenvalues:");
    R.print("calculated solution for eigenvectors:");
    return 0;
} //END:unit_known

int unit_known_arma() {
    /* This unit test will check the compute the eigenvalues and eigenvectors
     * of a 4x4-matrix using the integrated jacobi-solver.
     * Thereafter the values will be checked against armadillos eigsys-function.
     */
    double val, tol;
    int n = 3;
    tol = 1e-10;
    //vec matrix_values = randu(4,4);
    /*
     * mat A_test(4,4);
    mat A_arma(4,4);
    for (int i=0; i<4; i++){
        for (int j=i; j<4; j++){
            val = matrix_values(i,j);
            A_test(i,j) = val;
            A_test(j,i) = val;
            A_arma(i,j) = val;
            A_arma(j,i) = val;
        }
    }
    */
    mat A_test = {{1,2,3}, {2,3,4}, {3,4,1}};
    mat A_arma = A_test;
    mat R_test(3,3);
    mat R_arma(3,3);
    vec lambda_test(3);
    vec lambda_test_sorted(3);
    vec lambda_arma(3);

    solve_jacobi(A_test, R_test, lambda_test, 3);
    //sort eigenvalues
    uvec indeces_sorted = sort_index(lambda_test); // array of indeces of lambda when sorted
    for (int i=0; i<3; i++){
        int sort_i = indeces_sorted(i);
        lambda_test_sorted(i) = lambda_test(sort_i);
    }
    eig_sym(lambda_arma, R_arma, A_arma);
    lambda_arma.print();
    lambda_test.print();
    for (int i=0; i<3; i++){
        if (abs(lambda_arma(i) - lambda_test(i)) > tol){
            cout << "unit test 'unit_known_arma' failed" << endl;
            return 1;
        }
    }
    return 0;
}
