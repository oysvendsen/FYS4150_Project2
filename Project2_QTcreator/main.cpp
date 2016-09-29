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

void jacobi_rotation(mat& A, int k, int l){
    /* Take a matrix A and turn it to matrix B by jacobi_rotation. */
    // declare local variables
    // tau =
    // t = //t^2 + 2tau*t - 1 = 0
    // c =
    // s =
    // for (i=1;...
    // b_ik =
    // b_il =
    // b_kk =
    // b_ll =
    // overwrite A
} //END: jacobi_rotation
void find_largest(){
    /* Take matrix A and calculate the difference squared between the upper triangle
     * and the lower triangle, return difference.
     */
    //declare local variables
} //END: find_largest
double unit_symmetry(){} //END: unit_symmetry
double unit_orthogonality(){} //END: unit_orthogonality
double unit_known(){} //END:unit_known
