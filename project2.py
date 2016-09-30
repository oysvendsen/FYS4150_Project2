import pylab as pl
import astropy as ap
import sys
import os
import time

def rotation_mat(A, k, l):
    """
    function for rotation in an iteration
    """
    tau         = (A[l,l] - A[k,k])/float(2*A[k,l])
    if tau > 0:
        t       = -tau + (1 + tau**2)**0.5
    else:
        t       = -tau - (1 + tau**2)**0.5
    c           = 1./(1+t**2)
    s           = c*t

    B   =   A
    B[:,k]  = A[:,k]*c - A[:,l]*s
    B[k,:]  = B[:,k]
    B[:,l]  = A[:,l]*c + A[:,k]*s
    B[l,:]  = B[:,l]
    B[k,k]  = A[k,k]*(c**2)- 2*A[k,l]*c*s + A[l,l]*(s**2)
    B[l,l]  = A[l,l]*(c**2)- 2*A[k,l]*c*s + A[k,k]*(s**2)
    B[k,l]  = 0 
    B[l,k]  = 0
    return B

def find_args(A):
    """
    finds k and l as arguments for the maximum value of the matrix
    """
    if pl.shape(A)[0] == pl.shape(A)[1]:
        sz  = pl.shape(A)[0]
        print "sz=", sz
    else:
        sys.exit("Matrix dimensions don't match")

    eps     = A[0,1]
    
    for i in pl.arange(0, sz-1):
        for j in pl.arange(i+1, sz):
            A_val = abs(A[i,j])
            if A_val >= eps:
                max_i = i
                max_j = j
                eps = A_val

        # stores the maxima;tot.max arg val. -> tot.max row arg
        # stores args of col.s' maxes

    return max_i, max_j, eps

def make_A(n):
    """
    creates our specialized matrix given in the assignment
    """
    rho_0   = 0.
    rho_end = 10.
    rho     = pl.linspace(rho_0, rho_end, n)
    h       = (rho_end - rho_0)/float(n)

    eup     = -1./(h**2)*pl.ones(len(rho)-1)
    edown   = -1./(h**2)*pl.ones(len(rho))[1:]
    Vee     = rho**2
    di      = 2./(h**2) + Vee

    A = pl.diag(di, k=0) + pl.diag(eup, k=-1) + pl.diag(edown, k=1)
    return A

def run():
    """
    main body of the program
    """
    fasit_eig   = pl.array([3,7,11])
    drac        = 0
    epstol      = 1e-8

    n   = int(sys.argv[1]) + 1 # or something
    A   = make_A(n)
    k, l, epsloop = find_args(A)

    while abs(epsloop) > abs(epstol):
        A     = rotation_mat(A, k, l)
        k, l, epsloop = find_args(A)
        drac += 1

    print "Total number of iterations:", drac
    print pl.diag(A)

def unit_mirror(matrix):
    diff_tot = 0
    eps_tot = 0
    n = len(matrix[0]) #assume quadratic matrix
    for i in range(n-1): #rows of upper triangular
        for j in range(i+1,n):
            diff_tot += abs(matrix[i,j]-matrix[j,i])
            eps_tot += abs(matrix[i,j]-matrix[j,i])**2
def unit_known(matrix):
    None
def unit_orthogonal(matrix):
    None
run()




