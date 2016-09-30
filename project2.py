import pylab as pl
import astropy as ap
import sys
import os
import time

# def rotation(diag, eup, edown, k, l):
#     """
#     performs the iteration in question
#     """
#     new_diag        = diag
#     new_eup         = eup
#     new_edown       = edown

#     tau             = (diag[l,l] - diag[k,k])/float(2*eup[l])
#     t               = -tau + (1 + tau**2)**0.5
#     c               = 1./(1+t**2)

#     new_edown[k]    = edown[k]*pl.cos(theta) \
#                         - eup[l]*pl.sin(theta)
#     new_eup[l]      = eup[l]*pl.cos(theta) \
#                         + edown[k]*pl.sin(theta)
#     new_diag[k]     = diag[k]*(pl.cos(theta)**2) \
#                         - 2.*eup[l]*pl.cos(theta)*pl.sin(theta) \
#                         + diag[l]*(pl.sin(theta)**2)
#     new_diag[l]     = diag[l]*(pl.cos(theta)**2) \
#                         + 2.*eup[l]*pl.cos(theta)*pl.sin(theta) \
#                         + diag[k]*(pl.sin(theta)**2)
#     new_eup[l]

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
    B[:,l]  = A[:,l]*c + A[:,k]*s
    B[k,k]  = A[k,k]*(c**2)- 2*A[k,l]*c*s + A[l,l]*(s**2)
    B[l,l]  = A[l,l]*(c**2)- 2*A[k,l]*c*s + A[k,k]*(s**2)
    B[k,l]  = 0 # (A[k,k] - A[l,l])*pl.cos(theta)*pl.sin(theta) \
                # + A[k,l]*((pl.cos(theta)**2) - (pl.sin(theta)**2))
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

    maxes   = pl.zeros(sz)
    rows    = pl.zeros(sz)
    cols    = pl.zeros(sz)
    eps     = A[0,1]

    # for i in pl.arange(0, sz-1):
    #     print "i=", i
    #     maxes   = pl.amax(abs(A[i,i+1:]))
    #     print "A_i_row=", A[i,i+1:]
    #     print "index of max A_i_row=", maxes  
    #     cols[i] = pl.argmax(abs(A[i,i+1:])) + 1
    #     print "A_i_row=", 
    #     print "index of max A_i_row=", cols

    # epsloop = pl.amax(maxes)
    # k       = pl.argmax(maxes)
    # l       = pl.int64(cols[k])
    
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
        print A
        print drac
        print k, l
        raw_input()
        A     = rotation_mat(A, k, l)
        k, l, epsloop = find_args(A)
        drac += 1

    print "Total number of iterations:", drac
    print A

run()




