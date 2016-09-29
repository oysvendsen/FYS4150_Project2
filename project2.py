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
    tau             = (A[l,l] - A[k,k])/float(2*A[k,l])
    t               = -tau + (1 + tau**2)**0.5
    c               = 1./(1+t**2)
    s               = c*t

    B = A
    B[:,k]  = A[:,k]*c - A[:,l]*s
    B[:,l]  = A[:,l]*c + A[:,k]*s
    B[k,k]  = A[k,k]*(c**2)- 2*A[k,l]*c*s + A[l,l]*(s**2)
    B[l,l]  = A[l,l]*(c**2)- 2*A[k,l]*c*s + A[k,k]*(s**2)
    B[k,l]  = 0 #(A[k,k] - A[l,l])*pl.cos(theta)*pl.sin(theta) \
                # + A[k,l]*((pl.cos(theta)**2) - (pl.sin(theta)**2))
    B[l,k]  = 0

    return B









def run():

    for i in range(len(diag)):



l = 0
n = int(1e3) + 1 # or something

rho_0   = 0.
rho_end = 1000.
rho     = pl.linspace(rho_0, rho_end, n)
h       = (rho_end - rho_0)/float(n)

u       = pl.zeros(n)
eup     = -1./(h**2)*pl.ones(len(rho))
edown   = -1./(h**2)*pl.ones(len(rho))
Vee     = rho**2
diag    = 2./(h**2) + Vee




