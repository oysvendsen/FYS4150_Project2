import numpy as np
import matplotlib.pyplot as pl
import sys
import os
from matplotlib import rc
rc('font',**{'family':'serif'})

"""
Find Eigenvalues of Schrodinger equation for 
three values of rho_N, and several values of N.
When are they sufficently close to known solution?
"""

def get_arrays(filename):
    #fetch data from file
    with open(filename, 'r') as infile:
        data = infile.read()
    #rewrite data to nested list
    
    data_splat = data.split('\n')

    eigvals = []
    eigvecs = []
    #sort out columns as arrays in dictionary
    for line in np.range(len(data_splat)):
        if line % 2 == 0:
            eigvals.append(np.float64(data_splat[line][0]))
        else:
            eigvecs.append(pl.array(data_splat[line].split(', ')).astype(np.float64))
    data_dict = {'lambda': np.array(eigvals), 'eigvecs': np.array(eigvecs)}

    return data_dict

def write2file(outstring,
               filename=curdir+"/data/interacting.dat",
               append=True):
    """
    If 'append' is True:
    -open 'filename'.
    -append 'outstring' to end of file
    -close file
    If 'append' is False:
    -open a new file 'filename' (deleting the old one)
    -write 'outstring' to file
    -close file
    """
    outstring = str(outstring)
    if append:
        with open(filename,'a') as outfile:
            outfile.write(outstring + "\n")
    else:
        with open(filename,'w') as outfile:
            outfile.write(outstring + "\n")
    return outstring


#print get_arrays(filename="data/test.dat", length=3)

eig_tol     = 1e-4
eig_known   = np.array([3.0, 7.0, 11.0])
rho_n       = 5. # tbd
n           = 500 # tbd
write2file("New data-file", append=False)

rho0        = 0.
omega       = np.array([0.01, 0.5, 1., 5.]).

    # still need the non-interacting eigvectors, made in a previous run
old_data = get_arrays(filename="data/project2_noninteracting_rho0=0_rhoN=%d_N=%d_omega=%d.dat" % ( rho_n, n))
non_eigvecs = old_data['eigvecs'] # non-interacting eigvecs

# ./main.cpp n rho_N rho_N omega 'interact'
#run program "Project2_QTcreator" for these cases.
for omega in omegas:
    #run program with rho=0,..,rho_n with n steps for a non-interacting case.
    os.system("./build-Project2_QTcreator-Desktop_Qt_5_7_0_GCC_64bit-Release/Project2_QTcreator %d %f %f %f %s" % (n, rho0, rho_n, omega, "interact"))
    #fetch arrays written to "data/project2_noninteracting_rho0=0_rhoN=%d_n=%d.dat"%(rho_n,n)
    data = get_arrays( filename="data/project2_interacting_rho0=0_rhoN=%d_N=%d_omega=%d.dat" % (rho_n, n, (omega/float(1000)) ) )

    #is eigenvalues withing tolerance?
    eigvals = data['lambda']
    within_tol = []
    
    #test the value of the three lowest eigenvalues
    for i in np.arange(3):
        if abs(eigvals[i] - eig_known[i]) < eig_tol:
            within_tol.append(True)
        else:
            within_tol.append(False)

    #store result in data-file.
    print write2file("rho_0=%f, rho_N=%f, n=%f"%(rho_0, rho_n, n))
    if within_tol:
        print write2file("all three eigenvalues WITHIN tolerance")
        for i in np.arange(3):
            print write2file("lambda_%d = %e"%(i,eigvals[i]))
    else:
        print write2file("all three eigenvalues OUTSIDE tolerance")
        for i in np.arange(3):
            print write2file("lambda_%d = %e"%(i,eigvals[i]))

    # plots eigenvectors relevant to current omega
    pl.figure()
    rho = np.linspace(rho0, rho_n, n)
    eigvecs = data['eigvecs']
    for k in range(len(eigvecs)):
        pl.plot(rho, non_eigvecs[k], label=r'Non-interact. eigvec. $v_%d$' % k)
        pl.plot(rho, eigvecs[k], label=r'Interact. eigvec. $v_%d$' % k)

    pl.title(r'Non-interactive- vs Interactive Eigenvectors for $\omega_r$ = %.2f' % omega)
    pl.xlabel(r'$\rho$')
    pl.ylabel(r'$v$')
    pl.legend(loc='best')
    pl.savefig('interacting_eig_plot_omega=%.2f' % omega, dpi=300)

