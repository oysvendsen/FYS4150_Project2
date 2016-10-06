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
    for i in np.arange(len(data_splat) - 1):
        if i % 2 == 0:
            eigvals.append(np.float64(data_splat[i]))
        else:
            line_list = np.array(data_splat[i].split(', '))
            eigvecs.append(line_list.astype(np.float64))
    data_dict = {'lambda': np.array(eigvals), 'eigvecs': np.array(eigvecs)}

    return data_dict

def write2file(outstring,
               filename=os.getcwd()+"/pythondata/noninteracting.dat",
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

def test_rhon_n_dependency(rho_n, n):
    #run program with rho=0,..,rho_n with n steps for a non-interacting case.
    os.system("./build-Project2_QTcreator-Desktop_Qt_5_7_0_GCC_64bit-Release/Project2_QTcreator %d %f" % (n, rho_n))
    #fetch arrays written to "data/project2_noninteracting_rho0=0_rhoN=%d_n=%d.dat"%(rho_n,n)
    data = get_arrays(filename="data/project2_noninteracting_rho0=0_rhoN=%d_N=%d.dat"%(rho_n,n))

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
    print write2file("rho_N=%f, n=%f"%(rho_n, n))
    if any(within_tol):
        print write2file("all three eigenvalues WITHIN tolerance")
        for i in np.arange(3):
            print write2file("lambda_%d = %e"%(i,eigvals[i]))
    else:
        print write2file("all three eigenvalues OUTSIDE tolerance")
        for i in np.arange(3):
            print write2file("lambda_%d = %e"%(i,eigvals[i]))

#print get_arrays(filename="data/test.dat", length=3)

eig_tol = 1e-4
eig_known = np.array([3.0, 7.0, 11.0])
rho_N_range = np.array(range(1,15))
N_range = np.array([50, 100, 200, 300])
write2file("New data-file", append=False)

# ./main.cpp n rho_N
#run program "Project2_QTcreator" for these cases.
for rho_n in rho_N_range:
    test_rhon_n_dependency(rho_n, 200)
for n in N_range:
    test_rhon_n_dependency(4.0, n)

"""
# plots eigenvectors relevant to current omega
pl.figure()
rho = np.linspace(0, rho_n, n)
eigvecs = data['eigvecs']
for k in range(len(eigvecs)):
    pl.plot(rho, eigvecs[k], label=r'Non-interact. eigvec. $v_%d$' % k)

pl.title(r'Non-interactive eigenvectors for $\rho_N$ = %d, $N$=%d ' % (rho_n, n))
pl.xlabel(r'$\rho$')
pl.ylabel(r'$v$')
pl.legend(loc='best')
pl.savefig('non_interacting_eigvec_plot_rhoN=%d_N=%d' % (rho_n, n), dpi=300)
"""
