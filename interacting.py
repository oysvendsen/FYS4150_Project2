import numpy as np
import matplotlib.pyplot as pl
import sys
import os

"""
Find Eigenvalues of Schrodinger equation for 
three values of rho_N, and several values of N.
When are they sufficently close to known solution?
"""

def get_arrays(filename, length):
    #fetch data from file
    with open(filename, 'r') as infile:
        data = infile.read()
    #rewrite data to nested list
    
    data_splat = data.split('\n')

    eigvals = []
    eigvecs = []
    #sort out columns as arrays in dictionary
    for line in range(len(data_splat)):
        if line % 2 == 0:
            eigvals.append(np.float64(data_splat[line][0]))
        else:
            eigvecs.append(pl.array(data_splat[line].split(', ')).astype(np.float64))
    data_dict = {'lambdas': np.array(eigvals), 'eigvecs': np.array(eigvecs)}

    return data_dict

def write2file(outstring,
               filename=curdir+"/data/noninteracting.dat",
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

eig_tol = 1e-4
eig_known = np.array([3.0, 7.0, 11.0])
rho_N_range = np.array([1.0, 5.0, 10.0])
N_range = np.array([50, 100, 500, 1000])
write2file("New data-file", append=False)

# ./main.cpp n rho_N
#run program "Project2_QTcreator" for these cases.
for rho_n in rho_N_range:
    for n in N_range:
        #run program with rho=0,..,rho_n with n steps for a non-interacting case.
        os.system("./build-Project2_QTcreator-Desktop_Qt_5_7_0_GCC_64bit-Release/Project2_QTcreator")
        #fetch arrays written to "data/project2_noninteracting_rho0=0_rhoN=%d_n=%d.dat"%(rho_n,n)
        data = get_arrays(filename="data/project2_noninteracting_rho0=0_rhoN=%d_n=%d.dat"%(rho_n,n), length=n+1)
        #sort the lambda's in increasing order
        eig_sorted = []
        for i in range(len(eigenvalues)):
            lowest = eigenvalues[0]
            for j in range(len(eigenvalues)):
                if eigenvalues[j] <= lowest:
                    lowest = eigenvalues[j]
            eig_sorted[i] = lowest

        #is eigenvalues withing tolerance?
        eigenvalues = data['lambda']
        within_tol = []
        
        #test the value of the three lowest eigenvalues
        for i in range(3):
            if abs(eig_sorted[i] - eig_known[i]) < eig_tol:
                within_tol.append(True)
            else:
                within_tol.append(False)

        #store result in data-file.
        print write2file("rho_0=%f, rho_N=%f, n=%f"%(rho_0, rho_n, n))
        if within_tol:
            print write2file("all three eigenvalues WITHIN tolerance")
            for i in range(3):
                print write2file("lambda_%d = %e"%(i,eig_sorted[i]))
        else:
            print write2file("all three eigenvalues OUTSIDE tolerance")
            for i in range(3):
                print write2file("lambda_%d = %e"%(i,eig_sorted[i]))
