
"""
Reads in INTS/OVERLAP and SCF/MO1 from two bin files
"""

from re import A, search
import os
import sys
import subprocess
from urllib.parse import uses_fragment
import numpy as np
import h5py

def main() :

    if len(sys.argv) == 3:
      dataBin1, dataBin2 = sys.argv[1], sys.argv[2]
      threshold = 0.100000
    elif len(sys.argv) == 4:
      dataBin1, dataBin2, threshold = sys.argv[1], sys.argv[2], float(sys.argv[3])
    else:
      print("Too many or too few arguments")
      sys.exit(1)

    typeBin1, typeBin2 = determineMethod(dataBin1,dataBin2)

    mo1, mo2, overlap = prepareMO(dataBin1,dataBin2,typeBin1,typeBin2) 

    outputMOOverlap(mo1, mo2, overlap, threshold)



def determineMethod(dataBin1,dataBin2) :
    """
    takes in bin files and determines type 
    params:
        dataBin1/2 - .bin files
    returns:
        integer typeBin1 - REFTYPE of Bin 1 
        integer typeBin2 - REFTYPE of Bin 2 
    """
    typeBin1 = -1 
    typeBin2 = -1 

    bin1 = h5py.File(dataBin1)

    typeBin1 = bin1["/REF/REFTYPE"][0]

    print("RefType of Bin1 is ",typeBin1)

    bin2 = h5py.File(dataBin2)

    typeBin2 = bin2["/REF/REFTYPE"][0]

    print("RefType of Bin2 is ",typeBin2)

    if typeBin1 > typeBin2:
      print("larger ref type should come second")
      sys.exit(1)

    return typeBin1, typeBin2 


def prepareMO(dataBin1,dataBin2,typeBin1,typeBin2) :
    """
    Takes in bin file types, finds MOs, and then 
    prepares MOs to be same dimension 
    params:
        dataBin1/2 - bin files 
        typeBin1/2 - ref types 
    returns:
        np mo1 - mo entry of bin 1
        np mo2 - mo entry of bin 2
    """

    # assuming same basis set so grabbing overlap from bin1
    bin1 = h5py.File(dataBin1)

    mo1 = bin1["/SCF/MO1"][:]

    print("mo1 size",mo1.shape)

    #np.savetxt("mo1.txt ", mo1)

    # handling of overlap assumes 2c
    overlapAO = bin1["/INTS/OVERLAP"][:]

    print("overlap initial AO size ",overlapAO.shape)

    overlap = make_2c_block(overlapAO)

    print("overlap 2c AO size ",overlap.shape)

    #np.savetxt("overlapAO.txt", overlap)

    bin2 = h5py.File(dataBin2)

    mo2 = bin2["/SCF/MO1"][:]

    print("mo2 size",mo2.shape)

    if mo1.shape != mo2.shape:
      print("MOs should have same size. Please read the README for help.")
      sys.exit(1)

    #np.savetxt("mo2.txt", mo2)

    if typeBin1 == typeBin2:
      print("Ref type the same for both bin files, so no MO preparation is needed")
    else:
      print("Do not know how to convert between the two ref types. Please read README.")

    return mo1, mo2, overlap


def outputMOOverlap(mo1, mo2, overlap, threshold) :
    """
    takes in two MOs and overlap and calculate orbital projection
    params:
        mo1 - mo on bin file 1
        mo2 - mo on bin file 2
    """

    np.set_printoptions(precision=4)

    overlapMOcomp = mo1 @ overlap @ mo2.conj().transpose()
    overlapMO = np.abs(overlapMOcomp)

    i_range, j_range = overlapMO.shape

#    for i in range(i_range):
#      print(i, overlapMOcomp[:,i].conj().T @ overlapMOcomp[:,i])

    print("Printing out orbital contributions > ",threshold)

    for i in range(i_range) :
      print("**** Contributions from mo1 to orbital ",i+1," from mo2")
      for j in range(j_range) :
        if overlapMO[j,i]*overlapMO[j,i] > threshold :
          print(f"MO #{j+1:2d} : {overlapMO[j,i]*overlapMO[j,i]:6.4f}")
        
    #np.savetxt("overlapMO.txt", overlapMO)


# Written by Shiv Upadhyay
def make_2c_block(matrix):
    return np.block(
        [[matrix, np.zeros_like(matrix)], [np.zeros_like(matrix), matrix]]
    )

if __name__ == "__main__" :
    main()
