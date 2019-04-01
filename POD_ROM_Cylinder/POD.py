import numpy as np
from numpy import linalg as LA

def pod_mos(snapshot_matrix,data_tsteps):

    # Generate POD bases - method of snapshots
    new_mat = np.matmul(np.transpose(snapshot_matrix),snapshot_matrix)
    w,v = LA.eig(new_mat)

    # Bases calculation
    phi = np.real(np.matmul(snapshot_matrix,v))

    trange = np.arange(data_tsteps)
    phi[:,trange] = phi[:,trange]/np.sqrt(w[:])

    # Coefficients of the POD bases - Rows are mode numbers, columns are time
    coefficient_matrix = np.matmul(np.transpose(phi),snapshot_matrix)

    return phi, coefficient_matrix

def reconstruct(phi,coefficient_matrix,mean,xdof,ydof):
    rec_snap = np.matmul(phi,coefficient_matrix) + mean
    return np.reshape(rec_snap,newshape=(xdof,ydof))