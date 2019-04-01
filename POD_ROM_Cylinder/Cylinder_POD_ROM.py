import numpy as np
import matplotlib.pyplot as plt
from numpy import linalg as LA
from scipy.io import loadmat
import matplotlib.animation as animation
from mpl_toolkits.axes_grid1 import make_axes_locatable
import tensorflow as tf

from FC_NETS import resnet_for_dynamics, evaluate_rom_deployment
from LSTM_NETS import lstm_for_dynamics, evaluate_rom_deployment_lstm
from POD import pod_mos, reconstruct

np.random.seed(10)
tf.random.set_random_seed(10)

tsteps = int(input('Number of time steps to train with?\n'))
num_trunc = int(input('Number of modes to truncate to?\n'))

#--------------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------------
# Data reading section
#--------------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------------

# Loading flow around cylinder data
x = loadmat('VonKarmanStreet.mat') # This is a dictionary
vel = x['U'] # Converting to velocity

# print(x['DESCRIPTION'])
components, xdof, ydof, data_tsteps = np.shape(vel)[0], np.shape(vel)[1], np.shape(vel)[2], np.shape(vel)[3]

vis = input('Visualize velocities and POD basis? [0] - No, [1] - Yes\n')

if int(vis) == 1:
    fig,ax = plt.subplots(nrows=2)
    def animate(i):
           ax[0].clear()
           ax[0].imshow(vel[0,:,:,i])
           ax[0].set_title('U Snapshot no.'+'%03d'%(i)) 

           ax[1].clear()
           ax[1].imshow(vel[1,:,:,i])
           ax[1].set_title('V Snapshot no.'+'%03d'%(i)) 

    interval = 0.1#in seconds     
    ani = animation.FuncAnimation(fig,animate,interval=interval*1e+3,blit=False)
    plt.show()

#--------------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------------
# POD Method of Snapshots section
#--------------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------------

# Make snapshot matrix - remove mean - dimension 2 is time
snapshot_matrix = np.zeros((components,xdof*ydof,data_tsteps))
for t in range(data_tsteps):
    snapshot_matrix[0,:,t] = vel[0,:,:,t].flatten()
    snapshot_matrix[1,:,t] = vel[1,:,:,t].flatten()

snapshot_matrix_mean = np.mean(snapshot_matrix,axis=2)

snapshot_matrix[0,:,:] = snapshot_matrix[0,:,:]-snapshot_matrix_mean[0,:,None]
snapshot_matrix[1,:,:] = snapshot_matrix[1,:,:]-snapshot_matrix_mean[1,:,None]

phi_u, coefficient_matrix_u = pod_mos(snapshot_matrix[0,:,:],data_tsteps) # Coefficient matrix dimension 1 is time
phi_v, coefficient_matrix_v = pod_mos(snapshot_matrix[1,:,:],data_tsteps)

if int(vis) == 1:
    fig,ax = plt.subplots(nrows=2)
    def animate(i):
           ax[0].clear()
           ax[0].imshow(np.reshape(phi_u[:,i],newshape=(xdof,ydof))) # POD matrix dimension 1 is basis
           ax[0].set_title('U POD Basis no.'+'%03d'%(i)) 

           ax[1].clear()
           ax[1].imshow(np.reshape(phi_v[:,i],newshape=(xdof,ydof)))
           ax[1].set_title('V POD Basis no.'+'%03d'%(i)) 

    interval = 0.5#in seconds     
    ani = animation.FuncAnimation(fig,animate,interval=interval*1e+3,blit=False)
    plt.show()

# Setup truncated coefficient matrix - reduced DOF evolution in POD space
coeff_trunc = np.concatenate((coefficient_matrix_u[0:num_trunc,:],coefficient_matrix_v[0:num_trunc,:]),axis=0)
coeff_trunc_train = coeff_trunc[:,:tsteps]

# Need to scale the rows individually to ensure common scale
cft_max = np.amax(coeff_trunc_train,axis=1)
cft_min = np.amin(coeff_trunc_train,axis=1)

coeff_trunc_scaled = -1.0 + 2.0*(coeff_trunc_train[:,:]-cft_min[:,None])/(cft_max[:,None]-cft_min[:,None])

#--------------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------------
# Learning dynamics using ML section
#--------------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------------
# See the performance for different ML frameworks
valid_model = lstm_for_dynamics(coeff_trunc_scaled)
output_state_ml, output_tracker_ml = evaluate_rom_deployment_lstm(valid_model,coeff_trunc_scaled,num_trunc,data_tsteps)

output_state_ml = np.transpose(output_state_ml)
output_tracker_ml = np.transpose(output_tracker_ml)

# Rescale ML prediction
output_state_ml[:,0] = (1.0+output_state_ml[:,0])*(cft_max[:]-cft_min[:])/2.0+cft_min[:]
output_tracker_ml[:,:] = (1.0+output_tracker_ml[:,:])*(cft_max[:,None]-cft_min[:,None])/2.0+cft_min[:,None]

#Find final state values and compare
u_true = reconstruct(phi_u[:,0:num_trunc],coeff_trunc[0:num_trunc,-1],snapshot_matrix_mean[0,:],xdof,ydof)
u_ml = reconstruct(phi_u[:,0:num_trunc],output_state_ml[0:num_trunc,-1],snapshot_matrix_mean[0,:],xdof,ydof)

#--------------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------------
# Visualization and performance assessment section
#--------------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------------

# Visualization of performance in physical domain
fig, ax = plt.subplots(nrows=3,ncols=1)

im1 = ax[0].imshow(u_true)
ax[0].set_title('True')
divider = make_axes_locatable(ax[0])
cax1 = divider.append_axes('right', size='5%', pad=0.4)
fig.colorbar(im1, cax=cax1, orientation='vertical')

im2 = ax[1].imshow(u_ml)
ax[1].set_title('ML')
divider = make_axes_locatable(ax[1])
cax2 = divider.append_axes('right', size='5%', pad=0.4)
fig.colorbar(im2, cax=cax2, orientation='vertical')


im3 = ax[2].imshow(np.abs(u_true-u_ml))
ax[2].set_title('L1 Difference')
divider = make_axes_locatable(ax[2])
cax3 = divider.append_axes('right', size='5%', pad=0.4)
fig.colorbar(im3, cax=cax3, orientation='vertical')

plt.show()

# Visualization of modal evolution
fig, ax = plt.subplots(nrows=3,ncols=1)
ax[0].plot(coeff_trunc[0,:],label='True')
ax[0].plot(output_tracker_ml[0,:],label='ML')
ax[0].set_title('Mode 1')

ax[1].plot(coeff_trunc[1,:],label='True')
ax[1].plot(output_tracker_ml[1,:],label='ML')
ax[1].set_title('Mode 2')

ax[2].plot(coeff_trunc[2,:],label='True')
ax[2].plot(output_tracker_ml[2,:],label='ML')
ax[2].set_title('Mode 3')

plt.legend()
plt.show()


