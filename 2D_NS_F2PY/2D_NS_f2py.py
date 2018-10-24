import numpy as np
import matplotlib.pyplot as plt
import time
from numpy import f2py
from keras.models import Sequential, Model, load_model
from keras.layers import Dense, dot, Input
import keras.backend as K

import Fortran_functions
import Multigrid_solver
import Gauss_Siedel
import Spectral_Poisson
import Approximate_Deconvolution
import Fortran_Gaussian_Filter
import Fortran_Pade_Filter
import Ml_convolution
import Fortran_Dynamic

#-------------------------------------------------------------------------------------#
#-------------------------------------------------------------------------------------#
#Standard functions for solver and inputs
#-------------------------------------------------------------------------------------#
#-------------------------------------------------------------------------------------#
def init_domain():

    global nx, ny, kappa, Re_n, lx, ly, lt, nt, dt, dx, dy, gs_tol
    global problem
    global sigma#For AD-LES
    global closure_choice
    global n_samples #For ML-AD-LES framework
    
    nx = 256
    ny = 256
    kappa = 2.0
    Re_n = 32000.0
    kappa = 2.0
    lx = 2.0 * np.pi
    ly = 2.0 * np.pi
    dx = lx / float(nx)
    dy = ly / float(ny)#Not nx-1 since last point is to be assumed within the domain
    sigma = 1.0

    gs_tol = 1.0e-2#Keep low for GS, keep high for MG

    n_samples = 8000#For ML - AD - LES classification

    problem = 'HIT'
    closure_choice = 8#0 - uDNS, 1- AD LES, 2 - ML AD LES, 3 - RF LES, 4 - ML SGS, 5 - Smag SGS, 6 - Leith SGS, 7 - Nearest neighbor LES, 8 - ML_Convolution_AD, 9 - Feature SGS, 10 - Dyn Smag, 11 - Dyn Leith

    if problem == 'TGV':
        lt = 0.1
    else:
        lt = 4.0

    dt = 1.0e-3
    nt = int(lt/dt)

    omega = np.zeros(shape=(nx,ny),dtype='double', order='F')
    psi = np.zeros(shape=(nx, ny), dtype='double', order='F')

    return omega, psi

def initialize_ic_bc(omega, psi):
    if problem == 'TGV':
        for i in range(nx):
            for j in range(ny):
                x = float(i) * dx
                y = float(j) * dy
                omega[i, j] = 2.0 * kappa * np.cos(kappa * x) * np.cos(kappa * y)
                psi[i, j] = 1.0 / kappa * np.cos(kappa * x) * np.cos(kappa * y)

    elif problem == 'HIT':
        Fortran_functions.hit_init_cond(omega, dx, dy)

        Spectral_Poisson.solve_poisson(psi, -omega, dx, dy)
        # Gauss_Siedel.solve_poisson(psi,-omega, dx, dy, gs_tol)
        # Multigrid_solver.solve_poisson_periodic(psi,-omega, dx, dy, gs_tol)

    elif problem == 'DSL':

        delta = 0.05
        sig = 15.0 / np.pi

        for i in range(nx):
            for j in range(ny):
                x = float(i) * dx
                y = float(j) * dy

                if y <= np.pi:
                    omega[i, j] = delta * np.cos(x) - sig / (np.cosh(sig * (y - np.pi / 2.0)) ** 2)
                    # psi[i,j] = delta*np.cos(x) - np.log(np.cosh(np.pi*sig/2.0 - sig*y))/sig
                else:
                    omega[i, j] = delta * np.cos(x) + sig / (np.cosh(sig * (3.0 * np.pi / 2.0 - y)) ** 2)
                    # psi[i,j] = delta*np.cos(x) - np.log(np.cosh(3.0*np.pi*sig/2.0 - sig*y))/sig

        Spectral_Poisson.solve_poisson(psi, -omega, dx, dy)

def post_process(omega,psi):

    if problem == 'TGV':

        fig, ax = plt.subplots(nrows=3,ncols=1)

        levels = np.linspace(-4,4,10)

        ax[0].set_title("Numerical Solution - TGV")
        plot1 = ax[0].contourf(omega[:, :],levels=levels)
        plt.colorbar(plot1, format="%.2f",ax=ax[0])

        omega_true = np.zeros(shape=(nx,ny),dtype='double')

        for i in range(nx):
            for j in range(ny):

                x = float(i) * dx
                y = float(j) * dy

                omega_true[i, j] = 2.0 * kappa * np.cos(kappa * x) * np.cos(kappa * y)*np.exp(-2.0*kappa*kappa*lt/Re_n)

        ax[1].set_title("Exact Solution")
        plot2 = ax[1].contourf(omega_true[:, :],levels=levels)
        plt.colorbar(plot2, format="%.2f",ax=ax[1])


        #levels = np.linspace(0.0,1.0e-3,10)

        ax[2].set_title("L1 - Error")
        plot3 = ax[2].contourf(np.abs(omega_true[:, :]-omega[:,:]))
        plt.colorbar(plot3,ax=ax[2])

        plt.show()


    elif problem == 'HIT':

        fig, ax = plt.subplots(nrows=1,ncols=1)

        #levels = np.linspace(-50,50,10)

        ax.set_title("Numerical Solution - HIT")
        plot1 = ax.contourf(omega[:, :])#,levels=levels)
        plt.colorbar(plot1, format="%.2f")

        plt.show()

        arr_len = int(0.5*np.sqrt(float(nx*nx + ny*ny)))-1
        eplot = np.zeros(arr_len+1,dtype='double')
        kplot = np.arange(0,arr_len+1,1,dtype='double')

        Fortran_functions.spec(omega,eplot)

        scale_plot = np.array([[10,0.1],[100,1.0e-4]])

        plt.loglog(kplot,eplot)
        plt.loglog(scale_plot[:,0],scale_plot[:,1])
        plt.xlim([1,1.0e3])
        plt.ylim([1e-8,1])
        plt.xlabel('k')
        plt.ylabel('E(k)')
        plt.title('Angle averaged energy spectra')
        plt.show()

        np.save('Field.npy',omega)
        np.save('Spectra.npy',[kplot, eplot])

    elif problem == 'DSL':
        fig, ax = plt.subplots(nrows=2,ncols=1)

        #levels = np.linspace(-50,50,10)

        ax[0].set_title("Numerical Solution - DSL - Vorticity")
        plot1 = ax[0].contourf(omega[:, :])#,levels=levels)
        plt.colorbar(plot1, format="%.2f",ax=ax[0])

        ax[1].set_title("Numerical Solution - DSL - Streamfunction")
        plot2 = ax[1].contourf(psi[:, :])#,levels=levels)
        plt.colorbar(plot2, format="%.2f",ax=ax[1])

        plt.show()

#-------------------------------------------------------------------------------------#
#-------------------------------------------------------------------------------------#
#Machine Learning model function deployment
#-------------------------------------------------------------------------------------#
#-------------------------------------------------------------------------------------#
def load_pretrained_model():

    if closure_choice == 2:
        # Initialization - just load model directly from hd5
        model = load_model('ml_ad_model.hd5')
        ml_model = K.function([model.layers[0].input],
                                [model.layers[-1].output])  # One layer input one classification output
    elif closure_choice == 4:

        def coeff_determination(y_true, y_pred):
            SS_res = K.sum(K.square(y_true - y_pred))
            SS_tot = K.sum(K.square(y_true - K.mean(y_true)))
            return (1 - SS_res / (SS_tot + K.epsilon()))

        # Initialization - just load model directly from hd5
        model = load_model('ml_sgs_model.hd5',custom_objects={'coeff_determination': coeff_determination})
        ml_model = K.function([model.layers[0].input],
                                [model.layers[-1].output])  # One layer input, one regression output
    elif closure_choice == 7:

        def coeff_determination(y_true, y_pred):
            SS_res = K.sum(K.square(y_true - y_pred))
            SS_tot = K.sum(K.square(y_true - K.mean(y_true)))
            return (1 - SS_res / (SS_tot + K.epsilon()))

        # Initialization - just load model directly from hd5
        model = load_model('ml_sgs_model.hd5',custom_objects={'coeff_determination': coeff_determination})
        ml_model = K.function([model.layers[0].input],
                                [model.layers[-1].output])  # One layer input, one regression output

    elif closure_choice == 8:

        def coeff_determination(y_true, y_pred):
            SS_res = K.sum(K.square(y_true - y_pred))
            SS_tot = K.sum(K.square(y_true - K.mean(y_true)))
            return (1 - SS_res / (SS_tot + K.epsilon()))

        # Initialization - just load model directly from hd5
        model = load_model('Forward_Filter.hd5',custom_objects={'coeff_determination': coeff_determination})
        ml_model_forward = K.function([model.layers[0].input],
                                [model.layers[-1].output])  # One layer input, one regression output

        # Initialization - just load model directly from hd5
        model = load_model('Inverse_Filter.hd5', custom_objects={'coeff_determination': coeff_determination})
        ml_model_inverse = K.function([model.layers[0].input],
                                      [model.layers[-1].output])  # One layer input, one regression output

        return ml_model_forward, ml_model_inverse

    elif closure_choice == 9:

        def coeff_determination(y_true, y_pred):
            SS_res = K.sum(K.square(y_true - y_pred))
            SS_tot = K.sum(K.square(y_true - K.mean(y_true)))
            return (1 - SS_res / (SS_tot + K.epsilon()))

        # Initialization - just load model directly from hd5
        model = load_model('Feature_SGS_Model.hd5',custom_objects={'coeff_determination': coeff_determination})
        ml_model = K.function([model.layers[0].input],
                                [model.layers[-1].output])  # One layer input, one regression output

    return ml_model#Return the precompiled function

def deploy_model_classification(omega,ml_model):

    #Sampling the field randomly for omega stencils
    sampling_matrix = Fortran_functions.field_sampler(omega, n_samples)
    #Classify
    return_matrix = ml_model([sampling_matrix])[0]#Using precompiled ML network
    #Summing up different predictions for filter strength
    return_matrix = np.sum(return_matrix,axis=0)
    #Finding the filter strength with max predictions
    bucket = np.argmax(return_matrix,axis=0)

    #Modifying the global sigma
    global sigma

    #Note inverse relationship with training
    if bucket == 0:
        sigma = 2.0
    elif bucket == 1:
        sigma = 1.5
    elif bucket == 2:
        sigma = 1.0

def deploy_model_regression(omega,psi,ml_model):

    # Validated - samples the entire field for omega, psi stencil and 2 pointwise invariants (20 inputs per point)
    sampling_matrix = Fortran_functions.sgs_field_sampler(omega, psi, dx, dy)

    # Prediction from precompiled keras function - validated
    return_matrix = ml_model([sampling_matrix])[0]

    # Laplacian calculator in Fortran
    laplacian = Fortran_functions.laplacian_calculator(omega,dx,dy)

    # Placeholders incase we need normalization - superfluous now
    source_max = 100.0
    source_min = -100.0

    # Defining the sgs array
    sgs = np.zeros(shape=(nx,ny),dtype='double', order='F')

    # Reshaping keras prediction in sgs array - validated using visualization
    # laplacian used for postprocessing SGS output for stability
    Fortran_functions.sgs_reshape_non_sdev(return_matrix,sgs,laplacian,source_max,source_min)

    return sgs

def deploy_model_nearest_neighbor(omega,psi,ml_model):

    # Validated - samples the entire field for omega, psi stencil and 2 pointwise invariants (20 inputs per point)
    sampling_matrix = Fortran_functions.sgs_field_sampler(omega, psi, dx, dy)

    # Prediction from precompiled keras function - validated
    return_matrix = ml_model([sampling_matrix])[0]

    # Laplacian calculator in Fortran
    laplacian = Fortran_functions.laplacian_calculator(omega,dx,dy)

    # Placeholders incase we need normalization - superfluous now
    source_max = 100.0
    source_min = -100.0

    # Defining the sgs array
    sgs_ml = np.zeros(shape=(nx,ny),dtype='double', order='F')

    # Reshaping keras prediction in sgs array - validated using visualization
    # laplacian used for postprocessing SGS output for stability
    Fortran_functions.sgs_reshape_non_sdev(return_matrix,sgs_ml,laplacian,source_max,source_min)

    # Smagorinsky model estimate
    sgs_smag = Fortran_functions.smag_source_term(omega, psi, dx, dy)

    # Leith model estimate
    sgs_leith = Fortran_functions.leith_source_term(omega, psi, dx, dy)

    # AD estimate
    sfs_ad = Approximate_Deconvolution.approximate_deconvolution(omega,psi,dx,dy,sigma)

    # Final sgs estimate according to nearest neighbor classification
    sgs = np.zeros(shape=(nx, ny), dtype='double', order='F')
    Fortran_functions.sgs_classify(sgs,sgs_ml,sgs_smag,sgs_leith,sfs_ad)

    return sgs

def deploy_model_convolution(omega,psi,ml_model_forward,ml_model_inverse):

    #Number of iterative resubstitutions
    n_ad_iter = 3

    #Recording some parameters needed for normalization
    omega_mean = Ml_convolution.mean_calc(omega)
    psi_mean = Ml_convolution.mean_calc(psi)

    omega_std = Ml_convolution.std_calc(omega)
    psi_std = Ml_convolution.std_calc(psi)

    wtemp = np.copy(omega)
    stemp = np.copy(psi)

    #Calculating jacobian using grid-resolved variables
    jc_fc = Ml_convolution.jacobian_calc(omega, psi, dx, dy)

    jc_fc_mean = Ml_convolution.mean_calc(jc_fc)
    jc_fc_std = Ml_convolution.std_calc(jc_fc)

    Ml_convolution.normalize_field(wtemp,omega_mean,omega_std)
    Ml_convolution.normalize_field(stemp,psi_mean,psi_std)

    # Using inverse filter to find deconvolved fields
    omega_sampler = Ml_convolution.conv_field_sampler(wtemp)
    Ml_convolution.field_reshape(ml_model_inverse([omega_sampler])[0], wtemp)
    psi_sampler = Ml_convolution.conv_field_sampler(stemp)
    Ml_convolution.field_reshape(ml_model_inverse([psi_sampler])[0], stemp)


    Ml_convolution.denormalize_field(wtemp,omega_mean,omega_std)
    Ml_convolution.denormalize_field(stemp,psi_mean,psi_std)

    #Calculating Jacobian of deconvolved variables
    jc_ad = Ml_convolution.jacobian_calc(wtemp, stemp, dx, dy)
    Ml_convolution.normalize_field(jc_ad,jc_fc_mean,jc_fc_std)

    del (wtemp,stemp)

    #Convolving the Jacobian
    jc_ad_sampler = Ml_convolution.conv_field_sampler(jc_ad)
    jc_ad_f = np.copy(omega)
    Ml_convolution.field_reshape(ml_model_forward([jc_ad_sampler])[0], jc_ad_f)

    Ml_convolution.denormalize_field(jc_ad_f, jc_fc_mean, jc_fc_std)

    # Laplacian calculator in Fortran
    laplacian = Ml_convolution.laplacian_calculator(omega,dx,dy)
    # sfs = Ml_convolution.ml_sfs_calculator_weiss(omega,psi,jc_fc, jc_ad_f, laplacian, dx, dy)
    #sfs = Ml_convolution.ml_sfs_calculator_bs_1(jc_fc, jc_ad_f, laplacian)
    # sfs = Ml_convolution.ml_sfs_calculator_bs_2(jc_fc, jc_ad_f, laplacian)
    # sfs = Ml_convolution.ml_sfs_calculator_bs_3(jc_fc, jc_ad_f, laplacian)
    # sfs = Ml_convolution.ml_sfs_calculator_bs_4(jc_fc, jc_ad_f, laplacian)
    # sfs = Ml_convolution.ml_sfs_calculator_bs_5(jc_fc, jc_ad_f, laplacian)
    sfs = Ml_convolution.ml_sfs_calculator_bs_6(jc_fc, jc_ad_f, laplacian)

    return sfs

def deploy_model_feature_regression(omega,psi,ml_model):
    # Samples the entire field (pointwise) for features
    sampling_matrix = Fortran_functions.field_feature_sampler(omega, psi, dx, dy)
    # Prediction from precompiled keras function - validated
    return_matrix = ml_model([sampling_matrix])[0]

    # Laplacian calculator in Fortran
    laplacian = Fortran_functions.laplacian_calculator(omega, dx, dy)

    # Placeholders incase we need normalization - superfluous now
    source_max = 100.0
    source_min = -100.0

    # Defining the sgs array
    sgs = np.zeros(shape=(nx, ny), dtype='double', order='F')

    # Reshaping keras prediction in sgs array - validated using visualization
    # laplacian used for postprocessing SGS output for stability
    Fortran_functions.sgs_reshape_non_sdev(return_matrix, sgs, laplacian, source_max, source_min)

    return sgs

#-------------------------------------------------------------------------------------#
#-------------------------------------------------------------------------------------#
#Standard eddy-viscosity closures/ AD
#-------------------------------------------------------------------------------------#
#-------------------------------------------------------------------------------------#
def smag_sgs_calc(omega,psi):

    sgs = Fortran_functions.smag_source_term_div(omega,psi,dx,dy)
    return sgs

def leith_sgs_calc(omega,psi):

    sgs = Fortran_functions.leith_source_term_div(omega,psi,dx,dy)
    return sgs

def ad_sfs_calc(omega,psi):
    sfs = Approximate_Deconvolution.approximate_deconvolution(omega,psi,dx,dy,sigma)
    return sfs

#-------------------------------------------------------------------------------------#
#-------------------------------------------------------------------------------------#
#Dynamic eddy-viscosity closures
#-------------------------------------------------------------------------------------#
#-------------------------------------------------------------------------------------#
def dynamic_smag_calc(omega,psi):

    laplacian = Fortran_functions.laplacian_calculator(omega, dx, dy)
    sgs = Fortran_Dynamic.dynamic_smagorinsky(omega,psi,laplacian,dx,dy)

    return sgs

def dynamic_leith_calc(omega,psi):

    laplacian = Fortran_functions.laplacian_calculator(omega, dx, dy)
    sgs = Fortran_Dynamic.dynamic_leith(omega,psi,laplacian,dx,dy)

    return sgs

#-------------------------------------------------------------------------------------#
#-------------------------------------------------------------------------------------#
#Runge Kutta third-order for different closures
#-------------------------------------------------------------------------------------#
#-------------------------------------------------------------------------------------#
def tvdrk3_fortran_rf_les(omega,psi):

    Fortran_Pade_Filter.filter_pade(omega,sigma)

    oneth = 1.0 / 3.0
    twoth = 2.0 / 3.0

    # Step 1
    # Calculate RHS
    f = Fortran_functions.rhs_periodic(psi, omega, dx, dy, Re_n)
    omega_1 = omega + dt * (f)

    # Fortran update for Poisson Equation
    # Multigrid_solver.solve_poisson_periodic(psi, -omega_1, dx, dy, gs_tol)
    Spectral_Poisson.solve_poisson(psi,-omega_1, dx, dy)

    # Step 2
    # Calculate RHS
    f = Fortran_functions.rhs_periodic(psi, omega_1, dx, dy, Re_n)
    omega_2 = 0.75 * omega + 0.25 * omega_1 + 0.25 * dt * (f)

    # Fortran update for Poisson Equation
    Spectral_Poisson.solve_poisson(psi,-omega_2, dx, dy)

    # Step 2
    # Calculate RHS
    f = Fortran_functions.rhs_periodic(psi, omega_2, dx, dy, Re_n)

    for i in range(nx):
        for j in range(ny):
            omega[i, j] = oneth * omega[i, j] + twoth * omega_2[i, j] + twoth * dt * (f[i, j])


    # Fortran update for Poisson Equation
    Spectral_Poisson.solve_poisson(psi,-omega, dx, dy)

def tvdrk3_fortran_ml_ad_les(omega,psi,ml_model):

    #Need to add ML based sigma estimation - modifies global
    deploy_model_classification(omega,ml_model)

    oneth = 1.0 / 3.0
    twoth = 2.0 / 3.0

    # Step 1
    # Calculate RHS
    f = Fortran_functions.rhs_periodic(psi, omega, dx, dy, Re_n)
    sfs = Approximate_Deconvolution.approximate_deconvolution(omega, psi, dx, dy, sigma)
    omega_1 = omega + dt * (f + sfs)

    # Fortran update for Poisson Equation
    # Multigrid_solver.solve_poisson_periodic(psi, -omega_1, dx, dy, gs_tol)
    Spectral_Poisson.solve_poisson(psi,-omega_1, dx, dy)

    # Step 2
    # Calculate RHS
    f = Fortran_functions.rhs_periodic(psi, omega_1, dx, dy, Re_n)
    sfs = Approximate_Deconvolution.approximate_deconvolution(omega_1, psi, dx, dy, sigma)
    omega_2 = 0.75 * omega + 0.25 * omega_1 + 0.25 * dt * (f + sfs)

    # Fortran update for Poisson Equation
    Spectral_Poisson.solve_poisson(psi,-omega_2, dx, dy)

    # Step 2
    # Calculate RHS
    f = Fortran_functions.rhs_periodic(psi, omega_2, dx, dy, Re_n)
    sfs = Approximate_Deconvolution.approximate_deconvolution(omega_2, psi, dx, dy, sigma)

    for i in range(nx):
        for j in range(ny):
            omega[i, j] = oneth * omega[i, j] + twoth * omega_2[i, j] + twoth * dt * (f[i, j] + sfs[i, j])

    # Fortran update for Poisson Equation
    Spectral_Poisson.solve_poisson(psi,-omega, dx, dy)

def tvdrk3_fortran_ml_sgs(omega,psi,ml_model):

    oneth = 1.0 / 3.0
    twoth = 2.0 / 3.0

    # Step 1
    # Calculate RHS
    f = Fortran_functions.rhs_periodic(psi, omega, dx, dy, Re_n)
    #Need to add ML based sgs computation
    sgs = deploy_model_regression(omega,psi,ml_model)

    omega_1 = omega + dt * (f + sgs)

    # Fortran update for Poisson Equation
    #Multigrid_solver.solve_poisson_periodic(psi, -omega_1, dx, dy, gs_tol)
    Spectral_Poisson.solve_poisson(psi,-omega_1, dx, dy)

    # Step 2
    # Calculate RHS
    f = Fortran_functions.rhs_periodic(psi, omega_1, dx, dy, Re_n)
    #Need to add ML based sgs computation
    sgs = deploy_model_regression(omega_1,psi,ml_model)
    omega_2 = 0.75 * omega + 0.25 * omega_1 + 0.25 * dt * (f + sgs)

    # Fortran update for Poisson Equation
    Spectral_Poisson.solve_poisson(psi,-omega_2, dx, dy)

    # Step 2
    # Calculate RHS
    f = Fortran_functions.rhs_periodic(psi, omega_2, dx, dy, Re_n)
    #Need to add ML based sgs computation
    sgs = deploy_model_regression(omega_2,psi,ml_model)

    for i in range(nx):
        for j in range(ny):
            omega[i, j] = oneth * omega[i, j] + twoth * omega_2[i, j] + twoth * dt * (f[i, j] + sgs[i, j])

    # Fortran update for Poisson Equation
    Spectral_Poisson.solve_poisson(psi,-omega, dx, dy)

def tvdrk3_fortran_ml_nearest_neighbor(omega,psi,ml_model):

    oneth = 1.0 / 3.0
    twoth = 2.0 / 3.0

    # Step 1
    # Calculate RHS
    f = Fortran_functions.rhs_periodic(psi, omega, dx, dy, Re_n)
    #Need to add ML based sgs computation
    sgs = deploy_model_nearest_neighbor(omega,psi,ml_model)

    omega_1 = omega + dt * (f + sgs)

    # Fortran update for Poisson Equation
    #Multigrid_solver.solve_poisson_periodic(psi, -omega_1, dx, dy, gs_tol)
    Spectral_Poisson.solve_poisson(psi,-omega_1, dx, dy)

    # Step 2
    # Calculate RHS
    f = Fortran_functions.rhs_periodic(psi, omega_1, dx, dy, Re_n)
    #Need to add ML based sgs computation
    sgs = deploy_model_nearest_neighbor(omega_1,psi,ml_model)
    omega_2 = 0.75 * omega + 0.25 * omega_1 + 0.25 * dt * (f + sgs)

    # Fortran update for Poisson Equation
    Spectral_Poisson.solve_poisson(psi,-omega_2, dx, dy)

    # Step 2
    # Calculate RHS
    f = Fortran_functions.rhs_periodic(psi, omega_2, dx, dy, Re_n)
    #Need to add ML based sgs computation
    sgs = deploy_model_nearest_neighbor(omega_2,psi,ml_model)

    for i in range(nx):
        for j in range(ny):
            omega[i, j] = oneth * omega[i, j] + twoth * omega_2[i, j] + twoth * dt * (f[i, j] + sgs[i, j])

    # Fortran update for Poisson Equation
    Spectral_Poisson.solve_poisson(psi,-omega, dx, dy)

def tvdrk3_fortran_smag_sgs(omega,psi):

    oneth = 1.0 / 3.0
    twoth = 2.0 / 3.0

    # Step 1
    # Calculate RHS
    f = Fortran_functions.rhs_periodic(psi, omega, dx, dy, Re_n)
    #Standard Smag
    sgs = smag_sgs_calc(omega,psi)

    omega_1 = omega + dt * (f + sgs)

    # Fortran update for Poisson Equation
    #Multigrid_solver.solve_poisson_periodic(psi, -omega_1, dx, dy, gs_tol)
    Spectral_Poisson.solve_poisson(psi,-omega_1, dx, dy)

    # Step 2
    # Calculate RHS
    f = Fortran_functions.rhs_periodic(psi, omega_1, dx, dy, Re_n)
    # Standard Smag
    sgs = smag_sgs_calc(omega_1,psi)
    omega_2 = 0.75 * omega + 0.25 * omega_1 + 0.25 * dt * (f + sgs)

    # Fortran update for Poisson Equation
    Spectral_Poisson.solve_poisson(psi,-omega_2, dx, dy)

    # Step 2
    # Calculate RHS
    f = Fortran_functions.rhs_periodic(psi, omega_2, dx, dy, Re_n)
    # Standard Smag
    sgs = smag_sgs_calc(omega_2,psi)

    for i in range(nx):
        for j in range(ny):
            omega[i, j] = oneth * omega[i, j] + twoth * omega_2[i, j] + twoth * dt * (f[i, j] + sgs[i, j])

    # Fortran update for Poisson Equation
    Spectral_Poisson.solve_poisson(psi,-omega, dx, dy)

def tvdrk3_fortran_dyn_smag_sgs(omega,psi):

    oneth = 1.0 / 3.0
    twoth = 2.0 / 3.0

    # Step 1
    # Calculate RHS
    f = Fortran_functions.rhs_periodic(psi, omega, dx, dy, Re_n)
    #Dynamic Smagorinsky SGS
    sgs = dynamic_smag_calc(omega,psi)

    omega_1 = omega + dt * (f + sgs)

    # Fortran update for Poisson Equation
    #Multigrid_solver.solve_poisson_periodic(psi, -omega_1, dx, dy, gs_tol)
    Spectral_Poisson.solve_poisson(psi,-omega_1, dx, dy)

    # Step 2
    # Calculate RHS
    f = Fortran_functions.rhs_periodic(psi, omega_1, dx, dy, Re_n)
    # Dynamic Smagorinsky SGS
    sgs = dynamic_smag_calc(omega_1,psi)
    omega_2 = 0.75 * omega + 0.25 * omega_1 + 0.25 * dt * (f + sgs)

    # Fortran update for Poisson Equation
    Spectral_Poisson.solve_poisson(psi,-omega_2, dx, dy)

    # Step 2
    # Calculate RHS
    f = Fortran_functions.rhs_periodic(psi, omega_2, dx, dy, Re_n)
    # Dynamic Smagorinsky SGS
    sgs = dynamic_smag_calc(omega_2,psi)

    for i in range(nx):
        for j in range(ny):
            omega[i, j] = oneth * omega[i, j] + twoth * omega_2[i, j] + twoth * dt * (f[i, j] + sgs[i, j])

    # Fortran update for Poisson Equation
    Spectral_Poisson.solve_poisson(psi,-omega, dx, dy)

def tvdrk3_fortran_ml_convolution(omega,psi,ml_model_forward, ml_model_inverse):

    oneth = 1.0 / 3.0
    twoth = 2.0 / 3.0

    # Step 1
    # Calculate RHS
    f = Fortran_functions.rhs_periodic(psi, omega, dx, dy, Re_n)
    #Need to add ML based sgs computation
    sgs = deploy_model_convolution(omega,psi,ml_model_forward, ml_model_inverse)

    omega_1 = omega + dt * (f + sgs)

    # Fortran update for Poisson Equation
    #Multigrid_solver.solve_poisson_periodic(psi, -omega_1, dx, dy, gs_tol)
    Spectral_Poisson.solve_poisson(psi,-omega_1, dx, dy)

    # Step 2
    # Calculate RHS
    f = Fortran_functions.rhs_periodic(psi, omega_1, dx, dy, Re_n)
    #Need to add ML based sgs computation
    sgs = deploy_model_convolution(omega_1,psi,ml_model_forward, ml_model_inverse)
    omega_2 = 0.75 * omega + 0.25 * omega_1 + 0.25 * dt * (f + sgs)

    # Fortran update for Poisson Equation
    Spectral_Poisson.solve_poisson(psi,-omega_2, dx, dy)

    # Step 2
    # Calculate RHS
    f = Fortran_functions.rhs_periodic(psi, omega_2, dx, dy, Re_n)
    #Need to add ML based sgs computation
    sgs = deploy_model_convolution(omega_2,psi,ml_model_forward, ml_model_inverse)

    for i in range(nx):
        for j in range(ny):
            omega[i, j] = oneth * omega[i, j] + twoth * omega_2[i, j] + twoth * dt * (f[i, j] + sgs[i, j])

    # Fortran update for Poisson Equation
    Spectral_Poisson.solve_poisson(psi,-omega, dx, dy)

def tvdrk3_fortran_leith_sgs(omega,psi):

    oneth = 1.0 / 3.0
    twoth = 2.0 / 3.0

    # Step 1
    # Calculate RHS
    f = Fortran_functions.rhs_periodic(psi, omega, dx, dy, Re_n)
    # Standard Leith
    sgs = leith_sgs_calc(omega,psi)

    omega_1 = omega + dt * (f + sgs)

    # Fortran update for Poisson Equation
    #Multigrid_solver.solve_poisson_periodic(psi, -omega_1, dx, dy, gs_tol)
    Spectral_Poisson.solve_poisson(psi,-omega_1, dx, dy)

    # Step 2
    # Calculate RHS
    f = Fortran_functions.rhs_periodic(psi, omega_1, dx, dy, Re_n)
    # Standard Leith
    sgs = leith_sgs_calc(omega_1,psi)
    omega_2 = 0.75 * omega + 0.25 * omega_1 + 0.25 * dt * (f + sgs)

    # Fortran update for Poisson Equation
    Spectral_Poisson.solve_poisson(psi,-omega_2, dx, dy)

    # Step 2
    # Calculate RHS
    f = Fortran_functions.rhs_periodic(psi, omega_2, dx, dy, Re_n)
    # Standard Leith
    sgs = leith_sgs_calc(omega_2,psi)

    for i in range(nx):
        for j in range(ny):
            omega[i, j] = oneth * omega[i, j] + twoth * omega_2[i, j] + twoth * dt * (f[i, j] + sgs[i, j])

    # Fortran update for Poisson Equation
    Spectral_Poisson.solve_poisson(psi,-omega, dx, dy)

def tvdrk3_fortran_dyn_leith_sgs(omega,psi):

    oneth = 1.0 / 3.0
    twoth = 2.0 / 3.0

    # Step 1
    # Calculate RHS
    f = Fortran_functions.rhs_periodic(psi, omega, dx, dy, Re_n)
    #Dynamic Leith SGS
    sgs = dynamic_leith_calc(omega,psi)

    omega_1 = omega + dt * (f + sgs)

    # Fortran update for Poisson Equation
    #Multigrid_solver.solve_poisson_periodic(psi, -omega_1, dx, dy, gs_tol)
    Spectral_Poisson.solve_poisson(psi,-omega_1, dx, dy)

    # Step 2
    # Calculate RHS
    f = Fortran_functions.rhs_periodic(psi, omega_1, dx, dy, Re_n)
    # Dynamic Leith SGS
    sgs = dynamic_leith_calc(omega_1,psi)
    omega_2 = 0.75 * omega + 0.25 * omega_1 + 0.25 * dt * (f + sgs)

    # Fortran update for Poisson Equation
    Spectral_Poisson.solve_poisson(psi,-omega_2, dx, dy)

    # Step 2
    # Calculate RHS
    f = Fortran_functions.rhs_periodic(psi, omega_2, dx, dy, Re_n)
    # Dynamic Leith SGS
    sgs = dynamic_leith_calc(omega_2,psi)

    for i in range(nx):
        for j in range(ny):
            omega[i, j] = oneth * omega[i, j] + twoth * omega_2[i, j] + twoth * dt * (f[i, j] + sgs[i, j])

    # Fortran update for Poisson Equation
    Spectral_Poisson.solve_poisson(psi,-omega, dx, dy)

def tvdrk3_fortran_ml_feature_sgs(omega,psi,ml_model):

    oneth = 1.0 / 3.0
    twoth = 2.0 / 3.0

    # Step 1
    # Calculate RHS
    f = Fortran_functions.rhs_periodic(psi, omega, dx, dy, Re_n)
    #Need to add ML based sgs computation - pointwise
    sgs = deploy_model_feature_regression(omega,psi,ml_model)

    omega_1 = omega + dt * (f + sgs)

    # Fortran update for Poisson Equation
    #Multigrid_solver.solve_poisson_periodic(psi, -omega_1, dx, dy, gs_tol)
    Spectral_Poisson.solve_poisson(psi,-omega_1, dx, dy)

    # Step 2
    # Calculate RHS
    f = Fortran_functions.rhs_periodic(psi, omega_1, dx, dy, Re_n)
    # Need to add ML based sgs computation - pointwise
    sgs = deploy_model_feature_regression(omega_1,psi,ml_model)
    omega_2 = 0.75 * omega + 0.25 * omega_1 + 0.25 * dt * (f + sgs)

    # Fortran update for Poisson Equation
    Spectral_Poisson.solve_poisson(psi,-omega_2, dx, dy)

    # Step 2
    # Calculate RHS
    f = Fortran_functions.rhs_periodic(psi, omega_2, dx, dy, Re_n)
    # Need to add ML based sgs computation - pointwise
    sgs = deploy_model_feature_regression(omega_2,psi,ml_model)

    for i in range(nx):
        for j in range(ny):
            omega[i, j] = oneth * omega[i, j] + twoth * omega_2[i, j] + twoth * dt * (f[i, j] + sgs[i, j])

    # Fortran update for Poisson Equation
    Spectral_Poisson.solve_poisson(psi,-omega, dx, dy)

def tvdrk3_fortran(omega,psi):

    oneth = 1.0/3.0
    twoth = 2.0/3.0

    #Step 1
    #Calculate RHS
    f = Fortran_functions.rhs_periodic(psi,omega, dx, dy, Re_n)
    omega_1 = omega + dt*(f)

    #Fortran update for Poisson Equation
    #Multigrid_solver.solve_poisson_periodic(psi,-omega_1, dx, dy, gs_tol)
    Spectral_Poisson.solve_poisson(psi,-omega_1, dx, dy)

    #Step 2
    #Calculate RHS
    f = Fortran_functions.rhs_periodic(psi,omega_1, dx, dy, Re_n)
    omega_2 = 0.75*omega + 0.25*omega_1 + 0.25*dt*(f)

    #Fortran update for Poisson Equation
    Spectral_Poisson.solve_poisson(psi,-omega_2, dx, dy)

    #Step 2
    #Calculate RHS
    f = Fortran_functions.rhs_periodic(psi,omega_2, dx, dy, Re_n)

    for i in range(nx):
        for j in range(ny):
            omega[i,j] = oneth*omega[i,j] + twoth*omega_2[i,j] + twoth*dt*(f[i,j])

    #Fortran update for Poisson Equation
    Spectral_Poisson.solve_poisson(psi,-omega, dx, dy)

def tvdrk3_fortran_ad_les(omega,psi):

    oneth = 1.0/3.0
    twoth = 2.0/3.0

    #Step 1
    #Calculate RHS
    f = Fortran_functions.rhs_periodic(psi,omega, dx, dy, Re_n)
    sfs = ad_sfs_calc(omega,psi)
    omega_1 = omega + dt*(f+sfs)

    #Fortran update for Poisson Equation
    # Multigrid_solver.solve_poisson_periodic(psi,-omega_1, dx, dy, gs_tol)
    Spectral_Poisson.solve_poisson(psi,-omega_1, dx, dy)

    #Step 2
    #Calculate RHS
    f = Fortran_functions.rhs_periodic(psi,omega_1, dx, dy, Re_n)
    sfs = ad_sfs_calc(omega_1, psi)
    omega_2 = 0.75*omega + 0.25*omega_1 + 0.25*dt*(f+sfs)

    #Fortran update for Poisson Equation
    Spectral_Poisson.solve_poisson(psi,-omega_2, dx, dy)

    #Step 2
    #Calculate RHS
    f = Fortran_functions.rhs_periodic(psi,omega_2, dx, dy, Re_n)
    sfs = ad_sfs_calc(omega_2, psi)

    for i in range(nx):
        for j in range(ny):
            omega[i,j] = oneth*omega[i,j] + twoth*omega_2[i,j] + twoth*dt*(f[i,j]+sfs[i,j])

    #Fortran update for Poisson Equation
    Spectral_Poisson.solve_poisson(psi,-omega, dx, dy)

#-------------------------------------------------------------------------------------#
#-------------------------------------------------------------------------------------#
#Main time integrator
#-------------------------------------------------------------------------------------#
#-------------------------------------------------------------------------------------#

def main_func():

    #Initialize my domain and constants
    omega, psi = init_domain()

    initialize_ic_bc(omega,psi)

    if closure_choice == 2 or closure_choice == 4 or closure_choice == 7 or closure_choice==9:#For the ML SFS, ML SGS variants
        ml_model = load_pretrained_model()
    elif closure_choice == 8:
        ml_model_forward, ml_model_inverse = load_pretrained_model()

    t = 0.0

    clock_time_init = time.clock()
    # Defining numpy array for storing sigma
    sigma_history = np.asarray([sigma], dtype='double')

    for tstep in range(nt):

        t = t + dt
        print(t)

        #TVD - RK3 Fortran
        if closure_choice == 0:
            tvdrk3_fortran(omega,psi)
        elif closure_choice == 1:#AD LES
            tvdrk3_fortran_ad_les(omega,psi)
        elif closure_choice == 2:#AD ML LES - classification
            tvdrk3_fortran_ml_ad_les(omega,psi,ml_model)
            print ('Sigma value of :',sigma)
            sigma_history = np.concatenate((sigma_history, np.asarray([sigma])), axis=0)
        elif closure_choice == 3:#RF LES
            tvdrk3_fortran_rf_les(omega,psi)
        elif closure_choice == 4:#ML SGS
            tvdrk3_fortran_ml_sgs(omega,psi,ml_model)
        elif closure_choice == 5:#Smag SGS
            tvdrk3_fortran_smag_sgs(omega,psi)
        elif closure_choice == 6:#Leith SGS
            tvdrk3_fortran_leith_sgs(omega,psi)
        elif closure_choice == 7:#Nearest neighbor LES
            tvdrk3_fortran_ml_nearest_neighbor(omega,psi,ml_model)
        elif closure_choice == 8:#ML Based deconvolution
            tvdrk3_fortran_ml_convolution(omega,psi,ml_model_forward, ml_model_inverse)
        elif closure_choice == 9:#ML based feature SGS
            tvdrk3_fortran_ml_feature_sgs(omega, psi, ml_model)
        elif closure_choice == 10:#Dynamic Smagorinsky
            tvdrk3_fortran_dyn_smag_sgs(omega, psi)
        elif closure_choice == 11:#Dynamic Leith
            tvdrk3_fortran_dyn_leith_sgs(omega, psi)

        if np.isnan(np.sum(omega))==1:
            print('overflow')
            exit()


    total_clock_time = time.clock()-clock_time_init

    print('Total Clock Time = ',total_clock_time)

    post_process(omega,psi)

    if closure_choice == 2:
        np.savetxt('Sigma_history.txt',sigma_history)
        plt.figure()
        plt.plot(sigma_history[:])
        plt.show()

#-------------------------------------------------------------------------------------#
#-------------------------------------------------------------------------------------#
##### RUN HERE #####
#-------------------------------------------------------------------------------------#
#-------------------------------------------------------------------------------------#
main_func()

