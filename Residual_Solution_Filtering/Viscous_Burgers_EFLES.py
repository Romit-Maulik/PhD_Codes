import numpy as np
import matplotlib.pyplot as plt

# Viscosity
Re_n = 5000

# Spatial
nx = 1024
nr = 32
xlength = 2.0*np.pi
dx = xlength / float(nx)

# Temporal
ft = 1.0
nt = 10000
dt = ft/nt

#EFLES Choice
residual_smoothing = 0 # Activated if 1, make sure to switch solution smoothing to opposite
solution_smoothing = 1 # Activated if 1, make sure to switch residual smoothing to opposite
sigma = 1.0 # Filter strength

def plot_results(u,espec):
    fig, ax = plt.subplots(nrows=1,ncols=2)
    ax[0].set_title('VBE Solution field')
    ax[0].plot(u[:])

    ax[1].set_title('Energy spectrum')
    ax[1].loglog(espec[:])

    plt.show()

def initialize():  # Initializes in Fourier domain - one time cost so don't hate me for the for loops
    pi = np.pi

    dx1 = dx / float(nr)

    nx_new = nx * nr
    xlength_new = dx1 * float(nx_new)

    r1 = np.arange(0, nx_new // 2,step=1)
    r2 = np.asarray([0])
    r3 = np.arange(-nx_new // 2,0,step=1)

    wavevector = np.concatenate((r1,r2,r3),axis=0)
    kx = np.array([(2 * pi) * i / xlength_new for i in wavevector])  # Note no im heref
    acons = 2.0 / (10.0 ** (5.0)) / (3.0 * (pi ** (0.5)))

    array_hat = np.zeros(nx_new, dtype=np.complex)  # The array is of type complex
    phase = np.zeros(2 * nx_new, dtype='double')

    np.random.seed(0)
    rand = np.random.uniform(0.0, 1.0)
    phase[0] = np.cos(2.0 * pi * rand)
    phase[1] = 0.0
    phase[nx_new] = np.cos(2.0 * pi * rand)
    phase[nx_new + 1] = 0.0

    k = 3
    for i in range(1, nx_new // 2):
        rand = np.random.uniform(0.0, 1.0)
        phase[k - 1] = np.cos(2.0 * pi * rand)
        phase[k] = np.sin(2.0 * pi * rand)
        phase[2 * nx_new - k + 1] = np.cos(2.0 * pi * rand)
        phase[2 * nx_new - k + 2] = -np.sin(2.0 * pi * rand)
        k = k + 2

    k = 0
    for i in range(0, nx_new):
        espec_ip = np.exp(-(kx[i] / 10.0) ** (2.0))
        espec = acons * (kx[i] ** 4) * espec_ip
        array_hat[i] = nx_new * (np.sqrt(2.0 * espec) * (phase[k] + phase[k + 1]))
        k = k + 2

    temp_array = np.real(np.fft.ifft(array_hat))
    copy_array = np.zeros(nx, dtype='double')

    for i in range(0, nx):
        copy_array[i] = temp_array[i * nr]

    del temp_array
    return copy_array

def rhs_calc(u):
    '''
    :param u: The solution field at any given instant
    :return: the RHS for explicit time integration
    '''

    uud = 0.5*first_derivative(u*u)
    us = second_derivative(u)

    # plot_field(uud)
    # exit()

    rhs = - uud + 1.0/(Re_n)*us

    return rhs

def first_derivative(u):
    '''
    :param u: The solution field at any given instant
    :return: The first derivative for a periodic domain
    We shall utilize 6th order explicit finite difference calculations for the derivative
    '''
    utemp = np.zeros(shape=(nx+6),dtype='double')# 3 extra ghost points

    # Copy array
    utemp[3:nx+3] = u[:]

    #Update ghost points
    utemp[0:3] = u[nx-3:nx]
    utemp[nx+3:nx+6] = u[0:3]

    # Coefficients - first derivative
    am3 = -1.0/60.0
    am2 = 9.0/60.0
    am1 = -45.0/60.0
    ap1 = 45.0/60.0
    ap2 = -9.0/60.0
    ap3 = 1.0/60.0

    ud = 1.0 / dx * (
                am3 * utemp[0:nx] + am2 * utemp[1:nx + 1] + am1 * utemp[2:nx + 2] + ap1 * utemp[4:nx + 4] + ap2 * utemp[
                                                                                                                  5:nx + 5] + ap3 * utemp[
                                                                                                                                    6:nx + 6])

    del utemp
    return ud

def second_derivative(u):
    '''
    :param u: The solution field at any given instant
    :return: The second derivative for a periodic domain
    We shall utilize 6th order explicit finite difference calculations for the derivative
    '''
    utemp = np.zeros(shape=(nx+6),dtype='double')# 3 extra ghost points

    # Copy array
    utemp[3:nx+3] = u[:]

    #Update ghost points
    utemp[0:3] = u[nx-3:nx]
    utemp[nx+3:nx+6] = u[0:3]

    # Coefficients - second derivative
    am3 = 2.0/180.0
    am2 = -27.0/180.0
    am1 = 270.0/180.0
    a0 = -490.0/180.0
    ap1 = 270.0/180.0
    ap2 = -27.0/180.0
    ap3 = 2.0/180.0

    ud = 1.0 / (dx*dx) * (
                am3 * utemp[0:nx] + am2 * utemp[1:nx + 1] + am1 * utemp[2:nx + 2] + a0*utemp[3:nx+3] + ap1 * utemp[4:nx + 4] + ap2 * utemp[
                                                                                                                  5:nx + 5] + ap3 * utemp[
                                                                                                                                    6:nx + 6])

    del utemp
    return ud

def filter(u):
    '''
    :param u: Any field at any given instant
    :return: Filtered field
    Maulik and San, Fluids 2017 for reference of filter
    '''
    utemp = np.zeros(shape=(nx + 6), dtype='double')  # 3 extra ghost points

    # Copy array
    utemp[3:nx + 3] = u[:]

    # Update ghost points
    utemp[0:3] = u[nx - 3:nx]
    utemp[nx + 3:nx + 6] = u[0:3]

    # Filter coefficients
    f0 = 5.0/16.0
    f1 = -15.0/64.0
    f2 = 3.0/32.0
    f3 = -1.0/64.0

    uf = utemp[3:nx+3]-sigma*(f0*utemp[3:nx+3]+f1*(utemp[2:nx+2]+utemp[4:nx+4])+f2*(utemp[1:nx+1]+utemp[5:nx+5])+f3*(utemp[0:nx]+utemp[6:nx+6]))

    return uf

def spectra_calculation(u):
    # Take Fourier transform
    array_hat = np.fft.fft(u)

    # Normalizing data
    array_new = np.copy(array_hat / float(nx))
    # Energy Spectrum
    espec = 0.5 * np.absolute(array_new)**2
    # Angle Averaging
    eplot = np.zeros(nx // 2, dtype='double')
    for i in np.arange(1, nx // 2):
        eplot[i] = 0.5 * (espec[i] + espec[nx - i])

    return eplot

def time_integrator(u):
    '''
    :param u: Initial condition of u
    :return: Returns final u at a specified time
    '''

    t = 0
    while t < ft:
        t = t + dt

        # Substep 1
        rhs = rhs_calc(u)

        if residual_smoothing == 1:
            filter(rhs)

        u_1 = u + dt*rhs

        if solution_smoothing == 1:
            filter(u_1)

        #Substep 2
        rhs = rhs_calc(u_1)

        if residual_smoothing == 1:
            filter(rhs)

        u_1 = 0.75*u + 0.25*u_1 + 0.25*dt*rhs

        if solution_smoothing == 1:
            filter(u_1)

        #Substep 3
        rhs = rhs_calc(u_1)

        if residual_smoothing == 1:
            filter(rhs)

        u[:] = 1.0/3.0*u[:] + 2.0/3.0*u_1[:] + 2.0/3.0*dt*rhs[:]

        if solution_smoothing == 1:
            filter(u)

    return u


if __name__ == '__main__':

    #Initialize Burgers turbulence
    u = initialize()
    e_init = spectra_calculation(u)
    plot_results(u,e_init)

    # Time integration
    uf = time_integrator(u)

    # Final time results
    e_fin = spectra_calculation(u)
    plot_results(uf,e_fin)