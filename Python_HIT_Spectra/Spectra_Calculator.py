import numpy as np
import os
import matplotlib.pyplot as plt
from scipy.ndimage import gaussian_filter

#Function to output numpy 3D array in tecplot format
def tecplot_op(u,v,w,filename,labels,data_size):

    x = np.zeros((data_size+1),dtype='double')
    y = np.zeros((data_size+1), dtype='double')
    z = np.zeros((data_size+1), dtype='double')

    for i in range(1, data_size+1):
        x[i] = float(i) / data_size * 2.0 * np.pi
        y[i] = float(i) / data_size * 2.0 * np.pi
        z[i] = float(i) / data_size * 2.0 * np.pi


    f = open(filename,"w+")
    f.write('variables = "x", "y", "z",'+labels+'\n')
    f.write('zone T=Zone_ i=          65 j=          65 k=          65 f=point\n')

    output_array = np.zeros(((data_size+1)**3,6),dtype='double')
    iter = 0
    for k in range(data_size+1):
        for j in range(data_size+1):
            for i in range(data_size+1):
                output_array[iter, 0] = x[i]
                output_array[iter, 1] = y[j]
                output_array[iter, 2] = z[k]
                output_array[iter, 3] = u[i,j,k]
                output_array[iter, 4] = v[i,j,k]
                output_array[iter, 5] = w[i,j,k]
                iter = iter + 1

    np.savetxt(f,output_array,fmt='%10.7e')
    f.close()
    del output_array

def load_data(filename,data_size):
    # Loading data from Tecplot file (check number of skiprows)
    fine_data = np.loadtxt(filename, skiprows=2)
    x = np.zeros((data_size + 1), dtype='double')
    y = np.zeros((data_size + 1), dtype='double')
    z = np.zeros((data_size + 1), dtype='double')

    iter = 0
    for i in range(0, data_size):
        x[i] = fine_data[iter, 0]
        y[i] = fine_data[iter, 0]
        z[i] = fine_data[iter, 0]
        iter = iter + 1

    # Periodic BC
    x[data_size] = 2.0 * np.pi
    y[data_size] = 2.0 * np.pi
    z[data_size] = 2.0 * np.pi

    # Store in numpy array
    ufine = np.arange((data_size + 1) * (data_size + 1) * (data_size + 1), dtype='double').reshape(data_size + 1,
                                                                                                   data_size + 1,
                                                                                                   data_size + 1)
    vfine = np.arange((data_size + 1) * (data_size + 1) * (data_size + 1), dtype='double').reshape(data_size + 1,
                                                                                                   data_size + 1,
                                                                                                   data_size + 1)
    wfine = np.arange((data_size + 1) * (data_size + 1) * (data_size + 1), dtype='double').reshape(data_size + 1,
                                                                                                   data_size + 1,
                                                                                                   data_size + 1)

    iter = 0
    for k in range(0, data_size + 1):
        for j in range(0, data_size + 1):
            for i in range(0, data_size + 1):
                ufine[i, j, k] = fine_data[iter, 3]
                vfine[i, j, k] = fine_data[iter, 4]
                wfine[i, j, k] = fine_data[iter, 5]
                iter = iter + 1


    return ufine, vfine, wfine


def filter_operation(ufine,sigma):
    ufil = gaussian_filter(ufine,sigma)

    return ufil


def energy_spectra_calculation(u,v,w,data_size,labels):#3D energy spectra calculation
    #real domain
    ureal = np.arange((data_size) * (data_size) * (data_size), dtype='double').reshape(data_size,
                                                                                       data_size,
                                                                                       data_size)

    vreal = np.arange((data_size) * (data_size) * (data_size), dtype='double').reshape(data_size,
                                                                                       data_size,
                                                                                       data_size)

    wreal = np.arange((data_size) * (data_size) * (data_size), dtype='double').reshape(data_size,
                                                                                       data_size,
                                                                                       data_size)
    for k in range(0, data_size):
        for j in range(0, data_size):
            for i in range(0, data_size):
                ureal[i, j, k] = u[i, j, k]
                vreal[i, j, k] = v[i, j, k]
                wreal[i, j, k] = w[i, j, k]


    ufreqdom = np.fft.fftn(ureal)
    vfreqdom = np.fft.fftn(vreal)
    wfreqdom = np.fft.fftn(wreal)
    #Normalizing
    ufreqdomnorm = np.copy(ufreqdom / float(data_size ** 3))
    vfreqdomnorm = np.copy(vfreqdom / float(data_size ** 3))
    wfreqdomnorm = np.copy(wfreqdom / float(data_size ** 3))
    #Energy calculation
    espec = 0.5 * (np.absolute(ufreqdomnorm)**2 + np.absolute(vfreqdomnorm)**2 + np.absolute(wfreqdomnorm)**2)

    #Angle Averaging
    arr_len = int(np.sqrt(float(3*(data_size/2)**2)))
    eplot = np.zeros(arr_len+1,dtype='double')
    kplot = np.arange(0,arr_len+1,1,dtype='double')

    #Wavenumbers
    kx = np.array([i for i in
                   list(range(0, data_size // 2)) + list([0]) + list(range(-data_size // 2 + 1, 0))],dtype='double')
    ky = np.copy(kx)
    kz = np.copy(kx)

    print (kx)


    for p in range(1,arr_len+1):
        eplot[p]=0.0
        for k in range(0, data_size):
            for j in range(0, data_size):
                for i in range(0, data_size):
                    kr = np.sqrt(kx[i]*kx[i] + ky[j]*ky[j] + kz[k]*kz[k])

                    if kr>=float(p)-0.5 and kr<float(p)+0.5:
                        eplot[p] = eplot[p] + espec[i,j,k]


    plt.loglog(kplot,eplot,label=labels)



if __name__ == "__main__":
    #Data size
    data_size = 64

    u,v,w = load_data('TGV_1600_3_Coarse.plt',data_size)

    plt.figure()

    energy_spectra_calculation(u,v,w,data_size,'True Spectra')

    uf = filter_operation(u, 1.0)
    vf = filter_operation(v, 1.0)
    wf = filter_operation(w, 1.0)

    #tecplot_op(uf,vf,wf,'Filtered_Fields.plt','"ufil", "vfil", "wfil"', data_size)

    energy_spectra_calculation(uf,vf,wf, data_size, 'Filtered Spectra')

    x = np.array([2, 200])
    y = np.array([10, 0.0046415888])
    plt.loglog(x, y, 'r--',label='Scaling')
    plt.legend()
    plt.show()
