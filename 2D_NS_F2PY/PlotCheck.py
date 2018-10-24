import numpy as np
import matplotlib.pyplot as plt
import Fortran_functions

def load_field(filename,data_size):
    #Loading data from Tecplot file (check number of skiprows) - no tecplot header
    fine_data = np.loadtxt(filename)

    # Store in numpy array
    vort= np.arange((data_size + 1)* (data_size + 1), dtype='double').reshape(data_size + 1, data_size + 1)

    iter = 0
    for j in range(0, data_size + 1):
        for i in range(0, data_size + 1):
            vort[i, j] = fine_data[iter, 0]
            iter = iter + 1

    vortper = np.arange((data_size) * (data_size), dtype='double').reshape(data_size, data_size)

    for j in range(0, data_size):
        for i in range(0, data_size):
            vortper[i, j] = vort[i, j]

    return vortper


#Load first
kplot_udns  = np.load('Spectra_UDNS.npy')[0]
eplot_udns = np.load('Spectra_UDNS.npy')[1]
eplot_ml_dcd = np.load('Spectra_ML_DCD.npy')[1]

eplot_ml_bs_1 = np.load('Spectra_BS_1.npy')[1]
eplot_ml_bs_2 = np.load('Spectra_BS_2.npy')[1]
eplot_ml_bs_3 = np.load('Spectra_BS_3.npy')[1]
eplot_ml_bs_4 = np.load('Spectra_BS_4.npy')[1]
eplot_ml_bs_5 = np.load('Spectra_BS_5.npy')[1]
eplot_ml_bs_6 = np.load('Spectra_BS_6.npy')[1]

#DNS
kplot_DNS = np.load('DNS_32_Spectra.npy')[0]
eplot_DNS = np.load('DNS_32_Spectra.npy')[1]


#Fourier filtered DNS
# data_size = 256
# coarse_vort = load_field('Source.504_32k',data_size)
# arr_len = int(0.5*np.sqrt(float(data_size*data_size + data_size*data_size)))-1
# eplot_FDNS = np.zeros(arr_len+1,dtype='double')
# kplot_FDNS = np.arange(0,arr_len+1,1,dtype='double')
# Fortran_functions.spec(coarse_vort,eplot_FDNS)
# #

scale_plot = np.array([[10,0.1],[100,1.0e-4]])

plt.loglog(kplot_DNS,eplot_DNS,label=r'DNS')
plt.loglog(kplot_udns,eplot_udns,label=r'UNS',color='red')
plt.loglog(kplot_udns,eplot_ml_dcd,label='DCD No Backscatter',color='orange')

plt.loglog(kplot_udns,eplot_ml_bs_1,label='DCD BS 1')
plt.loglog(kplot_udns,eplot_ml_bs_2,label='DCD BS 2')
plt.loglog(kplot_udns,eplot_ml_bs_3,label='DCD BS 3')
plt.loglog(kplot_udns,eplot_ml_bs_4,label='DCD BS 4')
plt.loglog(kplot_udns,eplot_ml_bs_5,label='DCD BS 5')
plt.loglog(kplot_udns,eplot_ml_bs_6,label='DCD BS 6')

plt.loglog(scale_plot[:,0],scale_plot[:,1],color='black',linestyle='dashed',label=r'$k^{-3}$ scaling')
plt.xlim([1,500])
plt.ylim([1e-8,1])
plt.xlabel('k',fontsize='16')
plt.ylabel('E(k)',fontsize='16')

plt.tick_params(axis='both', labelsize='12')
plt.legend(fontsize='12',loc='upper right')
plt.show()

