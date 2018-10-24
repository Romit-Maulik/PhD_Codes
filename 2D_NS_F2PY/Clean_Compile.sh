find . -name "*.so" -exec rm {} \;

#OpenMP compiled
#f2py -m Fortran_functions --fcompiler=gfortran --f90flags='-fopenmp' -lgomp -c Fortran_functions.f95
#f2py -m Approximate_Deconvolution --fcompiler=gfortran --f90flags='-fopenmp' -lgomp -c Approximate_Deconvolution.f95
#f2py -m Multigrid_solver --fcompiler=gfortran --f90flags='-fopenmp' -lgomp -c Multigrid_solver.f95


#Serial compiled
f2py -c --fcompiler=gfortran -m Fortran_functions Fortran_functions.f95
f2py -c --fcompiler=gfortran -m Fortran_Dynamic Fortran_Dynamic.f95
f2py -c --fcompiler=gfortran -m Multigrid_solver Multigrid_solver.f95
f2py -c --fcompiler=gfortran -m Gauss_Siedel Gauss_Siedel.f95
f2py -c --fcompiler=gfortran -m Spectral_Poisson Spectral_Poisson.f95
f2py -c --fcompiler=gfortran -m Fortran_Gaussian_Filter Fortran_Gaussian_Filter.f95
f2py -c --fcompiler=gfortran -m Fortran_Pade_Filter Fortran_Pade_Filter.f95
f2py -c --fcompiler=gfortran -m Approximate_Deconvolution Approximate_Deconvolution.f95
f2py -c --fcompiler=gfortran -m Ml_convolution Ml_convolution.f95

#Change this according to your distribution of python
for x in *.cpython-36m-x86_64-linux-gnu.so; do mv "$x" "${x%.cpython-36m-x86_64-linux-gnu.so}.so"; done