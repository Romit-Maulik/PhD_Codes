!------------------------------------------------------------------
!Calculates right-hand-side for TVD RK3 implementations
!Second-order Arakawa Jacobian calculation
!Validated
!------------------------------------------------------------------
subroutine rhs_periodic(psi,omega,f,nx,ny,dx,dy,Re_N)
!$ use omp_lib
implicit none


integer, intent(in) :: nx, ny
integer :: i, j
double precision, intent(in) :: dx, dy, Re_N
double precision, dimension(0:nx-1,0:ny-1),intent(in)::psi,omega
double precision, dimension(0:nx-1,0:ny-1),intent(out)::f
double precision :: jj1, jj2, jj3, d2wdy2, d2wdx2

double precision, dimension(:,:), allocatable :: psi_new, omega_new

allocate(psi_new(-1:nx,-1:ny))
allocate(omega_new(-1:nx,-1:ny))

!$OMP PARALLEL DO
do j = 0,ny-1
do i = 0,nx-1
psi_new(i,j) = psi(i,j)
omega_new(i,j) = omega(i,j)
end do
end do
!$OMP END PARALLEL DO

call periodic_bc_update(nx,ny,psi_new)
call periodic_bc_update(nx,ny,omega_new)

!$OMP PARALLEL DO
do j = 0,ny-1
do i = 0,nx-1


jj1 = 1.0/(4.0*dx*dy) * ((omega_new(i+1,j)-omega_new(i-1,j)) * (psi_new(i,j+1) - psi_new(i,j-1)) &
			- (omega_new(i,j+1)-omega_new(i,j-1)) * (psi_new(i+1,j) - psi_new(i-1,j)))

jj2 = 1.0 / (4.0 * dx * dy) * (omega_new(i+1, j) * (psi_new(i+1, j+1) - psi_new(i+1, j-1)) &
                                         - omega_new(i-1, j) * (psi_new(i-1, j+1) - psi_new(i-1, j-1)) &
                                         - omega_new(i, j+1) * (psi_new(i+1, j+1) - psi_new(i-1, j+1)) &
                                         + omega_new(i, j-1) * (psi_new(i+1, j-1) - psi_new(i-1, j-1)) &
                                          )

jj3 = 1.0 / (4.0 * dx * dy) * (omega_new(i+1, j+1) * (psi_new(i, j+1) - psi_new(i+1, j)) &
                                        -  omega_new(i-1, j-1) * (psi_new(i-1, j) - psi_new(i, j-1)) &
                                        -  omega_new(i-1, j+1) * (psi_new(i, j+1) - psi_new(i-1, j)) &
                                        +  omega_new(i+1, j-1) * (psi_new(i+1, j) - psi_new(i, j-1)) &
                                          )

d2wdy2 = (omega_new(i, j+1) + omega_new(i, j-1) - 2.0 * omega_new(i, j)) / (dy * dy)
d2wdx2 = (omega_new(i+1, j) + omega_new(i-1, j) - 2.0 * omega_new(i, j)) / (dx * dx)

f(i, j) = (-(jj1 + jj2 + jj3)/3.0 + 1.0 / Re_n * (d2wdy2 + d2wdx2))

end do
end do
!$OMP END PARALLEL DO

end subroutine



!---------------------------------------------------------------------------!
!subroutine - boundary condition update
!Validated
!---------------------------------------------------------------------------!
subroutine periodic_bc_update(nx,ny,u)
implicit none

integer :: nx, ny, i, j
double precision, dimension(-1:nx,-1:ny) :: u

do i = 0,nx-1
u(i,-1) = u(i,ny-1)
u(i,ny) = u(i,0)
end do

do j = -1,ny
u(-1,j) = u(nx-1,j)
u(nx,j) = u(0,j)
end do

end subroutine



!---------------------------------------------------------------------------!
!Compute 2D vorticity field from the energy spectrum
!Periodic, equidistant grid
!Validated
!---------------------------------------------------------------------------!
subroutine hit_init_cond(w_org,nx,ny,dx,dy)
implicit none
integer,intent(inout)::nx,ny
double precision,intent(inout) ::w_org(0:nx-1,0:ny-1)
double precision, intent(in) :: dx, dy
integer :: ni,nj,ii,jj
double precision,dimension(:,:),allocatable ::w
double precision ::ran,pi,kk,E4
double precision,parameter:: tiny=1.0d-10
double precision,allocatable ::data1d(:),phase2d(:,:,:),ksi(:,:),eta(:,:)
double precision,allocatable ::kx(:),ky(:),ww(:,:)
integer::i,j,k,isign,ndim,nn(2),seed

allocate(w(-2:nx+2,-2:ny+2))

w = 0.0d0

seed = 19

!expand it to dns grid

ni = nx
nj = ny

nx = 2048
ny = 2048

ii = nx/ni
jj = ny/nj

ndim =2
nn(1)=nx
nn(2)=ny

allocate(kx(0:nx-1),ky(0:ny-1))
allocate(ksi(0:nx/2,0:ny/2),eta(0:nx/2,0:ny/2))
allocate(data1d(2*nx*ny))
allocate(phase2d(2,0:nx-1,0:ny-1))
allocate(ww(0:nx,0:ny))

!Set seed for the random number generator between [0,1]
CALL RANDOM_SEED(seed)

pi = 4.0d0*datan(1.0d0)


!Wave numbers 
do i=0,nx/2-1
kx(i)      = dfloat(i)
kx(i+nx/2) = dfloat(i-nx/2)
end do
kx(0) = tiny

do j=0,ny/2-1
ky(j)      = dfloat(j)
ky(j+ny/2) = dfloat(j-ny/2)
end do
ky(0) = tiny

!Random numbers in the first quadrant
do j=0,ny/2
do i=0,nx/2
CALL RANDOM_NUMBER(ran)
ksi(i,j) =2.0d0*pi*ran
end do
end do

do j=0,ny/2
do i=0,nx/2
CALL RANDOM_NUMBER(ran)
eta(i,j) =2.0d0*pi*ran
end do
end do

!Random phase
do j=0,ny-1
do i=0,nx-1
phase2d(1,i,j)       = 0.0d0
phase2d(2,i,j)       = 0.0d0
end do
end do
  
do j=1,ny/2-1
do i=1,nx/2-1
!I.st
phase2d(1,i,j)       = dcos(ksi(i,j)+eta(i,j)) 
phase2d(2,i,j)       = dsin(ksi(i,j)+eta(i,j)) 
!II.nd
phase2d(1,nx-i,j)    = dcos(-ksi(i,j)+eta(i,j)) 
phase2d(2,nx-i,j)    = dsin(-ksi(i,j)+eta(i,j)) 
!IV.th
phase2d(1,i,ny-j)    = dcos(ksi(i,j)-eta(i,j)) 
phase2d(2,i,ny-j)    = dsin(ksi(i,j)-eta(i,j)) 
!III.rd
phase2d(1,nx-i,ny-j) = dcos(-ksi(i,j)-eta(i,j)) 
phase2d(2,nx-i,ny-j) = dsin(-ksi(i,j)-eta(i,j)) 
end do
end do


!vorticity amplitudes in Fourier space 
k=1
do j=0,ny-1
do i=0,nx-1   
    kk = dsqrt(kx(i)*kx(i) + ky(j)*ky(j))
    data1d(k)   =  dsqrt(kk*E4(kk)/pi)*phase2d(1,i,j)
	data1d(k+1) =  dsqrt(kk*E4(kk)/pi)*phase2d(2,i,j)   
k = k + 2
end do
end do

!find the velocity in physical space
!forward fourier transform
isign= 1
call fourn(data1d,nn,ndim,isign)


k=1
do j=0,ny-1
do i=0,nx-1
ww(i,j)=data1d(k)
k=k+2
end do
end do


! periodicity
do j=0,ny-1
ww(nx,j)=ww(0,j)
end do
do i=0,nx
ww(i,ny)=ww(i,0)
end do

!back to the local grid
nx = ni
ny = nj
do j=0,ny
do i=0,nx
w(i,j)=ww(i*ii,j*jj)
end do
end do

deallocate(data1d,phase2d,ksi,eta,ww)


do j = 0,ny-1
do i = 0,nx-1
w_org(i,j) = w(i,j)
end do
end do

deallocate(w)


return
end


!---------------------------------------------------------------------------!
!Given energy spectrum
!Used for initial field calculation for 2D HIT
!Validated
!---------------------------------------------------------------------------!
double precision function E4(kr)
implicit none
double precision:: kr,pi,c,k0
k0 = 10.0d0
pi = 4.0d0*datan(1.0d0)
c = 4.0d0/(3.0d0*dsqrt(pi)*(k0**5))
!c = 1.0d0/(4.0d0*pi*(k0**6))
!c = 1.0d0/(2.0d0*pi*(k0**6))
E4 = c*(kr**4)*dexp(-(kr/k0)**2)
end



!-----------------------------------------------------------------!
! fft routine from numerical recipes
! ndim: dimension of the transform (i.e.; 2 for 2d problems)
! nn  : number of points in each direction
! data: one-dimensional array including real and imaginary part 
!-----------------------------------------------------------------!
subroutine fourn(data,nn,ndim,isign)
implicit none
integer:: ndim,isign
integer:: nn(ndim)
real*8:: data(*)
real*8:: wr,wi,wpr,wpi,wtemp,theta,tempr,tempi
integer::ntot,n,nrem,nprev,idim,ip1,ip2,ip3,i1,i2,i3
integer::i2rev,i3rev,ibit,ifp1,ifp2,k1,k2

      ntot=1
      do 11 idim=1,ndim
        ntot=ntot*nn(idim)
11    continue
      nprev=1
      do 18 idim=1,ndim
        n=nn(idim)
        nrem=ntot/(n*nprev)
        ip1=2*nprev
        ip2=ip1*n
        ip3=ip2*nrem
        i2rev=1
        do 14 i2=1,ip2,ip1
          if(i2.lt.i2rev)then
            do 13 i1=i2,i2+ip1-2,2
              do 12 i3=i1,ip3,ip2
                i3rev=i2rev+i3-i2
                tempr=data(i3)
                tempi=data(i3+1)
                data(i3)=data(i3rev)
                data(i3+1)=data(i3rev+1)
                data(i3rev)=tempr
                data(i3rev+1)=tempi
12            continue
13          continue
          endif
          ibit=ip2/2
1         if ((ibit.ge.ip1).and.(i2rev.gt.ibit)) then
            i2rev=i2rev-ibit
            ibit=ibit/2
          go to 1
          endif
          i2rev=i2rev+ibit
14      continue
        ifp1=ip1
2       if(ifp1.lt.ip2)then
          ifp2=2*ifp1
          theta=isign*6.28318530717959d0/(ifp2/ip1)
          wpr=-2.d0*dsin(0.5d0*theta)**2
          wpi=dsin(theta)
          wr=1.d0
          wi=0.d0
          do 17 i3=1,ifp1,ip1
            do 16 i1=i3,i3+ip1-2,2
              do 15 i2=i1,ip3,ifp2
                k1=i2
                k2=k1+ifp1
                tempr=sngl(wr)*data(k2)-sngl(wi)*data(k2+1)
                tempi=sngl(wr)*data(k2+1)+sngl(wi)*data(k2)
                data(k2)=data(k1)-tempr
                data(k2+1)=data(k1+1)-tempi
                data(k1)=data(k1)+tempr
                data(k1+1)=data(k1+1)+tempi
15            continue
16          continue
            wtemp=wr
            wr=wr*wpr-wi*wpi+wr
            wi=wi*wpr+wtemp*wpi+wi
17        continue
          ifp1=ifp2
        go to 2
        endif
        nprev=n*nprev
18    continue

return
end

!-----------------------------------------------------------------!
!Compute energy spectra and send back to python - postprocessing
!Validated
!-----------------------------------------------------------------!
subroutine spec(w_org,nx,ny,eplot,n)
implicit none
integer,intent(in) ::nx,ny,n
double precision, intent(in)::w_org(0:nx-1,0:ny-1)
double precision::pi
double precision,dimension(0:n),intent(inout) :: eplot
integer::i,j,k,ic
double precision::kx(0:nx-1),ky(0:ny-1),kk
double precision,parameter:: tiny=1.0d-10
double precision,dimension(:),allocatable:: data1d
double precision,dimension(:,:),allocatable::es
integer,parameter::ndim=2
integer::nn(ndim),isign

allocate(data1d(2*nx*ny))

pi = 4.0d0*datan(1.0d0)

nn(1)= nx
nn(2)= ny

!finding fourier coefficients of w 
!invese fourier transform
!find the vorticity in Fourier space
k=1
do j=0,ny-1  
do i=0,nx-1   
  data1d(k)   =  w_org(i,j)
  data1d(k+1) =  0.0d0    
k = k + 2
end do
end do
!normalize
do k=1,2*nx*ny
data1d(k)=data1d(k)/dfloat(nx*ny)
end do
!inverse fourier transform
isign= -1
call fourn(data1d,nn,ndim,isign)


!Wave numbers 
do i=0,nx/2-1
kx(i)      = dfloat(i)
kx(i+nx/2) = dfloat(i-nx/2)
end do
kx(0) = tiny

do j=0,ny/2-1
ky(j)      = dfloat(j)
ky(j+ny/2) = dfloat(j-ny/2)
end do
ky(0) = tiny

!Energy spectrum (for all wavenumbers)
allocate(es(0:nx-1,0:ny-1))
k=1
do j=0,ny-1
do i=0,nx-1 
kk = dsqrt(kx(i)*kx(i) + ky(j)*ky(j))
es(i,j) = pi*(data1d(k)*data1d(k) + data1d(k+1)*data1d(k+1))/kk
k = k + 2
end do
end do

!Plot angle averaged energy spectrum
do k=1,n
eplot(k) = 0.0d0
ic = 0
do j=1,ny-1
do i=1,nx-1
kk = dsqrt(kx(i)*kx(i) + ky(j)*ky(j))
    if(kk.ge.(dfloat(k)-0.5d0).and.kk.le.(dfloat(k)+0.5d0)) then
    ic = ic + 1
    eplot(k) = eplot(k) + es(i,j)
    end if
end do
end do
eplot(k) = eplot(k) / dfloat(ic)
end do

deallocate(data1d,es)

return
end 


!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!subroutine for making a random sampling matrix from omega field
!Used for ML-AD classification
!Validated
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
subroutine field_sampler(omega,nx,ny,sampling_matrix,n_samples)

!$ use omp_lib
implicit none

integer, intent(in) :: nx,ny,n_samples
integer :: i,j,sample,seed
double precision,dimension(0:nx-1,0:ny-1),intent(in) :: omega
double precision,dimension(0:n_samples-1,0:8),intent(out) :: sampling_matrix
double precision :: ran_real, omega_max, omega_min
double precision, dimension(:,:), allocatable :: omega_norm

!Set seed for the random number generator between [0,1]
CALL RANDOM_SEED(seed)


!Normalize omega
allocate(omega_norm(0:nx-1,0:ny-1))
omega_max = maxval(omega)
omega_min = minval(omega)

do j = 0,ny-1
  do i = 0,nx-1
    omega_norm(i,j) = (omega(i,j)-omega_min)/(omega_max-omega_min)
  end do
end do


!$OMP PARALLEL DO
do sample = 0,n_samples-1

  !Real random number - we avoid boundaries
  call RANDOM_NUMBER(ran_real)
  !Map to integers
  i = 1+FLOOR((nx-1)*ran_real)
  !Real random number
  call RANDOM_NUMBER(ran_real)
  !Map to integers
  j = 1+FLOOR((ny-1)*ran_real)


  sampling_matrix(sample,0) = omega_norm(i,j)
  sampling_matrix(sample,1) = omega_norm(i,j+1)
  sampling_matrix(sample,2) = omega_norm(i,j-1)
  sampling_matrix(sample,3) = omega_norm(i+1,j)
  sampling_matrix(sample,4) = omega_norm(i+1,j+1)
  sampling_matrix(sample,5) = omega_norm(i+1,j-1)
  sampling_matrix(sample,6) = omega_norm(i-1,j)
  sampling_matrix(sample,7) = omega_norm(i-1,j+1)
  sampling_matrix(sample,8) = omega_norm(i-1,j-1)

  ! sampling_matrix(sample,0) = omega(i,j)
  ! sampling_matrix(sample,1) = omega(i+1,j)
  ! sampling_matrix(sample,2) = omega(i-1,j)
  ! sampling_matrix(sample,3) = omega(i,j+1)
  ! sampling_matrix(sample,4) = omega(i+1,j+1)
  ! sampling_matrix(sample,5) = omega(i-1,j+1)
  ! sampling_matrix(sample,6) = omega(i,j-1)
  ! sampling_matrix(sample,7) = omega(i+1,j-1)
  ! sampling_matrix(sample,8) = omega(i-1,j-1)

end do
!$OMP END PARALLEL DO


deallocate(omega_norm)

return
end




!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!Subroutine for sampling across domain for SGS predictions using ML
!Used for ML-SGS prediction - input preparation for JFM Rapids
!Validated
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
subroutine sgs_field_sampler(omega,psi,sampling_matrix,nx,ny,dx,dy)

!$ use omp_lib
implicit none
integer, intent(in) :: nx,ny
integer :: i,j, sample
double precision,dimension(0:nx-1,0:ny-1),intent(in) :: omega, psi
double precision, intent(in) :: dx, dy
double precision,dimension(0:nx*ny,0:19),intent(out) :: sampling_matrix
double precision, dimension(:,:), allocatable :: omega_new, psi_new
double precision, dimension(:,:), allocatable :: dsdy, d2sdx, d2sdy, d2sdxdy,dwdx,dwdy
double precision,dimension(:,:,:),allocatable :: inv
!double precision :: vort_mean, stream_mean, inv1_mean, inv2_mean
!double precision :: vort_std, stream_std, inv1_std, inv2_std


allocate(omega_new(-1:nx,-1:ny))
allocate(psi_new(-1:nx,-1:ny))

do j = 0,ny-1
  do i = 0,nx-1
    omega_new(i,j) = (omega(i,j))
    psi_new(i,j) = psi(i,j)
  end do
end do

!BC Update
call periodic_bc_update(nx,ny,omega_new)
call periodic_bc_update(nx,ny,psi_new)

!Calculating Smag Turbulence model invariants
allocate(d2sdxdy(0:nx-1,0:ny-1))
allocate(d2sdx(0:nx-1,0:ny-1))
allocate(d2sdy(0:nx-1,0:ny-1))
allocate(dsdy(-1:nx,-1:ny))

do j = 0,ny-1
  do i = 0,nx-1
    dsdy(i,j) = (psi_new(i,j+1)-psi_new(i,j-1))/(2.0d0*dy)
    d2sdx(i,j) = (psi_new(i+1,j)+psi_new(i-1,j)-2.0d0*psi_new(i,j))/(dx*dx)
    d2sdy(i,j) = (psi_new(i,j+1)+psi_new(i,j-1)-2.0d0*psi_new(i,j))/(dy*dy)
  end do
end do

call periodic_bc_update(nx,ny,dsdy)!Need for a second derivative for this quantity

do j = 0,ny-1
  do i = 0,nx-1
    d2sdxdy(i,j) = (dsdy(i+1,j)-dsdy(i-1,j))/(2.0d0*dx)
  end do
end do

allocate(inv(0:nx,0:ny,0:1))
!Smag invariant
do j = 0,ny-1
  do i = 0,nx-1
    inv(i,j,0) = dsqrt(4.0d0*d2sdxdy(i,j)**2 + (d2sdx(i,j)-d2sdy(i,j))**2)
  end do
end do

deallocate(d2sdxdy,dsdy,d2sdx,d2sdy)

!Calculating Leith turbulence model invariants
allocate(dwdy(0:nx-1,0:ny-1))
allocate(dwdx(0:nx-1,0:ny-1))
do j = 0,ny-1
do i = 0,nx-1
dwdy(i,j) = (omega_new(i,j+1)-omega_new(i,j-1))/(2.0d0*dy)
dwdx(i,j) = (omega_new(i+1,j)-omega_new(i-1,j))/(2.0d0*dx)
end do
end do

!Leith invariant
do j = 0,ny-1
  do i = 0,nx-1
    inv(i,j,1) = dsqrt(dwdx(i,j)**2 + dwdy(i,j)**2)
  end do
end do

deallocate(dwdx,dwdy)

!!Calculating means for input vectors
!vort_mean = 0.0d0
!stream_mean = 0.0d0
!inv1_mean = 0.0d0
!inv2_mean = 0.0d0

!do j = 0,ny-1
!  do i = 0,nx-1
!  	vort_mean = vort_mean + omega_new(i,j)
!  	stream_mean = stream_mean + psi_new(i,j)
!  	inv1_mean = inv1_mean + inv(i,j,0)
!  	inv2_mean = inv2_mean + inv(i,j,1)
!  end do
!end do

!vort_mean = vort_mean/dfloat(nx*ny)
!stream_mean = stream_mean/dfloat(nx*ny)
!inv1_mean = inv1_mean/dfloat(nx*ny)
!inv2_mean = inv2_mean/dfloat(nx*ny)

!!Calculating sdevs for input vectors
!vort_std = 0.0d0
!stream_std = 0.0d0
!inv1_std = 0.0d0
!inv2_std = 0.0d0

!do j = 0,ny-1
!  do i = 0,nx-1
!  	vort_std = vort_std + (omega_new(i,j)-vort_mean)**2
!  	stream_std = stream_std + (psi_new(i,j)-stream_mean)**2
!  	inv1_std = inv1_std + (inv(i,j,0)-inv1_mean)**2
!  	inv2_std = inv2_std + (inv(i,j,1)-inv2_mean)**2
!  end do
!end do

!vort_std = dsqrt(vort_std/dfloat(nx*ny-1))
!stream_std = dsqrt(stream_std/dfloat(nx*ny-1))
!inv1_std = dsqrt(inv1_std/dfloat(nx*ny-1))
!inv2_std = dsqrt(inv2_std/dfloat(nx*ny-1))


!!Changing input data arrays - normalization
!do j = 0,ny-1
!  do i = 0,nx-1
!  	omega_new(i,j) = (omega_new(i,j) - vort_mean)/vort_std
!	psi_new(i,j) = (psi_new(i,j) - stream_mean)/stream_std
!	inv(i,j,0) = (inv(i,j,0) - inv1_mean)/inv1_std
!	inv(i,j,1) = (inv(i,j,1) - inv2_mean)/inv2_std
!  end do
!end do

sample = 0
!Preparing sample matrix
!$OMP PARALLEL DO
do j = 0,ny-1
  do i = 0,nx-1
    sampling_matrix(sample,0) = omega_new(i,j)
    sampling_matrix(sample,1) = omega_new(i,j+1)
    sampling_matrix(sample,2) = omega_new(i,j-1)
    sampling_matrix(sample,3) = omega_new(i+1,j)
    sampling_matrix(sample,4) = omega_new(i+1,j+1)
    sampling_matrix(sample,5) = omega_new(i+1,j-1)
    sampling_matrix(sample,6) = omega_new(i-1,j)
    sampling_matrix(sample,7) = omega_new(i-1,j+1)
    sampling_matrix(sample,8) = omega_new(i-1,j-1)

    sampling_matrix(sample,9) = psi_new(i,j)
    sampling_matrix(sample,10) = psi_new(i,j+1)
    sampling_matrix(sample,11) = psi_new(i,j-1)
    sampling_matrix(sample,12) = psi_new(i+1,j)
    sampling_matrix(sample,13) = psi_new(i+1,j+1)
    sampling_matrix(sample,14) = psi_new(i+1,j-1)
    sampling_matrix(sample,15) = psi_new(i-1,j)
    sampling_matrix(sample,16) = psi_new(i-1,j+1)
    sampling_matrix(sample,17) = psi_new(i-1,j-1)

    sampling_matrix(sample,18) = inv(i,j,0)
    sampling_matrix(sample,19) = inv(i,j,1)

    sample = sample + 1

  end do
end do
!$OMP END PARALLEL DO

deallocate(omega_new,psi_new,inv)

return
end


!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!Subroutine for sampling across domain for SGS predictions using ML
!Used for Feature-SGS prediction//No stencil, only pointwise
!Validated
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
subroutine field_feature_sampler(omega,psi,sampling_matrix,nx,ny,dx,dy)

!$ use omp_lib
implicit none
integer, intent(in) :: nx,ny
integer :: i,j, sample
double precision,dimension(0:nx-1,0:ny-1),intent(in) :: omega, psi
double precision, intent(in) :: dx, dy
double precision,dimension(0:nx*ny,0:8),intent(out) :: sampling_matrix
double precision, dimension(:,:), allocatable :: omega_new, psi_new, lapc
double precision, dimension(:,:), allocatable :: dsdy, d2sdx, d2sdy, d2sdxdy,dwdx,dwdy
double precision,dimension(:,:,:),allocatable :: inv
double precision :: inv1_mean, inv2_mean, inv3_mean, inv4_mean, lap_val
double precision :: inv1_std, inv2_std, inv3_std, inv4_std
double precision, dimension(:,:),allocatable :: s11, s12, s22, r12


allocate(omega_new(-1:nx,-1:ny))
allocate(psi_new(-1:nx,-1:ny))

do j = 0,ny-1
  do i = 0,nx-1
    omega_new(i,j) = (omega(i,j))
    psi_new(i,j) = psi(i,j)
  end do
end do

!BC Update
call periodic_bc_update(nx,ny,omega_new)
call periodic_bc_update(nx,ny,psi_new)

!Calculating Smag Turbulence model invariants
allocate(d2sdxdy(0:nx-1,0:ny-1))
allocate(d2sdx(0:nx-1,0:ny-1))
allocate(d2sdy(0:nx-1,0:ny-1))
allocate(dsdy(-1:nx,-1:ny))

do j = 0,ny-1
  do i = 0,nx-1
    dsdy(i,j) = (psi_new(i,j+1)-psi_new(i,j-1))/(2.0d0*dy)
    d2sdx(i,j) = (psi_new(i+1,j)+psi_new(i-1,j)-2.0d0*psi_new(i,j))/(dx*dx)
    d2sdy(i,j) = (psi_new(i,j+1)+psi_new(i,j-1)-2.0d0*psi_new(i,j))/(dy*dy)
  end do
end do

call periodic_bc_update(nx,ny,dsdy)!Need for a second derivative for this quantity

do j = 0,ny-1
  do i = 0,nx-1
    d2sdxdy(i,j) = (dsdy(i+1,j)-dsdy(i-1,j))/(2.0d0*dx)
  end do
end do

allocate(inv(0:nx,0:ny,0:3))
!Smag invariant
do j = 0,ny-1
  do i = 0,nx-1
    inv(i,j,0) = dsqrt(4.0d0*d2sdxdy(i,j)**2 + (d2sdx(i,j)-d2sdy(i,j))**2)
  end do
end do



!Calculating Leith turbulence model invariants
allocate(dwdy(0:nx-1,0:ny-1))
allocate(dwdx(0:nx-1,0:ny-1))
do j = 0,ny-1
do i = 0,nx-1
dwdy(i,j) = (omega_new(i,j+1)-omega_new(i,j-1))/(2.0d0*dy)
dwdx(i,j) = (omega_new(i+1,j)-omega_new(i-1,j))/(2.0d0*dx)
end do
end do

!Leith invariant
do j = 0,ny-1
  do i = 0,nx-1
    inv(i,j,1) = dsqrt(dwdx(i,j)**2 + dwdy(i,j)**2)
  end do
end do

!BL invariant
do j = 0,ny-1
  do i = 0,nx-1
    inv(i,j,2) = dabs(omega_new(i,j))
  end do
end do

!CS invariant
do j = 0,ny-1
  do i = 0,nx-1
    inv(i,j,3) = dsqrt(d2sdx(i,j)**2 + d2sdy(i,j)**2)
  end do
end do


deallocate(dwdx,dwdy)
deallocate(d2sdxdy,dsdy,d2sdx,d2sdy)

!Calculating means for input vectors
inv1_mean = 0.0d0
inv2_mean = 0.0d0
inv3_mean = 0.0d0
inv4_mean = 0.0d0

do j = 0,ny-1
  do i = 0,nx-1
  	inv1_mean = inv1_mean + inv(i,j,0)
  	inv2_mean = inv2_mean + inv(i,j,1)
	inv3_mean = inv3_mean + inv(i,j,2)
  	inv4_mean = inv4_mean + inv(i,j,3)
  end do
end do


inv1_mean = inv1_mean/dfloat(nx*ny)
inv2_mean = inv2_mean/dfloat(nx*ny)
inv3_mean = inv3_mean/dfloat(nx*ny)
inv4_mean = inv4_mean/dfloat(nx*ny)


!Calculating sdevs for input vectors

inv1_std = 0.0d0
inv2_std = 0.0d0
inv3_std = 0.0d0
inv4_std = 0.0d0

do j = 0,ny-1
  do i = 0,nx-1
  	inv1_std = inv1_std + (inv(i,j,0)-inv1_mean)**2
  	inv2_std = inv2_std + (inv(i,j,1)-inv2_mean)**2
  	inv3_std = inv3_std + (inv(i,j,2)-inv3_mean)**2
  	inv4_std = inv4_std + (inv(i,j,3)-inv4_mean)**2
  end do
end do


inv1_std = dsqrt(inv1_std/dfloat(nx*ny-1))
inv2_std = dsqrt(inv2_std/dfloat(nx*ny-1))
inv3_std = dsqrt(inv3_std/dfloat(nx*ny-1))
inv4_std = dsqrt(inv4_std/dfloat(nx*ny-1))


!Changing input data arrays - normalization
do j = 0,ny-1
  do i = 0,nx-1
	inv(i,j,0) = (inv(i,j,0) - inv1_mean)/inv1_std
	inv(i,j,1) = (inv(i,j,1) - inv2_mean)/inv2_std
	inv(i,j,2) = (inv(i,j,2) - inv3_mean)/inv3_std
	inv(i,j,3) = (inv(i,j,3) - inv4_mean)/inv4_std
  end do
end do


!Finding field features
allocate(s11(0:nx-1,0:ny-1))
allocate(s12(0:nx-1,0:ny-1))
allocate(s22(0:nx-1,0:ny-1))
allocate(r12(0:nx-1,0:ny-1))

call strain_rotation_tensor_calc(nx,ny,psi,s11,s12,s22,r12,dx,dy)

!Find laplacian
allocate(lapc(0:nx-1,0:ny-1))

!Calculating laplacian
do j = 0,ny-1
  do i = 0,nx-1
    lap_val = (omega_new(i+1,j)+omega_new(i-1,j)-2.0d0*omega_new(i,j))/(dx*dx)
    lapc(i,j) = lap_val + (omega_new(i,j+1)+omega_new(i,j-1)-2.0d0*omega_new(i,j))/(dy*dy)
  end do
end do

!Normalizing Laplacian to zero mean and unit variance
!Calculating means for laplacian - reusing inv1_mean, inv1_std
inv1_mean = 0.0d0
do j = 0,ny-1
  do i = 0,nx-1
  	inv1_mean = inv1_mean + lapc(i,j)
  end do
end do
inv1_mean = inv1_mean/dfloat(nx*ny)
!Calculating sdevs for input vectors
inv1_std = 0.0d0
do j = 0,ny-1
  do i = 0,nx-1
  	inv1_std = inv1_std + (lapc(i,j)-inv1_mean)**2
  end do
end do

inv1_std = dsqrt(inv1_std/dfloat(nx*ny-1))
!Changing input data arrays - normalization
do j = 0,ny-1
  do i = 0,nx-1
	lapc(i,j) = (lapc(i,j) - inv1_mean)/inv1_std
  end do
end do


sample = 0
!Preparing sample matrix
!$OMP PARALLEL DO
do j = 0,ny-1
  do i = 0,nx-1

    sampling_matrix(sample,0) = s11(i,j)
    sampling_matrix(sample,1) = s12(i,j)
    sampling_matrix(sample,2) = s22(i,j)
    sampling_matrix(sample,3) = r12(i,j)

    sampling_matrix(sample,4) = inv(i,j,0)
    sampling_matrix(sample,5) = inv(i,j,1)
    sampling_matrix(sample,6) = inv(i,j,2)
    sampling_matrix(sample,7) = inv(i,j,3)

    sampling_matrix(sample,8) = lapc(i,j)

    sample = sample + 1

  end do
end do
!$OMP END PARALLEL DO

deallocate(omega_new,psi_new,inv)
deallocate(s11,s12,s22,r12,lapc)

return
end


!------------------------------------------------------------------------------!
!------------------------------------------------------------------------------!
!subroutine to calculate strain and rotation tensors
!Outputs normalized to zero mean and unit variance
!------------------------------------------------------------------------------!
!------------------------------------------------------------------------------!
subroutine strain_rotation_tensor_calc(nxc,nyc,sc,s11,s12,s22,r12,dx,dy)
	implicit none

	integer :: nxc,nyc, i, j
	double precision, dimension(0:nxc-1,0:nyc-1) :: sc,s11,s12,s22,r12
	double precision :: dx, dy

	double precision, dimension(:,:), allocatable :: u, v, sc_temp
	double precision :: s11_m, s12_m, s22_m, r12_m
	double precision :: s11_s, s12_s, s22_s, r12_s

	allocate(u(-1:nxc,-1:nyc))
	allocate(v(-1:nxc,-1:nyc))
	allocate(sc_temp(-1:nxc,-1:nyc))

	do j=0,nyc-1
		do i = 0,nxc-1
			sc_temp(i,j)=sc(i,j)
		end do
	end do

	call periodic_bc_update(nxc,nyc,sc_temp)

	!Calculate u and v
	do j=0,nyc-1
		do i = 0,nxc-1
			u(i,j)=(sc_temp(i,j+1)-sc_temp(i,j-1))/(2.0d0*dy)
			v(i,j)=-(sc_temp(i+1,j)-sc_temp(i-1,j))/(2.0d0*dx)
		end do
	end do

	call periodic_bc_update(nxc,nyc,u)
	call periodic_bc_update(nxc,nyc,v)

	!Calculate strain and rotation components - 2d
	do j=0,nyc-1
		do i = 0,nxc-1
			s11(i,j) = (u(i+1,j)-u(i-1,j))/(2.0d0*dx)
			s22(i,j) = (v(i,j+1)-v(i,j-1))/(2.0d0*dy)
			s12(i,j) = 0.5d0*((u(i,j+1)-u(i,j-1))/(2.0d0*dy) + (v(i+1,j)-v(i-1,j))/(2.0d0*dx))
			r12(i,j) = 0.5d0*((u(i,j+1)-u(i,j-1))/(2.0d0*dy) - (v(i+1,j)-v(i-1,j))/(2.0d0*dx))
		end do
	end do

	!Normalize for deployment
	!Calculating means for input vectors
	s11_m = 0.0d0
	s12_m = 0.0d0
	s22_m = 0.0d0
	r12_m = 0.0d0

	do j = 0,nyc-1
	  do i = 0,nxc-1
	  	s11_m = s11_m + s11(i,j)
	  	s12_m = s12_m + s12(i,j)
	  	s22_m = s22_m + s22(i,j)
	  	r12_m = r12_m + r12(i,j)
	  end do
	end do


	s11_m = s11_m/dfloat(nxc*nyc)
	s12_m = s12_m/dfloat(nxc*nyc)
	s22_m = s22_m/dfloat(nxc*nyc)
	r12_m = r12_m/dfloat(nxc*nyc)

	!Calculating sdevs for input vectors

	s11_s = 0.0d0
	s12_s = 0.0d0
	s22_s = 0.0d0
	r12_s = 0.0d0

	do j = 0,nyc-1
	  do i = 0,nxc-1
	  	s11_s = s11_s + (s11(i,j)-s11_m)**2
	  	s12_s = s12_s + (s12(i,j)-s12_m)**2
	  	s22_s = s22_s + (s22(i,j)-s22_m)**2
	  	r12_s = r12_s + (r12(i,j)-r12_m)**2
	  end do
	end do


	s11_s = dsqrt(s11_s/dfloat(nxc*nyc-1))
	s12_s = dsqrt(s12_s/dfloat(nxc*nyc-1))
	s22_s = dsqrt(s22_s/dfloat(nxc*nyc-1))
	r12_s = dsqrt(r12_s/dfloat(nxc*nyc-1))

	!Changing input data arrays - normalization
	do j = 0,nyc-1
	  do i = 0,nxc-1
		s11(i,j) = (s11(i,j) - s11_m)/s11_s
		s12(i,j) = (s12(i,j) - s12_m)/s12_s
		s22(i,j) = (s22(i,j) - s22_m)/s22_s
		r12(i,j) = (r12(i,j) - r12_m)/r12_s
	  end do
	end do


	deallocate(u,v,sc_temp)

	return
end

!--------------------------------------------------------------------------
!Calculates laplacian for given 2D matrix input from python
!Validated
!--------------------------------------------------------------------------
subroutine laplacian_calculator(omega,laplacian,nx,ny,dx,dy)

!$ use omp_lib
implicit none
integer, intent(in) :: nx,ny
integer :: i,j
double precision,dimension(0:nx-1,0:ny-1),intent(in) :: omega
double precision,dimension(0:nx-1,0:ny-1),intent(out) :: laplacian
double precision,intent(in) :: dx, dy
double precision, dimension(:,:), allocatable :: omega_new
double precision :: lap_val

!Normalize omega
allocate(omega_new(-1:nx,-1:ny))

do j = 0,ny-1
  do i = 0,nx-1
    omega_new(i,j) = (omega(i,j))
  end do
end do

!BC Update
call periodic_bc_update(nx,ny,omega_new)

!Calculating laplacian
do j = 0,ny-1
  do i = 0,nx-1
    lap_val = (omega_new(i+1,j)+omega_new(i-1,j)-2.0d0*omega_new(i,j))/(dx*dx)
    laplacian(i,j) = lap_val + (omega_new(i,j+1)+omega_new(i,j-1)-2.0d0*omega_new(i,j))/(dy*dy)
  end do
end do

deallocate(omega_new)

return
end

!--------------------------------------------------------------------------
!Subroutine for reshaping Keras prediction into nx x ny shape
!The return matrix here contains values for source term (Pi from Keras)
!Validated
!--------------------------------------------------------------------------

subroutine sgs_reshape_non_sdev(return_matrix,nrows,ncols,sgs,laplacian,nx,ny,source_max,source_min)

!$ use omp_lib
implicit none

integer, intent(in) :: nx,ny,nrows,ncols
integer :: i,j, sample
double precision, dimension(0:nrows,ncols),intent(in) :: return_matrix
double precision, dimension(0:nx-1,0:ny-1),intent(inout) :: sgs
double precision, dimension(0:nx-1,0:ny-1),intent(in) :: laplacian
double precision, intent(in) :: source_max, source_min

sample = 0
!Reading return matrix
!$OMP PARALLEL DO
do j = 0,ny-1
  do i = 0,nx-1
	!This if statement to ensure no negative numerical viscosities
	!Key for JFM Rapids manuscript
  	if (laplacian(i,j)<0 .and. return_matrix(sample,1)<0) then
  		sgs(i,j) = return_matrix(sample,1)
  	else if (laplacian(i,j)>0 .and. return_matrix(sample,1)>0) then
  		sgs(i,j) = return_matrix(sample,1)
  	else
  		sgs(i,j) = 0.0d0
  	end if

    sample = sample + 1
  
  end do
end do
!$OMP END PARALLEL DO


return
end

!--------------------------------------------------------------------------
!Subroutine for reshaping Keras prediction into nx x ny shape
!The return matrix here contains values for source term (Pi from Keras)
!A form of local averaging to statistically represent backscatter in a stencil
!--------------------------------------------------------------------------

subroutine sgs_reshape_backscatter(return_matrix,nrows,ncols,sgs,laplacian,nx,ny,source_max,source_min)

!$ use omp_lib
implicit none

integer, intent(in) :: nx,ny,nrows,ncols
integer :: i,j, sample, k
double precision, dimension(0:nrows,ncols),intent(in) :: return_matrix
double precision, dimension(0:nx-1,0:ny-1),intent(inout) :: sgs
double precision, dimension(0:nx-1,0:ny-1),intent(in) :: laplacian
double precision, intent(in) :: source_max, source_min


!Temp variables in Fortran
double precision, dimension(:,:), allocatable :: sgs_temp, lap_temp
double precision, dimension(1:9) :: nu, lap
double precision :: sgs_val, lap_val, nu_av, tiny_val, lap_av


tiny_val = 1.0d-10

allocate(sgs_temp(-1:nx,-1:ny))
allocate(lap_temp(-1:nx,-1:ny))

sample = 0
!Reading return matrix
!$OMP PARALLEL DO
do j = 0,ny-1
  do i = 0,nx-1
	sgs_temp(i,j) = return_matrix(sample,1)
	lap_temp(i,j) = laplacian(i,j)
	sample = sample + 1
  end do
end do
!$OMP END PARALLEL DO

!Boundary conditions update
!BC Update
call periodic_bc_update(nx,ny,sgs_temp)
call periodic_bc_update(nx,ny,lap_temp)

!Local averaging
!$OMP PARALLEL DO
do j = 0,ny-1
  do i = 0,nx-1

	!Calculate eddy viscosities
	sgs_val = sgs_temp(i,j)
	lap_val = lap_temp(i,j)
	nu(1) = sgs_val/lap_val
	lap(1) = lap_val

	sgs_val = sgs_temp(i,j+1)
	lap_val = lap_temp(i,j+1)
	nu(2) = sgs_val/lap_val
	lap(2) = lap_val

	sgs_val = sgs_temp(i,j-1)
	lap_val = lap_temp(i,j-1)
	nu(3) = sgs_val/lap_val
	lap(3) = lap_val

	sgs_val = sgs_temp(i+1,j)
	lap_val = lap_temp(i+1,j)
	nu(4) = sgs_val/lap_val
	lap(4) = lap_val

	sgs_val = sgs_temp(i+1,j+1)
	lap_val = lap_temp(i+1,j+1)
	nu(5) = sgs_val/lap_val
	lap(5) = lap_val

	sgs_val = sgs_temp(i+1,j-1)
	lap_val = lap_temp(i+1,j-1)
	nu(6) = sgs_val/lap_val
	lap(6) = lap_val

	sgs_val = sgs_temp(i-1,j)
	lap_val = lap_temp(i-1,j)
	nu(7) = sgs_val/lap_val
	lap(7) = lap_val

	sgs_val = sgs_temp(i-1,j+1)
	lap_val = lap_temp(i-1,j+1)
	nu(8) = sgs_val/lap_val
	lap(8) = lap_val

	sgs_val = sgs_temp(i-1,j-1)
	lap_val = lap_temp(i-1,j-1)
	nu(9) = sgs_val/lap_val
	lap(9) = lap_val

	nu_av = 0.0d0
	lap_av = 0.0d0
	do k = 1,9
	    nu_av = nu_av + nu(k)
	    lap_av = lap_av + lap(k)
	end do
	nu_av = nu_av/9.0d0
	lap_av = lap_av/9.0d0


	if (nu(1)>tiny_val.and.nu_av>nu(1)) then
           sgs(i,j) = nu(1)*laplacian(i,j)
	else
	   sgs(i,j) = 0.0d0
	end if

 
  end do
end do
!$OMP END PARALLEL DO

deallocate(sgs_temp,lap_temp)


return
end

!--------------------------------------------------------------------------
!Subroutine for reshaping Keras prediction into nx x ny shape
!The return matrix here contains values for source term (Pi from Keras)
!Validated
!--------------------------------------------------------------------------

subroutine sgs_reshape(return_matrix,nrows,ncols,sgs,laplacian,nx,ny,source_max,source_min)

!$ use omp_lib
implicit none

integer, intent(in) :: nx,ny,nrows,ncols
integer :: i,j, sample
double precision, dimension(0:nrows,ncols),intent(in) :: return_matrix
double precision, dimension(0:nx-1,0:ny-1),intent(inout) :: sgs
double precision, dimension(0:nx-1,0:ny-1),intent(in) :: laplacian
double precision, intent(in) :: source_max, source_min
double precision :: mean, sdev

sample = 0
mean = 0.0d0
!Reading return matrix
!$OMP PARALLEL DO
do j = 0,ny-1
  do i = 0,nx-1
	sgs(i,j) = return_matrix(sample,1)
	mean = mean + sgs(i,j)
    	sample = sample + 1
  end do
end do

mean = mean /dfloat(nx*ny)

!Calculating standard deviation of prediction
sdev = 0.0d0
do j = 0,ny-1
  do i = 0,nx-1
	sgs(i,j) = return_matrix(sample,1)
	sdev = sdev + (sgs(i,j)-mean)**2
    	sample = sample + 1
  end do
end do

sdev = dsqrt(sdev/dfloat(nx*ny-1))

sample = 0
!Reading return matrix
!$OMP PARALLEL DO
do j = 0,ny-1
  do i = 0,nx-1
	!This if statement to ensure no negative numerical viscosities
	!Key for JFM Rapids manuscript
  	if (laplacian(i,j)<0 .and. return_matrix(sample,1)<0) then
  		sgs(i,j) = return_matrix(sample,1)
	  	sgs(i,j) = min(return_matrix(sample,1),sdev)
  	else if (laplacian(i,j)>0 .and. return_matrix(sample,1)>0) then
  		sgs(i,j) = return_matrix(sample,1)
		sgs(i,j) = min(return_matrix(sample,1),sdev)
	else
		sgs(i,j) = 0.0d0
 	end if


    sample = sample + 1
  
  end do
end do
!$OMP END PARALLEL DO


return
end


!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!Subroutine for smagorinsky based SGS prediction
!Validated
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
subroutine smag_source_term(omega,psi,sgs,nx,ny,dx,dy)

!$ use omp_lib
implicit none
integer, intent(in) :: nx,ny
integer :: i,j
double precision,dimension(0:nx-1,0:ny-1),intent(in) :: omega, psi
double precision, intent(in) :: dx, dy
double precision,dimension(0:nx-1,0:ny-1),intent(out) :: sgs
double precision, dimension(:,:), allocatable :: omega_new, psi_new
double precision, dimension(:,:), allocatable :: dsdy, d2sdx, d2sdy, d2sdxdy
double precision :: c_turb, lap_val

allocate(omega_new(-1:nx,-1:ny))
allocate(psi_new(-1:nx,-1:ny))

do j = 0,ny-1
  do i = 0,nx-1
    omega_new(i,j) = (omega(i,j))
    psi_new(i,j) = psi(i,j)
  end do
end do

!BC Update
call periodic_bc_update(nx,ny,omega_new)
call periodic_bc_update(nx,ny,psi_new)

!Calculating Smag Turbulence model invariants
allocate(d2sdxdy(0:nx-1,0:ny-1))
allocate(d2sdx(0:nx-1,0:ny-1))
allocate(d2sdy(0:nx-1,0:ny-1))
allocate(dsdy(-1:nx,-1:ny))

do j = 0,ny-1
  do i = 0,nx-1
    dsdy(i,j) = (psi_new(i,j+1)-psi_new(i,j-1))/(2.0d0*dy)
    d2sdx(i,j) = (psi_new(i+1,j)+psi_new(i-1,j)-2.0d0*psi_new(i,j))/(dx*dx)
    d2sdy(i,j) = (psi_new(i,j+1)+psi_new(i,j-1)-2.0d0*psi_new(i,j))/(dy*dy)
  end do
end do

call periodic_bc_update(nx,ny,dsdy)!Need for a second derivative for this quantity

do j = 0,ny-1
  do i = 0,nx-1
    d2sdxdy(i,j) = (dsdy(i+1,j)-dsdy(i-1,j))/(2.0d0*dx)
  end do
end do

!Smag invariant
do j = 0,ny-1
  do i = 0,nx-1
    sgs(i,j) = dsqrt(4.0d0*d2sdxdy(i,j)**2 + (d2sdx(i,j)-d2sdy(i,j))**2)
  end do
end do

deallocate(d2sdxdy,dsdy,d2sdx,d2sdy)

!Combining invariant with turbulence coefficients and laplacian
c_turb = 1.0d0*dx
do j = 0,ny-1
  do i = 0,nx-1
    lap_val = (omega_new(i+1,j)+omega_new(i-1,j)-2.0d0*omega_new(i,j))/(dx*dx)
    lap_val = lap_val + (omega_new(i,j+1)+omega_new(i,j-1)-2.0d0*omega_new(i,j))/(dy*dy)
    sgs(i,j) = c_turb*c_turb*sgs(i,j)*lap_val
  end do
end do

deallocate(omega_new,psi_new)

return
end


!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!Subroutine for smagorinsky based SGS prediction - using divergence of flux
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
subroutine smag_source_term_div(omega,psi,sgs,nx,ny,dx,dy)

!$ use omp_lib
implicit none
integer, intent(in) :: nx,ny
integer :: i,j
double precision,dimension(0:nx-1,0:ny-1),intent(in) :: omega, psi
double precision, intent(in) :: dx, dy
double precision,dimension(0:nx-1,0:ny-1),intent(out) :: sgs
double precision, dimension(:,:), allocatable :: omega_new, psi_new
double precision, dimension(:,:), allocatable :: dsdy, d2sdx, d2sdy, d2sdxdy
double precision :: c_turb, lap_val
double precision,dimension(:,:),allocatable :: temp_term, sgs_div

allocate(omega_new(-1:nx,-1:ny))
allocate(psi_new(-1:nx,-1:ny))

do j = 0,ny-1
  do i = 0,nx-1
    omega_new(i,j) = (omega(i,j))
    psi_new(i,j) = psi(i,j)
  end do
end do

!BC Update
call periodic_bc_update(nx,ny,omega_new)
call periodic_bc_update(nx,ny,psi_new)

!Calculating Smag Turbulence model invariants
allocate(d2sdxdy(0:nx-1,0:ny-1))
allocate(d2sdx(0:nx-1,0:ny-1))
allocate(d2sdy(0:nx-1,0:ny-1))
allocate(dsdy(-1:nx,-1:ny))

do j = 0,ny-1
  do i = 0,nx-1
    dsdy(i,j) = (psi_new(i,j+1)-psi_new(i,j-1))/(2.0d0*dy)
    d2sdx(i,j) = (psi_new(i+1,j)+psi_new(i-1,j)-2.0d0*psi_new(i,j))/(dx*dx)
    d2sdy(i,j) = (psi_new(i,j+1)+psi_new(i,j-1)-2.0d0*psi_new(i,j))/(dy*dy)
  end do
end do

call periodic_bc_update(nx,ny,dsdy)!Need for a second derivative for this quantity

do j = 0,ny-1
  do i = 0,nx-1
    d2sdxdy(i,j) = (dsdy(i+1,j)-dsdy(i-1,j))/(2.0d0*dx)
  end do
end do

!Smag invariant
do j = 0,ny-1
  do i = 0,nx-1
    sgs(i,j) = dsqrt(4.0d0*d2sdxdy(i,j)**2 + (d2sdx(i,j)-d2sdy(i,j))**2)!This is solely EV kernel now
  end do
end do

c_turb = 1.0d0*dx

!Calculating divergence term
allocate(temp_term(-1:nx,-1:ny))
allocate(sgs_div(0:nx-1,0:ny-1))
do j = 0,ny-1
do i = 0,nx-1
	temp_term(i,j) = sgs(i,j)
end do
end do

call periodic_bc_update(nx,ny,temp_term)
!calculating gradient of ev kernel
do j = 0,ny-1
do i = 0,nx-1
	sgs_div(i,j) = (temp_term(i+1,j)-temp_term(i-1,j))*(omega_new(i+1,j)-omega_new(i-1,j))/(4.0*dx*dx)
	sgs_div(i,j) = sgs_div(i,j) + (temp_term(i,j+1)-temp_term(i,j-1))*(omega_new(i,j+1)-omega_new(i,j-1))/(4.0*dy*dy)
	sgs_div(i,j) = c_turb*c_turb*sgs_div(i,j)
end do
end do


deallocate(d2sdxdy,dsdy,d2sdx,d2sdy)

!Combining invariant with turbulence coefficients and laplacian
do j = 0,ny-1
  do i = 0,nx-1
    lap_val = (omega_new(i+1,j)+omega_new(i-1,j)-2.0d0*omega_new(i,j))/(dx*dx)
    lap_val = lap_val + (omega_new(i,j+1)+omega_new(i,j-1)-2.0d0*omega_new(i,j))/(dy*dy)
    sgs(i,j) = c_turb*c_turb*sgs(i,j)*lap_val + sgs_div(i,j)
  end do
end do

deallocate(omega_new,psi_new,temp_term,sgs_div)

return
end


!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!Subroutine for Leith SGS calculation
!Validated
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
subroutine leith_source_term(omega,psi,sgs,nx,ny,dx,dy)

!$ use omp_lib
implicit none
integer, intent(in) :: nx,ny
integer :: i,j
double precision,dimension(0:nx-1,0:ny-1),intent(in) :: omega, psi
double precision, intent(in) :: dx, dy
double precision,dimension(0:nx-1,0:ny-1),intent(out) :: sgs
double precision, dimension(:,:), allocatable :: omega_new, psi_new
double precision, dimension(:,:), allocatable :: dwdx,dwdy
double precision :: c_turb, lap_val


allocate(omega_new(-1:nx,-1:ny))
allocate(psi_new(-1:nx,-1:ny))

do j = 0,ny-1
  do i = 0,nx-1
    omega_new(i,j) = (omega(i,j))
    psi_new(i,j) = psi(i,j)
  end do
end do

!BC Update
call periodic_bc_update(nx,ny,omega_new)
call periodic_bc_update(nx,ny,psi_new)

!Calculating Leith turbulence model invariants
allocate(dwdy(0:nx-1,0:ny-1))
allocate(dwdx(0:nx-1,0:ny-1))
do j = 0,ny-1
do i = 0,nx-1
dwdy(i,j) = (omega_new(i,j+1)-omega_new(i,j-1))/(2.0d0*dy)
dwdx(i,j) = (omega_new(i+1,j)-omega_new(i-1,j))/(2.0d0*dx)
end do
end do

!Leith invariant
do j = 0,ny-1
  do i = 0,nx-1
    sgs(i,j) = dsqrt(dwdx(i,j)**2 + dwdy(i,j)**2)
  end do
end do

deallocate(dwdx,dwdy)

!Combining invariants with turbulence coefficients and laplacian
c_turb = 1.0d0*dx
do j = 0,ny-1
  do i = 0,nx-1
    lap_val = (omega_new(i+1,j)+omega_new(i-1,j)-2.0d0*omega_new(i,j))/(dx*dx)
    lap_val = lap_val + (omega_new(i,j+1)+omega_new(i,j-1)-2.0d0*omega_new(i,j))/(dy*dy)
    sgs(i,j) = c_turb*c_turb*c_turb*sgs(i,j)*lap_val
  end do
end do

deallocate(omega_new,psi_new)

return
end


!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!Subroutine for Leith SGS calculation
!Validated
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
subroutine leith_source_term_div(omega,psi,sgs,nx,ny,dx,dy)

!$ use omp_lib
implicit none
integer, intent(in) :: nx,ny
integer :: i,j
double precision,dimension(0:nx-1,0:ny-1),intent(in) :: omega, psi
double precision, intent(in) :: dx, dy
double precision,dimension(0:nx-1,0:ny-1),intent(out) :: sgs
double precision, dimension(:,:), allocatable :: omega_new, psi_new
double precision, dimension(:,:), allocatable :: dwdx,dwdy
double precision :: c_turb, lap_val
double precision, dimension(:,:), allocatable :: sgs_div, temp_term

allocate(omega_new(-1:nx,-1:ny))
allocate(psi_new(-1:nx,-1:ny))

do j = 0,ny-1
  do i = 0,nx-1
    omega_new(i,j) = (omega(i,j))
    psi_new(i,j) = psi(i,j)
  end do
end do

!BC Update
call periodic_bc_update(nx,ny,omega_new)
call periodic_bc_update(nx,ny,psi_new)

!Calculating Leith turbulence model invariants
allocate(dwdy(0:nx-1,0:ny-1))
allocate(dwdx(0:nx-1,0:ny-1))
do j = 0,ny-1
do i = 0,nx-1
dwdy(i,j) = (omega_new(i,j+1)-omega_new(i,j-1))/(2.0d0*dy)
dwdx(i,j) = (omega_new(i+1,j)-omega_new(i-1,j))/(2.0d0*dx)
end do
end do

!Leith invariant
do j = 0,ny-1
  do i = 0,nx-1
    sgs(i,j) = dsqrt(dwdx(i,j)**2 + dwdy(i,j)**2)
  end do
end do

!---------------------------------------------------------------------------------

c_turb = 1.0d0*dx

!Calculating divergence term
allocate(temp_term(-1:nx,-1:ny))
allocate(sgs_div(0:nx-1,0:ny-1))
do j = 0,ny-1
do i = 0,nx-1
	temp_term(i,j) = sgs(i,j)
end do
end do

call periodic_bc_update(nx,ny,temp_term)
!calculating gradient of ev kernel
do j = 0,ny-1
do i = 0,nx-1
	sgs_div(i,j) = (temp_term(i+1,j)-temp_term(i-1,j))*(omega_new(i+1,j)-omega_new(i-1,j))/(4.0*dx*dx)
	sgs_div(i,j) = sgs_div(i,j) + (temp_term(i,j+1)-temp_term(i,j-1))*(omega_new(i,j+1)-omega_new(i,j-1))/(4.0*dy*dy)
	sgs_div(i,j) = c_turb*c_turb*c_turb*sgs_div(i,j)
end do
end do

!Combining invariants with turbulence coefficients and laplacian
do j = 0,ny-1
  do i = 0,nx-1
    lap_val = (omega_new(i+1,j)+omega_new(i-1,j)-2.0d0*omega_new(i,j))/(dx*dx)
    lap_val = lap_val + (omega_new(i,j+1)+omega_new(i,j-1)-2.0d0*omega_new(i,j))/(dy*dy)
    sgs(i,j) = c_turb*c_turb*c_turb*sgs(i,j)*lap_val + sgs_div(i,j)
  end do
end do

!---------------------------------------------------------------------------------

deallocate(dwdx,dwdy)
deallocate(omega_new,psi_new,temp_term,sgs_div)

return
end



!-----------------------------------------------------------------------------------------!
!Subroutine for nearest neighbor ML model classification
!-----------------------------------------------------------------------------------------!
subroutine sgs_classify(sgs,sgs_ml,sgs_smag,sgs_leith,sfs_ad,nx,ny)
implicit none

integer, intent(in) :: nx, ny
integer :: i, j
double precision, intent(in), dimension(0:nx-1,0:ny-1) :: sgs_ml, sgs_smag, sgs_leith,sfs_ad
double precision, intent(inout), dimension(0:nx-1,0:ny-1) :: sgs
double precision :: diff_sma, diff_lei, diff_ad


do j=0,ny-1
  do i = 0,nx-1
    
      diff_sma = dabs(sgs_ml(i,j)-sgs_smag(i,j))
      diff_lei = dabs(sgs_ml(i,j)-sgs_leith(i,j))
      diff_ad = dabs(sgs_ml(i,j)-sfs_ad(i,j))

      if (diff_sma < diff_lei .and. diff_sma < diff_ad) then
        sgs(i,j) = sgs_smag(i,j)
      else if (diff_lei < diff_sma .and. diff_lei < diff_ad) then
        sgs(i,j) = sgs_leith(i,j)
      else
        sgs(i,j) = sfs_ad(i,j)
      end if

  end do
end do


return
end
