!-----------------------------------------------------------------------------------!
!Implemented natively - Romit Maulik
!Filtering (Gaussian) - using definition in classical image processing literature
!Refer Wolfram Mathematica notebook in Kappa_SFS_Predictions for kernel development
!Validated
!-----------------------------------------------------------------------------------!
subroutine filter_gaussian(q_org,nx,ny,sigma)
!$ use omp_lib

implicit none
integer,intent(in)::nx,ny
double precision, intent(inout) :: q_org(0:nx-1,0:ny-1)
double precision, intent(in) :: sigma

integer ::i,j,k
double precision,allocatable:: q(:,:),u(:,:),v(:,:)
double precision :: sumval,pi
double precision,dimension(-3:3) :: kernel

pi = datan(1.0d0)*4.0d0

!Obtain kernel for convolution from classical Gaussian used in image processing literature
sumval = 0.0d0
do i = -3,3
kernel(i) = dexp(-(dfloat(i*i)/(2.0d0*sigma**2)))/dsqrt(2.0d0*(sigma**2)*pi)
sumval = sumval + kernel(i)
end do

!Normalize the kernel to unity
do i =-3,3
kernel(i) = kernel(i)/sumval
end do

allocate(q(-3:nx+2,-3:ny+2))
allocate(u(-3:nx+2,-3:ny+2))
allocate(v(-3:nx+2,-3:ny+2))

do j=0,ny-1
do i=0,nx-1
u(i,j) = q_org(i,j)
end do
end do

call periodic_bc_update_for_filter(nx,ny,u)

!$OMP PARALLEL DO
!filter in y
do j=0,ny-1
do i=0,nx-1

sumval = 0.0d0
do k = -3,3
sumval = sumval + kernel(k)*u(i,j+k)
end do
v(i,j) = sumval

end do
end do
!$OMP END PARALLEL DO

call periodic_bc_update_for_filter(nx,ny,v)

!filter in x
!$OMP PARALLEL DO
do j=0,ny-1
do i=0,nx-1


sumval = 0.0d0
do k = -3,3
sumval = sumval + kernel(k)*v(i,j+k)
end do
q(i,j) = sumval
              
end do
end do
!$OMP END PARALLEL DO


call periodic_bc_update_for_filter(nx,ny,q)

do j = 0,ny-1
     do i = 0,nx-1
          q_org(i,j) = q(i,j)
     end do
end do

deallocate(u,v,q)

return
end


!---------------------------------------------------------------------------!
!subroutine - boundary condition update for +-3
!Validated
!---------------------------------------------------------------------------!
subroutine periodic_bc_update_for_filter(nx,ny,u)
implicit none

integer :: nx, ny, i, j
double precision, dimension(-3:nx+2,-3:ny+2) :: u

do i = 0,nx-1
u(i,-1) = u(i,ny-1)
u(i,-2) = u(i,ny-2)
u(i,-3) = u(i,ny-3)


u(i,ny) = u(i,0)
u(i,ny+1) = u(i,1)
u(i,ny+2) = u(i,2)
end do

do j = -3,ny+2
u(-1,j) = u(nx-1,j)
u(-2,j) = u(nx-2,j)
u(-3,j) = u(nx-3,j)


u(nx,j) = u(0,j)
u(nx+1,j) = u(1,j)
u(nx+2,j) = u(2,j)

end do

end subroutine


!----------------------------------------------------------------------------!
!Subroutine for source term calculation using AD
!Validated - used for ML-AD-SFS predictions
!----------------------------------------------------------------------------!
subroutine approximate_deconvolution(omega,psi,source,nx,ny,dx,dy,sigma)
!$ use omp_lib

implicit none

integer, intent(in) :: nx, ny
double precision, intent(in) :: dx, dy, sigma
double precision, dimension(0:nx-1,0:ny-1), intent(in) :: omega, psi
double precision, dimension(0:nx-1,0:ny-1), intent(out) :: source
double precision, dimension(:,:),allocatable :: jcf,jcad,psi_ad,omega_ad
integer :: i,j

allocate(jcf(0:nx-1,0:ny-1))
allocate(jcad(0:nx-1,0:ny-1))

!Compute Jacobian of filtered variables
call jacobian_calc(omega,psi,nx,ny,jcf,dx,dy)

!AD process
allocate(psi_ad(0:nx-1,0:ny-1))
allocate(omega_ad(0:nx-1,0:ny-1))

call adm(nx,ny,psi,psi_ad, sigma)
call adm(nx,ny,omega,omega_ad, sigma)

!Compute Jacobian of deconvolved variables
call jacobian_calc(omega_ad,psi_ad,nx,ny,jcad,dx,dy)

!$OMP PARALLEL DO
do j = 0,ny-1
	do i = 0,nx-1
		source(i,j) = jcf(i,j)-jcad(i,j)
	end do
end do
!$OMP END PARALLEL DO

deallocate(jcf,jcad,omega_ad,psi_ad)

return
end



!-------------------------------------------------------------------------!
!Subroutine for calculation of Jacobian
!Validated
!-------------------------------------------------------------------------!
subroutine jacobian_calc(omega,psi,nx,ny,jc,dx,dy)
!$ use omp_lib
implicit none

integer :: nx, ny, i, j
double precision, dimension(0:nx-1,0:ny-1) :: omega, psi, jc
double precision, dimension(:,:), allocatable :: psi_new, omega_new
double precision :: jj1, jj2, jj3, dx, dy

allocate(psi_new(-1:nx,-1:ny))
allocate(omega_new(-1:nx,-1:ny))

do j = 0,ny-1
do i = 0,nx-1
psi_new(i,j) = psi(i,j)
omega_new(i,j) = omega(i,j)
end do
end do

call periodic_bc_update_for_jacobian(nx,ny,psi_new)
call periodic_bc_update_for_jacobian(nx,ny,omega_new)

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

jc(i, j) = (jj1 + jj2 + jj3)/3.0

end do
end do
!$OMP END PARALLEL DO

deallocate(psi_new,omega_new)

return
end

!-------------------------------------------------------------------------!
!Subroutine for approximate deconvolution of omega and psi variables
!Utilizes 3 iterative resubstitutions
!Validated
!-------------------------------------------------------------------------!
subroutine adm(nx,ny,uf,u_ad,sigma)
implicit none

integer :: nx, ny, k, ad_iter
double precision :: sigma
double precision,dimension(0:nx-1,0:ny-1) :: uf, u_ad
double precision, dimension(:,:),allocatable :: utemp

allocate(utemp(0:nx-1,0:ny-1))

ad_iter = 3

!Initialize as filtered variable
u_ad = uf

do k = 1,ad_iter

utemp = u_ad
call filter_gaussian(utemp,nx,ny,sigma)
u_ad = u_ad + (uf - utemp)

end do

deallocate(utemp)

return
end



!---------------------------------------------------------------------------!
!subroutine - boundary condition update
!Validated
!---------------------------------------------------------------------------!
subroutine periodic_bc_update_for_jacobian(nx,ny,u)
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
