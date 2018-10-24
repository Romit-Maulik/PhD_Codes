!-----------------------------------------------------------------------------------!
!Filtering (Gaussian) - using definition in classical image processing literature
!Same definition as in Approximate_Deconvolution.f95
!Validated
!-----------------------------------------------------------------------------------!
subroutine filter_gaussian(q_org,nx,ny,sigma)
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

call periodic_bc_update(nx,ny,u)

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

call periodic_bc_update(nx,ny,v)

!filter in x

do j=0,ny-1
do i=0,nx-1


sumval = 0.0d0
do k = -3,3
sumval = sumval + kernel(k)*v(i,j+k)
end do
q(i,j) = sumval



                           
end do
end do
 
call periodic_bc_update(nx,ny,q)

do j = 0,ny-1
     do i = 0,nx-1
          q_org(i,j) = q(i,j)
     end do
end do

deallocate(u,v,q)

return
end


!---------------------------------------------------------------------------!
!subroutine - boundary condition update
!Validated
!---------------------------------------------------------------------------!
subroutine periodic_bc_update(nx,ny,u)
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
