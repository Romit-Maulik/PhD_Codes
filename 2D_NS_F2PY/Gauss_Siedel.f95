!------------------------------------------------------------------
!Poisson solver using Gauss Siedel
!Not explicitly used in deployment - only for validation
!u-field (psi), f-source (-omega)
!------------------------------------------------------------------
subroutine solve_poisson(u,f,nx,ny,dx,dy,gs_tol)
implicit none

integer, intent(in) :: nx, ny
double precision, intent(in) :: dx, dy, gs_tol
double precision, intent(in), dimension(0:nx-1,0:ny-1) :: f
double precision, intent(inout), dimension(0:nx-1,0:ny-1) :: u

double precision, dimension(:,:), allocatable :: u_update_1, u_update_2
integer :: not_converged, i, j, gsiter

allocate(u_update_1(-1:nx,-1:ny))
allocate(u_update_2(-1:nx,-1:ny))

do j = 0, ny-1
do i = 0, nx-1

u_update_1(i,j) = u(i,j)
u_update_2(i,j) = u(i,j)

end do
end do

!BC update
do j = 0,ny-1
u_update_1(-1,j) = u_update_1(nx-1,j)
u_update_1(nx,j) = u_update_1(0,j)

u_update_2(-1,j) = u_update_2(nx-1,j)
u_update_2(nx,j) = u_update_2(0,j)
end do

do i = -1,nx
u_update_1(i,-1) = u_update_1(i,ny-1)
u_update_1(i,ny) = u_update_1(i,0)

u_update_2(i,-1) = u_update_2(i,ny-1)
u_update_2(i,ny) = u_update_2(i,0)
end do

!GS Iterations
not_converged = 1
gsiter = 0
do while (not_converged == 1)

do j = 0, ny - 1
do i = 0, nx - 1

u_update_2(i, j) = -0.25*(f(i, j) * dx * dx - (u_update_2(i+1, j) + u_update_2(i-1, j) &
			+ u_update_2(i, j+1) + u_update_2(i, j-1)))

end do
end do

gsiter = gsiter + 1


!Update BC
do j = 0,ny-1
u_update_2(-1,j) = u_update_2(nx-1,j)
u_update_2(nx,j) = u_update_2(0,j)
end do

do i = -1,nx
u_update_2(i,-1) = u_update_2(i,ny-1)
u_update_2(i,ny) = u_update_2(i,0)
end do

if (maxval(u_update_2-u_update_1) < gs_tol) then

do j = 0, ny-1
do i = 0, nx-1
u(i,j) = u_update_2(i,j)
end do
end do

!print*,'Converged: ',gsiter,' Iterations of GS , Residual: ',maxval(u_update_2-u_update_1)

deallocate(u_update_1,u_update_2)
exit

else


do j = -1,ny
do i = -1,nx
u_update_1(i,j) = u_update_2(i,j)
end do
end do

end if

end do


!Done
return
end
