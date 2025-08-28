program heat
implicit none

integer(4) :: n_x_steps, i, j, n_t_steps, plot_skips
real(8) :: dt, L, h, kappa, coeff
real(8), allocatable, dimension(:, :) :: temp

!read initial conditions




!initial conditions

kappa = 1.0d0 !diffusion constant
L = 1.0d0 !system goes from -L/2 to L/2
n_x_steps = 61 !last step index
h = L/(n_x_steps - 1) !spatial grid size
dt = 1.0d-4 !h**2/2.0/kappa-1.0d-6
n_t_steps = 1000
plot_skips =  10
coeff = kappa * dtx / h ** 2
allocate(temp(n_t_steps,n_x_steps))

!stability check
if (coeff < 0.5) then
	write(*, *) 'solution is expected to be stable'
else
	write(*, *) 'warning: solution is not expected to be stable'
endif

! set initial and boundary conditions
temp = 0.0d0 !initialize temp to be zero at all points

temp(1, ceiling(n_x_steps/2.0)) = 1.0d0 / h !initial temperature in center

!boundary conditions are t(:,1) = t(:,N) = 0
!end points are unchanged

!loop over desired number of time steps
do i = 1, (n_t_steps - 1)
	!loop over spatial grid point
	do j = 2, (n_x_steps - 1)
		temp(i+1, j) = temp(i, j) + coeff * (temp(i, j + 1) + temp(i, j - 1) - 2.0d0 * temp(i, j))
		!temp (i+1, j) = temp(i+1, j) + 0.001 * temp(i, j) !creation term (neutron diffusion)
	enddo
	!temp(i + 1, 1) = temp(i + 1, 2) !Neumann conditions
	!temp(i + 1, n_x_steps) = temp(i+1, n_x_steps - 1) !Neumann Conditions
enddo

open(unit = 100, file = 'temp.dat')

do i = 10, n_t_steps, plot_skips
	write(100, *) temp(i, :)
enddo

end program heat		