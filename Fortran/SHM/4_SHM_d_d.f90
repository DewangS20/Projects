program SHM_d
implicit none

real, allocatable, dimension(:) :: a, v, x, t, ke, u, e
real :: dt, k, m, pi, period, b, fd, wd
integer :: i, n, cycles

open(unit = 100, file = 'tvx_ddshm.dat')
open(unit = 200, file = 'tvv_ddshm.dat')
open(unit = 300, file = 'tva_ddshm.dat')
open(unit = 400, file = 'tvke_ddshm.dat')
open(unit = 500, file = 'tvu_ddshm.dat')
open(unit = 600, file = 'tve_ddshm.dat')

pi = 4.0 * atan(1.0)
k = 1.0 !N/m
m = 1.0 !kg
period = 2 * pi * sqrt(m/k) !s

write(*, *) 'dampening value'
read(*, *) b

write(*, *) 'how many cycles'
read(*, *) cycles
dt = 0.001 !use a smaller time interval to get a straighter line on total energy
n = floor(cycles * period / dt) !number of points on the array
allocate(x(n), a(n), u(n))
allocate(v(n-1), t(n -1), ke(n-1), e(n-1))

!initial conditions
t(1) = 0.0
x(1) = 5.0
v(1) = 0.0
a(1) = ((-1 * k * x(1)) - (b * v(1)) + (fd * cos(wd * t(1))))/ m
fd = 1.0
write(*, *) 'driven frequency?'
read(*, *) wd

ke(1) = 0.5 * m * (v(1) ** 2)
u(1) = 0.5 * k * (x(1) ** 2)
e(1) = ke(1) + u(1)

x(2) = x(1) + (v(1) * dt) + (0.5 * a(1) * (dt ** 2))
u(2) = 0.5 * k * x(2) * x(2)

!i loop
do i = 2, n - 1
	t(i) = t(i - 1) + dt
	a(i) = ((-1 * k * x(i-1)) - (b * v(i - 1)) + (fd * cos(wd * t(i))))/ m
	x(i + 1) = 2 * x(i) - x(i - 1) + ((dt ** 2) * a(i))
	v(i) = (x(i + 1) - x(i - 1)) / (2 * dt)
	ke(i) = 0.5 * m * (v(i) ** 2)
	u(i + 1) = 0.5 * k * (x(i + 1) ** 2)
	e(i) = ke(i) + u(i)
enddo

do i = 1, n-1
	
	write(100, *) t(i), x(i)
	write(200, *) t(i), v(i) 
	write(300, *) t(i), a(i)
	write(400, *) t(i), ke(i)
	write(500, *) t(i), u(i)
	write(600, *) t(i), e(i)

enddo

end program SHM_d
