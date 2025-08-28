program lorentz
implicit none

real(8), dimension(6) :: s, temp, k1, k2, k3, k4
real(8) :: t, dt, pi, q, m, a(3), k, u, ef(3), b(3)
integer :: i
!s(1:3) position s(4:6) velocity
!r is magnitude 

open(unit = 100, file = '7_position.dat')
!open(unit = 200, file = '6_velocity.dat')
!open(unit = 300, file = '6_acceleration.dat')
!open(unit = 150, file = '6_helix.dat')
!open(unit = 400, file = '6_ke.dat')
!open(unit = 500, file = '6_u.dat')
!open(unit = 600, file = '6_te.dat')



!parameters
pi = 4.0d0 * atan(1.0d0)
dt = 0.0001 

m = 1 
ef = 0
ef(2) = 0.1
ef(3) = (0.1)
b = 0.0
b(3) = 0.1
q = 1.0d0

!initial values
s = 0.0d0
s(4) = 1.0d0
t = 0.0d0

do i = 1, 100000, 1
	!rum kutta 
	a(1) = q * (ef(1) + (s(5) * b(3) - s(6) * b(2))) / m
	a(2) = q * (ef(2) + (s(6) * b(1) - s(4) * b(3))) / m
	a(3) = q * (ef(3) + (s(4) * b(2) - s(5) * b(1))) / m
	k1(1:3) = s(4:6)
	k1(4:6) = a(1:3)
	temp = s + 0.5 * dt * k1
	a(1) = q * (ef(1) + (temp(5) * b(3) - temp(6) * b(2))) / m
	a(2) = q * (ef(2) + (temp(6) * b(1) - temp(4) * b(3))) / m
	a(3) = q * (ef(3) + (temp(4) * b(2) - temp(5) * b(1))) / m
	k2(1:3) = temp(4:6)
	k2(4:6) = a(1:3)
	temp = s + 0.5 * dt * k2
	a(1) = q * (ef(1) + (temp(5) * b(3) - temp(6) * b(2))) / m
	a(2) = q * (ef(2) + (temp(6) * b(1) - temp(4) * b(3))) / m
	a(3) = q * (ef(3) + (temp(4) * b(2) - temp(5) * b(1))) / m
	k3(1:3) = temp(4:6)
	k3(4:6) = a(1:3)
	temp = s + dt * k3
	a(1) = q * (ef(1) + (temp(5) * b(3) - temp(6) * b(2))) / m
	a(2) = q * (ef(2) + (temp(6) * b(1) - temp(4) * b(3))) / m
	a(3) = q * (ef(3) + (temp(4) * b(2) - temp(5) * b(1))) / m
	k4(1:3) = temp(4:6)
	k4(4:6) = a(1:3)
	s = s + dt * (k1 + 2 * k2 + 2 * k3 + k4) / 6
	
	
	t = t + dt
	
	write(100, *) s(1:3)
end do
end program lorentz

