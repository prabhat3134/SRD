PROGRAM RP
implicit none
INTEGER, PARAMETER :: Lx = 2000, Ly = 100, NC =3, f_n = 200 
INTEGER :: i, j, k, np
REAL*8 :: avg(Lx,f_n,NC), temp(Ly, Lx)
character (len=100) :: fname, fmt_spec, variable

variable = '/RP_vx_vy_'
write(fmt_spec, '(a, I0, a)') '(', Lx, 'F10.5)'

avg = 0.0d0
DO j = 1,f_n
	temp = 0.0d0
	write(fname,'(a,I0,a)') 'results/RP_rho_',j,'.dat'
	open(unit=20, file = trim(fname), action = 'read')
	DO k = 1,Ly
		read(20,fmt = trim(fmt_spec)) temp(k,:)
	END DO
	close(20)
	avg(:,j,1) = sum( temp, 1)/( Ly*1.0d0)
end do
DO j = 1,f_n
	temp = 0.0d0
	write(fname,'(a,I0,a)') 'results/RP_temp_',j,'.dat'
	open(unit=20, file = trim(fname), action = 'read')
	DO k = 1,Ly
		read(20,fmt = trim(fmt_spec)) temp(k,:)
	END DO
	close(20)
	avg(:,j,2) = sum( temp, 1)/( Ly*1.0d0)
end do
DO j = 1,f_n
	temp = 0.0d0
	write(fname,'(a,I0,a)') 'results/RP_vx_vy_',j,'.dat'
	open(unit=20, file = trim(fname), action = 'read')
	DO k = 1,Ly
		read(20,fmt = trim(fmt_spec)) temp(k,:)
	END DO
	close(20)
	avg(:,j,3) = sum( temp, 1)/( Ly*1.0d0)
end do

write(fmt_spec, '(a, I0, a)') '(', f_n, 'F10.5)'

open(unit=20, file = './avg_rho.dat', action = 'write', status = 'replace')
DO k = 1,Lx
	write(20,fmt = trim(fmt_spec)) avg(k,:,1)
END DO
close(20)

open(unit=20, file = './avg_temp.dat', action = 'write', status = 'replace')
DO k = 1,Lx
	write(20,fmt = trim(fmt_spec)) avg(k,:,2)
END DO
close(20)

open(unit=20, file = './avg_vel.dat', action = 'write', status = 'replace')
DO k = 1,Lx
	write(20,fmt = trim(fmt_spec)) avg(k,:,3)
END DO
close(20)

END PROGRAM RP
