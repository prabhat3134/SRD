module SRD_lib

implicit none

INTEGER, PARAMETER :: dp = selected_real_kind(15, 307)
REAL(kind=dp), PARAMETER :: pi=4.D0*DATAN(1.D0), e = 2.71828

 contains
!*******************************************
!	To generate random seed
subroutine init_random_seed()

      INTEGER :: i, n, clock
      INTEGER, DIMENSION(:), ALLOCATABLE :: seed

      CALL RANDOM_SEED(size = n)
      ALLOCATE(seed(n))

      CALL SYSTEM_CLOCK(COUNT=clock)

      seed = clock + 37 * (/ (i - 1, i = 1, n) /)
      CALL RANDOM_SEED(PUT = seed)

      DEALLOCATE(seed)
end subroutine init_random_seed


! For eliminating particle within a box
subroutine box_eliminate(rx,ry,rx1,ry1,l0,b0,width)
implicit none

real(kind=dp), dimension(:) :: rx,ry
integer ::  l0, b0, width
logical :: l1(size(rx)), l2(size(rx))
real(kind=dp), dimension(:) , allocatable::rx1,ry1

l1 = ry<b0+width .and. ry>b0 .and. rx<l0+width .and. rx>l0
l2 = .not. l1
rx1 = pack(rx,l2)
ry1 = pack(ry,l2)

end subroutine box_eliminate

!******************************************
! Getting normal distribution from uniform distribution

subroutine random_normal(vel,np,av,std, v1, v2)
implicit none
integer  :: np, i
real(kind=dp) :: vel(np), v1(np), v2(np)	! v1,v2 are two sets of random number
real :: av,std

!v1(1) = ran(-2)
!v2(1) = ran()

!DO i=1,2*np
!v1(i) = ran()
!v2(i) = ran()
!END DO

! Box Muller Transformation
vel = std*sqrt(-2*log(v1))*cos(2*pi*v2)+av

end subroutine random_normal

!******************************************
! getting weibull distribution for velocity
subroutine random_weibull(vel,np,kbt)
implicit none
real(kind=dp) :: vel(np)
integer :: np,i
real :: a,b,v,Pv, vmax, vprob, Pmax, kbt, P
i=1
a = sqrt(2*kbt)
b = 2.0
vprob = sqrt(kbt)		! most probable speed
vmax = 5*vprob			! maximum speed possible
Pmax = 1/sqrt(kbt*e)		!maximum probability

do while (i<=np)
	
	call random_number(v)
	v = v*vmax
	Pv = ( b/(a**b) ) * v**(b-1) * EXP( -1*(v/a)**b )
	call random_number(P)
	P = P*Pmax
	
	if (Pv>P) then
	  vel(i)=v
	  i=i+1
	end if	
end do
end subroutine random_weibull

!******************************************
! Intialize the domain
subroutine initialize(rx,ry,vx,vy, Ly, Lx, np, av, std, randposx, randposy, v1, v2, v3, v4)
implicit none
INTEGER :: np, Ly, Lx, i
REAL :: av, std
REAL(kind=dp) :: rx(np), ry(np), vx(np), vy(np), randposx(np), randposy(np), v1(np), v2(np), v3(np), v4(np) 

rx = randposx*Lx
ry = randposy*Ly

call random_normal(vx,np,av,std,v1,v2)
call random_normal(vy,np,av,std,v3,v4)

end subroutine initialize

!*******************************************
! Cell partition for dividing particles into grid cells
subroutine partition(rx,ry,head,list,Lx,Ly,np)    	! assumed no grid shifting
implicit none
real(kind=dp), dimension(:) :: rx, ry
integer, dimension(:) :: list
integer, dimension(:,:)    ::   head
integer :: np, Lx, Ly, xindex, yindex,i
head = 0
list = 0

do i=1,np
	yindex = ceiling(ry(i))	
	xindex = ceiling(rx(i))		
	if (yindex==0) yindex=1 
	if (xindex==0) xindex=1
	! linked list algorithm
	if (head(yindex,xindex)==0) then
		head(yindex,xindex)=i
	else 
		list(i) = head(yindex,xindex)
		head(yindex,xindex)=i
	end if	
end do
end subroutine partition

!*******************************************
! Collision Scheme
subroutine collision(vx, vy, tempx, tempy, head, list, Lx, Ly, alpha, t, randcol)
implicit none
real(kind=dp), dimension(:) :: vx, vy
integer, dimension(:) :: list
integer, dimension(:,:)    ::   head
integer :: i,j,count1, ipar, jpar, Lx, Ly, out_unit, k, t
real	:: r, alpha, vx_temp, vy_temp, alpha1
real(kind=dp) :: vxcom(Ly,Lx), vycom(Ly,Lx), tempx(Ly,Lx), tempy(Ly,Lx), randcol(Lx*Ly)

k=0
out_unit=20
vxcom=0.0
vycom=0.0
tempx=0.0
tempy=0.0

!call init_random_seed()
 !open (unit=out_unit,file="coll_rand.dat",action="read",status="old")
!			do i=1,100
!				read(out_unit,'(f15.12)'), r(i)
!			end do
 !close(out_unit) 

 !open (unit=out_unit,file="../Codes/vxcom.dat",action="write",status="replace")
do i=1,Ly	
	do j=1,Lx
		count1 = 0
		k = k+1		
		r = randcol(k)
		
		if (r >=0.5) then
			alpha1 = -alpha
		else 
			alpha1 = alpha
		end if			
		
		ipar = head(i,j)					
		do while (ipar/=0)
			!if (t==209 .and. i==1 .and. j==6) write(*,*), ipar, vx(ipar), vx(197)
			vxcom(i,j) = vxcom(i,j) + vx(ipar)
			vycom(i,j) = vycom(i,j) + vy(ipar)
			count1 = count1 + 1
			ipar = list(ipar)
		end do
		if (count1/=0) then
			vxcom(i,j) = vxcom(i,j)/count1
			vycom(i,j) = vycom(i,j)/count1
		end if
		
		! Rotation of velocity
		ipar = head(i,j)
		do while (ipar/=0)
			vx_temp  = vx(ipar)-vxcom(i,j) 
			vy_temp  = vy(ipar)-vycom(i,j) 			
			vx(ipar) = vxcom(i,j) + cos(alpha1)*vx_temp + sin(alpha1)*vy_temp
			vy(ipar) = vycom(i,j) - sin(alpha1)*vx_temp + cos(alpha1)*vy_temp			
			ipar = list(ipar)
			tempx(i,j) = tempx(i,j) + ((vx(ipar)-vxcom(i,j))**2)/2 
			tempy(i,j) = tempy(i,j) + ((vy(ipar)-vycom(i,j))**2)/2 
		end do		
		tempx(i,j) = tempx(i,j)/(count1-1)
		tempy(i,j) = tempy(i,j)/(count1-1)		
	end do
	!write(out_unit,*), vxcom(i,:)
end do		
	!close(out_unit)
end subroutine collision

!********************************************
! Boundary wall boundary conditions
subroutine thermal_bc(rx, ry, rx1, ry1, vx, vy, Ly, dt_c, g, kbt)
implicit none
real(kind=dp), dimension(:) :: rx, ry, rx1, ry1, vx, vy
real :: g, dt_c, kbt, av=0.0, std, ydiff, t_app, t_dec, vxwall
integer :: Ly, i, np, np1
logical :: Wall (size(ry)), NWall(size(ry))
real(kind=dp), dimension(:) , allocatable :: rx2, ry2, vx1, vy1, rx0, ry0, vx0, vy0, rx3, ry3, vxthermal, vythermal

std   = sqrt(kbt)
Wall  = (ry1 >(Ly*1.0)) .or. (ry1<0.0)
NWall = .not. Wall
np    = count(Wall)
np1   = count(NWall)
write(*,*) np,np1
allocate(vxthermal(np), vythermal(np))

! Unchanged particles still inside domain
rx2 = pack(rx1, NWall)
ry2 = pack(ry1, Nwall)
vx1 = pack(vx, NWall)
vy1 = pack(vy, NWall)

! Changed particles crossing boundaries
rx0 = pack(rx, Wall)   ! old particle positions
ry0 = pack(ry, Wall)
vx0 = pack(vx, Wall)
vy0 = pack(vy, Wall)
rx3 = pack(rx1, Wall)   ! New particle positions
ry3 = pack(ry1, Wall)

!call random_normal(vxthermal,np,av,std)
!call random_weibull(vythermal,np,kbt)
do i=1,np
	ydiff = ry3(i)-ry0(i)
	
	if (ydiff<0) then
	      t_app = ry0(i)/abs(vy0(i))
	      t_dec = dt_c - t_app
	      vy0(i) = vythermal(i)
	      ry3(i) = 0.0+vy0(i)*t_dec	
	else
	      t_app = (Ly-ry0(i))/abs(vy0(i))
	      t_dec = dt_c - t_app
	      vy0(i) = -vythermal(i)
	      ry3(i) = Ly+ vy0(i)*t_dec	
	end if
	
	vx0(i)  = vx0(i)-g*(dt_c-t_app)
	rx3(i)  = rx0(i) + vx0(i)*t_app      
      vx0(i)  = vxthermal(i) + g*t_dec;
      rx3(i)  = rx1(i) + vx0(i)*t_dec
end do

rx1(1:np1) = rx2
rx1(np1+1:np+np1) = rx3
ry1(1:np1) = ry2
ry1(np1+1:np+np1) = ry3
vx(1:np1) = vx1
vx(np1+1:np+np1) = vx0
vy(1:np1) = vy1
vy(np1+1:np+np1) = vy0

Wall  = (ry1 >Ly*1.0) .or. (ry1<0.0)
NWall = .not. Wall
np    = count(Wall)
np1   = count(NWall)
write(*,*) np,np1
 
deallocate(rx2, ry2, vx1, vy1, rx0, ry0, vx0, vy0, rx3, ry3, vxthermal, vythermal)
end subroutine thermal_bc

!****************************************************************
! Another subroutine for thermal boundary conditions
subroutine thermal_wall(rx, ry, rx1, ry1, vx, vy, Ly, dt_c, g, kbt)
implicit none
integer, parameter :: np1=100		! Expected no. of particles crossing boundary
real(kind=dp), dimension(:) :: rx, ry, rx1, ry1, vx, vy
real :: g, dt_c, kbt, av=0.0, std, t_app, t_dec
integer :: Ly, i, np, j=0
logical :: Wall(size(ry))
real(kind=dp) :: vxwall(np1), vywall(np1)
std   = sqrt(kbt)
np = size(rx)

Wall = (ry1 > (Ly*1.0) .or. ry1 < 0.0)
!call random_normal(vxwall,np1,av,std)
!call random_weibull(vywall,np1,kbt)
do i = 1,np
	if (Wall(i)) then	
		if ((ry1(i)-ry(i))<0.0) then
			j = j+1
			t_app = ry(i)/abs(vy(i))
	      	t_dec = dt_c - t_app
	      	vy(i) = vywall(j)
	      	ry1(i) = 0.0 + vy(i)*t_dec	
		else 
			j = j+1
			t_app = (Ly-ry(i))/abs(vy(i))
	      	t_dec = dt_c - t_app
	      	vy(i) = -vywall(j)
	      	ry1(i) = Ly + vy(i)*t_dec	
		end if		
		vx(i)  = vx(i)-g*(dt_c-t_app)
		rx1(i)  = rx(i) + vx(i)*t_app      
      	vx(i)  = vxwall(j) + g*t_dec
      	rx1(i)  = rx1(i) + vx(i)*t_dec       	       				
	end if
	if (j==np1) then		! if more particles are crossing, generate new numbers
		!call random_normal(vxwall,np1,av,std)
		!call random_weibull(vywall,np1,kbt)
		j=0	
	end if
end do

end subroutine thermal_wall

!***********************************************************************
! Writing data to files for analysis
subroutine writing(rx,ry,vx,vy,head,list, density, Lx ,Ly)
implicit none
real(kind=dp), dimension(:) :: rx, ry, vx, vy
integer , dimension(:,:)    :: head
integer, dimension (:)      :: list
integer :: Lx, Ly, i, j, den=0, ipar
real(kind=dp) :: density(Ly,Lx), vx_com(Ly,Lx), vy_com(Ly,Lx)
density(Ly,Lx)=0.0
vx_com(Ly,Lx)=0.0
vy_com(Ly,Lx)=0.0

do j=1, Lx
	do i=1, Ly
		ipar = head(i,j)		
		do while (ipar/=0)
			vx_com(i,j) = vx_com(i,j) + vx(ipar)
			vy_com(i,j) = vy_com(i,j) + vy(ipar)
			den = den +1
			ipar = list(ipar)
		end do
		if (den/=0) then		
		density(i,j) = den
		vx_com(i,j)  = vx_com(i,j)/den	
		vy_com(i,j)  = vy_com(i,j)/den
		end if	
	end do
end do

end subroutine writing
!***********************************************************************
! Subroutine for bounce back boundary condition at wall
subroutine bounce_wall(rx, ry, rx1, ry1, vx, vy, Ly, dt_c, g, kbt,np)
implicit none
real(kind=dp), dimension(np)  :: rx, ry
real(kind=dp), dimension(np)  :: rx1, ry1, vx, vy
real :: g, dt_c, kbt, t_app, t_dec
integer :: Ly, i, np
logical :: Wall(size(ry))

Wall = (ry1 > (Ly*1.0) .or. ry1 < 0.0)  			  		

do i = 1,np
	if (Wall(i)) then			
		if ((ry1(i)-ry(i))<0.0) then			
			t_app = ry(i)/abs(vy(i))
	      	t_dec = dt_c - t_app
	      	vy(i) = -vy(i)
	      	ry1(i) = 0.0 + vy(i)*t_dec	
		else 			
			t_app = (Ly-ry(i))/abs(vy(i))
	      	t_dec = dt_c - t_app
	      	vy(i) = -vy(i)
	      	ry1(i) = Ly + vy(i)*t_dec	
		end if		
		vx(i)  = vx(i)-g*(dt_c-t_app)
		rx1(i)  = rx(i) + vx(i)*t_app      
      	vx(i)  = -vx(i) + g*t_dec
      	rx1(i)  = rx1(i) + vx(i)*t_dec  		
	end if	
end do

end subroutine bounce_wall



end module SRD_lib
