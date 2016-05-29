program testing
implicit none
real :: v1(1000)
integer :: i, f(4,3), j
logical :: l(1000)
character(len=1024) :: filename

call random_number(v1)
v1 = v1*10
l = v1>5.0
write(*,*) count(l)
do i = 1,1000
	if (l(i)) then
		v1(i) = 10-v1(i)
	end if
end do
l = v1>5.0
write(*,*) count(l)

do j = 1,3
do i = 1,4
f(i,j) = i+j-1
end do
write (filename, "(A5,I0.5)") "hello",j*10000

    print *, trim(filename)//'.dat'
 end do
 write(*,*) f
 write(*,*) 'trying to get array f row-wise'
 write(*,*) f(2,:)
 

    

    


end program testing
