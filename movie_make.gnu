#! gnuplot

set terminal jpeg
set pm3d map
set cbrange [0:40]
do for [i=1:200]{

	set output 'image/RP_rho_'.i.'.jpg'
	splot './RP_rho_'.i.'.dat'  every ::1650::1850 matrix notitle 
}



