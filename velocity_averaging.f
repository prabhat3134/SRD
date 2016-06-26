%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Velocity Averaging %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

! for averaging the values 
    			if (t>150000) then 
    				t_count = t_count+1
    				vxcom   = 0.0
    				vycom   = 0.0
    			 	do  i=1,Ly
    			  		p_count = 0
    			  		do j = 1,Lx    			  			
    			  			ipar = head(i,j)    			  		
    			  			do while (ipar/=0)
    			  				vxcom(i) = vxcom(i) + vx(ipar)
    			  				vycom(i) = vycom(i) + vy(ipar)
    			  				p_count  = p_count +1
    			  				ipar     = list(ipar)    			  		
    			  			end do    			  			
    			  		end do
    			  		vxcom(i) = vxcom(i)/p_count
    			  		vycom(i) = vycom(i)/p_count
    			  	end do 
    			  	vx1 = vx1 + vxcom
    			  	vy1 = vy1 + vycom
    			  	if (mod(t_count,10000)==0) then
    			  		vx1_temp = vx1/t_count
    			  		vy1_temp = vy1/t_count    			  		
    			  		write (fname, "(A16,I0.6)") "./data/BB_vel_no_MBS",t                             
    			  		open (unit=out_unit,file=trim(fname)//'.dat',action="write",status="replace")
    			  		do i = 1,Ly
    			  			write (out_unit,'(I2.2,F10.5, F10.5)') i, vx1_temp(i), vy1_temp(i)			  		
    			  		end do
    			  		close(out_unit)					
    			  	end if    			  	  		  
    			end if		

		close(out_unit)
	! For averaged velocity profile
  		vx1_temp = vx1/t_count
    	        vy1_temp = vy1/t_count    			  		
 		write (fname, "(A26,I0.6)") "../Codes/thermal_vel_poise",t                             
    		open (unit=out_unit,file=trim(fname)//'.dat',action="write",status="replace")
    		do i = 1,Ly
    			write (out_unit,'(I2.2,F10.5, F10.5)') i, vx1_temp(i), vy1_temp(i)			  		
    		end do
    		close(out_unit)	
 		 
