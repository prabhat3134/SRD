%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% File Writing Modules %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

---------------------------------------------------------------------------------------------------------------------------------------------------------------
! Writing the positions and velocities of all particles
if (mod(t,250)==0 .or. t==2099) then
		write(*,*) "about to get screwed"
		write (fname, "(A14,I0.6)") "./data/position_velocity",t                             
    		
		open (unit=out_unit2,file=trim(fname)//'.dat',action="write",status="replace")
			do i = 1,np
				write (out_unit2,'(F10.5, F10.5, F10.6, F10.6)') rx(i), ry(i), vx(i), vy(i)	
			end do
    			  		
			
    			  		
		close(out_unit2)	
end if

---------------------------------------------------------------------------------------------------------------------------------------------------------------
! Writing the temperature of all the cells

if (mod(t,250)==0 .or. t==2099) then
		write(*,*) "about to get screwed"
		write (fname, "(A14,I0.6)") "./data/position_velocity",t                             
    		
		open (unit=out_unit2,file=trim(fname)//'.dat',action="write",status="replace")
		do i = 1,Ly
    			  write (out_unit,*) (temp(i,j), j=1,Lx)	
    			  write (out_unit,'(F9.6)') vx(i)	  		
    		end do
    			  			  		
		close(out_unit2)	
end if

----------------------------------------------------------------------------------------------------------------------------------------------------------------
! Writing the difference between KE and momentum before and after collision

			if (mod(t,1000)==0) then
				
				momx_diff = momx_after_col - momx_before_col
				momy_diff = momy_after_col - momy_before_col
				ke_diff = ke_after_col - ke_before_col
																
				counter = counter + 1;
			
				write (fname, "(A14,I0.6)") "./data/momx_diff_t",t                             
    				open (unit=out_unit2,file=trim(fname)//'.dat',action="write",status="replace")
				do fc_i = 1,Ly
					write (out_unit2,*) (momx_diff(fc_i,fc_j), fc_j=1,Lx)	
				end do
    			  	close(out_unit2)	
					
				write (fname, "(A14,I0.6)") "./data/momy_diff_t",t                             
    			  	open (unit=out_unit2,file=trim(fname)//'.dat',action="write",status="replace")
				do fc_i = 1,Ly
					write (out_unit2,*) (momy_diff(fc_i,fc_j), fc_j=1,Lx)	
				end do
    				close(out_unit2)	
					
				write (fname, "(A14,I0.6)") "./data/ke_diff_t",t                             
    				open (unit=out_unit2,file=trim(fname)//'.dat',action="write",status="replace")
				do fc_i = 1,Ly
					write (out_unit2,*) (ke_diff(fc_i,fc_j), fc_j=1,Lx)	
				end do
    				close(out_unit2)	
					
			end if   		

----------------------------------------------------------------------------------------------------------------------------------------------------------------

