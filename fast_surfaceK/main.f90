program MAIN

use com_var
use read_input
use initialization
!use integrate
use write_output
use observables
use integrate
implicit none
!+++++++++++++	accessory variables	++++++++++++++++

logical :: exit_stat = .false.
integer(kind = i15) :: i
integer :: fr_number
integer, parameter :: general_set_ = 0, new_conf_ = 1, write_initial_ =0 , write_along_=1,control_=2,write_finished_ =3
integer, parameter :: get_hydro_=1

!+++++++++	observables	++++++++++++++++++++

double precision, parameter :: difference = 0.000

!----COMMENTS
!
!
! -------------Read input file and initialize-------------

call idate(yesterday) 
call itime(before)
call cpu_time(start_time) 
write(*,*) '+++++++++++++++++STARTING RUN++++++++++++++++'
write ( *, 1000 )  yesterday(2), yesterday(1), yesterday(3), before
1000 format ( 'Date ', i2.2, '/', i2.2, '/', i4.4, '; time ',i2.2, ':', i2.2, ':', i2.2 )

call read_input_file()
call Init(general_set_)
if (read_last .eqv. .false.) call Init(new_conf_)  

call write_main_out(write_initial_)

!					INTEGRATE

!--------------------------------------------------

print *, 'beginning integration'

call init_int()!initialize integration function

meanc_eq= Av_ceq(exit_stat)!Compute average concentration
mass0 =  get_mass()!compute total mass
mass = mass0
!fluxL = get_flux()
fr_number=0
crystal_vel = hydro(h0,pressure)
fluxL = get_flux()!flux at t=0
cumul_flux = 0.0d0


call write_init()! write on file initial data.
call write_hout(0)!write height profile,force and concentration
call write_3D(0)
!call write_crystal_vel(0_i15)
call write_main_out(write_along_,fr_number)


!+++++++++++++++++++++++++	MAIN LOOP	+++++++++++++++++++++++

!--------- CAREFULL: less computational efficient but better physically to check always concentration -------


i = 2

do while (i<(N+1))

	fluxL = get_flux()!this is the flux at the previous step
	cumul_flux = cumul_flux + fluxL*deltat
	mass = get_mass()!this is the mass at the previous step
		
	h(:) =  step(h)!INTEGRATION STEP
	
	!if (mod((i-1),freq_vel) == 0) then
	!	call write_crystal_vel(i-1)
	!endif
	
	if (mod((i-1),freq_step) == 0) then!Write output every freq_fr
		meanc_eq = Av_ceq(exit_stat)!also uploads c_eq
		if(exit_stat) stop 'An unphysical value of the mean concentration was obtained'
		fr_number=fr_number+1
		
		call write_main_out(write_along_,fr_number)
		call write_hout(fr_number)!write output on file every Nframe
		!print *, 2.0d0*Q(M)/L, crystal_vel
		
!# comment out for faster code -------
!		call write_3D(fr_number)
!# comment out for faster code -------
		
		
	endif


	i = i+1
enddo
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


call write_main_out(control_)

11 continue


h_BC = h(M)


call write_final()


call idate(today)   
call itime(now)
call cpu_time(stop_time)
call write_main_out(write_finished_)
deallocate(h)
deallocate(h0)
deallocate(c_eq)

if(allocated(Q)) then
	deallocate(Q)
	deallocate(mu)
endif


end program
