!##############################
! EXPANDING BOX ZERO VELOCITY (Sadwitch)
!
! the code includes matrix inversion for local growth rate determination
!The matrixes for inversion are always declared locally with current sizes.
!This to avoid problems with inversion or useless huge matrix inversion
!
!PROBLEMS
!	i) Mass conservation not so good. 
!		This because the calculation of Q_BC uses only first derivatives (cannot exploit curvature as for der_BC, 
!		I don't have the second derivative of mu_BC)
!	ii) Way integrals are computed not the most precise.. Could be improved considering for instance the average between previous and following point.. (This is a criticism valid for all versions of the code even the one at fixed box) 
!RECENT IMPROVEMENTS: 
! Curvature is time dependent and depends on derivative at the boundary and (half) pore size
!
!#############




program MAIN

use com_var
use read_input
use initialization
use write_output
use observables
use integrate
implicit none
!+++++++++++++	accessory variables	++++++++++++++++
real percent_completion
logical :: exit_stat = .false.
integer(kind = i15) :: i
integer :: fr_number
integer, parameter :: general_set_ = 0, new_conf_ = 1, write_initial_ =0 , write_along_=1,control_=2,write_finished_ =3
integer, parameter :: get_hydro_=1
double precision, dimension(:),allocatable:: h !profile and initial profile
double precision :: zeta_BC



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

allocate(h(M_max))
allocate(c_eq(M_max))


call read_input_file() !creates h0

call Init(general_set_)
!if (read_last .eqv. .false.) call Init(new_conf_)  Use configuration generated by constant size box

h = 0.0d0
h(1:M) = h0(:) !assign initial value from h0
zeta_BC = hs-h_BC


call init_int(h0)!initialize integration function

call write_main_out(write_initial_,0,h0)



!					INTEGRATE

!--------------------------------------------------

print *, 'beginning integration'




mean_conc= Av_ceq(exit_stat)!Compute average concentration
mass0 =  get_mass(h0)!compute total mass
mass = mass0
fr_number=0
!crystal_vel = hydro(h0,pressure) is in init
fluxL = get_flux()!flux at t=0
cumul_flux = 0.0d0


call write_init()! write on file initial data.
call write_hout(h0)!write height profile,force and concentration
!call write_3D(0,h0)



lateral_vel = -1.0d0/der_BC * nu*(1.0d0+sup-exp(mu_BC))

crystal_force = get_cryst_force(h0)



call write_crystal_force(h0,0_i15)
call write_main_out(write_along_,fr_number,h0)


!+++++++++++++++++++++++++	MAIN LOOP	+++++++++++++++++++++++

!--------- CAREFULL: less computational efficient but better physically to check always concentration -------


i = 2

do while (i<(N+1))

	cumul_flux = cumul_flux + fluxL*deltat

	h(1:M) =  step(h(1:M))!INTEGRATION STEP calculate also delta_BC and lateral vel

	mass = get_mass(h(1:M))
	fluxL = get_flux()
	der_BC = (h(M) - h(M-1) - &
	&(0.5d0*deltax*deltax+deltax*delta_BC)*curv_BC)/(deltax+0.5d0*deltax*deltax/L+deltax*delta_BC/L)

!	der_BC = (h(M) - h(M-1))*dx
	curv_BC = 1.0d0/(sqrt(1.0d0+der_BC**2)*(W2-zeta_BC))

	lateral_vel = -1.0d0/der_BC * nu*(1.0d0+sup-exp(mu_BC))


	delta_BC = delta_BC + lateral_vel*deltat

!##########		BOX EXPANSION	#############


if(delta_BC>deltax) then


	!print *, '---------- BOX EXPANSION ------'

	if ((M+1) > M_max) then
		print *, 'Contact size larger than simulation box. Exit'
		flush(0)
		stop 
	endif

	delta_BC = delta_BC-deltax

	h(M+1) =(h_BC-0.5d0*delta_BC*dx*(h(M-1)-4.0d0*h(M)))*deltax/(1.5d0*delta_BC+deltax)



	M = M+1
	L = (M-1)*deltax

	!print *, 'L = ', L
	!write(*,*) 'x_BC =', delta_BC + L,'h_bc =', h_bc, 'h(L) =',h(M), 'h(L-1)=',h(M-1)
	!print *, 'der_BC=', der_BC
	!print *, 'delta_BC=', delta_BC







elseif (delta_BC <0) then
	!print *, '---------- BOX CONTRACTION ------'

	delta_BC = deltax + delta_BC


	h(M-1) =(h_BC-0.5d0*delta_BC*dx*(h(M-3)-4*h(M-2)))*deltax/(1.5d0*delta_BC+deltax)

	M = M-1

	L = (M-1)*deltax


	!print *, 'L = ', L
	!print *, 'Delta_BC = ', delta_BC
	!write(*,*) 'x_BC =', delta_BC + L,'h_bc =', h_bc, 'h(L) =',h(M), 'h(L-1)=',h(M-1)


endif



!####################


	if (mod((i-1),freq_step) == 0) then!Write output every freq_fr

		mean_conc = Av_ceq(exit_stat)!also uploads c_eq
		if(exit_stat) then 
		flush(0)
		stop 'An unphysical value of the mean concentration was obtained'
		endif
		
		fr_number=fr_number+1

		crystal_force = get_cryst_force(h)

		call write_main_out(write_along_,fr_number)

		call write_hout(h(1:M))

		call write_crystal_force(h,i)

		if (isnan(h(1))) then
			print *, 'x_BC =', delta_BC + L,'h_bc =', h_bc, 'h(L) =',h(M), 'h(L-1)=',h(M-1)
			print *, 'der_BC=',(h_BC-h(M))/delta_BC
			flush(0)
			stop 'h0 is a NaN'
		endif
		percent_completion = real(i)/real(N+1) *100

		write (*,"(A,F4.1,A)") "Completed: ", percent_completion, "%"

	endif

!# comment out for faster execution -------
!		call write_3D(fr_number,h)
!# comment out for faster code -------
	! endif
		
	i = i+1

enddo
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


call write_main_out(control_,0,h(1:M))







call idate(today)   
call itime(now)
call cpu_time(stop_time)
call write_main_out(write_finished_)
call write_final(h(1:M))
deallocate(h)
deallocate(h0)
deallocate(c_eq)
deallocate(Q)
deallocate(mu)
deallocate(growth_vel)




end program
