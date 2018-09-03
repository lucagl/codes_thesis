!-------------------------------------------------------------------------------------------
!
!
!				CONTAINS ROUTINES COMPUTING OBSERVABLES TO BE CHECKED DURING INTEGRATION
!
!
!
!-------------------------------------------------------------------------------------------
MODULE observables

use com_var
use integrate, only : dx
implicit none

double precision new_2dterm

CONTAINS
	
	function Av_ceq(exit_stat) result(meanC)

		implicit none

		integer :: j
		double precision :: meanC
		logical ,intent (inout) ::	exit_stat
		
		meanC= 0.0d0
		
		do j= 1,M
			c_eq(j) = c0*exp(mu(j))
			meanC = meanC + c_eq(j)+ nui*growth_vel(j)
		enddo


		meanC = (meanC+mu_BC+ nui*growthVel_BC)/(M+1)

		if (meanC<0) then
			exit_stat = .true.
		endif

	end function Av_ceq


	function get_mass(z) result (mass)

		implicit none
		double precision, dimension(:), intent(in) :: z
		double precision	::	mass,r
		integer :: j

		mass = 0.0d0
		!x_BC = L+delta_BC

		do j=1,M-1
			r = (j-1) * deltax
			mass = mass + r *(hs-z(j)) * deltax
		enddo
		mass = 2.0*PI*(mass+(hs-z(M))*L*delta_BC)
	end function get_mass

	

	function get_flux() result(hydro_fluxL)
!+++++++ Problem: won't be very precise when expanding box. 
!		Would be more precise if these quantities calculated inside evolution function before box expansion


		implicit none
		!double precision new_2dterm
		integer j
		double precision hydro_fluxL,fluxLocal_BC,Q_BC,x_BC
		!double precision, dimension(:), intent(in) :: z

		x_BC = L+delta_BC
		!fluxLocal_BC = nui*(hs-h_BC)*(growthVel_BC - growth_vel(M))/delta_BC
		fluxLocal_BC = nui*(hs-h_BC)*(growth_vel(M) - growth_vel(M-1))*dx
		!Q_BC = (hs-h_BC)*(mu_BC - mu(M))/delta_BC
		 Q_BC = (hs-h_BC)*(mu(M) - mu(M-1))*dx !1st order and more stable of previous
		hydro_fluxL = -2.0d0*PI *x_BC*Q_BC -2.0d0*PI *x_BC*fluxLocal_BC + 2.0d0*PI*lateral_vel*x_BC*(hs-h_BC)

		!hydro_fluxL = -2.0d0*PI *L*Q(M) - PI*L**2 * crystal_vel -2.0d0*PI *L*F_local_gr 


	end function get_flux

	function get_average_height(z) result(mean)
		implicit none
		double precision :: mean
		double precision, dimension(:), intent(in) :: z
		integer :: j
	
		mean = 0.0d0
		do j =1, M
			mean = mean + z(j)
		enddo
		mean = mean/M
	end function get_average_height

	function get_cryst_force(z)
	implicit none

	double precision :: get_cryst_force
	double precision :: intforce,x
	double precision, dimension(:), intent(in) :: z

	integer :: i

	intforce = 0.0d0

	do i = 1, M-1
	x = (i-1)*deltax
	intforce = intforce +force(z(i)) * x 
	enddo

	intforce = intforce*deltax+ delta_BC * L*force(z(M))


	get_cryst_force =  2.0d0*PI*intforce


	end function get_cryst_force

	




end MODULE
