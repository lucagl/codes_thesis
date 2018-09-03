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
use integrate, only : sqrdx,dx,onesixth
implicit none

double precision new_2dterm

CONTAINS
	
	function Av_ceq(exit_stat) result(meanC)

		implicit none

		integer :: j
		!double precision, dimension (M) , intent(in) :: c_eq
		double precision :: meanC
		logical ,intent (inout) ::	exit_stat
		

	
		 

		meanC= 0.0d0
		do j=1, M
			c_eq(j) = c0*(1.0d0 + mu(j))!linear Gibbs Thomson relation
			meanC = meanC + c_eq(j)+ nui*growth_vel(j)
		enddo

		meanC = meanC/(M+1)

		if (meanC<0) then
		exit_stat = .true.
		endif

	end function Av_ceq


	function get_mass() result (mass)

		implicit none
		double precision	::	mass,r
		integer :: j

		mass = 0.0d0

		!do j =1, M
		!	mass = mass + h(j)*deltax!(((h(j) + h(j-1))/2.0d0))*deltax
		!enddo
		do j=1,M
			r = (j-1) * deltax
			mass = mass + r *(hs-h(j)) * deltax
		enddo
		mass = 2.0*PI*mass
	end function get_mass

	
	function get_flux() result(hydro_fluxL)

		implicit none
		!double precision new_2dterm
		integer j
		double precision hydro_fluxL,F_local_gr

		F_local_gr = 0.5d0*nui*(hs-h_BC)*(3.0d0*growth_vel(M) -4.0d0* growth_vel(M-1)+growth_vel(M-2))*dx
		hydro_fluxL = -2.0d0*PI *L*Q(M) - PI*L**2 * crystal_vel -2.0d0*PI *L*F_local_gr

	end function get_flux

	subroutine get_average_height(mean)
		implicit none
		double precision, intent(out) :: mean
		integer :: j
	
		mean = 0.0d0
		do j =1, M
			mean = mean + h(j)
		enddo
		mean = mean/M
	end subroutine get_average_height
	
	!function get_crystalisation() result(cryst_force)
	!	implicit none 
	!	double precision :: cryst_force
	!	cryst_force = intforce
		
	!	end function get_crystalisation
	




end MODULE
