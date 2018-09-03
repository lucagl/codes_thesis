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
	
	function Av_ceq(exit_stat) result(meanc_eq)

		implicit none

		integer :: j
		!double precision, dimension (M) , intent(in) :: c_eq
		double precision :: meanc_eq,r1
		logical ,intent (inout) ::	exit_stat
		

		!c0 factorize and gamma in front of the curvature
		
		!c_eq(1) = c0*(1.0d0-gamma*2.0d0*(h(2) + h(2) - 2.0d0* h(1)) * sqrdx  + force(hs-h(j)))!using De L'Hopital relation + mirror boundary
		do j= 1,M-1
			!r1 = 1.0d0/((j-1) * deltax)
			!c_eq(j) = c0*(1.0d0- gamma*((h(j+1) + h(j-1) - 2.0d0* h(j))*sqrdx-r1*0.5d0*(h(j+1)-h(j-1))*dx) + force(hs-h(j))) 
			c_eq(j) = c0*(1.0d0 + mu(j))
		enddo
		
!		Tricky, if constant derivative h change and then external concentration vary
!		Therefore better to integrate mean conc without the bundary
		c_eq(M) = c_inf

		meanc_eq = 0.0d0
		do j=1, M-1
			meanc_eq = meanc_eq + c_eq(j)
		enddo

		meanc_eq = meanc_eq/(M-1)

		if (meanc_eq<0) then
				print*, 'error'
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
		double precision r1
		!double precision new_2dterm
		integer j
		double precision hydro_fluxL


		!hydro_fluxL = (hs-h(M))*(P(M)-(-(h(M)+h(M-2)-2*h(M-1))*sqrdx+force(hs-h(M-1))))*dx + crystal_vel * L
		! new_2dterm = 0.0d0
		 !do j = 2, M
		!	r= 1.0d0/((j-1)*deltax)
		!	new_2dterm = new_2dterm + Q(j)*r1 * deltax
		!enddo
		!new_2dterm = new_2dterm + onesixth*(8.0d0*Q(2)-Q(3))!adding contribution in 0
		!hydro_fluxL = Q(M) + new_2dterm + crystal_vel * L !a meno di 2pi
		
		hydro_fluxL = -2.0d0*PI *L*Q(M) - PI*L**2 * crystal_vel
		
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
