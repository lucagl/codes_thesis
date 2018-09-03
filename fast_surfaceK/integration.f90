!-------------------------------------------------------------------------------------------
!
!
!				CONTAINS THE SPECIFIC FUNCTIONS RELATED TO THE INTEGRATION
!
!
!-------------------------------------------------------------------------------------------


MODULE integrate

use com_var

implicit none


double precision, save :: sqrdx, dx,onethrd,onesixth
double precision, dimension(:),allocatable, private, save :: temp_force

public init_int,sqrdx,dx



!private P ! Q might be used to print flux


CONTAINS

    function step(z)
        implicit none

        double precision, dimension(:), intent(in) :: z
        integer ::      j
        double precision, dimension (size(z)) ::  step !step is a M dimensional array (only spatial indexes)
		double precision :: r1,derQ1

		crystal_vel = hydro(z,pressure) !returns also the pressure 

		mu(1) = 2.0d0*onethrd*gamma*(0.5d0*z(3)+7.5d0*z(1)-8.0d0*z(2))*sqrdx + temp_force(1)!using De L'Hopital relation
        do j = 2, M-1
			r1 = 1.0d0/((j-1) * deltax)
            mu(j) = gamma*(-(z(j+1)+z(j-1)-2.0d0*z(j)) * sqrdx - r1*0.5d0*(z(j+1)-z(j-1))*dx) + temp_force(j)
        enddo

       
		!Q(2) = (hs-z(2))*0.5d0*(4.0d0*mu(3)-mu(4)-3.0d0*mu(2))*dx!to avoid using mu(1)
        do j=2,M-1
            Q(j) = 0.5d0*(hs - z(j)) * (mu(j+1) - mu(j-1)) * dx !Flux of mass
        enddo
		
		Q(M) = (hs-h_BC)*(mu(M) - mu(M-1))*dx
		!Q(M) = 0.5d0*(hs-h_BC)*(3.0d0*mu(M) -4.0d0* mu(M-1)+mu(M-2))*dx
		derQ1 = onesixth*(8.0d0*Q(2)-Q(3))*dx
	
		step(1) = z(1) + deltat * (2.0d0*derQ1 + crystal_vel)!using De L'Hopital relation

		!step(1) = onethrd * (4.0d0*z(2)-z(3))
	
		!HERE IS THE TIME STEP
        do j = 2, M-1
			r1 = 1.0d0/((j-1) * deltax)
            step(j)= z(j) + deltat * (0.5d0*(Q(j+1)-Q(j-1))* dx + r1*Q(j)  + crystal_vel)
        enddo

		step(M) = h_BC
		

    end function step

	function hydro(z,pres)
		implicit none
		double precision, dimension(M), intent(in) :: z
		double precision, dimension(M), intent(inout) :: pres
		double precision :: hydro
		!double precision :: friction
		double precision :: s,zita,x
		double precision, dimension(M) :: a,c !dummy variable
		integer :: i


		intforce = 0.0d0
		friction = 0.0d0
		s = 0.0d0


		
		
		a(1) = 0
		do i = 2, M
		x = (i-2)*deltax
		zita = (hs - z(i-1))
		s = s + (6.0d0 *viscosity*(x/(zita*zita*zita)))
		a(i) = s*deltax
		enddo

		c(:)=a(M)-a(:) !is the pressure minus a constant and divided by the crystal velocity
		!because int_r^x_bc = int_0^x_BC-int_0^r
		!and because int_0^x dr f(r) = \sum_{r=0}^{x-1} f(r) * dr. This is why I compute s in i-1

		do i = 1, M-1
		temp_force(i) = force(z(i))
		x = (i-1)*deltax
		intforce = intforce +temp_force(i) * x * deltax
		friction = friction + deltax*c(i)*x
		enddo

		hydro = (-load -2.0d0*PI*intforce)/(2.0d0*PI*friction)
		pres = c(:) * hydro
		
	end function hydro
!-------------------------------

    subroutine init_int

		implicit none

		integer :: j
		double precision :: r1
        allocate(Q(M))!flux of mass
        allocate(mu(M))!chemical potential
		allocate(temp_force(M))
		allocate(pressure(M))
		

		sqrdx = 1.0d0/(deltax*deltax)
		dx = 1.0d0/deltax
		onethrd = 1.0d0/(3.0d0)
		onesixth = 1.0d0/(6.0d0)
		mu(1) = 2.0d0*onethrd*gamma*(0.5d0*h0(3)+7.5d0*h0(1)-8.0d0*h0(2))*sqrdx + force(h0(1))!using De L'Hopital relation
		!print *, mu(1)
		do j = 2, M-1
			r1 = 1.0d0/((j-1) * deltax)
			mu(j) = -(h0(j+1)+h0(j-1)-2.0d0*h0(j)) * sqrdx - r1*0.5d0*(h0(j+1)-h0(j-1))*dx + force(h0(j))
		enddo
		mu(M) = (c_inf-1.0d0)
       
		Q(1) = 0.0d0 !no source or sink in the middle
		do j=2,M-1
			Q(j) = 0.5d0*(hs - h0(j)) * (mu(j+1) - mu(j-1)) * dx !Flux of mass
		enddo
		!Q(M) = 0.5d0*(hs-h_BC)*(3.0d0*mu(M) -4.0d0* mu(M-1)+mu(M-2))*dx
		Q(M) = (hs-h_BC)*(mu(M)-mu(M-1))*dx
		
		!crystal_vel = hydro(h0,pressure)

    end subroutine init_int



end MODULE






