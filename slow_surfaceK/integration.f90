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
double precision, dimension(:),allocatable, save :: temp_force,d
double precision, dimension(:,:),allocatable, save :: A,Ainv
!double precision, dimension(M), save :: temp_force,growth_vel !dynamical local array (automatically erased when exiting routine)
public init_int,sqrdx,dx



!private P ! Q might be used to print flux


CONTAINS

    function step(z)
        implicit none

        double precision, dimension(M), intent(in) :: z
        integer ::      j
        double precision, dimension (M) ::  step !step is a M dimensional array (only spatial indexes)
		double precision :: h_der
	

		crystal_vel = hydro(z,pressure) !returns also the pressure 
		growth_vel = get_local_growth(z)


		




!HERE IS THE TIME STEP #############

		!step(1) = z(1) + deltat * (2.0d0*derQ1 + crystal_vel+2.0d0*derF1)!using De L'Hopital relation


!----------

        do j = 1, M-1
			!r1 = 1.0d0/((j-1) * deltax)
			!h_der_1 = 0.5d0*(Q(j+1)-Q(j-1))* dx + r1*Q(j)  + crystal_vel + 0.5d0*(F(j+1)-F(j-1))* dx + r1*F(j) 
			h_der= growth_vel(j) + crystal_vel

			!profile evolution
            step(j)= z(j) + deltat * h_der
        enddo
		
		step(M) = h_BC !changes here when growing box..


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

    subroutine init_int()

		implicit none
		integer :: j,i
		double precision :: r1,zita,der_z,derQ1
	!	double precision, dimension(M,M) :: A, Ainv !local dynamical array (used just once and deleted)

        allocate(Q(M))!diffusion flux of mass
		!allocate(F(M))!attachment flux of mass
        allocate(mu(M))!chemical potential
		allocate(temp_force(M))
		allocate(pressure(M))
		allocate(growth_vel(M))!\partial_t h -u_cz
		allocate(A(M,M),Ainv(M,M))
		allocate(d(M))
		


		sqrdx = 1.0d0/(deltax*deltax)
		dx = 1.0d0/deltax
		onethrd = 1.0d0/(3.0d0)
		onesixth = 1.0d0/(6.0d0)
		nui = 1.0d0/nu

		crystal_vel = hydro(h0,pressure)

		mu(1) = 2.0d0*onethrd*(0.5d0*h0(3)+7.5d0*h0(1)-8.0d0*h0(2))*sqrdx + force(h0(1))!using De L'Hopital relation
		do j = 2, M-1
			r1 = 1.0d0/((j-1) * deltax)
			mu(j) = -(h0(j+1)+h0(j-1)-2.0d0*h0(j)) * sqrdx - r1*0.5d0*(h0(j+1)-h0(j-1))*dx + force(h0(j))
		enddo
		!mu(M) = (c_inf-1.0d0) !boundary supersaturation
		mu(M) = sup+nui*crystal_vel
		Q(1) = 0.0d0 !no source or sink in the middle
		do j=2,M-1
			Q(j) = 0.5d0*(hs - h0(j)) * (mu(j+1) - mu(j-1)) * dx !Flux of diffusion initial conf
		enddo
		!Q(M) = (hs-h_BC)*(mu(M) - mu(M-1))*dx
		Q(M) = 0.5d0*(hs-h_BC)*(3.0d0*mu(M) -4.0d0* mu(M-1)+mu(M-2))*dx
!INITIALIZATION OF local Crystallization vel #############
	!option 1: solve linear system
	
!build matrix #############
	A =0.0d0
		do j=1,M 
			do i = 2,M-1
				r1 = 1.0d0/((i-1) * deltax)
				zita = (hs-h0(i)) 
				der_z = -0.5d0*dx*(h0(i+1)-h0(i-1))!- because is zita
				d(i) = 0.5d0*(Q(i+1)-Q(i-1))* dx + r1*Q(i) 
				if(j==i) then
					A(i,j) =  1.0d0+ nui*2.0d0*zita *sqrdx
				elseif(j==(i+1)) then 
					A(i,j) = - nui*(zita*(0.5d0*dx*r1+sqrdx)+0.5d0*dx*der_z)
				elseif(j==(i-1)) then 
					A(i,j) = nui*(zita*(0.5d0*dx*r1-sqrdx)+0.5d0*dx*der_z)
				endif
			enddo
		enddo
	!boundary conditions
	derQ1 = onesixth*(8.0d0*Q(2)-Q(3))*dx
	d(1) = 2.0d0*derQ1 
	d(M)= -crystal_vel

	A(1,1) = 1.0d0+nui*5.0d0*(hs-h0(1))*sqrdx
	A(1,2) = -nui*16.0d0/3.0d0*(hs-h0(1))*sqrdx
	A(1,3) = +nui*(hs-h0(1))/3.0d0*sqrdx
	!A(1,4:)=0.0d0
	A(M,M) = 1
!boundary conditions missing for A(M,M-1) ? not sure



	!invert matrix
	Ainv = invert(A)
	growth_vel = 0
	!compute initial value of local growth_vel
	do j=1,M
		do i = 1,M
			growth_vel(i)=growth_vel(i)+ Ainv(i,j)*d(j)
		enddo
	enddo

	!growth_vel = 0.0d0 just to try..


!+++++++ CHECK +++++++++


!open(unit = 8, file = 'test_d.txt',action = 'write')
!
!d_test = 0
!
!do j=1,M
!	do i = 1,M
!		d_test(i) = d_test(i) + A(i,j) * growth_vel(j)
!	enddo
!enddo
!
!do j=1,M
!	write(8,*) (j-1)*deltax,d(j),d_test(j)
!enddo
!close(8)
!+++++++++++++++++++

end subroutine init_int

function get_local_growth(z) result(local_vel)

	implicit none
	double precision, dimension(:), intent(in) :: z
	double precision, dimension (size(z)) ::  local_vel
	double precision r1,zita,der_z,derQ1
	integer :: i,j

	mu(1) = 2.0d0*onethrd*(0.5d0*z(3)+7.5d0*z(1)-8.0d0*z(2))*sqrdx + temp_force(1)!using De L'Hopital relation

	do j = 2, M-1
		r1 = 1.0d0/((j-1) * deltax)
		mu(j) = (-(z(j+1)+z(j-1)-2.0d0*z(j)) * sqrdx - r1*0.5d0*(z(j+1)-z(j-1))*dx) + temp_force(j)
	enddo
	mu(M) = sup+nui*crystal_vel
	
	Q(1) = 0.0d0
	do j=2,M-1
		Q(j) = 0.5d0*(hs - z(j)) * (mu(j+1) - mu(j-1)) * dx 
	enddo
	!Q(M) = (hs-h_BC)*(mu(M) - mu(M-1))*dx
	Q(M) = 0.5d0*(hs-h_BC)*(3.0d0*mu(M) -4.0d0 *mu(M-1)+mu(M-2))*dx
	derQ1 = onesixth*(8.0d0*Q(2)-Q(3))*dx
	!build matrix #############
	A =0.0d0
	do j=1,M 
		do i = 2,M-1
			r1 = 1.0d0/((i-1) * deltax)
			zita = (hs-z(i)) 
			der_z = -0.5d0*dx*(z(i+1)-z(i-1))!- because is (hs-h)
			d(i) = 0.5d0*(Q(i+1)-Q(i-1))* dx + r1*Q(i) 
			if(j==i) then
				A(i,j) =  1.0d0+ nui*2.0d0*zita *sqrdx
			elseif(j==(i+1)) then 
				A(i,j) = - nui*(zita*(0.5d0*dx*r1+sqrdx)+0.5d0*dx*der_z)
			elseif(j==(i-1)) then 
				A(i,j) = nui*(zita*(0.5d0*dx*r1-sqrdx)+0.5d0*dx*der_z)
			endif
		enddo
	enddo
	!boundary conditions

	d(1) = 2.0d0*derQ1 
	!crystal_vel = hydro(z,pressure)
	d(M)= -crystal_vel 
	A(1,1) = 1.0d0+nui*5.0d0*(hs-z(1))*sqrdx
	A(1,2) = -nui*16.0d0/3.0d0*(hs-z(1))*sqrdx
	A(1,3) = +nui*(hs-z(1))/3.0d0*sqrdx
	!A(1,4:)=0.0d0
	A(M,M) = 1
	!boundary conditions missing for A(M,M-1) ? not sure



	!invert matrix
	Ainv = invert(A)
	local_vel = 0
	do j=1,M
		do i = 1,M
			local_vel(i) =local_vel(i)+ Ainv(i,j)*d(j)
		enddo
	enddo



end function get_local_growth



function invert(A) result(Ainv)



double precision, dimension(:,:), intent(in) :: A
double precision, dimension(size(A,1),size(A,2)) :: Ainv

double precision, dimension(size(A,1)) :: work  ! work array for LAPACK
integer, dimension(size(A,1)) :: ipiv   ! pivot indices
integer :: n, info

! External procedures defined in LAPACK
external DGETRF
external DGETRI

! Store A in Ainv to prevent it from being overwritten by LAPACK
Ainv = A
n = size(A,1)

! DGETRF computes an LU factorization of a general M-by-N matrix A
! using partial pivoting with row interchanges.
call DGETRF(n, n, Ainv, n, ipiv, info)

if (info /= 0) then
stop 'Matrix is numerically singular!'
end if

! DGETRI computes the inverse of a matrix using the LU factorization
! computed by DGETRF.
call DGETRI(n, Ainv, n, ipiv, work, n, info)

if (info /= 0) then
stop 'Matrix inversion failed!'
end if
end function invert



end MODULE






