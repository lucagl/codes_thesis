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

double precision, save :: sqrdx,onethrd,onesixth,dx
!double precision, save :: der_BC



CONTAINS

    function step(z)
        implicit none

        double precision, dimension(M), intent(in) :: z
        integer ::      j
        double precision, dimension (M) ::  step 
		double precision :: h_der
		double precision :: curv_L


		growth_vel = get_local_growth(z)

        do j = 1, M
            step(j)= z(j) + deltat * growth_vel(j)
        enddo


		!step(M) = (h_BC-0.5d0*delta_BC*dx*(step(M-2)-4.0d0*step(M-1)))*deltax/(1.5d0*delta_BC+deltax)

		!der_BC = (step(M) - step(M-1))*dx



    end function step


!-------------------------------

    subroutine init_int(z0)

		implicit none
		integer :: j,i
		double precision, dimension (M) :: z0
		double precision :: r1,zita,der_z,derQ1
		double precision, dimension (M) :: d
		double precision, dimension (M,M) :: A,Ainv
        allocate(Q(M_max))!diffusion flux of mass
        allocate(mu(M_max))!chemical potential
		allocate(growth_vel(M_max))!\partial_t h -u_cz
	


		dx = 1.0d0/deltax
		sqrdx = 1.0d0/(deltax*deltax)
		onethrd = 1.0d0/(3.0d0)
		onesixth = 1.0d0/(6.0d0)
		nui = 1.0d0/nu
		
	

		mu_BC = curv_BC+force(h_BC)


		der_BC = (z0(M) - z0(M-1) - &
		&(0.5*deltax*deltax + deltax*delta_BC)*curv_BC)/(deltax+0.5*deltax*deltax/L + deltax*Delta_BC/L)

		lateral_vel = -1.0d0/der_BC * nu*(1.0d0+sup-exp(mu_BC))

		


		mu(1) = 2.0d0*onethrd*gamma*(0.5d0*z0(3)+7.5d0*z0(1)-8.0d0*z0(2))*sqrdx + force(z0(1))!using De L'Hopital relation
		do j = 2, M-1
			r1 = 1.0d0/((j-1) * deltax)
			mu(j) = -(z0(j+1)+z0(j-1)-2.0d0*z0(j)) * sqrdx - r1*0.5d0*(z0(j+1)-z0(j-1))*dx + force(z0(j))
		enddo
		!mu(M) = (mu_BC + mu(M-1)*delta_BC*dx)*deltax/(deltax+delta_BC)
		mu(M) =(mu_BC-0.5d0*delta_BC*dx*(mu(M-2)-4.0d0*mu(M-1)))*deltax/(1.5d0*delta_BC+deltax)
		Q(1) = 0.0d0 !no source or sink in the middle
		do j=2,M-1
			Q(j) = 0.5d0*(hs - z0(j)) * exp(mu(j))*(mu(j+1) - mu(j-1)) * dx !Flux of diffusion initial conf
		enddo
		Q(M) = 0.5d0*(hs-z0(M))*exp(mu(M))*(3.0d0*mu(M) -4.0d0* mu(M-1)+mu(M-2))*dx
		
!INITIALIZATION OF local Crystallization vel #############
	
!build matrix #############
	A =0.0d0
	do j=1,M 
		do i = 2,M-1
			r1 = 1.0d0/((i-1) * deltax)
			zita = (hs-z0(i)) 
			der_z = -0.5d0*dx*(z0(i+1)-z0(i-1))!- because is zita
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
	
	growthVel_BC = nu*(1.0d0+sup-exp(mu_BC))!= -crystal_vel - lateral_vel*der_BC

	d(1) = 2.0d0*derQ1 
	d(M) = growthVel_BC


	A(1,1) = 1.0d0+nui*5.0d0*(hs-z0(1))*sqrdx
	A(1,2) = -nui*16.0d0/3.0d0*(hs-z0(1))*sqrdx
	A(1,3) = +nui*(hs-z0(1))/3.0d0*sqrdx
	!A(M,M-1) = -delta_BC*dx
	!A(M,M) = 1.0d0+delta_BC*dx
	A(M,M-2) = 0.5d0*Delta_BC*dx
	A(M,M-1) = -2.0d0*Delta_BC*dx
	A(M,M) = 1.0d0+1.5d0*Delta_BC*dx
	!invert matrix
	Ainv = invert(A)
	growth_vel = 0
	!compute initial value of local growth_vel
	do j=1,M
		do i = 1,M
			growth_vel(i) =growth_vel(i)+ Ainv(i,j)*d(j)
		enddo
	enddo

	!growth_vel = 0.0d0 just to try..



end subroutine init_int


function get_local_growth(z) result(local_vel)

	implicit none
	double precision, dimension(:), intent(in) :: z
	double precision, dimension (size(z)) ::  local_vel
	double precision r1,zita,der_z,derQ1
	double precision, dimension (M) :: d
	double precision, dimension (M,M) :: A,Ainv !need local array to avoid problems with matrix inversion
	integer :: i,j



	mu(1) = 2.0d0*onethrd*gamma*(0.5d0*z(3)+7.5d0*z(1)-8.0d0*z(2))*sqrdx + force(z(1))!using De L'Hopital relation



	do j = 2, M-1
		r1 = 1.0d0/((j-1) * deltax)
		mu(j) = gamma*(-(z(j+1)+z(j-1)-2.0d0*z(j)) * sqrdx - r1*0.5d0*(z(j+1)-z(j-1))*dx) + force(z(j))
	enddo
	!mu(M) = (mu_BC + mu(M-1)*delta_BC*dx)*deltax/(deltax+delta_BC)
	mu(M) = (mu_BC-0.5d0*delta_BC*dx*(mu(M-2)-4.0d0*mu(M-1)))*deltax/(1.5d0*delta_BC+deltax)

	do j=2,M-1
		Q(j) = 0.5d0*(hs - z(j)) * exp(mu(j))*(mu(j+1) - mu(j-1)) * dx 
	enddo
	Q(M) = 0.5d0*(hs-z(M))*exp(mu(M))*(3.0d0*mu(M) -4.0d0* mu(M-1)+mu(M-2))*dx

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

	growthVel_BC = nu*(1.0d0+sup-exp(mu_BC))

	d(1) = 2.0d0*derQ1 
	d(M) = growthVel_BC

	A(1,1) = 1.0d0+nui*5.0d0*(hs-z(1))*sqrdx
	A(1,2) = -nui*16.0d0/3.0d0*(hs-z(1))*sqrdx
	A(1,3) = +nui*(hs-z(1))/3.0d0*sqrdx
!A(M,M-1) = -delta_BC*dx
!A(M,M) = 1.0d0+delta_BC*dx
	A(M,M-2) = 0.5d0*Delta_BC*dx
	A(M,M-1) = -2.0d0*Delta_BC*dx
	A(M,M) = 1.0d0+1.5d0*Delta_BC*dx
	



	!invert matrix
	Ainv = invert(A)
	local_vel = 0

	!compute initial value of local growth_vel
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






