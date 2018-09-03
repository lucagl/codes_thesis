MODULE read_input
use com_var
implicit none
logical, private :: check, slow_attachment!,check_integration
!logical, save :: read_last
integer, private ::  initial_generator
character (len=30) :: read_name,read_lenght,read_sup,read_load
integer, save :: read_potential
double precision, save :: hs_overwrite,d_read
double precision, dimension(:),allocatable, save :: h0
!integer, dimension (:), allocatable, save :: mask

CONTAINS

	subroutine read_input_file()
		inquire(file ='input.in',EXIST = check)
		if (check .eqv. .false.) then
		stop 'Input file not found.\n Remember it must be called input.in'
		else
		open (unit = 1,file = 'input.in', status = 'old', action = 'read')
		endif
		read (1,*) read_name
			call get_run_name(read_name)
		read(1,*) N
		read(1,*) deltat
		c0 =1.0d0
		gamma = 1.0d0
		read(1,*) read_potential
		if (read_potential==1) then
			read(1,*) d_read !Debyle lenght
		endif
		if (read_potential==3) then
			read(1,*) hs_overwrite
		endif
		if (read_potential==4) then
			read(1,*) d_read !Debyle lenght
			read(1,*) B !Potential strenght
			read(1,*) hs_overwrite
			read(1,*) hc
			!read(1,*) hc
		endif

		read(1,*) nu


		read(1,*) freq_step
		!read(1,*) freq_vel

		
		!read_last = iflast
		!if (iflast) then
		
		read(1,*) sup
		write(read_sup,'(F0.2)') sup
		call get_sup_name(read_sup)
		
!		read(1,*) iflast
!		print *, iflast
!		if(iflast .eqv. .false.) then
!			read(1,*) curv_BC !for cryst force read sandwith height and deduce mu_BC
!			!print *, 'curv_BC=', curv_BC
!		endif
		read(1,*) W2 !half width of the pore
		call read_old(run_name)
		close(1)


		write(read_lenght,'(F0.1)'), L
		
		call get_lenght_name(read_lenght)

	end subroutine read_input_file


	subroutine read_old(s)

		implicit none
		integer j
		character (len = *) :: s
		!integer :: ifPot,j
		double precision :: dummy,r1,max_curv
		integer :: sub_M,max_pos,sub_M_2
		double precision, dimension(3) :: line
		double precision, dimension(:), allocatable :: curv_0,der_0
		
		!if(iflast .eqv. .true.) then
			print*, '+++++++ Reading from EXPANDING BOX LC file +++++++++++'

			inquire(file = run_name//'_exp.lc',EXIST = check)
			if (check .eqv. .false.) stop 'Last configuration file does not exist'

			open (unit = 3,file=s//'_exp.lc', status = 'old', action = 'read')
			read(3,*) L, deltax,delta_BC
			print *,L,deltax,delta_BC
			M = int(L/deltax) +1
			print *, M
			allocate (h0(M))
			allocate (curv_0(M))
			allocate (der_0(M))

			read(3,*) hs
			read(3,*) h_BC
			read(3,*) (h0(j), j=1,M)

			der_BC = (h0(M) - h0(M-1))/deltax!first order
			curv_BC = 1.0d0/(sqrt(1.0d0+der_BC**2)*(W2-(hs-h_BC)))

			print*, 'Delta_BC=',  delta_BC + L
			print*, 'xbc = ', delta_BC + L
			print*, 'h_bc=', h_BC, 'zeta_BC=', (hs-h_BC)
			print* ,'der_BC = ', der_BC
			print* ,'curv_BC = ', curv_BC
			print*, 'W', 2*W2
			
			if(W2<(hs-h_BC)) stop 'Zeta_BC must be smaller than the half width. STOP'

			open(unit =9, file = 'initial_curvature.txt', action = 'write')	
			curv_0(1) = 2.0d0*1.0/3.0*(0.5d0*h0(3)+7.5d0*h0(1)-8.0d0*h0(2))/(deltax**2)
			write(9,*) 0,curv_0(1)
			do j = 2, M-1
			r1 = 1.0d0/((j-1) * deltax)
			curv_0(j) = -(h0(j+1)+h0(j-1)-2.0d0*h0(j))/(deltax**2) - r1*0.5d0*(h0(j+1)-h0(j-1))/(deltax)
			write(9,*) (j-1)*deltax,curv_0(j)
			enddo
			curv_0(M) =  -(h0(M-2)-2*h0(M-1)+h0(M))/(deltax**2) - 1.0d0/L*(h0(M)-h0(M-1))/(deltax)
			write(9,*) L,curv_0(M)
			
			open(unit =9, file = 'initial_der.txt', action = 'write')	
			der_0(1) = 0
			write(9,*) 0,0
			do j = 2, M-1
			der_0(j) = 0.5d0*(h0(j+1)-h0(j-1))/deltax
			write(9,*) (j-1)*deltax,der_0(j)
			enddo
			der_0(M) =  (h0(M) - h0(M-1))/deltax
			write(9,*) L,der_0(M)

			close(9)

			deallocate (curv_0)
			deallocate (der_0)

!		else
!			print*, '+++++++ Reading from FIXED BOX LC file +++++++++++'
!			inquire(file = run_name//'.lc',EXIST = check)
!			if (check .eqv. .false.) stop 'last configuration file does not exist'
!			open (unit = 3,file=s//'.lc', status = 'old', action = 'read')
!
!			read(3,*) L, deltax
!			print *,L,deltax
!			M = int(L/deltax) +1
!			allocate (h0(M))
!			allocate (curv_0(M))
!			!allocate (mask(M))
!			read(3,*) hs
!			read(3,*) dummy !to account for old way of printing
!			read(3,*) (h0(j), j=1,M)
!			!print *, h0(1)
!			close(3)
!			
!			
!
!			der_BC = (h0(M) - h0(M-1))/deltax
!
!			!curv_BC = 1.0d0/(sqrt(1+der_BC**2)*(W/2-h_BC))
!			print *, 'Der_BC =', der_BC
!			!print *, 'curv_BC=', curv_BC
!			print *,'W=', W
!			!++++++++ INITIALIZATION OF X_bc and H_bc	
!
!			print*, 'INITIALIZING x_BC and h_BC'
!
!
!			! Compute curvature everywhere of initial profile then compare to chosen curvature
!			!
!			open(unit =9, file = 'initial_curvature', action = 'write')	
!			curv_0(1) = 2.0d0*1.0/3.0*(0.5d0*h0(3)+7.5d0*h0(1)-8.0d0*h0(2))/(deltax**2)
!			write(9,*) 0,curv_0(1)
!			do j = 2, M-1
!				r1 = 1.0d0/((j-1) * deltax)
!				curv_0(j) = -(h0(j+1)+h0(j-1)-2.0d0*h0(j))/(deltax**2) - r1*0.5d0*(h0(j+1)-h0(j-1))/(deltax)
!				write(9,*) (j-1)*deltax,curv_0(j)
!			enddo
!			curv_0(M) =  -(h0(M-2)-2*h0(M-1)+h0(M))/(deltax**2) - 1.0d0/L*(h0(M)-h0(M-1))/(deltax)
!			write(9,*) L,curv_0(M)
!
!			close(9)
!			
!			open(unit =9, file = 'initial_der', action = 'write')	
!			der_0(1) = 0
!			write(9,*) 0,0
!			do j = 2, M-1
!			der_0(j) = 0.5d0*(h0(j+1)-h0(j-1))*dx
!			write(9,*) (j-1)*deltax,der_0(j)
!			enddo
!			der_0(M) =  (h0(M) - h0(M-1))/deltax
!			write(9,*) L,curv_0(M)
!
!			close(9)
!
!
!
!			max_curv = maxval(curv_0)
!			
!			if(max_curv<curv_BC) stop 'Cannot realise the boundary curvature'
!
!
!			print *, M-minloc(abs(curv_0(M:1:-1)-curv_BC),1), (M-minloc(abs(curv_0(M:1:-1)-curv_BC),1))*deltax
!			sub_M = (M-minloc(abs(curv_0(M:1:-1)-curv_BC),1))+1
!
!			max_pos = maxloc(curv_0,1)
!
!			print *, M-minloc(abs(curv_0(M:max_pos:-1)-curv_BC),1), (M-minloc(abs(curv_0(M:max_pos:-1)-curv_BC),1))*deltax
!			sub_M_2 = M-minloc(abs(curv_0(M:max_pos:-1)-curv_BC),1)
!
!			if(sub_M_2>sub_M) sub_M = sub_M_2
!
!			!y = h0(sub_M-2:sub_M+2)
!			!print *,curv_0(sub_M-2:sub_M+2)
!			line = linreg(curv_0(sub_M-1:sub_M+1),5)
!
!			delta_BC=(curv_BC-line(2))/line(1)
!
!			print *, delta_BC
!			print *, 'x_BC = ',((sub_M-1)-1)*deltax + delta_BC
!			
!			L = ((sub_M-1)-1)*deltax+int(delta_BC/deltax) * deltax
!
!			M = sub_M-1+int(delta_BC/deltax)
!
!			delta_BC = DMOD(delta_BC,deltax)
!
!			
!
!			!print * ,'M=', M,'L=',L 
!			!sub_M = minloc(abs(curv_0-mu_BC),1)
!			
!			h_BC = h0(M) + delta_BC * 0.5d0*(h0(M-2)-4.0d0*h0(M-1)+3.0d0*h0(M))/deltax
!			!h_BC = h0(M) + delta_BC * (h0(M)-h0(M-1))/deltax
!
!			print *, 'Check: '
!			print *, 'Delta_BC =', delta_BC, 'M =', M , 'L =', L, '=',(M-1)*deltax
!
!			!print *, 'der_BC = ', der_BC
!
!
!			!++++++++++++++++
!
!
!		endif


		
	end subroutine read_old


subroutine get_run_name(input_string)

implicit none 

character(len=*), intent(in) :: input_string
!character(len=:), allocatable, intent(out) :: run_name

allocate(character(len=LEN(trim(input_string))):: run_name)
run_name= trim(input_string)

end subroutine get_run_name



subroutine get_lenght_name(input_string)

implicit none 

character(len=*), intent(in) :: input_string
!character(len=:), allocatable, intent(out) :: run_name

allocate(character(len=LEN(trim(input_string))+1):: lenght)

lenght= trim('_initial_L'//trim(input_string))

end subroutine get_lenght_name


subroutine get_sup_name(input_string)

implicit none 

character(len=*), intent(in) :: input_string
!character(len=:), allocatable, intent(out) :: run_name

allocate(character(len=LEN(trim(input_string))+5):: supersaturation)

supersaturation= trim('cinf_'//trim(input_string))

end subroutine get_sup_name


subroutine get_load_name(input_string)

implicit none 

character(len=*), intent(in) :: input_string
!character(len=:), allocatable, intent(out) :: run_name

allocate(character(len=LEN(trim(input_string))+5):: ext_force)

ext_force= trim('load_'//trim(input_string))

end subroutine get_load_name



function linreg (y,N) result(l)

implicit none 

double precision, dimension (:) :: y
double precision, dimension (3) :: l !result
integer :: N,j
double precision :: x

double precision         ::  b                                                        ! y-intercept of least-squares best fit line
double precision           ::  m                                                        ! slope of least-squares best fit line
double precision           ::  s = 0.0d0                                                ! number of data points
double precision          ::  r                                                        ! squared correlation coefficient
character (len=80)  ::  str                                                      ! input string
double precision           ::  sumx  = 0.0d0                                            ! sum of x
double precision         ::  sumx2 = 0.0d0                                            ! sum of x**2
double precision        ::  sumxy = 0.0d0                                            ! sum of x * y
double precision          ::  sumy  = 0.0d0                                            ! sum of y
double precision          ::  sumy2 = 0.0d0                                            ! sum of y**2
! input y data

!write (unit=*, fmt="(a)") " LINREG - Perform linear regression"                  ! print introductory message
!write (unit=*, fmt="(a/)") "   (Enter END to stop data entry and compute"//  &
!" linear regression.)"

do j = 1, N                                                                     ! loop for all data points
x = (j-1)*deltax
s = s + 1.0d0                                                                 ! increment number of data points by 1
sumx  = sumx + x                                                            ! compute sum of x
sumx2 = sumx2 + x*x                                                           ! compute sum of x**2
sumxy = sumxy + x * y(j)                                                            ! compute sum of x * y
sumy  = sumy + y(j)                                                              ! compute sum of y
sumy2 = sumy2 + y(j) * y(j)                                                         ! compute sum of y**2
end do

m = (s * sumxy  -  sumx * sumy) / (s * sumx2 - sumx**2)                          ! compute slope
b = (sumy * sumx2  -  sumx * sumxy) / (s * sumx2  -  sumx**2)                    ! compute y-intercept
r = (sumxy - sumx * sumy / s) /                                     &            ! compute correlation coefficient
sqrt((sumx2 - sumx**2/s) * (sumy2 - sumy**2/s))

write (unit=*, fmt="(/a,es15.6)") " Slope        m = ", m                        ! print results
write (unit=*, fmt="(a, es15.6)") " y-intercept  b = ", b
write (unit=*, fmt="(a, es15.6)") " Correlation  r = ", r

l(1) = m
l(2) = b
l(3) = r


end function linreg




!function linreg (x,y,N) result(line)
!
!implicit none 
!
!double precision, dimension (:) :: x,y
!double precision, dimension (3) :: line !result
!integer :: N,j
!
!double precision         ::  b                                                        ! y-intercept of least-squares best fit line
!double precision           ::  m                                                        ! slope of least-squares best fit line
!double precision           ::  s = 0.0d0                                                ! number of data points
!double precision          ::  r                                                        ! squared correlation coefficient
!character (len=80)  ::  str                                                      ! input string
!double precision           ::  sumx  = 0.0d0                                            ! sum of x
!double precision         ::  sumx2 = 0.0d0                                            ! sum of x**2
!double precision        ::  sumxy = 0.0d0                                            ! sum of x * y
!double precision          ::  sumy  = 0.0d0                                            ! sum of y
!double precision          ::  sumy2 = 0.0d0                                            ! sum of y**2
!                                                     ! input y data
!
!!write (unit=*, fmt="(a)") " LINREG - Perform linear regression"                  ! print introductory message
!!write (unit=*, fmt="(a/)") "   (Enter END to stop data entry and compute"//  &
!!" linear regression.)"
!
!do j = 1, N                                                                     ! loop for all data points
!
!s = s + 1.0d0                                                                 ! increment number of data points by 1
!sumx  = sumx + x(j)                                                              ! compute sum of x
!sumx2 = sumx2 + x(j)*x(j)                                                            ! compute sum of x**2
!sumxy = sumxy + x(j) * y(j)                                                            ! compute sum of x * y
!sumy  = sumy + y(j)                                                              ! compute sum of y
!sumy2 = sumy2 + y(j) * y(j)                                                         ! compute sum of y**2
!end do
!
!m = (s * sumxy  -  sumx * sumy) / (s * sumx2 - sumx**2)                          ! compute slope
!b = (sumy * sumx2  -  sumx * sumxy) / (s * sumx2  -  sumx**2)                    ! compute y-intercept
!r = (sumxy - sumx * sumy / s) /                                     &            ! compute correlation coefficient
!sqrt((sumx2 - sumx**2/s) * (sumy2 - sumy**2/s))
!
!write (unit=*, fmt="(/a,es15.6)") " Slope        m = ", m                        ! print results
!write (unit=*, fmt="(a, es15.6)") " y-intercept  b = ", b
!write (unit=*, fmt="(a, es15.6)") " Correlation  r = ", r
!
!line(1) = m
!line(2) = b
!line(3) = r
!
!
!end function linreg
!






end MODULE
