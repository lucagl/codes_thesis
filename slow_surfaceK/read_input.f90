MODULE read_input
use com_var
implicit none
logical, private :: iflast, check, slow_attachment!,check_integration
logical, save :: read_last
integer, private ::  initial_generator
character (len=30) :: read_name,read_lenght,read_sup,read_load
integer, save :: read_boundary,read_potential
double precision, save :: hs_overwrite,d_read
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
		!read(1,*) c0
		!read(1,*) gamma
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

		read(1,*) iflast
		read_last = iflast
		if (iflast) then
		inquire(file = run_name//'.lc',EXIST = check)
			if (check .eqv. .false.) stop 'last configuration file does not exist'
			call read_old(run_name)
			read(1,*) sup
			write(read_sup,'(F0.2)') sup
			call get_sup_name(read_sup)
		go to 10
		else
			if (.not.(read_potential==2)) stop 'For generating new configuration needs attractive potential' !this could be avoided, is more an historical problem..
		endif
		read(1,*) L
		read(1,*) deltax
		M = int(L/deltax) +1
		print *, 'seeking equilibrium conf'
		print *, 'curv_BC = 0, der_BC = equilibrium inf crystal value'
		read(1,*) h_BC
		read(1,*) hs
		10 continue
		read(1,*) load
			write(read_load,'(I0)') int(load)
			call get_load_name(read_load)
		read(1,*) viscosity
		close(1)
		write(read_lenght,'(F0.1)'), L
		
		call get_lenght_name(read_lenght)

	end subroutine read_input_file


	subroutine read_old(s)

		implicit none
		integer j
		character (len = *) :: s
		!integer :: ifPot,j

		print*, '+++++++ Reading from lc file +++++++++++'

		open (unit = 3,file=s//'.lc', status = 'old', action = 'read')

		read(3,*) L, deltax
		print *,L,deltax
		M = int(L/deltax) +1
		allocate (h0(M))
		allocate(h(M))
		allocate(c_eq(M))
		read(3,*) hs
		print *,hs
		read(3,*) der_BC
		read(3,*) (h0(j), j=1,M)
		!print *, h0(1)
		h_BC = h0(M)
		h(:) = h0(:) !assign first value to h0
		
		close(3)

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

lenght= trim('L'//trim(input_string))

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

end MODULE
