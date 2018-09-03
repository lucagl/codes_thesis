MODULE  write_output!contains also time information

use com_var

!use ifport

implicit none
double precision, private :: curv_L,dri2,L_contact
CHARACTER*160, private ::  fileplace
!double precision, parameter, private :: target_value = 0.000001

CONTAINS

	subroutine write_main_out(s,frame,z)
		implicit none 
		integer, intent(in) :: s
		integer,optional :: frame	
		double precision :: h_max
		double precision,optional, dimension(:), intent(in) :: z
		!if(present(frame)) fr_number = frame
		if (s == 0) then
			open(unit=2,file='OUTPUT'//lenght//supersaturation//'.out',action ='write')
			write (2,*) ''
			write (2,*) ''
			write ( 2, 1000 )  yesterday(1), yesterday(2), yesterday(3), before
			1000 format ( 'Date ', i2.2, '/', i2.2, '/', i4.4, '; time ',i2.2, ':', i2.2, ':', i2.2 )
			write (2,*) ''
			write (2,*) ''
			write (2,'(A,A)') trim(run_name),  '       OUTPUT      '
			write (2,'(A,3L)') 'Reading from last configuration'
			write (2,*) 'Maximum number of time steps  ', N
			write (2,*) 'Delta t  ', deltat
			write (2,*) 'Frame output frequency  ', freq_step
			write (2,*) 'Initial Lenght of the box   ', L
			write (2,*) 'Delta x  ', deltax
			write (2,*) 'Kinetic constant ', nu
			write (2,*) 'Height of the substrate hs = ', hs

			
			write (2,*) 'Initial boundary derivative', ((h_BC-z(M))/delta_BC)
			
			curv_L=-(z(M-2) - 2.0d0*z(M-1) + z(M))/(deltax**2)- 1.0d0/L*(z(M)-z(M-1))/deltax !2d curvature
			write(2,*) 'Initial curvature in L(last bin) =',  curv_L
			write (2,*) 'Initial boundary height and position: x_bc',L+delta_BC,'h_BC= ', h_BC
			write (2,*) 'Boundary supersaturation', sup, 'force at the boundary = ', force(h_BC)
			write (2,*) 'Boundary curvature', curv_BC, 'boundary chemical potential= ', mu_BC

			write (*,*) 'Initial boundary position and height : x_bc',L+delta_BC,'h_BC= ', h_BC
			write (*,*) 'Force at the boundary', force(h_BC)
			write (*,*) 'Boundary curvature', curv_BC, 'boundary chemical potential= ', mu_BC
			
			write(2,*) ''
			write(2,*) ''
			write(2,*) 'frame          mean_conc               mass	   	      flux	       M(0) + cumul_flux(L)	size		force'
		
		elseif(s==1) then
			!open(unit=2,file='OUTPUT_'//lenght//ext_force//supersaturation//'.out',action ='write',status = 'old')
			write(2,1001) frame,mean_conc,mass,fluxL,mass0+cumul_flux,L,crystal_force
			1001 format(I0,'	',F16.8,'	',F16.8,'	',F16.8,'	',F16.8,'	',F5.1,'	',F16.8)
			
		elseif(s==2) then
			write(*,*) 'Integration finished'
			write(*,*) '******** CHECK *********'
			write(2,*) 'Integration finished'
			write(2,*) '******** CHECK *********'
			curv_L=-(z(M-2) - 2.0d0*z(M-1) + z(M))/(deltax**2) - 1.0d0/L*(z(M)-z(M-1))/deltax
			write (*,*) 'Boundary curvature ~',mu_BC - force(h_BC),'Boundary curvature from profile=', curv_L
			write (2,*) ' Boundary curvature imposed ~',mu_BC- force(h_BC),' boundary curvature from profile=', curv_L
!			write (2,*) 'Boundary curvature from sup ~',sup-nui*growth_vel(M)-force(h_BC),&
!			&'Boundary curvature from profile=', curv_BC
!			write(*,*)'Check: bondary supersaturation =',(curv_BC + force(h_BC))+nui*growth_vel(M),mu(M)+nui*growth_vel(M), sup
!			write (2,*) 'Check: bondary supersaturation =', (curv_BC + force(h_BC))+nui*growth_vel(M), sup
			




			write (2,*) 'h0=', z(1)
			h_max = MAXVAL(z)
			write (2,*) 'hmax-h0=', h_max-z(1)
			write(*,*) 'x_BC =', delta_BC + L
			write(2,*) 'x_BC =', delta_BC + L
			write(2,*) 'Boundary concentration =',c0*(sup+1), &
&'Mean concetration =', mean_conc, 'Conc. in L(last bin)',1+mu(M)+ nui*growth_vel(M) 
			write(*,*) 'Boundary concentrationm =',c0*(sup+1), &
&'Mean concetration =', mean_conc, 'Conc. in L(last bin)',1+mu(M)+ nui*growth_vel(M) 
			write (*,*) 'h0=', z(1)
			write (*,*) '*********************'
			write (2,*) '*********************'

		elseif(s==3) then
			!open(unit=2,file='OUTPUT_'//lenght//ext_force//supersaturation//'.out',action ='write',status='old')
			

			write(2,*) 'CPU running time: ',stop_time-start_time, 'seconds'
			write(*,*) 'CPU running time: ',stop_time-start_time, 'seconds'
			write ( *, 1000 )  today(1), today(2), today(3), now
			write ( 2, 1000 )  today(1), today(2), today(3), now
			write(*,*) 'Elapsed human readable time',today(1)-yesterday(1),'days',now(1)-before(1),'hours',now(2)-before(2),'min'&
						&,now(3)-before(3),'sec'
			write(2,*) 'Elapsed human readable time',today(1)-yesterday(1),'days',now(1)-before(1),'hours',now(2)-before(2),'min'&
			&,now(3)-before(3),'sec'
		endif


	end subroutine write_main_out

    subroutine write_hout(z)

        implicit none
        integer :: u, j
		!integer, intent(in) :: r
		double precision, dimension(M), intent(in) :: z
		character(10) ::  number
		character(15) :: frame
		character(45) :: name
		double precision, dimension (2*M-1) :: dummy

!        write(number,'(I0)'), r
!		number = '_'//number
!		frame = trim(number)
		!f_name = index(frame,' ') - 1

        dummy(1:M) = z(M:1:-1)
        dummy(M:2*M-1) = z(1:M)
		!name = trim('_'//run_name//trim(frame))


		
		write(102,*) -(L+delta_BC),h_BC
		do j = 1, 2*M-1
			write(102, *) ((j-M)*deltax), dummy(j)
		enddo
		write(102,*) (L+delta_BC),h_BC
		write(102,*) '#'

		dummy(1:M) = growth_vel(M:1:-1)
		dummy(M:2*M-1) =  growth_vel(1:M)

		write(103,*) -(L+delta_BC), growthVel_BC
		do j=1,2*M-1	
		write(103,*) ((j-M)*deltax), dummy(j)
		enddo
		write(103,*) (L+delta_BC), growthVel_BC
		write(103,*) '#'

		dummy(1:M) = c_eq(M:1:-1) + nui*growth_vel(M:1:-1)
		dummy(M:2*M-1) = c_eq(1:M)+ nui*growth_vel(1:M)

		write(104,*) -(L+delta_BC),(1.0d0 + mu_BC+nui*growthVel_BC)
		do j=1,2*M-1	
		write(104,*) ((j-M)*deltax), dummy(j)
		enddo
		write(104,*) (L+delta_BC),(1.0d0 + mu_BC+nui*growthVel_BC)
		write(104,*) '#'
		


    end subroutine write_hout
!-------------------------------------



    subroutine write_init()

		implicit none
		integer :: j
		
		double precision, dimension (2*M-1) :: dummy
		INTEGER*4 getcwd, status
		character*160 dirname
		character*10 dir
		logical ex
		
		dri2 = 1./(deltax*deltax)
		
		status = getcwd( dirname )
		fileplace = trim(dirname)//'/plots/'
		print *, fileplace
		
		inquire(file ='plots',DIRECT = dir, EXIST = ex)
		print *, ex
		if (ex .eqv. .true.) then
			print *, 'A folder containing plots already existed, it will be deleted'
			call system('rm -r '//fileplace)
		endif
		call system('mkdir '//fileplace)


		
		
		open(unit = 102, file =trim(fileplace)//'heights.out', action ='write')
		write(102,*) hs,W2
		open(unit =104, file = trim(fileplace)//'concs.out', action = 'write')
		open(unit=103,file=trim(fileplace)//'lGrowthrates.out',action = 'write')


		open (15,file = 'crystal_pressure'//ext_force//supersaturation//'.out',action = 'write')
		write (15,*) '# Contact size estimated seeking for maximum second derivative'
		write (15,*) '# t				contact_size		pressure'

!		open (18,file = 'crystal_pressure_alt'//ext_force//supersaturation//'.out',action = 'write')
!		write (18,*) '# Contact size estimated seeking for maximum third derivative'
!		write (18,*) '# t				contact_size		pressure'
		
!		open (19,file = 'crystal_pressure_alt'//ext_force//supersaturation//'.out',action = 'write')
!		write (19,*) '# Contact size estimated using interaction cut-off'
!		write (19,*) '# t				contact_size		pressure'

		open (16,file = 'crystal_force'//ext_force//supersaturation//'.out',action = 'write')
		write (16,*) '# t				total_size		force'

		open(unit=17,file = 'lateral_vel.out',action = 'write')
		write(17,*) '# t				x_BC		lateral_vel		boundary_der	boundary curvature'

	end subroutine write_init

	subroutine write_final(z)

	implicit none
	integer :: j,k,Ntheta
	double precision :: r1
	double precision, dimension (2*M-1) :: dummy
	double precision, dimension(:), intent(in) :: z
	double precision, dimension (M) :: curv
	
	close(2)
	close(102)
	close(103)
	close(104)
	close(15)
	close(16)
	close(17)
!	close(18)
!	close(19)
	

!	open(unit = 3,file = 'final_2nd_der.out', action = 'write') 
!
!	curv(1) = 1.0/3.0*(0.5d0*z(3)+7.5d0*z(1)-8.0d0*z(2))/(deltax**2)
!	write(3,*) 0,curv(1)
!	do j = 2, M-1
!	r1 = 1.0d0/((j-1) * deltax)
!	curv(j) = -(z(j+1)+z(j-1)-2.0d0*z(j))/(deltax**2)
!	write(3,*) (j-1)*deltax,curv(j)
!	enddo
!	curv(M) =  -(z(M-2)-2*z(M-1)+z(M))/(deltax**2)
!	write(3,*) L,curv(M)
!
!	close(3)
open(unit = 3,file = 'final_profile.out', action = 'write') 

dummy(1:M) = z(M:1:-1)
dummy(M:2*M-1) = z(1:M)
!name = trim('_'//run_name//trim(frame))



write(3,*) -(L+delta_BC),h_BC
do j = 1, 2*M-1
write(3, *) ((j-M)*deltax), dummy(j)
enddo
write(3,*) (L+delta_BC),h_BC

close(3)
open(unit = 3,file = 'final_interaction.out', action = 'write') 

!name = trim('_'//run_name//trim(frame))



write(3,*) -(L+delta_BC),force(h_BC)
do j = 1, 2*M-1
write(3, *) ((j-M)*deltax), force(dummy(j))
enddo
write(3,*) (L+delta_BC),force(h_BC)

close(3)


open(unit = 3,file = 'final_conc.out', action = 'write') 

dummy(1:M) = c_eq(M:1:-1)+ nui*growth_vel(M:1:-1)
dummy(M:2*M-1) = c_eq(1:M)+ nui*growth_vel(1:M)
!name = trim('_'//run_name//trim(frame))



write(3,*) -(L+delta_BC),exp(mu_BC) + nui*growthVel_BC
do j = 1, 2*M-1
write(3, *) ((j-M)*deltax), dummy(j)
enddo
write(3,*) (L+delta_BC),exp(mu_BC) + nui*growthVel_BC

close(3)
!	Ntheta = 50
!	dtheta = 2*PI/Ntheta
!	write(3,*) '#	x			y			z'
!!	write(3,*) 0d0,	0d0,	h(1)
!	do j = 1, M,4
!		write (3,*)
!		r = (j-1)*deltax
!		do k = 1,Ntheta +1
!			theta = (k-1) * dtheta
!			x = r*cos(theta)
!			y = r*sin(theta)
!			write (3,*) x,y,z(j)
!		enddo
!	enddo	
!	close(3)
	
	open(unit = 5, file = run_name//'_exp.lc', action = 'write')!Final configuration to be used as input
	write(5,*) L, deltax,delta_BC
	!write(5,*) curv_BC
	write(5,*) hs
	write(5,*) h_BC
	write(5,*) (z(j),j=1,M)

	close(5)

	

	end subroutine write_final

	subroutine write_crystal_force(z,s)

	double precision, dimension(M), intent(in) :: z
	integer(kind=i15), intent(in) :: s
	double precision, dimension (M) :: curv
	
	
	integer j

	curv(1) = 1.0/3.0*(0.5d0*z(3)+7.5d0*z(1)-8.0d0*z(2))*dri2
	do j = 2, M-1
		curv(j) = -(z(j+1)+z(j-1)-2.0d0*z(j))*dri2
	enddo
	curv(M) =  -(z(M-2)-2*z(M-1)+z(M))*dri2
!	do j = 3, M-2
!		third_der(j) = (0.5*(z(j+2)-z(j-2))-(z(j+1)-z(j-1)))*dri3
!	enddo
!	do j=1,M
!		interaction(j) = abs(force(z(j)))
!	enddo

	L_contact = (maxloc(curv,1)-1)*deltax
	!L_contact_2 = (maxloc(abs(third_der),1)-1)*deltax

	!L_contact_3 = (minloc(abs(interaction - target_value), 1)-1) *deltax
	
	write(15,*)	s*deltat,L_contact,crystal_force/((L_contact)*(L_contact)*PI)
!ALTERNATIVE WAYS OF COMPUTING PRESSURE. NO BETTER PERFORMANCE NOTICED..
	!write(18,*)	s*deltat,L_contact_2,crystal_force/((L_contact_2)*(L_contact_2)*PI)
!	write(19,*)	s*deltat,L_contact_3,crystal_force/((L_contact_3)*(L_contact_3)*PI)

	write(16,*)	s*deltat,(L+delta_BC),crystal_force
	
	write(17,*) s*deltat,(L+delta_BC),lateral_vel,der_BC,curv_BC



	end subroutine write_crystal_force
	
	
subroutine write_3D(s,z)

implicit none
integer, intent(in) :: s
character(10) ::  number
integer :: j,k,Ntheta
double precision theta,r,x,y,dtheta
double precision, dimension(M), intent(in) :: z

write(number,'(I0)'), s
number = '_'//number

open(unit = 115,file =trim(fileplace)//'3d_height'//trim(number)//'.frame', action = 'write')

!do j = 1, 2*M-1
!write (3,*) ((j-M)*deltax), dummy(j) , hs
!enddo
Ntheta = 50
dtheta = 2*PI/Ntheta
write(115,*) '#	x			y			z'
!	write(3,*) 0d0,	0d0,	h(1)
do j = 1, M,4
write (115,*)
r = (j-1)*deltax
do k = 1,Ntheta +1
theta = (k-1) * dtheta
x = r*cos(theta)
y = r*sin(theta)
write (115,*) x,y,z(j)
enddo
enddo	
close(115)





end subroutine write_3D



end MODULE
