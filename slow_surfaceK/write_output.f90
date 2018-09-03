MODULE  write_output!contains also time information

use com_var
use read_input, only : read_boundary,read_last
!use ifport

implicit none
double precision, private :: curv_BC
CHARACTER*160, private ::  fileplace


CONTAINS

	subroutine write_main_out(s,frame)
		implicit none 
		integer, intent(in) :: s
		integer,optional :: frame	
		double precision :: h_max
		!if(present(frame)) fr_number = frame
		if (s == 0) then
			open(unit=2,file='OUTPUT_'//lenght//ext_force//supersaturation//'.out',action ='write')
			write (2,*) ''
			write (2,*) ''
			write ( 2, 1000 )  yesterday(1), yesterday(2), yesterday(3), before
			1000 format ( 'Date ', i2.2, '/', i2.2, '/', i4.4, '; time ',i2.2, ':', i2.2, ':', i2.2 )
			write (2,*) ''
			write (2,*) ''
			write (2,'(A,A)') trim(run_name),  '       OUTPUT      '
			write (2,'(A,3L)') 'Reading from last configuration   ', read_last
			write (2,*) 'Maximum number of time steps  ', N
			write (2,*) 'Load =', load
			write (2,*) 'Viscosity =', viscosity
			write (2,*) 'Delta t  ', deltat
			write (2,*) 'Frame output frequency  ', freq_step
			write (2,*) 'Lenght of the box   ', L
			write (2,*) 'Delta x  ', deltax
			write (2,*) 'Kinetic constant ', nu


			curv_BC=-(h(M-2) - 2.0d0*h(M-1) + h(M))/(deltax**2)- 1.0d0/L*(h(M)-h(M-1))/deltax !2d curvature
			write (2,*) 'Initial boundary height and curvature:   ','h_BC= ', h_BC,'curv_BC= ', curv_BC
			write (2,*) 'Boundary supersaturation', sup, 'force at the boundary = ', force(h_BC)
			write (2,*) 'Boundary derivative', -der_BC, +((h(M)-h(M-1))/deltax) !- 0.5* curv_BC *deltax)
			write (2,*) 'Height of the substrate hs = ', hs
			write(2,*) ''
			write(2,*) ''
			write(2,*) 'frame          mean_conc               mass	   	      flux(L)	       M(0) + cumul_flux(L)		 velocity		cryst_force'
		
		elseif(s==1) then
			!open(unit=2,file='OUTPUT_'//lenght//ext_force//supersaturation//'.out',action ='write',status = 'old')
			write(2,1001) frame,mean_conc,mass,-fluxL,mass0+cumul_flux,crystal_vel,2.0d0*PI*intforce
			1001 format(I0,'	',F16.8,'	',F16.8,'	',F16.8,'	',F16.8,'	',F16.8,'	',F16.8)
			
		elseif(s==2) then
	
			write(*,*) 'Integration finished'
			write(*,*) '******** CHECK *********'
			write(2,*) 'Integration finished'
			write(2,*) '******** CHECK *********'
			write(2,*)  'Final friction term:', friction			
			curv_BC=-(h(M-2) - 2.0d0*h(M-1) + h(M))/(deltax**2) - 1.0d0/L*(h(M)-h(M-1))/deltax
			write (*,*) ' boundary curvature from sup ~',sup-nui*growth_vel(M)-force(h_BC),mu(M),' boundary curvature from profile=', curv_BC
			write (2,*) ' boundary curvature from sup ~',sup-nui*growth_vel(M)-force(h_BC),' boundary curvature from profile=', curv_BC
			!write (2,*) 'initial boundary curvature=',curv_BC,'final boundary curvature=', -(h(M-2) - 2.0d0*h(M-1) + h(M))/(deltax**2)
			 !not the most refinite scheme for the boundaries..
			write(*,*)'Check: bondary supersaturation =',(curv_BC + force(h_BC))+nui*growth_vel(M),mu(M)+nui*growth_vel(M), sup
			write (2,*) 'Check: bondary supersaturation =', (curv_BC + force(h_BC))+nui*growth_vel(M), sup
			!write(2,*) 'Final step=', i-1
			write (2,*) 'h0=', h(1)
			h_max = MAXVAL(h)
			write (2,*) 'hmax-h0=', h_max-h(1)
			!write (2,*) 'Force in the center', force(hs-h(1)),'Boundary curvature + force', curv_BC + force(hs-h_BC)
			!write (2,*) 'Force at the boundary = ', force(hs-h_BC)

			write(2,*) 'Boundary concentration =',c0*(sup+1), 'Mean concentration =', mean_conc, 'difference=', abs(c0*(sup+1)-mean_conc)
			!write(2,*) 'Boundary derivative = ', der_BC, 'equilibrium infinite crystal value =',sqrt(-2.0d0/gamma *U(hs - h00))

			write(*,*) 'Boundary concentration =',c0*(sup+1), 'Mean concentration =', mean_conc, 'difference=', abs(c0*(sup+1)-mean_conc)
			!write(*,*) 'Final step=', i-1
			write (*,*) 'h0=', h(1)
			write (*,*) '*********************'
			write (2,*) '*********************'

		elseif(s==3) then
			

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

    subroutine write_hout()! h is the current height,concentration,interaction force and pressure

        implicit none
        integer :: u, j
!		integer, intent(in) :: r
		character(10) ::  number
		character(15) :: frame
		character(45) :: name
		double precision, dimension (2*M-1) :: dummy

!        write(number,'(I0)'), r
!		number = '_'//number
!		frame = lenght//trim(number)
!		!f_name = index(frame,' ') - 1


		!name = trim('_'//run_name//'_'//trim(frame))
		dummy(1:M) = h(M:1:-1)
		dummy(M:2*M-1) = h(1:M)

		do j = 1, 2*M-1
		write(102, *) ((j-M)*deltax), dummy(j)
		enddo
		write(102,*) '#'

		dummy(1:M) = growth_vel(M:1:-1)
		dummy(M:2*M-1) =  growth_vel(1:M)

		do j=1,2*M-1	
		write(103,*) ((j-M)*deltax), dummy(j)
		enddo
		write(103,*) '#'

		dummy(1:M) = c_eq(M:1:-1) + nui*growth_vel(M:1:-1)
		dummy(M:2*M-1) = c_eq(1:M)+ nui*growth_vel(1:M)

		do j=1,2*M-1	
		write(104,*) ((j-M)*deltax), dummy(j)
		enddo
		write(104,*) '#'




		dummy(1:M) = pressure(M:1:-1)
		dummy(M:2*M-1) = pressure(1:M)

		do j=1,2*M-1	
		write(105,*) ((j-M)*deltax), dummy(j)
		enddo
		write(105,*) '#'

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
		
	
		status = getcwd( dirname )
		fileplace = trim(dirname)//'/plots/'
		!print *, dirname
		print *, fileplace
		
		inquire(file ='plots',DIRECT = dir, EXIST = ex)
		print *, ex
		!print *,dir
		if (ex .eqv. .true.) then
			print *, 'A folder containing plots already existed, it will be deleted'
			call system('rm -r '//fileplace)
		endif
		call system('mkdir '//fileplace)

		open (15,file = 'crystal_velocity'//ext_force//lenght//supersaturation//'.out',action = 'write')
		write (15,*) '# t		vel'
		open(unit=16,file = 'zita0'//ext_force//'.out',action = 'write')
		write (16,*) '# t		zita0'


		open(unit = 102, file =trim(fileplace)//'heights.out', action ='write')
		open(unit =105, file = trim(fileplace)//'pressures.out', action = 'write')
		open(unit =104, file = trim(fileplace)//'concs.out', action = 'write')
		open(unit=103,file=trim(fileplace)//'lGrowthrates.out',action = 'write')

	end subroutine write_init

	subroutine write_final()

	implicit none
	integer :: j,k,Ntheta
	double precision theta,r,x,y,dtheta
	double precision, dimension (2*M-1) :: dummy

	close(2)
	close(102)
	close(103)
	close(104)
	close(105)
	close(15)
	close(16)
	close(17)

	
	open(unit = 3,file = 'final_height.out', action = 'write')!write final surface in 3d!!
		dummy(1:M) = h(M:1:-1)
		dummy(M:2*M-1) = h(1:M)

		do j = 1, 2*M-1
		write(3, *) ((j-M)*deltax), dummy(j)
		enddo

	close(3)
	
	open(unit = 3,file = 'final_conc.out', action = 'write')!write final surface in 3d!!
	dummy(1:M) = c_eq(M:1:-1) + nui*growth_vel(M:1:-1)
	dummy(M:2*M-1) = c_eq(1:M)+ nui*growth_vel(1:M)

	do j = 1, 2*M-1
	write(3, *) ((j-M)*deltax), dummy(j)
	enddo

	close(3)
	open(unit = 5, file = run_name//'.lc', action = 'write')!Final configuration to be used as input
	write(5,*) L, deltax
	write(5,*) hs
	der_BC = -((h(M)-h(M-1))/deltax)
	write(5,*) der_BC
	write(5,*) (h(j),j=1,M)
	
	close(5)

	

	end subroutine write_final

	subroutine write_crystal_vel(s)

	integer(kind=i15), intent(in) :: s

	

	write(15,*)	s*deltat,crystal_vel
	
	write(16,*) s*deltat,(hs-h(1))
	

	end subroutine write_crystal_vel
	
	
subroutine write_3D(s)

implicit none
integer, intent(in) :: s
character(10) ::  number
integer :: j,k,Ntheta
double precision theta,r,x,y,dtheta

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
write (115,*) x,y,h(j)
enddo
enddo	
close(115)





end subroutine write_3D



end MODULE
