!-------------------------------------------------------------------------------------------
!
!
!				CONTAINS ALL GENERAL INITIALISATION ROUTINES
!
!
!
!-------------------------------------------------------------------------------------------

MODULE initialization
use com_var
use read_input, only : read_potential,hs_overwrite,d_read,get_sup_name,h0
use get_force
!use integrate

implicit none

double precision, private :: h00


CONTAINS

!--------------------

subroutine Init(s)
implicit none
integer, intent(in) ::s
integer j,i
character (len=30) :: dummy_string

if(s==0) then
	open(unit=2,file='OUTPUT'//lenght//supersaturation//'.out',action ='write')
	write(*,*) '+++++++++ Initialising force ++++++++'
	write(2,*) '+++++++++ Initialising force ++++++++'

	write(2,*) ''
	if (read_potential==1) then
		write(2,'(A)') ('Potential used LDVO (only repulsion)')
		write(*,'(A)') ('Potential used LDVO (only repulsion)')
		write (*,*) 'Debyle lenght constant =', d_read
		write (2,*) 'Debyle lenght constant =', d_read

		call init_DLVO_RP() !Initialize potential constant and hs
		
		elseif (read_potential==2) then
		write(2,'(A)') ('Potential used: single well')
		write(*,'(A)') ('Potential used: single well')
		call init_PL()
		
		elseif (read_potential ==3) then
			write(2,'(A)') ('Pure power law repulsion witj substrate')
			write(*,'(A)') ('Pure power law repulsion witj substrate')	
			call init_PL_RP()
			write (*,*) 'Position of interface =', hs
			write (2,*) 'Position of interface =', hs
		
		elseif (read_potential ==4) then
			write(2,'(A)') ('Effective interaction')
			write(*,'(A)') ('Effective interaction')
			call init_eff()
			write (*,*) 'Position of substrate =', hs
			write (2,*) 'Position of substrate =', hs
			write (*,*) 'Position of dust grain =', hc
			write (2,*) 'Position of dust grain =', hc
			write (*,*) 'Debyle lenght constant =', d_read
			write (2,*) 'Debyle lenght constant =', d_read
			write (*,*) 'Potential amplitude =', B
			write (2,*) 'Potential amplitude =', B
		elseif (read_potential ==5) then
			write(2,'(A)') ('No interaction with substrate')
			write(*,'(A)') ('No interaction with substrate')
			force => f_0
			U =>U_0
			h00 = 1.0d0
		endif
	
!!!! NOT IMPLEMENTED HERE !!!!!
elseif (s==1) then
	allocate(h0(M))
	allocate(c_eq(M_max))
	open(unit=2,file='OUTPUT'//lenght//supersaturation//'.out',action ='write')

	print *, 'Using zero curvature eq der hyperbolic tangent'
	write (2,*) 'Using zero curvature eq der hyperbolic tangent'

	

	 

	!c_inf = (1+0.0d0) + force(h_BC)!(1+gamma*curv_BC) + force(hs-h_BC)
	sup = 0!zero supersaturation 

!	do j =1, M
!		h0(j) = Initfunc(j-1) 
!	enddo

	
!	print *, 'h_BC =', h_BC, 'fmax=', force(h_BC)
!	print *,'der_BC =', -der_BC, ((h0(M)-h0(M-1))/deltax)
!	print *, 'h_BC=', h0(M)
!	print *, 'hs=', hs
	
	write(dummy_string,'(F0.2)') sup
	call get_sup_name(dummy_string)

	
!	print *,  'A new initial configuration was generated'
!	write (2,*) 'A new initial configuration was generated'


endif


end subroutine Init

!--------------------------------------------


!function Initfunc(i) result(Initfunc_hyp)
!
!implicit none
!integer, intent(in) ::  i
!double precision :: Initfunc_hyp
!double precision :: lambda,x0,w!parameters used for generating initial configuration
!
!der_BC =  sqrt(-2.0d0 *U(h00))
!x0 = ((h_BC-h00)/der_BC) + L
!w = 1.0d0!sqrt(-1.0d0/force_1der(h00))
!Initfunc_hyp= h00 +&
!&0.5d0*((1.0d0-dexp(2.0d0*(i*deltax -x0)/w))/(1.0d0+dexp(2.0d0*(i*deltax-x0)/w))-1)*((i*deltax -L)*der_BC+h00-h_BC)
!
!
!end function Initfunc



!-------------------------


subroutine init_eff()
implicit none
dinv = 1.0/d_read
!hc = hs  !not clean, but just to use value hs from old equilibrated simulations as value of hc
hs = hs_overwrite
!no need to define h00 since this is used to generate initial configuration which is not implemented for repulsive potentials
U => U_eff
force => f_eff!assign correct pointer to function

end subroutine init_eff

subroutine init_DLVO_RP()
implicit none
dinv = 1.0/d_read
h00 =1.0d0
U => U_DLVO
force => f_DLVO!assign correct pointer to function

end subroutine init_DLVO_RP

!-----------------------------
!-----------------------------

subroutine init_PL()
implicit none

U => U_PL
force=> fPL!assign correct pointer to function
B = 1.0d0 
h00 = hs - B!put central distance at minimum of potential if generating new configuration


end subroutine init_PL

subroutine init_PL_RP()
implicit none

U => U_PL_RP
force=> fPL_RP!assign correct pointer to function
B = 1.0d0 
hs = hs_overwrite


end subroutine init_PL_RP



!-----------------------------
!-----------------------------



end MODULE
