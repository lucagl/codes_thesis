MODULE com_var

implicit none


INTERFACE
function func(x)
double precision, intent(in) ::x
double precision :: func
end function func
end INTERFACE

INTERFACE
function int_func(x)
integer, intent(in) ::x
double precision :: int_func
end function int_func
end INTERFACE

!INTERFACE
!function integration(z)
!double precision, dimension(:), intent(in) :: z
!double precision, dimension(size(z)):: integration
!end function integration
!end INTERFACE

!INTERFACE
!subroutine initialize_integration()
!end subroutine
!end INTERFACE




double precision, parameter :: PI = 3.14159265358979
integer, parameter :: i15 = selected_int_kind(15)
integer(kind= i15) :: N,freq_step,freq_vel
integer ::  M !time, M space
integer, parameter :: M_max = 10000

double precision :: lateral_vel,L,h_BC,delta_BC,curv_BC

double precision ::  deltax, deltat, sup,c0, gamma,hs
character(len=:), allocatable :: run_name,lenght,supersaturation,ext_force
double precision ::  B, dinv,hc,nu,nui
double precision, dimension(:),allocatable::  c_eq
double precision, dimension(:),allocatable:: mu, Q,growth_vel

double precision :: mean_conc,mass0,mass,fluxL,cumul_flux
double precision :: start_time, stop_time
integer *4 :: yesterday(4),today(4), before(3),now(3)
double precision :: growthVel_BC,mu_BC
double precision :: crystal_force

double precision :: der_BC

procedure (func), pointer :: force, force_1der, force_2der,U



SAVE

end MODULE
