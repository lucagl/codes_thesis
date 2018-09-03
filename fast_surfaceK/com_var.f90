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
integer :: L, M !time, M space
double precision ::  deltax, deltat, h_BC, c_inf, der_BC,c0, gamma,hs
character(len=:), allocatable :: run_name,lenght,concentration,ext_force
double precision ::  B, dinv,hc
double precision, dimension(:),allocatable:: h0, c_eq
double precision, dimension(:),allocatable:: h, mu, Q,pressure
double precision :: load, viscosity
double precision :: meanc_eq,mass0,mass,fluxL, cumul_flux,crystal_vel,intforce,friction
double precision :: start_time, stop_time
integer *4 :: yesterday(4),today(4), before(3),now(3)

procedure (func), pointer :: force, force_1der, force_2der,U



SAVE

end MODULE
