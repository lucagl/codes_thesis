
!-------------------------------------------------------------------------------------------
!
!
!				CONTAINS ANALYTICAL FORCES FROM GIVEN POTENTIAL
!
!
!-------------------------------------------------------------------------------------------



MODULE get_force

use com_var, only : B,dinv,hc,hs
!use initialization, only : B,d,dinv,hc
implicit none


CONTAINS

	function U_eff(x)
	implicit none
	double precision, intent(in) :: x
	double precision:: U_eff

		U_eff = B * (dexp(-(hc-x)*dinv)/(hc-x))

	end function U_eff


    function f_eff(x)
        implicit none
        double precision, intent(in) :: x
        double precision :: f_eff
		double precision :: zita_ci, zita_ci2
		
		zita_ci = 1.0d0/(hc-x)
		zita_ci2 = 1.0d0/((hc-x)*(hc-x))
		
        f_eff =  B * (dinv*zita_ci + zita_ci2)*dexp(-(hc-x)*dinv)

    end function f_eff

function U_DLVO(x)
implicit none
double precision, intent(in) :: x
double precision:: U_DLVO

U_DLVO = dexp(-(hs-x)*dinv)

end function U_DLVO


function f_DLVO(x)
implicit none
double precision, intent(in) :: x
double precision :: f_DLVO


f_DLVO =  dinv * dexp(-(hs-x)*dinv)

end function f_DLVO


function U_PL(x)
implicit none
double precision, intent(in) :: x
double precision:: U_PL

U_PL = B/(3.0d0 * (hs-x)*(hs-x)*(hs-x)) - 1.0d0/(2.0d0*(hs-x)*(hs-x))

end function U_PL

function fPL(x)
implicit none
double precision, intent(in) :: x
double precision:: fPL,x3,x4

x4 = 1.0d0/((hs-x)*(hs-x)*(hs-x)*(hs-x))
x3 = 1.0d0/((hs-x)*(hs-x)*(hs-x))
fPL = B*x4 - x3

end function fPL

function U_PL_RP(x)
implicit none
double precision, intent(in) :: x
double precision:: U_PL_RP

U_PL_RP = B/(3.0d0 * (hs-x)*(hs-x)*(hs-x))

end function U_PL_RP

function fPL_RP(x)
implicit none
double precision, intent(in) :: x
double precision:: fPL_RP,x4

x4 = 1.0d0/((hs-x)*(hs-x)*(hs-x)*(hs-x))
fPL_RP = B*x4

end function fPL_RP

function f_0(x)
implicit none
double precision, intent(in) :: x
double precision :: f_0
f_0 = 0.0d0

end function f_0

function U_0(x)
implicit none
double precision, intent(in) :: x
double precision :: U_0
U_0 = 0.0d0

end function U_0


end MODULE

