subroutine ex_lda_sr(rho_a,rho_b,ex,vx_a,vx_b)
 include 'constants.include.F'
 implicit none 
 double precision, intent(out) ::  ex
 double precision, intent(out) ::  vx_a,vx_b
 double precision, intent(in)  ::  rho_a,rho_b


 double precision ::  rho_a_2,rho_b_2
 double precision :: z0,z1,z2,z3,z4,z6,z8,z16,z24,z96,z12
 double precision :: ex_a,ex_b

 double precision :: f12,f13,f14,f32,f23,f43,f16
 double precision :: ckf 
 double precision :: a, akf,a2, a3

 z0  = 0.D0
 z1  = 1.D0
 z2  = 2.D0
 z3  = 3.D0
 z4  = 4.D0
 z6  = 6.D0
 z8  = 8.D0
 z12 = 12.D0
 z16 = 16.D0
 z24 = 24.D0
 z96 = 96.D0
 f12 = z1/z2
 f13 = z1/z3
 f14 = z1/z4
 f32 = z3/z2
 f23 = z2/z3
 f43 = z4/z3
 f16 = z1/z6
 ckf = (z3*pi*pi)**f13

!Density and kF
 rho_a_2=rho_a*2.D0
 akf = ckf*(rho_a_2**f13)
 a = mu_erf/(z2*akf)
 a2 = a*a
 a3 = a2*a

!Test on the value of a

!Limit for small a (expansion not so important as for large a)
 if (a.lt.1.d-9) then
   ex_a = -z3/z8*rho_a_2*(z24*rho_a_2/pi)**f13
   vx_a = - ((z3/pi)*rho_a_2)**f13

!Intermediate values of a
 elseif (a.le.100d0) then
   ex_a = - (rho_a_2*(z24*rho_a_2/pi)**f13) * (z3/z8-a*(sqpi*derf(f12/a)+(z2*a-z4*a3)*dexp(-f14/a2)-z3*a+z4*a3))
   vx_a =  -(z3*rho_a_2/pi)**f13 + z2*a*mu_erf/pi*(dexp(-f14/a2)-z1)+mu_erf/sqpi * derf(f12/a)


!Expansion for large a
 elseif (a.lt.1.d+9) then
   ex_a = -(rho_a_2*(z24*rho_a_2/pi)**f13) * z1/(z96*a2)
   vx_a = -pi*rho_a_2/(z2*mu_erf*mu_erf)

!Limit for large a
 else
   ex_a = 0.d0
   vx_a = 0.d0
 end if

!Density and kF
 rho_b_2= rho_b * 2.d0
 akf = ckf*(rho_b_2**f13)
 a = mu_erf/(z2*akf)
 a2 = a*a
 a3 = a2*a

!Test on the value of a

!Limit for small a (expansion not so important as for large a)
 if (a.lt.1.d-9) then
   ex_b = -z3/z8*rho_b_2*(z24*rho_b_2/pi)**f13
   vx_b = - ((z3/pi)*rho_b_2)**f13

!Intermediate values of a
 elseif (a.le.100d0) then
   ex_b = - (rho_b_2*(z24*rho_b_2/pi)**f13)*(z3/z8-a*(sqpi*derf(f12/a)+(z2*a-z4*a3)*dexp(-f14/a2)-z3*a+z4*a3))
   vx_b = -(z3*rho_b_2/pi)**f13+ z2*a*mu_erf/pi*(dexp(-f14/a2)-z1)+mu_erf/sqpi* derf(f12/a)

!Expansion for large a
 elseif (a.lt.1.d+9) then
   ex_b = - (rho_b_2*(z24*rho_b_2/pi)**f13) *z1/(z96*a2)
   vx_b = - pi*rho_b_2/(z2*mu_erf*mu_erf)

!Limit for large a
 else
   ex_b = z0
   vx_b = 0.d0
 end if

 ex = (ex_a+ex_b) * 0.5d0



end
