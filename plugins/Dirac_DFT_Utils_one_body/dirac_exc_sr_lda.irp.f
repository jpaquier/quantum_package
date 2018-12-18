 subroutine dirac_ex_lda(rho,ex,vx)
 include 'constants.include.F'
 implicit none
 double precision, intent(in) :: rho
 double precision, intent(out) :: ex,vx
 double precision :: tmp
 tmp = rho**(c_1_3)
 ex = cx_lda * tmp*tmp*tmp*tmp
 vx = cx_lda * c_4_3 * tmp
 end


!subroutine dirac_ec_lda(rho_a,rho_b,ec,vc_a,vc_b)
!     implicit none
!include 'constants.include.F'
!     double precision, intent(out) ::  ec
!     double precision, intent(out) ::  vc_a,vc_b
!     double precision, intent(in)  ::  rho_a,rho_b
!! Double precision numbers
!     double precision :: rsfac,rho,rs,rhoa,rhob,z
!     double precision :: eccoul, ecd, ecz, ecdd, eczd
!     double precision :: vcup,vcdown
!     rsfac = (3.0d0/(4.0d0*pi))**c_1_3
!! Test on density
!     rho = rho_a + rho_b
!     if (dabs(rho).ge.1.d-10) then
!     rs=rsfac/(rho**c_1_3)
!     rhoa=max(rho_a,1.0d-15)
!     rhob=max(rho_b,1.0d-15)
!     z=(rhoa-rhob)/(rhoa+rhob)
!     call ecPW(rs,z,eccoul,ecd,ecz,ecdd,eczd)
!     ec=(eccoul)*rho
!     vcup=eccoul-rs/3.d0*ecd-(z-1.d0)*ecz
!     vcdown=eccoul-rs/3.d0*ecd-(z+1.d0)*ecz
!     vc_a = vcup
!     vc_b = vcdown
!     else
!      ec = 1.d-15
!      vc_a = 1.d-15
!      vc_b = 1.d-15
!     endif
!end

!subroutine dirac_ec_lda_sr(mu,rho_a,rho_b,ec,vc_a,vc_b)
!     implicit none
!include 'constants.include.F'
!     double precision, intent(out) ::  ec
!     double precision, intent(out) ::  vc_a,vc_b
!     double precision, intent(in)  ::  mu,rho_a,rho_b
!! Double precision numbers
!     double precision :: rsfac,rho,rs,rhoa,rhob,z
!     double precision :: eccoul, ecd, ecz, ecdd, eczd
!     double precision :: eclr,vcup,vcdown,vclrup,vclrdown,vclrupd,vclrdownd
!     rsfac = (3.0d0/(4.0d0*pi))**c_1_3
!     ec = 0.d0
!     vc_a = 0.d0
!     vc_b = 0.d0
!! Test on density
!     rho = rho_a + rho_b
!     if (dabs(rho).ge.1.d-12) then
!     rs=rsfac/(rho**c_1_3)
!     rhoa=max(rho_a,1.0d-15)
!     rhob=max(rho_b,1.0d-15)
!     z=(rhoa-rhob)/(rhoa+rhob)
!     call ecPW(rs,z,eccoul,ecd,ecz,ecdd,eczd)
!     call ecorrlr(rs,z,mu,eclr)
!     ec=(eccoul-eclr)*rho
!     vcup=eccoul-rs/3.d0*ecd-(z-1.d0)*ecz
!     vcdown=eccoul-rs/3.d0*ecd-(z+1.d0)*ecz
!     call vcorrlr(rs,z,mu,vclrup,vclrdown,vclrupd,vclrdownd)
!     vc_a = vcup-vclrup
!     vc_b = vcdown-vclrdown
!     else
!      ec = 1.d-15
!      vc_a = 1.d-15
!      vc_b = 1.d-15
!     endif
!end

!subroutine dirac_ex_lda_sr(mu,rho_a,rho_b,ex,vx_a,vx_b)
!include 'constants.include.F'
!implicit none 
!double precision, intent(out) ::  ex
!double precision, intent(out) ::  vx_a,vx_b
!double precision, intent(in)  ::  rho_a,rho_b,mu
!double precision ::  rho_a_2,rho_b_2
!double precision :: z0,z1,z2,z3,z4,z6,z8,z16,z24,z96,z12
!double precision :: ex_a,ex_b
!double precision :: f12,f13,f14,f32,f23,f43,f16
!double precision :: ckf 
!double precision :: a, akf,a2, a3
!z0  = 0.D0
!z1  = 1.D0
!z2  = 2.D0
!z3  = 3.D0
!z4  = 4.D0
!z6  = 6.D0
!z8  = 8.D0
!z12 = 12.D0
!z16 = 16.D0
!z24 = 24.D0
!z96 = 96.D0
!f12 = 0.5d0
!f13 = 0.3333333333333333d0
!f14 = 0.25d0
!f32 = 1.5d0
!f23 = 0.6666666666666666d0
!f43 = 1.3333333333333333d0
!f16 = 0.16666666666666666d0
!ckf = 3.0936677262801355d0
!!Density and kF
!rho_a_2=rho_a*2.D0
!akf = ckf*(rho_a_2**f13)
!a = mu/(z2*akf)
!a2 = a*a
!a3 = a2*a
!!Test on the value of a
!!Limit for small a (expansion not so important as for large a)
!if (a.lt.1.d-9) then
!  ex_a = -z3/z8*rho_a_2*(z24*rho_a_2/pi)**f13
!  vx_a = - ((z3/pi)*rho_a_2)**f13
!!Intermediate values of a
!elseif (a.le.100d0) then
!  ex_a = - (rho_a_2*(z24*rho_a_2/pi)**f13) * (z3/z8-a*(sqpi*derf(f12/a)+(z2*a-z4*a3)*dexp(-f14/a2)-z3*a+z4*a3))
!  vx_a =  -(z3*rho_a_2/pi)**f13 + z2*a*mu/pi*(dexp(-f14/a2)-z1)+mu/sqpi * derf(f12/a)
!!Expansion for large a
!elseif (a.lt.1.d+9) then
!  ex_a = -(rho_a_2*(z24*rho_a_2/pi)**f13) * z1/(z96*a2)
!  vx_a = -pi*rho_a_2/(z2*mu*mu)
!!Limit for large a
!else
!  ex_a = 0.d0
!  vx_a = 0.d0
!end if
!!Density and kF
!rho_b_2= rho_b * 2.d0
!akf = ckf*(rho_b_2**f13)
!a = mu/(z2*akf)
!a2 = a*a
!a3 = a2*a
!!Test on the value of a
!!Limit for small a (expansion not so important as for large a)
!if (a.lt.1.d-9) then
!  ex_b = -z3/z8*rho_b_2*(z24*rho_b_2/pi)**f13
!  vx_b = - ((z3/pi)*rho_b_2)**f13
!!Intermediate values of a
!elseif (a.le.100d0) then
!  ex_b = - (rho_b_2*(z24*rho_b_2/pi)**f13)*(z3/z8-a*(sqpi*derf(f12/a)+(z2*a-z4*a3)*dexp(-f14/a2)-z3*a+z4*a3))
!  vx_b = -(z3*rho_b_2/pi)**f13+ z2*a*mu/pi*(dexp(-f14/a2)-z1)+mu/sqpi* derf(f12/a)
!!Expansion for large a
!elseif (a.lt.1.d+9) then
!  ex_b = - (rho_b_2*(z24*rho_b_2/pi)**f13) *z1/(z96*a2)
!  vx_b = - pi*rho_b_2/(z2*mu*mu)
!!Limit for large a
!else
!  ex_b = z0
!  vx_b = 0.d0
!end if
! ex = (ex_a+ex_b) * 0.5d0
!end
