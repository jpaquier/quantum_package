subroutine ex_lda_sr(mu,rho_a,rho_b,ex,vx_a,vx_b)
 include 'constants.include.F'
 implicit none 
 double precision, intent(out) ::  ex
 double precision, intent(out) ::  vx_a,vx_b
 double precision, intent(in)  ::  rho_a,rho_b,mu


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
 f12 = 0.5d0
 f13 = 0.3333333333333333d0
 f14 = 0.25d0
 f32 = 1.5d0
 f23 = 0.6666666666666666d0
 f43 = 1.3333333333333333d0
 f16 = 0.16666666666666666d0
 ckf = 3.0936677262801355d0

!Density and kF
 rho_a_2=rho_a*2.D0
 akf = ckf*(rho_a_2**f13)
 a = mu/(z2*akf)
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
   vx_a =  -(z3*rho_a_2/pi)**f13 + z2*a*mu/pi*(dexp(-f14/a2)-z1)+mu/sqpi * derf(f12/a)


!Expansion for large a
 elseif (a.lt.1.d+9) then
   ex_a = -(rho_a_2*(z24*rho_a_2/pi)**f13) * z1/(z96*a2)
   vx_a = -pi*rho_a_2/(z2*mu*mu)

!Limit for large a
 else
   ex_a = 0.d0
   vx_a = 0.d0
 end if

!Density and kF
 rho_b_2= rho_b * 2.d0
 akf = ckf*(rho_b_2**f13)
 a = mu/(z2*akf)
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
   vx_b = -(z3*rho_b_2/pi)**f13+ z2*a*mu/pi*(dexp(-f14/a2)-z1)+mu/sqpi* derf(f12/a)

!Expansion for large a
 elseif (a.lt.1.d+9) then
   ex_b = - (rho_b_2*(z24*rho_b_2/pi)**f13) *z1/(z96*a2)
   vx_b = - pi*rho_b_2/(z2*mu*mu)

!Limit for large a
 else
   ex_b = z0
   vx_b = 0.d0
 end if

 ex = (ex_a+ex_b) * 0.5d0

end

subroutine ex_pbe_sr(mu,rho_a,rho_b,grd_rho_a_2,grd_rho_b_2,grd_rho_a_b,ex,vx_rho_a,vx_rho_b,vx_grd_rho_a_2,vx_grd_rho_b_2,vx_grd_rho_a_b)
BEGIN_DOC
!rho_a = density alpha
!rho_b = density beta
!grd_rho_a_2 = (gradient rho_a)^2
!grd_rho_b_2 = (gradient rho_b)^2
!grd_rho_a_b = (gradient rho_a).(gradient rho_b)
!ex = exchange energy density at point r
!vx_rho_a = d ex / d rho_a
!vx_rho_b = d ex / d rho_b
!vx_grd_rho_a_2 = d ex / d grd_rho_a_2
!vx_grd_rho_b_2 = d ex / d grd_rho_b_2
!vx_grd_rho_a_b = d ex / d grd_rho_a_b
END_DOC

 implicit none

! input
 double precision, intent(in) :: mu,rho_a, rho_b
 double precision, intent(in) :: grd_rho_a_2, grd_rho_b_2, grd_rho_a_b

! output
 double precision, intent(out) :: ex
 double precision, intent(out) :: vx_rho_a, vx_rho_b
 double precision, intent(out) :: vx_grd_rho_a_2, vx_grd_rho_b_2, vx_grd_rho_a_b

! function
  double precision berf
  double precision dberfda

! local
  double precision, parameter :: tol=1d-12
  double precision, parameter :: f13=0.333333333333333d0

  double precision exerflda,vxerflda_a,vxerflda_b
  double precision dexerfldadrho
  double precision exerfpbe_a, exerfpbe_b
  double precision dexerfpbedrho_a, dexerfpbedrho_b
  double precision dexerfpbeddrho2_a, dexerfpbeddrho2_b

  double precision rho,drho2
  double precision rho_a_2, rho_b_2
  double precision t1,t2,t3,t4
  double precision kappa,sq,sqs,sqss,fx,fxs,ksig

! Parameter of the modified interaction

! initialization
  ex=0.d0
  vx_rho_a=0.d0
  vx_rho_b=0.d0
  vx_grd_rho_a_2=0.d0
  vx_grd_rho_b_2=0.d0
  vx_grd_rho_a_b=0.d0

  
! spin scaling relation Ex[rho_a,rho_b] = (1/2) (Ex[2rho_a,2rho_a] + Ex[2rho_b,2rho_b])

! two times spin alpha density
  rho = max(rho_a,tol)*2.d0

! test on density
  if (rho >= tol) then

!  call srLDA Ex[2*rho_a,2*rho_a]
   call ex_lda_sr(mu,rho_a,rho_a,exerflda,vxerflda_a,vxerflda_b)
   dexerfldadrho = (vxerflda_a + vxerflda_b)*0.5d0

!  square of two times spin alpha density gradient
   drho2=max(grd_rho_a_2,0d0)*4.0d0

   kappa=0.804d0
   sq=drho2*2.6121172985233599567768d-2*rho**(-8d0/3d0)
   fx=1d0+kappa-kappa/(1d0+berf(1.616204596739954813d-1*mu*rho**(-f13))*sq/kappa)
   exerfpbe_a=exerflda*fx

!  Derivatives
   sqs=-8d0*sq/(3d0*rho)
   fxs=kappa**2*(-1.616204596739954813d-1*mu*rho**(-4d0*f13)/3d0*dberfda(1.616204596739954813d-1*mu*rho**(-f13))*sq+berf(1.616204596739954813d-1*mu*rho**(-f13))*sqs)/(kappa+berf(1.616204596739954813d-1*mu*rho**(-f13))*sq)**2
   dexerfpbedrho_a=dexerfldadrho*fx+exerflda*fxs
   sqss=2.6121172985233599567768d-2*rho**(-8d0/3d0)
   dexerfpbeddrho2_a=exerflda*berf(1.616204596739954813d-1*mu*rho**(-1.d0/3.d0))*sqss*kappa**2/(kappa+berf(1.616204596739954813d-1*mu*rho**(-1.d0/3.d0))*sq)**2

 endif
   

! two times spin beta density
  rho = max(rho_b,tol)*2.d0

! test on density
  if (rho >= tol) then

!  call srLDA Ex[2*rho_b,2*rho_b]
   call ex_lda_sr(mu,rho_b,rho_b,exerflda,vxerflda_a,vxerflda_b)
   dexerfldadrho = (vxerflda_a + vxerflda_b)*0.5d0

!  square of two times spin beta density gradient
   drho2=max(grd_rho_b_2,0d0)*4.0d0

   kappa=0.804d0
   sq=drho2*2.6121172985233599567768d-2*rho**(-8d0/3d0)
   fx=1d0+kappa-kappa/(1d0+berf(1.616204596739954813d-1*mu*rho**(-f13))*sq/kappa)
   exerfpbe_b=exerflda*fx

!  Derivatives
   sqs=-8d0*sq/(3d0*rho)
   fxs=kappa**2*(-1.616204596739954813d-1*mu*rho**(-4d0*f13)/3d0*dberfda(1.616204596739954813d-1*mu*rho**(-f13))*sq+berf(1.616204596739954813d-1*mu*rho**(-f13))*sqs)/(kappa+berf(1.616204596739954813d-1*mu*rho**(-f13))*sq)**2
   dexerfpbedrho_b=dexerfldadrho*fx+exerflda*fxs
   sqss=2.6121172985233599567768d-2*rho**(-8d0/3d0)
   dexerfpbeddrho2_b=exerflda*berf(1.616204596739954813d-1*mu*rho**(-1.d0/3.d0))*sqss*kappa**2/(kappa+berf(1.616204596739954813d-1*mu*rho**(-1.d0/3.d0))*sq)**2

  endif


  ex = (exerfpbe_a+exerfpbe_b)*0.5d0
  vx_rho_a = dexerfpbedrho_a 
  vx_rho_b = dexerfpbedrho_a 
  vx_grd_rho_a_2 =   2.d0*dexerfpbeddrho2_a
  vx_grd_rho_b_2 =   2.d0*dexerfpbeddrho2_b
  vx_grd_rho_a_b = 0.d0

  end

!-------------------------------------------
      function berf(a)
!-------------------------------------------
!  Second-order exchange gradient expansion coefficient for erf
!  interaction
!  a = mu/(2*kF)
!
!  Author : J. Toulouse
!  Date   : 10-03-04
!-------------------------------------------
      implicit none
      include 'constants.include.F'

      double precision a
      double precision eta,fak,berf,berf_dexp

! function
      double precision derf

      eta=19.0d0
      fak=2.540118935556d0*dexp(-eta*a*a)

      if(a .lt. 0.075d0) then
!      expansion for small mu to avoid numerical problems
!      denominator becomes zero for a approximately 0.4845801308
!      (and for one negative and two complex values of a)
       berf = (-7d0+72.d0*a*a)/(27.d0*(-3d0-24.d0*a*a+32.d0*a**4+8d0*dsqrt(pi)*a))

      else if(a .gt. 50.d0) then
       berf = 1.d0/(72.d0*a*a)-1.d0/(17280.d0*a**4)- 23.d0/(358400.d0*a**6)

      else


!      Code generated by Mathematica
       berf_dexp=dexp(2.5d-1/a**2)
       berf = (1.851851851851851851851852d-2*(-1.d0 + 1.44d2*a**4*(-1.d0  &
         + berf_dexp) - 2.d0*a**2*(1.1d1 + 7.d0*berf_dexp                 &
        )))/(a**2*(3.2d1*a**4*(-1.d0 + berf_dexp) - 3.d0*berf_dexp        &
        + 1.417963080724412821838534d1*a*derf(5.d-1/a)*berf_dexp         &
        - 8.d0*a**2*(-2.d0 + 3.d0*berf_dexp)))

      end if

      berf=berf*fak

      return
      end

!-------------------------------------------
      function dberfda(a)
!-------------------------------------------
!  Derivative of second-order exchange gradient
!  expansion coefficient for erf interaction
!  a = mu/(2*kF)
!
!  Author : J. Toulouse
!  Date   : 10-03-04
!-------------------------------------------
      implicit none
      include 'constants.include.F'

      double precision a
      double precision eta,fak,dfakda,berf,dberfda,berf_dexp
      double precision t1,t2,tdexp,t3,t4,t5

      eta=19.0d0
      fak=2.540118935556d0*dexp(-eta*a*a)
      dfakda=-2.0d0*eta*a*fak

      if(a .lt. 0.075d0) then
!      expansion for small mu to avoid numerical problems
!      denominator becomes zero for a approximately 0.4845801308
!      (and for one negative and two complex values of a)
       berf = (-7d0+72.d0*a*a)/(27.d0*(-3d0-24.d0*a*a+32.d0*a**4+8d0*dsqrt(pi)*a))
       dberfda = (8d0*(-96.d0*a + 112.d0*a**3 - 576.d0*a**5      &
        + 7d0*dsqrt(pi) + 72.d0*a**2*dsqrt(pi)))/                &
        (27.d0*(3d0 + 24.d0*a**2 - 32.d0*a**4 - 8d0*a*dsqrt(pi))**2)

      else if(a .gt. 50.d0) then
       berf = 1.d0/(72.d0*a*a)-1.d0/(17280.d0*a**4)- 23.d0/(358400.d0*a**6)
       dberfda = - 1.d0/(36.d0*a**3) +  1.d0/(4320.d0*a**5)+ 69.d0/(179200.d0*a**7)


      else

!      Code generated by Mathematica
       berf_dexp=dexp(2.5d-1/a**2)

       berf = (1.851851851851851851851852d-2*(-1.d0 + 1.44d2*a**4*(-1.d0  + berf_dexp) - 2.d0*a**2*(1.1d1 + 7.d0*berf_dexp )))/(a**2*(3.2d1*a**4*(-1.d0 + berf_dexp) - 3.d0*berf_dexp + 1.417963080724412821838534d1*a*derf(5.d-1/a)*berf_dexp - 8.d0*a**2*(-2.d0 + 3.d0*berf_dexp)))

       tdexp=dexp(2.5d-1/a**2)
       t1 = (1.851851851851851851851852d-2*(5.76d2*a**3*(-1.d0 + tdexp ) + (7.d0*tdexp)/a - 7.2d1*a*tdexp - 4.d0*a*(1.1d1 + 7.d0*tdexp)))/(a**2*(3.2d1*a**4*(-1.d0 + tdexp) - 3.d0*tdexp + 1.417963080724412821838534d1*a*derf(5.d-1/a)*tdexp - 8.d0*a**2*(-2.d0 + 3.d0*tdexp)))
       t2 = -1.851851851851851851851852d-2/a**2
       t3 = -8.d0/a + 1.28d2*a**3*(-1.d0 + tdexp) + (1.5d0*tdexp)/a**3 + (1.2d1*tdexp)/a - 1.6d1*a* tdexp + 1.417963080724412821838534d1*derf(5.d-1/a)*tdexp - (7.08981540362206410919267d0*derf(5.d-1/a)*tdexp)/a**2 - 1.6d1*a*(-2.d0 + 3.d0*tdexp)
       t4 = (-1.d0 + 1.44d2*a**4*(-1.d0 + tdexp) - 2.d0*a**2*(1.1d1 + 7.d0*tdexp))/(3.2d1*a**4*(-1.d0 + tdexp) - 3.d0*tdexp + 1.417963080724412821838534d1*a*derf(5.d-1/a)*tdexp - 8.d0*a**2*(-2.d0 + 3.d0*tdexp))**2
       t5 = (-3.703703703703703703703704d-2*(-1.d0 + 1.44d2*a**4*(-1.d0 + tdexp) - 2.d0*a**2*(1.1d1 + 7.d0*tdexp )))/(a**3*(3.2d1*a**4*(-1.d0 + tdexp) - 3.d0*tdexp+ 1.417963080724412821838534d1*a*derf(5.d-1/a)*tdexp- 8.d0*a**2*(-2.d0 + 3.d0*tdexp)))
       dberfda = t1 + t2*t3*t4 + t5

      end if

      dberfda=dberfda*fak+berf*dfakda

      return
      end
