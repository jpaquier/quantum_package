 subroutine dirac_ec_lda_sr(mu,rho,e_c,v_c)
      implicit none
 include 'constants.include.F'
      double precision, intent(out) ::  e_c
      double precision, intent(out) ::  v_c
      double precision, intent(in)  ::  mu,rho
     e_c = 0.d0 
     v_c = 0.d0
 end

 subroutine dirac_ex_lda_sr(mu,rho,e_x,v_x)
 include 'constants.include.F'
 implicit none 
 double precision, intent(out) ::  e_x
 double precision, intent(out) ::  v_x
 double precision, intent(in)  ::  rho,mu
 double precision :: kF, c, tmp_c, tmp_mu
 
!Limit for small a (expansion not so important as for large a)
!if (a.lt.1.d-9) then
!  e_x = -z3/z8*rho_a_2*(z24*rho_a_2/pi)**f13
!  v_x = - ((z3/pi)*rho_a_2)**f13

!Intermediate values of a
!elseif (a.le.100d0) then
!  e_x = - (rho_a_2*(z24*rho_a_2/pi)**f13) * (z3/z8-a*(sqpi*derf(f12/a)+(z2*a-z4*a3)*dexp(-f14/a2)-z3*a+z4*a3))
!  v_x =  -(z3*rho_a_2/pi)**f13 + z2*a*mu/pi*(dexp(-f14/a2)-z1)+mu/sqpi * derf(f12/a)


!Expansion for large a
!elseif (a.lt.1.d+9) then
!  e_x = -(rho_a_2*(z24*rho_a_2/pi)**f13) * z1/(z96*a2)
!  v_x = -pi*rho_a_2/(z2*mu*mu)

!Limit for large a
!else
!  e_x = 0.d0
!  v_x = 0.d0
!end if

!double precision, intent(out) ::  e_x
!double precision, intent(out) ::  v_x
!double precision, intent(in)  ::  rho,mu
!double precision :: kF, c, tmp_c, tmp_mu
 kF = (3* pi**2 *rho)**(c_1_3)
 c = speed_of_light
 tmp_c = c/kF
 tmp_mu = mu/kF 
 
 ! Non-relativistic equations
 if (tmp_c .gt. 500.d0) then

  ! Linear range-separation for very low values of tmp_mu
  if (tmp_mu .lt. 1.d-9) then

   e_x = -0.008062883608299872d0*kF**4 + 0.01905478546791209d0*kF**4*tmp_mu  

   v_x = -0.3183098861837907d0*kF + 0.5641895835477563d0*kF*tmp_mu

  ! Medium values of tmp_mu
  elseif (tmp_mu .le. 100) then
 
  e_x =  -0.008062883608299872d0*kF**4*(1.d0 + 0.6666666666666667d0*tmp_mu**2*(3.d0 - 1.d0*tmp_mu**2 + (-2.d0 + tmp_mu**2)/dexp(1.d0/tmp_mu**2)) -              &
        2.363271801207355d0*tmp_mu*derf(1.d0/tmp_mu)) 

  v_x = 0.3183098861837907d0*kF*(-1.d0 + (-1.d0 + dexp(-1.d0/tmp_mu**2))*tmp_mu**2 + 1.772453850905516d0*tmp_mu*derf(1.d0/tmp_mu)) 
  
  ! For very large values of tmp_mu
  elseif (tmp_mu .lt. 1.d+9) then

   e_x = (-0.0008958759564777636d0*kF**4)/tmp_mu**2 

   v_x = (-0.05305164769729845d0*kF)/tmp_mu**2

  ! Limit for large tmp_mu 
  else

   e_x = 0.d0
   v_x = 0.d0

  endif

 ! Relativistic equations
 else
 
  ! For the Coulomb ee interaction
  if (dirac_interaction == "Coulomb") then
 
   ! Linear range-separation for very low values of tmp_mu
   if (tmp_mu .lt. 1.d-9) then 

    e_x =  0.001343813934716645d0*kF**4*(-2.d0 + tmp_c**2 + 2.d0*(1.d0 + tmp_c**2)**2*Log(1.d0 + tmp_c**(-2)) - 2.d0*Sqrt(1.d0 + tmp_c**2)*(2.d0 +               &
        3.d0*tmp_c**2)*Log(Sqrt(1.d0 + tmp_c**(-2)) + 1.d0/tmp_c) + 3.d0*tmp_c**4*Log(Sqrt(1.d0 + tmp_c**(-2)) + 1.d0/tmp_c)**2) +                               &
        0.01905478546791209*kF**4*tmp_mu                                                                                
 
    v_x =   (0.1061032953945969d0*kF*(-1.d0 - 1.d0*tmp_c**2 + (1.d0 + tmp_c**2)**2*Log(1.d0 + tmp_c**(-2)) +  Sqrt(1.d0 + tmp_c**2)*(2.d0 +                      &
        3.d0*tmp_c**2)*Log(tmp_c/(1.d0 + Sqrt(1.d0 + tmp_c**2)))))/(1.d0 + tmp_c**2) +                                                                          &
        0.5641895835477563*kF*tmp_mu
 
   ! Medium values of tmp_mu
   elseif (tmp_mu .le. 100) then

    if (dirac_exchange_functional == "dirac_short_range_LDA_P2") then
     e_x = (0.00004479379782388818d0*kF**4*(60.d0*(-3.d0 + 2.d0*tmp_mu**2*(-3.d0 + tmp_mu**2 - (1.d0*(-2.d0 + tmp_mu**2))*dexp(-1.d0/tmp_mu**2)) +                       &               
        7.089815403622064d0*tmp_mu*derf(1.d0/tmp_mu)) + (20.d0*(1.d0 + 6.d0*tmp_mu**4*(3.d0 - 2.d0*tmp_mu**2 + (-1.d0 + 2.d0*tmp_mu**2)*dexp(-1.d0/tmp_mu**2)) -        &
        10.6347231054331d0*tmp_mu**3*derf(1.d0/tmp_mu))**2 - (3.d0*(-2.d0*tmp_mu**2*(-2.d0 + tmp_mu**2) + dexp(1.d0/tmp_mu**2)*(-3.d0 - 6.d0*tmp_mu**2 +                &
        2.d0*tmp_mu**4 + 7.089815403622064d0*tmp_mu*derf(1.d0/tmp_mu)))*(24.d0*tmp_mu**4*(4.d0 + 13.d0*tmp_mu**2 - 27.d0*tmp_mu**4) + dexp(1.d0/tmp_mu**2)*(-13.d0 -    &
        540.d0*tmp_mu**4 - 960.d0*tmp_mu**6 + 648.d0*tmp_mu**8 + 21.26944621086619d0*tmp_mu**3*(8.d0 + 45.d0*tmp_mu**2)*derf(1.d0/tmp_mu))))*dexp(-2.d0/tmp_mu**2))/    &
        (tmp_c**2*(1.d0 + 6.d0*tmp_mu**4*(3.d0- 2.d0*tmp_mu**2 + (-1.d0 + 2.d0*tmp_mu**2)*dexp(-1.d0/tmp_mu**2)) - 10.6347231054331*tmp_mu**3*derf(1.d0/tmp_mu)))))/    &
        (1.d0 + (0.05d0*(24.d0*tmp_mu**4*(4.d0 + 13.d0*tmp_mu**2 - 27.d0*tmp_mu**4) + dexp(1.d0/tmp_mu**2)*(-13.d0 - 540.d0*tmp_mu**4 - 960.d0*tmp_mu**6 +             & 
        648.d0*tmp_mu**8 + 21.26944621086619d0*tmp_mu**3*(8.d0 + 45.d0*tmp_mu**2)*derf(1.d0/tmp_mu))))/(tmp_c**2*(6.d0*tmp_mu**4 - 12.d0*tmp_mu**6 +                    &
        dexp(1.d0/tmp_mu**2)*(-1.d0 - 18.d0*tmp_mu**4 + 12.d0*tmp_mu**6 + 10.6347231054331d0*tmp_mu**3*derf(1.d0/tmp_mu))))) 
  
     v_x = (0.0353677651315323d0*kF*(864.d0*tmp_mu**10*(96.d0 + 559.d0*tmp_mu**2 - 882.d0*tmp_mu**4 - 1932.d0*tmp_mu**6 + 2734.d0*tmp_mu**8 + 150.d0*tmp_c**4*(1.d0 -    &
        2.d0*tmp_mu**2)**2 + 30.d0*tmp_c**2*(8.d0 + 5.d0*tmp_mu**2 - 86.d0*tmp_mu**4 + 88.d0*tmp_mu**6)) - 432.d0*dexp(1.d0/tmp_mu**2)*tmp_mu**6*(52.d0 +              & 
        296.d0*tmp_mu**2 + 100.d0*tmp_c**4*(1.d0 + tmp_mu**2 + 9.d0*tmp_mu**4 - 48.d0*tmp_mu**6 + 36.d0*tmp_mu**8) + tmp_mu**4*(2879.d0 + 4.d0*tmp_mu**2*(2791.d0 -    &
        5520.d0*tmp_mu**2 - 2895.d0*tmp_mu**4 + 4101.d0*tmp_mu**6)) + 5.*tmp_c**2*(29.d0 + 4.d0*tmp_mu**2*(23.d0 + 254.d0*tmp_mu**2 - 247.d0*tmp_mu**4 -               &
        978.d0*tmp_mu**6 + 792.d0*tmp_mu**8))) + 9.d0*dexp(2.d0/tmp_mu**2)*tmp_mu**2*(169.d0 + 400.d0*tmp_c**4*(-1.d0 - 18.d0*tmp_mu**4 + 12.d0*tmp_mu**6)*(-1.d0 -    &
        12.d0*tmp_mu**2 - 6.d0*tmp_mu**4 + 36.d0*tmp_mu**6) + 8.d0*tmp_mu**2*(247.d0 + tmp_mu**2*(2911.d0 + 6.d0*tmp_mu**2*(2018.d0 + 16613.d0*tmp_mu**2 +             &
        28718.d0*tmp_mu**4 - 55296.d0*tmp_mu**6 - 11568.d0*tmp_mu**8 + 16404.d0*tmp_mu**10))) + 40.d0*tmp_c**2*(13.d0 + 8.d0*tmp_mu**2*(18.d0 + tmp_mu**2*(143.d0 +    &
        3.d0*tmp_mu**2*(185.d0 + 774.d0*tmp_mu**2 - 923.d0*tmp_mu**4 - 1182.d0*tmp_mu**6 + 792.d0*tmp_mu**8))))) - 1.d0*dexp(3.d0/tmp_mu**2)*(1261.d0 +                &
        3600.d0*tmp_c**4*(1.d0 + tmp_mu**2)*(1.d0 + 18.d0*tmp_mu**4 - 12.d0*tmp_mu**6)**2 + 120.d0*tmp_c**2*(-1.d0 - 18.d0*tmp_mu**4 + 12.d0*tmp_mu**6)*(-34.d0 +      &
        3.d0*tmp_mu**2*(-13.d0 - 500.d0*tmp_mu**2 - 1520.d0*tmp_mu**4 - 132.d0*tmp_mu**6 + 528.d0*tmp_mu**8)) + 9.d0*tmp_mu**2*(169.d0 + 8.d0*tmp_mu**2*(1455.d0 +     &
        tmp_mu**2*(4165.d0 + 6.d0*tmp_mu**2*(5629.d0 + 24834.d0*tmp_mu**2 + 29620.d0*tmp_mu**4 - 34980.d0*tmp_mu**6 - 3852.d0*tmp_mu**8 + 5468.d0*tmp_mu**10))))) +    &
        15.95208465814964d0*dexp(1.d0/tmp_mu**2)*tmp_mu*derf(1.d0/tmp_mu)* (144.d0*tmp_mu**8*(12.d0*(4.d0 + 5.d0*tmp_c**2)**2 + 2.d0*(726.d0 + 525.d0*tmp_c**2 -        &
        400.d0*tmp_c**4)*tmp_mu**2 + (113.d0 + 400.d0*tmp_c**2*(-11.d0 + tmp_c**2))*tmp_mu**4 + 4.d0*(-1357.d0 + 490.d0*tmp_c**2)*tmp_mu**6 + 2736.d0*tmp_mu**8) -     &
        48.d0*dexp(1.d0/tmp_mu**2)*tmp_mu**4*(104.d0 + 738.d0*tmp_mu**2 + 200.d0*tmp_c**4*(1.d0 + 2.d0*tmp_mu**2 + 15.d0*tmp_mu**4 - 36.d0*tmp_mu**6 +                 &
        12.d0*tmp_mu**8) + tmp_mu**4*(6881.d0 + 6.d0*tmp_mu**2*(5646.d0 - 1731.d0*tmp_mu**2 - 8164.d0*tmp_mu**4 + 2736.d0*tmp_mu**6)) + 10.d0*tmp_c**2*(29.d0 +        &
        3.d0*tmp_mu**2*(43.d0 + 410.d0*tmp_mu**2 + 214.d0*tmp_mu**4 - 1272.d0*tmp_mu**6 + 392.d0*tmp_mu**8)) + 5.317361552716548d0*tmp_mu**3*(-12.d0*(4.d0 +           &
        5.d0*tmp_c**2)**2 + 2.d0*(-893.d0 - 900.d0*tmp_c**2 + 200.d0*tmp_c**4)*tmp_mu**2 + (-2501.d0 + 2680.d0*tmp_c**2)*tmp_mu**4 +                                   &
        4320.d0*tmp_mu**6)*derf(1.d0/tmp_mu)) + dexp(2.d0/tmp_mu**2)*(169.d0 + 400.d0*tmp_c**4*(-1.d0 - 18.d0*tmp_mu**4 + 12.d0*tmp_mu**6)*(-1.d0 +                     &
        6.d0*tmp_mu**2*(-2.d0 - 5.d0*tmp_mu**2 + 2.d0*tmp_mu**4)) + 4.d0*tmp_mu**2*(509.d0 + 3.d0*tmp_mu**2*(2323.d0 + 4.d0*tmp_mu**2*(2575.d0 + 21734.d0*tmp_mu**2 +  &
        56868.d0*tmp_mu**4 - 2517.d0*tmp_mu**6 - 32700.d0*tmp_mu**8 + 8208.d0*tmp_mu**10))) + 40.d0*tmp_c**2*(13.d0 + tmp_mu**2*(149.d0 + 12.d0*tmp_mu**2*(124.d0 +    &
        3.d0*tmp_mu**2*(157.d0 + 798.d0*tmp_mu**2 + 305.d0*tmp_mu**4 - 832.d0*tmp_mu**6 + 196.d0*tmp_mu**8)))) - 85.07778484346477d0*tmp_mu**3*derf(1.d0/tmp_mu)*(52.d0+&
        442.d0*tmp_mu**2 + 3.d0*tmp_mu**4*(1354.d0 + 7925.d0*tmp_mu**2 + 6821.d0*tmp_mu**4 - 4320.d0*tmp_mu**6) - 100.d0*tmp_c**4*(-1.d0 - 3.d0*tmp_mu**2 -            &
        21.d0*tmp_mu**4 + 12.d0*tmp_mu**6) + 5.d0*tmp_c**2*(29.d0 + 2.d0*tmp_mu**2*(83.d0 + 732.d0*tmp_mu**2 + 1344.d0*tmp_mu**4 - 804.d0*tmp_mu**6)) -                &
        3.544907701811032d0*tmp_mu**3*(6.d0*(4.d0 + 5.d0*tmp_c**2)**2 + 5.d0*(212.d0 + 255.d0*tmp_c**2)*tmp_mu**2 + 2700.d0*tmp_mu**4)*derf(1.d0/tmp_mu))))))/          &
        (dexp(1.d0/tmp_mu**2)*(-24.d0*tmp_mu**4*(-4.d0 - 13.d0*tmp_mu**2 + 27.d0*tmp_mu**4 + 5.d0*tmp_c**2*(-1.d0 + 2.d0*tmp_mu**2)) + dexp(1.d0/tmp_mu**2)*(-13.d0 -  & 
        20.d0*tmp_c**2 - 180.d0*(3.d0 + 2.d0*tmp_c**2)*tmp_mu**4 + 240.d0*(-4.d0 + tmp_c**2)*tmp_mu**6 + 648.d0*tmp_mu**8 + 21.26944621086619d0*tmp_mu**3*(8.d0 +      &
        10.d0*tmp_c**2 + 45.d0*tmp_mu**2)*derf(1.d0/tmp_mu)))**2)
    else
     print*, 'Exchange functional required does not exist ...'
     print*,'dirac_exchange_functional',dirac_exchange_functional
     stop
    endif

   ! For very large values of tmp_mu
   elseif (tmp_mu .lt. 1.d+9) then

    e_x =  (-0.0001119844945597204d0*kF**4*(4.d0 + 9.d0*(tmp_c**2 + tmp_c**4) + 9.d0*tmp_c**4*Log(Sqrt(1.d0 + tmp_c**(-2)) + 1.d0/tmp_c)*(-2.d0*Sqrt(1.d0 + tmp_c**2) +  & 
        tmp_c**2*Log(Sqrt(1.d0 + tmp_c**(-2)) + 1.d0/tmp_c))))/tmp_mu**2

    v_x =  (-0.01326291192432461d0*kF*(2.d0 + 5.d0*tmp_c**2 + 3.d0*tmp_c**4 + 3.d0*tmp_c**4*Sqrt(1.d0 + tmp_c**2)*Log(tmp_c/(1.d0 + Sqrt(1.d0 + tmp_c**2)))))/((1.d0 +    & 
        tmp_c**2)*tmp_mu**2) 

   ! Limit for large tmp_mu 
   else
  
    e_x = 0.d0
    v_x = 0.d0 
    
   endif 
  endif 
 endif
 end
