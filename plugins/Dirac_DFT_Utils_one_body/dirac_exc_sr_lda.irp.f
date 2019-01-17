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
 double precision :: kF_4, tmp_c_2,tmp_c_4,tmp_c_m2
 double precision :: tmp_mu_2,tmp_mu_3,tmp_mu_4,tmp_mu_6,tmp_mu_8 
 double precision :: z0,z1,z2,z3,z4,z5,z6,z8,z12,z13,z18,z20,z24,z27
 double precision :: z45,z60,z96,z540,z648,z960
 double precision :: f13,f120,f23,f43
 double precision :: ckf 
 double precision :: c1,c2,c2bis,c2ter,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12,c13,c14,c15

 z1  = 1.d0
 z2  = 2.d0
 z3 = 3.d0
 z4 = 4.d0
 z5 = 5.d0
 z6 = 6.d0
 z8 = 8.d0
 z12 = 12.d0
 z13 = 13.d0
 z18 = 18.d0
 z20 = 20.d0
 z24 = 24.d0
 z27 = 27.d0
 z45 = 45.d0
 z60 = 60.d0
 z540 = 540.d0
 z648 = 648.d0
 z960 = 960.d0
 f13 = 0.3333333333333333d0
 f120 = 0.05d0
 f23 = 0.6666666666666667d0
 f43 = 1.3333333333333333d0
 ckf = 3.0936677262801355d0
 c1  = 0.008062883608299872d0
 c2  = 0.01905478546791209d0
 c2bis = 0.005375255738866582d0
 c2ter = 3.544907701811032d0
 c3  = 0.3183098861837907d0
 c4  = 0.5641895835477563d0
 c5  = 2.363271801207355d0
 c6  = 0.3183098861837907d0
 c7  = 1.772453850905516d0
 c8 = 0.0008958759564777636d0
 c9 = 0.05305164769729845d0
 c10 = 0.001343813934716645d0
 c11 = 0.1061032953945969d0
 c12 = 0.00004479379782388818d0
 c13 = 7.089815403622064d0
 c14 = 10.6347231054331d0
 c15 = 21.26944621086619d0
 kF = ckf*(rho**f13)
 kF_4 = kF*kF*kF*kF
 c = speed_of_light
 tmp_c = c/kF
 tmp_c_2 = tmp_c*tmp_c
 tmp_c_4 = tmp_c_2*tmp_c_2
 tmp_c_m2 = z1/tmp_c_2
 tmp_mu = mu/kF
 tmp_mu_2 = tmp_mu*tmp_mu
 tmp_mu_3 = tmp_mu_2*tmp_mu
 tmp_mu_4 = tmp_mu_3*tmp_mu
 tmp_mu_6 = tmp_mu_4*tmp_mu_2
 tmp_mu_8 = tmp_mu_6*tmp_mu_2



 ! Non-relativistic equations
 if (tmp_c .gt. 500.d0) then

  ! Linear/quadratic range-separation for very low values of tmp_mu
  if (tmp_mu .lt. 1.d-3) then

  !e_x = (-c1 + c2*tmp_mu)*(kF_4)
   e_x = (-c1 + c2bis*(c2ter - z3*tmp_mu)*tmp_mu)*(kF_4)

  !v_x = (-c3 + c4*tmp_mu)*kF
   v_x = (-c3 + c3*(c7 - tmp_mu)*tmp_mu)*kF  
  
  ! Medium values of tmp_mu
  elseif (tmp_mu .le. 1.d+1) then
 
   e_x =  -c1*kF_4*(z1 + f23*tmp_mu_2*(z3 - z1*tmp_mu_2 + (-z2 + tmp_mu_2)/dexp(z1/tmp_mu_2)) -              &
        c5*tmp_mu*derf(z1/tmp_mu)) 

   v_x = c6*kF*(-z1 + (-z1 + dexp(-z1/tmp_mu_2))*tmp_mu_2 + c7*tmp_mu*derf(z1/tmp_mu)) 
  
  ! For very large values of tmp_mu
  elseif (tmp_mu .lt. 1.d+9) then
  
   e_x = (-c8*kF_4)/tmp_mu_2 

   v_x = (-c9*kF)/tmp_mu_2

  ! Limit for large tmp_mu 
  else

   e_x = 0.d0
   v_x = 0.d0

  endif

 ! Relativistic equations
 else
 
  ! For the Coulomb ee interaction
  if (dirac_interaction == "Coulomb") then
 
   ! Linear/quadratic range-separation for very low values of tmp_mu
   if (tmp_mu .lt. 1.d-1) then 

  ! e_x = kF_4*(c10*(-z2 + tmp_c_2 + z2*(z1 + tmp_c_2)**2*dlog(z1 + z1/tmp_c_2) - z2*dsqrt(z1 + tmp_c_2)*(z2 + z3*tmp_c_2)*dlog(dsqrt(z1 + tmp_c_m2)    & 
  !     + z1/tmp_c) + z3*tmp_c_4*dlog(dsqrt(z1 + tmp_c_m2) + z1/tmp_c)**2) + c2*tmp_mu)                                                                                
    e_x =  kF_4*(c10*(-z2 + tmp_c_2 + z2*(z1 + tmp_c_2)**2*dlog(z1 + z1/tmp_c_2) - z2*dsqrt(z1 + tmp_c_2)*(z2 + z3*tmp_c_2)*dlog(dsqrt(z1 + tmp_c_m2)    & 
        + z1/tmp_c) + z3*tmp_c_4*dlog(dsqrt(z1 + tmp_c_m2) + z1/tmp_c)**2) +c2bis*(c2ter - z3*tmp_mu)*tmp_mu)

  ! v_x =  kF*((c11*(-z1 - z1*tmp_c_2 + (z1 + tmp_c_2)**2*dlog(z1 + z1/tmp_c_2) +  dsqrt(z1 + tmp_c_2)*(z2 + z3*tmp_c_2)*dlog(tmp_c/(z1 +                 & 
  !     dsqrt(z1 + tmp_c_2)))))/(z1 + tmp_c_2) + c4*tmp_mu)
    v_x =  kF*((c11*(-z1 - z1*tmp_c_2 + (z1 + tmp_c_2)**2*dlog(z1 + z1/tmp_c_2) +  dsqrt(z1 + tmp_c_2)*(z2 + z3*tmp_c_2)*dlog(tmp_c/(z1 +                 & 
        dsqrt(z1 + tmp_c_2)))))/(z1 + tmp_c_2) + c3*(c7 - tmp_mu)*tmp_mu)

   ! Medium values of tmp_mu
   elseif (tmp_mu .le. 1d+1) then

    if (dirac_exchange_functional == "dirac_short_range_LDA_P2") then
 
     e_x = (c12*kF_4*(z60*(-z3 + z2*tmp_mu_2*(-z3 + tmp_mu_2 - (z1*(-z2 + tmp_mu_2))*dexp(-z1/tmp_mu_2)) +                       &               
        c13*tmp_mu*derf(z1/tmp_mu)) + (z20*(z1 + z6*tmp_mu_4*(z3 - z2*tmp_mu_2 + (-z1 + z2*tmp_mu_2)*dexp(-z1/tmp_mu_2)) -        &
        c14*tmp_mu_3*derf(z1/tmp_mu))**2 - (z3*(-z2*tmp_mu_2*(-z2 + tmp_mu_2) + dexp(z1/tmp_mu_2)*(-z3 - z6*tmp_mu_2 +                &
        z2*tmp_mu_4 + c13*tmp_mu*derf(z1/tmp_mu)))*(z24*tmp_mu_4*(z4 + z13*tmp_mu_2 - z27*tmp_mu_4) + dexp(z1/tmp_mu_2)*(-z13 -    &
        z540*tmp_mu_4 - z960*tmp_mu_6 + z648*tmp_mu_8 + c15*tmp_mu_3*(z8 + z45*tmp_mu_2)*derf(z1/tmp_mu))))*dexp(-z2/tmp_mu_2))/    &
        (tmp_c_2*(z1 + z6*tmp_mu_4*(z3- z2*tmp_mu_2 + (-z1 + z2*tmp_mu_2)*dexp(-z1/tmp_mu_2)) - c14*tmp_mu_3*derf(z1/tmp_mu)))))/    &
        (z1 + (f120*(z24*tmp_mu_4*(z4 + z13*tmp_mu_2 - z27*tmp_mu_4) + dexp(z1/tmp_mu_2)*(-z13 - z540*tmp_mu_4 - z960*tmp_mu_6 +             & 
        z648*tmp_mu_8 + c15*tmp_mu_3*(z8 + z45*tmp_mu_2)*derf(z1/tmp_mu))))/(tmp_c_2*(z6*tmp_mu_4 - z12*tmp_mu_6 +                    &
        dexp(z1/tmp_mu_2)*(-z1 - z18*tmp_mu_4 + z12*tmp_mu_6 + c14*tmp_mu_3*derf(z1/tmp_mu))))) 
  
     v_x = (0.0353677651315323d0*kF*(864.d0*tmp_mu**10*(96.d0 + 559.d0*tmp_mu_2 - 882.d0*tmp_mu_4 - 1932.d0*tmp_mu_6 + 2734.d0*tmp_mu_8 + 150.d0*tmp_c_4*(z1 -    &
        z2*tmp_mu_2)**2 + 30.d0*tmp_c_2*(z8 + z5*tmp_mu_2 - 86.d0*tmp_mu_4 + 88.d0*tmp_mu_6)) - 432.d0*dexp(z1/tmp_mu_2)*tmp_mu_6*(52.d0 +              & 
        296.d0*tmp_mu_2 + 100.d0*tmp_c_4*(z1 + tmp_mu_2 + 9.d0*tmp_mu_4 - 48.d0*tmp_mu_6 + 36.d0*tmp_mu_8) + tmp_mu_4*(2879.d0 + 4.d0*tmp_mu_2*(2791.d0 -    &
        5520.d0*tmp_mu_2 - 2895.d0*tmp_mu_4 + 4101.d0*tmp_mu_6)) + 5.*tmp_c_2*(29.d0 + 4.d0*tmp_mu_2*(23.d0 + 254.d0*tmp_mu_2 - 247.d0*tmp_mu_4 -               &
        978.d0*tmp_mu_6 + 792.d0*tmp_mu_8))) + 9.d0*dexp(z2/tmp_mu_2)*tmp_mu_2*(169.d0 + 400.d0*tmp_c_4*(-z1 - z18*tmp_mu_4 + z12*tmp_mu_6)*(-z1 -    &
        z12*tmp_mu_2 - z6*tmp_mu_4 + 36.d0*tmp_mu_6) + z8*tmp_mu_2*(247.d0 + tmp_mu_2*(2911.d0 + z6*tmp_mu_2*(2018.d0 + 16613.d0*tmp_mu_2 +             &
        28718.d0*tmp_mu_4 - 55296.d0*tmp_mu_6 - 11568.d0*tmp_mu_8 + 16404.d0*tmp_mu**10))) + 40.d0*tmp_c_2*(z13 + z8*tmp_mu_2*(18.d0 + tmp_mu_2*(143.d0 +    &
        3.d0*tmp_mu_2*(185.d0 + 774.d0*tmp_mu_2 - 923.d0*tmp_mu_4 - 1182.d0*tmp_mu_6 + 792.d0*tmp_mu_8))))) - z1*dexp(z3/tmp_mu_2)*(1261.d0 +                &
        3600.d0*tmp_c_4*(z1 + tmp_mu_2)*(z1 + z18*tmp_mu_4 - z12*tmp_mu_6)**2 + 120.d0*tmp_c_2*(-z1 - z18*tmp_mu_4 + z12*tmp_mu_6)*(-34.d0 +      &
        3.d0*tmp_mu_2*(-z13 - 500.d0*tmp_mu_2 - 1520.d0*tmp_mu_4 - 132.d0*tmp_mu_6 + 528.d0*tmp_mu_8)) + 9.d0*tmp_mu_2*(169.d0 + z8*tmp_mu_2*(1455.d0 +     &
        tmp_mu_2*(4165.d0 + z6*tmp_mu_2*(5629.d0 + 24834.d0*tmp_mu_2 + 29620.d0*tmp_mu_4 - 34980.d0*tmp_mu_6 - 3852.d0*tmp_mu_8 + 5468.d0*tmp_mu**10))))) +    &
        15.95208465814964d0*dexp(z1/tmp_mu_2)*tmp_mu*derf(z1/tmp_mu)* (144.d0*tmp_mu_8*(12.d0*(z4 + z5*tmp_c_2)**2 + z2*(726.d0 + 525.d0*tmp_c_2 -        &
        400.d0*tmp_c_4)*tmp_mu_2 + (113.d0 + 400.d0*tmp_c_2*(-11.d0 + tmp_c_2))*tmp_mu_4 + 4.d0*(-1357.d0 + 490.d0*tmp_c_2)*tmp_mu_6 + 2736.d0*tmp_mu_8) -     &
        48.d0*dexp(z1/tmp_mu_2)*tmp_mu_4*(104.d0 + 738.d0*tmp_mu_2 + 200.d0*tmp_c_4*(z1 + z2*tmp_mu_2 + 15.d0*tmp_mu_4 - 36.d0*tmp_mu_6 +                 &
        z12*tmp_mu_8) + tmp_mu_4*(6881.d0 + z6*tmp_mu_2*(5646.d0 - 1731.d0*tmp_mu_2 - 8164.d0*tmp_mu_4 + 2736.d0*tmp_mu_6)) + 10.d0*tmp_c_2*(29.d0 +        &
        3.d0*tmp_mu_2*(43.d0 + 410.d0*tmp_mu_2 + 214.d0*tmp_mu_4 - 1272.d0*tmp_mu_6 + 392.d0*tmp_mu_8)) + 5.317361552716548d0*tmp_mu_3*(-12.d0*(z4 +           &
        z5*tmp_c_2)**2 + z2*(-893.d0 - 900.d0*tmp_c_2 + 200.d0*tmp_c_4)*tmp_mu_2 + (-2501.d0 + 2680.d0*tmp_c_2)*tmp_mu_4 +                                   &
        4320.d0*tmp_mu_6)*derf(z1/tmp_mu)) + dexp(z2/tmp_mu_2)*(169.d0 + 400.d0*tmp_c_4*(-z1 - z18*tmp_mu_4 + z12*tmp_mu_6)*(-z1 +                     &
        z6*tmp_mu_2*(-z2 - z5*tmp_mu_2 + z2*tmp_mu_4)) + 4.d0*tmp_mu_2*(509.d0 + 3.d0*tmp_mu_2*(2323.d0 + 4.d0*tmp_mu_2*(2575.d0 + 21734.d0*tmp_mu_2 +  &
        56868.d0*tmp_mu_4 - 2517.d0*tmp_mu_6 - 32700.d0*tmp_mu_8 + 8208.d0*tmp_mu**10))) + 40.d0*tmp_c_2*(z13 + tmp_mu_2*(149.d0 + z12*tmp_mu_2*(124.d0 +    &
        3.d0*tmp_mu_2*(157.d0 + 798.d0*tmp_mu_2 + 305.d0*tmp_mu_4 - 832.d0*tmp_mu_6 + 196.d0*tmp_mu_8)))) - 85.07778484346477d0*tmp_mu_3*derf(z1/tmp_mu)*(52.d0+&
        442.d0*tmp_mu_2 + 3.d0*tmp_mu_4*(1354.d0 + 7925.d0*tmp_mu_2 + 6821.d0*tmp_mu_4 - 4320.d0*tmp_mu_6) - 100.d0*tmp_c_4*(-z1 - 3.d0*tmp_mu_2 -            &
        21.d0*tmp_mu_4 + z12*tmp_mu_6) + z5*tmp_c_2*(29.d0 + z2*tmp_mu_2*(83.d0 + 732.d0*tmp_mu_2 + 1344.d0*tmp_mu_4 - 804.d0*tmp_mu_6)) -                &
        3.544907701811032d0*tmp_mu_3*(z6*(z4 + z5*tmp_c_2)**2 + z5*(212.d0 + 255.d0*tmp_c_2)*tmp_mu_2 + 2700.d0*tmp_mu**4)*derf(z1/tmp_mu))))))/          &
        (dexp(z1/tmp_mu_2)*(-z24*tmp_mu**4*(-4.d0 - z13*tmp_mu_2 + z27*tmp_mu**4 + z5*tmp_c_2*(-z1 + z2*tmp_mu_2)) + dexp(z1/tmp_mu_2)*(-z13 -  & 
        20.d0*tmp_c_2 - 180.d0*(z3 + z2*tmp_c_2)*tmp_mu**4 + 240.d0*(-4.d0 + tmp_c_2)*tmp_mu_6 + z648*tmp_mu_8 + c15*tmp_mu_3*(z8 +      &
        10.d0*tmp_c_2 + z45*tmp_mu**2)*derf(z1/tmp_mu)))**2)




































!    e_x = (0.00004479379782388818d0*kF**4*(60.d0*(-3.d0 + 2.d0*tmp_mu**2*(-3.d0 + tmp_mu**2 - (1.d0*(-2.d0 + tmp_mu**2))*dexp(-1.d0/tmp_mu**2)) +                       & 
!       7.089815403622064d0*tmp_mu*derf(1.d0/tmp_mu)) + (20.d0*(1.d0 + 6.d0*tmp_mu**4*(3.d0 - 2.d0*tmp_mu**2 + (-1.d0 + 2.d0*tmp_mu**2)*dexp(-1.d0/tmp_mu**2)) -        &
!       10.6347231054331d0*tmp_mu**3*derf(1.d0/tmp_mu))**2 - (3.d0*(-2.d0*tmp_mu**2*(-2.d0 + tmp_mu**2) + dexp(1.d0/tmp_mu**2)*(-3.d0 - 6.d0*tmp_mu**2 +                &
!       2.d0*tmp_mu**4 + 7.089815403622064d0*tmp_mu*derf(1.d0/tmp_mu)))*(24.d0*tmp_mu**4*(4.d0 + 13.d0*tmp_mu**2 - 27.d0*tmp_mu**4) + dexp(1.d0/tmp_mu**2)*(-13.d0 -    &
!       540.d0*tmp_mu**4 - 960.d0*tmp_mu**6 + 648.d0*tmp_mu**8 + 21.26944621086619d0*tmp_mu**3*(8.d0 + 45.d0*tmp_mu**2)*derf(1.d0/tmp_mu))))*dexp(-2.d0/tmp_mu**2))/    &
!       (tmp_c**2*(1.d0 + 6.d0*tmp_mu**4*(3.d0- 2.d0*tmp_mu**2 + (-1.d0 + 2.d0*tmp_mu**2)*dexp(-1.d0/tmp_mu**2)) - 10.6347231054331*tmp_mu**3*derf(1.d0/tmp_mu)))))/    &
!       (1.d0 + (0.05d0*(24.d0*tmp_mu**4*(4.d0 + 13.d0*tmp_mu**2 - 27.d0*tmp_mu**4) + dexp(1.d0/tmp_mu**2)*(-13.d0 - 540.d0*tmp_mu**4 - 960.d0*tmp_mu**6 +             & 
!       648.d0*tmp_mu**8 + 21.26944621086619d0*tmp_mu**3*(8.d0 + 45.d0*tmp_mu**2)*derf(1.d0/tmp_mu))))/(tmp_c**2*(6.d0*tmp_mu**4 - 12.d0*tmp_mu**6 +                    &
!       dexp(1.d0/tmp_mu**2)*(-1.d0 - 18.d0*tmp_mu**4 + 12.d0*tmp_mu**6 + 10.6347231054331d0*tmp_mu**3*derf(1.d0/tmp_mu))))) 
  
!    v_x = (0.0353677651315323d0*kF*(864.d0*tmp_mu**10*(96.d0 + 559.d0*tmp_mu**2 - 882.d0*tmp_mu**4 - 1932.d0*tmp_mu**6 + 2734.d0*tmp_mu**8 + 150.d0*tmp_c**4*(1.d0 -    &
!       2.d0*tmp_mu**2)**2 + 30.d0*tmp_c**2*(8.d0 + 5.d0*tmp_mu**2 - 86.d0*tmp_mu**4 + 88.d0*tmp_mu**6)) - 432.d0*dexp(1.d0/tmp_mu**2)*tmp_mu**6*(52.d0 +              & 
!       296.d0*tmp_mu**2 + 100.d0*tmp_c**4*(1.d0 + tmp_mu**2 + 9.d0*tmp_mu**4 - 48.d0*tmp_mu**6 + 36.d0*tmp_mu**8) + tmp_mu**4*(2879.d0 + 4.d0*tmp_mu**2*(2791.d0 -    &
!       5520.d0*tmp_mu**2 - 2895.d0*tmp_mu**4 + 4101.d0*tmp_mu**6)) + 5.*tmp_c**2*(29.d0 + 4.d0*tmp_mu**2*(23.d0 + 254.d0*tmp_mu**2 - 247.d0*tmp_mu**4 -               &
!       978.d0*tmp_mu**6 + 792.d0*tmp_mu**8))) + 9.d0*dexp(2.d0/tmp_mu**2)*tmp_mu**2*(169.d0 + 400.d0*tmp_c**4*(-1.d0 - 18.d0*tmp_mu**4 + 12.d0*tmp_mu**6)*(-1.d0 -    &
!       12.d0*tmp_mu**2 - 6.d0*tmp_mu**4 + 36.d0*tmp_mu**6) + 8.d0*tmp_mu**2*(247.d0 + tmp_mu**2*(2911.d0 + 6.d0*tmp_mu**2*(2018.d0 + 16613.d0*tmp_mu**2 +             &
!       28718.d0*tmp_mu**4 - 55296.d0*tmp_mu**6 - 11568.d0*tmp_mu**8 + 16404.d0*tmp_mu**10))) + 40.d0*tmp_c**2*(13.d0 + 8.d0*tmp_mu**2*(18.d0 + tmp_mu**2*(143.d0 +    &
!       3.d0*tmp_mu**2*(185.d0 + 774.d0*tmp_mu**2 - 923.d0*tmp_mu**4 - 1182.d0*tmp_mu**6 + 792.d0*tmp_mu**8))))) - 1.d0*dexp(3.d0/tmp_mu**2)*(1261.d0 +                &
!       3600.d0*tmp_c**4*(1.d0 + tmp_mu**2)*(1.d0 + 18.d0*tmp_mu**4 - 12.d0*tmp_mu**6)**2 + 120.d0*tmp_c**2*(-1.d0 - 18.d0*tmp_mu**4 + 12.d0*tmp_mu**6)*(-34.d0 +      &
!       3.d0*tmp_mu**2*(-13.d0 - 500.d0*tmp_mu**2 - 1520.d0*tmp_mu**4 - 132.d0*tmp_mu**6 + 528.d0*tmp_mu**8)) + 9.d0*tmp_mu**2*(169.d0 + 8.d0*tmp_mu**2*(1455.d0 +     &
!       tmp_mu**2*(4165.d0 + 6.d0*tmp_mu**2*(5629.d0 + 24834.d0*tmp_mu**2 + 29620.d0*tmp_mu**4 - 34980.d0*tmp_mu**6 - 3852.d0*tmp_mu**8 + 5468.d0*tmp_mu**10))))) +    &
!       15.95208465814964d0*dexp(1.d0/tmp_mu**2)*tmp_mu*derf(1.d0/tmp_mu)* (144.d0*tmp_mu**8*(12.d0*(4.d0 + 5.d0*tmp_c**2)**2 + 2.d0*(726.d0 + 525.d0*tmp_c**2 -        &
!       400.d0*tmp_c**4)*tmp_mu**2 + (113.d0 + 400.d0*tmp_c**2*(-11.d0 + tmp_c**2))*tmp_mu**4 + 4.d0*(-1357.d0 + 490.d0*tmp_c**2)*tmp_mu**6 + 2736.d0*tmp_mu**8) -     &
!       48.d0*dexp(1.d0/tmp_mu**2)*tmp_mu**4*(104.d0 + 738.d0*tmp_mu**2 + 200.d0*tmp_c**4*(1.d0 + 2.d0*tmp_mu**2 + 15.d0*tmp_mu**4 - 36.d0*tmp_mu**6 +                 &
!       12.d0*tmp_mu**8) + tmp_mu**4*(6881.d0 + 6.d0*tmp_mu**2*(5646.d0 - 1731.d0*tmp_mu**2 - 8164.d0*tmp_mu**4 + 2736.d0*tmp_mu**6)) + 10.d0*tmp_c**2*(29.d0 +        &
!       3.d0*tmp_mu**2*(43.d0 + 410.d0*tmp_mu**2 + 214.d0*tmp_mu**4 - 1272.d0*tmp_mu**6 + 392.d0*tmp_mu**8)) + 5.317361552716548d0*tmp_mu**3*(-12.d0*(4.d0 +           &
!       5.d0*tmp_c**2)**2 + 2.d0*(-893.d0 - 900.d0*tmp_c**2 + 200.d0*tmp_c**4)*tmp_mu**2 + (-2501.d0 + 2680.d0*tmp_c**2)*tmp_mu**4 +                                   &
!       4320.d0*tmp_mu**6)*derf(1.d0/tmp_mu)) + dexp(2.d0/tmp_mu**2)*(169.d0 + 400.d0*tmp_c**4*(-1.d0 - 18.d0*tmp_mu**4 + 12.d0*tmp_mu**6)*(-1.d0 +                     &
!       6.d0*tmp_mu**2*(-2.d0 - 5.d0*tmp_mu**2 + 2.d0*tmp_mu**4)) + 4.d0*tmp_mu**2*(509.d0 + 3.d0*tmp_mu**2*(2323.d0 + 4.d0*tmp_mu**2*(2575.d0 + 21734.d0*tmp_mu**2 +  &
!       56868.d0*tmp_mu**4 - 2517.d0*tmp_mu**6 - 32700.d0*tmp_mu**8 + 8208.d0*tmp_mu**10))) + 40.d0*tmp_c**2*(13.d0 + tmp_mu**2*(149.d0 + 12.d0*tmp_mu**2*(124.d0 +    &
!       3.d0*tmp_mu**2*(157.d0 + 798.d0*tmp_mu**2 + 305.d0*tmp_mu**4 - 832.d0*tmp_mu**6 + 196.d0*tmp_mu**8)))) - 85.07778484346477d0*tmp_mu**3*derf(1.d0/tmp_mu)*(52.d0+&
!       442.d0*tmp_mu**2 + 3.d0*tmp_mu**4*(1354.d0 + 7925.d0*tmp_mu**2 + 6821.d0*tmp_mu**4 - 4320.d0*tmp_mu**6) - 100.d0*tmp_c**4*(-1.d0 - 3.d0*tmp_mu**2 -            &
!       21.d0*tmp_mu**4 + 12.d0*tmp_mu**6) + 5.d0*tmp_c**2*(29.d0 + 2.d0*tmp_mu**2*(83.d0 + 732.d0*tmp_mu**2 + 1344.d0*tmp_mu**4 - 804.d0*tmp_mu**6)) -                &
!       3.544907701811032d0*tmp_mu**3*(6.d0*(4.d0 + 5.d0*tmp_c**2)**2 + 5.d0*(212.d0 + 255.d0*tmp_c**2)*tmp_mu**2 + 2700.d0*tmp_mu**4)*derf(1.d0/tmp_mu))))))/          &
!       (dexp(1.d0/tmp_mu**2)*(-24.d0*tmp_mu**4*(-4.d0 - 13.d0*tmp_mu**2 + 27.d0*tmp_mu**4 + 5.d0*tmp_c**2*(-1.d0 + 2.d0*tmp_mu**2)) + dexp(1.d0/tmp_mu**2)*(-13.d0 -  & 
!       20.d0*tmp_c**2 - 180.d0*(3.d0 + 2.d0*tmp_c**2)*tmp_mu**4 + 240.d0*(-4.d0 + tmp_c**2)*tmp_mu**6 + 648.d0*tmp_mu**8 + 21.26944621086619d0*tmp_mu**3*(8.d0 +      &
!       10.d0*tmp_c**2 + 45.d0*tmp_mu**2)*derf(1.d0/tmp_mu)))**2)
    else
     print*, 'Exchange functional required does not exist ...'
     print*,'dirac_exchange_functional',dirac_exchange_functional
     stop
    endif

   ! For very large values of tmp_mu
   elseif (tmp_mu .lt. 1.d+9) then

    e_x =  (-0.0001119844945597204d0*kF**4*(4.d0 + 9.d0*(tmp_c**2 + tmp_c**4) + 9.d0*tmp_c**4*dlog(dsqrt(1.d0 + tmp_c**(-2)) + 1.d0/tmp_c)*(-2.d0*dsqrt(1.d0 + tmp_c**2) +  & 
        tmp_c**2*dlog(dsqrt(1.d0 + tmp_c**(-2)) + 1.d0/tmp_c))))/tmp_mu**2

    v_x =  (-0.01326291192432461d0*kF*(2.d0 + 5.d0*tmp_c**2 + 3.d0*tmp_c**4 + 3.d0*tmp_c**4*dsqrt(1.d0 + tmp_c**2)*dlog(tmp_c/(1.d0 + dsqrt(1.d0 + tmp_c**2)))))/((1.d0 +    & 
        tmp_c**2)*tmp_mu**2) 

   ! Limit for large tmp_mu 
   else
  
    e_x = 0.d0
    v_x = 0.d0 
    
   endif 
  endif 
 endif
 end
