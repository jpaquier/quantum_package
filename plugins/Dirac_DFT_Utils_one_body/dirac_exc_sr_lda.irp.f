 subroutine dirac_ex_lda(rho,ex,vx)
 include 'constants.include.F'
 implicit none
 double precision, intent(in) :: rho
 double precision, intent(out) :: ex,vx
 double precision :: kF,tmp_c
 kF = (3* pi**2 *rho)**(c_1_3)
 tmp_c = speed_of_light/kF
 ex =     -0.008062883608299872*kF**4*(0.8333333333333333 + 0.3333333333333333*tmp_c**2 +  & 
          0.6666666666666667*Sqrt(1.d0 + tmp_c**2)*ArcSinh(1.d0/tmp_c) -                   &  
          0.5*(Sqrt(1.d0 + tmp_c**2) - tmp_c**2*ArcSinh(1.d0/tmp_c))**2 -                  &  
          0.3333333333333333*(1.d0 + tmp_c**2)**2*Log(1.d0 + tmp_c**(-2)))                   

 vx =    -2.583856390024985*kF*(2.d0 + tmp_c*(2.d0*Sqrt(1.d0 + tmp_c**(-2))*Sqrt(1.d0 + tmp_c**2) +     &   
          tmp_c*(5.d0 + 3.d0*tmp_c*(tmp_c - Sqrt(1.d0 + tmp_c**(-2))*Sqrt(1.d0 + tmp_c**2)))) +         &
          tmp_c*(8.d0*Sqrt(1.d0 + tmp_c**(-2)) + 3.d0*tmp_c**2*(4.d0*Sqrt(1.d0 + tmp_c**(-2)) +         &
          Sqrt(1.d0 + tmp_c**(-2))*tmp_c**2 - 1.d0*tmp_c*Sqrt(1.d0 + tmp_c**2)))*ArcCsch(tmp_c)         &
          - 4.d0*Sqrt(1.d0 + tmp_c**(-2))*tmp_c*(1.d0 + tmp_c**2)**1.5*Log(1.d0 + tmp_c**(-2)))         &
                /(Sqrt(1.d0 + tmp_c**(-2))*tmp_c*Sqrt(1.d0 + tmp_c**2))      
 end


 subroutine dirac_ec_lda(rho,ec,vc)
      implicit none
 include 'constants.include.F'
      double precision, intent(out) ::  ec
      double precision, intent(out) ::  vc
      double precision, intent(in)  ::  rho
  ec = 0.d0
  vc = 0.d0
 end

 subroutine dirac_ec_lda_sr(mu,rho,ec,vc)
      implicit none
 include 'constants.include.F'
      double precision, intent(out) ::  ec
      double precision, intent(out) ::  vd0c
      double precision, intent(in)  ::  mu,rho
     ec = 0.d0 
     vc = 0.d0
 end

 subroutine dirac_ex_lda_sr(mu,rho,ex,vx)
 include 'constants.include.F'
 implicit none 
 double precision, intent(out) ::  ex
 double precision, intent(out) ::  vx
 double precision, intent(in)  ::  rho,mu
 double precision :: kF, tmp_c, tmp_mu
 kF = (3* pi**2 *rho)**(c_1_3)
 tmp_c = speed_of_light/kF
 tmp_mu = mu_erf/kF
 if (dirac_exchange_functional == dirac_short_range_LDA_P2) then
  ex =  ((0.07957747154594767*kF*(-2.d0*tmp_mu**2*(-2.d0 + tmp_mu**2) + 2.718281828459045**(1.d0/tmp_mu**2)*(-3.d0 - 6.d0*tmp_mu**2 + 2.d0*tmp_mu**4) + &
        7.089815403622064*2.718281828459045**(1.d0/tmp_mu**2)*tmp_mu*Erf(1.d0/tmp_mu)))/E**(1.d0/tmp_mu**2) + (0.02652582384864922*kF*(1.d0 +           &   
        (18.d0 - 6.d0/E**(1.d0/tmp_mu**2))*tmp_mu**4 + 12.d0*(-1.d0 + E**(-1.d0/tmp_mu**2))*tmp_mu**6 - 10.6347231054331*tmp_mu**3*Erf(1.d0/tmp_mu)) -  &
        (0.003978873577297383*kF*(-2.d0*tmp_mu**2*(-2.d0 + tmp_mu**2) + 2.718281828459045**(1.d0/tmp_mu**2)*(-3.d0 - 6.d0*tmp_mu**2 + 2.d0*tmp_mu**4) + &
        7.089815403622064*2.718281828459045**(1.d0/tmp_mu**2)*tmp_mu*Erf(1.d0/tmp_mu))*(96.d0*tmp_mu**4 + 312.d0*tmp_mu**6 - 648.d0*tmp_mu**8 +         &
        2.718281828459045**(1.d0/tmp_mu**2)*(-13.d0 - 540.d0*tmp_mu**4 - 960.d0*tmp_mu**6 + 648.d0*tmp_mu**8)  +                                        & 
        957.12507948898*2.718281828459045**(1.d0/tmp_mu**2)*tmp_mu**3*(0.177777777777778 + 1.d0*tmp_mu**2)*Erf(1.d0/tmp_mu)))                           &
        /(E**(2.d0/tmp_mu**2)*(1.d0 + (18.d0 - 6.d0/E**(1.d0/tmp_mu**2))*tmp_mu**4 + 12.d0*(-1.d0 + E**(-1.d0/tmp_mu**2))*tmp_mu**6 -                   &
        10.6347231054331*tmp_mu**3*Erf(1.d0/tmp_mu))))/tmp_c**2)                                                                                        &
        /(1.d0 - (0.05d0*(96.d0*tmp_mu**4 + 312.d0*tmp_mu**6 - 648.d0*tmp_mu**8 + 2.718281828459045**(1/tmp_mu**2)*(-13.d0 -                            &
        540.d0*tmp_mu**4 - 960.d0*tmp_mu**6 + 648.d0*tmp_mu**8) + 957.12507948898*2.718281828459045**(1.d0/tmp_mu**2)*tmp_mu**3*(0.177777777777778      &
        + 1.d0*tmp_mu**2)*Erf(1.d0/tmp_mu)))/(E**(1.d0/tmp_mu**2)*tmp_c**2*(1.d0 + (18.d0 - 6.d0/E**(1.d0/tmp_mu**2))*tmp_mu**4 + 12.d0*(-1.d0          &
        + E**(-1.d0/tmp_mu**2))*tmp_mu**6 - 10.6347231054331*tmp_mu**3*Erf(1.d0/tmp_mu)))) 
  
  vx =  
  
 endif 
 end
