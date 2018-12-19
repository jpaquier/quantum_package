 subroutine dirac_ex_lda(rho,ex,vx)
 include 'constants.include.F'
 implicit none
 double precision, intent(in) :: rho
 double precision, intent(out) :: ex,vx
 double precision :: kF,tmp_c
 kF = (3* pi**2 *rho)**(c_1_3)
 tmp_c = speed_of_light/kF
 ex =    0.001343813934716645*kF**4*(-2.d0 + tmpc**2 - 2.d0*Sqrt(1.d0 + tmpc**2)*(2.d0 + 3.d0*tmpc**2)*ArcCsch(tmpc)    &
         + 3.d0*tmpc**4*ArcCsch(tmpc)**2 + 2.d0*(1.d0 + tmpc**2)**2*Log(1.d0 + tmpc**(-2)))

 vx =   (0.1061032953945969*kF*(-1.d0*Sqrt(1.d0 + tmpc**2)*(2.d0 + 3.d0*tmpc**2)*ArcCsch(tmpc) +                        &
        (1.d0 + tmpc**2)*(-1.d0 + (1.d0 + tmpc**2)*Log(1.d0 + tmpc**(-2)))))/(1.d0 + tmpc**2)
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
  ex =  (0.00004479379782388818*kF**4*(60.d0*(-3.d0 + 2.d0*tmpmu**2*(-3.d0 + tmpmu**2 - (1.d0*(-2.d0 + tmpmu**2))*dexp(-1.d0/tmpmu**2)) +
        7.089815403622064*tmpmu*Erf(1.d0/tmpmu)) + (20.d0*(1.d0 + 6.d0*tmpmu**4*(3.d0 - 2.d0*tmpmu**2 + (-1.d0 + 2.d0*tmpmu**2)*dexp(-1.d0/tmpmu**2)) -
        10.6347231054331*tmpmu**3*Erf(1.d0/tmpmu))**2 - (3.d0*(-2.d0*tmpmu**2*(-2.d0 + tmpmu**2) + dexp(1.d0/tmpmu**2)*(-3.d0 - 6.d0*tmpmu**2 + 
        2.d0*tmpmu**4 + 7.089815403622064*tmpmu*Erf(1.d0/tmpmu)))*(24.d0*tmpmu**4*(4.d0 + 13.d0*tmpmu**2 - 27.d0*tmpmu**4) + dexp(1.d0/tmpmu**2)*(-13.d0 - 
        540.d0*tmpmu**4 - 960.d0*tmpmu**6 + 648.d0*tmpmu**8 + 21.26944621086619*tmpmu**3*(8.d0 + 45.d0*tmpmu**2)*Erf(1.d0/tmpmu))))*dexp(-2.d0/tmpmu**2))/
        (tmpc**2*(1.d0 + 6.d0*tmpmu**4*(3.d0- 2.d0*tmpmu**2 + (-1.d0 + 2.d0*tmpmu**2)*dexp(-1.d0/tmpmu**2)) - 10.6347231054331*tmpmu**3*Erf(1.d0/tmpmu)))))/
        (1.d0 + (0.05d0*(24.d0*tmpmu**4*(4.d0 + 13.d0*tmpmu**2 - 27.d0*tmpmu**4) + dexp(1.d0/tmpmu**2)*(-13.d0 - 540.d0*tmpmu**4 - 960.d0*tmpmu**6 +
        648.d0*tmpmu**8 + 21.26944621086619*tmpmu**3*(8.d0 + 45.d0*tmpmu**2)*Erf(1.d0/tmpmu))))/(tmpc**2*(6.d0*tmpmu**4 - 12.d0*tmpmu**6 +      
        dexp(1.d0/tmpmu**2)*(-1.d0 - 18.d0*tmpmu**4 + 12.d0*tmpmu**6 + 10.6347231054331*tmpmu**3*Erf(1.d0/tmpmu))))) 
  
  vx =  
  
 endif 
 end
