subroutine Ec_sr_PBE(mu,rhoc,rhoo,sigmacc,sigmaco,sigmaoo,ec)
!************************************************************************
!     Short-range PBE correlation energy functional for erf interaction
!
!
!************************************************************************
include 'constants.include.F'
      implicit none
! input
      double precision, intent(in) ::  rhoc,rhoo,mu
      double precision, intent(in) ::  sigmacc,sigmaco,sigmaoo
! output
      double precision, intent(out) ::  ec
! local
      double precision tol
      parameter(tol=1d-12)

      character(len=30) namedummy

      double precision eccerflda
      double precision vrhoccerflda
      double precision vrhoocerflda

      double precision ecclda
      double precision vrhocclda
      double precision vrhooclda

      integer i,igrad
      double precision rho,drho2,rhoa,rhob
      double precision ecerflda
      double precision eclda,decldadrho
      double precision ecerfpbe
      double precision arglog,alpha,beta,gamma
      double precision Aa,Ab,Ac,tq
      double precision zeta,phi,phi2,phi3,phi4
      double precision, parameter :: f13=0.333333333333333d0


! Parameter of the modified interaction

      ec = 0.d0

! First-type gradient functional
     igrad=1

     alpha=2.78d0
     gamma=3.1091d-2

!    test on density
     if (dabs(rhoc).lt.tol) return
     double precision :: vc_a,vc_b
!    Spin polarisation
     rhoa=max((rhoc+rhoo)*.5d0,1.0d-15)
     rhob=max((rhoc-rhoo)*.5d0,1.0d-15)

     call ec_lda_sr(mu,rhoa,rhob,eccerflda,vc_a,vc_b)
     ecerflda = eccerflda
     vrhoccerflda = 0.5d0 * (vc_a + vc_b)
     vrhoocerflda = 0.5d0 * (vc_a - vc_b)

!    Density
     rho = rhoc
     rho = max(rho,1.d-10)

!    Square of density gradient
     drho2 = sigmacc

     zeta = (rhoa-rhob)/(rhoa+rhob)
     zeta = max(zeta,1.d-10)

!    LDA energy density
     double precision :: vc_a_lda,vc_b_lda
     call ec_lda(rhoa,rhob,ecclda,vc_a_lda,vc_b_lda)
     eclda = ecclda
     decldadrho = 0.5d0 * (vc_a_lda+vc_b_lda)
     decldadrho = 0.5d0 * (vc_a_lda-vc_b_lda)

     if ((ecerflda/eclda).le.0d0) then
        beta=0d0
     else
        beta=6.6725d-2*(ecerflda/eclda)**alpha
     endif
     phi=((1d0+zeta)**(2d0/3d0)+(1d0-zeta)**(2d0/3d0))/2d0
     phi2=phi*phi
     phi3=phi2*phi
     phi4=phi3*phi
     tq=drho2*6.346820607d-2*rho**(-7d0/3d0)/phi2
!    tq=drho2*6.346820607d0-2*rho*rho*rho**(0.3333333333333d0)
!    tq=drho2*6.346820607d0-2*rho*rho!*rho**(0.3333333333333d0)
     Ab=dexp(-ecerflda/(rho*gamma*phi3))-1d0
     if (dabs(Ab).le.dabs(beta*tol)) then
        ecerfpbe=ecerflda
     else
        Aa=beta/(gamma*Ab)
        Ac=1d0+Aa*tq+Aa**2*tq**2
        if (Aa.lt.tol) Aa=tol
        arglog=1d0+beta*(1d0-1d0/Ac)/(gamma*Aa)
        arglog=max(arglog,1.d-10)
        ecerfpbe=ecerflda+rho*phi3*gamma*dlog(arglog)
     end if

     ec = ecerfpbe

end

subroutine numerical_derivative_of_sr_pbe_correlation (mu,rhoc,rhoo,sigmacc,sigmaco,sigmaoo,vrhoc,vrhoo,vsigmacc,vsigmaco,vsigmaoo)
implicit none
double precision, intent (in) :: rhoc,rhoo,sigmacc,sigmaco,sigmaoo,mu
double precision, intent (out) :: vrhoc,vrhoo,vsigmacc,vsigmaco,vsigmaoo
double precision :: delta_rhoc,delta_rhoo,delta_sigmacc,delta_sigmaco,delta_sigmaoo
double precision :: e_plus,e_minus 
delta_rhoc =    max(1d-15,1d-5*rhoc)
delta_rhoo =    max(1d-15,1d-5*rhoo)
delta_sigmacc = max(1d-15,1d-5*sigmacc)
delta_sigmaco = max(1d-15,1d-5*sigmaco)
delta_sigmaoo = max(1d-15,1d-5*sigmaoo)

e_plus=0d0
e_minus=0d0
call Ec_sr_PBE(mu,rhoc+delta_rhoc,rhoo,sigmacc,sigmaco,sigmaoo,e_plus)
call Ec_sr_PBE(mu,rhoc-delta_rhoc,rhoo,sigmacc,sigmaco,sigmaoo,e_minus)
vrhoc = 0.5d0*(e_plus - e_minus)/delta_rhoc

e_plus=0d0
e_minus=0d0
call Ec_sr_PBE(mu,rhoc,rhoo+delta_rhoo,sigmacc,sigmaco,sigmaoo,e_plus)
call Ec_sr_PBE(mu,rhoc,rhoo-delta_rhoo,sigmacc,sigmaco,sigmaoo,e_minus)
vrhoo = 0.5d0*(e_plus - e_minus)/delta_rhoo

e_plus=0d0
e_minus=0d0
call Ec_sr_PBE(mu,rhoc,rhoo,sigmacc+delta_sigmacc,sigmaco,sigmaoo,e_plus)
call Ec_sr_PBE(mu,rhoc,rhoo,sigmacc-delta_sigmacc,sigmaco,sigmaoo,e_minus)
vsigmacc = 0.5d0*(e_plus - e_minus)/delta_sigmacc

e_plus=0d0
e_minus=0d0
call Ec_sr_PBE(mu,rhoc,rhoo,sigmacc,sigmaco+delta_sigmaco,sigmaoo,e_plus)
call Ec_sr_PBE(mu,rhoc,rhoo,sigmacc,sigmaco-delta_sigmaco,sigmaoo,e_minus)
vsigmaco = 0.5d0*(e_plus - e_minus)/delta_sigmaco

e_plus=0d0
e_minus=0d0
call Ec_sr_PBE(mu,rhoc,rhoo,sigmacc,sigmaco,sigmaoo+delta_sigmaoo,e_plus)
call Ec_sr_PBE(mu,rhoc,rhoo,sigmacc,sigmaco,sigmaoo-delta_sigmaoo,e_minus)
vsigmaoo = 0.5d0*(e_plus - e_minus)/delta_sigmaoo

end
