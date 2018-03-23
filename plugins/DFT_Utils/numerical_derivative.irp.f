subroutine E_sr_PBE(rhoc,rhoo,sigmacc,sigmaco,sigmaoo,ec)
!************************************************************************
!     Short-range PBE correlation energy functional for erf interaction
!
!
!************************************************************************
include 'constants.include.F'
      implicit none
! input
      double precision, intent(in) ::  rhoc,rhoo
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
      double precision mu
      double precision rho,drho2,rhoa,rhob
      double precision ecerflda
      double precision eclda,decldadrho
      double precision ecerfpbe
      double precision arglog,alpha,beta,gamma
      double precision Aa,Ab,Ac,tq
      double precision zeta,phi,phi2,phi3,phi4
      double precision, parameter :: f13=0.333333333333333d0


! Parameter of the modified interaction
      mu = mu_erf

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

     call ec_lda_sr(rhoa,rhob,eccerflda,vc_a,vc_b)
     ecerflda = eccerflda
     vrhoccerflda = 0.5d0 * (vc_a + vc_b)
     vrhoocerflda = 0.5d0 * (vc_a - vc_b)

!    Density
     rho = rhoc

!    Square of density gradient
     drho2 = sigmacc

     zeta = (rhoa-rhob)/(rhoa+rhob)

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
     Ab=dexp(-ecerflda/(rho*gamma*phi3))-1d0
     if (dabs(Ab).le.dabs(beta*tol)) then
        ecerfpbe=ecerflda
     else
        Aa=beta/(gamma*Ab)
        Ac=1d0+Aa*tq+Aa**2*tq**2
        if (Aa.lt.tol) Aa=tol
        arglog=1d0+beta*(1d0-1d0/Ac)/(gamma*Aa)
        ecerfpbe=ecerflda+rho*phi3*gamma*dlog(arglog)
     end if

     ec = ecerfpbe

end

subroutine numerical_derivative_of_sr_pbe_correlation (rhoc,rhoo,sigmacc,sigmaco,sigmaoo,vrhoc,vrhoo,vsigmacc,vsigmaco,vsigmaoo)
implicit none
double precision, intent (in) :: rhoc,rhoo,sigmacc,sigmaco,sigmaoo
double precision, intent (out) :: vrhoc,vrhoo,vsigmacc,vsigmaco,vsigmaoo
double precision :: delta_rho 
double precision :: e_plus,e_minus 
delta_rho = 1d-5

e_plus=0d0
e_minus=0d0
call E_sr_PBE(rhoc+delta_rho,rhoo,sigmacc,sigmaco,sigmaoo,e_plus)
call E_sr_PBE(rhoc-delta_rho,rhoo,sigmacc,sigmaco,sigmaoo,e_minus)
vrhoc = 0.5d0*(e_plus - e_minus)/delta_rho

e_plus=0d0
e_minus=0d0
call E_sr_PBE(rhoc,rhoo+delta_rho,sigmacc,sigmaco,sigmaoo,e_plus)
call E_sr_PBE(rhoc,rhoo-delta_rho,sigmacc,sigmaco,sigmaoo,e_minus)
vrhoo = 0.5d0*(e_plus - e_minus)/delta_rho

e_plus=0d0
e_minus=0d0
call E_sr_PBE(rhoc,rhoo,sigmacc+delta_rho,sigmaco,sigmaoo,e_plus)
call E_sr_PBE(rhoc,rhoo,sigmacc-delta_rho,sigmaco,sigmaoo,e_minus)
vsigmacc = 0.5d0*(e_plus - e_minus)/delta_rho

e_plus=0d0
e_minus=0d0
call E_sr_PBE(rhoc,rhoo,sigmacc,sigmaco+delta_rho,sigmaoo,e_plus)
call E_sr_PBE(rhoc,rhoo,sigmacc,sigmaco-delta_rho,sigmaoo,e_minus)
vsigmaco = 0.5d0*(e_plus - e_minus)/delta_rho

e_plus=0d0
e_minus=0d0
call E_sr_PBE(rhoc,rhoo,sigmacc,sigmaco,sigmaoo+delta_rho,e_plus)
call E_sr_PBE(rhoc,rhoo,sigmacc,sigmaco,sigmaoo-delta_rho,e_minus)
vsigmaoo = 0.5d0*(e_plus - e_minus)/delta_rho

print*, vrhoc,vrhoo,vsigmacc,vsigmaco,vsigmaoo
end
