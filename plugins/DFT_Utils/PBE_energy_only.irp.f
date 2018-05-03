subroutine Energy_x_pbe_sr(mu,rho_a,rho_b,grd_rho_a_2,grd_rho_b_2,grd_rho_a_b,ex)
BEGIN_DOC
!rho_a = density alpha
!rho_b = density beta
!grd_rho_a_2 = (gradient rho_a)^2
!grd_rho_b_2 = (gradient rho_b)^2
!grd_rho_a_b = (gradient rho_a).(gradient rho_b)
!ex = exchange energy density at point r
END_DOC

 implicit none

! input
 double precision, intent(in) :: mu,rho_a, rho_b
 double precision, intent(in) :: grd_rho_a_2, grd_rho_b_2, grd_rho_a_b

! output
 double precision, intent(out) :: ex

! function
  double precision berf

! local
  double precision, parameter :: tol=1d-12
  double precision, parameter :: f13=0.333333333333333d0

  double precision exerflda,vxerflda_a,vxerflda_b
  double precision exerfpbe_a, exerfpbe_b

  double precision rho,drho2
  double precision kappa,sq,fx


! initialization
  ex=0.d0

  
! spin scaling relation Ex[rho_a,rho_b] = (1/2) (Ex[2rho_a,2rho_a] + Ex[2rho_b,2rho_b])

! two times spin alpha density
  rho = max(rho_a,tol)*2.d0

! test on density
  if (rho >= tol) then

!  call srLDA Ex[2*rho_a,2*rho_a]
   call ex_lda_sr(mu,rho_a,rho_a,exerflda,vxerflda_a,vxerflda_b)

!  square of two times spin alpha density gradient
   drho2=max(grd_rho_a_2,0d0)*4.0d0

   kappa=0.804d0
   sq=drho2*2.6121172985233599567768d-2*rho**(-8d0/3d0)
   fx=1d0+kappa-kappa/(1d0+berf(1.616204596739954813d-1*mu*rho**(-f13))*sq/kappa)
   exerfpbe_a=exerflda*fx

 endif
   

! two times spin beta density
  rho = max(rho_b,tol)*2.d0

! test on density
  if (rho >= tol) then

!  call srLDA Ex[2*rho_b,2*rho_b]
   call ex_lda_sr(mu,rho_b,rho_b,exerflda,vxerflda_a,vxerflda_b)

!  square of two times spin beta density gradient
   drho2=max(grd_rho_b_2,0d0)*4.0d0

   kappa=0.804d0
   sq=drho2*2.6121172985233599567768d-2*rho**(-8d0/3d0)
   fx=1d0+kappa-kappa/(1d0+berf(1.616204596739954813d-1*mu*rho**(-f13))*sq/kappa)
   exerfpbe_b=exerflda*fx

  endif

  ex = (exerfpbe_a+exerfpbe_b)*0.5d0

 end



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
