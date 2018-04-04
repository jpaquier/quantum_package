subroutine dftfun_ecerfpbe(mu,rhoc,rhoo,sigmacc,sigmaco,sigmaoo,ec,vrhoc,vrhoo,vsigmacc,vsigmaco,vsigmaoo)
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
      double precision, intent(out) ::  vrhoc,vrhoo
      double precision, intent(out) ::  vsigmacc,vsigmaco,vsigmaoo
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
      double precision ecerflda,decerfldadrho
      double precision eclda,decldadrho
      double precision ecerfpbe,decerfpbedrho,decerfpbedrhoo
      double precision decerfpbeddrho2
      double precision arglog,arglogs,arglogss,alpha,beta,betas,gamma
      double precision Aa,Ab,Ac,Aas,tq,tqs,tqss,decerfpur,decpur
      double precision t1,t2,t3,t4,t5,t6,t7,t8,t9,t10
      double precision t11,t12,t13,t14,t15,t16,t17,t18,t19
      double precision zeta,phi,phi2,phi3,phi4,phis,arglogsc
      double precision dlogarglog
      double precision, parameter :: f13=0.333333333333333d0


! Parameter of the modified interaction

      ec = 0.d0
      vrhoc  = 0.d0
      vrhoo = 0.d0
      vsigmacc = 0.d0
      vsigmaco = 0.d0
      vsigmaoo = 0.d0

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

!    Square of density gradient
     drho2 = sigmacc

     zeta = (rhoa-rhob)/(rhoa+rhob)

!    LDA energy density
     double precision :: vc_a_lda,vc_b_lda
     call ec_lda(rhoa,rhob,ecclda,vc_a_lda,vc_b_lda)
     eclda = ecclda
    !decldadrho = 0.5d0 * (vc_a_lda+vc_b_lda)
    !decldadrho = 0.5d0 * (vc_a_lda-vc_b_lda)

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

!     if(ldebug) write(*,*)"ecerfpbe=",ecerfpbe

! Derive


!    LDA energy density derivative
     decerfldadrho = vrhoccerflda
     decldadrho = 0.5d0 * (vc_a_lda+vc_b_lda)

     decerfpur=(decerfldadrho-ecerflda/rho)/rho
     decpur=(decldadrho-eclda/rho)/rho
     betas=alpha*beta*(decerfpur*rho/ecerflda-decpur*rho/eclda)
     phis=((rhoa - rhob)*((rhoa/(rhoa + rhob))**f13 - (rhob/(rhoa + rhob))**f13))/(3d0*2d0**f13*(rhoa/(rhoa + rhob))**f13*(rhob/(rhoa + rhob))**f13*(rhoa + rhob)**2)
     if (dabs(Ab).le.dabs(beta*tol)) then
        decerfpbedrho=decerfldadrho
     else
        Aas=betas/(gamma*Ab)+Aa*(1d0+1d0/Ab)*(decerfpur/phi3-3d0*phis*ecerflda/(rho*phi4))/gamma
        tqs=-7d0*tq/(3d0*rho)-2d0*tq*phis/phi
        arglogs=betas*tq*(1d0+Aa*tq)/(Ac*gamma)+beta*tqs*(1d0+Aa*tq)/(Ac*gamma)-beta*tq*Aa*tq*(Aas*tq+Aa*tqs)*(2d0+Aa*tq)/(Ac**2*gamma)
        dlogarglog=dlog(arglog)
        decerfpbedrho=decerfldadrho+gamma*(phi3*dlogarglog+3d0*rho*phis*phi2*dlogarglog+rho*phi3*arglogs/arglog)
     end if

     if (dabs(Ab).le.dabs(beta*tol)) then
        decerfpbeddrho2=0.0d0
     else
        arglogsc=Ab*(Aa+2d0*Aa*Aa*tq)/(Ac*Ac)
        tqss=6.346820607d-2*rho**(-7d0/3d0)/phi2
        arglogss=tqss*arglogsc
        decerfpbeddrho2=rho*gamma*phi3*arglogss/arglog
     end if

!    LDA energy density derivative
     decerfldadrho = vrhoocerflda
     decldadrho = 0.5d0 * (vc_a_lda-vc_b_lda)

     decerfpur=decerfldadrho/rho
     decpur=decldadrho/rho
     betas=alpha*beta*(decerfpur*rho/ecerflda-decpur*rho/eclda)
     phis=(rhob*(rhoa/(rhoa + rhob))**(2d0*f13)-rhoa*(rhob/(rhoa + rhob))**(2d0*f13))/(3d0*2d0**f13*rhoa*rhob)

     if (dabs(Ab).le.dabs(beta*tol)) then
        decerfpbedrhoo=decerfldadrho
     else
        Aas=betas/(gamma*Ab)+Aa*(1d0+1d0/Ab)*(decerfpur/phi3-3d0*phis*ecerflda/(rho*phi4))/gamma
        tqs=-2d0*tq*phis/phi
        arglogs=betas*tq*(1d0+Aa*tq)/(Ac*gamma)+beta*tqs*(1d0+Aa*tq)/(Ac*gamma)-beta*tq*Aa*tq*(Aas*tq+Aa*tqs)*(2d0+Aa*tq)/(Ac**2*gamma)
        decerfpbedrhoo=decerfldadrho+gamma*(3d0*rho*phis*phi2*dlog(arglog)+rho*phi3*arglogs/arglog)
     end if

!    derivatives
     vrhoc = vrhoc + decerfpbedrho
     vrhoo = vrhoo + decerfpbedrhoo
     vsigmacc = vsigmacc + decerfpbeddrho2


end
