
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
