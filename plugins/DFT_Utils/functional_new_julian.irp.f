 BEGIN_PROVIDER [double precision, energy_Hxc, (N_states)]
 implicit none
 include 'Utils/constants.include.F'
 integer :: istate
 do istate = 1, N_states
  energy_Hxc(istate) =  pi/(2.d0 * mu_erf**2) * E_cor_tot(istate) 
 enddo
 END_PROVIDER
