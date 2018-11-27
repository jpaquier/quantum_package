 BEGIN_PROVIDER [double precision, psi_kinetic_energy, (N_states) ]
&BEGIN_PROVIDER [double precision, psi_nuclear_elec_energy, (N_states) ]
 implicit none
 integer :: i,j,istate
 psi_kinetic_energy = 0.d0
 psi_nuclear_elec_energy = 0.d0 
 do istate = 1, N_states 
  do i = 1, mo_tot_num
   do j = 1, mo_tot_num
    psi_kinetic_energy(istate)      += ( one_body_dm_mo_alpha(j,i,istate)+one_body_dm_mo_beta(j,i,istate)) * mo_kinetic_integral(j,i) 
    psi_nuclear_elec_energy(istate) += ( one_body_dm_mo_alpha(j,i,istate)+one_body_dm_mo_beta(j,i,istate)) * mo_nucl_elec_integral(j,i) 
   enddo
  enddo
 enddo
END_PROVIDER 
