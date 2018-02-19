BEGIN_PROVIDER [double precision, elec_energy_dft_pure_corr_funct, (N_states)]
 implicit none
 BEGIN_DOC
! pure_corr_mr_dft_energy = <\Psi| H |\Psi> + E_{c,md,sr}
 END_DOC
 elec_energy_dft_pure_corr_funct = psi_energy + Energy_c_md

END_PROVIDER 
