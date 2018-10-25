program write_integrals_for_dft_ecmd_lda
 implicit none
 read_wf = .true.
 touch read_wf
 disk_access_mo_one_integrals = "None"
 touch disk_access_only_mo_one_integrals
 disk_access_mo_integrals = "None"
 touch disk_access_mo_integrals
 disk_access_ao_integrals = "None"
 touch disk_access_ao_integrals
 call write_all_integrals_for_mrdft_ecmd_lda
 if(projected_wft_for_dft)then
  call print_projected_energy_dft
 else 
  call print_variational_energy_dft_mu_of_r
 endif

end


