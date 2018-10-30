program write_integrals_for_dft
 implicit none
 read_wf = .true.
 touch read_wf
 disk_access_mo_one_integrals = "None"
 touch disk_access_mo_one_integrals
 disk_access_mo_integrals = "None"
 touch disk_access_mo_integrals
 disk_access_ao_integrals = "None"
 touch disk_access_ao_integrals
 call write_all_integrals_for_mrdft
 density_for_dft = "WFT" 
 touch density_for_dft 
 if(projected_wft_for_dft)then
  call print_projected_energy_dft
 else 
  call print_variational_energy_dft
 endif
 call ezfio_set_data_energy_and_density_data_one_body_alpha_dm_mo(one_body_dm_mo_alpha)
 call ezfio_set_data_energy_and_density_data_one_body_beta_dm_mo(one_body_dm_mo_beta)


end


