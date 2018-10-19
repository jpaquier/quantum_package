program save_auxiliary_mos
  implicit none
  read_wf = .True. 
  touch read_wf
  call routine

end


subroutine routine
 implicit none
 call ezfio_set_data_energy_and_density_data_one_body_alpha_dm_mo(dm_alpha_on_natorb_basis)
 call ezfio_set_data_energy_and_density_data_one_body_beta_dm_mo(dm_beta_on_natorb_basis)
end
