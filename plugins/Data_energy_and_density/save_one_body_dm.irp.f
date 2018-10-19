program save_one_body_density
  implicit none
  read_wf = .True.
  touch read_wf
  call routine

end

subroutine routine
 
 call ezfio_set_data_energy_and_density_data_one_body_alpha_dm_mo(one_body_dm_mo_alpha)
 call ezfio_set_data_energy_and_density_data_one_body_beta_dm_mo(one_body_dm_mo_beta)
end
