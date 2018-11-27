program write_effective_RSDFT_hamiltonian
 implicit none
 BEGIN_DOC
 ! This programs writes the effective RS-DFT Hamiltonian into the EZFIO folder. 
 ! The next programs that will run unto the EZFIO folder will, by default, have the one- and two-body integrals loaded from the EZFIO data. 
 END_DOC
 read_wf = .true.
 touch read_wf
 disk_access_mo_one_integrals = "None"
 touch disk_access_only_mo_one_integrals
 disk_access_mo_integrals = "None"
 touch disk_access_mo_integrals
 disk_access_ao_integrals = "None"
 touch disk_access_ao_integrals
 call routines_write_int
 call routines_compute_energy
end

subroutine routines_write_int
 implicit none
 call write_all_integrals_for_mrdft
 density_for_dft = "WFT" 
 touch density_for_dft 
end 
subroutine routines_compute_energy
 implicit none
 if(projected_wft_for_dft)then
  call print_projected_energy_dft
 else 
  call print_variational_energy_dft
 endif
 call ezfio_set_data_energy_and_density_data_one_body_alpha_dm_mo(one_body_dm_mo_alpha)
 call ezfio_set_data_energy_and_density_data_one_body_beta_dm_mo(one_body_dm_mo_beta)

end


