program pouet
 implicit none
 read_wf = .true.
 touch read_wf
 disk_access_mo_one_integrals = "None"
 touch disk_access_only_mo_one_integrals
 disk_access_mo_integrals = "None"
 touch disk_access_mo_integrals
 disk_access_ao_integrals = "None"
 touch disk_access_ao_integrals
 density_for_dft = "WFT" 
 touch density_for_dft 
 if(projected_wft_for_dft)then
  call print_projected_energy_dft
 else 
  call print_variational_energy_dft
 endif


end
