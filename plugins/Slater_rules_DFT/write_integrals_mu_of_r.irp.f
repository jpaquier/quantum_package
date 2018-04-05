program write_integrals_for_dft
 implicit none
 read_wf = .true.
 touch read_wf
 disk_access_mo_one_integrals = "None"
 touch disk_access_only_mo_one_integrals
 disk_access_mo_integrals = "None"
 touch disk_access_mo_integrals
 disk_access_ao_integrals = "None"
 touch disk_access_ao_integrals
 mu_erf = 1.d+10
 touch mu_erf 
 exchange_functional =  "basis_set_short_range_LDA"
 touch exchange_functional
 correlation_functional =  "basis_set_short_range_LDA"
 touch correlation_functional
 md_correlation_functional =  "basis_set_short_range_LDA"
 touch md_correlation_functional
 call write_all_integrals_for_mrdft
 call print_variational_energy_dft

end


