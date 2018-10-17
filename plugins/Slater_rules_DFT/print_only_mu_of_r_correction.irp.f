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
!md_correlation_functional =  "basis_set_short_range_LDA"
!touch md_correlation_functional
 call print_contribution_dft_mu_of_r

end


