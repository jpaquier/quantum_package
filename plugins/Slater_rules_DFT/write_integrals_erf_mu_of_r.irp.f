program write_integrals
 implicit none
 disk_access_mo_one_integrals = "None"
 touch disk_access_mo_one_integrals
 disk_access_mo_integrals = "None"
 touch disk_access_mo_integrals
 disk_access_ao_integrals = "None"
 touch disk_access_ao_integrals
 call routine

end

subroutine routine
 implicit none
 call save_erf_mu_of_r_bi_elec_integrals_ao
 call save_erf_mu_of_r_bi_elec_integrals_mo

end

