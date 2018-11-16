program write_integrals
 implicit none
 BEGIN_DOC
 ! This program saves the bielec erf integrals into the EZFIO folder
 END_DOC 
 disk_access_mo_integrals = "None"
 touch disk_access_mo_integrals
 disk_access_ao_integrals = "None"
 touch disk_access_ao_integrals
 call routine

end

subroutine routine
 implicit none
 call save_erf_bi_elec_integrals_ao
 call save_erf_bi_elec_integrals_mo

end

