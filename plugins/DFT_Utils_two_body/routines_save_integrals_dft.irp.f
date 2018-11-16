
subroutine save_one_e_effective_potential  
 implicit none
 BEGIN_DOC 
! used to save the effective_one_e_potential into the one-body integrals in the ezfio folder
! this effective_one_e_potential is computed with the current density 
! and will couple the WFT with DFT for the next regular WFT calculation
 END_DOC
 call write_one_e_integrals('mo_ne_integral', effective_one_e_potential_without_kin, &
      size(effective_one_e_potential_without_kin,1), size(effective_one_e_potential_without_kin,2))
 call write_one_e_integrals('mo_kinetic_integral',mo_kinetic_integral ,&
      size(mo_kinetic_integral,1), size(mo_kinetic_integral,2))

 print *,  'Effective DFT potential is written on disk on the mo_ne_integral integrals'
 call ezfio_set_integrals_monoelec_disk_access_mo_one_integrals("Read")

end

subroutine write_all_integrals_for_mrdft
 implicit none
 BEGIN_DOC
 ! saves all integrals needed for RS-DFT-MRCI calculation: one-body effective potential and two-elec erf integrals
 END_DOC
 call save_one_e_effective_potential
 call save_erf_bi_elec_integrals_ao_into_regular_integrals_ao
 call save_erf_bi_elec_integrals_mo_into_regular_integrals_mo
end

