
subroutine save_one_e_effective_potential  
 implicit none
 BEGIN_DOC 
! used to save the effective_one_e_potential into the one-body integrals in the ezfio folder
! this effective_one_e_potential is computed with the current density 
! and will couple the WFT with DFT for the next regular WFT calculation
 END_DOC
 double precision, allocatable :: tmp(:,:),tmp_2(:,:)
 allocate(tmp(size(effective_one_e_potential,1),size(effective_one_e_potential,2)),tmp_2(mo_tot_num,mo_tot_num))
 integer :: i,j
 do i = 1, mo_tot_num
  do j = 1, mo_tot_num
   tmp(i,j) = effective_one_e_potential(i,j,1)
  enddo
 enddo
 call write_one_e_integrals('mo_ne_integral', tmp, &
      size(tmp,1), size(tmp,2))
 print *,  'Effective DFT potential is written on disk on the mo_ne_integral integrals'
 call ezfio_set_integrals_monoelec_disk_access_mo_one_integrals("Read")

 tmp_2 = 0.d0
 call write_one_e_integrals('mo_kinetic_integral', tmp_2,&
      size(tmp_2,1), size(tmp_2,2))
 print *,  'MO kinetic integrals written and set to ZERO'


 deallocate(tmp,tmp_2)

end

subroutine save_erf_bi_elec_integrals_mo
 implicit none
 integer :: i,j,k,l
 PROVIDE mo_bielec_integrals_erf_in_map
 call ezfio_set_work_empty(.False.)
 call map_save_to_disk(trim(ezfio_filename)//'/work/mo_ints',mo_integrals_erf_map)
 call ezfio_set_integrals_bielec_disk_access_mo_integrals("Read")
end

subroutine save_erf_bi_elec_integrals_ao
 implicit none
 integer :: i,j,k,l
 PROVIDE ao_bielec_integrals_erf_in_map
 call ezfio_set_work_empty(.False.)
 call map_save_to_disk(trim(ezfio_filename)//'/work/ao_ints',ao_integrals_erf_map)
 call ezfio_set_integrals_bielec_disk_access_ao_integrals("Read")
end

