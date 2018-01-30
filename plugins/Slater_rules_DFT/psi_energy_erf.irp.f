BEGIN_PROVIDER [ double precision, psi_energy_erf, (N_states) ]
&BEGIN_PROVIDER [ double precision, psi_energy_core, (N_states) ]
&BEGIN_PROVIDER [ double precision, psi_energy_hartree, (N_states) ]
&BEGIN_PROVIDER [ double precision, psi_energy_core_and_hartree, (N_states) ]
&BEGIN_PROVIDER [ double precision, total_electronic_energy, (N_states) ]
  implicit none
  integer :: i,j
  double precision :: array(mo_tot_num,mo_tot_num),average
  BEGIN_DOC
! Energy of the current wave function
  END_DOC
! array = mo_nucl_elec_integral + mo_kinetic_integral
! call get_average(array,one_body_dm_mo,average)
! psi_energy_core = average
! do i = 1, mo_tot_num
!  do j = 1, mo_tot_num
!   array(i,j) = short_range_Hartree_operator(i,j)
!  enddo
! enddo
! call get_average(array,one_body_dm_mo,average)
! psi_energy_hartree = 0.5d0 * average


! psi_energy_hartree = 0.d0
! psi_energy_core = 0.d0
! psi_energy_core_and_hartree = psi_energy_monoelec_dft 
! call u_0_H_u_0_erf(psi_energy_erf,psi_coef,N_det,psi_det,N_int,N_states,psi_det_size)
  double precision :: hij_core, hij_hartree, hij_erf, contrib
  psi_energy_hartree = 0.d0
  psi_energy_core = 0.d0
  psi_energy_erf = 0.d0
  do i = 1, N_Det
   do j = 1, N_det
    contrib = psi_coef(i,1) * psi_coef(j,1)
    call i_H_j_dft_general(psi_det(1,1,i),psi_det(1,1,j),N_int,hij_core, hij_hartree, hij_erf)
    psi_energy_erf += contrib * hij_erf
    psi_energy_core += contrib * hij_core 
    psi_energy_hartree += contrib * hij_hartree
!   print*, psi_energy_erf,psi_energy_core,psi_energy_hartree
   enddo
  enddo
  psi_energy_core_and_hartree = psi_energy_core + psi_energy_hartree
! total_electronic_energy = psi_energy_hartree + psi_energy_core + psi_energy_erf + energy_x + energy_c
  total_electronic_energy = psi_energy_core_and_hartree + psi_energy_erf + energy_x + energy_c
END_PROVIDER

 BEGIN_PROVIDER [ double precision, ref_bitmask_energy_erf ]
&BEGIN_PROVIDER [ double precision, bi_elec_ref_bitmask_energy_erf ]
  use bitmasks
  implicit none
  BEGIN_DOC
  ! Energy of the reference bitmask used in Slater rules
  END_DOC
  
  integer                        :: occ(N_int*bit_kind_size,2)
  integer                        :: i,j
  
  call bitstring_to_list(ref_bitmask(1,1), occ(1,1), i, N_int)
  call bitstring_to_list(ref_bitmask(1,2), occ(1,2), i, N_int)
  
  
  ref_bitmask_energy_erf = 0.d0
  bi_elec_ref_bitmask_energy_erf = 0.d0
  
  do j= 1, elec_alpha_num
    do i = j+1, elec_alpha_num
      bi_elec_ref_bitmask_energy_erf += mo_bielec_integral_erf_jj_anti(occ(i,1),occ(j,1))
      ref_bitmask_energy_erf += mo_bielec_integral_erf_jj_anti(occ(i,1),occ(j,1))
    enddo
  enddo
  
  do j= 1, elec_beta_num
    do i = j+1, elec_beta_num
      bi_elec_ref_bitmask_energy += mo_bielec_integral_erf_jj_anti(occ(i,2),occ(j,2))
      ref_bitmask_energy += mo_bielec_integral_erf_jj_anti(occ(i,2),occ(j,2))
    enddo
    do i= 1, elec_alpha_num
      bi_elec_ref_bitmask_energy += mo_bielec_integral_erf_jj(occ(i,1),occ(j,2))
      ref_bitmask_energy += mo_bielec_integral_erf_jj(occ(i,1),occ(j,2))
    enddo
  enddo
  
END_PROVIDER

