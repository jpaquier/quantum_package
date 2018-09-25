
BEGIN_PROVIDER [ double precision, psi_energy_erf_mu_of_r, (N_states) ]
  use bitmasks
  implicit none
  BEGIN_DOC
  ! Computes e_0 = <Psi|W_{ee}^{lr}|Psi>/<Psi|Psi>
  !
  END_DOC
  integer :: i
  call u_0_H_u_0_erf_mu_of_r(psi_energy_erf_mu_of_r,psi_coef,N_det,psi_det,N_int,N_states,psi_det_size)
  do i=N_det+1,N_states
    psi_energy_erf_mu_of_r(i) = 0.d0
  enddo
END_PROVIDER


BEGIN_PROVIDER [ double precision, psi_energy_core_and_sr_mu_of_r_hartree, (N_states) ]
  implicit none
  BEGIN_DOC
! psi_energy_core                = <Psi| h_{core} + v_{H}^{sr}|Psi>
  END_DOC
  psi_energy_core_and_sr_mu_of_r_hartree = psi_energy_core + short_range_Hartree
END_PROVIDER


BEGIN_PROVIDER [ double precision, total_range_separated_electronic_energy_mu_of_r, (N_states) ]
  implicit none
  BEGIN_DOC
! total_range_separated_electronic_energy = <Psi| h_{core} |Psi> + (1/2) <Psi| v_{H}^{sr} |Psi> + <i|W_{ee}^{lr}|i> + E_{x} + E_{c}
  END_DOC
  total_range_separated_electronic_energy = psi_energy_core + short_range_Hartree_mu_of_r + psi_energy_erf_mu_of_r + energy_x + energy_c
END_PROVIDER


BEGIN_PROVIDER [ double precision, two_elec_energy_mu_of_r_dft, (N_states) ]
  implicit none
  BEGIN_DOC
! two_elec_energy_dft = (1/2) <Psi| v_{H}^{sr} |Psi> + <i|W_{ee}^{lr}|i> 
  END_DOC
   two_elec_energy_dft = short_range_Hartree + psi_energy_erf_mu_of_r
END_PROVIDER


 BEGIN_PROVIDER [ double precision, ref_bitmask_energy_erf_mu_of_r ]
&BEGIN_PROVIDER [ double precision, bi_elec_ref_bitmask_energy_erf_mu_of_r ]
  use bitmasks
  implicit none
  BEGIN_DOC
  ! Energy of the reference bitmask used in Slater rules
  END_DOC
  
  integer                        :: occ(N_int*bit_kind_size,2)
  integer                        :: i,j
  
  call bitstring_to_list(ref_bitmask(1,1), occ(1,1), i, N_int)
  call bitstring_to_list(ref_bitmask(1,2), occ(1,2), i, N_int)
  
  
  ref_bitmask_energy_erf_mu_of_r = 0.d0
  bi_elec_ref_bitmask_energy_erf_mu_of_r = 0.d0
  
  do j= 1, elec_alpha_num
    do i = j+1, elec_alpha_num
      bi_elec_ref_bitmask_energy_erf_mu_of_r += mo_bielec_integral_erf_mu_of_r_jj_anti(occ(i,1),occ(j,1))
      ref_bitmask_energy_erf_mu_of_r += mo_bielec_integral_erf_mu_of_r_jj_anti(occ(i,1),occ(j,1))
    enddo
  enddo
  
  do j= 1, elec_beta_num
    do i = j+1, elec_beta_num
      bi_elec_ref_bitmask_energy += mo_bielec_integral_erf_mu_of_r_jj_anti(occ(i,2),occ(j,2))
      ref_bitmask_energy += mo_bielec_integral_erf_mu_of_r_jj_anti(occ(i,2),occ(j,2))
    enddo
    do i= 1, elec_alpha_num
      bi_elec_ref_bitmask_energy += mo_bielec_integral_erf_mu_of_r_jj(occ(i,1),occ(j,2))
      ref_bitmask_energy += mo_bielec_integral_erf_mu_of_r_jj(occ(i,1),occ(j,2))
    enddo
  enddo
  
END_PROVIDER

