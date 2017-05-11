BEGIN_PROVIDER [ double precision, psi_energy_erf, (N_states) ]
&BEGIN_PROVIDER [ double precision, psi_energy_core, (N_states) ]
&BEGIN_PROVIDER [ double precision, psi_energy_hartree, (N_states) ]
&BEGIN_PROVIDER [ double precision, psi_energy_core_and_hartree, (N_states) ]
&BEGIN_PROVIDER [ double precision, total_electronic_energy, (N_states) ]
  implicit none
  integer :: i,j
  double precision :: array(mo_tot_num_align,mo_tot_num),average
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
  psi_energy_hartree = 0.d0
  psi_energy_core = 0.d0
  psi_energy_core_and_hartree = psi_energy_monoelec_dft 
  call u_0_H_u_0_erf(psi_energy_erf,psi_coef,N_det,psi_det,N_int,N_states,psi_det_size)
! total_electronic_energy = psi_energy_hartree + psi_energy_core + psi_energy_erf + energy_x + energy_c
  total_electronic_energy = psi_energy_core_and_hartree + psi_energy_erf + energy_x + energy_c
END_PROVIDER

