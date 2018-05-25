subroutine opt_orb
 implicit none
 if(n_det_generators_restart.gt.1)then
  print*, 'USING MCSCF ORBITAL OPTIMIZATION'
  double precision, allocatable :: psi_coef_normalized(:,:)
  double precision :: accu(N_states)
  integer :: i,j
  allocate(psi_coef_normalized(N_det_generators_restart, N_States))
  do i = 1, N_states
   accu(i) = 0.d0
   do j = 1, N_det_generators_restart
    accu(i) += psi_coef_generators_restart(j,i) ** 2
   enddo
   accu(i) = 1.d0/dsqrt(accu(i))
   do j = 1, N_det_generators_restart
    psi_coef_normalized(j,i) = psi_coef_generators_restart(j,i) * accu(i)
   enddo
  enddo
  call set_psi_det_as_input_psi(N_det_generators_restart,psi_det_generators_restart,psi_coef_normalized)  
  touch psi_det psi_coef 
  call diagonalize_CI 
  call clear_mo_map
  call casscf_routine
! call diag_inactive_virt_and_update_mos_MR_Fock  
 deallocate(psi_coef_normalized)
 else 
  call initialize_mo_coef_begin_iteration

  call damping_SCF
  print*, 'USING ROHF-LIKE  ORBITAL OPTIMIZATION'
  call diag_inactive_virt_and_update_mos_SR_Fock  
 endif
 call save_mos
end
