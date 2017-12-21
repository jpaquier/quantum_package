program MRPT_Utils
  implicit none
  read_wf = .True.
  touch read_wf
! call routine
! call routine_2
  call routine_3
end

subroutine routine
 implicit none
 integer :: i
 do i = 1, mo_tot_num
  print*,i ,'',mo_class(i)
 enddo
 print*, 'core'
 do i = 1, n_core_orb
  print*, list_core(i)
 enddo
 call debug_det(core_bitmask, N_int)
 print*, 'inact'
 do i = 1, n_inact_orb
  print*, list_inact(i)
 enddo
 call debug_det(inact_bitmask, N_int)
 print*, 'act'
 do i = 1, n_act_orb
  print*, list_act(i)
 enddo
 call debug_det(act_bitmask, N_int)
 print*, 'virt'
 do i = 1, n_virt_orb
  print*, list_virt(i)
 enddo
 call debug_det(virt_bitmask, N_int)

end

subroutine routine_3
 implicit none
!provide fock_virt_total_spin_trace
 provide delta_ij_mrpt 
 
 print *,  'N_det    = ', N_det
 print *,  'N_states = ', N_states
 print *,  'PT2      = ', second_order_pt_new(1)
 print *,  'E        = ', CI_energy(1)
 print *,  'E+PT2    = ', CI_energy(1)+second_order_pt_new(1)
 print *,'****** DIAGONALIZATION OF DRESSED MATRIX ******'
 print *,  'E dressed= ', CI_dressed_pt2_new_energy(1)
 integer :: i
 do i = 1, N_det_ref
  write(*, '(2(F10.7,X))')psi_coef(i,1),CI_dressed_pt2_new_eigenvectors(i,1)
 enddo

end

subroutine routine_2
 implicit none
 integer :: i
 do i = 1, n_core_inact_orb
  print*,fock_core_inactive_total(i,1,1),fock_core_inactive(i)
 enddo
 double precision :: accu
 accu = 0.d0
 do i = 1, n_act_orb
  integer :: j_act_orb
  j_act_orb = list_act(i)
  accu += one_body_dm_mo_alpha(j_act_orb,j_act_orb,1)
  print*,one_body_dm_mo_alpha(j_act_orb,j_act_orb,1),one_body_dm_mo_beta(j_act_orb,j_act_orb,1) 
 enddo
 print*,'accu = ',accu

end


