program print_1h2p
 implicit none
 read_wf = .True.
 touch read_wf
 call routine_2h1p
end

subroutine routine_2h1p
 implicit none
 double precision,allocatable :: matrix_1h2p(:,:,:) 
 allocate (matrix_1h2p(N_det_ref,N_det_ref,N_states))
 integer :: i,j,istate
 double precision :: accu
 double precision :: accu_bis(N_states)
 do i = 1, N_det_ref
  do j = 1, N_det_ref
   do istate = 1, N_states
    matrix_1h2p(i,j,istate) = 0.d0
   enddo
  enddo
 enddo
 double precision :: pt2(N_states)
 call contrib_2h1p_dm_based(pt2)
 print*, 'pt2 = ',pt2
 matrix_1h2p = 0.d0
 call H_apply_mrpt_2h1p(matrix_1h2p,N_det_ref)
 do istate = 1, N_states
 accu = 0.d0
 do i = 1, N_det_ref
  do j = 1, N_det_ref 
   accu += matrix_1h2p(i,j,istate) * psi_ref_coef(i,istate) * psi_ref_coef(j,istate)
  enddo
 enddo
 print*, 'second order ', accu
 print*, 'H apply - density matrix = ',accu - pt2(istate)
 enddo
 
 integer :: a,b,ispin
 accu_bis = 0.d0
 do istate = 1, N_states
  accu_bis(istate) += effective_fock_operator_pure_diag_2h1p(istate) 
  do ispin = 1, 2
  do i = 1, n_act_orb
   a = list_act(i)
   do j = 1, n_act_orb
    b = list_act(j)
    accu_bis(istate) += effective_fock_operator_2h1p(i,j,ispin,istate) * one_body_dm_mo_spin_index(b,a,istate,ispin)
   enddo
  enddo
  enddo
  print*, 'accu_bis = ',accu_bis(istate)
  print*, 'H apply - density matrix = ',accu_bis(istate) - pt2(istate)
 enddo
end

subroutine routine_1h2p
 implicit none
 double precision,allocatable :: matrix_1h2p(:,:,:) 
 allocate (matrix_1h2p(N_det_ref,N_det_ref,N_states))
 integer :: i,j,istate
 double precision :: accu
 double precision :: accu_bis(N_states)
 do i = 1, N_det_ref
  do j = 1, N_det_ref
   do istate = 1, N_states
    matrix_1h2p(i,j,istate) = 0.d0
   enddo
  enddo
 enddo

 call contrib_1h2p_dm_based(accu_bis)
 print*, 'accu_bis',accu_bis
 print*, '1h2p    ',contribution_1h2p
 matrix_1h2p = 0.d0
 call H_apply_mrpt_1h2p(matrix_1h2p,N_det_ref)
 do istate = 1, N_states
 accu = 0.d0
 do i = 1, N_det_ref
  do j = 1, N_det_ref 
   accu += matrix_1h2p(i,j,istate) * psi_ref_coef(i,istate) * psi_ref_coef(j,istate)
  enddo
 enddo
 print*, 'second order ', accu
 print*, 'H apply - density matrix = ',accu - accu_bis(istate)
 enddo
 accu = 0.d0
 integer :: ispin,a,b
 do istate = 1, N_states
  do ispin = 1, 2
   do i = 1, n_act_orb
    a = list_act(i)
    do j = 1, n_act_orb 
     b = list_act(j)
     accu += effective_fock_operator_1h2p(i,j,ispin,istate) * one_body_dm_mo_spin_index(b,a,istate,ispin)
    enddo
   enddo
  enddo
  print*, 'accu = ',accu
 print*, 'H apply - density matrix = ',accu - accu_bis(istate)
 enddo


 deallocate (matrix_1h2p)
end
