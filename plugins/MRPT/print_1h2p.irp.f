program print_1h2p
 implicit none
 read_wf = .True.
 touch read_wf
 call routine
end

subroutine routine
 implicit none
 double precision,allocatable :: matrix_1h2p(:,:,:) 
 allocate (matrix_1h2p(N_det,N_det,N_states))
 integer :: i,j,istate
 double precision :: accu
 double precision :: accu_bis(N_states)
 do i = 1, N_det
  do j = 1, N_det
   do istate = 1, N_states
    matrix_1h2p(i,j,istate) = 0.d0
   enddo
  enddo
 enddo
 if(.True.)then
 call contrib_1h2p_dm_based(accu_bis)
 print*, 'accu_bis     ', accu_bis(1)
 provide effective_fock_operator_1h2p
 call matrix_1h2p_dm_based(matrix_1h2p)
 accu = 0.d0
 do i = 1, N_det
  do j = 1, N_det 
   accu += matrix_1h2p(i,j,1) * psi_coef(i,1) * psi_coef(j,1)
  enddo
 enddo
 print*, 'second order ', accu
 matrix_1h2p = 0.d0
 call give_1h2p_contrib(matrix_1h2p)
 accu = 0.d0
 do i = 1, N_det
  do j = 1, N_det 
   accu += matrix_1h2p(i,j,1) * psi_coef(i,1) * psi_coef(j,1)
  enddo
 enddo
 print*, 'second order ', accu
!matrix_1h2p = 0.d0
!call H_apply_mrpt_1h2p(matrix_1h2p,N_det)
!accu = 0.d0
!do i = 1, N_det
! do j = 1, N_det 
!  accu += matrix_1h2p(i,j,1) * psi_coef(i,1) * psi_coef(j,1)
! enddo
!enddo
 endif


 deallocate (matrix_1h2p)
end
