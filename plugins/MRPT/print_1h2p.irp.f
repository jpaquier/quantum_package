program print_1h2p
 implicit none
 read_wf = .True.
 touch read_wf
 call routine
end

subroutine routine
 implicit none
 double precision,allocatable :: matrix_1h2p(:,:,:) 
 double precision,allocatable :: matrix_2h1p(:,:,:) 
 allocate (matrix_1h2p(N_det,N_det,N_states))
 allocate (matrix_2h1p(N_det,N_det,N_states))
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
 provide effective_fock_operator_1h2p
 call give_2h1p_contrib(matrix_2h1p)
 call contrib_2h1p_dm_based(accu_bis)
 print*, 'accu_bis     ', accu_bis(1)
 accu = 0.d0
 do i = 1, N_det
  do j = 1, N_det 
   accu += matrix_2h1p(i,j,1) * psi_coef(i,1) * psi_coef(j,1)
  enddo
 enddo
 print*, 'second order ', accu
 endif


 deallocate (matrix_2h1p)
end
