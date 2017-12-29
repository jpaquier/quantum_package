program print_1h2p
 implicit none
 read_wf = .True.
 touch read_wf
 call routine 
end

subroutine routine3
 implicit none
 integer :: i
 provide fock_core_inactive
 provide one_anhil
!do i = 1, mo_tot_num
! print*, 'fock_core_inactive',fock_core_inactive(i)
!enddo

end

subroutine routine4
 implicit none
 double precision,allocatable :: matrix_1h2p(:,:,:) 
 allocate (matrix_1h2p(N_det,N_det,N_states))
 call give_1h2p_contrib(matrix_1h2p)


end

subroutine routine
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

subroutine routine2
  use bitmasks
 implicit none
 integer :: i,j,m
 integer :: idet,jdet,hole,part
 integer :: i_state,ispin
 double precision :: accu(N_states)
 double precision, allocatable :: delta_ij_tmp(:,:,:)
 integer           :: degree(N_det)
 integer           :: idx(0:N_det)
 integer :: exc(0:2,2,2)
 double precision :: phase
 integer :: occ(N_int*bit_kind_size,2)
 integer  :: n_elec_tmp(2)
 integer :: iorb_a
 allocate (delta_ij_tmp(N_det,N_det,N_states))

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 ! 1h2p   
 delta_ij_tmp = 0.d0
!call give_1h2p_contrib(delta_ij_tmp)
!call H_apply_mrpt_1h2p(delta_ij_tmp,N_det)
 accu = 0.d0
 do i_state = 1, N_states
  do idet = 1, N_det
   !!! Diagonal element 
   call bitstring_to_list_ab(psi_active(1,1,idet), occ, n_elec_tmp, N_int)
   do ispin = 1, 2
    do i = 1, n_elec_tmp(ispin)
     iorb_a = list_act_reverse(occ(i,ispin))
     delta_ij_tmp(idet,idet,i_state) += effective_fock_operator_1h2p(iorb_a,iorb_a,ispin,i_state)
    enddo
   enddo
   accu(i_state) += delta_ij_tmp(idet,idet,i_State) * psi_coef(idet,i_state)**2
   !!! Extra diagonal elements 
   call get_excitation_degree_vector_mono(psi_det,psi_det(1,1,idet),degree,N_int,N_det,idx)
   do jdet = 1, idx(0)
    if(idx(jdet)==idet)cycle
    call get_mono_excitation(psi_det(1,1,idet),psi_det(1,1,idx(jdet)),exc,phase,N_int)
    if (exc(0,1,1) == 1) then
       ! Mono alpha
       hole = list_act_reverse(exc(1,1,1))   !!!  a_a
       part = list_act_reverse(exc(1,2,1))   !!!  a^{\dagger}_{b}
       ispin =  1
    else
       ! Mono beta
       hole = list_act_reverse(exc(1,1,2))   !!!  a_a
       part = list_act_reverse(exc(1,2,2))   !!!  a^{\dagger}_{b}
       ispin =  2
    endif
    delta_ij_tmp(idet,idx(jdet),i_state) += effective_fock_operator_1h2p(hole,part,ispin,i_state) * phase
    accu(i_state) += effective_fock_operator_1h2p(hole,part,ispin,i_state) * phase * psi_coef(idet,i_state) * psi_coef(idx(jdet),i_state)
   enddo
  enddo
 enddo
 print*, '1h2p              = ',accu
 print*, 'contribution_1h2p = ',contribution_1h2p
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
stop



end
