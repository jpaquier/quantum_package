
use bitmasks
subroutine contrib_1h2p_dm_based(accu)
 double precision, intent(out) :: accu(N_states)
 implicit none
 integer :: i_i,i_r,i_v,i_a,i_b
 integer :: i,r,v,a,b
 integer :: ispin,jspin
 integer :: istate
 double precision :: active_int(n_act_orb,2)
 double precision :: delta_e(n_act_orb,2,N_states)
 double precision :: get_mo_bielec_integral
 accu = 0.d0
 do i_i = 1, n_inact_orb
  i = list_inact(i_i)
  do i_r = 1, n_virt_orb
   r = list_virt(i_r)
   do i_v = 1, n_virt_orb
    v = list_virt(i_v)
    do i_a = 1, n_act_orb
     a = list_act(i_a)
     active_int(i_a,1) = get_mo_bielec_integral(i,a,r,v,mo_integrals_map) ! direct
     active_int(i_a,2) = get_mo_bielec_integral(i,a,v,r,mo_integrals_map) ! exchange
     do istate = 1, N_states
      do jspin=1, 2
       delta_e(i_a,jspin,istate) = one_anhil(i_a,jspin,istate)                        &
                                 - fock_virt_total_spin_trace(r,istate)               & 
                                 - fock_virt_total_spin_trace(v,istate)               & 
                                 + fock_core_inactive_total_spin_trace(i,istate)        
       delta_e(i_a,jspin,istate) = 1.d0/delta_e(i_a,jspin,istate)  
      enddo
     enddo
    enddo
    do i_a = 1, n_act_orb
     a = list_act(i_a)
     do i_b = 1, n_act_orb
      b = list_act(i_b)
      do ispin = 1, 2 ! spin of (i --> r)
       do jspin = 1, 2 ! spin of (a --> v)
        if(ispin == jspin .and. r.le.v)cycle ! condition not to double count 
        do istate = 1, N_states
         if(ispin == jspin)then
          accu(istate) += (active_int(i_a,1) - active_int(i_a,2)) * one_body_dm_mo_spin_index(a,b,istate,ispin)   &
                                                          * (active_int(i_b,1) - active_int(i_b,2)) & 
                                                          * delta_e(i_a,jspin,istate)
         else 
          accu(istate) += active_int(i_a,1)  * one_body_dm_mo_spin_index(a,b,istate,ispin) * delta_e(i_a,ispin,istate) & 
                        * active_int(i_b,1) 
         endif
        enddo
       enddo
      enddo
     enddo
    enddo
   enddo
  enddo
 enddo


end


subroutine matrix_1h2p_dm_based(matrix)
 double precision, intent(inout) :: matrix(N_det,N_det,N_states)
 implicit none
 integer :: i_i,i_r,i_v,i_a,i_b
 integer :: i,r,v,a,b
 integer :: ispin,jspin
 integer :: istate
 double precision :: active_int(n_act_orb,2)
 double precision :: delta_e(n_act_orb,2,N_states)
 double precision :: get_mo_bielec_integral
 double precision :: accu(N_states)
 accu = 0.d0
 print*, 'in the matrix  ....'
 do i_i = 1, n_inact_orb
  i = list_inact(i_i)
  do i_r = 1, n_virt_orb
   r = list_virt(i_r)
   do i_v = 1, n_virt_orb
    v = list_virt(i_v)
    do i_a = 1, n_act_orb
     a = list_act(i_a)
     active_int(i_a,1) = get_mo_bielec_integral(i,a,r,v,mo_integrals_map) ! direct
     active_int(i_a,2) = get_mo_bielec_integral(i,a,v,r,mo_integrals_map) ! exchange
     do istate = 1, N_states
      do jspin=1, 2
       delta_e(i_a,jspin,istate) = one_anhil(i_a,jspin,istate)                        &
                                 - fock_virt_total_spin_trace(r,istate)               & 
                                 - fock_virt_total_spin_trace(v,istate)               & 
                                 + fock_core_inactive_total_spin_trace(i,istate)        
       delta_e(i_a,jspin,istate) = 1.d0/delta_e(i_a,jspin,istate)  
      enddo
     enddo
    enddo
    do i_a = 1, n_act_orb
     a = list_act(i_a)
     do i_b = 1, n_act_orb
      b = list_act(i_b)
      do ispin = 1, 2 ! spin of (i --> r)
       do jspin = 1, 2 ! spin of (a --> v)
        if(ispin == jspin .and. r.le.v)cycle ! condition not to double count 
        do istate = 1, N_states
         if(ispin == jspin)then
          accu(istate) += (active_int(i_a,1) - active_int(i_a,2)) * one_body_dm_mo_spin_index(a,b,istate,ispin)   &
                                                          * (active_int(i_b,1) - active_int(i_b,2)) & 
                                                          * delta_e(i_a,jspin,istate)
         else 
          accu(istate) += active_int(i_a,1)  * one_body_dm_mo_spin_index(a,b,istate,ispin) * delta_e(i_a,ispin,istate) & 
                        * active_int(i_b,1) 
         endif
        enddo
       enddo
      enddo
     enddo
    enddo
   enddo
  enddo
 enddo
 print*, 'accu1',accu

 double precision :: accu_2(n_act_orb,n_act_orb,2,N_states)
 accu_2 = 0.d0
 do i_a = 1, n_act_orb
  a = list_act(i_a)
    do jspin=1, 2
     do istate = 1, N_states
      do i_i = 1, n_inact_orb
       i = list_inact(i_i)
       do i_r = 1, n_virt_orb
        r = list_virt(i_r)
        do i_v = 1, n_virt_orb
         v = list_virt(i_v)
         do i_b = 1, n_act_orb 
          b = list_act(i_b)
          active_int(i_b,1) = get_mo_bielec_integral(i,b,r,v,mo_integrals_map) ! direct
          active_int(i_b,2) = get_mo_bielec_integral(i,b,v,r,mo_integrals_map) ! exchange
         enddo
          delta_e(i_a,jspin,istate) = one_anhil(i_a,jspin,istate)                        &
                                    - fock_virt_total_spin_trace(r,istate)               & 
                                    - fock_virt_total_spin_trace(v,istate)               & 
                                    + fock_core_inactive_total_spin_trace(i,istate)        
          delta_e(i_a,jspin,istate) = 1.d0/delta_e(i_a,jspin,istate)  
         do i_b = 1, n_act_orb 
          b = list_act(i_b)
          do ispin = 1, 2 ! spin of (i --> r)
           if(ispin == jspin .and. r.le.v)cycle ! condition not to double count 
            if(ispin == jspin)then
             accu_2(i_a,i_b,jspin,istate) += (active_int(i_a,1) - active_int(i_a,2)) * (active_int(i_b,1) - active_int(i_b,2)) * delta_e(i_a,jspin,istate)
                                                             
            else 
             accu_2(i_a,i_b,jspin,istate) += active_int(i_a,1)   * delta_e(i_a,jspin,istate) * active_int(i_b,1)  
            endif
          enddo
         enddo
        enddo
       enddo
      enddo
    enddo
   enddo
 enddo
 accu = 0.d0
 do istate = 1, N_states
  do jspin = 1, 2 
   do i_a = 1,n_act_orb
    a = list_act(i_a)
    do i_b = 1, n_act_orb 
     b = list_act(i_b)
     accu(istate) += accu_2(i_a,i_b,jspin,istate) * one_body_dm_mo_spin_index(a,b,istate,jspin)
    enddo
   enddo
  enddo
 enddo
 print*, 'accu2',accu

end



 BEGIN_PROVIDER [double precision, effective_fock_operator_1h2p, (n_act_orb,n_act_orb,2,N_states)]
&BEGIN_PROVIDER [double precision, contribution_1h2p, (N_states)]
 implicit none
 integer :: i_i,i_r,i_v,i_a,i_b
 integer :: i,r,v,a,b
 integer :: ispin,jspin
 integer :: istate
 double precision :: active_int(n_act_orb,2)
 double precision :: delta_e(n_act_orb,2,N_states)
 double precision :: get_mo_bielec_integral

 effective_fock_operator_1h2p = 0.d0
 do i_a = 1, n_act_orb
  a = list_act(i_a)
    do jspin=1, 2
     do istate = 1, N_states
      do i_i = 1, n_inact_orb
       i = list_inact(i_i)
       do i_r = 1, n_virt_orb
        r = list_virt(i_r)
        do i_v = 1, n_virt_orb
         v = list_virt(i_v)
         do i_b = 1, n_act_orb 
          b = list_act(i_b)
          active_int(i_b,1) = get_mo_bielec_integral(i,b,r,v,mo_integrals_map) ! direct
          active_int(i_b,2) = get_mo_bielec_integral(i,b,v,r,mo_integrals_map) ! exchange
         enddo
          delta_e(i_a,jspin,istate) = one_anhil(i_a,jspin,istate)                        &
                                    - fock_virt_total_spin_trace(r,istate)               & 
                                    - fock_virt_total_spin_trace(v,istate)               & 
                                    + fock_core_inactive_total_spin_trace(i,istate)        
          delta_e(i_a,jspin,istate) = 1.d0/delta_e(i_a,jspin,istate)  
         do i_b = 1, n_act_orb 
          b = list_act(i_b)
          do ispin = 1, 2 ! spin of (i --> r)
           if(ispin == jspin .and. r.le.v)cycle ! condition not to double count 
            if(ispin == jspin)then
             effective_fock_operator_1h2p(i_a,i_b,jspin,istate) += (active_int(i_a,1) - active_int(i_a,2)) * (active_int(i_b,1) - active_int(i_b,2)) * delta_e(i_a,jspin,istate)
                                                             
            else 
             effective_fock_operator_1h2p(i_a,i_b,jspin,istate) += active_int(i_a,1)   * delta_e(i_a,jspin,istate) * active_int(i_b,1)  
            endif
          enddo
         enddo
        enddo
       enddo
      enddo
    enddo
   enddo
 enddo
 contribution_1h2p = 0.d0
 do istate = 1, N_states
  do jspin = 1, 2 
   do i_a = 1,n_act_orb
    a = list_act(i_a)
    do i_b = 1, n_act_orb 
     b = list_act(i_b)
     contribution_1h2p(istate) += effective_fock_operator_1h2p(i_a,i_b,jspin,istate) * one_body_dm_mo_spin_index(a,b,istate,jspin)
    enddo
   enddo
  enddo
 enddo
 print*, 'contribution_1h2p2',contribution_1h2p

END_PROVIDER 





 BEGIN_PROVIDER [integer, list_det_single_exc, (n_act_orb,n_act_orb,2,2)]
&BEGIN_PROVIDER [logical, single_exc_ok, (n_act_orb,n_act_orb,2)]
&BEGIN_PROVIDER [double precision, phase_single_exc, (n_act_orb,n_act_orb,2)]
 implicit none
 BEGIN_DOC
! list_det_single_exc(iorb_a,iorb_b,ispin,1) == index of a determinant who is connected to 
!                                         list_det_single_exc(iorb_a,iorb_b,ispin,2) by a single excitation a^{\dagger}_(iorb_b,ipsin) a_{(iorb_a,ipsin)}
 END_DOC

 integer :: idet,jdet,k,l,hole,part
 integer :: exc(0:2,2,2)
 double precision :: phase
 integer           :: degree(N_det)
 integer           :: idx(0:N_det)
 single_exc_ok = .True.
 do idet = 1, N_det
   call get_excitation_degree_vector_mono(psi_det,psi_det(1,1,idet),degree,N_int,N_det,idx)
   do jdet = 1, idx(0)
     if(idx(jdet).ne.idet)then
      call get_mono_excitation(psi_det(1,1,idet),psi_det(1,1,idx(jdet)),exc,phase,N_int)
      if (exc(0,1,1) == 1) then
        ! Mono alpha
        hole = list_act_reverse(exc(1,1,1))   !!!  a_a
        part = list_act_reverse(exc(1,2,1))   !!!  a^{\dagger}_{b}
        list_det_single_exc(hole,part,1,1) = idet
        list_det_single_exc(hole,part,1,2) = idx(jdet)
        single_exc_ok(hole,part,1) = .True.
        phase_single_exc(hole,part,1) = phase
      else
        ! Mono beta
        hole = list_act_reverse(exc(1,1,2))   !!!  a_a
        part = list_act_reverse(exc(1,2,2))   !!!  a^{\dagger}_{b}
        list_det_single_exc(hole,part,2,1) = idet
        list_det_single_exc(hole,part,2,2) = idx(jdet)
        single_exc_ok(hole,part,2) = .False.
        phase_single_exc(hole,part,2) = phase
      endif
     endif
   enddo
 enddo

END_PROVIDER 


subroutine contrib_2h1p_dm_based(accu)
 implicit none
 integer :: i_i,i_j,i_v,i_a,i_b
 integer :: i,j,v,a,b
 integer :: ispin,jspin
 integer :: istate
 double precision, intent(out) :: accu(N_states)
 double precision :: active_int(n_act_orb,2)
 double precision :: delta_e(n_act_orb,2,N_states)
 double precision :: get_mo_bielec_integral
 accu = 0.d0
 do i_i = 1, n_inact_orb
  i = list_inact(i_i)
  do i_j = 1, n_inact_orb
   j = list_inact(i_j)
   do i_v = 1, n_virt_orb
    v = list_virt(i_v)
    do i_a = 1, n_act_orb
     a = list_act(i_a)
     active_int(i_a,1) = get_mo_bielec_integral(i,j,v,a,mo_integrals_map) ! direct
     active_int(i_a,2) = get_mo_bielec_integral(i,j,a,v,mo_integrals_map) ! exchange
     do istate = 1, N_states
      do jspin=1, 2
       delta_e(i_a,jspin,istate) = one_creat(i_a,jspin,istate)  - fock_virt_total_spin_trace(v,istate)  &
                                 + fock_core_inactive_total_spin_trace(i,istate)       &
                                 + fock_core_inactive_total_spin_trace(j,istate)        
       delta_e(i_a,jspin,istate) = 1.d0/delta_e(i_a,jspin,istate)  
      enddo
     enddo
    enddo
    do i_a = 1, n_act_orb
     a = list_act(i_a)
     do i_b = 1, n_act_orb
      b = list_act(i_b)
      do ispin = 1, 2 ! spin of (i --> v)
       do jspin = 1, 2 ! spin of (j --> a)
        if(ispin == jspin .and. i.le.j)cycle ! condition not to double count 
!!!!!!!!!!!!!!!
!!!!! TEST FOR THE SAME SPIN
!!!!!!!!!!!!!!!
!       if(ispin == jspin)cycle ! condition not to double count 
        do istate = 1, N_states
         if(ispin == jspin)then
          accu(istate) += (active_int(i_a,1) - active_int(i_a,2)) * one_body_dm_dagger_mo_spin_index(a,b,istate,ispin)   &
                                                          * (active_int(i_b,1) - active_int(i_b,2)) & 
                                                          * delta_e(i_a,jspin,istate)
         else 
          accu(istate) += active_int(i_a,1)  * one_body_dm_dagger_mo_spin_index(a,b,istate,ispin) * delta_e(i_a,ispin,istate) & 
                        * active_int(i_b,1) 
         endif
        enddo
       enddo
      enddo
     enddo
    enddo
   enddo
  enddo
 enddo


end


!subroutine contrib_2p_dm_based(accu)
!implicit none
!integer :: i_r,i_v,i_a,i_b,i_c,i_d
!integer :: r,v,a,b,c,d
!integer :: ispin,jspin
!integer :: istate
!double precision, intent(out) :: accu(N_states)
!double precision :: active_int(n_act_orb,n_act_orb,2)
!double precision :: delta_e(n_act_orb,n_act_orb,2,2,N_states)
!double precision :: get_mo_bielec_integral
!accu = 0.d0
!do i_r = 1, n_virt_orb
! r = list_virt(i_r)
! do i_v = 1, n_virt_orb
!   v = list_virt(i_v)
!   do i_a = 1, n_act_orb
!    a = list_act(i_a)
!    do i_b = 1, n_act_orb
!     b = list_act(i_b)
!     active_int(i_a,i_b,1) = get_mo_bielec_integral(a,b,r,v,mo_integrals_map) ! direct
!     active_int(i_a,i_b,2) = get_mo_bielec_integral(a,b,v,r,mo_integrals_map) ! direct
!     do istate = 1, N_states
!      do jspin=1, 2 ! spin of i_a
!       do ispin = 1, 2 ! spin of i_b
!        delta_e(i_a,i_b,jspin,ispin,istate) = two_anhil(i_a,i_b,jspin,ispin,istate)    &
!                                  - fock_virt_total_spin_trace(r,istate)               & 
!                                  - fock_virt_total_spin_trace(v,istate)                 
!        delta_e(i_a,i_b,jspin,ispin,istate) = 1.d0/delta_e(i_a,i_b,jspin,ispin,istate)  
!       enddo
!      enddo
!     enddo
!    enddo
!   enddo
!   ! diagonal terms 
!   do i_a = 1, n_act_orb
!    a = list_act(i_a)
!    do i_b = 1, n_act_orb
!     b = list_act(i_b)
!     do ispin = 1, 2 ! spin of (a --> r)
!      do jspin = 1, 2 ! spin of (b --> v)
!       if(ispin == jspin .and. r.le.v)cycle ! condition not to double count 
!       if(ispin == jspin .and. a.le.b)cycle ! condition not to double count 
!       do istate = 1, N_states
!        if(ispin == jspin)then
!         double precision :: contrib_spin
!         if(ispin == 1)then
!          contrib_spin = two_body_dm_aa_diag_act(i_a,i_b)
!         else
!          contrib_spin = two_body_dm_bb_diag_act(i_a,i_b)
!         endif
!         accu(istate) += (active_int(i_a,i_b,1) - active_int(i_a,i_b,2)) * contrib_spin  &
!                       * (active_int(i_a,i_b,1) - active_int(i_a,i_b,2)) & 
!                       * delta_e(i_a,i_b,ispin,jspin,istate)
!        else 
!         accu(istate) += 0.5d0 * active_int(i_a,i_b,1)  * two_body_dm_ab_diag_act(i_a,i_b) * delta_e(i_a,i_b,ispin,jspin,istate) & 
!                               * active_int(i_a,i_b,1) 
!        endif
!       enddo
!      enddo
!     enddo
!    enddo
!   enddo
!  enddo
! enddo


!end

