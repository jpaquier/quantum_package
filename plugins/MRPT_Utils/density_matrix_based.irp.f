


 BEGIN_PROVIDER [double precision, effective_fock_operator_2h1p, (n_act_orb,n_act_orb,2,N_states)]
&BEGIN_PROVIDER [double precision, contribution_2h1p, (N_states)]
 implicit none
 integer :: i_i,i_v,i_j,i_a,i_b
 integer :: i,v,j,a,b
 integer :: ispin,jspin
 integer :: istate
 double precision :: active_int(n_act_orb,2)
 double precision :: delta_e(n_act_orb,2,N_states)
 double precision :: get_mo_bielec_integral

 effective_fock_operator_2h1p = 0.d0
 do i_i = 1, n_inact_orb
  i = list_inact(i_i)
  do i_j = 1, n_inact_orb
   j = list_inact(i_j)
   do i_v = 1, n_virt_orb
    v = list_virt(i_v)
     do istate = 1, N_states
      do i_b = 1, n_act_orb 
       b = list_act(i_b)
       active_int(i_b,1) = get_mo_bielec_integral(i,j,v,b,mo_integrals_map) ! direct
       active_int(i_b,2) = get_mo_bielec_integral(i,j,b,v,mo_integrals_map) ! exchange
       do jspin=1, 2
        delta_e(i_b,jspin,istate) = one_creat(i_b,jspin,istate)                        &
                                  - fock_virt_total_spin_trace(v,istate)               & 
                                  + fock_core_inactive_total_spin_trace(j,istate)               & 
                                  + fock_core_inactive_total_spin_trace(i,istate)        
        delta_e(i_b,jspin,istate) = 1.d0/delta_e(i_b,jspin,istate)  
       enddo
      enddo
     enddo
     do i_a = 1, n_act_orb
      a = list_act(i_a)
      do i_b = 1, n_act_orb 
       b = list_act(i_b)
       do jspin=1, 2
        do istate = 1, N_states
          do ispin = 1, 2 ! spin of (i --> r)
           if(ispin == jspin .and. i.le.j)cycle ! condition not to double count 
            if(ispin == jspin)then
             effective_fock_operator_2h1p(i_a,i_b,jspin,istate) += (active_int(i_a,1) - active_int(i_a,2)) * (active_int(i_b,1) - active_int(i_b,2)) * delta_e(i_a,jspin,istate)
                                                             
            else 
             effective_fock_operator_2h1p(i_a,i_b,jspin,istate) += active_int(i_a,1)   * delta_e(i_a,jspin,istate) * active_int(i_b,1)  
            endif
          enddo
         enddo
        enddo
       enddo
      enddo
    enddo
   enddo
 enddo


 contribution_2h1p = 0.d0
 do istate = 1, N_states
  do jspin = 1, 2 
   do i_a = 1,n_act_orb
    a = list_act(i_a)
    do i_b = 1, n_act_orb 
     b = list_act(i_b)
     contribution_2h1p(istate) += effective_fock_operator_2h1p(i_a,i_b,jspin,istate) * one_body_dm_dagger_mo_spin_index(a,b,istate,jspin)
    enddo
   enddo
  enddo
 enddo
 print*, 'contribution_2h1p',contribution_2h1p

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




subroutine contrib_1h2p_a_b_dm_based(pt2)
BEGIN_DOC
 ! gives the contribution of the 1h2p excitations of opposite spin
END_DOC
 implicit none
 double precision, intent(out) :: pt2(N_states)
 
 integer :: i_i,i_t,i_v,i_a,i_b
 integer :: i,t,v,a,b
 integer :: istate
 double precision :: delta_e(n_act_orb,2,N_states)
 double precision :: get_mo_bielec_integral
 double precision :: active_int(n_act_orb,2)
 double precision :: accu_direct
 integer :: ispin
 pt2 = 0.d0
 do i_i = 1, n_inact_orb
  i = list_inact(i_i)
  do i_t = 1, n_virt_orb
   t = list_virt(i_t)
   do i_v = 1, n_virt_orb
    v = list_virt(i_v)
    do i_a = 1, n_act_orb
     a = list_act(i_a)
     active_int(i_a,1) = get_mo_bielec_integral(i,a,t,v,mo_integrals_map) ! direct
     active_int(i_a,2) = get_mo_bielec_integral(i,a,v,t,mo_integrals_map) ! exchange
     do istate = 1, N_states
      do ispin = 1, 2
       delta_e(i_a,ispin,istate) = one_anhil(i_a,ispin,istate)                        &
                                 - fock_virt_total_spin_trace(t,istate)       & 
                                 - fock_virt_total_spin_trace(v,istate)       & 
                                 + fock_core_inactive_total_spin_trace(i,istate)        
       delta_e(i_a,ispin,istate) = 1.d0/delta_e(i_a,ispin,istate)  
      enddo
     enddo
    enddo
    do ispin = 1, 2
     do istate = 1, N_states 
      do i_a = 1, n_act_orb
       a = list_act(i_a)
       do i_b = 1, n_act_orb
        b = list_act(i_b)
        accu_direct = active_int(i_a,1) * active_int(i_b,1) 
        pt2(istate) += accu_direct  * delta_e(i_a,ispin,istate) * one_body_dm_mo_spin_index(b,a,istate,ispin)
       enddo
      enddo
     enddo
    enddo
   enddo
  enddo
 enddo


end

subroutine contrib_1h2p_parallel_spin_dm_based(pt2)
BEGIN_DOC
 ! gives the contribution of the 1h2p excitations of parallel spin
END_DOC
 implicit none
 double precision, intent(out) :: pt2(N_states)
 
 integer :: i_i,i_t,i_v,i_a,i_b
 integer :: i,t,v,a,b
 integer :: istate
 double precision :: delta_e(n_act_orb,2,N_states)
 double precision :: get_mo_bielec_integral
 double precision :: active_int(n_act_orb,2)
 double precision :: accu_direct
 integer :: ispin
 pt2 = 0.d0
 do i_i = 1, n_inact_orb
  i = list_inact(i_i)
  do i_v = 1, n_virt_orb
   v = list_virt(i_v)
   do i_t = i_v+1, n_virt_orb
    t = list_virt(i_t)
    do i_a = 1, n_act_orb
     a = list_act(i_a)
     active_int(i_a,1) = get_mo_bielec_integral(i,a,t,v,mo_integrals_map) ! direct
     active_int(i_a,2) = get_mo_bielec_integral(i,a,v,t,mo_integrals_map) ! exchange
     do istate = 1, N_states
      do ispin = 1, 2
       delta_e(i_a,ispin,istate) = one_anhil(i_a,ispin,istate)                        &
                                 - fock_virt_total_spin_trace(t,istate)       & 
                                 - fock_virt_total_spin_trace(v,istate)       & 
                                 + fock_core_inactive_total_spin_trace(i,istate)        
       delta_e(i_a,ispin,istate) = 1.d0/delta_e(i_a,ispin,istate)  
      enddo
     enddo
    enddo
    do ispin = 1, 2
     do istate = 1, N_states 
      do i_a = 1, n_act_orb
       a = list_act(i_a)
       do i_b = 1, n_act_orb
        b = list_act(i_b)
        accu_direct = (active_int(i_a,1) - active_int(i_a,2)) * (active_int(i_b,1) - active_int(i_b,2))
        pt2(istate) += accu_direct  * delta_e(i_a,ispin,istate) * one_body_dm_mo_spin_index(b,a,istate,ispin)
!       if(one_body_dm_mo_spin_index(a,b,istate,ispin).gt.0.d0)then
!       print*, active_int(i_a,1),1.d0/delta_e(i_a,1,istate),one_body_dm_mo_spin_index(a,b,istate,ispin)
!       endif
       enddo
      enddo
     enddo
    enddo
   enddo
  enddo
 enddo


end

subroutine contrib_1h2p_dm_based(pt2)
 implicit none
 double precision, intent(out) :: pt2(N_states)
 
 integer :: i_i,i_t,i_v,i_a,i_b
 integer :: i,t,v,a,b
 integer :: istate
 double precision :: delta_e(n_act_orb,2,N_states)
 double precision :: get_mo_bielec_integral
 double precision :: active_int(n_act_orb,2)
 double precision :: accu_direct
 integer :: ispin
 pt2 = 0.d0
 do i_i = 1, n_inact_orb
  i = list_inact(i_i)
  do i_v = 1, n_virt_orb
   v = list_virt(i_v)
   do i_t = 1,i_v
    t = list_virt(i_t)
    do i_a = 1, n_act_orb
     a = list_act(i_a)
     active_int(i_a,1) = get_mo_bielec_integral(i,a,t,v,mo_integrals_map) ! direct
     do istate = 1, N_states
      do ispin = 1, 2
       delta_e(i_a,ispin,istate) = one_anhil(i_a,ispin,istate)                        &
                                 - fock_virt_total_spin_trace(t,istate)       & 
                                 - fock_virt_total_spin_trace(v,istate)       & 
                                 + fock_core_inactive_total_spin_trace(i,istate)        
       delta_e(i_a,ispin,istate) = 1.d0/delta_e(i_a,ispin,istate)  
      enddo
     enddo
    enddo
    do ispin = 1, 2
     do istate = 1, N_states 
      do i_a = 1, n_act_orb
       a = list_act(i_a)
       do i_b = 1, n_act_orb
        b = list_act(i_b)
        accu_direct = active_int(i_a,1) * active_int(i_b,1)
        pt2(istate) += accu_direct  * delta_e(i_a,ispin,istate) * one_body_dm_mo_spin_index(b,a,istate,ispin)
       enddo
      enddo
     enddo
    enddo
   enddo
   do i_t = i_v+1, n_virt_orb
    t = list_virt(i_t)
    do i_a = 1, n_act_orb
     a = list_act(i_a)
     active_int(i_a,1) = get_mo_bielec_integral(i,a,t,v,mo_integrals_map) ! direct
     active_int(i_a,2) = get_mo_bielec_integral(i,a,v,t,mo_integrals_map) ! exchange
     do istate = 1, N_states
      do ispin = 1, 2
       delta_e(i_a,ispin,istate) = one_anhil(i_a,ispin,istate)                        &
                                 - fock_virt_total_spin_trace(t,istate)       & 
                                 - fock_virt_total_spin_trace(v,istate)       & 
                                 + fock_core_inactive_total_spin_trace(i,istate)        
       delta_e(i_a,ispin,istate) = 1.d0/delta_e(i_a,ispin,istate)  
      enddo
     enddo
    enddo
    do ispin = 1, 2
     do istate = 1, N_states 
      do i_a = 1, n_act_orb
       a = list_act(i_a)
       do i_b = 1, n_act_orb
        b = list_act(i_b)
        accu_direct = (active_int(i_a,1) - active_int(i_a,2)) * (active_int(i_b,1) - active_int(i_b,2)) + active_int(i_a,1) * active_int(i_b,1)
        pt2(istate) += accu_direct  * delta_e(i_a,ispin,istate) * one_body_dm_mo_spin_index(b,a,istate,ispin)
       enddo
      enddo
     enddo
    enddo
   enddo
  enddo
 enddo


end

 BEGIN_PROVIDER [double precision, effective_fock_operator_1h2p, (n_act_orb,n_act_orb,2,N_states)]
&BEGIN_PROVIDER [double precision, contribution_1h2p, (N_states)]

 implicit none
 
 integer :: i_i,i_t,i_v,i_a,i_b
 integer :: i,t,v,a,b
 integer :: istate
 double precision :: delta_e(n_act_orb,2,N_states)
 double precision :: get_mo_bielec_integral
 double precision :: active_int(n_act_orb,2)
 double precision :: accu_direct
 integer :: ispin
 contribution_1h2p = 0.d0
 effective_fock_operator_1h2p = 0.d0
 do i_i = 1, n_inact_orb
  i = list_inact(i_i)
  do i_v = 1, n_virt_orb
   v = list_virt(i_v)
   do i_t = 1,i_v
    t = list_virt(i_t)
    do i_a = 1, n_act_orb
     a = list_act(i_a)
     active_int(i_a,1) = get_mo_bielec_integral(i,a,t,v,mo_integrals_map) ! direct
     do istate = 1, N_states
      do ispin = 1, 2
       delta_e(i_a,ispin,istate) = one_anhil(i_a,ispin,istate)                        &
                                 - fock_virt_total_spin_trace(t,istate)       & 
                                 - fock_virt_total_spin_trace(v,istate)       & 
                                 + fock_core_inactive_total_spin_trace(i,istate)        
       delta_e(i_a,ispin,istate) = 1.d0/delta_e(i_a,ispin,istate)  
      enddo
     enddo
    enddo
    do ispin = 1, 2
     do istate = 1, N_states 
      do i_a = 1, n_act_orb
       a = list_act(i_a)
       do i_b = 1, n_act_orb
        b = list_act(i_b)
        accu_direct = active_int(i_a,1) * active_int(i_b,1)
        contribution_1h2p(istate) += accu_direct  * delta_e(i_a,ispin,istate) * one_body_dm_mo_spin_index(b,a,istate,ispin)
        effective_fock_operator_1h2p(i_a,i_b,ispin,istate) += accu_direct  * delta_e(i_a,ispin,istate) 
       enddo
      enddo
     enddo
    enddo
   enddo
   do i_t = i_v+1, n_virt_orb
    t = list_virt(i_t)
    do i_a = 1, n_act_orb
     a = list_act(i_a)
     active_int(i_a,1) = get_mo_bielec_integral(i,a,t,v,mo_integrals_map) ! direct
     active_int(i_a,2) = get_mo_bielec_integral(i,a,v,t,mo_integrals_map) ! exchange
     do istate = 1, N_states
      do ispin = 1, 2
       delta_e(i_a,ispin,istate) = one_anhil(i_a,ispin,istate)                        &
                                 - fock_virt_total_spin_trace(t,istate)       & 
                                 - fock_virt_total_spin_trace(v,istate)       & 
                                 + fock_core_inactive_total_spin_trace(i,istate)        
       delta_e(i_a,ispin,istate) = 1.d0/delta_e(i_a,ispin,istate)  
      enddo
     enddo
    enddo
    do ispin = 1, 2
     do istate = 1, N_states 
      do i_a = 1, n_act_orb
       a = list_act(i_a)
       do i_b = 1, n_act_orb
        b = list_act(i_b)
        accu_direct = (active_int(i_a,1) - active_int(i_a,2)) * (active_int(i_b,1) - active_int(i_b,2)) + active_int(i_a,1) * active_int(i_b,1)
        contribution_1h2p(istate) += accu_direct  * delta_e(i_a,ispin,istate) * one_body_dm_mo_spin_index(b,a,istate,ispin)
        effective_fock_operator_1h2p(i_a,i_b,ispin,istate) += accu_direct  * delta_e(i_a,ispin,istate) 
       enddo
      enddo
     enddo
    enddo
   enddo
  enddo
 enddo


END_PROVIDER

