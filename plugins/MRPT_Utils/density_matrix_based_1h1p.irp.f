subroutine test_1h1p(pt2)
 implicit none
 double precision, intent(out) :: pt2(N_states)
 integer :: i_i,i_v,i_j,i_k,i_a,i_b,ispin,jspin
 integer :: i,v,j,k,a,b
 integer :: istate
 double precision :: delta_e(N_states)
 double precision :: get_mo_bielec_integral
 double precision :: active_int(n_act_orb,2)
 double precision :: core_inactive_int(n_core_inact_orb,2)
 double precision :: accu(N_states)
 pt2 = 0.d0
 do i_i = 1, n_inact_orb
  i = list_inact(i_i)
  do i_v = 1, n_virt_orb
   v = list_virt(i_v)
   delta_e = 0.d0
   accu = 0.d0
   do istate = 1, N_states
    delta_e(istate)  = fock_core_inactive_total_spin_trace(i,istate) & 
                     - fock_virt_total_spin_trace(v,istate)       &
                     + one_anhil_one_creat_inact_virt(i_i,i_v,istate) 
    delta_e(istate) = 1.d0/delta_e(istate)
   enddo
   do istate = 1, N_states
    accu(istate) += mo_mono_elec_integral(i,v) !* mo_mono_elec_integral(i,v) * delta_e(istate)
   enddo
   do i_j = 1, n_core_inact_orb
     j = list_core_inact(i_j)
     core_inactive_int(i_j,1) = get_mo_bielec_integral(i,j,v,j,mo_integrals_map) ! direct
     core_inactive_int(i_j,2) = get_mo_bielec_integral(i,j,j,v,mo_integrals_map) ! exchange
     do istate = 1, N_states
      accu(istate) += (2.d0 * core_inactive_int(i_j,1) - core_inactive_int(i_j,2)) 
     enddo
     do i_k = 1, n_core_inact_orb
 
     enddo
   enddo
   do istate = 1, N_states
    pt2(istate) += 2.d0 * accu(istate) * accu(istate) * delta_e(istate)
   enddo
  enddo
 enddo
end

 BEGIN_PROVIDER [double precision, scalar_core_inact_contrib_1h1p, (N_states)]
&BEGIN_PROVIDER [double precision, array_core_inact_contrib_1h1p, (n_virt_orb, n_inact_orb)]

 implicit none
 integer :: i_i,i_v,i_j,i_k,i_a,i_b,ispin,jspin
 integer :: i,v,j,k,a,b
 integer :: istate
 double precision :: delta_e(N_states)
 double precision :: get_mo_bielec_integral
 double precision :: active_int(n_act_orb,2)
 double precision :: core_inactive_int(n_core_inact_orb,2)
 double precision :: accu,test
 test = 0.d0
 scalar_core_inact_contrib_1h1p = 0.d0
 array_core_inact_contrib_1h1p = 0.d0
 do i_i = 1, n_inact_orb
  i = list_inact(i_i)
  do i_v = 1, n_virt_orb
   v = list_virt(i_v)
   delta_e = 0.d0
   accu = 0.d0
   do istate = 1, N_states
    delta_e(istate)  = fock_core_inactive_total_spin_trace(i,istate) & 
                     - fock_virt_total_spin_trace(v,istate)       &
                     + one_anhil_one_creat_inact_virt(i_i,i_v,istate) 
    delta_e(istate) = 1.d0/delta_e(istate)
   enddo
   accu += mo_mono_elec_integral(i,v) 
   do i_j = 1, n_core_inact_orb
     j = list_core_inact(i_j)
     core_inactive_int(i_j,1) = get_mo_bielec_integral(i,j,v,j,mo_integrals_map) ! direct
     core_inactive_int(i_j,2) = get_mo_bielec_integral(i,j,j,v,mo_integrals_map) ! exchange
     accu += (2.d0 * core_inactive_int(i_j,1) - core_inactive_int(i_j,2)) 
   enddo
   array_core_inact_contrib_1h1p(i_v,i_i) =  accu
   do istate = 1, N_states
    scalar_core_inact_contrib_1h1p(istate) += 2.d0 * accu * accu * delta_e(istate)
   enddo
   test += 2.d0 * array_core_inact_contrib_1h1p(i_v,i_i) **2 * delta_e(1)
  enddo
 enddo
 print*, 'test',test,scalar_core_inact_contrib_1h1p(1)
END_PROVIDER 

 BEGIN_PROVIDER [double precision, effective_active_energies_1h1p, (n_act_orb,N_states)]
&BEGIN_PROVIDER [double precision, effective_coulomb_1h1hp, (n_act_orb,n_act_orb,2,2,N_states)]
 implicit none
 integer :: i_i,i_v,i_j,i_k,i_a,i_b,ispin,jspin
 integer :: i,v,j,k,a,b
 integer :: istate
 double precision :: delta_e(N_states)
 double precision :: get_mo_bielec_integral
 double precision :: active_int(n_act_orb,2)
 double precision :: core_inactive_int(n_core_inact_orb,2)
 double precision :: accu(n_act_orb), accu_double(n_act_orb,n_act_orb,2,2)
 effective_active_energies_1h1p = 0.d0
 effective_coulomb_1h1hp = 0.d0
 do i_i = 1, n_inact_orb
  i = list_inact(i_i) 
  do i_v = 1, n_virt_orb
   v = list_virt(i_v)
   do istate = 1, N_states
    delta_e(istate)  = fock_core_inactive_total_spin_trace(i,istate) & 
                     - fock_virt_total_spin_trace(v,istate)       &
                     + one_anhil_one_creat_inact_virt(i_i,i_v,istate) 
    delta_e(istate) = 1.d0/delta_e(istate)
   enddo
   accu = 0.d0
   do i_a = 1, n_act_orb
    a = list_act(i_a)
    active_int(i_a,1) = get_mo_bielec_integral(i,a,v,a,mo_integrals_map) ! direct
    active_int(i_a,2) = get_mo_bielec_integral(i,a,a,v,mo_integrals_map) ! exchange
    accu(i_a) +=  2.d0 * array_core_inact_contrib_1h1p(i_v,i_i) * (2.d0 * active_int(i_a,1) - active_int(i_a,2))
   enddo
   do istate = 1, N_states
    do i_a = 1, n_act_orb
     effective_active_energies_1h1p(i_a,istate)  += accu(i_a) * delta_e(istate)
    enddo
   enddo
   accu_double = 0.d0
   do jspin = 1,2 
    do ispin = 1, 2
     do i_b = 1, n_act_orb
      do i_a = 1, n_act_orb
       accu_double(i_a,i_b,ispin,jspin) += 2.d0 * active_int(i_a,1) * active_int(i_b,1) & 
                                           - active_int(i_a,1) * active_int(i_b,2)      & 
                                           - active_int(i_a,2) * active_int(i_b,1)        
      enddo
     enddo
    enddo
   enddo
   do ispin = 1, 2
    do i_b = 1, n_act_orb
     do i_a = 1, n_act_orb
      accu_double(i_a,i_b,ispin,ispin) +=  active_int(i_a,2) * active_int(i_b,2)   
     enddo
    enddo
   enddo
   do istate = 1, N_states
    do jspin = 1,2 
     do ispin = 1, 2
      do i_b = 1, n_act_orb
       do i_a = 1, n_act_orb
        effective_coulomb_1h1hp(i_a,i_b,ispin,jspin,istate) += accu_double(i_a,i_b,ispin,jspin) * delta_e(istate)
       enddo
      enddo
     enddo
    enddo
   enddo
  enddo 
 enddo

END_PROVIDER 
