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
!   print*, 'elta_e(istate) print', delta_e(istate)
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
   enddo
   do istate = 1, N_states
    pt2(istate) += 2.d0 * accu(istate) * accu(istate) * delta_e(istate)
   enddo
  enddo
 enddo
end

 BEGIN_PROVIDER [double precision, effective_fock_operator_pure_diag_1h1p, (N_states)]

 implicit none
 double precision, intent(out) :: effective_fock_operator_pure_diag_1h1p(N_states)
 integer :: i_i,i_v,i_j,i_k,i_a,i_b,ispin,jspin
 integer :: i,v,j,k,a,b
 integer :: istate
 double precision :: delta_e(N_states)
 double precision :: get_mo_bielec_integral
 double precision :: active_int(n_act_orb,2)
 double precision :: core_inactive_int(n_core_inact_orb,2)
 double precision :: accu(N_states)
 effective_fock_operator_pure_diag_1h1p = 0.d0
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
!   print*, 'elta_e(istate) print', delta_e(istate)
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
   enddo
   do istate = 1, N_states
    effective_fock_operator_pure_diag_1h1p(istate) += 2.d0 * accu(istate) * accu(istate) * delta_e(istate)
   enddo
  enddo
 enddo
END_PROVIDER 
