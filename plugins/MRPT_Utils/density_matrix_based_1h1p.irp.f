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
 BEGIN_DOC
! scalar_core_inact_contrib_1h1p = total contribution of the 1h1p single excitation from all the core-inactive block
! array_core_inact_contrib_1h1p(v,i)= contribution of ONE SINLE EXCITATION 1h1p i-->v from all the core-inactive block
 END_DOC
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
&BEGIN_PROVIDER [double precision, effective_Fock_1h1hp, (n_act_orb,n_act_orb,2,N_states)]
&BEGIN_PROVIDER [double precision, effective_pseudo_Fock_1h1hp, (n_act_orb,2,n_act_orb,n_act_orb,2,N_states)]
&BEGIN_PROVIDER [double precision, effective_pseudo_bielec_1h1hp, (n_act_orb,n_act_orb,2,n_act_orb,n_act_orb,2,N_states)]
&BEGIN_PROVIDER [double precision, effective_active_energies_double_1h1p, (n_act_orb,2,N_states)]
&BEGIN_PROVIDER [double precision, effective_coulomb_double_1h1hp, (n_act_orb,n_act_orb,2,N_states)]
&BEGIN_PROVIDER [double precision, effective_Fock_1h1hp_double, (n_act_orb,n_act_orb,2,N_states)]
&BEGIN_PROVIDER [double precision, effective_pseudo_Fock_double_1h1hp, (n_act_orb,n_act_orb,n_act_orb,2,N_states)]
  BEGIN_DOC
! effective_active_energies_1h1p(i_a) = effective one-body energy coming from ALL the single excitations 1h1p
! effective_coulomb_1h1hp(i_a,i_b,ispin,jspin) = effective two coulomb operator of type \hat{n}_{a,ispin} \hat{n}_{b,ispin} coming from ALL the single excitations 1h1p
! effective_Fock_1h1hp(i_a,i_b) = effective Fock operator of type a^{\dagger}_{b} a_{a}  (does not depend on the spin) coming from ALL the single excitations 1h1p
! effective_pseudo_Fock_1h1hp(i_a,ispin,i_b,i_c,jspin) = effective pseudo Fock operator of type a^{\dagger}_{c,ispin} a_{b,jspin} \hat{n}_{a,ispin}
! effective_pseudo_bielec_1h1hp(i_a,i_b,ispin,i_c,i_d,jspin) = a^{\dagger}_{i_b,ispin} a^{\dagger}_{i_d,jspin} a_{i_c,jspin} a_{i_a,ispin}
  END_DOC
 implicit none
 integer :: i_i,i_v,i_j,i_k,i_a,i_b,ispin,jspin,i_c,i_d
 integer :: i,v,j,k,a,b,c,d
 integer :: istate
 double precision :: delta_e(N_states),delta_e_ab(n_act_orb,n_act_orb,2,N_states),delta_e_ab_bis(n_act_orb,2,n_act_orb,2,N_states)
 double precision :: get_mo_bielec_integral
 double precision :: active_int(n_act_orb,2),active_int_double(n_act_orb,n_act_orb,2)
 double precision :: core_inactive_int(n_core_inact_orb,2)
 double precision :: accu_2(n_act_orb,n_act_orb),accu_3(n_act_orb,2,n_act_orb,n_act_orb,2)
 integer :: other_spin(2)
 other_spin(1) = 2
 other_spin(2) = 1


 effective_active_energies_1h1p = 0.d0
 effective_coulomb_1h1hp = 0.d0
 effective_Fock_1h1hp = 0.d0
 effective_Fock_1h1hp_double = 0.d0
 effective_pseudo_bielec_1h1hp = 0.d0
 effective_active_energies_double_1h1p = 0.d0
 effective_coulomb_double_1h1hp = 0.d0
 effective_pseudo_Fock_double_1h1hp = 0.d0
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
   accu_2 = 0.d0
   do i_a = 1, n_act_orb
    a = list_act(i_a)
    active_int(i_a,1) = get_mo_bielec_integral(i,a,v,a,mo_integrals_map) ! direct
    active_int(i_a,2) = get_mo_bielec_integral(i,a,a,v,mo_integrals_map) ! exchange
    do i_b = 1, n_act_orb
     b = list_act(i_b)
     active_int_double(i_b,i_a,1)     = get_mo_bielec_integral(i,a,v,b,mo_integrals_map) ! V_{ia}^{vb}
     active_int_double(i_b,i_a,2)     = get_mo_bielec_integral(i,a,b,v,mo_integrals_map) ! V_{ia}^{bv}
     accu_2(i_b,i_a) = (2.d0 * active_int_double(i_b,i_a,1) - active_int_double(i_b,i_a,2))
     do ispin = 1,2
      do istate = 1, N_states
       delta_e_ab(i_b,i_a,ispin,istate) = fock_core_inactive_total_spin_trace(i,istate) &  
                                        - fock_virt_total_spin_trace(v,istate)          &
                                        + one_anhil_one_creat(i_b,i_a,ispin,ispin,istate)! a^{\dagger}_{b,ispin} a_{a,ispin}
       delta_e_ab(i_b,i_a,ispin,istate) = 1.d0/delta_e_ab(i_b,i_a,ispin,istate)  
       do jspin = 1, 2
         delta_e_ab_bis(i_b,jspin,i_a,ispin,istate) = fock_core_inactive_total_spin_trace(i,istate)   & 
                                                - fock_virt_total_spin_trace(v,istate)            &
                                                + one_anhil_one_creat(i_b,i_a,jspin,ispin,istate) 
                                                        ! a^{\dagger}_{b,jspin} a_{a,ispin}
       enddo
      enddo
     enddo
    enddo
   enddo
   double precision :: contrib,contrib_test
!!!!!!!!!!!!!!!////////////////// EFFECTIVE ENERGIES and DIAGONAL BIELEC FROM THE SINGLES 
   do istate = 1, N_states
    do i_a = 1, n_act_orb
     effective_active_energies_1h1p(i_a,istate)  += 2.d0 * array_core_inact_contrib_1h1p(i_v,i_i) & 
                                                 * (2.d0 * active_int(i_a,1) - active_int(i_a,2)) * delta_e(istate)
     do ispin = 1, 2
      do jspin = 1, 2
       do i_b = 1, n_act_orb
        effective_coulomb_1h1hp(i_b,i_a,jspin,ispin,istate) +=  (2.d0 * active_int(i_a,1) * active_int(i_b,1) & 
                              - active_int(i_a,2) * active_int(i_b,1) - active_int(i_a,1) * active_int(i_b,2) ) * delta_e(istate)
       enddo
      enddo
      jspin = ispin
      do i_b = 1, n_act_orb
       effective_coulomb_1h1hp(i_b,i_a,jspin,ispin,istate) +=  active_int(i_a,2) * active_int(i_b,2) * delta_e(istate)
      enddo
     enddo
    enddo


!!!!!!!!!!!!!!!////////////////// EFFECTIVE ENERGIES and DIAGONAL BIELEC FROM THE DOUBLES V_{ib}^{va}
    do ispin = 1, 2
     do i_a = 1, n_act_orb
      do i_b = 1, n_act_orb
       if(i_a==i_b)cycle
       contrib =delta_e_ab(i_b,i_a,ispin,istate)                                      &
               *(2.d0 * active_int_double(i_b,i_a,1)  * active_int_double(i_b,i_a,1)  & 
                      - active_int_double(i_b,i_a,2)  * active_int_double(i_b,i_a,1)  &  
                      - active_int_double(i_b,i_a,1)  * active_int_double(i_b,i_a,2)  &  
                      + active_int_double(i_b,i_a,2)  * active_int_double(i_b,i_a,2))   
       
       effective_coulomb_double_1h1hp(i_b,i_a,ispin,istate) += contrib
       effective_active_energies_double_1h1p(i_a,ispin,istate) += contrib
      enddo
     enddo
    enddo

    do ispin = 1, 2
     do i_a = 1, n_act_orb
      do i_b = 1, n_act_orb
       effective_Fock_1h1hp(i_b,i_a,ispin,istate) += accu_2(i_b,i_a) * array_core_inact_contrib_1h1p(i_v,i_i) & 
                                                  *  (delta_e(istate) + delta_e_ab(i_b,i_a,ispin,istate))
        do i_c = 1, n_act_orb
         if(i_c == i_a)cycle
         if(i_c == i_b)cycle
         contrib = delta_e_ab(i_c,i_a,ispin,istate)                                       &
                   *(2.d0 * active_int_double(i_c,i_a,1)  * active_int_double(i_c,i_b,1)  & 
                          - active_int_double(i_c,i_a,2)  * active_int_double(i_c,i_b,1)  &  
                          - active_int_double(i_c,i_a,1)  * active_int_double(i_c,i_b,2)  &  
                          + active_int_double(i_c,i_a,2)  * active_int_double(i_c,i_b,2))   
         effective_Fock_1h1hp_double(i_b,i_a,ispin,istate) += contrib
         effective_pseudo_Fock_double_1h1hp(i_c,i_b,i_a,ispin,istate)  += contrib
         contrib = delta_e_ab(i_a,i_c,ispin,istate)                                       &
                   *(2.d0 * active_int_double(i_a,i_c,1)  * active_int_double(i_b,i_c,1)  & 
                          - active_int_double(i_a,i_c,2)  * active_int_double(i_b,i_c,1)  &  
                          - active_int_double(i_a,i_c,1)  * active_int_double(i_b,i_c,2)  &  
                          + active_int_double(i_a,i_c,2)  * active_int_double(i_b,i_c,2))   
         effective_pseudo_Fock_double_1h1hp(i_c,i_b,i_a,ispin,istate)  += contrib
        enddo
      enddo
     enddo
    enddo
   enddo
   accu_3 = 0.d0
   do jspin = 1,2 
    do ispin = 1, 2
     do i_a = 1, n_act_orb
      do i_b = 1, n_act_orb
       do i_c = 1, n_act_orb
        accu_3(i_c,ispin,i_b,i_a,jspin) += 2.d0 * active_int(i_c,1) * active_int_double(i_b,i_a,1) & 
                                         -        active_int(i_c,2) * active_int_double(i_b,i_a,1) &
                                         -        active_int(i_c,1) * active_int_double(i_b,i_a,2)  
       enddo
      enddo
     enddo
    enddo
   enddo
   do ispin = 1, 2
    do i_a = 1, n_act_orb
     do i_b = 1, n_act_orb
      do i_c = 1, n_act_orb
       accu_3(i_c,ispin,i_b,i_a,ispin) +=  active_int(i_c,2) * active_int_double(i_b,i_a,2)   
      enddo
     enddo
    enddo
   enddo
   do istate = 1, N_states
    do jspin = 1,2 
     do i_b = 1, n_act_orb
      do i_a = 1, n_act_orb
       do ispin = 1, 2
        do i_d = 1, n_act_orb
         effective_pseudo_Fock_1h1hp(i_d,ispin,i_a,i_b,jspin,istate) += accu_3(i_d,ispin,i_a,i_b,jspin) & 
                                                                     * (delta_e(istate) + delta_e_ab(i_a,i_b,jspin,istate) )
         do i_c = 1, n_act_orb
          effective_pseudo_bielec_1h1hp(i_c,i_d,ispin,i_a,i_b,jspin,istate) +=  delta_e_ab(i_a,i_b,jspin,istate)        *                              &
                                                                                 (2.d0 * active_int_double(i_a,i_b,1)  * active_int_double(i_d,i_c,1)  & 
                                                                                       - active_int_double(i_a,i_b,2)  * active_int_double(i_d,i_c,1)  & 
                                                                                       - active_int_double(i_a,i_b,1)  * active_int_double(i_d,i_c,2))   
         enddo
        enddo
       enddo
      enddo
     enddo
    enddo
    do jspin = 1,2 
     do i_b = 1, n_act_orb
      do i_a = 1, n_act_orb
       do ispin = jspin,jspin
        do i_d = 1, n_act_orb
         do i_c = 1, n_act_orb
          effective_pseudo_bielec_1h1hp(i_c,i_d,ispin,i_a,i_b,jspin,istate) += delta_e_ab(i_a,i_b,jspin,istate) &
                                                                                 * active_int_double(i_a,i_b,2) * active_int_double(i_d,i_c,2)   


         enddo
        enddo
       enddo
      enddo
     enddo
    enddo
   enddo

  enddo 
 enddo

END_PROVIDER 


 BEGIN_PROVIDER [double precision, effective_active_energies_double_bis_1h1p, (n_act_orb,2,N_states)]
&BEGIN_PROVIDER [double precision, effective_coulomb_double_bis_1h1hp,  (n_act_orb,n_act_orb,2,2,N_states)]
&BEGIN_PROVIDER [double precision, effective_Fock_1h1hp_double_bis, (n_act_orb,n_act_orb,2,N_states)]
&BEGIN_PROVIDER [double precision, effective_pseudo_Fock_double_bis_1h1hp, (n_act_orb,2,n_act_orb,n_act_orb,2,N_states)]
&BEGIN_PROVIDER [double precision, effective_pseudo_bielec_double_bis_1h1hp, (n_act_orb,n_act_orb,2,n_act_orb,n_act_orb,2,N_states)]
 implicit none
  BEGIN_DOC
! effective_active_energies_1h1p(i_a) = effective one-body energy coming from ALL the single excitations 1h1p
! effective_coulomb_1h1hp(i_a,i_b,ispin,jspin) = effective two coulomb operator of type \hat{n}_{a,ispin} \hat{n}_{b,ispin} coming from ALL the single excitations 1h1p
! effective_Fock_1h1hp(i_a,i_b) = effective Fock operator of type a^{\dagger}_{b} a_{a}  (does not depend on the spin) coming from ALL the single excitations 1h1p
! effective_pseudo_Fock_1h1hp(i_a,ispin,i_b,i_c,jspin) = effective pseudo Fock operator of type a^{\dagger}_{c,ispin} a_{b,jspin} \hat{n}_{a,ispin}
! effective_pseudo_bielec_1h1hp(i_a,i_b,ispin,i_c,i_d,jspin) = a^{\dagger}_{i_b,ispin} a^{\dagger}_{i_d,jspin} a_{i_c,jspin} a_{i_a,ispin}
  END_DOC

 integer :: i_i,i_v,i_j,i_k,i_a,i_b,ispin,jspin,i_c,i_d
 integer :: i,v,j,k,a,b,c,d
 integer :: istate
 double precision :: delta_e_ab_bis(n_act_orb,2,n_act_orb,2,N_states)
 double precision :: get_mo_bielec_integral
 double precision :: active_int_double_bis(n_act_orb,n_act_orb,2)
 logical :: test_1,test_2
 integer :: other_spin(2)
 other_spin(1) = 2
 other_spin(2) = 1

 effective_active_energies_double_bis_1h1p = 0.d0
 effective_coulomb_double_bis_1h1hp = 0.d0
 effective_Fock_1h1hp_double_bis = 0.d0
 effective_pseudo_Fock_double_bis_1h1hp = 0.d0
 effective_pseudo_bielec_double_bis_1h1hp = 0.d0

 do i_i = 1, n_inact_orb
  i = list_inact(i_i) 
  do i_v = 1, n_virt_orb
   v = list_virt(i_v)
   do i_a = 1, n_act_orb
    a = list_act(i_a)
    do i_b = 1, n_act_orb
     b = list_act(i_b)
     active_int_double_bis(i_b,i_a,1) = get_mo_bielec_integral(i,a,b,v,mo_integrals_map)  ! V_{ia}^{bv}
     active_int_double_bis(i_b,i_a,2) = get_mo_bielec_integral(i,a,v,b,mo_integrals_map)  ! V_{ia}^{vb}
     do istate = 1, N_states
      do ispin = 1, 2
       do jspin = 1, 2
         delta_e_ab_bis(i_b,jspin,i_a,ispin,istate) = fock_core_inactive_total_spin_trace(i,istate)   & 
                                                - fock_virt_total_spin_trace(v,istate)            &
                                                + one_anhil_one_creat(i_b,i_a,jspin,ispin,istate) 
                                                        ! a^{\dagger}_{b,jspin} a_{a,ispin}
         delta_e_ab_bis(i_b,jspin,i_a,ispin,istate) = 1.d0/delta_e_ab_bis(i_b,jspin,i_a,ispin,istate)
       enddo
      enddo 
     enddo
    enddo
   enddo
   do istate = 1, N_states 
    do i_a = 1, n_act_orb
     do ispin = 1, 2
!!!!!!!!!!! /////////// ACTIVE ENERGIES
      do i_b = 1, n_act_orb   ! opposite spin contribution i-->b + a-->v
        effective_active_energies_double_bis_1h1p(i_a,ispin,istate) += active_int_double_bis(i_b,i_a,1) * active_int_double_bis(i_b,i_a,1) * delta_e_ab_bis(i_b,other_spin(ispin),i_a,ispin,istate)
      enddo
      do i_b = 1, n_act_orb   ! same spin contribution i-->b + a-->v
     ! if(i_a==i_b)cycle
     !  effective_active_energies_double_bis_1h1p(i_a,ispin,istate) +=  (active_int_double_bis(i_b,i_a,1) - active_int_double_bis(i_b,i_a,2)) * (active_int_double_bis(i_b,i_a,1) -active_int_double_bis(i_b,i_a,2) ) * delta_e_ab_bis(i_b,ispin,i_a,ispin,istate)
     ! 
      enddo
!!!!!!!!!!! /////////// ACTIVE COULOMB OPERATOR   n_{b,\sigma} n_{a,\sigma '}
      do i_b = 1, n_act_orb ! opposite spin contribution 
       effective_coulomb_double_bis_1h1hp(i_a,i_b,ispin,other_spin(ispin),istate) +=  delta_e_ab_bis(i_b,other_spin(ispin),i_a,ispin,istate) & 
                                            *  active_int_double_bis(i_b,i_a,1) *  active_int_double_bis(i_b,i_a,1)
      enddo
      do i_b = 1, n_act_orb ! same spin contribution 
     ! if(i_a==i_b)cycle
     ! effective_coulomb_double_bis_1h1hp(i_a,i_b,ispin,ispin,istate) +=  delta_e_ab_bis(i_b,ispin,i_a,ispin,istate) & 
     !                                      *  (active_int_double_bis(i_b,i_a,1) - active_int_double_bis(i_b,i_a,2)) *  (active_int_double_bis(i_b,i_a,1) - active_int_double_bis(i_b,i_a,2))
      enddo
     enddo
    enddo
!!!!!!!!!!! /////////// ACTIVE FOCK OPERATOR   a^{\dagger}_{b,ispin}  a_{a,ispin}
    do i_c = 1, n_act_orb
     do ispin = 1, 2
      do i_b = 1, n_act_orb   
       do i_a = 1, n_act_orb ! opposite spin contribution i-->a + b-->v
        effective_Fock_1h1hp_double_bis(i_b,i_c,ispin,istate) += active_int_double_bis(i_a,i_b,1) * active_int_double_bis(i_a,i_c,1) * delta_e_ab_bis(i_a,other_spin(ispin),i_b,ispin,istate) 
       enddo
       do i_a = 1, n_act_orb ! same spin contribution i-->c + b-->v
     !  if(i_a==i_b)cycle
     !  if(i_a==i_c)cycle
     !  effective_Fock_1h1hp_double_bis(i_b,i_c,ispin,istate) += (active_int_double_bis(i_a,i_b,1) - active_int_double_bis(i_a,i_b,2)) * (active_int_double_bis(i_a,i_c,1) - active_int_double_bis(i_a,i_c,2)) * delta_e_ab_bis(i_a,ispin,i_b,ispin,istate) 
       enddo
      enddo
     enddo
    enddo
!!!!!!!!!!! /////////// ACTIVE PSEUDO FOCK OPERATOR a^{\dagger}_{a,ispin} a_{c,ispin} \hat{n}_{b,jspin}
    do i_a = 1, n_act_orb
     do ispin = 1, 2
      do i_c = 1, n_act_orb ! opposite spin contribution 
       do i_b = 1, n_act_orb
        effective_pseudo_Fock_double_bis_1h1hp(i_b,other_spin(ispin),i_a,i_c,ispin,istate) += (active_int_double_bis(i_a,i_b,1) * active_int_double_bis(i_c,i_b,1) * delta_e_ab_bis(i_a,ispin,i_b,other_spin(ispin),istate)) & 
                                                                                            + (active_int_double_bis(i_b,i_a,1) * active_int_double_bis(i_b,i_c,1) * delta_e_ab_bis(i_b,other_spin(ispin),i_c,ispin,istate))
       enddo
      enddo
     enddo
    enddo

    do i_a = 1, n_act_orb
     do ispin = 1, 2
      do i_c = 1, n_act_orb ! same spin contribution 
       if(i_a==i_b)cycle
       do i_b = 1, n_act_orb
        if(i_b==i_a)cycle
        if(i_b==i_c)cycle
    !   effective_pseudo_Fock_double_bis_1h1hp(i_b,ispin,i_a,i_c,ispin,istate) += & 
    !  (active_int_double_bis(i_a,i_b,1) - active_int_double_bis(i_a,i_b,2)) * (active_int_double_bis(i_c,i_b,1) - active_int_double_bis(i_c,i_b,2)) & 
    ! * delta_e_ab_bis(i_a,ispin,i_b,ispin,istate) & 
    !+ (active_int_double_bis(i_b,i_a,1) - active_int_double_bis(i_b,i_a,2)) * (active_int_double_bis(i_b,i_c,1) - active_int_double_bis(i_b,i_c,2)) & 
    ! * delta_e_ab_bis(i_b,ispin,i_c,ispin,istate) 
       enddo
      enddo
     enddo
    enddo

!!!!!!!!!!! /////////// ACTIVE EXTRA DIAGONAL COULOMB OPERATOR a^{\dagger}_{b,ispin} a^{\dager}_{d,jspin} a_{c,jspin} a_{a,ispin}
    do ispin = 1, 2
     do i_b = 1, n_act_orb
      do i_a = 1, n_act_orb 
       test_1 = ((i_a == 2.and.i_b==1.and.ispin==1))
       if(i_a==i_b)cycle
       jspin = other_spin(ispin)
       do i_d = 1, n_act_orb
        do i_c = 1, n_act_orb
        if(i_c==i_d)cycle
         effective_pseudo_bielec_double_bis_1h1hp(i_c,i_d,jspin,i_a,i_b,ispin,istate) += active_int_double_bis(i_c,i_b,1) * active_int_double_bis(i_d,i_a,1) * delta_e_ab_bis(i_d,jspin,i_a,ispin,istate)! V_{ib}^{cv} * V_{ia}^{dv}/delta_a^d
        enddo
       enddo
      enddo
     enddo
    enddo 

   !do ispin = 1, 2
   ! do i_b = 1, n_act_orb
   !  do i_a = 1, n_act_orb 
   !   test_1 = ((i_a == 2.and.i_b==1.and.ispin==1))
   !   if(i_a==i_b)cycle
   !   jspin = ispin
   !   do i_d = 1, n_act_orb
   !    do i_c = 1, n_act_orb
   !    if(i_c==i_d)cycle
   !    if(i_c==i_a)cycle
   !    if(i_c==i_b)cycle
   !     effective_pseudo_bielec_double_bis_1h1hp(i_c,i_d,jspin,i_a,i_b,ispin,istate) += & 
   !     (active_int_double_bis(i_c,i_b,1) - active_int_double_bis(i_c,i_b,2)) * &  !  V_{ib}^{cv} - V_{ib}^{vc} 
   !     (active_int_double_bis(i_d,i_a,1) - active_int_double_bis(i_d,i_a,2)) * &  !* V_{ia}^{dv} - V_{ia}^{vd}
   !     delta_e_ab_bis(i_d,jspin,i_a,ispin,istate)                                 ! /delta_a^d
   !    enddo
   !   enddo
   !  enddo
   ! enddo
   !enddo 


   enddo
  enddo
 enddo


 END_PROVIDER



 BEGIN_PROVIDER [double precision, effective_coulomb_diag_1h1p_total,  (n_act_orb,n_act_orb,2,2,N_states)]
&BEGIN_PROVIDER [double precision, effective_energies_1h1p_total,  (n_act_orb,2,N_states)]
 implicit none
 integer :: a,b,ispin,jspin,i_a,i_b,i_c,c,i_d,d
 integer :: i,j,istate,k
 double precision :: contrib_diag_one_bod 
 double precision :: contrib_diag_two_bod 
 effective_coulomb_diag_1h1p_total = 0.d0
 effective_energies_1h1p_total = 0.d0
 do istate = 1, N_states
  do ispin = 1, 2
   do i_a = 1, n_act_orb
    contrib_diag_one_bod =  effective_active_energies_double_1h1p(i_a,ispin,istate)        & ! (i-->v)(a-->b)
                          + effective_active_energies_1h1p(i_a,istate)                     & ! (i-->v)(a<->a)
                          + effective_active_energies_double_bis_1h1p(i_a,ispin,istate)      ! (i-->a)(b-->v)
    effective_energies_1h1p_total(i_a,ispin,istate) = contrib_diag_one_bod 
    do jspin = 1,2
     do i_b = 1, n_act_orb
      effective_coulomb_diag_1h1p_total(i_b,i_a,jspin,ispin,istate) += effective_coulomb_1h1hp(i_a,i_b,jspin,ispin,istate)  & 
                                                                     - effective_coulomb_double_bis_1h1hp(i_b,i_a,jspin,ispin,istate) 
     enddo
    enddo
    do i_b = 1, n_act_orb
     effective_coulomb_diag_1h1p_total(i_b,i_a,ispin,ispin,istate) -= effective_coulomb_double_1h1hp(i_b,i_a,ispin,istate) 
    enddo
   enddo
  enddo
 enddo



 END_PROVIDER
