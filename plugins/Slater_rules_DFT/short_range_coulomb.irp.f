 BEGIN_PROVIDER [double precision, short_range_Hartree_operator, (mo_tot_num,mo_tot_num)]
&BEGIN_PROVIDER [double precision, short_range_Hartree]
 implicit none
 BEGIN_DOC
! short_range_Hartree_operator(i,j) = \int dr i(r)j(r) \int r' \rho(r') W_{ee}^{sr}
! short_range_Hartree = 0.5 * \sum_{i,j} \rho_{ij} short_range_Hartree_operator(i,j) 
!                     = 0.5 * \int dr \int r' \rho(r) \rho(r') W_{ee}^{sr}
 END_DOC
 integer :: i,j,k,l,m,n
 double precision :: get_mo_bielec_integral,get_mo_bielec_integral_erf
 double precision :: integral, integral_erf, contrib
 short_range_Hartree_operator = 0.d0
 short_range_Hartree = 0.d0
 do i = 1, mo_tot_num
  do j = 1, mo_tot_num
   if(dabs(one_body_dm_mo(i,j)).le.1.d-10)cycle
   do k = 1, mo_tot_num
    do l = 1, mo_tot_num
     integral = get_mo_bielec_integral(i,k,j,l,mo_integrals_map) ! <ik|jl> = (ij|kl)
     integral_erf = get_mo_bielec_integral_erf(i,k,j,l,mo_integrals_erf_map)
     contrib = one_body_dm_mo(i,j) * (integral  - integral_erf)
     short_range_Hartree_operator(l,k) += contrib 
     short_range_Hartree += contrib * one_body_dm_mo(k,l) 
    enddo
   enddo
  enddo
 enddo
 short_range_Hartree = short_range_Hartree * 0.5d0
 print*, 'short_range_Hartree',short_range_Hartree
END_PROVIDER


 BEGIN_PROVIDER [double precision, effective_one_e_potential, (mo_tot_num, mo_tot_num,N_states)]
&BEGIN_PROVIDER [double precision, effective_one_e_potential_without_kin, (mo_tot_num, mo_tot_num,N_states)]
 implicit none
 integer :: i,j,i_state
 effective_one_e_potential = 0.d0
 BEGIN_DOC 
! effective_one_e_potential(i,j) = <i| v_{H}^{sr} |j> + <i| h_{core} |j> + <i| v_{xc} |j>
! Taking the expectation value does not provide any energy
! but effective_one_e_potential(i,j) is the potential coupling DFT and WFT part to be used in any WFT calculation
 END_DOC
 do i_state = 1, N_states
  do i = 1, mo_tot_num
   do j = 1, mo_tot_num
    effective_one_e_potential(i,j,i_state) = short_range_Hartree_operator(i,j) + mo_nucl_elec_integral(i,j) + mo_kinetic_integral(i,j) & 
                                   + 0.5d0 * (potential_x_alpha_mo(i,j,i_state) + potential_c_alpha_mo(i,j,i_state)                               &
                                   +          potential_x_beta_mo(i,j,i_state) + potential_c_beta_mo(i,j,i_state)   )
    effective_one_e_potential_without_kin(i,j,i_state) = short_range_Hartree_operator(i,j) + mo_nucl_elec_integral(i,j)  & 
                                   + 0.5d0 * (potential_x_alpha_mo(i,j,i_state) + potential_c_alpha_mo(i,j,i_state)                               &
                                   +          potential_x_beta_mo(i,j,i_state) + potential_c_beta_mo(i,j,i_state)   )
   enddo
  enddo
 enddo
END_PROVIDER 

BEGIN_PROVIDER [double precision, one_e_energy_potential, (mo_tot_num, mo_tot_num)]
 implicit none
 integer :: i,j,i_state
 BEGIN_DOC 
! one_e_energy_potential(i,j) = <i|h_{core}|j> + \int dr i(r)j(r) \int r' \rho(r') W_{ee}^{sr}
! If one take the expectation value over Psi, one gets the total one body energy
 END_DOC
  do i = 1, mo_tot_num
   do j = 1, mo_tot_num
    one_e_energy_potential(i,j) = mo_nucl_elec_integral(i,j) + mo_kinetic_integral(i,j) + short_range_Hartree_operator(i,j) * 0.5d0
   enddo
  enddo

END_PROVIDER 

BEGIN_PROVIDER [double precision, Fock_matrix_expectation_value]
 implicit none
  call get_average(effective_one_e_potential,one_body_dm_mo,Fock_matrix_expectation_value)

END_PROVIDER 

BEGIN_PROVIDER [double precision, Trace_v_xc]
 implicit none
 integer :: i,j
 double precision :: tmp(mo_tot_num,mo_tot_num)
  tmp = 0.d0
  do i = 1, mo_tot_num
   do j = 1, mo_tot_num
     tmp(i,j) =   + 0.5d0 * (potential_x_alpha_mo(i,j,1) + potential_c_alpha_mo(i,j,1)&
                  +     potential_x_beta_mo(i,j,1) + potential_c_beta_mo(i,j,1)   )
   enddo
  enddo
  call get_average(tmp,one_body_dm_mo,Trace_v_xc)

END_PROVIDER 
