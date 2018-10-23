program test_dm
 implicit none
 read_wf = .True.
 touch read_wf 
 call routine
!call test_expectation_value_2
end

subroutine routine
 implicit none
 integer :: i,j,k,l,i_state
 i_state = 1
 do i = 1, mo_tot_num
  do j = 1, mo_tot_num
   do k = 1, mo_tot_num
    do l = 1, mo_tot_num
     if(dabs(two_bod_alpha_beta_mo_contracted_serial(l,k,j,i,i_state) - two_bod_alpha_beta_mo_contracted(l,k,j,i,i_state)).gt.1.d-10)then
      print*,i,j,k,l
      print*,dabs(two_bod_alpha_beta_mo_contracted_serial(l,k,j,i,i_state) - two_bod_alpha_beta_mo_contracted(l,k,j,i,i_state)),two_bod_alpha_beta_mo_contracted_serial(l,k,j,i,i_state) , two_bod_alpha_beta_mo_contracted(l,k,j,i,i_state)
     endif
    enddo
   enddo
  enddo
 enddo
end

subroutine test_expectation_value
 implicit none
 double precision:: test_1,test_2,test
 integer :: j,k,l,istate
 double precision :: r(3),rho2
 double precision :: two_dm_in_r
 double precision :: wall_1, wall_2,test_1_contract_integral

do istate = 1, N_states
 test = 0.d0
 r(1) = 0.d0
 r(2) = 0.d0
 r(3) = 0.d0
 call wall_time(wall_1)
  do j = 1, nucl_num
   do k = 1, n_points_radial_grid  -1
    do l = 1, n_points_integration_angular
     r(1) = grid_points_per_atom(1,l,k,j)
     r(2) = grid_points_per_atom(2,l,k,j)
     r(3) = grid_points_per_atom(3,l,k,j)
     test = test_1_contract_integral(r)
     test_1 += test * final_weight_functions_at_grid_points(l,k,j)
    enddo
   enddo
  enddo
 enddo
 call wall_time(wall_2)
 print*,'wall time bourrin num inte = ',wall_2 - wall_1
 print*,'test_1            = ',test_1
 print*,'psi_energy_bielec = ',psi_energy_bielec
end

subroutine test_expectation_value_2
 implicit none
 double precision:: integral_1,integral_2,test_1,test_2,integral_3
 integer :: j1,k1,l1
 integer :: j2,k2,l2,istate
 double precision :: r1(3),r2(3)
 double precision :: two_dm_in_r
 double precision :: wall_1, wall_2,test_1_contract,test_1_contract_integral

 integral_1 = 0.d0
 integral_2 = 0.d0
 istate = 1
 do j1 = 1, nucl_num
  do k1 = 1, n_points_radial_grid  -1
   do l1 = 1, n_points_integration_angular
    r1(1) = grid_points_per_atom(1,l1,k1,j1)
    r1(2) = grid_points_per_atom(2,l1,k1,j1)
    r1(3) = grid_points_per_atom(3,l1,k1,j1)
    test_1 = test_1_contract_integral(r1)
    integral_1 += test_1 * final_weight_functions_at_grid_points(l1,k1,j1)
    print*,k1
    integral_3 = 0.d0
    do j2 = 1, nucl_num
     do k2 = 1, n_points_radial_grid  -1
      do l2 = 1, n_points_integration_angular
       r2(1) = grid_points_per_atom(1,l2,k2,j2)
       r2(2) = grid_points_per_atom(2,l2,k2,j2)
       r2(3) = grid_points_per_atom(3,l2,k2,j2)
       test_2 = test_1_contract(r1,r2)
       integral_2 += test_2 * final_weight_functions_at_grid_points(l1,k1,j1) * final_weight_functions_at_grid_points(l2,k2,j2)
       integral_3 += test_2 * final_weight_functions_at_grid_points(l2,k2,j2)
      enddo
     enddo
    enddo
    print*,integral_3,test_1

   enddo
  enddo
 enddo
 print*,'integral_1        = ',integral_1
 print*,'integral_2        = ',integral_2
 print*,'psi_energy_bielec = ',psi_energy_bielec
end

