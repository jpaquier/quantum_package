program test_dm
 implicit none
 read_wf = .True.
 touch read_wf 
 call routine
end

subroutine routine
 implicit none
 integer :: i,j,k,l,i_state
 i_state = 1
 do i = 1, mo_tot_num
  do j = 1, mo_tot_num
   do k = 1, mo_tot_num
    do l = 1, mo_tot_num
     if(dabs(two_bod_alpha_beta_mo_contracted(l,k,j,i,i_state) - two_bod_alpha_beta_mo_contracted_parallel(l,k,j,i,i_state)).gt.1.d-10)then
      print*,i,j,k,l
      print*,dabs(two_bod_alpha_beta_mo_contracted(l,k,j,i,i_state) - two_bod_alpha_beta_mo_contracted_parallel(l,k,j,i,i_state)),two_bod_alpha_beta_mo_contracted(l,k,j,i,i_state) , two_bod_alpha_beta_mo_contracted_parallel(l,k,j,i,i_state)
     endif
    enddo
   enddo
  enddo
 enddo
end

subroutine test_expectation_value
 implicit none
 double precision:: test_1,test_2
 integer :: j,k,l,istate
 double precision :: r(3),rho2
 double precision :: two_dm_in_r
 double precision :: wall_1, wall_2

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
     
     test += rho2 * final_weight_functions_at_grid_points(l,k,j)
    enddo
   enddo
  enddo
 enddo
 call wall_time(wall_2)
 print*,'wall time bourrin num inte = ',wall_2 - wall_1
end

