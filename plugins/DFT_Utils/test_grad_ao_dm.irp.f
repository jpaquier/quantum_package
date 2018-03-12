program pouet
 read_wf = .True.
 touch read_wf
!call test_grad_density
!call test_v_corel_old
!call test_v_corel_new
 call test_v_corel
end

subroutine test_grad_ao
 implicit none
 integer :: i,j,k,l,m,n
 double precision :: r(3),rdx_plus(3),rdx_minus(3),accu(3)
 double precision :: grad_aos_array(3,ao_num),grad_aos_array_bis(3,ao_num)
 double precision :: aos_array(ao_num),aos_array_plus(ao_num),aos_array_minus(ao_num)
 double precision :: dr

 do n = 1, 16
 dr = 10d0**(-n)
 r = 0d0
 accu = 0d0
  do j = 1,nucl_num
   do k = 1, n_points_radial_grid -1
    do l = 1 , n_points_integration_angular
     r(1) = grid_points_per_atom(1,l,k,j)
     r(2) = grid_points_per_atom(2,l,k,j)
     r(3) = grid_points_per_atom(3,l,k,j)
     call give_all_aos_and_grad_at_r_new(r,aos_array,grad_aos_array)
     do m = 1,3
      rdx_plus = r
      rdx_plus(m) = r(m) + dr
      rdx_minus = r
      rdx_minus(m) = r(m) - dr
      call give_all_aos_at_r(rdx_plus,aos_array_plus)
      call give_all_aos_at_r(rdx_minus,aos_array_minus)
       
      do i = 1, ao_num
        grad_aos_array_bis(m,i) = (aos_array_plus(i) - aos_array_minus(i)) / (2.d0 * dr)
        accu(m) += dabs(grad_aos_array_bis(m,i) - grad_aos_array(m,i)) *  final_weight_functions_at_grid_points(l,k,j)
      enddo
     enddo
    enddo
   enddo
  enddo
  print*,dr,accu(1)
 enddo

end


subroutine test_grad_density
 implicit none
 integer :: i,j,k,l,m,n
 double precision :: r(3),rdx_plus(3),rdx_minus(3),accu_a(3),accu_b(3)
 double precision :: grad_aos_array(3,ao_num),grad_aos_array_bis(3,ao_num)
 double precision :: aos_array(ao_num),aos_array_plus(ao_num),aos_array_minus(ao_num)
 double precision :: dr, dm_a_plus,dm_b_plus,dm_a_minus,dm_b_minus,dm_a,dm_b,grad_dm_a(3),grad_dm_b(3),grad_dm_a_bis(3),grad_dm_b_bis(3)

 do n = 1, 16
 dr = 10d0**(-n)
 r = 0d0
 accu_a = 0d0
 accu_b = 0d0
 
  do j = 1,nucl_num
   do k = 1, n_points_radial_grid -1
    do l = 1 , n_points_integration_angular
     r(1) = grid_points_per_atom(1,l,k,j)
     r(2) = grid_points_per_atom(2,l,k,j)
     r(3) = grid_points_per_atom(3,l,k,j)
     call density_and_grad_alpha_beta_and_all_aos_and_grad_aos_at_r(r,dm_a,dm_b, grad_dm_a, grad_dm_b, aos_array, grad_aos_array)
     call give_all_aos_and_grad_at_r_new(r,aos_array,grad_aos_array)
     do m = 1,3
      rdx_plus = r
      rdx_plus(m) = r(m) + dr
      rdx_minus = r
      rdx_minus(m) = r(m) - dr
      call dm_dft_alpha_beta_at_r(rdx_plus,dm_a_plus,dm_b_plus)
      call dm_dft_alpha_beta_at_r(rdx_minus,dm_a_minus,dm_b_minus)
       
      grad_dm_a_bis(m) = (dm_a_plus - dm_a_minus) / (2.d0 * dr)
      grad_dm_b_bis(m) = (dm_b_plus - dm_b_minus) / (2.d0 * dr)
      accu_a(m) += dabs(grad_dm_a_bis(m) - grad_dm_a(m)) *  final_weight_functions_at_grid_points(l,k,j)
      accu_b(m) += dabs(grad_dm_b_bis(m) - grad_dm_b(m)) *  final_weight_functions_at_grid_points(l,k,j)
     enddo
    enddo
   enddo
  enddo
  print*,dr,accu_a(1), accu_b(1)
 enddo

end


subroutine test_v_corel_old
implicit none
integer :: i,j
double precision :: wall_0,wall_1
 call wall_time(wall_0)
 provide potential_c_beta_ao
 call wall_time(wall_1)
 print*,'time to provide old potential =',dabs(wall_1-wall_0)
end

subroutine test_v_corel_new
implicit none
integer :: i,j
double precision :: wall_0,wall_1
 call wall_time(wall_0)
 provide potential_c_beta_ao_new
 call wall_time(wall_1)
 print*,'time to provide new potential =',dabs(wall_1-wall_0)
end

subroutine test_v_corel
implicit none
integer :: i,j
 print*,'energy_c    = ',energy_c
 print*,'energy_c_new= ',energy_c_new
 print*,'energy_x    = ',energy_x
 print*,'energy_x_new= ',energy_x_new
!provide potential_c_beta_ao
 double precision :: accu
 accu = 0.d0
 do i = 1, ao_num
  do j = 1, ao_num
   if(dabs(potential_c_beta_ao(j,i,1)-potential_c_beta_ao_new(j,i,1)).gt.1.d-10)then
    print*,'AHAHA' 
    print*,i,j,dabs(potential_c_beta_ao(j,i,1)-potential_c_beta_ao_new(j,i,1))
    print*,potential_c_beta_ao(j,i,1),potential_c_beta_ao_new(j,i,1)
    accu += dabs(potential_c_beta_ao(j,i,1)-potential_c_beta_ao_new(j,i,1))
   endif
  enddo
 enddo
 print*,'accu   = ',accu
 print*,'accu/n = ',accu/(dble(ao_num*ao_num))
end
