program pouet
 read_wf = .True.
 touch read_wf
! call test_grad_ao
 !call test_grad_lapl_ao
 call test_grad_lapl_mo
!call test_grad_lapl_ao
! call test_grad_density
!call test_v_corel_old
!call test_v_corel_new
!call test_v_corel
!call test_int
end

subroutine test_grad_lapl_ao
 implicit none
 integer :: i,j,k,l,m,n
 double precision :: r(3),rdx_plus(3),rdx_minus(3),accu(3),accu_2(3)
 double precision :: grad_aos_array(ao_num,3),grad_aos_array_bis(ao_num,3)
 double precision :: aos_array(ao_num),aos_array_plus(ao_num),aos_array_minus(ao_num)
 double precision :: lapl_aos_array(ao_num,3),lapl_aos_array_bis(ao_num,3) 
 double precision :: grad_aos_array_plus(3,ao_num),grad_aos_array_minus(3,ao_num) 
 double precision :: dr
print*,'dr,error grad, error lapl'
 do n = 1, 16
 dr = 10d0**(-n)
 r = 0d0
 accu = 0d0
 accu_2 = 0d0
  do j = 1,nucl_num
   do k = 1, n_points_radial_grid -1
    do l = 1 , n_points_integration_angular
     r(1) = grid_points_per_atom(1,l,k,j)
     r(2) = grid_points_per_atom(2,l,k,j)
     r(3) = grid_points_per_atom(3,l,k,j)
     call give_all_aos_and_grad_and_lapl_at_r(r,aos_array,grad_aos_array,lapl_aos_array)
     do m = 1,3
      rdx_plus = r
      rdx_plus(m) = r(m) + dr
      rdx_minus = r
      rdx_minus(m) = r(m) - dr

      call give_all_aos_and_grad_at_r(rdx_plus,aos_array_plus,grad_aos_array_plus)
      call give_all_aos_and_grad_at_r(rdx_minus,aos_array_minus,grad_aos_array_minus) 

      do i = 1, ao_num
        grad_aos_array_bis(i,m) = (aos_array_plus(i) - aos_array_minus(i)) /(2.d0 * dr)
        accu(m) += dabs(grad_aos_array_bis(i,m) - grad_aos_array(i,m)) * final_weight_functions_at_grid_points(l,k,j)

        lapl_aos_array_bis(i,m) = (grad_aos_array_plus(m,i) - grad_aos_array_minus(m,i))/(2.d0 * dr)
        accu_2(m) += dabs(lapl_aos_array_bis(i,m) - lapl_aos_array(i,m)) *final_weight_functions_at_grid_points(l,k,j)
      enddo
     enddo
    enddo
   enddo
  enddo
  write(*,'(100(F16.10,X))')dr,accu(1),accu_2(1),accu_2(2), accu_2(3)
 enddo

end

subroutine test_grad_lapl_mo
 implicit none
 integer :: i,j,k,l,m,n
 double precision :: r(3),rdx_plus(3),rdx_minus(3),accu(3),accu_2(3)
 double precision :: grad_mos_array_bis(ao_num,3)
 double precision :: mos_array_plus(ao_num),mos_array_minus(ao_num)
 double precision :: lapl_mos_array_bis(mo_tot_num,3)
 double precision :: grad_mos_array_plus(3,mo_tot_num),grad_mos_array_minus(3,mo_tot_num)
 double precision :: dr
 print*,'\\\\\\\\\\\\\\\\\'
 print*,' '
 print*,'Test MO'
 print*,'dr,error grad, error lapl'
 do n = 1, 16
  dr = 10d0**(-n)
  r = 0d0
  accu = 0d0
  accu_2= 0d0
  do i = 1, n_points_final_grid
   r(1) = final_grid_points(1,i)
   r(2) = final_grid_points(2,i)
   r(3) = final_grid_points(3,i)
   do m = 1,3
    rdx_plus = r
    rdx_plus(m) = r(m) + dr
    rdx_minus = r
    rdx_minus(m) = r(m) - dr 
    call give_all_mos_and_grad_at_r(rdx_plus,mos_array_plus,grad_mos_array_plus)
    call give_all_mos_and_grad_at_r(rdx_minus,mos_array_minus,grad_mos_array_minus)    
    do j = 1, mo_tot_num
     grad_mos_array_bis(j,m) = (mos_array_plus(j) - mos_array_minus(j))/(2.d0 * dr)
     accu(m) += dabs(grad_mos_array_bis(j,m) - mos_grad_in_r_array(j,i,m)) * final_weight_functions_at_final_grid_points(i)

     lapl_mos_array_bis(j,m) = (grad_mos_array_plus(m,j) - grad_mos_array_minus(m,j))/(2.d0 * dr)
     accu_2(m) += dabs(lapl_mos_array_bis(j,m) - mos_lapl_in_r_array(j,i,m)) * final_weight_functions_at_final_grid_points(i)
    enddo
   enddo
  enddo
  write(33,'(100(F16.10,X))'),dr,accu(1),accu(2),accu(3),accu_2(1),accu_2(2),accu_2(3)
 enddo
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
     call give_all_aos_and_grad_at_r(r,aos_array,grad_aos_array)
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
     call give_all_aos_and_grad_at_r(r,aos_array,grad_aos_array)
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
!provide potential_c_beta_ao_new
 call wall_time(wall_1)
 print*,'time to provide new potential =',dabs(wall_1-wall_0)
end

subroutine test_v_corel
implicit none
integer :: i,j
 print*,'energy_c    = ',energy_c
!print*,'energy_c_new= ',energy_c_new
 print*,'energy_x    = ',energy_x
!print*,'energy_x_new= ',energy_x_new
!provide potential_c_beta_ao
 double precision :: accu
 accu = 0.d0
 do i = 1, ao_num
  do j = 1, ao_num
  !if(dabs(potential_c_beta_ao(j,i,1)-potential_c_beta_ao_new(j,i,1)).gt.1.d-10)then
  ! print*,'AHAHA' 
  ! print*,i,j,dabs(potential_c_beta_ao(j,i,1)-potential_c_beta_ao_new(j,i,1))
    print*,potential_c_beta_ao(j,i,1),potential_x_beta_ao(j,i,1)
  ! accu += dabs(potential_c_beta_ao(j,i,1)-potential_c_beta_ao_new(j,i,1))
  !endif
  enddo
 enddo
 print*,'accu   = ',accu
 print*,'accu/n = ',accu/(dble(ao_num*ao_num))
end

subroutine test_int
 implicit none
 integer            :: power_A(3),power_B(3),power_C(3),power_D(3)
 double precision   :: alpha,beta,delta,gama
 double precision   :: A_center(3), B_center(3), C_center(3), D_center(3)
 double precision :: integral, ao_bielec_integral
  include 'Utils/constants.include.F'

 double precision               :: P_new(0:max_dim,3),P_center(3),fact_p,pp
 double precision               :: Q_new(0:max_dim,3),Q_center(3),fact_q,qq
 integer                        :: iorder_p(3), iorder_q(3)
 double precision               :: general_primitive_integral
 double precision               :: p_inv,q_inv
 double precision               :: ERI

 alpha = 0.5d0
 beta  = 0.5d0
 delta = 0.5d0
 gama  = 0.5d0
 A_center = 0.d0
 B_center = 0.d0
 C_center = 0.d0
 D_center = 0.d0
 power_A = 0
 power_B = 0
 power_C = 0
 power_D = 0
 
 call give_explicit_poly_and_gaussian(P_new,P_center,pp,fact_p,iorder_p,&
     alpha,beta,                 &
     power_A,power_B,A_center,B_center,n_pt_max_integrals)
 p_inv = 1.d0/pp

 call give_explicit_poly_and_gaussian(Q_new,Q_center,qq,fact_q,iorder_q,&
     delta,gama,                 &
     power_C,power_D,C_center,D_center,n_pt_max_integrals)
 q_inv = 1.d0/qq

 print*,'fact_p = ',fact_p
 print*,'fact_q = ',fact_q
 integral = general_primitive_integral(n_pt_max_integrals,              &
     P_new,P_center,fact_p,pp,p_inv,iorder_p,             &
     Q_new,Q_center,fact_q,qq,q_inv,iorder_q)
 print*,'integral = ',integral
 integral = ERI(alpha,beta,delta,gama,power_A(1),power_B(1),power_C(1),power_D(1),power_A(2),power_B(2),power_C(2),power_D(2),power_A(3),power_B(3),power_C(3),power_D(3))
 print*,'integral = ',integral
end

