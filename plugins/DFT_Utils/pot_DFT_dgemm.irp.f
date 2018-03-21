 BEGIN_PROVIDER [double precision, energy_x_new, (N_states)]
&BEGIN_PROVIDER [double precision, energy_c_new, (N_states)]
&BEGIN_PROVIDER [double precision, potential_x_alpha_ao_new,(ao_num,ao_num,N_states)]
&BEGIN_PROVIDER [double precision, potential_x_beta_ao_new,(ao_num,ao_num,N_states)]
&BEGIN_PROVIDER [double precision, potential_c_alpha_ao_new,(ao_num,ao_num,N_states)]
&BEGIN_PROVIDER [double precision, potential_c_beta_ao_new,(ao_num,ao_num,N_states)]

 implicit none
 integer :: i,j,k,l,istate
 integer :: m,n,p,q
 double precision, allocatable :: aos_array(:)
 double precision, allocatable :: r(:)
 double precision :: rho_a(N_states),rho_b(N_states),ex(N_states),ec(N_states)
 double precision :: vx_rho_a(N_states),vx_rho_b(N_states),vc_rho_a(N_states),vc_rho_b(N_states)
 double precision, allocatable :: tmp_c_a(:,:,:),tmp_c_b(:,:,:),tmp_x_a(:,:,:),tmp_x_b(:,:,:)
 double precision :: weight
 double precision :: dvx_rho_a(3,N_states), dvx_rho_b(3,N_states)
 double precision :: dvc_rho_a(3,N_states), dvc_rho_b(3,N_states)
 double precision :: vx_grad_rho_a_2(N_states),vx_grad_rho_b_2(N_states),vx_grad_rho_a_b(N_states)
 double precision :: vc_grad_rho_a_2(N_states),vc_grad_rho_b_2(N_states),vc_grad_rho_a_b(N_states)
 double precision :: grad_rho_a_2(N_states),grad_rho_b_2(N_states),grad_rho_a_b(N_states)
 double precision :: grad_rho_a(3,N_states),grad_rho_b(3,N_states)
 double precision :: grad_aos_array(3,ao_num)
 double precision :: grad_aos_array_transpose(ao_num,3)
 double precision :: dtmp_x_a(3,N_states),dtmp_x_b(3,N_states)
 double precision :: dtmp_c_a(3,N_states),dtmp_c_b(3,N_states)
 double precision :: ao_matrix(ao_num,n_points_integration_angular) 
 double precision :: ao_matrix_vc_a(ao_num,n_points_integration_angular,N_states) 
 double precision :: ao_matrix_vx_a(ao_num,n_points_integration_angular,N_states) 
 double precision :: ao_matrix_vc_b(ao_num,n_points_integration_angular,N_states) 
 double precision :: ao_matrix_vx_b(ao_num,n_points_integration_angular,N_states) 
 double precision ::  contrib_grad_xa(3,N_states)
 double precision ::  contrib_grad_xb(3,N_states)
 double precision ::  contrib_grad_ca(3,N_states)
 double precision ::  contrib_grad_cb(3,N_states)
 double precision :: ao_matrix_dvc_a(ao_num,n_points_integration_angular,3,N_states) 
 double precision :: ao_matrix_dvx_a(ao_num,n_points_integration_angular,3,N_states) 
 double precision :: ao_matrix_dvc_b(ao_num,n_points_integration_angular,3,N_states) 
 double precision :: ao_matrix_dvx_b(ao_num,n_points_integration_angular,3,N_states) 
 double precision :: grad_ao_matrix(ao_num,n_points_integration_angular,3) 
 double precision :: grad_ao_matrix_dvc_a(ao_num,n_points_integration_angular,3,N_states) 
 double precision :: grad_ao_matrix_dvx_a(ao_num,n_points_integration_angular,3,N_states) 
 double precision :: grad_ao_matrix_dvc_b(ao_num,n_points_integration_angular,3,N_states) 
 double precision :: grad_ao_matrix_dvx_b(ao_num,n_points_integration_angular,3,N_states) 
  energy_x_new = 0.d0
  energy_c_new = 0.d0
  potential_c_alpha_ao_new = 0.d0
  potential_x_alpha_ao_new = 0.d0
  potential_c_beta_ao_new = 0.d0
  potential_x_beta_ao_new = 0.d0

  do j = 1, nucl_num
   do k = 1, n_points_radial_grid  -1
    allocate(tmp_c_a(ao_num,ao_num,N_states),tmp_c_b(ao_num,ao_num,N_states),tmp_x_a(ao_num,ao_num,N_states),tmp_x_b(ao_num,ao_num,N_states),aos_array(ao_num),r(3))
    tmp_c_a = 0.d0
    tmp_c_b = 0.d0
    tmp_x_a = 0.d0
    tmp_x_b = 0.d0

    if(DFT_TYPE=="LDA")then
     call create_matrix_for_LDA_like_potential_update_energy(k,j,ao_matrix,ao_matrix_vc_a,ao_matrix_vc_b,ao_matrix_vx_a,ao_matrix_vx_b,energy_x_new,energy_c_new)
    do l = 1, n_points_integration_angular 
   !  r(1) = grid_points_per_atom(1,l,k,j)
   !  r(2) = grid_points_per_atom(2,l,k,j)
   !  r(3) = grid_points_per_atom(3,l,k,j)
   !  weight = final_weight_functions_at_grid_points(l,k,j)
   !  call dm_dft_alpha_beta_and_all_aos_at_r(r,rho_a,rho_b,aos_array)
   !  call LDA_type_functional(rho_a,rho_b,vx_rho_a,vx_rho_b,vc_rho_a,vc_rho_b,ex,ec)

   !  do istate = 1, N_states
   !   energy_x_new(istate) += weight *  ex(istate) 
   !   energy_c_new(istate) += weight *  ec(istate)
   !   vx_rho_a(istate) *= weight
   !   vc_rho_a(istate) *= weight
   !   vx_rho_b(istate) *= weight
   !   vc_rho_b(istate) *= weight
   !   do i = 1 , ao_num 
   !    ao_matrix(i,l)             = aos_array(i)
   !    ao_matrix_vc_a(i,l,istate) = vc_rho_a(istate) * aos_array(i)
   !    ao_matrix_vc_b(i,l,istate) = vc_rho_b(istate) * aos_array(i)
   !    ao_matrix_vx_a(i,l,istate) = vx_rho_a(istate) * aos_array(i)
   !    ao_matrix_vx_b(i,l,istate) = vx_rho_b(istate) * aos_array(i)
   !   enddo
   !  enddo
    enddo 
     call update_pots_scalar_dgemm(ao_matrix_vc_a,ao_matrix_vx_a,ao_matrix_vc_b,ao_matrix_vx_b,ao_matrix,potential_c_alpha_ao_new,potential_x_alpha_ao_new,potential_c_beta_ao_new,potential_x_beta_ao_new)

    else if (DFT_TYPE=="GGA")then
     do l = 1, n_points_integration_angular 
      r(1) = grid_points_per_atom(1,l,k,j)
      r(2) = grid_points_per_atom(2,l,k,j)
      r(3) = grid_points_per_atom(3,l,k,j)
      call density_and_grad_alpha_beta_and_all_aos_and_grad_aos_at_r(r,rho_a,rho_b, grad_rho_a, grad_rho_b, aos_array, grad_aos_array)
      grad_rho_a_2 = 0.d0
      grad_rho_b_2 = 0.d0
      grad_rho_a_b = 0.d0
      do istate = 1, N_states
       do m = 1, 3
        grad_rho_a_2(istate) += grad_rho_a(m,istate) * grad_rho_a(m,istate)
        grad_rho_b_2(istate) += grad_rho_b(m,istate) * grad_rho_b(m,istate)
        grad_rho_a_b(istate) += grad_rho_a(m,istate) * grad_rho_b(m,istate)
       enddo
      enddo
      call GGA_type_functionals(rho_a,rho_b,grad_rho_a_2,grad_rho_b_2,grad_rho_a_b,ex,vx_rho_a,vx_rho_b,vx_grad_rho_a_2,vx_grad_rho_b_2,vx_grad_rho_a_b, &  
                                                                                   ec,vc_rho_a,vc_rho_b,vc_grad_rho_a_2,vc_grad_rho_b_2,vc_grad_rho_a_b )
      do i = 1, ao_num
       ao_matrix(i,l) = aos_array(i)
       do m = 1,3
        grad_aos_array_transpose(i,m)=grad_aos_array(m,i)
       enddo
      enddo
      do m = 1,3
       do i = 1,ao_num
        grad_ao_matrix(i,l,m) = grad_aos_array_transpose(i,m)
       enddo
      enddo
      weight = final_weight_functions_at_grid_points(l,k,j) 
      do istate=1,N_states
       energy_x_new(istate) += weight *  ex(istate)
       energy_c_new(istate) += weight *  ec(istate)
       vx_rho_a(istate) *= weight
       vc_rho_a(istate) *= weight
       vx_rho_b(istate) *= weight
       vc_rho_b(istate) *= weight

       do m= 1,3
        contrib_grad_xa(m,istate) = weight * (2.d0 * vx_grad_rho_a_2(istate) *  grad_rho_a(m,istate) + vx_grad_rho_a_b(istate)  * grad_rho_b(m,istate)) 
        contrib_grad_xb(m,istate) = weight * (2.d0 * vx_grad_rho_b_2(istate) *  grad_rho_b(m,istate) + vx_grad_rho_a_b(istate)  * grad_rho_a(m,istate)) 
        contrib_grad_ca(m,istate) = weight * (2.d0 * vc_grad_rho_a_2(istate) *  grad_rho_a(m,istate) + vc_grad_rho_a_b(istate)  * grad_rho_b(m,istate)) 
        contrib_grad_cb(m,istate) = weight * (2.d0 * vc_grad_rho_b_2(istate) *  grad_rho_b(m,istate) + vc_grad_rho_a_b(istate)  * grad_rho_a(m,istate)) 
       enddo 
       do i = 1, ao_num
        ao_matrix_vc_a(i,l,istate) = vc_rho_a(istate) * aos_array(i)
        ao_matrix_vc_b(i,l,istate) = vc_rho_b(istate) * aos_array(i)
        ao_matrix_vx_a(i,l,istate) = vx_rho_a(istate) * aos_array(i)
        ao_matrix_vx_b(i,l,istate) = vx_rho_b(istate) * aos_array(i)
        do m = 1,3
         ao_matrix_dvc_a(i,l,m,istate) = contrib_grad_ca(m,istate) * aos_array(i) 
         ao_matrix_dvc_b(i,l,m,istate) = contrib_grad_cb(m,istate) * aos_array(i) 
         ao_matrix_dvx_a(i,l,m,istate) = contrib_grad_xa(m,istate) * aos_array(i) 
         ao_matrix_dvx_b(i,l,m,istate) = contrib_grad_xb(m,istate) * aos_array(i) 
         grad_ao_matrix_dvc_a(i,l,m,istate) = contrib_grad_ca(m,istate) * grad_aos_array_transpose(i,m) 
         grad_ao_matrix_dvc_b(i,l,m,istate) = contrib_grad_cb(m,istate) * grad_aos_array_transpose(i,m) 
         grad_ao_matrix_dvx_a(i,l,m,istate) = contrib_grad_xa(m,istate) * grad_aos_array_transpose(i,m) 
         grad_ao_matrix_dvx_b(i,l,m,istate) = contrib_grad_xb(m,istate) * grad_aos_array_transpose(i,m) 
        enddo
       enddo
      enddo
     enddo
     call update_pots_scalar_dgemm(ao_matrix_vc_a,ao_matrix_vx_a,ao_matrix_vc_b,ao_matrix_vx_b,ao_matrix,potential_c_alpha_ao_new,potential_x_alpha_ao_new,potential_c_beta_ao_new,potential_x_beta_ao_new)
     call update_pots_gradient_dgemm(ao_matrix,ao_matrix_dvc_a,ao_matrix_dvc_b,ao_matrix_dvx_a,ao_matrix_dvx_b,grad_ao_matrix,grad_ao_matrix_dvc_a,grad_ao_matrix_dvc_b,grad_ao_matrix_dvx_a,grad_ao_matrix_dvx_b,potential_c_alpha_ao_new,potential_x_alpha_ao_new,potential_c_beta_ao_new,potential_x_beta_ao_new )
    endif
    deallocate(tmp_x_a,tmp_x_b,tmp_c_a,tmp_c_b,aos_array,r)
   enddo
  enddo

END_PROVIDER 


 subroutine update_pots_scalar_dgemm(ao_matrix_vc_a,ao_matrix_vx_a,ao_matrix_vc_b,ao_matrix_vx_b,ao_matrix,pot_c_alpha_ao_new,pot_x_alpha_ao_new,pot_c_beta_ao_new,pot_x_beta_ao_new)
  implicit none
  double precision, intent(in) :: ao_matrix(ao_num,n_points_integration_angular), ao_matrix_vc_a(ao_num,n_points_integration_angular,N_states), ao_matrix_vc_b(ao_num,n_points_integration_angular,N_states), ao_matrix_vx_a(ao_num,n_points_integration_angular,N_states), ao_matrix_vx_b(ao_num,n_points_integration_angular,N_states)
  double precision, intent(inout):: pot_c_alpha_ao_new(ao_num,ao_num,N_states)
  double precision, intent(inout):: pot_x_alpha_ao_new(ao_num,ao_num,N_states)
  double precision, intent(inout):: pot_c_beta_ao_new(ao_num,ao_num,N_states)
  double precision, intent(inout):: pot_x_beta_ao_new(ao_num,ao_num,N_states)
  integer :: istate
  do istate = 1, N_states
   call dgemm('N','T',ao_num,ao_num,n_points_integration_angular,1.d0,ao_matrix_vc_a(1,1,istate),size(ao_matrix_vc_a,1),ao_matrix,size(ao_matrix,1),1.d0,pot_c_alpha_ao_new(1,1,istate),size(pot_c_alpha_ao_new,1))
   call dgemm('N','T',ao_num,ao_num,n_points_integration_angular,1.d0,ao_matrix_vx_a(1,1,istate),size(ao_matrix_vx_a,1),ao_matrix,size(ao_matrix,1),1.d0,pot_x_alpha_ao_new(1,1,istate),size(pot_x_alpha_ao_new,1))
   call dgemm('N','T',ao_num,ao_num,n_points_integration_angular,1.d0,ao_matrix_vc_b(1,1,istate),size(ao_matrix_vc_b,1),ao_matrix,size(ao_matrix,1),1.d0,pot_c_beta_ao_new(1,1,istate),size(pot_c_beta_ao_new,1))
   call dgemm('N','T',ao_num,ao_num,n_points_integration_angular,1.d0,ao_matrix_vx_b(1,1,istate),size(ao_matrix_vx_b,1),ao_matrix,size(ao_matrix,1),1.d0,pot_x_beta_ao_new(1,1,istate),size(pot_x_beta_ao_new,1))
  enddo
 end
 
 subroutine update_pots_gradient_dgemm(ao_matrix,ao_matrix_dvc_a,ao_matrix_dvc_b,ao_matrix_dvx_a,ao_matrix_dvx_b,grad_ao_matrix,grad_ao_matrix_dvc_a,grad_ao_matrix_dvc_b,grad_ao_matrix_dvx_a,grad_ao_matrix_dvx_b,pot_c_alpha_ao_new,pot_x_alpha_ao_new,pot_c_beta_ao_new,pot_x_beta_ao_new )
  implicit none
  double precision, intent(in) :: ao_matrix(ao_num,n_points_integration_angular)
  double precision, intent(in) :: ao_matrix_dvc_a(ao_num,n_points_integration_angular,3,N_states)
  double precision, intent(in) :: ao_matrix_dvc_b(ao_num,n_points_integration_angular,3,N_states)
  double precision, intent(in) :: ao_matrix_dvx_a(ao_num,n_points_integration_angular,3,N_states)
  double precision, intent(in) :: ao_matrix_dvx_b(ao_num,n_points_integration_angular,3,N_states)
  double precision, intent(in) :: grad_ao_matrix(ao_num,n_points_integration_angular,3) 
  double precision, intent(in) :: grad_ao_matrix_dvc_a(ao_num,n_points_integration_angular,3,N_states)
  double precision, intent(in) :: grad_ao_matrix_dvc_b(ao_num,n_points_integration_angular,3,N_states)
  double precision, intent(in) :: grad_ao_matrix_dvx_a(ao_num,n_points_integration_angular,3,N_states)
  double precision, intent(in) :: grad_ao_matrix_dvx_b(ao_num,n_points_integration_angular,3,N_states)
  double precision, intent(inout):: pot_c_alpha_ao_new(ao_num,ao_num,N_states)
  double precision, intent(inout):: pot_x_alpha_ao_new(ao_num,ao_num,N_states)
  double precision, intent(inout):: pot_c_beta_ao_new(ao_num,ao_num,N_states)
  double precision, intent(inout):: pot_x_beta_ao_new(ao_num,ao_num,N_states)
  integer :: istate,m
  do istate = 1, N_states
   do m= 1,3
    call dgemm('N','T',ao_num,ao_num,n_points_integration_angular,1.d0,ao_matrix_dvc_a(1,1,m,istate),size(ao_matrix_dvc_a,1),grad_ao_matrix(1,1,m),size(grad_ao_matrix,1),1.d0,pot_c_alpha_ao_new(1,1,istate),size(pot_c_alpha_ao_new,1))
    call dgemm('N','T',ao_num,ao_num,n_points_integration_angular,1.d0,ao_matrix_dvx_a(1,1,m,istate),size(ao_matrix_dvx_a,1),grad_ao_matrix(1,1,m),size(grad_ao_matrix,1),1.d0,pot_x_alpha_ao_new(1,1,istate),size(pot_x_alpha_ao_new,1))
    call dgemm('N','T',ao_num,ao_num,n_points_integration_angular,1.d0,ao_matrix_dvc_b(1,1,m,istate),size(ao_matrix_dvc_b,1),grad_ao_matrix(1,1,m),size(grad_ao_matrix,1),1.d0,pot_c_beta_ao_new(1,1,istate),size(pot_c_beta_ao_new,1))
    call dgemm('N','T',ao_num,ao_num,n_points_integration_angular,1.d0,ao_matrix_dvx_b(1,1,m,istate),size(ao_matrix_dvx_b,1),grad_ao_matrix(1,1,m),size(grad_ao_matrix,1),1.d0,pot_x_beta_ao_new(1,1,istate),size(pot_x_beta_ao_new,1))
    call dgemm('N','T',ao_num,ao_num,n_points_integration_angular,1.d0,grad_ao_matrix_dvc_a(1,1,m,istate),size(grad_ao_matrix_dvc_a,1),ao_matrix,size(ao_matrix,1),1.d0,pot_c_alpha_ao_new(1,1,istate),size(pot_c_alpha_ao_new,1))
    call dgemm('N','T',ao_num,ao_num,n_points_integration_angular,1.d0,grad_ao_matrix_dvx_a(1,1,m,istate),size(grad_ao_matrix_dvx_a,1),ao_matrix,size(ao_matrix,1),1.d0,pot_x_alpha_ao_new(1,1,istate),size(pot_x_alpha_ao_new,1))
    call dgemm('N','T',ao_num,ao_num,n_points_integration_angular,1.d0,grad_ao_matrix_dvc_b(1,1,m,istate),size(grad_ao_matrix_dvc_b,1),ao_matrix,size(ao_matrix,1),1.d0,pot_c_beta_ao_new(1,1,istate),size(pot_c_beta_ao_new,1))
    call dgemm('N','T',ao_num,ao_num,n_points_integration_angular,1.d0,grad_ao_matrix_dvx_b(1,1,m,istate),size(grad_ao_matrix_dvx_b,1),ao_matrix,size(ao_matrix,1),1.d0,pot_x_beta_ao_new(1,1,istate),size(pot_x_beta_ao_new,1))
   enddo
  enddo
 end

 subroutine create_matrix_for_LDA_like_potential_update_energy(k_radial,j_nucl,ao_matrix,ao_matrix_vc_a,ao_matrix_vc_b,ao_matrix_vx_a,ao_matrix_vx_b,ex_tot,ec_tot)
 implicit none
 integer, intent(in) :: j_nucl,k_radial,l_angular
 double precision, intent(inout) :: ex_tot(N_states), ec_tot(N_states) 
 double precision, intent(inout) :: ao_matrix(ao_num,n_points_integration_angular)
 double precision, intent(inout) :: ao_matrix_vc_a(ao_num,n_points_integration_angular,N_states)
 double precision, intent(inout) :: ao_matrix_vc_b(ao_num,n_points_integration_angular,N_states)
 double precision, intent(inout) :: ao_matrix_vx_a(ao_num,n_points_integration_angular,N_states)
 double precision, intent(inout) :: ao_matrix_vx_b(ao_num,n_points_integration_angular,N_states)
 double precision :: weight,rho_a(N_states), rho_b(N_states), ex(N_states), ec(N_states), vc_rho_a(N_states), vc_rho_b(N_states), vx_rho_a(N_states), vx_rho_b(N_states), aos_array(ao_num), r(3)
 integer :: i,istate

  r(1) = grid_points_per_atom(1,l_angular,k_radial,j_nucl)
  r(2) = grid_points_per_atom(2,l_angular,k_radial,j_nucl)
  r(3) = grid_points_per_atom(3,l_angular,k_radial,j_nucl)
  weight = final_weight_functions_at_grid_points(l_angular,k_radial,j_nucl)
  call dm_dft_alpha_beta_and_all_aos_at_r(r,rho_a,rho_b,aos_array)
  call LDA_type_functional(rho_a,rho_b,vx_rho_a,vx_rho_b,vc_rho_a,vc_rho_b,ex,ec)

  do i = 1 , ao_num 
   ao_matrix(i,l_angular) = aos_array(i)
  enddo

  do istate = 1, N_states
   ex_tot(istate) += weight *  ex(istate) 
   ec_tot(istate) += weight *  ec(istate)
   vx_rho_a(istate) *= weight
   vc_rho_a(istate) *= weight
   vx_rho_b(istate) *= weight
   vc_rho_b(istate) *= weight
   do i = 1 , ao_num 
    ao_matrix_vc_a(i,l_angular,istate) = vc_rho_a(istate) * aos_array(i)
    ao_matrix_vc_b(i,l_angular,istate) = vc_rho_b(istate) * aos_array(i)
    ao_matrix_vx_a(i,l_angular,istate) = vx_rho_a(istate) * aos_array(i)
    ao_matrix_vx_b(i,l_angular,istate) = vx_rho_b(istate) * aos_array(i)
   enddo
  enddo
 end


