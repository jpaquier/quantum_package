 BEGIN_PROVIDER [double precision, energy_x_new, (N_states)]
&BEGIN_PROVIDER [double precision, energy_c_new, (N_states)]
&BEGIN_PROVIDER [double precision, potential_x_alpha_ao_new,(ao_num,ao_num,N_states)]
&BEGIN_PROVIDER [double precision, potential_x_beta_ao_new,(ao_num,ao_num,N_states)]
&BEGIN_PROVIDER [double precision, potential_c_alpha_ao_new,(ao_num,ao_num,N_states)]
&BEGIN_PROVIDER [double precision, potential_c_beta_ao_new,(ao_num,ao_num,N_states)]

 implicit none
 integer :: j,k,l,istate
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
  energy_x_new = 0.d0
  energy_c_new = 0.d0
  potential_c_alpha_ao_new = 0.d0
  potential_c_beta_ao_new = 0.d0
  potential_x_alpha_ao_new_new = 0.d0
  potential_x_beta_ao_new = 0.d0
  do j = 1, nucl_num
   do k = 1, n_points_radial_grid  -1
    allocate(tmp_c_a(ao_num,ao_num,N_states),tmp_c_b(ao_num,ao_num,N_states),tmp_x_a(ao_num,ao_num,N_states),tmp_x_b(ao_num,ao_num,N_states),aos_array(ao_num),r(3))
    tmp_c_a = 0.d0
    tmp_c_b = 0.d0
    tmp_x_a = 0.d0
    tmp_x_b = 0.d0
    do l = 1, n_points_integration_angular 
     r(1) = grid_points_per_atom(1,l,k,j)
     r(2) = grid_points_per_atom(2,l,k,j)
     r(3) = grid_points_per_atom(3,l,k,j)
     weight = final_weight_functions_at_grid_points(l,k,j)
     if(DFT_TYPE=="LDA")then
      call dm_dft_alpha_beta_and_all_aos_at_r(r,rho_a,rho_b,aos_array)
      call LDA_type_functional(rho_a,rho_b,vx_rho_a,vx_rho_b,vc_rho_a,vc_rho_b,ex,ec)
      do istate = 1, N_states
       energy_x_new(istate) += weight *  ex(istate) 
       energy_c_new(istate) += weight *  ec(istate)
       vx_rho_a(istate) *= weight
       vc_rho_a(istate) *= weight
       vx_rho_b(istate) *= weight
       vc_rho_b(istate) *= weight
      enddo
      call update_potentials_scalar_dger(tmp_c_a,tmp_c_b,tmp_x_a,tmp_x_b,vc_rho_a,vc_rho_b,vx_rho_a,vx_rho_b,aos_array)

     else if (DFT_TYPE=="GGA")then
      call density_and_grad_alpha_beta_and_all_aos_and_grad_aos_at_r(r,rho_a,rho_b, grad_rho_a, grad_rho_b, aos_array, grad_aos_array)
      grad_rho_a_2 = 0.d0
      grad_rho_b_2 = 0.d0
      grad_rho_a_b = 0.d0
      do istate = 1, N_states
       do m = 1, 3
        grad_rho_a_2(istate) += grad_rho_a(m,istate)*grad_rho_a(m,istate)
        grad_rho_b_2(istate) += grad_rho_b(m,istate)*grad_rho_b(m,istate)
        grad_rho_a_b(istate) += grad_rho_a(m,istate)*grad_rho_b(m,istate)
       enddo
      enddo
      call GGA_type_functionals(rho_a,rho_b,grad_rho_a_2,grad_rho_b_2,grad_rho_a_b,ex,vx_rho_a,vx_rho_b,vx_grad_rho_a_2,vx_grad_rho_b_2,vx_grad_rho_a_b, &  
                                                                                   ec,vc_rho_a,vc_rho_b,vc_grad_rho_a_2,vc_grad_rho_b_2,vc_grad_rho_a_b )
      do istate = 1, N_states
       vx_rho_a(istate) *= weight
       vc_rho_a(istate) *= weight
       vx_rho_b(istate) *= weight
       vc_rho_b(istate) *= weight
       energy_x_new(istate) += weight *  ex(istate) 
       energy_c_new(istate) += weight *  ec(istate)
       do m = 1, 3
        dtmp_x_a(m,istate) = (2.d0 * vx_grad_rho_a_2(istate) *  grad_rho_a(m,istate) + vx_grad_rho_a_b(istate)  * grad_rho_b(m,istate)) * weight
        dtmp_x_b(m,istate) = (2.d0 * vx_grad_rho_b_2(istate) *  grad_rho_b(m,istate) + vx_grad_rho_a_b(istate)  * grad_rho_a(m,istate)) * weight
       enddo
      enddo
      call update_potentials_scalar_dger(tmp_c_a,tmp_c_b,tmp_x_a,tmp_x_b,vc_rho_a,vc_rho_b,vx_rho_a,vx_rho_b,aos_array)
      call dtranspose(grad_aos_array,3,grad_aos_array_transpose,ao_num,3,ao_num)
      call update_potentials_gradient_dger(dtmp_x_a,dtmp_x_b,dtmp_c_a,dtmp_c_b,aos_array,grad_aos_array_transpose,tmp_x_a,tmp_x_b,tmp_c_a,tmp_c_b)
     endif
    enddo ! angular points 

    do istate = 1,N_states
     potential_c_alpha_ao_new_new(:,:,istate) = potential_c_alpha_ao_new_new(:,:,istate) + tmp_c_a(:,:,istate)
     potential_x_alpha_ao_new_new(:,:,istate) = potential_x_alpha_ao_new_new(:,:,istate) + tmp_x_a(:,:,istate)
     potential_c_beta_ao_new(:,:,istate)  =  potential_c_beta_ao_new(:,:,istate) + tmp_c_b(:,:,istate)
     potential_x_beta_ao_new(:,:,istate)  =  potential_x_beta_ao_new(:,:,istate) + tmp_x_b(:,:,istate)
    enddo
    deallocate(tmp_x_a,tmp_x_b,tmp_c_a,tmp_c_b,aos_array,r)
   enddo
  enddo

END_PROVIDER 

