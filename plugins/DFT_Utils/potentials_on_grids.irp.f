 BEGIN_PROVIDER[double precision, aos_vc_alpha_lda_w, (n_points_final_grid,ao_num,N_states)]
&BEGIN_PROVIDER[double precision, aos_vc_beta_lda_w, (n_points_final_grid,ao_num,N_states)]
&BEGIN_PROVIDER[double precision, aos_vx_alpha_lda_w, (n_points_final_grid,ao_num,N_states)]
&BEGIN_PROVIDER[double precision, aos_vx_beta_lda_w, (n_points_final_grid,ao_num,N_states)]
&BEGIN_PROVIDER[double precision, energy_x_LDA, (N_states) ]
&BEGIN_PROVIDER[double precision, energy_c_LDA, (N_states) ]
 implicit none
 BEGIN_DOC
! aos_vxc_alpha_lda_w(j,i) = ao_i(r_j) * (v^x_alpha(r_j) + v^c_alpha(r_j)) * W(r_j)
 END_DOC
 integer :: istate,i,j
 double precision :: r(3)
 double precision :: mu,weight
 double precision :: e_c,vc_a,vc_b,e_x,vx_a,vx_b
 double precision, allocatable :: rhoa(:),rhob(:)
 allocate(rhoa(N_states), rhob(N_states))
 energy_x_LDA = 0.d0
 energy_c_LDA = 0.d0
 do istate = 1, N_states
  do i = 1, n_points_final_grid
   r(1) = final_grid_points(1,i)
   r(2) = final_grid_points(2,i)
   r(3) = final_grid_points(3,i)
   weight=final_weight_functions_at_final_grid_points(i)
   call dm_dft_alpha_beta_at_r(r,rhoa(istate),rhob(istate))
   call ec_lda_sr(mu_erf,rhoa(istate),rhob(istate),e_c,vc_a,vc_b)
   call ex_lda_sr(mu_erf,rhoa(istate),rhob(istate),e_x,vx_a,vx_b)
   energy_x_LDA(istate) += weight * e_x
   energy_c_LDA(istate) += weight * e_c
   do j =1, ao_num
    aos_vc_alpha_lda_w(i,j,istate) = vc_a * aos_in_r_array(j,i)*weight
    aos_vc_beta_lda_w(i,j,istate)  = vc_b * aos_in_r_array(j,i)*weight
    aos_vx_alpha_lda_w(i,j,istate) = vx_a * aos_in_r_array(j,i)*weight
    aos_vx_beta_lda_w(i,j,istate)  = vx_b * aos_in_r_array(j,i)*weight
   enddo
  enddo
 enddo

 END_PROVIDER 


 BEGIN_PROVIDER [double precision, potential_x_alpha_ao_LDA,(ao_num,ao_num,N_states)]
&BEGIN_PROVIDER [double precision, potential_x_beta_ao_LDA,(ao_num,ao_num,N_states)]
&BEGIN_PROVIDER [double precision, potential_c_alpha_ao_LDA,(ao_num,ao_num,N_states)]
&BEGIN_PROVIDER [double precision, potential_c_beta_ao_LDA,(ao_num,ao_num,N_states)]
 implicit none
 integer :: istate
 double precision :: wall_1,wall_2
 call wall_time(wall_1)
 do istate = 1, N_states 
  call dgemm('N','N',ao_num,ao_num,n_points_final_grid,1.d0,aos_in_r_array,ao_num,aos_vc_alpha_lda_w(1,1,istate),n_points_final_grid,0.d0,potential_c_alpha_ao_LDA(1,1,istate),ao_num)
  call dgemm('N','N',ao_num,ao_num,n_points_final_grid,1.d0,aos_in_r_array,ao_num,aos_vc_beta_lda_w(1,1,istate) ,n_points_final_grid,0.d0,potential_c_beta_ao_LDA(1,1,istate),ao_num)
  call dgemm('N','N',ao_num,ao_num,n_points_final_grid,1.d0,aos_in_r_array,ao_num,aos_vx_alpha_lda_w(1,1,istate),n_points_final_grid,0.d0,potential_x_alpha_ao_LDA(1,1,istate),ao_num)
  call dgemm('N','N',ao_num,ao_num,n_points_final_grid,1.d0,aos_in_r_array,ao_num,aos_vx_beta_lda_w(1,1,istate) ,n_points_final_grid,0.d0,potential_x_beta_ao_LDA(1,1,istate),ao_num)
 enddo
 call wall_time(wall_2)
 print*,'time to provide potential_x/c_alpha/beta_ao_LDA = ',wall_2 - wall_1

 END_PROVIDER 

