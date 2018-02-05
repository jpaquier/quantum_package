
 BEGIN_PROVIDER [double precision, energy_x, (N_states)]
&BEGIN_PROVIDER [double precision, energy_c, (N_states)]
&BEGIN_PROVIDER [double precision, potential_x_alpha_ao,(ao_num,ao_num,N_states)]
&BEGIN_PROVIDER [double precision, potential_x_beta_ao,(ao_num,ao_num,N_states)]
&BEGIN_PROVIDER [double precision, potential_c_alpha_ao,(ao_num,ao_num,N_states)]
&BEGIN_PROVIDER [double precision, potential_c_beta_ao,(ao_num,ao_num,N_states)]

 implicit none
 integer :: j,k,l,istate
 integer :: m,n
 double precision, allocatable :: aos_array(:)
 double precision, allocatable :: r(:)
 double precision :: rho_a(N_states),rho_b(N_states),ex(N_states),ec(N_states)
 double precision :: vx_a(N_states),vx_b(N_states),vc_a(N_states),vc_b(N_states)
 potential_c_alpha_ao = 0.d0
 potential_c_beta_ao = 0.d0
 potential_x_alpha_ao = 0.d0
 potential_x_beta_ao = 0.d0
 double precision, allocatable :: tmp_c_a(:,:,:),tmp_c_b(:,:,:),tmp_x_a(:,:,:),tmp_x_b(:,:,:)
 print*,'providing the potentials ...'
  energy_x = 0.d0
  energy_c = 0.d0
  do j = 1, nucl_num
   do k = 1, n_points_radial_grid  -1
   !!$OMP PARALLEL DEFAULT(NONE)                                                                                                                   &
   !!OMP PRIVATE(r,l,rho_a,rho_b,ex,ec,vx_a,vx_b,vc_a,vc_b,aos_array,contrib_xa,contrib_xb,contrib_ca, contrib_cb,tmp_c_a,tmp_c_b,tmp_x_a,tmp_x_b) & 
   !!OMP SHARED(i,j,k,ao_num,n_points_integration_angular,exchange_functional,correlation_functional,final_weight_functions_at_grid_points,potential_x_alpha_ao,potential_x_beta_ao,potential_c_alpha_ao,potential_c_beta_ao)                &
   !!$OMP REDUCTION (+:energy_x)       &
   !!$OMP REDUCTION (+:energy_c)        
    allocate(tmp_c_a(ao_num,ao_num,N_states),tmp_c_b(ao_num,ao_num,N_states),tmp_x_a(ao_num,ao_num,N_states),tmp_x_b(ao_num,ao_num,N_states),aos_array(ao_num),r(3))
    
    tmp_c_a = 0.d0
    tmp_c_b = 0.d0
    tmp_x_a = 0.d0
    tmp_x_b = 0.d0
   !!$OMP DO SCHEDULE(static)
    do l = 1, n_points_integration_angular 
     
     r(1) = grid_points_per_atom(1,l,k,j)
     r(2) = grid_points_per_atom(2,l,k,j)
     r(3) = grid_points_per_atom(3,l,k,j)
     call dm_dft_alpha_beta_and_all_aos_at_r(r,rho_a,rho_b,aos_array)

!!!!!!!!!!! EXCHANGE PART
     do istate = 1, N_states
      if(exchange_functional.EQ."short_range_LDA")then
       call ex_lda_sr(rho_a(istate),rho_b(istate),ex(istate),vx_a(istate),vx_b(istate))
      else if(exchange_functional.EQ."LDA")then
       call ex_lda(rho_a(istate),rho_b(istate),ex(istate),vx_a(istate),vx_b(istate))
      else if(exchange_functional.EQ."None")then
       ex = 0.d0
       vx_a = 0.d0
       vx_b = 0.d0
      else
       print*, 'Exchange functional required does not exist ...'
       print*,'exchange_functional',exchange_functional
       stop
      endif

!!!!!!!!!!!! CORRELATION PART
      if(correlation_functional.EQ."short_range_LDA")then
       call ec_lda_sr(rho_a(istate),rho_b(istate),ec(istate),vc_a(istate),vc_b(istate))
      else if(correlation_functional.EQ."LDA")then
       call ec_lda(rho_a(istate),rho_b(istate),ec(istate),vc_a(istate),vc_b(istate))
      else if(correlation_functional.EQ."None")then
       ec = 0.d0
       vc_a = 0.d0
       vc_b = 0.d0
      else
       print*, 'Correlation functional required does not exist ...'
       print*, 'correlation_functional',correlation_functional
       stop
      endif
     enddo
     double precision :: contrib_xa,contrib_xb,contrib_ca, contrib_cb
     do istate = 1, N_states
      energy_x(istate) += final_weight_functions_at_grid_points(l,k,j) * ex(istate) 
      energy_c(istate) += final_weight_functions_at_grid_points(l,k,j) * ec(istate)

      contrib_xa =  vx_a(istate) * final_weight_functions_at_grid_points(l,k,j)
      contrib_ca =  vc_a(istate) * final_weight_functions_at_grid_points(l,k,j)
      contrib_xb =  vx_b(istate) * final_weight_functions_at_grid_points(l,k,j)
      contrib_cb =  vc_b(istate) * final_weight_functions_at_grid_points(l,k,j)
      call dger(ao_num,ao_num,contrib_xa,aos_array,1,aos_array,1,tmp_x_a,size(tmp_x_a,1))
      call dger(ao_num,ao_num,contrib_ca,aos_array,1,aos_array,1,tmp_c_a,size(tmp_c_a,1))
      call dger(ao_num,ao_num,contrib_xb,aos_array,1,aos_array,1,tmp_x_b,size(tmp_x_b,1))
      call dger(ao_num,ao_num,contrib_cb,aos_array,1,aos_array,1,tmp_c_b,size(tmp_c_b,1))
     enddo

    enddo
  !!$OMP END DO 
    
  !!$OMP CRITICAL
   
    do istate = 1,N_states
     potential_c_alpha_ao(:,:,istate) = potential_c_alpha_ao(:,:,istate) + tmp_c_a(:,:,istate)
     potential_x_alpha_ao(:,:,istate) = potential_x_alpha_ao(:,:,istate) + tmp_x_a(:,:,istate)
     potential_c_beta_ao(:,:,istate)  =  potential_c_beta_ao(:,:,istate) + tmp_c_b(:,:,istate)
     potential_x_beta_ao(:,:,istate)  =  potential_x_beta_ao(:,:,istate) + tmp_x_b(:,:,istate)
    enddo
  !!$OMP END CRITICAL
    deallocate(tmp_x_a,tmp_x_b,tmp_c_a,tmp_c_b,aos_array,r)
  !!$OMP END PARALLEL
   enddo
  enddo
 print*,'potentials provided !' 

END_PROVIDER 

 BEGIN_PROVIDER [double precision, potential_x_alpha_mo,(mo_tot_num,mo_tot_num,N_states)]
&BEGIN_PROVIDER [double precision, potential_x_beta_mo,(mo_tot_num,mo_tot_num,N_states)]
&BEGIN_PROVIDER [double precision, potential_c_alpha_mo,(mo_tot_num,mo_tot_num,N_states)]
&BEGIN_PROVIDER [double precision, potential_c_beta_mo,(mo_tot_num,mo_tot_num,N_states)]
 implicit none
 integer :: istate 
 do istate = 1, N_states
    call ao_to_mo(                                                   &
        potential_x_alpha_ao(1,1,istate),                                 &
        size(potential_x_alpha_ao,1),                                &
        potential_x_alpha_mo(1,1,istate),                                 &
        size(potential_x_alpha_mo,1)                                 &
        )

    call ao_to_mo(                                                   &
        potential_x_beta_ao(1,1,istate),                                  &
        size(potential_x_beta_ao,1),                                 &
        potential_x_beta_mo(1,1,istate),                                  &
        size(potential_x_beta_mo,1)                                  &
        )


    call ao_to_mo(                                                   &
        potential_c_alpha_ao(1,1,istate),                                 &
        size(potential_c_alpha_ao,1),                                &
        potential_c_alpha_mo(1,1,istate),                                 &
        size(potential_c_alpha_mo,1)                                 &
        )

    call ao_to_mo(                                                   &
        potential_c_beta_ao(1,1,istate),                                  &
        size(potential_c_beta_ao,1),                                 &
        potential_c_beta_mo(1,1,istate),                                  &
        size(potential_c_beta_mo,1)                                  &
        )

 enddo

END_PROVIDER 
