subroutine ex_lda(rho_a,rho_b,ex,vx_a,vx_b)
 include 'constants.include.F'
 implicit none
 double precision, intent(in) :: rho_a,rho_b
 double precision, intent(out) :: ex,vx_a,vx_b
 double precision :: tmp_a,tmp_b
 tmp_a = rho_a**(c_1_3)
 tmp_b = rho_b**(c_1_3)
 ex = cst_lda * (tmp_a*tmp_a*tmp_a*tmp_a + tmp_b*tmp_b*tmp_b*tmp_b)
 vx_a = cst_lda * c_4_3 * tmp_a
 vx_b = cst_lda * c_4_3 * tmp_b

end

 BEGIN_PROVIDER [double precision, lda_exchange, (N_states)]
&BEGIN_PROVIDER [double precision, lda_ex_potential_alpha_ao,(ao_num,ao_num,N_states)]
&BEGIN_PROVIDER [double precision, lda_ex_potential_beta_ao,(ao_num,ao_num,N_states)]

 implicit none
 integer :: i,j,k,l
 integer :: m,n
 double precision :: aos_array(ao_num)
 double precision :: r(3)
 lda_ex_potential_alpha_ao = 0.d0
 lda_ex_potential_beta_ao = 0.d0
 do l = 1, N_states
  lda_exchange(l) = 0.d0
  do j = 1, nucl_num
   do i = 1, n_points_radial_grid  -1 
    do k = 1, n_points_integration_angular 
     double precision :: rho_a,rho_b,ex
     double precision :: vx_a,vx_b
     rho_a = one_body_dm_mo_alpha_at_grid_points(k,i,j,l)
     rho_b = one_body_dm_mo_beta_at_grid_points(k,i,j,l)
     call ex_lda(rho_a,rho_b,ex,vx_a,vx_b) 
     lda_exchange(l) += final_weight_functions_at_grid_points(k,i,j) * ex
     r(1) = grid_points_per_atom(1,k,i,j) 
     r(2) = grid_points_per_atom(2,k,i,j) 
     r(3) = grid_points_per_atom(3,k,i,j) 
     call give_all_aos_at_r(r,aos_array)
     do m = 1, ao_num
      do n = 1, ao_num
       lda_ex_potential_alpha_ao(m,n,l) += (vx_a ) * aos_array(m)*aos_array(n) * final_weight_functions_at_grid_points(k,i,j)
       lda_ex_potential_beta_ao(m,n,l) += (vx_b) * aos_array(m)*aos_array(n) * final_weight_functions_at_grid_points(k,i,j)
      enddo
     enddo
    enddo
   enddo
  enddo
 enddo

END_PROVIDER 


 BEGIN_PROVIDER [double precision, energy_x, (N_states)]
&BEGIN_PROVIDER [double precision, energy_c, (N_states)]
&BEGIN_PROVIDER [double precision, potential_x_alpha_ao,(ao_num,ao_num,N_states)]
&BEGIN_PROVIDER [double precision, potential_x_beta_ao,(ao_num,ao_num,N_states)]
&BEGIN_PROVIDER [double precision, potential_c_alpha_ao,(ao_num,ao_num,N_states)]
&BEGIN_PROVIDER [double precision, potential_c_beta_ao,(ao_num,ao_num,N_states)]

 implicit none
 integer :: i,j,k,l
 integer :: m,n
 double precision, allocatable :: aos_array(:)
 double precision, allocatable :: r(:)
 double precision :: rho_a,rho_b,ex,ec
 double precision :: vx_a,vx_b,vc_a,vc_b
 potential_c_alpha_ao = 0.d0
 potential_c_beta_ao = 0.d0
 potential_x_alpha_ao = 0.d0
 potential_x_beta_ao = 0.d0
 double precision, allocatable :: tmp_c_a(:,:),tmp_c_b(:,:),tmp_x_a(:,:),tmp_x_b(:,:)
 print*,'providing the potentials ...'
 do i = 1, N_states
  energy_x(i) = 0.d0
  energy_c(i) = 0.d0
  do j = 1, nucl_num
   do k = 1, n_points_radial_grid  -1
   !!$OMP PARALLEL DEFAULT(NONE)                                                                                                                   &
   !!OMP PRIVATE(r,l,rho_a,rho_b,ex,ec,vx_a,vx_b,vc_a,vc_b,aos_array,contrib_xa,contrib_xb,contrib_ca, contrib_cb,tmp_c_a,tmp_c_b,tmp_x_a,tmp_x_b) & 
   !!OMP SHARED(i,j,k,ao_num,n_points_integration_angular,exchange_functional,correlation_functional,final_weight_functions_at_grid_points,potential_x_alpha_ao,potential_x_beta_ao,potential_c_alpha_ao,potential_c_beta_ao)                &
   !!$OMP REDUCTION (+:energy_x)       &
   !!$OMP REDUCTION (+:energy_c)        
    allocate(tmp_c_a(ao_num,ao_num),tmp_c_b(ao_num,ao_num),tmp_x_a(ao_num,ao_num),tmp_x_b(ao_num,ao_num),aos_array(ao_num),r(3))
    
    tmp_c_a = 0.d0
    tmp_c_b = 0.d0
    tmp_x_a = 0.d0
    tmp_x_b = 0.d0
   !!$OMP DO SCHEDULE(static)
    do l = 1, n_points_integration_angular 
     
     r(1) = grid_points_per_atom(1,l,k,j)
     r(2) = grid_points_per_atom(2,l,k,j)
     r(3) = grid_points_per_atom(3,l,k,j)
     call density_matrices_alpha_beta_and_all_aos_at_r(r,rho_a,rho_b,aos_array)

!!!!!!!!!!! EXCHANGE PART
     if(exchange_functional.EQ."short_range_LDA")then
      call ex_lda_sr(rho_a,rho_b,ex,vx_a,vx_b)
     else if(exchange_functional.EQ."LDA")then
      call ex_lda(rho_a,rho_b,ex,vx_a,vx_b) 
     else if(exchange_functional.EQ."None")then
      ex = 0.d0
      vx_a = 0.d0
      vx_b = 0.d0
     else
      print*, 'Exchange functional required does not exist ...'
      print*,'exchange_functional',exchange_functional
      stop
     endif

!!!!!!!!!!! CORRELATION PART
     if(correlation_functional.EQ."short_range_LDA")then
      call ec_lda_sr(rho_a,rho_b,ec,vc_a,vc_b)
     else if(correlation_functional.EQ."LDA")then
      call ec_lda(rho_a,rho_b,ec,vc_a,vc_b)
     else if(correlation_functional.EQ."None")then
      ec = 0.d0
      vc_a = 0.d0
      vc_b = 0.d0
     else
      print*, 'Correlation functional required does not exist ...'
      print*, 'correlation_functional',correlation_functional
      stop
     endif
     energy_x(i) += final_weight_functions_at_grid_points(l,k,j) * ex 
     energy_c(i) += final_weight_functions_at_grid_points(l,k,j) * ec

     double precision :: contrib_xa,contrib_xb,contrib_ca, contrib_cb
     contrib_xa =  vx_a * final_weight_functions_at_grid_points(l,k,j)
     contrib_ca =  vc_a * final_weight_functions_at_grid_points(l,k,j)
     contrib_xb =  vx_b * final_weight_functions_at_grid_points(l,k,j)
     contrib_cb =  vc_b * final_weight_functions_at_grid_points(l,k,j)
     call dger(ao_num,ao_num,contrib_xa,aos_array,1,aos_array,1,tmp_x_a,size(tmp_x_a,1))
     call dger(ao_num,ao_num,contrib_ca,aos_array,1,aos_array,1,tmp_c_a,size(tmp_c_a,1))
     call dger(ao_num,ao_num,contrib_xb,aos_array,1,aos_array,1,tmp_x_b,size(tmp_x_b,1))
     call dger(ao_num,ao_num,contrib_cb,aos_array,1,aos_array,1,tmp_c_b,size(tmp_c_b,1))

    enddo
  !!$OMP END DO 
    
  !!$OMP CRITICAL
    potential_c_alpha_ao(:,:,i) = potential_c_alpha_ao(:,:,i) + tmp_c_a(:,:)
    potential_x_alpha_ao(:,:,i) = potential_x_alpha_ao(:,:,i) + tmp_x_a(:,:)
    potential_c_beta_ao(:,:,i)  =  potential_c_beta_ao(:,:,i) + tmp_c_b(:,:)
    potential_x_beta_ao(:,:,i)  =  potential_x_beta_ao(:,:,i) + tmp_x_b(:,:)
  !!$OMP END CRITICAL
    deallocate(tmp_x_a,tmp_x_b,tmp_c_a,tmp_c_b,aos_array,r)
  !!$OMP END PARALLEL
   enddo
  enddo
 enddo
 print*,'potentials provided !' 

END_PROVIDER 

 BEGIN_PROVIDER [double precision, potential_x_alpha_mo,(mo_tot_num,mo_tot_num,N_states)]
&BEGIN_PROVIDER [double precision, potential_x_beta_mo,(mo_tot_num,mo_tot_num,N_states)]
&BEGIN_PROVIDER [double precision, potential_c_alpha_mo,(mo_tot_num,mo_tot_num,N_states)]
&BEGIN_PROVIDER [double precision, potential_c_beta_mo,(mo_tot_num,mo_tot_num,N_states)]
 implicit none
 
    call ao_to_mo(                                                   &
        potential_x_alpha_ao(1,1,1),                                 &
        size(potential_x_alpha_ao,1),                                &
        potential_x_alpha_mo(1,1,1),                                 &
        size(potential_x_alpha_mo,1)                                 &
        )

    call ao_to_mo(                                                   &
        potential_x_beta_ao(1,1,1),                                  &
        size(potential_x_beta_ao,1),                                 &
        potential_x_beta_mo(1,1,1),                                  &
        size(potential_x_beta_mo,1)                                  &
        )


    call ao_to_mo(                                                   &
        potential_c_alpha_ao(1,1,1),                                 &
        size(potential_c_alpha_ao,1),                                &
        potential_c_alpha_mo(1,1,1),                                 &
        size(potential_c_alpha_mo,1)                                 &
        )

    call ao_to_mo(                                                   &
        potential_c_beta_ao(1,1,1),                                  &
        size(potential_c_beta_ao,1),                                 &
        potential_c_beta_mo(1,1,1),                                  &
        size(potential_c_beta_mo,1)                                  &
        )



END_PROVIDER 
