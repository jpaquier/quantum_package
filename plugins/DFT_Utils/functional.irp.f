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
 double precision :: aos_array(ao_num)
 double precision :: r(3)
 double precision :: rho_a,rho_b,ex,ec
 double precision :: vx_a,vx_b,vc_a,vc_b
 potential_c_alpha_ao = 0.d0
 potential_c_beta_ao = 0.d0
 potential_x_alpha_ao = 0.d0
 potential_x_beta_ao = 0.d0
 print*,'providing the potentials ...'
      print*,'exchange_functional',exchange_functional
 do i = 1, N_states
  energy_x(i) = 0.d0
  energy_c(i) = 0.d0
  do j = 1, nucl_num
   do k = 1, n_points_radial_grid  -1
    do l = 1, n_points_integration_angular 
     rho_a = one_body_dm_mo_alpha_at_grid_points(l,k,j,i)
     rho_b =  one_body_dm_mo_beta_at_grid_points(l,k,j,i)

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
     r(1) = grid_points_per_atom(1,l,k,j) 
     r(2) = grid_points_per_atom(2,l,k,j) 
     r(3) = grid_points_per_atom(3,l,k,j) 
     call give_all_aos_at_r(r,aos_array)
     do m = 1, ao_num
      do n = 1, ao_num
       potential_x_alpha_ao(m,n,i) += (vx_a) * aos_array(m)*aos_array(n)  * final_weight_functions_at_grid_points(l,k,j)
       potential_x_beta_ao(m,n,i)  += (vx_b) * aos_array(m)*aos_array(n)  * final_weight_functions_at_grid_points(l,k,j)
       potential_c_alpha_ao(m,n,i) += (vc_a) * aos_array(m)*aos_array(n)  * final_weight_functions_at_grid_points(l,k,j)
       potential_c_beta_ao(m,n,i)  += (vc_b) * aos_array(m)*aos_array(n)  * final_weight_functions_at_grid_points(l,k,j)
      enddo
     enddo
    enddo
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
