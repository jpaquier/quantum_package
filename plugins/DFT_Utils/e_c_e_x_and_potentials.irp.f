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
    allocate(tmp_c_a(ao_num,ao_num,N_states),tmp_c_b(ao_num,ao_num,N_states),tmp_x_a(ao_num,ao_num,N_states),tmp_x_b(ao_num,ao_num,N_states),aos_array(ao_num),r(3))
    tmp_c_a = 0.d0
    tmp_c_b = 0.d0
    tmp_x_a = 0.d0
    tmp_x_b = 0.d0
    do l = 1, n_points_integration_angular 
     r(1) = grid_points_per_atom(1,l,k,j)
     r(2) = grid_points_per_atom(2,l,k,j)
     r(3) = grid_points_per_atom(3,l,k,j)
     double precision :: weight
     weight = final_weight_functions_at_grid_points(l,k,j)
     call exchange_correlation_functional(r,weight,ec,ex,tmp_c_a,tmp_c_b,tmp_x_a,tmp_x_b)
     do istate = 1, N_states
      energy_x(istate) +=  ex(istate) 
      energy_c(istate) +=  ec(istate) 
     enddo
    enddo
    do istate = 1,N_states
     potential_c_alpha_ao(:,:,istate) = potential_c_alpha_ao(:,:,istate) + tmp_c_a(:,:,istate)
     potential_x_alpha_ao(:,:,istate) = potential_x_alpha_ao(:,:,istate) + tmp_x_a(:,:,istate)
     potential_c_beta_ao(:,:,istate)  =  potential_c_beta_ao(:,:,istate) + tmp_c_b(:,:,istate)
     potential_x_beta_ao(:,:,istate)  =  potential_x_beta_ao(:,:,istate) + tmp_x_b(:,:,istate)
    enddo
    deallocate(tmp_x_a,tmp_x_b,tmp_c_a,tmp_c_b,aos_array,r)
   enddo
  enddo
 print*,'potentials provided !' 

END_PROVIDER 


 subroutine exchange_correlation_functional(r,weight,ec,ex,vc_a_array,vc_b_array,vx_a_array,vx_b_array)
  implicit none
  double precision, intent(in)  :: r(3), weight
  double precision, intent(out) :: ec(N_states),ex(N_states)
  double precision, intent(inout) :: vx_a_array(ao_num,ao_num,N_states),vx_b_array(ao_num,ao_num,N_states)
  double precision, intent(inout) :: vc_a_array(ao_num,ao_num,N_states),vc_b_array(ao_num,ao_num,N_states)

  double precision :: rho_a(N_states),rho_b(N_states)
  double precision :: aos_array(ao_num)
  double precision :: grad_rho_a(3,N_states),grad_rho_b(3,N_states)
  double precision :: grad_aos_array(3,ao_num),grad_aos_array_transpose(ao_num,3)
  double precision :: vx_a(N_states), vx_b(N_states), vc_a(N_states), vc_b(N_states)
  double precision :: contrib_grad_ca, contrib_grad_cb
  double precision :: ao_ao_grad_matrix(ao_num,ao_num)
  double precision :: ao_ao_grad_matrix_transpose(ao_num,ao_num)
  integer :: istate,i,j,k
  if(DFT_TYPE.EQ."LDA")then
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

!!!!!!!!!! CORRELATION PART
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

    double precision :: contrib_xa,contrib_xb,contrib_ca, contrib_cb
    contrib_xa = weight * vx_a(istate) 
    contrib_ca = weight * vc_a(istate) 
    contrib_xb = weight * vx_b(istate) 
    contrib_cb = weight * vc_b(istate) 
    call dger(ao_num,ao_num,contrib_xa,aos_array,1,aos_array,1,vx_a_array(1,1,istate),size(vx_a_array,1))
    call dger(ao_num,ao_num,contrib_ca,aos_array,1,aos_array,1,vc_a_array(1,1,istate),size(vc_a_array,1))
    call dger(ao_num,ao_num,contrib_xb,aos_array,1,aos_array,1,vx_b_array(1,1,istate),size(vx_b_array,1))
    call dger(ao_num,ao_num,contrib_cb,aos_array,1,aos_array,1,vc_b_array(1,1,istate),size(vc_b_array,1))
   enddo

  else if(DFT_TYPE.EQ."GGA")then
   call density_and_grad_alpha_beta_and_all_aos_and_grad_aos_at_r(r,rho_a,rho_b, grad_rho_a, grad_rho_b, aos_array, grad_aos_array)
   do istate = 1, N_states
    double precision :: dvc_a(3,N_states), dvc_b(3,N_states)
    call routine_gga_correlation(r,rho_a(istate),rho_b(istate),grad_rho_a(1,istate),grad_rho_b(1,istate),ec(istate),vc_a(istate),dvc_a(1,istate),vc_b(istate),dvc_b(1,istate))

    contrib_ca=weight*vc_a(istate)
    contrib_cb=weight*vc_b(istate)

!   call dger(ao_num,ao_num,contrib_ca,aos_array,1,aos_array,1,vc_a_array(1,1,istate),size(vc_a_array,1))
!   call dger(ao_num,ao_num,contrib_cb,aos_array,1,aos_array,1,vc_b_array(1,1,istate),size(vc_b_array,1))
    ao_ao_grad_matrix= 0.d0
    call dger(ao_num,ao_num,1.d0,aos_array,1,aos_array,1,ao_ao_grad_matrix(1,1),size(ao_ao_grad_matrix,1))
    vc_a_array(:,:,istate) += contrib_ca * ao_ao_grad_matrix(:,:)
    vc_b_array(:,:,istate) += contrib_cb * ao_ao_grad_matrix(:,:)
    call dtranspose(grad_aos_array,3,grad_aos_array_transpose,ao_num,3,ao_num)
    do k= 1,3
      contrib_grad_ca=weight*dvc_a(k,istate)
      contrib_grad_cb=weight*dvc_b(k,istate)
      ao_ao_grad_matrix= 0.d0
      call dger(ao_num,ao_num,1.D0,aos_array,1,grad_aos_array_transpose(1,k),1,ao_ao_grad_matrix(1,1),size(ao_ao_grad_matrix,1))
      call dtranspose(ao_ao_grad_matrix,ao_num,ao_ao_grad_matrix_transpose,ao_num,ao_num,ao_num)
      ao_ao_grad_matrix += ao_ao_grad_matrix_transpose
      vc_a_array(:,:,istate) += ao_ao_grad_matrix(:,:) * contrib_grad_ca
      vc_b_array(:,:,istate) += ao_ao_grad_matrix(:,:) * contrib_grad_cb
!     call dger(ao_num,ao_num,contrib_grad_ca,grad_aos_array_transpose(1,k),1,aos_array,1,vc_a_array(1,1,istate),size(vc_a_array,1))
!     call dger(ao_num,ao_num,contrib_grad_cb,aos_array,1,grad_aos_array_transpose(1,k),1,vc_b_array(1,1,istate),size(vc_b_array,1))
!     call dger(ao_num,ao_num,contrib_grad_cb,grad_aos_array_transpose(1,k),1,aos_array,1,vc_b_array(1,1,istate),size(vc_b_array,1))
    enddo
   enddo
  endif
  ex = ex * weight
  ec = ec * weight

 end

 subroutine routine_gga_correlation(r,rho_a,rho_b,grad_rho_a,grad_rho_b,ec,vc_a,dvc_a,vc_b,dvc_b)
  implicit none
  double precision, intent(in) :: r(3)
  double precision, intent(in) :: rho_a,rho_b,grad_rho_a(3),grad_rho_b(3)
  double precision, intent(out) :: dvc_b(3),dvc_a(3),vc_a,vc_b,ec 
  ec = rho_a + rho_b
  vc_a = rho_a
  vc_b = rho_b
  dvc_a = grad_rho_a(:)
  dvc_b = grad_rho_b(:)
 end

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


 BEGIN_PROVIDER [double precision, Energy_c_md, (N_states)]
 implicit none
 BEGIN_DOC
 ! Corelation energy for the multi determinent short range LDA. PRB 73 155111 2006
 END_DOC
 integer :: j,k,l,istate 
 double precision, allocatable :: aos_array(:), r(:), rho_a(:), rho_b(:), ec(:)
 logical :: dospin
 dospin = .false.
 allocate(aos_array(ao_num),r(3), rho_a(N_states), rho_b(N_states), ec(N_states))
  do j = 1, nucl_num
   do k = 1, n_points_radial_grid  -1
    do l = 1, n_points_integration_angular 
     
     r(1) = grid_points_per_atom(1,l,k,j)
     r(2) = grid_points_per_atom(2,l,k,j)
     r(3) = grid_points_per_atom(3,l,k,j)
     call dm_dft_alpha_beta_and_all_aos_at_r(r,rho_a,rho_b,aos_array)

     do istate = 1, N_states
!!!!!!!!!!!! CORRELATION PART
      if(md_correlation_functional.EQ."short_range_LDA")then
        call ESRC_MD_LDAERF (rho_a(istate),rho_b(istate),dospin,ec(istate))
      else if(md_correlation_functional.EQ."None")then
       ec = 0.d0
      else
       print*, 'Multi determinant correlation functional required does not exist ...'
       print*, 'md_correlation_functional',md_correlation_functional
       stop
      endif
     enddo
     do istate = 1, N_states
      energy_c_md(istate) += final_weight_functions_at_grid_points(l,k,j) * ec(istate)
     enddo

    enddo
   enddo
  enddo
 deallocate(aos_array,r,rho_a,rho_b, ec)
END_PROVIDER
