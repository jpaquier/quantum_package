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
  energy_x = 0.d0
  energy_c = 0.d0

!print*,'providing the potentials ...'
!double precision :: wall0,wall1
!call wall_time(wall0)
!double precision :: cpu0,cpu1
!call cpu_time(cpu0)
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
     call exchange_correlation_functional_at_r(r,weight,ec,ex,tmp_c_a,tmp_c_b,tmp_x_a,tmp_x_b)
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
!print*,'potentials provided !' 
!call wall_time(wall1)
!print*,'Wall time to provide old potentials   = ',dabs(wall1-wall0)
!call cpu_time(cpu1)
!print*,'CPU  time to provide old potentials   = ',dabs(cpu1-cpu0)

END_PROVIDER 


 subroutine exchange_correlation_functional_at_r(r,weight,ec,ex,vc_a_array,vc_b_array,vx_a_array,vx_b_array)
  implicit none
  double precision, intent(in)  :: r(3), weight
  double precision, intent(out) :: ec(N_states),ex(N_states)
  double precision, intent(inout) :: vx_a_array(ao_num,ao_num,N_states),vx_b_array(ao_num,ao_num,N_states)
  double precision, intent(inout) :: vc_a_array(ao_num,ao_num,N_states),vc_b_array(ao_num,ao_num,N_states)

  double precision :: rho_a(N_states),rho_b(N_states)
  double precision :: aos_array(ao_num)
  double precision :: grad_rho_a(3,N_states),grad_rho_b(3,N_states)
  double precision :: grad_aos_array(3,ao_num)
  double precision :: vx_a(N_states), vx_b(N_states), vc_a(N_states), vc_b(N_states)
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

    double precision :: grad_aos_array_transpose(ao_num,3)
    call dger(ao_num,ao_num,weight*vc_a(istate),aos_array,1,aos_array,1,vc_a_array(1,1,istate),size(vc_a_array,1))
    call dger(ao_num,ao_num,weight*vc_b(istate),aos_array,1,aos_array,1,vc_b_array(1,1,istate),size(vc_b_array,1))
    call dtranspose(grad_aos_array,3,grad_aos_array_transpose,ao_num,3,ao_num)
    do k= 1,3
     call dger(ao_num,ao_num,weight*dvc_a(k,istate),aos_array,1,grad_aos_array_transpose(1,k),1,vc_a_array(1,1,istate),size(vc_a_array,1))
     call dger(ao_num,ao_num,weight*dvc_a(k,istate),grad_aos_array_transpose(1,k),1,aos_array,1,vc_a_array(1,1,istate),size(vc_a_array,1))
     call dger(ao_num,ao_num,weight*dvc_b(k,istate),aos_array,1,grad_aos_array_transpose(1,k),1,vc_b_array(1,1,istate),size(vc_b_array,1))
     call dger(ao_num,ao_num,weight*dvc_b(k,istate),grad_aos_array_transpose(1,k),1,aos_array,1,vc_b_array(1,1,istate),size(vc_b_array,1))
    enddo
   !do j = 1, ao_num
   ! do i = 1, ao_num
   !  vc_a_array(i,j,istate) += weight*vc_a(istate)*aos_array(i)*aos_array(j) 
   !  vc_b_array(i,j,istate) += weight*vc_b(istate)*aos_array(i)*aos_array(j) 
   !  do k= 1,3
   !    vc_a_array(i,j,istate) += weight*dvc_a(k,istate)*(aos_array(i)*grad_aos_array(k,j)+grad_aos_array(k,i)*aos_array(j))
   !    vc_b_array(i,j,istate) += weight*dvc_b(k,istate)*(aos_array(i)*grad_aos_array(k,j)+grad_aos_array(k,i)*aos_array(j))
   !  enddo
   ! enddo  
   !enddo

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


 BEGIN_PROVIDER [double precision, energy_x_new, (N_states)]
&BEGIN_PROVIDER [double precision, energy_c_new, (N_states)]
&BEGIN_PROVIDER [double precision, potential_x_alpha_ao_new,(ao_num,ao_num,N_states)]
&BEGIN_PROVIDER [double precision, potential_x_beta_ao_new,(ao_num,ao_num,N_states)]
&BEGIN_PROVIDER [double precision, potential_c_alpha_ao_new,(ao_num,ao_num,N_states)]
&BEGIN_PROVIDER [double precision, potential_c_beta_ao_new,(ao_num,ao_num,N_states)]

 implicit none
 integer :: j,k,l,istate
 integer :: m,n
 double precision, allocatable :: aos_array(:)
 double precision, allocatable :: r(:)
 double precision :: rho_a(N_states),rho_b(N_states),ex(N_states),ec(N_states)
 double precision :: vx_a(N_states),vx_b(N_states),vc_a(N_states),vc_b(N_states)
 potential_c_alpha_ao_new = 0.d0
 potential_c_beta_ao_new = 0.d0
 potential_x_alpha_ao_new = 0.d0
 potential_x_beta_ao_new = 0.d0
 double precision, allocatable :: tmp_c_a(:,:,:),tmp_c_b(:,:,:),tmp_x_a(:,:,:),tmp_x_b(:,:,:)
 double precision, allocatable :: ao_matrix(:,:),diag_matrix(:,:,:)
  energy_x_new = 0.d0
  energy_c_new = 0.d0
!print*,'providing the potentials ...'
!double precision :: wall0,wall1
!call wall_time(wall0)
!double precision :: cpu0,cpu1
!call cpu_time(cpu0)
  allocate(ao_matrix(ao_num,n_points_integration_angular),diag_matrix(n_points_integration_angular,n_points_integration_angular,N_states))
  do j = 1, nucl_num
   do k = 1, n_points_radial_grid  -1
    allocate(tmp_c_a(ao_num,ao_num,N_states),tmp_c_b(ao_num,ao_num,N_states),tmp_x_a(ao_num,ao_num,N_states),tmp_x_b(ao_num,ao_num,N_states),aos_array(ao_num),r(3))
    call exchange_correlation_functional_at_radial_point(j,k,ec,ex,tmp_c_a,tmp_c_b,tmp_x_a,tmp_x_b)
    do istate = 1, N_states
     energy_x_new(istate) +=  ex(istate) 
     energy_c_new(istate) +=  ec(istate) 
    enddo
    do istate = 1,N_states
     potential_c_alpha_ao_new(:,:,istate) = potential_c_alpha_ao_new(:,:,istate) + tmp_c_a(:,:,istate)
     potential_x_alpha_ao_new(:,:,istate) = potential_x_alpha_ao_new(:,:,istate) + tmp_x_a(:,:,istate)
     potential_c_beta_ao_new(:,:,istate)  =  potential_c_beta_ao_new(:,:,istate) + tmp_c_b(:,:,istate)
     potential_x_beta_ao_new(:,:,istate)  =  potential_x_beta_ao_new(:,:,istate) + tmp_x_b(:,:,istate)
    enddo
    deallocate(tmp_x_a,tmp_x_b,tmp_c_a,tmp_c_b,aos_array,r)
   enddo
  enddo
!print*,'potentials provided !' 
!call wall_time(wall1)
!print*,'Wall time to provide new potentials   = ',dabs(wall1-wall0)
!call cpu_time(cpu1)
!print*,'CPU  time to provide new potentials   = ',dabs(cpu1-cpu0)

END_PROVIDER 

subroutine exchange_correlation_functional_at_radial_point(i_nucl,j_rad,ec,ex,v_array_c_a,v_array_c_b,v_array_x_a,v_array_x_b)
 implicit none
 BEGIN_DOC
! computes the exchange/correlation contributions together with all potentials for a a given radial point j_rad on a given nuclei i_nucl
! i_nucl = index of the current nuclei
! j_rad = index of the current radial point
! ec = correlation contribution
! ex = exchange contribution
! v_array_c_a = contribution to the alpha correlation potential 
! v_array_c_b = contribution to the beta correlation potential 
! v_array_x_a = contribution to the alpha exchange potential 
! v_array_x_b = contribution to the beta exchange potential 
 END_DOC
 integer, intent(in) :: i_nucl,j_rad
 double precision, intent(out) :: ec(N_states),ex(N_states)
 double precision, intent(out) :: v_array_c_a(ao_num,ao_num,N_states)
 double precision, intent(out) :: v_array_c_b(ao_num,ao_num,N_states)
 double precision, intent(out) :: v_array_x_a(ao_num,ao_num,N_states)
 double precision, intent(out) :: v_array_x_b(ao_num,ao_num,N_states)

 integer :: l,i,istate
 double precision :: r(3),diag_matrix(n_points_integration_angular,n_points_integration_angular,N_states)
 double precision :: ao_matrix(ao_num,n_points_integration_angular)
 double precision :: aos_array(ao_num)
 double precision :: ao_matrix_vc_a(ao_num,n_points_integration_angular,N_states)
 double precision :: ao_matrix_vc_b(ao_num,n_points_integration_angular,N_states)
 double precision :: ao_matrix_vx_a(ao_num,n_points_integration_angular,N_states)
 double precision :: ao_matrix_vx_b(ao_num,n_points_integration_angular,N_states)
 double precision :: ex_tmp(N_states),ec_tmp(N_states),rho_a(N_states),rho_b(N_states)
 double precision :: vx_a(N_states),vx_b(N_states),vc_a(N_states),vc_b(N_states)
 double precision :: weight
 double precision :: contrib_xa,contrib_xb,contrib_ca, contrib_cb
 ex = 0.d0
 ec = 0.d0

 if(DFT_TYPE.EQ."LDA")then
  do l = 1, n_points_integration_angular 
   r(1) = grid_points_per_atom(1,l,j_rad,i_nucl)
   r(2) = grid_points_per_atom(2,l,j_rad,i_nucl)
   r(3) = grid_points_per_atom(3,l,j_rad,i_nucl)
   call dm_dft_alpha_beta_and_all_aos_at_r(r,rho_a,rho_b,aos_array)
   do i = 1, ao_num
    ao_matrix(i,l) = aos_array(i)
   enddo
!!!!!!!!!!! EXCHANGE PART
   do istate = 1, N_states
    if(exchange_functional.EQ."short_range_LDA")then
     call ex_lda_sr(rho_a(istate),rho_b(istate),ex_tmp(istate),vx_a(istate),vx_b(istate))
    else if(exchange_functional.EQ."LDA")then
     call ex_lda(rho_a(istate),rho_b(istate),ex_tmp(istate),vx_a(istate),vx_b(istate))
    else if(exchange_functional.EQ."None")then
     ex_tmp = 0.d0
     vx_a = 0.d0
     vx_b = 0.d0
    else
     print*, 'Exchange functional required does not ex_tmpist ...'
     print*,'exchange_functional',exchange_functional
     stop
    endif
!!!!!!!!!! CORRELATION PART
    if(correlation_functional.EQ."short_range_LDA")then
     call ec_lda_sr(rho_a(istate),rho_b(istate),ec_tmp(istate),vc_a(istate),vc_b(istate))
    else if(correlation_functional.EQ."LDA")then
     call ec_lda(rho_a(istate),rho_b(istate),ec_tmp(istate),vc_a(istate),vc_b(istate))
    else if(correlation_functional.EQ."None")then
     ec_tmp = 0.d0
     vc_a = 0.d0
     vc_b = 0.d0
    else
     print*, 'Correlation functional required does not ex_tmpist ...'
     print*, 'correlation_functional',correlation_functional
     stop
    endif
    weight = final_weight_functions_at_grid_points(l,j_rad,i_nucl)
    ex(istate) += weight * ex_tmp(istate) 
    ec(istate) += weight * ec_tmp(istate) 

    contrib_xa = weight * vx_a(istate) 
    contrib_ca = weight * vc_a(istate) 
    contrib_xb = weight * vx_b(istate) 
    contrib_cb = weight * vc_b(istate) 
    do i = 1, ao_num
     ao_matrix_vc_a(i,l,istate) = contrib_ca * aos_array(i)
     ao_matrix_vc_b(i,l,istate) = contrib_cb * aos_array(i)
     ao_matrix_vx_a(i,l,istate) = contrib_xa * aos_array(i)
     ao_matrix_vx_b(i,l,istate) = contrib_xb * aos_array(i)
    enddo
   enddo
  enddo
  do istate = 1, N_states
  ! matrix product ! 
   call dgemm('N','T',ao_num,ao_num,n_points_integration_angular,1.d0,ao_matrix_vc_a(1,1,istate),size(ao_matrix_vc_a,1),ao_matrix,size(ao_matrix,1),0.d0,v_array_c_a(1,1,istate),size(v_array_c_a,1))
   call dgemm('N','T',ao_num,ao_num,n_points_integration_angular,1.d0,ao_matrix_vc_b(1,1,istate),size(ao_matrix_vc_a,1),ao_matrix,size(ao_matrix,1),0.d0,v_array_c_b(1,1,istate),size(v_array_c_b,1))
   call dgemm('N','T',ao_num,ao_num,n_points_integration_angular,1.d0,ao_matrix_vx_a(1,1,istate),size(ao_matrix_vx_a,1),ao_matrix,size(ao_matrix,1),0.d0,v_array_x_a(1,1,istate),size(v_array_x_a,1))
   call dgemm('N','T',ao_num,ao_num,n_points_integration_angular,1.d0,ao_matrix_vx_b(1,1,istate),size(ao_matrix_vx_a,1),ao_matrix,size(ao_matrix,1),0.d0,v_array_x_b(1,1,istate),size(v_array_x_b,1))
  enddo
 else if(DFT_TYPE.EQ."GGA")then
   !call density_and_grad_alpha_beta_and_all_aos_and_grad_aos_at_r(r,rho_a,rho_b, grad_rho_a, grad_rho_b, aos_array, grad_aos_array)
   !do istate = 1, N_states
   ! double precision :: dvc_a(3,N_states), dvc_b(3,N_states)
   ! call routine_gga_correlation(r,rho_a(istate),rho_b(istate),grad_rho_a(1,istate),grad_rho_b(1,istate),ec(istate),vc_a(istate),dvc_a(1,istate),vc_b(istate),dvc_b(1,istate))
 
   ! double precision :: grad_aos_array_transpose(ao_num,3)
   ! call dger(ao_num,ao_num,weight*vc_a(istate),aos_array,1,aos_array,1,vc_a_array(1,1,istate),size(vc_a_array,1))
   ! call dger(ao_num,ao_num,weight*vc_b(istate),aos_array,1,aos_array,1,vc_b_array(1,1,istate),size(vc_b_array,1))
   ! call dtranspose(grad_aos_array,3,grad_aos_array_transpose,ao_num,3,ao_num)
   ! do k= 1,3
   !  call dger(ao_num,ao_num,weight*dvc_a(k,istate),aos_array,1,grad_aos_array_transpose(1,k),1,vc_a_array(1,1,istate),size(vc_a_array,1))
   !  call dger(ao_num,ao_num,weight*dvc_a(k,istate),grad_aos_array_transpose(1,k),1,aos_array,1,vc_a_array(1,1,istate),size(vc_a_array,1))
   !  call dger(ao_num,ao_num,weight*dvc_b(k,istate),aos_array,1,grad_aos_array_transpose(1,k),1,vc_b_array(1,1,istate),size(vc_b_array,1))
   !  call dger(ao_num,ao_num,weight*dvc_b(k,istate),grad_aos_array_transpose(1,k),1,aos_array,1,vc_b_array(1,1,istate),size(vc_b_array,1))
   ! enddo
   !!do j = 1, ao_num
   !! do i = 1, ao_num
   !!  vc_a_array(i,j,istate) += weight*vc_a(istate)*aos_array(i)*aos_array(j) 
   !!  vc_b_array(i,j,istate) += weight*vc_b(istate)*aos_array(i)*aos_array(j) 
   !!  do k= 1,3
   !!    vc_a_array(i,j,istate) += weight*dvc_a(k,istate)*(aos_array(i)*grad_aos_array(k,j)+grad_aos_array(k,i)*aos_array(j))
   !!    vc_b_array(i,j,istate) += weight*dvc_b(k,istate)*(aos_array(i)*grad_aos_array(k,j)+grad_aos_array(k,i)*aos_array(j))
   !!  enddo
   !! enddo  
   !!enddo
 
   !enddo
 endif
  


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
