 BEGIN_PROVIDER [double precision, energy_x, (N_states)]
&BEGIN_PROVIDER [double precision, energy_c, (N_states)]
&BEGIN_PROVIDER [double precision, potential_x_alpha_ao,(ao_num,ao_num,N_states)]
&BEGIN_PROVIDER [double precision, potential_x_beta_ao,(ao_num,ao_num,N_states)]
&BEGIN_PROVIDER [double precision, potential_c_alpha_ao,(ao_num,ao_num,N_states)]
&BEGIN_PROVIDER [double precision, potential_c_beta_ao,(ao_num,ao_num,N_states)]

 implicit none
 integer :: j,k,l,istate
 integer :: m,n,p,q
 double precision, allocatable :: aos_array(:)
 double precision, allocatable :: r(:)
 double precision :: rho_a(N_states),rho_b(N_states),ex(N_states),ec(N_states)
 double precision :: vx_rho_a(N_states),vx_rho_b(N_states),vc_rho_a(N_states),vc_rho_b(N_states)
 potential_c_alpha_ao = 0.d0
 potential_c_beta_ao = 0.d0
 potential_x_alpha_ao = 0.d0
 potential_x_beta_ao = 0.d0
 double precision, allocatable :: tmp_c_a(:,:,:),tmp_c_b(:,:,:),tmp_x_a(:,:,:),tmp_x_b(:,:,:)
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
     if(DFT_TYPE=="LDA")then
      call dm_dft_alpha_beta_and_all_aos_at_r(r,rho_a,rho_b,aos_array)
      call LDA_type_functional(rho_a,rho_b,vx_rho_a,vx_rho_b,vc_rho_a,vc_rho_b,ex,ec)
      do istate = 1, N_states
       energy_x(istate) += weight *  ex(istate) 
       energy_c(istate) += weight *  ec(istate)
       vx_rho_a(istate) *= weight
       vc_rho_a(istate) *= weight
       vx_rho_b(istate) *= weight
       vc_rho_b(istate) *= weight
      enddo
      call update_potentials_scalar_dger(tmp_c_a,tmp_c_b,tmp_x_a,tmp_x_b,vc_rho_a,vc_rho_b,vx_rho_a,vx_rho_b,aos_array)

     else if (DFT_TYPE=="GGA")then
      double precision :: dvx_rho_a(3,N_states), dvx_rho_b(3,N_states)
      double precision :: dvc_rho_a(3,N_states), dvc_rho_b(3,N_states)
      double precision :: vx_grad_rho_a_2(N_states),vx_grad_rho_b_2(N_states),vx_grad_rho_a_b(N_states)
      double precision :: vc_grad_rho_a_2(N_states),vc_grad_rho_b_2(N_states),vc_grad_rho_a_b(N_states)
      double precision :: grad_rho_a_2(N_states),grad_rho_b_2(N_states),grad_rho_a_b(N_states)
      double precision :: grad_rho_a(3,N_states),grad_rho_b(3,N_states)
      double precision :: grad_aos_array(3,ao_num)
      call density_and_grad_alpha_beta_and_all_aos_and_grad_aos_at_r(r,rho_a,rho_b, grad_rho_a, grad_rho_b, aos_array, grad_aos_array)
      do istate = 1, N_states
       grad_rho_a_2 = 0.d0
       grad_rho_b_2 = 0.d0
       grad_rho_a_b = 0.d0
       do m = 1, 3
        grad_rho_a_2 += grad_rho_a(m,istate)*grad_rho_a(m,istate)
        grad_rho_b_2 += grad_rho_b(m,istate)*grad_rho_b(m,istate)
        grad_rho_a_b += grad_rho_a(m,istate)*grad_rho_b(m,istate)
       enddo
       if(exchange_functional.EQ."short_range_PBE")then
        call ex_pbe_sr(rho_a(istate),rho_b(istate),grad_rho_a_2(istate),grad_rho_b_2(istate),grad_rho_a_b(istate),ex(istate),vx_rho_a(istate),vx_rho_b(istate),vx_grad_rho_a_2(istate),vx_grad_rho_b_2(istate),vx_grad_rho_a_b(istate))
       else if(exchange_functional.EQ."None")then
        ex = 0.d0
        vx_rho_a = 0.d0
        vx_rho_b = 0.d0
        vx_grad_rho_a_2 = 0.d0
        vx_grad_rho_a_b = 0.d0
        vx_grad_rho_b_2 = 0.d0
       else 
        print*, 'Exchange functional required does not exist ...'
        print*,'exchange_functional',exchange_functional
        stop
       endif
        
       if(correlation_functional.EQ."None")then
        ec = 0.d0
        vc_rho_a = 0.d0
        vc_rho_b = 0.d0
        vc_grad_rho_a_2 = 0.d0
        vc_grad_rho_a_b = 0.d0
        vc_grad_rho_b_2 = 0.d0
       else 
        print*, 'Correlation functional required does not exist ...'
        print*, 'correlation_functional',correlation_functional
        stop
       endif
       
       double precision :: grad_aos_array_transpose(ao_num,3)
       double precision :: dtmp_x_a(3,N_states),dtmp_x_b(3,N_states)
       double precision :: dtmp_c_a(3,N_states),dtmp_c_b(3,N_states)
       call dger(ao_num,ao_num,weight*vx_rho_a(istate),aos_array,1,aos_array,1,tmp_x_a(1,1,istate),size(tmp_x_a,1))
       call dger(ao_num,ao_num,weight*vx_rho_b(istate),aos_array,1,aos_array,1,tmp_x_b(1,1,istate),size(tmp_x_b,1))
       call dtranspose(grad_aos_array,3,grad_aos_array_transpose,ao_num,3,ao_num)
       do m = 1, 3
        dtmp_x_a(m,istate) = 2.d0 * vx_grad_rho_a_2(istate) *  grad_rho_a(m,istate) + vx_grad_rho_a_b(istate)  * grad_rho_b(m,istate)
        dtmp_x_b(m,istate) = 2.d0 * vx_grad_rho_b_2(istate) *  grad_rho_b(m,istate) + vx_grad_rho_a_b(istate)  * grad_rho_a(m,istate)
       enddo
      !do m = 1, 3
      ! do n = 1, ao_num
      !  do p =1, ao_num
      !   tmp_x_a(p,n,istate) += dtmp_x_a(m,istate) * (aos_array(p) * grad_aos_array_transpose(n,m) + aos_array(n) * grad_aos_array_transpose(p,m))
      !   tmp_x_b(p,n,istate) += dtmp_x_b(m,istate) * (aos_array(p) * grad_aos_array_transpose(n,m) + aos_array(n) * grad_aos_array_transpose(p,m))
      !  enddo
      ! enddo
      !enddo
       do m= 1,3
        call dger(ao_num,ao_num,weight*dtmp_x_a(m,istate),aos_array,1,grad_aos_array_transpose(1,m),1,tmp_x_a(1,1,istate),size(tmp_x_a,1))
        call dger(ao_num,ao_num,weight*dtmp_x_a(m,istate),grad_aos_array_transpose(1,m),1,aos_array,1,tmp_x_a(1,1,istate),size(tmp_x_a,1))
        call dger(ao_num,ao_num,weight*dtmp_x_b(m,istate),aos_array,1,grad_aos_array_transpose(1,m),1,tmp_x_b(1,1,istate),size(tmp_x_b,1))
        call dger(ao_num,ao_num,weight*dtmp_x_b(m,istate),grad_aos_array_transpose(1,m),1,aos_array,1,tmp_x_b(1,1,istate),size(tmp_x_b,1))
       enddo
      enddo
     endif
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

END_PROVIDER 

subroutine LDA_type_functional(rho_a,rho_b,vx_a,vx_b,vc_a,vc_b,ex,ec)
 implicit none
 double precision, intent(in)  :: rho_a(N_states),rho_b(N_states)
 double precision, intent(out) :: ex(N_states),ec(N_states)
 double precision, intent(out) :: vx_a(N_states),vx_b(N_states),vc_a(N_states),vc_b(N_states)
 integer          :: istate
!!!!!!! EXCHANGE PART
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
   print*,'Exchange functional required does not exist ...'
   print*,'exchange_functional ',exchange_functional
   stop
  endif
!!!!!!! CORRELATION PART
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

end

subroutine update_potentials_scalar_dger(vc_a_array,vc_b_array,vx_a_array,vx_b_array,vc_a,vc_b,vx_a,vx_b,aos_array)
 implicit none
 double precision, intent(inout) :: vc_a_array(ao_num,ao_num,N_states), vc_b_array(ao_num,ao_num,N_states), vx_a_array(ao_num,ao_num,N_states), vx_b_array(ao_num,ao_num,N_states)
 double precision, intent(in)    :: vc_a(N_states), vc_b(N_states), vx_a(N_states), vx_b(N_states), aos_array(ao_num) 
 integer :: istate
 do istate = 1, N_states
  call dger(ao_num,ao_num,vx_a(istate),aos_array,1,aos_array,1,vx_a_array(1,1,istate),size(vx_a_array,1))
  call dger(ao_num,ao_num,vx_b(istate),aos_array,1,aos_array,1,vx_b_array(1,1,istate),size(vx_b_array,1))
  call dger(ao_num,ao_num,vc_a(istate),aos_array,1,aos_array,1,vc_a_array(1,1,istate),size(vc_a_array,1))
  call dger(ao_num,ao_num,vc_b(istate),aos_array,1,aos_array,1,vc_b_array(1,1,istate),size(vc_b_array,1))
 enddo
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
