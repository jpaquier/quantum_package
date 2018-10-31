!BEGIN_PROVIDER [double precision, dirac_energy_x, (N_states)]
!&BEGIN_PROVIDER [double precision, dirac_energy_c, (N_states)]
!&BEGIN_PROVIDER [double precision, dirac_potential_x_ao,(2*dirac_ao_num,2*dirac_ao_num,N_states)]
!&BEGIN_PROVIDER [double precision, dirac_potential_c_ao,(2*dirac_ao_num,2*dirac_ao_num,N_states)]
!implicit none
!integer :: j,k,l,istate
!integer :: m,n,p,q
!double precision, allocatable :: aos_array(:)
!double precision, allocatable :: r(:)
!double precision :: rho(N_states),ex(N_states),ec(N_states)
!double precision :: vx_rho(N_states),vc_rho(N_states)
!double precision, allocatable :: tmp_c(:,:,:),tmp_x(:,:,:)
!double precision :: weight
!double precision :: dvx_rho(3,N_states)
!double precision :: dvc_rho(3,N_states)
!double precision :: vx_grad_rho_2(N_states)
!double precision :: vc_grad_rho_2(N_states)
!double precision :: grad_rho_2(N_states)
!double precision :: grad_rho(3,N_states)
!double precision :: grad_aos_array(3,dirac_ao_num)
!double precision :: grad_aos_array_transpose(dirac_ao_num,3)
!double precision :: dtmp_x(3,N_states)
!double precision :: dtmp_c(3,N_states)
! dirac_energy_x = 0.d0
! dirac_energy_c = 0.d0
! dirac_potential_c_ao = 0.d0
! dirac_potential_x_ao = 0.d0
! do j = 1, nucl_num
!  do k = 1, n_points_radial_grid  -1
!   allocate(tmp_c(2*dirac_ao_num,2*dirac_ao_num,N_states),tmp_x(2*dirac_ao_num,2*dirac_ao_num,N_states),aos_array(dirac_ao_num),r(3))
!   tmp_c = 0.d0
!   tmp_x = 0.d0
!   do l = 1, n_points_integration_angular 
!    r(1) = grid_points_per_atom(1,l,k,j)
!    r(2) = grid_points_per_atom(2,l,k,j)
!    r(3) = grid_points_per_atom(3,l,k,j)
!    weight = final_weight_functions_at_grid_points(l,k,j)
!    if(DIRAC_DFT_TYPE=="LDA")then
!     call dm_dft_and_all_dirac_aos_at_r(r,rho,aos_array)
!     call LDA_type_dirac_functional(rho_a,rho_b,vx_rho_a,vx_rho_b,vc_rho_a,vc_rho_b,ex,ec)
!     do istate = 1, N_states
!      energy_x(istate) += weight *  ex(istate) 
!      energy_c(istate) += weight *  ec(istate)
!      vx_rho_a(istate) *= weight
!      vc_rho_a(istate) *= weight
!      vx_rho_b(istate) *= weight
!      vc_rho_b(istate) *= weight
!     enddo
!  !  call update_potentials_scalar_dger(tmp_c_a,tmp_c_b,tmp_x_a,tmp_x_b,vc_rho_a,vc_rho_b,vx_rho_a,vx_rho_b,aos_array)

!    else if (DFT_TYPE=="GGA")then
!  !  call density_and_grad_alpha_beta_and_all_aos_and_grad_aos_at_r(r,rho_a,rho_b, grad_rho_a, grad_rho_b, aos_array, grad_aos_array)
!     grad_rho_a_2 = 0.d0
!     grad_rho_b_2 = 0.d0
!     grad_rho_a_b = 0.d0
!     do istate = 1, N_states
!      do m = 1, 3
!       grad_rho_a_2(istate) += grad_rho_a(m,istate)*grad_rho_a(m,istate)
!       grad_rho_b_2(istate) += grad_rho_b(m,istate)*grad_rho_b(m,istate)
!       grad_rho_a_b(istate) += grad_rho_a(m,istate)*grad_rho_b(m,istate)
!      enddo
!     enddo
!     call GGA_type_dirac_functionals(rho_a,rho_b,grad_rho_a_2,grad_rho_b_2,grad_rho_a_b,ex,vx_rho_a,vx_rho_b,vx_grad_rho_a_2,vx_grad_rho_b_2,vx_grad_rho_a_b, &  
!                                                                                  ec,vc_rho_a,vc_rho_b,vc_grad_rho_a_2,vc_grad_rho_b_2,vc_grad_rho_a_b )
!     do istate = 1, N_states
!      energy_x(istate) += weight *  ex(istate) 
!      energy_c(istate) += weight *  ec(istate)
!      vx_rho_a(istate) *= weight
!      vc_rho_a(istate) *= weight
!      vx_rho_b(istate) *= weight
!      vc_rho_b(istate) *= weight
!      do m = 1, 3
!       dtmp_x_a(m,istate) = (2.d0 * vx_grad_rho_a_2(istate) *  grad_rho_a(m,istate) + vx_grad_rho_a_b(istate)  * grad_rho_b(m,istate)) * weight
!       dtmp_x_b(m,istate) = (2.d0 * vx_grad_rho_b_2(istate) *  grad_rho_b(m,istate) + vx_grad_rho_a_b(istate)  * grad_rho_a(m,istate)) * weight
!       dtmp_c_a(m,istate) = (2.d0 * vc_grad_rho_a_2(istate) *  grad_rho_a(m,istate) + vc_grad_rho_a_b(istate)  * grad_rho_b(m,istate)) * weight
!       dtmp_c_b(m,istate) = (2.d0 * vc_grad_rho_b_2(istate) *  grad_rho_b(m,istate) + vc_grad_rho_a_b(istate)  * grad_rho_a(m,istate)) * weight
!      enddo
!     enddo
!     call update_potentials_scalar_dger(tmp_c_a,tmp_c_b,tmp_x_a,tmp_x_b,vc_rho_a,vc_rho_b,vx_rho_a,vx_rho_b,aos_array)
!     call dtranspose(grad_aos_array,3,grad_aos_array_transpose,ao_num,3,ao_num)
!     call update_potentials_gradient_dger(dtmp_x_a,dtmp_x_b,dtmp_c_a,dtmp_c_b,aos_array,grad_aos_array_transpose,tmp_x_a,tmp_x_b,tmp_c_a,tmp_c_b)
!    endif
!   enddo ! angular points 

!   do istate = 1,N_states
!    potential_c_ao(:,:,istate) += tmp_c_a(:,:,istate)
!    potential_x_ao(:,:,istate) += tmp_x_a(:,:,istate)
!   enddo
!   deallocate(tmp_x,tmp_c,aos_array,r)
!  enddo
! enddo
!END_PROVIDER 

!subroutine LDA_type_dirac_functional(rho_a,rho_b,vx_a,vx_b,vc_a,vc_b,ex,ec)
! implicit none
! double precision, intent(in)  :: rho_a(N_states),rho_b(N_states)
! double precision, intent(out) :: ex(N_states),ec(N_states)
! double precision, intent(out) :: vx_a(N_states),vx_b(N_states),vc_a(N_states),vc_b(N_states)
! integer          :: istate
!!!!!!!! EXCHANGE PART
! do istate = 1, N_states
!  if(dirac_exchange_functional.EQ."dirac_short_range_LDA")then
! ! call ex_lda_sr(rho_a(istate),rho_b(istate),ex(istate),vx_a(istate),vx_b(istate))
!  else if(dirac_exchange_functional.EQ."None")then
!   ex = 0.d0
!   vx_a = 0.d0
!   vx_b = 0.d0
!  else
!   print*,'Exchange functional required does not exist ...'
!   print*,'dirac_exchange_functional ', dirac_exchange_functional
!   stop
!  endif
!!!!!!!! CORRELATION PART
!  if(dirac_correlation_functional.EQ."dirac_short_range_LDA")then
! ! call ec_lda_sr(rho_a(istate),rho_b(istate),ec(istate),vc_a(istate),vc_b(istate))
!  else if(dirac_correlation_functional.EQ."None")then
!   ec = 0.d0
!   vc_a = 0.d0
!   vc_b = 0.d0
!  else
!   print*, 'Correlation functional required does not exist ...'
!   print*, 'dirac_correlation_functional', dirac_correlation_functional
!   stop
!  endif
! enddo
!end

!subroutine GGA_type_dirac_functionals(rho_a,rho_b,grad_rho_a_2,grad_rho_b_2,grad_rho_a_b, & 
!                               ex,vx_rho_a,vx_rho_b,vx_grad_rho_a_2,vx_grad_rho_b_2,vx_grad_rho_a_b, &  
!                               ec,vc_rho_a,vc_rho_b,vc_grad_rho_a_2,vc_grad_rho_b_2,vc_grad_rho_a_b )
! implicit none
! double precision, intent(in)  :: rho_a(N_states),rho_b(N_states),grad_rho_a_2(N_states),grad_rho_b_2(N_states),grad_rho_a_b(N_states)
! double precision, intent(out) :: ex(N_states),vx_rho_a(N_states),vx_rho_b(N_states),vx_grad_rho_a_2(N_states),vx_grad_rho_b_2(N_states),vx_grad_rho_a_b(N_states)
! double precision, intent(out) :: ec(N_states),vc_rho_a(N_states),vc_rho_b(N_states),vc_grad_rho_a_2(N_states),vc_grad_rho_b_2(N_states),vc_grad_rho_a_b(N_states)
! integer          :: istate
! do istate = 1, N_states
!  if(dirac_exchange_functional.EQ."dirac_short_range_PBE")then
! ! call ex_pbe_sr(rho_a(istate),rho_b(istate),grad_rho_a_2(istate),grad_rho_b_2(istate),grad_rho_a_b(istate),ex(istate),vx_rho_a(istate),vx_rho_b(istate),vx_grad_rho_a_2(istate),vx_grad_rho_b_2(istate),vx_grad_rho_a_b(istate))
!  else if(dirac_exchange_functional.EQ."None")then
!   ex = 0.d0
!   vx_rho_a = 0.d0
!   vx_rho_b = 0.d0
!   vx_grad_rho_a_2 = 0.d0
!   vx_grad_rho_a_b = 0.d0
!   vx_grad_rho_b_2 = 0.d0
!  else 
!   print*, 'Exchange functional required does not exist ...'
!   print*,'dirac_exchange_functional', dirac_exchange_functional
!   stop
!  endif
!  if(dirac_correlation_functional.EQ."None")then
!   ec = 0.d0
!   vc_rho_a = 0.d0
!   vc_rho_b = 0.d0
!   vc_grad_rho_a_2 = 0.d0
!   vc_grad_rho_a_b = 0.d0
!   vc_grad_rho_b_2 = 0.d0
!  else 
!   print*, 'Correlation functional required does not exist ...'
!   print*, 'dirac_correlation_functional', dirac_correlation_functional
!   stop
!  endif
! enddo
!end


!subroutine update_potentials_scalar_dger(vc_a_array,vc_b_array,vx_a_array,vx_b_array,vc_a,vc_b,vx_a,vx_b,aos_array)
! implicit none
! double precision, intent(in)    :: vc_a(N_states), vc_b(N_states), vx_a(N_states), vx_b(N_states), aos_array(ao_num) 
! double precision, intent(inout) :: vc_a_array(ao_num,ao_num,N_states), vc_b_array(ao_num,ao_num,N_states), vx_a_array(ao_num,ao_num,N_states), vx_b_array(ao_num,ao_num,N_states)
! integer :: istate
! do istate = 1, N_states
!  call dger(ao_num,ao_num,vx_a(istate),aos_array,1,aos_array,1,vx_a_array(1,1,istate),size(vx_a_array,1))
!  call dger(ao_num,ao_num,vx_b(istate),aos_array,1,aos_array,1,vx_b_array(1,1,istate),size(vx_b_array,1))
!  call dger(ao_num,ao_num,vc_a(istate),aos_array,1,aos_array,1,vc_a_array(1,1,istate),size(vc_a_array,1))
!  call dger(ao_num,ao_num,vc_b(istate),aos_array,1,aos_array,1,vc_b_array(1,1,istate),size(vc_b_array,1))
! enddo
!end

!subroutine update_potentials_gradient_dger(dtmp_x_a,dtmp_x_b,dtmp_c_a,dtmp_c_b,aos_array,grad_aos_array_transpose,tmp_x_a,tmp_x_b,tmp_c_a,tmp_c_b)
! implicit none
! double precision, intent(in) :: dtmp_x_a(3,N_states),dtmp_x_b(3,N_states)
! double precision, intent(in) :: dtmp_c_a(3,N_states),dtmp_c_b(3,N_states)
! double precision, intent(in) :: aos_array(ao_num),grad_aos_array_transpose(ao_num,3)
! double precision, intent(inout):: tmp_x_a(ao_num,ao_num,N_states), tmp_x_b(ao_num,ao_num,N_states)
! double precision, intent(inout):: tmp_c_a(ao_num,ao_num,N_states), tmp_c_b(ao_num,ao_num,N_states)
! integer :: istate,m
! do istate = 1, N_states
!  do m= 1,3
!   call dger(ao_num,ao_num,dtmp_x_a(m,istate),aos_array,1,grad_aos_array_transpose(1,m),1,tmp_x_a(1,1,istate),size(tmp_x_a,1))
!   call dger(ao_num,ao_num,dtmp_x_a(m,istate),grad_aos_array_transpose(1,m),1,aos_array,1,tmp_x_a(1,1,istate),size(tmp_x_a,1))
!   call dger(ao_num,ao_num,dtmp_x_b(m,istate),aos_array,1,grad_aos_array_transpose(1,m),1,tmp_x_b(1,1,istate),size(tmp_x_b,1))
!   call dger(ao_num,ao_num,dtmp_x_b(m,istate),grad_aos_array_transpose(1,m),1,aos_array,1,tmp_x_b(1,1,istate),size(tmp_x_b,1))

!   call dger(ao_num,ao_num,dtmp_c_a(m,istate),aos_array,1,grad_aos_array_transpose(1,m),1,tmp_c_a(1,1,istate),size(tmp_c_a,1))
!   call dger(ao_num,ao_num,dtmp_c_a(m,istate),grad_aos_array_transpose(1,m),1,aos_array,1,tmp_c_a(1,1,istate),size(tmp_c_a,1))
!   call dger(ao_num,ao_num,dtmp_c_b(m,istate),aos_array,1,grad_aos_array_transpose(1,m),1,tmp_c_b(1,1,istate),size(tmp_c_b,1))
!   call dger(ao_num,ao_num,dtmp_c_b(m,istate),grad_aos_array_transpose(1,m),1,aos_array,1,tmp_c_b(1,1,istate),size(tmp_c_b,1))
!  enddo
! enddo
!end



!BEGIN_PROVIDER [double precision, dirac_ao_potential_xc, (2*dirac_ao_num, 2*dirac_ao_num)] 
!implicit none
!integer :: i,j,k,l
!dirac_ao_potential_xc = 0.d0
! do i = 1, 2*dirac_ao_num
!  do j = 1, 2*dirac_ao_num
!   dirac_ao_potential_xc(i,j)  =  dirac_potential_c_ao(i,j,1)  + dirac_potential_x_ao(i,j,1)
!  enddo
! enddo
!END_PROVIDER 
! 
!BEGIN_PROVIDER [double precision, dirac_e_exchange_dft]
! implicit none
! dirac_e_exchange_dft = dirac_energy_x(1)
!END_PROVIDER 

!BEGIN_PROVIDER [double precision, dirac_e_correlation_dft]
! implicit none
! dirac_e_correlation_dft = dirac_energy_c(1)
!END_PROVIDER 
