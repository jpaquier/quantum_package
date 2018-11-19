
 BEGIN_PROVIDER [double precision, mu_of_r_hf_coalescence_vector, (n_points_final_grid) ]
&BEGIN_PROVIDER [double precision, mu_of_r_hf_coalescence, (n_points_integration_angular,n_points_radial_grid,nucl_num) ]
 implicit none 
 BEGIN_DOC
 ! mu_of_r and mu_average computation 
 END_DOC
 integer :: i_point,k,i,j
 double precision :: r(3)
 double precision :: cpu0,cpu1,local_potential
 print*,'providing the mu_of_r_hf_coalescence_vector ...'
 call wall_time(cpu0)
 r = 0.d0
 call local_r12_operator_on_hf(r,r,local_potential)
 !$OMP PARALLEL DO &
 !$OMP DEFAULT (NONE)  &
 !$OMP PRIVATE (i_point,r,local_potential) & 
 !$OMP ShARED (n_points_final_grid,final_grid_points,mu_of_r_hf_coalescence_vector) 
 do i_point = 1, n_points_final_grid
  r(1) = final_grid_points(1,i_point)
  r(2) = final_grid_points(2,i_point)
  r(3) = final_grid_points(3,i_point)
  call local_r12_operator_on_hf(r,r,local_potential)
  mu_of_r_hf_coalescence_vector(i_point) =  local_potential * dsqrt(dacos(-1.d0)) * 0.5d0
 enddo
 !$OMP END PARALLEL DO
 do i_point = 1, n_points_final_grid
  k = index_final_points(1,i_point)
  i = index_final_points(2,i_point)
  j = index_final_points(3,i_point)
  mu_of_r_hf_coalescence(k,i,j) = mu_of_r_hf_coalescence_vector(i_point)
 enddo
 call wall_time(cpu1)
 print*,'Time to provide mu_of_r_hf_coalescence_vector = ',cpu1-cpu0
 END_PROVIDER 

 BEGIN_PROVIDER [double precision, mu_of_r_psi_coalescence_vector, (n_points_final_grid) ]
&BEGIN_PROVIDER [double precision, mu_of_r_psi_coalescence, (n_points_integration_angular,n_points_radial_grid,nucl_num) ]
 implicit none 
 BEGIN_DOC
 ! mu_of_r and mu_average computation 
 END_DOC
 integer :: i_point,k,i,j
 double precision :: r(3)
 double precision :: cpu0,cpu1,local_potential,two_body_dm
 print*,'providing the mu_of_r_psi_coalescence_vector ...'
 call wall_time(cpu0)
!r = 0.d0
!call local_r12_operator_on_hf(r,r,local_potential)
!!$OMP PARALLEL DO &
!!$OMP DEFAULT (NONE)  &
!!$OMP PRIVATE (i_point,r,local_potential,two_body_dm) & 
!!$OMP shARED (n_points_final_grid,final_grid_points,mu_of_r_psi_coalescence_vector) 
!do i_point = 1, n_points_final_grid
! r(1) = final_grid_points(1,i_point)
! r(2) = final_grid_points(2,i_point)
! r(3) = final_grid_points(3,i_point)
! call expectation_value_in_real_space(r,r,local_potential,two_body_dm)
! mu_of_r_psi_coalescence_vector(i_point) =  local_potential * dsqrt(dacos(-1.d0)) * 0.5d0
!enddo
!!$OMP END PARALLEL DO

 if(.True.)then
  provide on_top_of_r_vector_simple 
  provide f_psi_B
 endif
 !$OMP PARALLEL DO &
 !$OMP DEFAULT (NONE)  &
 !$OMP PRIVATE (i_point,r,local_potential,two_body_dm) & 
 !$OMP shARED (n_points_final_grid,final_grid_points,mu_of_r_psi_coalescence_vector,f_psi_B,on_top_of_r_vector_simple) 
 do i_point = 1, n_points_final_grid
  local_potential = f_psi_B(i_point) / on_top_of_r_vector_simple(i_point,1)
  if(on_top_of_r_vector_simple(i_point,1).gt.1.d-12.and.f_psi_B(i_point).gt.1.d-12)then
   local_potential = f_psi_B(i_point)/on_top_of_r_vector_simple(i_point,1)
  else 
   local_potential = 1.d-10
  endif
  mu_of_r_psi_coalescence_vector(i_point) =  local_potential * dsqrt(dacos(-1.d0)) * 0.5d0
 enddo
 !$OMP END PARALLEL DO


 do i_point = 1, n_points_final_grid
  k = index_final_points(1,i_point)
  i = index_final_points(2,i_point)
  j = index_final_points(3,i_point)
  mu_of_r_psi_coalescence(k,i,j) = mu_of_r_psi_coalescence_vector(i_point)
 enddo
 call wall_time(cpu1)
 print*,'Time to provide mu_of_r_psi_coalescence_vector = ',cpu1-cpu0
 END_PROVIDER 


