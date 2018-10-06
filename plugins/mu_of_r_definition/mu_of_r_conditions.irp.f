
 BEGIN_PROVIDER [double precision, mu_of_r, (n_points_integration_angular,n_points_radial_grid,nucl_num) ]
&BEGIN_PROVIDER [double precision, mu_of_r_vector, (n_points_final_grid) ]
 implicit none 
 BEGIN_DOC
 ! mu_of_r and mu_average computation 
 END_DOC
 integer :: j,k,l
 double precision, allocatable :: r(:)
 double precision :: local_potential,two_body_dm
 double precision :: cpu0,cpu1,integral_f,mu_integral,spherical_average
 print*,'providing the mu_of_r ...'
 call cpu_time(cpu0)
 allocate(r(3))
 do j = 1, nucl_num
  do k = 1, n_points_radial_grid  -1
   spherical_average = 0.d0
   do l = 1, n_points_integration_angular 
    r(1) = grid_points_per_atom(1,l,k,j)
    r(2) = grid_points_per_atom(2,l,k,j)
    r(3) = grid_points_per_atom(3,l,k,j)
    if(mu_of_r_potential.EQ."cusp_condition")then
     mu_of_r(l,k,j) = mu_of_r_cusp_condition(l,k,j) 
    else if(mu_of_r_potential.EQ."hf_coallescence")then
     call local_r12_operator_on_hf(r,r,local_potential)
     mu_of_r(l,k,j) =  local_potential * dsqrt(dacos(-1.d0)) * 0.5d0
    else if(mu_of_r_potential.EQ."psi_coallescence")then
     call expectation_value_in_real_space(r,r,local_potential,two_body_dm)
     mu_of_r(l,k,j) =  local_potential * dsqrt(dacos(-1.d0)) * 0.5d0
    else if(mu_of_r_potential.EQ."integral_hf")then
      call integral_of_f_12_on_hf(r,integral_f)
      mu_of_r(l,k,j) = mu_integral(integral_f,r)
    else 
      print*,'you requested the following mu_of_r_potential'
      print*,mu_of_r_potential
      print*,'which does not correspond to any of the options for such keyword'
      stop
    endif
    spherical_average += final_weight_functions_at_grid_points(l,k,j)*mu_of_r(l,k,j)
    mu_of_r_vector(index_final_points_reverse(l,k,j)) = mu_of_r(l,k,j)
   enddo
   double precision :: distance
   spherical_average = 1.d0/(4.d0 * dacos(-1.d0)) * spherical_average
   distance = dsqrt( (grid_points_per_atom(1,1,k,j) - nucl_coord(1,j))**2  + (grid_points_per_atom(2,1,k,j) - nucl_coord(2,j))**2  + (grid_points_per_atom(3,1,k,j) - nucl_coord(3,j))**2 )
!  write(33,'(100(F16.10,X))')distance,spherical_average,(one_body_dm_mo_alpha_at_grid_points(l,k,j,1)+one_body_dm_mo_beta_at_grid_points(l,k,j,1))
  enddo
 enddo
 call cpu_time(cpu1)
 print*,'Time to provide mu_of_r = ',cpu1-cpu0
 deallocate(r)
 END_PROVIDER 

 BEGIN_PROVIDER [double precision, mu_of_r_cusp_condition, (n_points_integration_angular,n_points_radial_grid,nucl_num) ]
 implicit none 
 BEGIN_DOC
 ! mu_of_r and mu_average computation 
 END_DOC
 integer :: i,j,k,i_point
 double precision :: r(3),two_dm,two_dm_laplacian,total_dm
 do i_point = 1, n_points_final_grid
  k = index_final_points(1,i_point)
  i = index_final_points(2,i_point)
  j = index_final_points(3,i_point)
  mu_of_r_cusp_condition(k,i,j) = mu_of_r_cusp_condition_vector(i_point,1)
 enddo

 END_PROVIDER 


