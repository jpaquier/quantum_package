
BEGIN_PROVIDER [integer, n_points_final_grid]
 BEGIN_DOC
! number of points which are non zero 
 END_DOC
 integer :: i,j,k,l
 double precision :: density, mos_array(mo_tot_num),r(3)
 n_points_final_grid = 0 
  do j = 1, nucl_num  
   do i = 1, n_points_radial_grid -1 
    do k = 1, n_points_integration_angular  
     r(1) = grid_points_per_atom(1,k,i,j)
     r(2) = grid_points_per_atom(2,k,i,j)
     r(3) = grid_points_per_atom(3,k,i,j)
     call give_all_mos_at_r(r,mos_array)
     density = 0.d0
     do l = 1, elec_alpha_num 
      density += 2.d0 * mos_array(l)**2
     enddo
     density = density * 0.5d0 /dble(elec_alpha_num)
!    if(dabs(final_weight_functions_at_grid_points(k,i,j) * density ).gt.threshold_grid_dft)then
      n_points_final_grid += 1
!    endif
    enddo
   enddo
  enddo
 print*,'n_points_final_grid = ',n_points_final_grid
 print*,'n max point         = ',n_points_integration_angular*n_points_radial_grid*nucl_num
END_PROVIDER 

 BEGIN_PROVIDER [double precision, final_grid_points, (3,n_points_final_grid)]
&BEGIN_PROVIDER [double precision, final_weight_functions_at_final_grid_points, (n_points_final_grid) ]
&BEGIN_PROVIDER [integer, index_final_points, (3,n_points_final_grid) ]
 implicit none
 integer :: i,j,k,l,i_count
 double precision :: density, mos_array(mo_tot_num),r(3)
 i_count = 0
  do j = 1, nucl_num  
   do i = 1, n_points_radial_grid -1 
    do k = 1, n_points_integration_angular  
     r(1) = grid_points_per_atom(1,k,i,j)
     r(2) = grid_points_per_atom(2,k,i,j)
     r(3) = grid_points_per_atom(3,k,i,j)
     call give_all_mos_at_r(r,mos_array)
     density = 0.d0
     do l = 1, elec_alpha_num 
      density += 2.d0 * mos_array(l)**2
     enddo
     density = density * 0.5d0 /dble(elec_alpha_num)
    !if(dabs(final_weight_functions_at_grid_points(k,i,j) * density ).gt.threshold_grid_dft)then
      i_count += 1
      final_grid_points(1,i_count) = grid_points_per_atom(1,k,i,j)
      final_grid_points(2,i_count) = grid_points_per_atom(2,k,i,j)
      final_grid_points(3,i_count) = grid_points_per_atom(3,k,i,j)
      final_weight_functions_at_final_grid_points(i_count) = final_weight_functions_at_grid_points(k,i,j)
      index_final_points(1,i_count) = k
      index_final_points(2,i_count) = i
      index_final_points(3,i_count) = j
    !endif
    enddo
   enddo
  enddo

 END_PROVIDER 
