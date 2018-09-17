 BEGIN_PROVIDER[double precision, aos_in_r_array, (ao_num,n_points_final_grid)]
&BEGIN_PROVIDER[double precision, aos_in_r_array_transp, (n_points_final_grid,ao_num)]
 implicit none
 integer :: i,j
 double precision :: aos_array(ao_num), r(3)
 do i = 1, n_points_final_grid
  r(1) = final_grid_points(1,i) 
  r(2) = final_grid_points(2,i) 
  r(3) = final_grid_points(3,i) 
  call give_all_aos_at_r(r,aos_array)
  do j = 1, ao_num
   aos_in_r_array(j,i) = aos_array(j)
   aos_in_r_array_transp(i,j) = aos_array(j)
  enddo
 enddo

 END_PROVIDER 
