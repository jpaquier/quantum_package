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

 BEGIN_PROVIDER[double precision, aos_grad_in_r_array, (ao_num,n_points_final_grid,3)]
&BEGIN_PROVIDER[double precision, aos_grad_in_r_array_transp, (n_points_final_grid,ao_num,3)]
 implicit none
 integer :: i,j,m
 double precision :: aos_array(ao_num), r(3)
 double precision :: aos_grad_array(3,ao_num)
 do i = 1, n_points_final_grid
  r(1) = final_grid_points(1,i)
  r(2) = final_grid_points(2,i)
  r(3) = final_grid_points(3,i)
  call give_all_aos_and_grad_at_r(r,aos_array,aos_grad_array)
  do m = 1, 3
   do j = 1, ao_num
    aos_grad_in_r_array(j,i,m) = aos_grad_array(m,j)
    aos_grad_in_r_array_transp(i,j,m) = aos_grad_array(m,j)
   enddo
  enddo
 enddo
 END_PROVIDER

 BEGIN_PROVIDER[double precision, aos_lapl_in_r_array, (ao_num,n_points_final_grid,3)]
&BEGIN_PROVIDER[double precision, aos_lapl_in_r_array_transp, (n_points_final_grid,ao_num,3)]
 implicit none
 integer :: i,j,m
 double precision :: aos_array(ao_num), r(3)
 double precision :: aos_grad_array(ao_num,3)
 double precision :: aos_lapl_array(ao_num,3)
 do m = 1, 3
  do i = 1, n_points_final_grid
   r(1) = final_grid_points(1,i)
   r(2) = final_grid_points(2,i)
   r(3) = final_grid_points(3,i)
   call give_all_aos_and_grad_and_lapl_at_r(r,aos_array,aos_grad_array,aos_lapl_array)
   do j = 1, ao_num
    aos_lapl_in_r_array(j,i,m) = aos_lapl_array(j,m)
    aos_lapl_in_r_array_transp(i,j,m) = aos_lapl_array(j,m)
   enddo
  enddo
 enddo
 END_PROVIDER

