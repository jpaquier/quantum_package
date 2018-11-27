 BEGIN_PROVIDER[double precision, mos_in_r_array, (mo_tot_num,n_points_final_grid)]
&BEGIN_PROVIDER[double precision, mos_in_r_array_transp,(n_points_final_grid,mo_tot_num)]
 implicit none
 integer :: i,j
 double precision :: mos_array(mo_tot_num), r(3)
 do i = 1, n_points_final_grid
  r(1) = final_grid_points(1,i)
  r(2) = final_grid_points(2,i)
  r(3) = final_grid_points(3,i)
  call give_all_mos_at_r(r,mos_array)
  do j = 1, mo_tot_num
   mos_in_r_array(j,i) = mos_array(j)
   mos_in_r_array_transp(i,j) = mos_array(j)
  enddo
 enddo
 END_PROVIDER


 BEGIN_PROVIDER[double precision, mos_grad_in_r_array,(mo_tot_num,n_points_final_grid,3)]
 implicit none
 integer :: m
 mos_grad_in_r_array = 0.d0
 do m=1,3
  call dgemm('N','N',mo_tot_num,n_points_final_grid,ao_num,1.d0,mo_coef_transp,mo_tot_num,aos_grad_in_r_array(1,1,m),ao_num,0.d0,mos_grad_in_r_array(1,1,m),mo_tot_num)
 enddo
 END_PROVIDER

 BEGIN_PROVIDER[double precision, mos_lapl_in_r_array,(mo_tot_num,n_points_final_grid,3)]
 implicit none
 integer :: m
 mos_lapl_in_r_array = 0.d0
 do m=1,3
  call dgemm('N','N',mo_tot_num,n_points_final_grid,ao_num,1.d0,mo_coef_transp,mo_tot_num,aos_lapl_in_r_array(1,1,m),ao_num,0.d0,mos_lapl_in_r_array(1,1,m),mo_tot_num)
 enddo
 END_PROVIDER


