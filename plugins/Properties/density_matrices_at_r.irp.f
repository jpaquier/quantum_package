
subroutine density_matrices_alpha_beta_at_r(r,dm_a,dm_b)
 implicit none
 double precision, intent(in) :: r(3)
 double precision, intent(out) :: dm_a,dm_b

 integer :: i,j
 double precision  :: aos_array(ao_num),aos_array_bis(ao_num),u_dot_v

 call give_all_aos_at_r(r,aos_array)


 aos_array_bis = aos_array
 ! alpha density
 call dgemv('N',ao_num,ao_num,1.d0,one_body_dm_ao_alpha,ao_num,aos_array,1,0.d0,aos_array_bis,1)
 dm_a = u_dot_v(aos_array,aos_array_bis,ao_num)

 ! beta density
 aos_array_bis = aos_array
 call dgemv('N',ao_num,ao_num,1.d0,one_body_dm_ao_beta,ao_num,aos_array,1,0.d0,aos_array_bis,1)
 dm_b = u_dot_v(aos_array,aos_array_bis,ao_num)

end


subroutine density_matrices_alpha_beta_and_all_aos_at_r(r,dm_a,dm_b,aos_array)
 implicit none
 double precision, intent(in) :: r(3)
 double precision, intent(out) :: dm_a,dm_b
 double precision, intent(out) :: aos_array(ao_num)

 integer :: i,j
 double precision  :: aos_array_bis(ao_num),u_dot_v

 call give_all_aos_at_r(r,aos_array)

 dm_a = 1.d0
 dm_b = 1.d0

 aos_array_bis = aos_array
 ! alpha density
!call dgemv('N',ao_num,ao_num,1.d0,one_body_dm_ao_alpha,ao_num,aos_array,1,0.d0,aos_array_bis,1)
 call dsymv('U',ao_num,1.d0,one_body_dm_ao_alpha,size(one_body_dm_ao_alpha,1),aos_array,1,0.d0,aos_array_bis,1)
 dm_a = u_dot_v(aos_array,aos_array_bis,ao_num)

 ! beta density
 aos_array_bis = aos_array
!call dgemv('N',ao_num,ao_num,1.d0,one_body_dm_ao_beta,ao_num,aos_array,1,0.d0,aos_array_bis,1)
 call dsymv('U',ao_num,1.d0,one_body_dm_ao_beta,size(one_body_dm_ao_beta,1),aos_array,1,0.d0,aos_array_bis,1)
 dm_b = u_dot_v(aos_array,aos_array_bis,ao_num)

end



subroutine density_matrices_and_gradients_alpha_beta_at_r(r,dm_a,dm_b, dm_a_grad, dm_b_grad)
 implicit none
 double precision, intent(in) :: r(3)
 double precision, intent(out) :: dm_a,dm_b
 double precision, intent(out) :: dm_a_grad(3),dm_b_grad(3)

 integer :: i,j
 double precision  :: aos_array(ao_num),aos_array_bis(ao_num),u_dot_v
 double precision  :: aos_grad_array_tmp(3,ao_num),aos_grad_array(ao_num,3)

 call give_all_aos_and_grad_at_r(r,aos_array,aos_grad_array_tmp)
 do i = 1, ao_num
  do j = 1, 3
   aos_grad_array(i,j) = aos_grad_array_tmp(j,i)
  enddo
 enddo

 aos_array_bis = aos_array


 ! alpha density
 call dgemv('N',ao_num,ao_num,1.d0,one_body_dm_ao_alpha,ao_num,aos_array,1,0.d0,aos_array_bis,1)
 dm_a = u_dot_v(aos_array,aos_array_bis,ao_num)
!double precision  :: accu
!accu = 0.d0
!do i = 1, ao_num
! do j = 1, ao_num 
!  accu += aos_array(i) * aos_array(j) * (one_body_dm_ao_alpha(i,j)) 
! enddo
!enddo
!if(dabs(accu - dm_a).gt.1.d-10)then
! print*,r 
! print*,accu,dm_a,dabs(accu-dm_a)
! pause
!endif
 dm_a_grad(1) = u_dot_v(aos_grad_array(1,1),aos_array_bis,ao_num)
 dm_a_grad(2) = u_dot_v(aos_grad_array(1,2),aos_array_bis,ao_num)
 dm_a_grad(3) = u_dot_v(aos_grad_array(1,3),aos_array_bis,ao_num)

 aos_array_bis = aos_array
 ! beta density
 call dgemv('N',ao_num,ao_num,1.d0,one_body_dm_ao_beta,ao_num,aos_array,1,0.d0,aos_array_bis,1)
 dm_b = u_dot_v(aos_array,aos_array_bis,ao_num)
!if(dm_b .ne.dm_a)then
! print*,'AHAH' 
! print*,dm_a,dm_b,dm_a-dm_b
! pause
!endif
 dm_b_grad(1) = u_dot_v(aos_grad_array(1,1),aos_array_bis,ao_num)
 dm_b_grad(2) = u_dot_v(aos_grad_array(1,2),aos_array_bis,ao_num)
 dm_b_grad(3) = u_dot_v(aos_grad_array(1,3),aos_array_bis,ao_num)

end

