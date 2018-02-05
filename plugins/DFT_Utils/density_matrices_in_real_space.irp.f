
subroutine dm_dft_alpha_beta_at_r(r,dm_a,dm_b)
 implicit none
 double precision, intent(in) :: r(3)
 double precision, intent(out) :: dm_a(N_states),dm_b(N_states)

 integer :: istate
 double precision  :: aos_array(ao_num),aos_array_bis(ao_num),u_dot_v

 call give_all_aos_at_r(r,aos_array)


 do istate = 1, N_states
  aos_array_bis = aos_array
  ! alpha density
  call dgemv('N',ao_num,ao_num,1.d0,one_body_dm_alpha_ao_for_dft(1,1,istate),ao_num,aos_array,1,0.d0,aos_array_bis,1)
  dm_a(istate) = u_dot_v(aos_array,aos_array_bis,ao_num)
 
  ! beta density
  aos_array_bis = aos_array
  call dgemv('N',ao_num,ao_num,1.d0,one_body_dm_beta_ao_for_dft(1,1,istate),ao_num,aos_array,1,0.d0,aos_array_bis,1)
  dm_b(istate) = u_dot_v(aos_array,aos_array_bis,ao_num)
 enddo

end


subroutine dm_dft_alpha_beta_and_all_aos_at_r(r,dm_a,dm_b,aos_array)
 implicit none
 double precision, intent(in) :: r(3)
 double precision, intent(out) :: dm_a(N_states),dm_b(N_states)
 double precision, intent(out) :: aos_array(ao_num)

 integer :: istate
 double precision  :: aos_array_bis(ao_num),u_dot_v

 call give_all_aos_at_r(r,aos_array)

 do istate = 1, N_states
  aos_array_bis = aos_array
  ! alpha density
  call dsymv('U',ao_num,1.d0,one_body_dm_alpha_ao_for_dft(1,1,istate),size(one_body_dm_alpha_ao_for_dft,1),aos_array,1,0.d0,aos_array_bis,1)
  dm_a(istate) = u_dot_v(aos_array,aos_array_bis,ao_num)
 
  ! beta density
  aos_array_bis = aos_array
  call dsymv('U',ao_num,1.d0,one_body_dm_beta_ao_for_dft(1,1,istate),size(one_body_dm_beta_ao_for_dft,1),aos_array,1,0.d0,aos_array_bis,1)
  dm_b(istate) = u_dot_v(aos_array,aos_array_bis,ao_num)
 enddo

end



subroutine dm_and_gradients_dft_alpha_beta_at_r(r,dm_a,dm_b, dm_a_grad, dm_b_grad)
 implicit none
 double precision, intent(in) :: r(3)
 double precision, intent(out) :: dm_a(N_states),dm_b(N_states)
 double precision, intent(out) :: dm_a_grad(3,N_states),dm_b_grad(3,N_states)

 integer :: i,j,istate
 double precision  :: aos_array(ao_num),aos_array_bis(ao_num),u_dot_v
 double precision  :: aos_grad_array_tmp(3,ao_num),aos_grad_array(ao_num,3)

 call give_all_aos_and_grad_at_r_new(r,aos_array,aos_grad_array_tmp)
 do i = 1, ao_num
  do j = 1, 3
   aos_grad_array(i,j) = aos_grad_array_tmp(j,i)
  enddo
 enddo
 
 do istate = 1, N_states

   aos_array_bis = aos_array
   ! alpha density
   call dgemv('N',ao_num,ao_num,1.d0,one_body_dm_alpha_ao_for_dft(1,1,istate),ao_num,aos_array,1,0.d0,aos_array_bis,1)
   dm_a(istate) = u_dot_v(aos_array,aos_array_bis,ao_num)
!!!double precision  :: accu
!!!accu = 0.d0
!!!do i = 1, ao_num
!!! do j = 1, ao_num 
!!!  accu += aos_array(i) * aos_array(j) * (one_body_dm_alpha_ao_for_dft(i,j)) 
!!! enddo
!!!enddo
!!!if(dabs(accu - dm_a).gt.1.d-10)then
!!! print*,r 
!!! print*,accu,dm_a,dabs(accu-dm_a)
!!! pause
!!!endif
   dm_a_grad(1,istate) = u_dot_v(aos_grad_array(1,1),aos_array_bis,ao_num)
   dm_a_grad(2,istate) = u_dot_v(aos_grad_array(1,2),aos_array_bis,ao_num)
   dm_a_grad(3,istate) = u_dot_v(aos_grad_array(1,3),aos_array_bis,ao_num)
  
   aos_array_bis = aos_array
   ! beta density
   call dgemv('N',ao_num,ao_num,1.d0,one_body_dm_beta_ao_for_dft(1,1,istate),ao_num,aos_array,1,0.d0,aos_array_bis,1)
   dm_b(istate) = u_dot_v(aos_array,aos_array_bis,ao_num)
!!!if(dm_b .ne.dm_a)then
!!! print*,'AHAH' 
!!! print*,dm_a,dm_b,dm_a-dm_b
!!! pause
!!!endif
   dm_b_grad(1,istate) = u_dot_v(aos_grad_array(1,1),aos_array_bis,ao_num)
   dm_b_grad(2,istate) = u_dot_v(aos_grad_array(1,2),aos_array_bis,ao_num)
   dm_b_grad(3,istate) = u_dot_v(aos_grad_array(1,3),aos_array_bis,ao_num)
 enddo
  
end

