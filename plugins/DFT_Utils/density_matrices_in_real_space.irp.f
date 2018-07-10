
subroutine dm_dft_alpha_beta_at_r(r,dm_a,dm_b)
 implicit none
 BEGIN_DOC
! input: r(1) ==> r(1) = x, r(2) = y, r(3) = z
! output : dm_a = alpha density evaluated at r(3)
! output : dm_b = beta  density evaluated at r(3)
 END_DOC
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
 BEGIN_DOC
! input: r(1) ==> r(1) = x, r(2) = y, r(3) = z
! output : dm_a = alpha density evaluated at r
! output : dm_b = beta  density evaluated at r
! output : aos_array(i) = ao(i) evaluated at r
 END_DOC
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


subroutine density_and_grad_alpha_beta_and_all_aos_and_grad_aos_at_r(r,dm_a,dm_b, grad_dm_a, grad_dm_b, aos_array, grad_aos_array)
 implicit none
 BEGIN_DOC
! input      : r(1) ==> r(1) = x, r(2) = y, r(3) = z
! output     : dm_a = alpha density evaluated at r
!            : dm_b = beta  density evaluated at r
!            : aos_array(i) = ao(i) evaluated at r
!            : grad_dm_a(1) = X gradient of the alpha density evaluated in r
!            : grad_dm_a(1) = X gradient of the beta  density evaluated in r
!            : grad_aos_array(1) = X gradient of the aos(i) evaluated at r
 END_DOC
 double precision, intent(in)  :: r(3)
 double precision, intent(out) :: dm_a(N_states),dm_b(N_states)
 double precision, intent(out) :: grad_dm_a(3,N_states),grad_dm_b(3,N_states)
 double precision, intent(out) :: grad_aos_array(3,ao_num)

 integer :: i,j,istate
 double precision  :: aos_array(ao_num),aos_array_bis(ao_num),u_dot_v
 double precision  :: aos_grad_array(ao_num,3), aos_grad_array_bis(ao_num,3)

 call give_all_aos_and_grad_at_r(r,aos_array,grad_aos_array)
 do i = 1, ao_num
  do j = 1, 3
   aos_grad_array(i,j) = grad_aos_array(j,i)
  enddo
 enddo
 
 do istate = 1, N_states

 ! aos_array_bis = aos_array
   ! alpha density
   ! aos_array_bis = \rho_ao * aos_array
   call dsymv('U',ao_num,1.d0,one_body_dm_alpha_ao_for_dft(1,1,istate),size(one_body_dm_beta_ao_for_dft,1),aos_array,1,0.d0,aos_array_bis,1)
   dm_a(istate) = u_dot_v(aos_array,aos_array_bis,ao_num)
   
   ! grad_dm(1) = \sum_i aos_grad_array(i,1) * aos_array_bis(i)
   grad_dm_a(1,istate) = u_dot_v(aos_grad_array(1,1),aos_array_bis,ao_num)
   grad_dm_a(2,istate) = u_dot_v(aos_grad_array(1,2),aos_array_bis,ao_num)
   grad_dm_a(3,istate) = u_dot_v(aos_grad_array(1,3),aos_array_bis,ao_num)
   grad_dm_a *= 2.d0
   ! aos_grad_array_bis = \rho_ao * aos_grad_array
 ! call dsymv('U',ao_num,1.d0,one_body_dm_alpha_ao_for_dft(1,1,istate),size(one_body_dm_beta_ao_for_dft,1),aos_grad_array,1,0.d0,aos_grad_array_bis,1)
 ! ! grad_dm(1) += \sum_i aos_grad_array(i,1) * aos_array_bis(i)
 ! grad_dm_a(1,istate)+= u_dot_v(aos_array(1),aos_grad_array_bis,ao_num)
 ! grad_dm_a(2,istate)+= u_dot_v(aos_array(1),aos_grad_array_bis,ao_num)
 ! grad_dm_a(3,istate)+= u_dot_v(aos_array(1),aos_grad_array_bis,ao_num)
  
  !aos_array_bis = aos_array
   ! beta density
   call dsymv('U',ao_num,1.d0,one_body_dm_beta_ao_for_dft(1,1,istate),size(one_body_dm_beta_ao_for_dft,1),aos_array,1,0.d0,aos_array_bis,1)
   dm_b(istate) = u_dot_v(aos_array,aos_array_bis,ao_num)

   ! grad_dm(1) = \sum_i aos_grad_array(i,1) * aos_array_bis(i)
   grad_dm_b(1,istate) = u_dot_v(aos_grad_array(1,1),aos_array_bis,ao_num)
   grad_dm_b(2,istate) = u_dot_v(aos_grad_array(1,2),aos_array_bis,ao_num)
   grad_dm_b(3,istate) = u_dot_v(aos_grad_array(1,3),aos_array_bis,ao_num)
   grad_dm_b *= 2.d0
   ! aos_grad_array_bis = \rho_ao * aos_grad_array
 ! call dsymv('U',ao_num,1.d0,one_body_dm_beta_ao_for_dft(1,1,istate),size(one_body_dm_beta_ao_for_dft,1),aos_grad_array,1,0.d0,aos_grad_array_bis,1)
 ! ! grad_dm(1) += \sum_i aos_grad_array(i,1) * aos_array_bis(i)
 ! grad_dm_b(1,istate)+= u_dot_v(aos_array(1),aos_grad_array_bis,ao_num)
 ! grad_dm_b(2,istate)+= u_dot_v(aos_array(1),aos_grad_array_bis,ao_num)
 ! grad_dm_b(3,istate)+= u_dot_v(aos_array(1),aos_grad_array_bis,ao_num)
 enddo
  
end


subroutine on_top_pair_density_in_real_space(r,rho2)
 implicit none
 double precision, intent(in)  :: r(3)
 double precision, intent(out) :: rho2(N_states)
 rho2 = 0.d0
 double precision :: mos_array(mo_tot_num)
 double precision :: aos_array(mo_tot_num)
 
 call give_all_mos_at_r(r,mos_array) 
!call give_all_aos_at_r(r,mos_array)
 double precision :: threshold,c1,c2,c3
 double precision :: a1,a2,a3
 integer :: i,j,k,l,istate
 threshold = 1.d-11 
 rho2 = 0.d0
 do istate = 1, N_states
  do i = 1, mo_tot_num
   a1 = mos_array(i)
   c1 = dabs(a1)
!  if(c1.le.threshold)cycle
   do j = 1, mo_tot_num 
    a2 = a1 * mos_array(j)
    c2 = dabs(a2)
!   if(dabs(one_body_dm_mo_for_dft(j,i,istate)).le.threshold)cycle
!   if(c2.le.threshold)cycle
    do k = 1, mo_tot_num 
     a3 = a2 * mos_array(k)
     c3 = dabs(a3)
!    if(c3.le.threshold)cycle
     do l = 1, mo_tot_num
      rho2 += a3 * mos_array(l) * two_bod_alpha_beta_mo(l,k,j,i,istate)
!     print*,a3,mos_array(l),two_bod_alpha_beta_mo(l,k,j,i,istate)
     enddo
    enddo
   enddo
  enddo
 enddo
end

subroutine on_top_pair_density_in_real_space_from_ao(r,rho2)
 implicit none
 double precision, intent(in)  :: r(3)
 double precision, intent(out) :: rho2(N_states)
 rho2 = 0.d0
 double precision :: aos_array(ao_num)
 
 call give_all_aos_at_r(r,aos_array)

 double precision :: threshold,c1,c2,c3
 double precision :: a1,a2,a3
 integer :: i,j,k,l,istate
 threshold = 1.d-8 
 rho2 = 0.d0
 do istate = 1, N_states
  do i = 1, ao_num
   a1 = aos_array(i)
   c1 = dabs(a1)
   if(c1.le.threshold)cycle
   do j = 1, ao_num 
    a2 = a1 * aos_array(j)
    c2 = dabs(a2)
    if(c2.le.threshold)cycle
    do k = 1, ao_num 
     a3 = a2 * aos_array(k)
     c3 = dabs(a3)
     if(c3.le.threshold)cycle
     do l = 1, ao_num
      rho2 += a3 * aos_array(l) * two_bod_alpha_beta_ao(l,k,j,i,istate)
     enddo
    enddo
   enddo
  enddo
 enddo
end
