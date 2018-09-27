subroutine spherical_averaged_two_dm_at_second_order(r1,r12,istate,two_dm,two_dm_laplacian,total_dm,total)
 implicit none
 integer, intent(in) :: istate
 double precision, intent(in) :: r1(3),r12
 double precision, intent(out) :: two_dm,two_dm_laplacian,total_dm,total
 double precision, allocatable :: mos_array_r1(:)
 double precision, allocatable :: mos_grad_array_r1(:,:)
 double precision, allocatable :: mos_lapl_array_r1(:,:)
 integer :: i,j,k,l,ix
 double precision :: lapl_k,lapl_l,grad_k,grad_l,scalar_grad
 allocate(mos_array_r1(mo_tot_num),mos_grad_array_r1(mo_tot_num,3),mos_lapl_array_r1(mo_tot_num,3))
 call give_all_mos_and_grad_and_lapl_at_r(r1,mos_array_r1,mos_grad_array_r1,mos_lapl_array_r1)
 two_dm = 0.d0
 two_dm_laplacian = 0.d0
 do i = 1, mo_tot_num
  do j = 1, mo_tot_num
   do k = 1, mo_tot_num
    lapl_k =mos_lapl_array_r1(k,1)+mos_lapl_array_r1(k,2)+mos_lapl_array_r1(k,3)
    do l = 1, mo_tot_num
     lapl_l =mos_lapl_array_r1(l,1)+mos_lapl_array_r1(l,2)+mos_lapl_array_r1(l,3)
     grad_l =mos_grad_array_r1(l,1)+mos_grad_array_r1(l,2)+mos_grad_array_r1(l,3)
     scalar_grad = 0.d0
     do ix = 1, 3
      scalar_grad += mos_grad_array_r1(l,ix) * mos_grad_array_r1(k,ix)
     enddo
     two_dm += two_bod_alpha_beta_mo_transposed(l,k,j,i,istate) * mos_array_r1(i) * mos_array_r1(j) * mos_array_r1(k) * mos_array_r1(l)
     two_dm_laplacian += two_bod_alpha_beta_mo_transposed(l,k,j,i,istate) * mos_array_r1(i) * mos_array_r1(j) * (lapl_k * mos_array_r1(l) + 2.d0 * scalar_grad + lapl_l * mos_array_r1(k))
    enddo
   enddo
  enddo
 enddo
 two_dm_laplacian = 1.d0/6.d0 * two_dm_laplacian
 
 total = two_dm + two_dm_laplacian * r12*r12 
 
end



