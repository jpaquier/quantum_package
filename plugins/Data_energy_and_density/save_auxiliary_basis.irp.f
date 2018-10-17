program save_auxiliary_mos
  implicit none
  BEGIN_DOC
! TODO
  END_DOC
 integer :: i,j,k
 double precision, allocatable  :: mo_coef_new(:,:), R(:,:),eigvalues(:), A(:,:), C(:,:)
 allocate(mo_coef_new(ao_num,mo_tot_num),R(mo_tot_num,mo_tot_num),eigvalues(mo_tot_num,A(mo_tot_num,mo_tot_num),C(mo_tot_num,mo_tot_num))
 do i =1, mo_tot_num
  do j = 1, mo_tot_num 
   A(j,i) = - one_body_dm_mo_alpha_average(j,i) - one_body_dm_mo_beta_average(j,i)
  enddo
 enddo
 call lapack_diag(eigvalues,R,A,mo_tot_num,mo_tot_num)
 ! R(i,j) = <i_MO|j_natorb>
 ! |j_natorb> = \sum_{i_MO} R(i,j) |i_MO>
 call get_pseudo_inverse(R,mo_tot_num,mo_tot_num,mo_tot_num,C,LDC)
 ! C(i,j) = <i_natorb|j_MO> 
 ! |j_MO> = \sum_{i_natorb} C(i,j) |i_natorb>
 
 
end
