 BEGIN_PROVIDER [double precision, natorb_coef_on_mo_basis, (mo_tot_num,mo_tot_num)]
&BEGIN_PROVIDER [double precision, natural_occ_numbers, (mo_tot_num) ]
 implicit none
 BEGIN_DOC 
 ! natorb_coef_on_mo_basis(i,j) = <i_MO|j_natorb>
 ! |j_natorb> = \sum_{i_MO} R(i,j) |i_MO>
 END_DOC
 integer :: i,j
 double precision, allocatable  ::  A(:,:)
 allocate(A(mo_tot_num,mo_tot_num))
 do i =1, mo_tot_num
  do j = 1, mo_tot_num 
   A(j,i) = - one_body_dm_mo_alpha_average(j,i) - one_body_dm_mo_beta_average(j,i)
  enddo
 enddo
 call lapack_diag(natural_occ_numbers,natorb_coef_on_mo_basis,A,mo_tot_num,mo_tot_num)
 do i = 1, mo_tot_num
  natural_occ_numbers(i) = - natural_occ_numbers(i)
 enddo
END_PROVIDER 

 BEGIN_PROVIDER [double precision, dm_alpha_on_natorb_basis, (mo_tot_num, mo_tot_num)]
&BEGIN_PROVIDER [double precision, dm_beta_on_natorb_basis, (mo_tot_num, mo_tot_num)]
 implicit none
 double precision, allocatable  ::  C(:,:)
 allocate(C(mo_tot_num,mo_tot_num))
 ! C = rho_alpha_mo * natorb_coef_on_mo_basis
 call dgemm('N','N',mo_tot_num,mo_tot_num,mo_tot_num,1.d0,one_body_dm_mo_alpha_average, mo_tot_num,natorb_coef_on_mo_basis, mo_tot_num,0.d0, C,mo_tot_num)
 ! dm_alpha_on_natorb_basis =  natorb_coef_on_mo_basis^T C 
 call dgemm('T','N',mo_tot_num,mo_tot_num,mo_tot_num,1.d0,natorb_coef_on_mo_basis, mo_tot_num,C, mo_tot_num,0.d0, dm_alpha_on_natorb_basis,mo_tot_num)

 ! C = rho_beta_mo * natorb_coef_on_mo_basis
 call dgemm('N','N',mo_tot_num,mo_tot_num,mo_tot_num,1.d0,one_body_dm_mo_beta_average, mo_tot_num,natorb_coef_on_mo_basis, mo_tot_num,0.d0, C,mo_tot_num)
 ! dm_beta_on_natorb_basis =  natorb_coef_on_mo_basis^T C 
 call dgemm('T','N',mo_tot_num,mo_tot_num,mo_tot_num,1.d0,natorb_coef_on_mo_basis, mo_tot_num,C, mo_tot_num,0.d0, dm_beta_on_natorb_basis,mo_tot_num)

END_PROVIDER 
 

