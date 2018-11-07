BEGIN_PROVIDER [double precision, mo_coef_spherical, (n_spherical_AOs_in_basis,mo_tot_num)]
 implicit none
 BEGIN_DOC
! MO coefficients on normalized spherical harmonics  
 END_DOC
 integer :: i,j,k
 mo_coef_spherical = 0.d0
 do i = 1, mo_tot_num
  do k =1, n_spherical_AOs_in_basis
   do j = 1, ao_num 
   mo_coef_spherical(k,i) += cartesian_to_spherical_matrix(j,k) * mo_coef(j,i)/ao_coef_normalization_factor(j)
   enddo
  enddo
 enddo

END_PROVIDER 
