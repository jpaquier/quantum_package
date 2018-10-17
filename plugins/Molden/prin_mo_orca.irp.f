program pouet
 implicit none
 call routine
end

subroutine test
 implicit none
 integer :: i,j,k,l
 double precision, allocatable :: mo_matrix(:,:)
 allocate(mo_matrix(ao_num,mo_tot_num))
 double precision, allocatable :: mo_matrix_carth(:,:)
 double precision, allocatable :: overlap_cart_to_sphe_2(:,:)
 allocate(mo_matrix_carth(ao_num,mo_tot_num), overlap_cart_to_sphe_2(6,5))
 do i = 1, ao_num
  do j = 1, mo_tot_num 
   mo_matrix(i,j)  = mo_coef(i,j)/ao_coef_normalization_factor(i)
  enddo
 enddo
 do i = 1, 6
  do k = 1, 5
   do l = 1, 6
    overlap_cart_to_sphe_2(i,k) += cart_to_sphe_2(l,k) * ao_overlap(i,l)
   enddo
  enddo
 enddo
 print*,''
 print*,'MO matrix inital '
 do i = 1, ao_num 
  write(*,'(100(F16.10,X))')mo_matrix(i,:)
 enddo
!print*,''
!print*,'overlap_cart_to_sphe_2'
!do i = 1, 6
! write(*,'(100(F16.10,X))')overlap_cart_to_sphe_2(i,:)
!enddo
 
 
 ! copy at hand
 mo_matrix_carth = 0.d0
 ! the d 
 do i = 1, 6
  do j = 1, mo_tot_num
   do k = 1, 5
    mo_matrix_carth(k,j) += mo_matrix(i,j) * overlap_cart_to_sphe_2(i,k) 
   enddo 
  enddo
 enddo
 do j = 1, mo_tot_num
  do k = 6, 9
   mo_matrix_carth(k,j) = mo_matrix(k+1,j)
  enddo
 enddo
 print*,''
 print*,'MO coef total '
 do i = 1, ao_num 
  write(*,'(100(F16.10,X))')mo_matrix_carth(i,:)
 enddo
 
 print*,''
 print*,'cartesian_to_spherical_matrix'
 do i = 1, ao_num
  write(*,'(100(F16.10,X))')cartesian_to_spherical_matrix(i,:)  
 enddo
 mo_matrix = 0.d0
 print*,'n_spherical_AOs_in_basis = ',n_spherical_AOs_in_basis
 print*,'mo_tot_num               = ',mo_tot_num
 do i = 1, mo_tot_num
  do k =1, n_spherical_AOs_in_basis
   do j = 1, ao_num 
   mo_matrix(k,i) += cartesian_to_spherical_matrix(j,k) * mo_coef(j,i)
   enddo
  enddo
 enddo
 print*,''
 print*,'MO coef beforelqskqdvjlkvjlkh'
 do i = 1, n_spherical_AOs_in_basis
  write(*,'(100(F16.10,X))')mo_matrix(i,:)
 enddo

 do i = 1, ao_num
  
 enddo


 double precision, allocatable :: SC(:,:)
 allocate(SC(ao_num,ao_num))
!print*,'overlap'
!do i = 1, ao_num
! write(*,'(100(F16.10,X))')ao_overlap(i,:)
!enddo
!print*,''
!print*,'cart to spher'
!do i = 1, 6
! write(*,'(100(F16.10,X))')cart_to_sphe_2(i,:)
!enddo
!print*,''

 
 print*,'norm factor'
 do i = 1, ao_num
  print*,i,ao_coef_normalization_libint_factor(i),ao_coef_normalization_factor(i) 
 enddo
 print*,''
 print*,'MO coef'
 do i = 1, ao_num 
  write(*,'(I3,X,100(F16.10,X))')i,mo_matrix_carth(i,:)
 enddo
!print*,'ao to carth'
!do i = 1, 6
! write(*,'(100(F16.10,X))')cart_to_sphe_2(i,:) 
!enddo
!print*,'ao to carth'
!do i = 1, 6
! write(*,'(100(F16.10,X))')cart_to_sphe_2(i,:) 
!enddo

end

subroutine routine
 implicit none
 integer :: i
 print*,'norm factor'
 do i = 1, ao_num
  print*,i,ao_coef_normalization_libint_factor(i),ao_coef_normalization_factor(i) 
 enddo
 print*,''
 print*,''
 print*,''
 do i = 1, n_spherical_AOs_in_basis
  write(*,'(I3,X,1000(F16.10,X))')i,mo_coef_spherical(i,46:50)  
 enddo

end
