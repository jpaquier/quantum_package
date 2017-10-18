subroutine diag_matrix_mo(matrix, ld_mo, list_orb, n_orb, ld_ao,mo_coef_new)
 implicit none
 integer, intent(in) :: n_orb,ld_mo,ld_ao
 integer, intent(in) :: list_orb(n_orb)
 double precision, intent(in) :: matrix(ld_mo,mo_tot_num)
 double precision, intent(inout) :: mo_coef_new(ld_ao,mo_tot_num)
 
 integer :: i,j,iorb,jorb,k,l
 double precision :: matrix_tmp(n_orb,n_orb)
 double precision :: mo_coef_tmp(ld_ao,n_orb)

 do i = 1, n_orb
  iorb = list_orb(i)
  do j = 1, n_orb
   jorb = list_orb(j)
   accu = 0.d0 
   do k = 1, ao_num
    do l = 1, ao_num
     accu += mo_coef_new(l,iorb) * mo_coef_new(k,jorb) * ao_overlap(l,k)
    enddo
   enddo
!  if(iorb == jorb.and. dabs(accu-1.d0).gt.1.d-12)then
!   print*, 'Warning !! In diag_matrix_mo routine, the old MOs might not be orthogonal ! '
!   print*, iorb,jorb,accu
!  endif
!  if(iorb .ne. jorb.and. dabs(accu).gt.1.d-12)then
!   print*, 'Warning !! In diag_matrix_mo routine, the old MOs might not be orthogonal ! '
!   print*, iorb,jorb,accu
!  endif
  enddo
 
 enddo


 do i = 1, n_orb
  iorb = list_orb(i)
  do j = 1, ao_num
   mo_coef_tmp(j,i) = mo_coef_new(j,iorb)
  enddo
  do j = 1, n_orb
   jorb = list_orb(j)
   matrix_tmp(j,i) = matrix(jorb,iorb)
  enddo
 enddo
 double precision :: eigvalues(n_orb),eigvectors(n_orb,n_orb)
 
!print*, 'matrix to be diagonalized '
!print*, 'n_orb',n_orb
 do i = 1, n_orb
! write(*,'(100(F10.5,X))')matrix_tmp(i,:)
  do j = 1, n_orb
   if(dabs(matrix_tmp(i,j) - matrix_tmp(j,i)).gt.1.d-10)then
    print*, 'Warning !! the matrix you want to diagonalized is not symmetric!'
    print*, i,j
    print*, matrix_tmp(j,i),matrix_tmp(i,j),dabs(matrix_tmp(i,j) - matrix_tmp(j,i))
   endif
  enddo
 enddo
 call lapack_diag(eigvalues,eigvectors,matrix_tmp,n_orb,n_orb)
 
 do i = 1, n_orb
! print*, 'eigenvector ',i
! write(*,'(100(F10.5,X))')eigvectors(:,i)
  do j = 1, n_orb 
   accu = 0.d0
   do  k = 1, n_orb
     accu += eigvectors(k,i) * eigvectors(k,j)
   enddo
 ! if(i==j)then
 !  if(dabs(accu)-1.d0.gt.1.d-12)then
 !   print*, 'Warning !! In diag_matrix_mo routine, the new eigenvectors might not be orthogonal ! '
 !   print*, i,j,accu
 !  endif
 ! endif
 ! if(i.ne.j)then
 !  if(dabs(accu).gt.1.d-12)then
 !   print*, 'Warning !! In diag_matrix_mo routine, the new eigenvectors might not be orthogonal ! '
 !   print*, i,j,accu
 !  endif
 ! endif
  enddo
 enddo

!!!!! copy the new eigenvectors
 do i = 1, n_orb 
  iorb = list_orb(i)
  do k = 1, ao_num
   mo_coef_new(k,iorb) = 0.d0
  enddo
  do j = 1, n_orb
   do k = 1, ao_num
    mo_coef_new(k,iorb) += eigvectors(j,i) * mo_coef_tmp(k,j)
   enddo
  enddo
 enddo

 double precision :: accu
 do i = 1, n_orb
  iorb = list_orb(i)
  do j = 1, n_orb
   jorb = list_orb(j)
   accu = 0.d0 
   do k = 1, ao_num
    do l = 1, ao_num
     accu += mo_coef_new(l,iorb) * mo_coef_new(k,jorb) * ao_overlap(l,k)
    enddo
   enddo
  !if(iorb == jorb.and. dabs(accu-1.d0).gt.1.d-12)then
  ! print*, 'Warning !! In diag_matrix_mo routine, the new MOs might not be orthogonal ! '
  ! print*, iorb,jorb,accu
  !endif
  !if(iorb .ne. jorb.and. dabs(accu).gt.1.d-12)then
  ! print*, 'Warning !! In diag_matrix_mo routine, the new MOs might not be orthogonal ! '
  ! print*, iorb,jorb,accu
  !endif
  enddo
 
 enddo

end
