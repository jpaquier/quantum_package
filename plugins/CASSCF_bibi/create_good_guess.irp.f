subroutine create_good_guess(u,e,i_st)
 implicit none
 integer, intent(in)  :: i_st
 double precision, intent(out) :: u(size_super_ci),e
 integer :: i,j
 double precision :: coef_pert,i_H_SCI_j
 integer :: list_good(size_super_ci)
 integer :: n_good
 integer :: ihole, jpart,iorb,jorb
 double precision :: dsqrt_2
 dsqrt_2 = dsqrt(2.d0)
 u = 0.d0
 n_good = 1
 list_good(n_good) = 1
 do i = 2, size_super_ci
  ihole = index_rotation_CI_reverse(i,1)
  iorb = list_core_inact(ihole)
  jpart = index_rotation_CI_reverse(i,2)
  jorb = list_virt(jpart)
  coef_pert = dabs(dsqrt_2 * Fock_matrix_alpha_beta_spin_average_mo(iorb,jorb,i_st) / diagonal_superci_matrix(i,i_st))
! if(diagonal_superci_matrix(i).lt.0.d0.or.coef_pert.gt.0.3d0.and.n_good.lt.10000)then
  if(coef_pert.gt.0.003d0.and.n_good.lt.10000)then
   n_good += 1 
   list_good(n_good) = i
  endif
 enddo
 print*, 'n_good = ',n_good
 double precision, allocatable :: matrix(:,:),eigenvectors(:,:),eigenvalues(:)
 allocate(matrix(n_good,n_good),eigenvectors(n_good,n_good),eigenvalues(n_good))
 do i = 1, n_good
  do j = 1, n_good
   matrix(i,j) = i_H_SCI_j(list_good(i),list_good(j),i_st)
  enddo
 enddo

 
 call lapack_diag(eigenvalues,eigenvectors,matrix,n_good,n_good)
 e = eigenvalues(1)
 print*, 'e = ',e
 do i = 1, n_good
  u(list_good(i)) = eigenvectors(i,1)
  print*, list_good(i),eigenvectors(i,1)
 enddo
 
 deallocate(matrix,eigenvectors,eigenvalues)

end
