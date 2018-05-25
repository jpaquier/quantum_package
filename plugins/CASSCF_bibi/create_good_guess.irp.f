subroutine create_good_guess(u,e,i_st)
 implicit none
 integer, intent(in)  :: i_st
 double precision, intent(out) :: u(size_super_ci),e
 integer :: i,j
 double precision :: coef_pert,i_H_SCI_j_state_specific
 integer :: list_good(size_super_ci)
 integer :: n_good
 integer :: ihole, jpart,iorb,jorb
 double precision :: dsqrt_2

 u = 0.d0

 integer, allocatable :: iorder(:)
 double precision, allocatable :: vec_tmp(:)
 allocate(vec_tmp(size_super_ci),iorder(size_super_ci))
 vec_tmp(1) = -1.d10
 iorder(1) = 1
 dsqrt_2 = dsqrt(2.d0)
 do i = 2, size_super_ci
  ihole = index_rotation_CI_reverse(i,1)
  iorb = list_core_inact(ihole)
  jpart = index_rotation_CI_reverse(i,2)
  jorb = list_virt(jpart)
  vec_tmp(i) = -dabs(dsqrt_2 * MR_Fock_matrix_alpha_beta_spin_average_mo(iorb,jorb,i_st) / diagonal_superci_matrix(i,i_st))
  iorder(i) = i
 enddo
 call dsort(vec_tmp,iorder,size_super_ci)

 n_good = 1
 list_good(n_good) = 1
 do i = 2, min(size_super_ci,n_det_max_jacobi)
  n_good += 1 
  list_good(n_good) = iorder(i)
 enddo
 print*, 'n_good = ',n_good
 double precision, allocatable :: matrix(:,:),eigenvectors(:,:),eigenvalues(:)
 allocate(matrix(n_good,n_good),eigenvectors(n_good,n_good),eigenvalues(n_good))
 do i = 1, n_good
  do j = 1, n_good
   matrix(i,j) = i_H_SCI_j_state_specific(list_good(i),list_good(j),i_st)
  enddo
 enddo

 
 call lapack_diag(eigenvalues,eigenvectors,matrix,n_good,n_good)
 e = eigenvalues(1)
 print*, 'e = ',e
 do i = 1, n_good
  u(list_good(i)) = eigenvectors(i,1)
  print*, list_good(i),eigenvectors(i,1)
 enddo
 
 deallocate(matrix,eigenvectors,eigenvalues,vec_tmp, iorder)


end

subroutine create_good_guess_state_average(u,e)
 implicit none
 double precision, intent(out) :: u(size_super_ci),e
 integer :: i,j
 double precision :: coef_pert,i_H_SCI_j_state_average
 integer :: list_good(size_super_ci)
 integer :: n_good
 integer :: ihole, jpart,iorb,jorb
 double precision :: dsqrt_2

 u = 0.d0

 integer, allocatable :: iorder(:)
 double precision, allocatable :: vec_tmp(:)
 allocate(vec_tmp(size_super_ci),iorder(size_super_ci))
 vec_tmp(1) = -1.d10
 iorder(1) = 1
 dsqrt_2 = dsqrt(2.d0)
 do i = 2, size_super_ci
  ihole = index_rotation_CI_reverse(i,1)
  iorb = list_core_inact(ihole)
  jpart = index_rotation_CI_reverse(i,2)
  jorb = list_virt(jpart)
  vec_tmp(i) = -dabs(dsqrt_2 * MR_Fock_matrix_spin_and_state_average_mo(iorb,jorb) / diagonal_superci_matrix_state_average(i))
  iorder(i) = i
 enddo
 call dsort(vec_tmp,iorder,size_super_ci)

 n_good = 1
 list_good(n_good) = 1
 do i = 2, min(size_super_ci,n_det_max_jacobi)
  n_good += 1 
  list_good(n_good) = iorder(i)
 enddo
 print*, 'n_good = ',n_good
 double precision, allocatable :: matrix(:,:),eigenvectors(:,:),eigenvalues(:)
 allocate(matrix(n_good,n_good),eigenvectors(n_good,n_good),eigenvalues(n_good))
 do i = 1, n_good
  do j = 1, n_good
   matrix(i,j) = i_H_SCI_j_state_average(list_good(i),list_good(j))
  enddo
 enddo

 
 call lapack_diag(eigenvalues,eigenvectors,matrix,n_good,n_good)
 e = eigenvalues(1)
 print*, 'e = ',e
 do i = 1, n_good
  u(list_good(i)) = eigenvectors(i,1)
 enddo
 
 deallocate(matrix,eigenvectors,eigenvalues,vec_tmp, iorder)
 

end

