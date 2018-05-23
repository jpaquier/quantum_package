program print_H_matrix_restart
 implicit none
 read_wf = .True.
 touch read_wf
 call routine

end

subroutine routine
 use bitmasks
 implicit none
 integer :: i,j
 integer, allocatable :: H_matrix_degree(:,:)
 double precision, allocatable :: H_matrix_phase(:,:)
 integer :: degree
 integer(bit_kind), allocatable :: keys_tmp(:,:,:)
 allocate(keys_tmp(N_int,2,N_det))
 do i = 1, N_det
  print*,''
  call debug_det(psi_det(1,1,i),N_int)
  do j = 1, N_int
   keys_tmp(j,1,i) = psi_det(j,1,i)
   keys_tmp(j,2,i) = psi_det(j,2,i)
  enddo
 enddo
 if(N_det.ge.10000)then
  print*,'Warning !!!'
  print*,'Number of determinants is ',N_det
  print*,'It means that the H matrix will be enormous !'
  print*,'stoppping ..'
  stop
 endif
 print*,''
 print*,'Determinants '
 do i = 1, N_det
 enddo
 allocate(H_matrix_degree(N_det,N_det),H_matrix_phase(N_det,N_det))
 integer         :: exc(0:2,2,2)
 double precision  :: phase
 do i = 1, N_det
  do j = i, N_det 
   call get_excitation_degree(psi_det(1,1,i),psi_det(1,1,j),degree,N_int)
   H_matrix_degree(i,j) = degree
   H_matrix_degree(j,i) = degree
   phase = 0.d0
   if(degree==1.or.degree==2)then
    call get_excitation(psi_det(1,1,i),psi_det(1,1,j),exc,degree,phase,N_int)
   endif
   H_matrix_phase(i,j) = phase
   H_matrix_phase(j,i) = phase
  enddo
 enddo
 print*,'H matrix '
 double precision :: ref_h_matrix,s2
 ref_h_matrix = H_matrix_all_dets(1,1)
 print*,'HF like determinant energy = ',ref_bitmask_energy+nuclear_repulsion
 print*,'Ref element of H_matrix    = ',ref_h_matrix+nuclear_repulsion
 print*,'Printing the H matrix ...'
 print*,''
 print*,''
!do i = 1, N_det
! H_matrix_all_dets(i,i) -= ref_h_matrix
!enddo

 do i = 1, N_det
  H_matrix_all_dets(i,i) += nuclear_repulsion
 enddo

!do i = 5,N_det
! H_matrix_all_dets(i,3) = 0.d0
! H_matrix_all_dets(3,i) = 0.d0
! H_matrix_all_dets(i,4) = 0.d0
! H_matrix_all_dets(4,i) = 0.d0
!enddo




 
 do i = 1, N_det
  write(*,'(I3,X,A3,1000(F16.7))')i,' | ',H_matrix_all_dets(i,:)
 enddo

 print*,''
 print*,''
 print*,''
 print*,'Printing the degree of excitations within the H matrix'
 print*,''
 print*,''
 do i = 1, N_det
  write(*,'(I3,X,A3,X,1000(I1,X))')i,' | ',H_matrix_degree(i,:)
 enddo


 print*,''
 print*,''
 print*,'Printing the phase of the Hamiltonian matrix elements '
 print*,''
 print*,''
 do i = 1, N_det
  write(*,'(I3,X,A3,X,1000(F3.0,X))')i,' | ',H_matrix_phase(i,:)
 enddo
 print*,''


 double precision, allocatable  :: eigenvectors(:,:), eigenvalues(:)
 double precision, allocatable  :: s2_eigvalues(:)
 allocate (eigenvectors(size(H_matrix_all_dets,1),N_det))
 allocate (eigenvalues(N_det),s2_eigvalues(N_det))
 call lapack_diag(eigenvalues,eigenvectors,                       &
     H_matrix_all_dets,size(H_matrix_all_dets,1),N_det)
 print*,'Two first eigenvectors '
 call u_0_S2_u_0(s2_eigvalues,eigenvectors,n_det,keys_tmp,N_int,N_det,size(eigenvectors,1))
 do j =1, N_states
   print*,'s2 = ',s2_eigvalues(j)
   print*,'e  = ',eigenvalues(j)
   print*,'coefs : '
   do i = 1, N_det
    print*,'i = ',i,eigenvectors(i,j)
   enddo
   if(j>1)then
    print*,'Delta E(H)  = ',eigenvalues(1) - eigenvalues(j)
    print*,'Delta E(eV) = ',(eigenvalues(1) - eigenvalues(j))*27.2114d0
   endif
 enddo
end
