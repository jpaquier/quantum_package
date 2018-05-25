
BEGIN_PROVIDER [integer, n_couple_ao_ao ]
 implicit none
 n_couple_ao_ao = ao_num * (ao_num  + 1) /2 
END_PROVIDER 

 BEGIN_PROVIDER [integer, index_ao_ao, (ao_num, ao_num)]
&BEGIN_PROVIDER [integer, index_ao_ao_reverse, (n_couple_ao_ao,2)]
 implicit none
 integer :: i,j,n_couple
 n_couple = 0
 do i = 1, ao_num
  do j = i, ao_num
   n_couple +=1 
   index_ao_ao(i,j) = n_couple
   index_ao_ao(j,i) = n_couple
   index_ao_ao_reverse(n_couple,1) = i
   index_ao_ao_reverse(n_couple,2) = j
  enddo
 enddo
 if(n_couple .ne.n_couple_ao_ao)then
  print*, 'PB !!! n_couple .ne.n_couple_ao_ao' 
 endif
 
END_PROVIDER 

BEGIN_PROVIDER [integer, n_couple_core_inact_virt]
 implicit none
 n_couple_core_inact_virt = n_core_inact_orb * n_virt_orb
END_PROVIDER 

 BEGIN_PROVIDER [integer, index_core_inact_virt, (n_core_orb, n_virt_orb)]
&BEGIN_PROVIDER [integer, index_core_inact_virt_reverse, (n_couple_core_inact_virt,2)]
 implicit none
 integer :: i,j,n_couple
 n_couple = 0
 do j = 1, n_virt_orb 
  do i = 1, n_core_orb
   n_couple +=1 
   index_core_inact_virt(i,j) = n_couple
   index_core_inact_virt(j,i) = n_couple
  enddo
 enddo
 if(n_couple .ne.n_couple_core_inact_virt)then
  print*, 'PB !!! n_couple .ne.n_couple_core_inact_virt' 
 endif
 
END_PROVIDER 

 BEGIN_PROVIDER [integer, n_occ_virt_rotations]
&BEGIN_PROVIDER [integer, size_super_ci]
 implicit none
 n_occ_virt_rotations = n_core_inact_orb*n_virt_orb
 size_super_ci = n_occ_virt_rotations + 1
 print*, 'size_super_ci = ',size_super_ci

END_PROVIDER 


 BEGIN_PROVIDER [integer, index_rotation_CI, (n_core_inact_orb, n_virt_orb)]
&BEGIN_PROVIDER [integer, index_rotation_CI_reverse, (size_super_ci,2)]
 implicit none
 BEGIN_DOC
! index_rotation_CI(i_core,j_virt) = index of the excitation i_core --> j_virt in the SUPER CI expansion
! index_rotation_CI_reverse(i,1) = i_core  
! index_rotation_CI_reverse(i,2) = j_virt 
 END_DOC
 integer :: i,j,k
 k = 1 
 do j = 1, n_virt_orb
  do i = 1, n_core_inact_orb
   k+=1
   index_rotation_CI(i,j) = k
   index_rotation_CI_reverse(k,1) = i  ! core_inact
   index_rotation_CI_reverse(k,2) = j  ! virt 
  enddo
 enddo
END_PROVIDER 


BEGIN_PROVIDER [double precision, diagonal_superci_matrix, (size_super_ci,N_states)]
 implicit none
 integer :: i,iorb,j,jorb,k,korb,l,lorb,m
 diagonal_superci_matrix = 0.d0
 print*, 'type_of_superci ',type_of_superci
 do m = 1, N_States
  do i = 1, n_core_inact_orb
   iorb = list_core_inact(i)
   do j = 1, n_virt_orb
    jorb = list_virt(j)
    if(type_of_superci == 0)then
     diagonal_superci_matrix(index_rotation_CI(i,j),m) = - MR_Fock_matrix_alpha_beta_spin_average_mo(iorb,iorb,m) + MR_Fock_matrix_alpha_beta_spin_average_mo(jorb,jorb,m) 
    else 
     diagonal_superci_matrix(index_rotation_CI(i,j),m) = - MR_Fock_matrix_alpha_beta_spin_average_mo(iorb,iorb,m) + MR_Fock_matrix_alpha_beta_spin_average_mo(jorb,jorb,m) & 
                                                       - transformed_occ1_virt2_virt2(i,j,j)          + 2.d0 * transformed_occ1_virt1_occ2_virt2(i,j,i,j)
    endif
   enddo
  enddo
 enddo


END_PROVIDER 



BEGIN_PROVIDER [double precision, superci_matrix, (size_super_ci,size_super_ci,N_states)]
 implicit none
 integer :: i,iorb,j,jorb,k,korb,l,lorb,m
 double precision :: dsqrt_2
 dsqrt_2 = dsqrt(2.d0)
 superci_matrix = 0.d0
 print*, 'Providing the superci matrix'
 
 if(type_of_superci == 0)then
  do m = 1, N_states
   do i = 1, n_core_inact_orb
    iorb = list_core_inact(i)
    do j = 1, n_virt_orb
     jorb = list_virt(j)
     ! Diagonal and Brillouin matrix elements 
     superci_matrix(index_rotation_CI(i,j),index_rotation_CI(i,j),m) =  - MR_Fock_matrix_alpha_beta_spin_average_mo(iorb,iorb,m) + MR_Fock_matrix_alpha_beta_spin_average_mo(jorb,jorb,m) 
     superci_matrix(1,index_rotation_CI(i,j),m) = dsqrt_2 * MR_Fock_matrix_alpha_beta_spin_average_mo(iorb,jorb,m)
     superci_matrix(index_rotation_CI(i,j),1,m) = dsqrt_2 * MR_Fock_matrix_alpha_beta_spin_average_mo(jorb,iorb,m)
     ! Interaction through the virt-virt MR_Fock operator
     do k = j+1, n_virt_orb
      korb = list_virt(k)
      superci_matrix(index_rotation_CI(i,j),index_rotation_CI(i,k),m) = MR_Fock_matrix_alpha_beta_spin_average_mo(jorb,korb,m)
      superci_matrix(index_rotation_CI(i,k),index_rotation_CI(i,j),m) = MR_Fock_matrix_alpha_beta_spin_average_mo(jorb,korb,m)
     enddo
     ! Interaction through the core-core MR_Fock operator
     do k = i+1, n_core_inact_orb
      korb = list_core_inact(k)
      superci_matrix(index_rotation_CI(i,j),index_rotation_CI(k,j),m) = - MR_Fock_matrix_alpha_beta_spin_average_mo(iorb,korb,m)
      superci_matrix(index_rotation_CI(k,j),index_rotation_CI(i,j),m) = - MR_Fock_matrix_alpha_beta_spin_average_mo(iorb,korb,m)
     enddo
    enddo
   enddo
  enddo

 else if(type_of_superci == 1)then


 do m = 1, N_states
  do i = 1, n_core_inact_orb
   iorb = list_core_inact(i)
   do j = 1, n_virt_orb
    jorb = list_virt(j)
    ! Diagonal and Brillouin matrix elements 
    superci_matrix(index_rotation_CI(i,j),index_rotation_CI(i,j),m) = -MR_Fock_matrix_alpha_beta_spin_average_mo(iorb,iorb,m) + MR_Fock_matrix_alpha_beta_spin_average_mo(jorb,jorb,m) & 
                                                                   - transformed_occ1_virt2_virt2(i,j,j)          + 2.d0 * transformed_occ1_virt1_occ2_virt2(i,j,i,j)
    superci_matrix(1,index_rotation_CI(i,j),m) = dsqrt_2 * MR_Fock_matrix_alpha_beta_spin_average_mo(iorb,jorb,m)
    superci_matrix(index_rotation_CI(i,j),1,m) = dsqrt_2 * MR_Fock_matrix_alpha_beta_spin_average_mo(jorb,iorb,m)
    ! Interaction through the virt-virt MR_Fock operator
    do k = j+1, n_virt_orb
     korb = list_virt(k)
     superci_matrix(index_rotation_CI(i,j),index_rotation_CI(i,k),m) = MR_Fock_matrix_alpha_beta_spin_average_mo(jorb,korb,m)  & 
                                                                   - transformed_occ1_virt2_virt2(i,k,j) + 2.d0 * transformed_occ1_virt1_occ2_virt2(i,k,i,j)
     superci_matrix(index_rotation_CI(i,k),index_rotation_CI(i,j),m) = superci_matrix(index_rotation_CI(i,j),index_rotation_CI(i,k),m)
    enddo
    ! Interaction through the core-core MR_Fock operator
    do k = i+1, n_core_inact_orb
     korb = list_core_inact(k)
     superci_matrix(index_rotation_CI(i,j),index_rotation_CI(k,j),m) = - MR_Fock_matrix_alpha_beta_spin_average_mo(iorb,korb,m) & 
                                                                     - transformed_virt1_occ2_occ2(j,k,i) + 2.d0 * transformed_occ1_virt1_occ2_virt2(i,j,k,j)
     superci_matrix(index_rotation_CI(k,j),index_rotation_CI(i,j),m) = superci_matrix(index_rotation_CI(i,j),index_rotation_CI(k,j),m)
    enddo
    
   enddo
  enddo
 enddo
 else if(type_of_superci == 2)then

 do m = 1, N_states
  do i = 1, n_core_inact_orb
   iorb = list_core_inact(i)
   do j = 1, n_virt_orb
    jorb = list_virt(j)
    ! Diagonal and Brillouin matrix elements 
    superci_matrix(index_rotation_CI(i,j),index_rotation_CI(i,j),m) = -MR_Fock_matrix_alpha_beta_spin_average_mo(iorb,iorb,m) + MR_Fock_matrix_alpha_beta_spin_average_mo(jorb,jorb,m) & 
                                                                   - transformed_occ1_virt2_virt2(i,j,j)          + 2.d0 * transformed_occ1_virt1_occ2_virt2(i,j,i,j)
    superci_matrix(1,index_rotation_CI(i,j),m) = dsqrt_2 * MR_Fock_matrix_alpha_beta_spin_average_mo(iorb,jorb,m)
    superci_matrix(index_rotation_CI(i,j),1,m) = dsqrt_2 * MR_Fock_matrix_alpha_beta_spin_average_mo(jorb,iorb,m)
    ! Interaction through the virt-virt MR_Fock operator
    do k = j+1, n_virt_orb
     korb = list_virt(k)
     superci_matrix(index_rotation_CI(i,j),index_rotation_CI(i,k),m) = MR_Fock_matrix_alpha_beta_spin_average_mo(jorb,korb,m)  & 
                                                                   - transformed_occ1_virt2_virt2(i,k,j) + 2.d0 * transformed_occ1_virt1_occ2_virt2(i,k,i,j)
     superci_matrix(index_rotation_CI(i,k),index_rotation_CI(i,j),m) = superci_matrix(index_rotation_CI(i,j),index_rotation_CI(i,k),m)
    enddo
    ! Interaction through the core-core MR_Fock operator
    do k = i+1, n_core_inact_orb
     korb = list_core_inact(k)
     superci_matrix(index_rotation_CI(i,j),index_rotation_CI(k,j),m) = - MR_Fock_matrix_alpha_beta_spin_average_mo(iorb,korb,m) & 
                                                                     - transformed_virt1_occ2_occ2(j,k,i) + 2.d0 * transformed_occ1_virt1_occ2_virt2(i,j,k,j)
     superci_matrix(index_rotation_CI(k,j),index_rotation_CI(i,j),m) = superci_matrix(index_rotation_CI(i,j),index_rotation_CI(k,j),m)
    enddo
    ! Hole-particle interaction 
    do l = 1, j-1
     do k = 1, i-1
      superci_matrix(index_rotation_CI(i,j),index_rotation_CI(k,l),m) = 2.d0 * transformed_occ1_virt1_occ2_virt2(k,l,i,j) - transformed_virt1_virt1_occ2_occ2(l,j,k,i)
     enddo
    enddo

    do l = 1, j-1
     do k = i+1, n_core_inact_orb
      superci_matrix(index_rotation_CI(i,j),index_rotation_CI(k,l),m) = 2.d0 * transformed_occ1_virt1_occ2_virt2(k,l,i,j) - transformed_virt1_virt1_occ2_occ2(l,j,k,i)
     enddo
    enddo

    do l = j+1, n_virt_orb
     do k = 1, i-1
      superci_matrix(index_rotation_CI(i,j),index_rotation_CI(k,l),m) = 2.d0 * transformed_occ1_virt1_occ2_virt2(k,l,i,j) - transformed_virt1_virt1_occ2_occ2(l,j,k,i)
     enddo
    enddo

    do l = j+1, n_virt_orb
     do k = i+1, n_core_inact_orb
      superci_matrix(index_rotation_CI(i,j),index_rotation_CI(k,l),m) = 2.d0 * transformed_occ1_virt1_occ2_virt2(k,l,i,j) - transformed_virt1_virt1_occ2_occ2(l,j,k,i)
     enddo
    enddo
    
   enddo
  enddo
 enddo

 endif

END_PROVIDER 


 BEGIN_PROVIDER [double precision, eigenvectors_sci, (size_super_ci,N_states)]
&BEGIN_PROVIDER [double precision, eigenvalues_sci, (N_states)]
 implicit none
 integer :: i,j,iter,m
 double precision :: e_guess
 double precision, allocatable :: grd_st_eigenvec(:),eigenvectors(:,:),eigenvalues(:),matrix_tmp(:,:)
 double precision :: u_dot_v
 double precision, allocatable :: delta_H_array(:)
 allocate(delta_H_array(size_super_ci))

 if(size_super_ci.le.n_det_max_jacobi)then
  allocate(eigenvectors(size_super_ci,size_super_ci),eigenvalues(size_super_ci),matrix_tmp(size_super_ci,size_super_ci))
  allocate(vec_tmp(size_super_ci),iorder(size_super_ci))
    print*, ''
    print*, '\\\\\\\\\\\\\\\\\ SUPERCI DIAGONALIZATION /////////////////'
    print*, ''
  do m = 1, N_states
    print*, ''
    print*, 'DIAGONALIZING THE SUPER CI MATRIX FOR STATE ',m
    matrix_tmp(:,:) = superci_matrix(:,:,m)
    call lapack_diag(eigenvalues,eigenvectors,matrix_tmp,size_super_ci,size_super_ci)
    eigenvectors_sci(:,m) = eigenvectors(:,1)
    eigenvalues_sci(m) = eigenvalues(1)
    print*, 'SCI eigenvalue = ',eigenvalues_sci(m)
!   do i = 1, size_super_ci
!    write(*,'(100(F16.10,X))')matrix_tmp(i,:)
!   enddo
  
!!!!!!!!!!!!! SC2 loop 
    logical :: do_sc2
    do_sc2 = .True.
    double precision :: amplitude_max
    integer :: imax
    integer, allocatable :: iorder(:)
    double precision, allocatable :: vec_tmp(:)
    print*, '***'
    do i = 1, size_super_ci
     vec_tmp(i) = -dabs(eigenvectors(i,1)/eigenvectors(1,1))
     iorder(i) = i
    enddo
    call dsort(vec_tmp,iorder,size_super_ci)
    write(*,'(A20,2X,F10.5)')'Largest amplitude = ',vec_tmp(2)
    if(iorder(1).ne.1)then
     do_sc2 = .False.
    endif
    if(dabs(vec_tmp(2)).gt.0.4d0)then
     do_sc2 = .False.
    endif
    if (.not. do_sc2)then
     print*, 'THE SUPERCI WAVE FUNCTION IS TOO MULTI-REFERENCE TO DO SC2 ...'
    endif
    if (do_sc2)then
     print*, 'DOING a few SC2 iterations ...'
     double precision :: e_before
     e_before = eigenvalues_sci(m)
     do iter = 1, 5
      call give_superci_sc2_dressing(delta_H_array,eigenvectors_sci(1,m),m)
      do i = 1, size_super_ci
        matrix_tmp(i,i) = diagonal_superci_matrix(i,m) + delta_H_array(i)  
      enddo
      call lapack_diag(eigenvalues,eigenvectors,matrix_tmp,size_super_ci,size_super_ci)
      eigenvectors_sci(:,m) = eigenvectors(:,1)
      eigenvalues_sci(m) = eigenvalues(1)
      print*, 'SCI(SC)2 energy= ',eigenvalues_sci(m)
      if(dabs(eigenvalues_sci(m) - e_before).lt.1.d-5)then
       exit
      endif 
     enddo
    endif
  enddo

!!!!!!!!!!!!!  
 else 
  allocate(grd_st_eigenvec(size_super_ci))
  do m = 1, N_states
    call create_good_guess(grd_st_eigenvec,e_guess,m)
    call davidson_diag_general_state_specific(grd_st_eigenvec,e_guess,size_super_ci,size_super_ci,m,1,N_int,6)
    eigenvectors_sci(:,m) = grd_st_eigenvec(:)
    eigenvalues_sci(m) = e_guess
  enddo
  deallocate(grd_st_eigenvec)
 endif
 
 deallocate(delta_H_array)
END_PROVIDER 

subroutine give_superci_sc2_dressing(delta_H_array,eigenvector,i_state)
 implicit none
 integer, intent(in)           :: i_state
 double precision, intent(in)  :: eigenvector(size_super_ci)
 double precision, intent(out) :: delta_H_array(size_super_ci)
 integer :: i
 integer :: iorb,jorb,aorb,borb,j,b,index_cj
 double precision :: dsqrt_2,inv_c0,coef_i
 dsqrt_2 = dsqrt(2.d0)
 inv_c0 = 1.d0/eigenvector(1)
 delta_H_array = 0.d0
 do i = 2, size_super_ci
  iorb = index_rotation_CI_reverse(i,1)
  aorb = index_rotation_CI_reverse(i,2)
  coef_i = eigenvector(i) * inv_c0
  do b = 1, n_virt_orb
   if(b== aorb)cycle
   borb = list_virt(b)
   do j = 1, n_core_inact_orb
    if(j== iorb)cycle
    jorb = list_core_inact(j)
    index_cj = index_rotation_CI(j,b)
    delta_H_array(i) += dsqrt_2 * MR_Fock_matrix_alpha_beta_spin_average_mo(jorb,borb,i_state) * eigenvector(index_cj)*inv_c0
!   delta_H_array(1) += eigenvector(index_cj)*inv_c0 * eigenvector(i) * inv_c0 * transformed_occ1_virt1_occ2_virt2(j,b,iorb,aorb)
   enddo
  enddo
 enddo
end 


BEGIN_PROVIDER [double precision, diagonal_superci_matrix_state_average, (size_super_ci)]
 implicit none
 integer :: i,j
 diagonal_superci_matrix_state_average = 0.d0
 do i = 1, N_states
  do j= 1, size_super_ci
   diagonal_superci_matrix_state_average(j) += diagonal_superci_matrix(j,i) * state_average_weight(i)
  enddo
 enddo

END_PROVIDER 

BEGIN_PROVIDER [double precision, superci_matrix_state_average, (size_super_ci,size_super_ci)]
 implicit none
 integer :: i,j,k
 superci_matrix_state_average = 0.d0
 print*, 'PROVIDING THE superci_matrix_state_average'
 do i = 1, N_states
   do j = 1, size_super_ci
    do k = 1, size_super_ci
     superci_matrix_state_average(k,j) += superci_matrix(k,j,i)  * state_average_weight(i)
    enddo
   enddo
 enddo
 ! do k = 1, size_super_ci
 !  write(*,'(100(F16.10,X))')superci_matrix_state_average(k,:)
 ! enddo
END_PROVIDER 


 BEGIN_PROVIDER [double precision, eigenvectors_sci_state_average, (size_super_ci)]
&BEGIN_PROVIDER [double precision, eigenvalues_sci_state_average]
 implicit none
 integer :: i,j,iter
 double precision :: e_guess
 double precision, allocatable :: grd_st_eigenvec(:),eigenvectors(:,:),eigenvalues(:),matrix_tmp(:,:)
 double precision :: u_dot_v
 double precision, allocatable :: delta_H_array(:)
 allocate(delta_H_array(size_super_ci))

 if(size_super_ci.le.n_det_max_jacobi)then
  allocate(eigenvectors(size_super_ci,size_super_ci),eigenvalues(size_super_ci),matrix_tmp(size_super_ci,size_super_ci))
  allocate(vec_tmp(size_super_ci),iorder(size_super_ci))
    print*, ''
    print*, '\\\\\\\\\\\\\\\\\ SUPERCI DIAGONALIZATION /////////////////'
    print*, ''
    print*, ''
    matrix_tmp(:,:) = superci_matrix_state_average(:,:)
    call lapack_diag(eigenvalues,eigenvectors,matrix_tmp,size_super_ci,size_super_ci)
    eigenvectors_sci_state_average(:) = eigenvectors(:,1)
    eigenvalues_sci_state_average = eigenvalues(1)
    print*, 'SCI eigenvalue = ',eigenvalues_sci_state_average
 !  do i = 1, size_super_ci
 !   write(*,'(100(F16.10,X))')matrix_tmp(i,:)
 !  enddo
  
!!!!!!!!!!!!! SC2 loop 
    logical :: do_sc2
    do_sc2 = .True.
    double precision :: amplitude_max
    integer :: imax
    integer, allocatable :: iorder(:)
    double precision, allocatable :: vec_tmp(:)
    print*, '***'
    do i = 1, size_super_ci
     vec_tmp(i) = -dabs(eigenvectors(i,1)/eigenvectors(1,1))
     iorder(i) = i
    enddo
    call dsort(vec_tmp,iorder,size_super_ci)
    write(*,'(A20,2X,F10.5)')'Largest amplitude = ',vec_tmp(2)
    if(iorder(1).ne.1)then
     do_sc2 = .False.
    endif
    if(dabs(vec_tmp(2)).gt.0.4d0)then
     do_sc2 = .False.
    endif
    if (.not. do_sc2)then
     print*, 'THE SUPERCI WAVE FUNCTION IS TOO MULTI-REFERENCE TO DO SC2 ...'
    endif
    if (do_sc2)then
     print*, 'DOING a few SC2 iterations ...'
     double precision :: e_before
     e_before = eigenvalues_sci_state_average
     do iter = 1, 5
      call give_superci_sc2_dressing_state_average(delta_H_array,eigenvectors_sci_state_average)
      do i = 1, size_super_ci
!       print*, diagonal_superci_matrix_state_average(i), delta_H_array(i),eigenvectors_sci_state_average(i)
        matrix_tmp(i,i) = diagonal_superci_matrix_state_average(i) + delta_H_array(i)  
      enddo
      call lapack_diag(eigenvalues,eigenvectors,matrix_tmp,size_super_ci,size_super_ci)
      eigenvectors_sci_state_average(:) = eigenvectors(:,1)
      eigenvalues_sci_state_average = eigenvalues(1)
      print*, 'SCI(SC)2 energy= ',eigenvalues_sci_state_average
      if(dabs(eigenvalues_sci_state_average - e_before).lt.1.d-5)then
       exit
      endif 
     enddo
    endif
!!!!!!!!!!!!!  
 else 
  allocate(grd_st_eigenvec(size_super_ci))
    call create_good_guess_state_average(grd_st_eigenvec,e_guess)
    call davidson_diag_general_state_average(grd_st_eigenvec,e_guess,size_super_ci,size_super_ci,1,N_int,6)
    eigenvectors_sci_state_average(:) = grd_st_eigenvec(:)
    eigenvalues_sci_state_average = e_guess
  deallocate(grd_st_eigenvec)
 endif
 
 deallocate(delta_H_array)
END_PROVIDER 

subroutine give_superci_sc2_dressing_state_average(delta_H_array,eigenvector)
 implicit none
 double precision, intent(in)  :: eigenvector(size_super_ci)
 double precision, intent(out) :: delta_H_array(size_super_ci)
 integer :: i
 integer :: iorb,jorb,aorb,borb,j,b,index_cj
 double precision :: dsqrt_2,inv_c0,coef_i
 dsqrt_2 = dsqrt(2.d0)
 inv_c0 = 1.d0/eigenvector(1)
 delta_H_array = 0.d0
 do i = 2, size_super_ci
  iorb = index_rotation_CI_reverse(i,1)
  aorb = index_rotation_CI_reverse(i,2)
  coef_i = eigenvector(i) * inv_c0
! print*, MR_Fock_matrix_spin_and_state_average_mo(list_core_inact(iorb),list_virt(aorb))
  do b = 1, n_virt_orb
   if(b== aorb)cycle
   borb = list_virt(b)
   do j = 1, n_core_inact_orb
    if(j== iorb)cycle
    jorb = list_core_inact(j)
    index_cj = index_rotation_CI(j,b)
    delta_H_array(i) += dsqrt_2 * MR_Fock_matrix_spin_and_state_average_mo(jorb,borb) * eigenvector(index_cj)*inv_c0
!   delta_H_array(1) += eigenvector(index_cj)*inv_c0 * eigenvector(i) * inv_c0 * transformed_occ1_virt1_occ2_virt2(j,b,iorb,aorb)
   enddo
  enddo
 enddo
end 
