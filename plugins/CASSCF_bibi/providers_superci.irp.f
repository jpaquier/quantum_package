
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


BEGIN_PROVIDER [double precision, diagonal_superci_matrix, (size_super_ci)]
 implicit none
 integer :: i,iorb,j,jorb,k,korb,l,lorb
 diagonal_superci_matrix = 0.d0
 print*, 'type_of_superci ',type_of_superci
 do i = 1, n_core_inact_orb
  iorb = list_core_inact(i)
  do j = 1, n_virt_orb
   jorb = list_virt(j)
   if(type_of_superci == 0)then
    diagonal_superci_matrix(index_rotation_CI(i,j)) = - Fock_matrix_alpha_beta_average_mo(iorb,iorb) + Fock_matrix_alpha_beta_average_mo(jorb,jorb) 
   else 
    diagonal_superci_matrix(index_rotation_CI(i,j)) = - Fock_matrix_alpha_beta_average_mo(iorb,iorb) + Fock_matrix_alpha_beta_average_mo(jorb,jorb) & 
                                                      - transformed_occ1_virt2_virt2(i,j,j)          + 2.d0 * transformed_occ1_virt1_occ2_virt2(i,j,i,j)
   endif
  enddo
 enddo


END_PROVIDER 


BEGIN_PROVIDER [integer, type_of_superci]
 implicit none
 type_of_superci = 2


END_PROVIDER 

BEGIN_PROVIDER [double precision, superci_matrix, (size_super_ci,size_super_ci)]
 implicit none
 integer :: i,iorb,j,jorb,k,korb,l,lorb
 double precision :: dsqrt_2
 dsqrt_2 = dsqrt(2.d0)
 superci_matrix = 0.d0
 print*, 'Providing the superci matrix'
 
 if(type_of_superci == 0)then
   do i = 1, n_core_inact_orb
    iorb = list_core_inact(i)
    do j = 1, n_virt_orb
     jorb = list_virt(j)
     ! Diagonal and Brillouin matrix elements 
     superci_matrix(index_rotation_CI(i,j),index_rotation_CI(i,j)) =  - Fock_matrix_alpha_beta_average_mo(iorb,iorb) + Fock_matrix_alpha_beta_average_mo(jorb,jorb) 
     superci_matrix(1,index_rotation_CI(i,j)) = dsqrt_2 * Fock_matrix_alpha_beta_average_mo(iorb,jorb)
     superci_matrix(index_rotation_CI(i,j),1) = dsqrt_2 * Fock_matrix_alpha_beta_average_mo(jorb,iorb)
     ! Interaction through the virt-virt Fock operator
     do k = j+1, n_virt_orb
      korb = list_virt(k)
      superci_matrix(index_rotation_CI(i,j),index_rotation_CI(i,k)) = Fock_matrix_alpha_beta_average_mo(jorb,korb)
      superci_matrix(index_rotation_CI(i,k),index_rotation_CI(i,j)) = Fock_matrix_alpha_beta_average_mo(jorb,korb)
     enddo
     ! Interaction through the core-core Fock operator
     do k = i+1, n_core_inact_orb
      korb = list_core_inact(k)
      superci_matrix(index_rotation_CI(i,j),index_rotation_CI(k,j)) = - Fock_matrix_alpha_beta_average_mo(iorb,korb)
      superci_matrix(index_rotation_CI(k,j),index_rotation_CI(i,j)) = - Fock_matrix_alpha_beta_average_mo(iorb,korb)
     enddo
    enddo
   enddo

 else if(type_of_superci == 1)then


  do i = 1, n_core_inact_orb
   iorb = list_core_inact(i)
   do j = 1, n_virt_orb
    jorb = list_virt(j)
    ! Diagonal and Brillouin matrix elements 
    superci_matrix(index_rotation_CI(i,j),index_rotation_CI(i,j)) = -Fock_matrix_alpha_beta_average_mo(iorb,iorb) + Fock_matrix_alpha_beta_average_mo(jorb,jorb) & 
                                                                   - transformed_occ1_virt2_virt2(i,j,j)          + 2.d0 * transformed_occ1_virt1_occ2_virt2(i,j,i,j)
    superci_matrix(1,index_rotation_CI(i,j)) = dsqrt_2 * Fock_matrix_alpha_beta_average_mo(iorb,jorb)
    superci_matrix(index_rotation_CI(i,j),1) = dsqrt_2 * Fock_matrix_alpha_beta_average_mo(jorb,iorb)
    ! Interaction through the virt-virt Fock operator
    do k = j+1, n_virt_orb
     korb = list_virt(k)
     superci_matrix(index_rotation_CI(i,j),index_rotation_CI(i,k)) = Fock_matrix_alpha_beta_average_mo(jorb,korb)  & 
                                                                   - transformed_occ1_virt2_virt2(i,k,j) + 2.d0 * transformed_occ1_virt1_occ2_virt2(i,k,i,j)
     superci_matrix(index_rotation_CI(i,k),index_rotation_CI(i,j)) = superci_matrix(index_rotation_CI(i,j),index_rotation_CI(i,k))
    enddo
    ! Interaction through the core-core Fock operator
    do k = i+1, n_core_inact_orb
     korb = list_core_inact(k)
     superci_matrix(index_rotation_CI(i,j),index_rotation_CI(k,j)) = - Fock_matrix_alpha_beta_average_mo(iorb,korb) & 
                                                                     - transformed_virt1_occ2_occ2(j,k,i) + 2.d0 * transformed_occ1_virt1_occ2_virt2(i,j,k,j)
     superci_matrix(index_rotation_CI(k,j),index_rotation_CI(i,j)) = superci_matrix(index_rotation_CI(i,j),index_rotation_CI(k,j))
    enddo
    
   enddo
  enddo
 else if(type_of_superci == 2)then

  do i = 1, n_core_inact_orb
   iorb = list_core_inact(i)
   do j = 1, n_virt_orb
    jorb = list_virt(j)
    ! Diagonal and Brillouin matrix elements 
    superci_matrix(index_rotation_CI(i,j),index_rotation_CI(i,j)) = -Fock_matrix_alpha_beta_average_mo(iorb,iorb) + Fock_matrix_alpha_beta_average_mo(jorb,jorb) & 
                                                                   - transformed_occ1_virt2_virt2(i,j,j)          + 2.d0 * transformed_occ1_virt1_occ2_virt2(i,j,i,j)
    superci_matrix(1,index_rotation_CI(i,j)) = dsqrt_2 * Fock_matrix_alpha_beta_average_mo(iorb,jorb)
    superci_matrix(index_rotation_CI(i,j),1) = dsqrt_2 * Fock_matrix_alpha_beta_average_mo(jorb,iorb)
    ! Interaction through the virt-virt Fock operator
    do k = j+1, n_virt_orb
     korb = list_virt(k)
     superci_matrix(index_rotation_CI(i,j),index_rotation_CI(i,k)) = Fock_matrix_alpha_beta_average_mo(jorb,korb)  & 
                                                                   - transformed_occ1_virt2_virt2(i,k,j) + 2.d0 * transformed_occ1_virt1_occ2_virt2(i,k,i,j)
     superci_matrix(index_rotation_CI(i,k),index_rotation_CI(i,j)) = superci_matrix(index_rotation_CI(i,j),index_rotation_CI(i,k))
    enddo
    ! Interaction through the core-core Fock operator
    do k = i+1, n_core_inact_orb
     korb = list_core_inact(k)
     superci_matrix(index_rotation_CI(i,j),index_rotation_CI(k,j)) = - Fock_matrix_alpha_beta_average_mo(iorb,korb) & 
                                                                     - transformed_virt1_occ2_occ2(j,k,i) + 2.d0 * transformed_occ1_virt1_occ2_virt2(i,j,k,j)
     superci_matrix(index_rotation_CI(k,j),index_rotation_CI(i,j)) = superci_matrix(index_rotation_CI(i,j),index_rotation_CI(k,j))
    enddo
    ! Hole-particle interaction 
    do l = 1, j-1
     do k = 1, i-1
      superci_matrix(index_rotation_CI(i,j),index_rotation_CI(k,l)) = 2.d0 * transformed_occ1_virt1_occ2_virt2(k,l,i,j) - transformed_virt1_virt1_occ2_occ2(l,j,k,i)
     enddo
    enddo

    do l = 1, j-1
     do k = i+1, n_core_inact_orb
      superci_matrix(index_rotation_CI(i,j),index_rotation_CI(k,l)) = 2.d0 * transformed_occ1_virt1_occ2_virt2(k,l,i,j) - transformed_virt1_virt1_occ2_occ2(l,j,k,i)
     enddo
    enddo

    do l = j+1, n_virt_orb
     do k = 1, i-1
      superci_matrix(index_rotation_CI(i,j),index_rotation_CI(k,l)) = 2.d0 * transformed_occ1_virt1_occ2_virt2(k,l,i,j) - transformed_virt1_virt1_occ2_occ2(l,j,k,i)
     enddo
    enddo

    do l = j+1, n_virt_orb
     do k = i+1, n_core_inact_orb
      superci_matrix(index_rotation_CI(i,j),index_rotation_CI(k,l)) = 2.d0 * transformed_occ1_virt1_occ2_virt2(k,l,i,j) - transformed_virt1_virt1_occ2_occ2(l,j,k,i)
     enddo
    enddo
    
   enddo
  enddo

 endif

END_PROVIDER 


 BEGIN_PROVIDER [double precision, eigenvectors_sci, (size_super_ci,1)]
&BEGIN_PROVIDER [double precision, eigenvalues_sci, (1)]
 implicit none
 integer :: i,j
 double precision :: e_guess
 double precision, allocatable :: grd_st_eigenvec(:),eigenvectors(:,:),eigenvalues(:)
 double precision :: u_dot_v

 if(size_super_ci.le.n_det_max_jacobi)then
  provide superci_matrix
  allocate(eigenvectors(size_super_ci,size_super_ci),eigenvalues(size_super_ci))
  call lapack_diag(eigenvalues,eigenvectors,superci_matrix,size_super_ci,size_super_ci)
  eigenvectors_sci(:,1) = eigenvectors(:,1)
  eigenvalues_sci(1) = eigenvalues(1)
 else 
  allocate(grd_st_eigenvec(size_super_ci))
  call create_good_guess(grd_st_eigenvec,e_guess)
  print*, 'e_guess = ',e_guess
  call davidson_diag_general(grd_st_eigenvec,e_guess,size_super_ci,size_super_ci,1,1,N_int,output_determinants)
  eigenvectors_sci(:,1) = grd_st_eigenvec(:)
  eigenvalues_sci(1) = e_guess
  deallocate(grd_st_eigenvec)
 endif
END_PROVIDER 

