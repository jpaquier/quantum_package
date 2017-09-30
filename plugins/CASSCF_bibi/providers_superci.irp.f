
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
 do i = 1, n_core_orb
  do j = 1, n_virt_orb 
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
 do i = 1, n_core_inact_orb
  do j = 1, n_virt_orb
   k+=1
   index_rotation_CI(i,j) = k
   index_rotation_CI_reverse(k,1) = i  ! core_inact
   index_rotation_CI_reverse(k,2) = j  ! virt 
  enddo
 enddo
END_PROVIDER 



BEGIN_PROVIDER [double precision, superci_matrix, (size_super_ci,size_super_ci)]
 implicit none
 integer :: i,iorb,j,jorb,k,korb,l,lorb
 double precision :: dsqrt_2
 dsqrt_2 = dsqrt(2.d0)
 
 logical  :: simple_fock_matrix 
 simple_fock_matrix = .False.
 if(simple_fock_matrix)then
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
!  do i = 1, size_super_ci
!   print*, superci_matrix(i,i),superci_matrix(1,i)
!  enddo

 else 

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



 BEGIN_PROVIDER [double precision, eigenvectors_sci, (size_super_ci,size_super_ci)]
&BEGIN_PROVIDER [double precision, eigenvalues_sci, (size_super_ci)]
 implicit none
 integer :: i,j

 call lapack_diag(eigenvalues_sci,eigenvectors_sci,superci_matrix,size_super_ci,size_super_ci)
!do i = 1, size_super_ci
! print*, eigenvalues_sci(i)
!enddo


END_PROVIDER 
