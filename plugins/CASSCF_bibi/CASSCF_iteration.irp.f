subroutine iteration_scf
 implicit none
 integer :: i,j,k,l,iorb,jorb,korb
 integer :: n_orb_rot
 integer :: n_ci
 double precision,allocatable  :: SCI_matrix(:,:), eigvalues(:),eigvectors(:,:)
 integer :: index_rotation_CI(n_core_inact_orb,n_virt_orb)
 integer, allocatable :: index_rotation_CI_reverse(:,:)
 n_orb_rot = n_core_inact_orb*n_virt_orb
 n_ci = n_orb_rot + 1
 allocate(SCI_matrix(n_ci,n_ci), eigvalues(n_ci),eigvectors(n_ci,n_ci),index_rotation_CI_reverse(n_ci,2))
 k = 1
 do i = 1, n_core_inact_orb
  do j = 1, n_virt_orb
   k+=1
   index_rotation_CI(i,j) = k
   index_rotation_CI_reverse(k,1) = i  ! core_inact
   index_rotation_CI_reverse(k,2) = j  ! virt 
  enddo
 enddo
 print*, k
 if(k.ne.n_orb_rot+1)then
  print*, 'PB !! k.ne.n_orb_rot'
  print*, k,n_orb_rot
 endif
 SCI_matrix = 0.d0
 ! Diagonal and Brillouin matrix elements 
 do i = 1, n_core_inact_orb
  iorb = list_core_inact(i)
  do j = 1, n_virt_orb
   jorb = list_virt(j)
   SCI_matrix(index_rotation_CI(i,j),index_rotation_CI(i,j)) = - Fock_matrix_mo(iorb,iorb) + Fock_matrix_mo(jorb,jorb)
   SCI_matrix(1,index_rotation_CI(i,j)) = Fock_matrix_mo(iorb,jorb)
   SCI_matrix(index_rotation_CI(i,j),1) = Fock_matrix_mo(jorb,iorb)
   do k = j+1, n_virt_orb
    korb = list_virt(k)
    SCI_matrix(index_rotation_CI(i,j),index_rotation_CI(i,k)) = Fock_matrix_mo(jorb,korb)
    SCI_matrix(index_rotation_CI(i,k),index_rotation_CI(i,j)) = Fock_matrix_mo(jorb,korb)
   enddo
   do k = i+1, n_core_inact_orb
    korb = list_core_inact(k)
    SCI_matrix(index_rotation_CI(i,j),index_rotation_CI(k,j)) = Fock_matrix_mo(iorb,korb)
    SCI_matrix(index_rotation_CI(k,j),index_rotation_CI(i,j)) = Fock_matrix_mo(iorb,korb)
   enddo
  enddo
 enddo
 do i = 1, n_ci
  print*, SCI_matrix(i,i),SCI_matrix(1,i)
 enddo


end
