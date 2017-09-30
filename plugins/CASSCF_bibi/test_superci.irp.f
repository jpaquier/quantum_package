program pouet
 implicit none
 read_wf = .True.
 touch read_wf
 call test_fock
end 

subroutine test_fock
 implicit none
 integer :: i,j
 print*, 'Reference energy = ',reference_energy_superci+nuclear_repulsion
  print*, eigenvalues_sci_brutal(1),eigenvalues_sci(1)
 do i = 1, size_super_ci
  do j = 1, size_super_ci 
!  if(dabs(superci_matrix(i,j)) - dabs(total_matrix_superci_brutal(i,j)).gt.1.d-10)then
   if(dabs(superci_matrix(i,j)-total_matrix_superci_brutal(i,j)).gt.1.d-10.and.i.ne.1)then
    print*, ' ***** ',i,j
    print*, index_rotation_CI_reverse(i,1),index_rotation_CI_reverse(i,2)
    print*, index_rotation_CI_reverse(j,1),index_rotation_CI_reverse(j,2)
    print*, superci_matrix(i,j),total_matrix_superci_brutal(i,j)
    print*, '***'
 
   endif
  enddo
 enddo

 

end 
