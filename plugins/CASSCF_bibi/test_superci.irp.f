program pouet
 implicit none
 read_wf = .True.
 touch read_wf
 call test_H_superci
end 

subroutine test_fock
 implicit none
 integer :: i,j
 print*, 'Reference energy = ',reference_energy_superci+nuclear_repulsion
  print*, eigenvalues_sci_brutal(1),eigenvalues_sci(1)
 do i = 1, size_super_ci
  do j = 1, size_super_ci 
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


subroutine test_H_superci
 implicit none
 double precision, allocatable :: u1(:),u0(:)
 allocate(u1(size_super_ci),u0(size_super_ci))
 double precision, allocatable :: v1(:),v0(:)
 allocate(v1(size_super_ci),v0(size_super_ci))
 integer :: i,j
 u0 = 0.d0
 u1 = 0.d0
 u0(1) = 1.d0

 v0 = 0.d0
 v1 = 0.d0
 v0(1) = 1.d0

 do i = 1, size_super_ci
  u1(i) += u0(i) * superci_matrix(i,i) 
  do j = i+1, size_super_ci
   u1(i) += u0(j) * superci_matrix(j,i) 
   u1(j) += u0(i) * superci_matrix(j,i)
  enddo
 enddo
 

 call apply_H_superci_to_vector(v0,v1)
 
 do i = 1, size_super_ci
  if(dabs(v1(i) - u1(i)).gt.1.d-12)then
   print*, 'PB !!! '
   print*, i,v1(i),u1(i)
  endif
 enddo

 print*, 'passed the first one '
 pause

 u0 = u1
 u1 = 0.d0
 v0 = v1
 v1 = 0.d0

 do i = 1, size_super_ci
  u1(i) += u0(i) * superci_matrix(i,i) 
  do j = i+1, size_super_ci
   u1(i) += u0(j) * superci_matrix(j,i) 
   u1(j) += u0(i) * superci_matrix(j,i)
  enddo
 enddo
 

 call apply_H_superci_to_vector(v0,v1)
 
 do i = 1, size_super_ci
  if(dabs(v1(i) - u1(i)).gt.1.d-12)then
   print*, 'PB !!! '
   print*, i,v1(i),u1(i)
  endif
 enddo

 deallocate(u1,u0,v0,v1)



end
