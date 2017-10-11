program pouet
 implicit none
 read_wf = .True.
 touch read_wf
!call test_MR_Fock_matrix_alpha_mo
!call test_superci_eigenvectors
!call test_natural_orbitals 
 call test_ref_energy
!call test_guess_superci 
!call test_fock_spin_average
!call test_H_superci
!call test_fock_from_act

end 

subroutine test_superci_eigenvectors
 implicit none
 print*, eigenvalues_sci_brutal(1),eigenvalues_sci(1)
end

subroutine test_fock_spin_average
 implicit none
!provide MR_Fock_matrix_alpha_beta_spin_average_mo
!provide superci_matrix
 provide transformed_occ1_virt2_virt2

end

subroutine test_fock
 implicit none
 integer :: i,j,iorb,jorb,m
 do m = 1, N_states
  print*, '\\\\\\\\\\\\\\\\\\\'
  print*, 'State ',m
 !print*, 'Reference energy = ',reference_energy_superci(m)
 !print*, 'eigenvalues',eigenvalues_sci(m),eigenvalues_sci_brutal(m)
  do i = 1, size_super_ci
   do j = 1, size_super_ci 
    if(dabs(superci_matrix(i,j,m)-total_matrix_superci_brutal(i,j,m)).gt.1.d-10)then
     print*, ' ***** ',i,j
     print*, index_rotation_CI_reverse(i,1),index_rotation_CI_reverse(i,2)
     print*, index_rotation_CI_reverse(j,1),index_rotation_CI_reverse(j,2)
     print*, superci_matrix(i,j,m),total_matrix_superci_brutal(i,j,m),dabs(superci_matrix(i,j,m)-total_matrix_superci_brutal(i,j,m))
     print*, '***'
    endif
   enddo
  enddo
  print*, '///////////////////'
 enddo

 

end 


subroutine test_H_superci
 implicit none
 double precision, allocatable :: u1(:),u0(:)
 allocate(u1(size_super_ci),u0(size_super_ci))
 double precision, allocatable :: v1(:),v0(:)
 allocate(v1(size_super_ci),v0(size_super_ci))
 integer :: i,j,m
 do m = 1, N_states
   u0 = 0.d0
   u1 = 0.d0
   u0(1) = 1.d0
  
   v0 = 0.d0
   v1 = 0.d0
   v0(1) = 1.d0
  
   do i = 1, size_super_ci
    u1(i) += u0(i) * superci_matrix(i,i,m) 
    do j = i+1, size_super_ci
     u1(i) += u0(j) * superci_matrix(j,i,m) 
     u1(j) += u0(i) * superci_matrix(j,i,m)
    enddo
   enddo
   
  
   call apply_H_superci_state_specific_to_vector(v0,v1,m)
   
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
    u1(i) += u0(i) * superci_matrix(i,i,m) 
    do j = i+1, size_super_ci
     u1(i) += u0(j) * superci_matrix(j,i,m) 
     u1(j) += u0(i) * superci_matrix(j,i,m)
    enddo
   enddo
   
  
   call apply_H_superci_state_specific_to_vector(v0,v1,m)
   
   do i = 1, size_super_ci
    if(dabs(v1(i) - u1(i)).gt.1.d-12)then
     print*, 'PB !!! '
     print*, i,v1(i),u1(i)
    endif
   enddo
 enddo

 deallocate(u1,u0,v0,v1)
end


subroutine test_guess_superci
 implicit none
 provide eigenvalues_sci 
end

subroutine test_ref_energy
 implicit none
 integer :: m
 do m = 1, N_states
  print*,reference_energy_superci(m),psi_energy(m)+nuclear_repulsion
 enddo
end 

subroutine test_natural_orbitals
 implicit none
 character*(64) :: label
 label = "Natural"
 integer :: i
 double precision :: accu
 accu = 0.d0 
 print*, ''
 print*, ''
 print*, 'SUPER CI DM'
 print*, ''
 do i = 1, mo_tot_num
 !write(*,'(100(F10.7,X))')super_ci_state_average_density_matrix_mo(i,:)
  accu += super_ci_state_average_density_matrix_mo(i,i)
 enddo
 print*, 'Number of electrons ',accu
 print*, 'Actual n elec       ',elec_alpha_num + elec_beta_num
 call mo_as_svd_vectors_of_mo_matrix(super_ci_state_average_density_matrix_mo,size(super_ci_state_average_density_matrix_mo,1),mo_tot_num,mo_tot_num,label)
 print*, ''
 integer :: m
 do m = 1, N_states
 print*, ''
 print*, ''
 print*, 'State ',m
  accu = 0.d0
  do i = 1, mo_tot_num
 ! write(*,'(100(F10.7,X))')density_matrix_mo_act(i,:,1)
   accu += density_matrix_mo_act(i,i,m)
  enddo
 print*, ''
 print*, 'Number of electrons ',accu
 print*, ''
 print*, ''
 enddo
 


end

subroutine test_fock_from_act
 implicit none
 integer :: i,j,iorb,jorb
 integer :: m

 do m = 1, N_states
  do i = 1, n_core_inact_orb
   iorb = list_core_inact(i)
   do j = 1, n_virt_orb
    jorb = list_virt(j)
    if(dabs(MR_Fock_matrix_alpha_from_act_mo_bis(jorb,iorb,m) - MR_Fock_matrix_beta_from_act_mo(iorb,jorb,m)).gt.1.d-10)then
     print*, iorb,jorb
     print*, MR_Fock_matrix_alpha_from_act_mo_bis(jorb,iorb,m),MR_Fock_matrix_beta_from_act_mo(iorb,jorb,m),dabs(MR_Fock_matrix_alpha_from_act_mo_bis(jorb,iorb,m) - MR_Fock_matrix_beta_from_act_mo(iorb,jorb,m))
    endif
   enddo
  enddo
 enddo



end

subroutine test_MR_Fock_matrix_alpha_mo
 implicit none
 provide MR_Fock_matrix_alpha_mo


end
