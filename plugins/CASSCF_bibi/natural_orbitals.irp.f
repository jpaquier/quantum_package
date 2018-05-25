BEGIN_PROVIDER [double precision, super_ci_state_specific_density_matrix_mo, (mo_tot_num,mo_tot_num,N_states)]
 implicit none
 integer :: i,j,iorb,jorb,index_ci,index_cj,k,korb,l,lorb,m
 double precision :: dsqrt_2
 dsqrt_2 = dsqrt(2.d0)

 super_ci_state_specific_density_matrix_mo = 0.d0

 do m = 1, N_states
   do i = 1,n_core_inact_orb
    iorb = list_core_inact(i)
    super_ci_state_specific_density_matrix_mo(iorb,iorb,m) = 2.d0
    do j = 1, n_virt_orb
     jorb = list_virt(j)
     index_ci = index_rotation_CI(i,j)
    ! DIAGONAL PART OF THE OCCUPIED 
     super_ci_state_specific_density_matrix_mo(iorb,iorb,m) -= eigenvectors_sci(index_ci,m) **2 
    ! OCC-VIRT PART
     super_ci_state_specific_density_matrix_mo(iorb,jorb,m) += eigenvectors_sci(1,m) * eigenvectors_sci(index_ci,m) * dsqrt_2
     super_ci_state_specific_density_matrix_mo(jorb,iorb,m) += eigenvectors_sci(1,m) * eigenvectors_sci(index_ci,m) * dsqrt_2
    enddo
   enddo
  
   ! DIAGONAL PART OF THE VIRTUAL 
   do i = 1, n_virt_orb
    iorb = list_virt(i)
    do j = 1, n_core_inact_orb
     index_ci = index_rotation_CI(j,i)
     super_ci_state_specific_density_matrix_mo(iorb,iorb,m) += eigenvectors_sci(index_ci,m) **2
    enddo
   enddo
  
   ! EXTRA DIAGONAL PART OF THE OCCUPIED 
   do i = 1,n_core_inact_orb
    iorb = list_core_inact(i)
    do j = i+1, n_core_inact_orb
     jorb = list_core_inact(j)
     do k = 1, n_virt_orb
      index_ci = index_rotation_CI(i,k)
      index_cj = index_rotation_CI(j,k)
      super_ci_state_specific_density_matrix_mo(iorb,jorb,m) -= eigenvectors_sci(index_ci,m) * eigenvectors_sci(index_cj,m)
      super_ci_state_specific_density_matrix_mo(jorb,iorb,m) -= eigenvectors_sci(index_ci,m) * eigenvectors_sci(index_cj,m)
     enddo
    enddo
   enddo
  
   ! EXTRA DIAGONAL PART OF THE VIRTUAL 
   do i = 1, n_virt_orb
    iorb = list_virt(i)
    do j = i+1, n_virt_orb
     jorb = list_virt(j)
     do k = 1, n_core_inact_orb
      index_ci = index_rotation_CI(k,i)
      index_cj = index_rotation_CI(k,j)
      super_ci_state_specific_density_matrix_mo(iorb,jorb,m) += eigenvectors_sci(index_ci,m) * eigenvectors_sci(index_cj,m)
      super_ci_state_specific_density_matrix_mo(jorb,iorb,m) += eigenvectors_sci(index_ci,m) * eigenvectors_sci(index_cj,m)
     enddo
    enddo 
   enddo
  
   ! ACTIVE PART OF THE DENSITY MATRIX
   do i = 1, n_act_orb
    iorb = list_act(i)
    do j = 1, n_act_orb
     jorb = list_act(j)
     super_ci_state_specific_density_matrix_mo(iorb,jorb,m) = one_body_dm_mo_alpha(iorb,jorb,m) + one_body_dm_mo_beta(iorb,jorb,m)
    enddo
   enddo

  enddo

END_PROVIDER 


BEGIN_PROVIDER [double precision, averaged_super_ci_state_specific_density_matrix_mo, (mo_tot_num,mo_tot_num)]
  implicit none
  integer                        :: i
  double precision :: n_elec_states(N_states)
  double precision :: n_elec_average
  integer :: j,k
  averaged_super_ci_state_specific_density_matrix_mo = 0.d0
  do i = 1, N_states
   do j = 1, mo_tot_num
    do k = 1, mo_tot_num
     averaged_super_ci_state_specific_density_matrix_mo(k,j) += super_ci_state_specific_density_matrix_mo(k,j,i) * state_average_weight(i)
    enddo
   enddo
  enddo
  n_elec_states = 0.d0
  n_elec_average = 0.d0
  do j = 1, mo_tot_num
   n_elec_average += averaged_super_ci_state_specific_density_matrix_mo(j,j)
   do i = 1, N_states
    n_elec_states(i) += super_ci_state_specific_density_matrix_mo(j,j,i)
   enddo
  enddo
  print*, 'Number of electrons using state specific DM'
  do i = 1,N_States
   print*, 'state  |  N_elec '
   write(*,'(I3,7X,F10.5)')i,n_elec_states(i)
  enddo
  print*, 'Number of electrons using state average  DM'
   print*, n_elec_average
END_PROVIDER 




BEGIN_PROVIDER [double precision, super_ci_state_average_density_matrix_mo, (mo_tot_num,mo_tot_num)]
 implicit none
 integer :: i,j,iorb,jorb,index_ci,index_cj,k,korb,l,lorb
 double precision :: dsqrt_2,n_elec_average
 dsqrt_2 = dsqrt(2.d0)

 super_ci_state_average_density_matrix_mo = 0.d0

   do i = 1,n_core_inact_orb
    iorb = list_core_inact(i)
    super_ci_state_average_density_matrix_mo(iorb,iorb) = 2.d0
    do j = 1, n_virt_orb
     jorb = list_virt(j)
     index_ci = index_rotation_CI(i,j)
    ! DIAGONAL PART OF THE OCCUPIED 
     super_ci_state_average_density_matrix_mo(iorb,iorb) -= eigenvectors_sci_state_average(index_ci) **2 
    ! OCC-VIRT PART
     super_ci_state_average_density_matrix_mo(iorb,jorb) += eigenvectors_sci_state_average(1) * eigenvectors_sci_state_average(index_ci) * dsqrt_2
     super_ci_state_average_density_matrix_mo(jorb,iorb) += eigenvectors_sci_state_average(1) * eigenvectors_sci_state_average(index_ci) * dsqrt_2
    enddo
   enddo
  
   ! DIAGONAL PART OF THE VIRTUAL 
   do i = 1, n_virt_orb
    iorb = list_virt(i)
    do j = 1, n_core_inact_orb
     index_ci = index_rotation_CI(j,i)
     super_ci_state_average_density_matrix_mo(iorb,iorb) += eigenvectors_sci_state_average(index_ci) **2
    enddo
   enddo
  
   ! EXTRA DIAGONAL PART OF THE OCCUPIED 
   do i = 1,n_core_inact_orb
    iorb = list_core_inact(i)
    do j = i+1, n_core_inact_orb
     jorb = list_core_inact(j)
     do k = 1, n_virt_orb
      index_ci = index_rotation_CI(i,k)
      index_cj = index_rotation_CI(j,k)
      super_ci_state_average_density_matrix_mo(iorb,jorb) -= eigenvectors_sci_state_average(index_ci) * eigenvectors_sci_state_average(index_cj)
      super_ci_state_average_density_matrix_mo(jorb,iorb) -= eigenvectors_sci_state_average(index_ci) * eigenvectors_sci_state_average(index_cj)
     enddo
    enddo
   enddo
  
   ! EXTRA DIAGONAL PART OF THE VIRTUAL 
   do i = 1, n_virt_orb
    iorb = list_virt(i)
    do j = i+1, n_virt_orb
     jorb = list_virt(j)
     do k = 1, n_core_inact_orb
      index_ci = index_rotation_CI(k,i)
      index_cj = index_rotation_CI(k,j)
      super_ci_state_average_density_matrix_mo(iorb,jorb) += eigenvectors_sci_state_average(index_ci) * eigenvectors_sci_state_average(index_cj)
      super_ci_state_average_density_matrix_mo(jorb,iorb) += eigenvectors_sci_state_average(index_ci) * eigenvectors_sci_state_average(index_cj)
     enddo
    enddo 
   enddo
  
   ! ACTIVE PART OF THE DENSITY MATRIX
   do i = 1, n_act_orb
    iorb = list_act(i)
    do j = 1, n_act_orb
     jorb = list_act(j)
     super_ci_state_average_density_matrix_mo(iorb,jorb) = one_body_dm_mo_alpha_average(iorb,jorb) + one_body_dm_mo_beta_average(iorb,jorb)
    enddo
   enddo

  n_elec_average = 0.d0
  do i = 1, mo_tot_num
   n_elec_average += super_ci_state_average_density_matrix_mo(i,i)
  enddo


  print*, 'Number of electrons using state average  DM'
  print*, n_elec_average

END_PROVIDER 



subroutine set_superci_natural_mos_state_specific
 implicit none
 BEGIN_DOC
 ! Set natural orbitals, obtained by diagonalization of the one-body density matrix in the MO basis of the SUPERCI wave function
 END_DOC
 character*(64) :: label
 
 label = "Natural"
 call mo_as_svd_vectors_of_mo_matrix(averaged_super_ci_state_specific_density_matrix_mo,size(averaged_super_ci_state_specific_density_matrix_mo,1),mo_tot_num,mo_tot_num,label)

end

subroutine set_superci_natural_mos_state_average
 implicit none
 BEGIN_DOC
 ! Set natural orbitals, obtained by diagonalization of the one-body density matrix in the MO basis of the SUPERCI wave function
 END_DOC
 character*(64) :: label
 
 label = "Natural"
 call mo_as_svd_vectors_of_mo_matrix(super_ci_state_average_density_matrix_mo,size(super_ci_state_average_density_matrix_mo,1),mo_tot_num,mo_tot_num,label)

end

