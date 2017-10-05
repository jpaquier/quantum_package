BEGIN_PROVIDER [double precision, super_ci_density_matrix_mo, (mo_tot_num,mo_tot_num)]
 implicit none
 integer :: i,j,iorb,jorb,index_ci,index_cj,k,korb,l,lorb
 double precision :: dsqrt_2
 dsqrt_2 = dsqrt(2.d0)

 super_ci_density_matrix_mo = 0.d0

 do i = 1,n_core_inact_orb
  iorb = list_core_inact(i)
  super_ci_density_matrix_mo(iorb,iorb) = 2.d0
  do j = 1, n_virt_orb
   jorb = list_virt(j)
   index_ci = index_rotation_CI(i,j)
  ! DIAGONAL PART OF THE OCCUPIED 
   super_ci_density_matrix_mo(iorb,iorb) -= eigenvectors_sci(index_ci,1) **2 
  ! OCC-VIRT PART
   super_ci_density_matrix_mo(iorb,jorb) += eigenvectors_sci(1,1) * eigenvectors_sci(index_ci,1) * dsqrt_2
   super_ci_density_matrix_mo(jorb,iorb) += eigenvectors_sci(1,1) * eigenvectors_sci(index_ci,1) * dsqrt_2
  enddo
 enddo

 ! DIAGONAL PART OF THE VIRTUAL 
 do i = 1, n_virt_orb
  iorb = list_virt(i)
  do j = 1, n_core_inact_orb
   index_ci = index_rotation_CI(j,i)
   super_ci_density_matrix_mo(iorb,iorb) += eigenvectors_sci(index_ci,1) **2
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
    super_ci_density_matrix_mo(iorb,jorb) -= eigenvectors_sci(index_ci,1) * eigenvectors_sci(index_cj,1)
    super_ci_density_matrix_mo(jorb,iorb) -= eigenvectors_sci(index_ci,1) * eigenvectors_sci(index_cj,1)
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
    super_ci_density_matrix_mo(iorb,jorb) += eigenvectors_sci(index_ci,1) * eigenvectors_sci(index_cj,1)
    super_ci_density_matrix_mo(jorb,iorb) += eigenvectors_sci(index_ci,1) * eigenvectors_sci(index_cj,1)
   enddo
  enddo 
 enddo

 ! ACTIVE PART OF THE DENSITY MATRIX
 do i = 1, n_act_orb
  iorb = list_act(i)
  do j = 1, n_act_orb
   jorb = list_act(j)
   super_ci_density_matrix_mo(iorb,jorb) = one_body_dm_mo_alpha_average(iorb,jorb) + one_body_dm_mo_beta_average(iorb,jorb)
  enddo
 enddo

END_PROVIDER 

subroutine set_superci_natural_mos
 implicit none
 BEGIN_DOC
 ! Set natural orbitals, obtained by diagonalization of the one-body density matrix in the MO basis of the SUPERCI wave function
 END_DOC
 character*(64) :: label
 
 label = "Natural"
 call mo_as_svd_vectors_of_mo_matrix(super_ci_density_matrix_mo,size(one_body_dm_mo,1),mo_tot_num,mo_tot_num,label)

end
