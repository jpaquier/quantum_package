 ! For tests only, this is the density matrix calculated with the first 2
 ! columns, easy to compare the L_alpha L_alpha bloc with the non-relativistic 
 ! one but meaningless when using real dirac_mo_coef because the first columns
 ! are eigenvectors corresponding to the continuum
!BEGIN_PROVIDER [complex*16, dirac_SCF_density_matrix_ao, (2*dirac_ao_num,2*dirac_ao_num) ]
!  implicit none
!  BEGIN_DOC
!  ! S^{-1}.P.S^{-1}  where P = C.C^t
!  END_DOC
!  call zgemm('N','C',2*dirac_ao_num,2*dirac_ao_num,elec_num,(1.d0,0.d0), &
!       dirac_mo_coef, size(dirac_mo_coef,1), &
!       dirac_mo_coef, size(dirac_mo_coef,1), (0.d0,0.d0), &
!       dirac_SCF_density_matrix_ao, size(dirac_SCF_density_matrix_ao,1))
!END_PROVIDER


 BEGIN_PROVIDER [complex*16,dirac_mo_coef_electronic, (2*dirac_ao_num,elec_num)]
 implicit none
 BEGIN_DOC
 ! This is the matrix of the eigenvectors of the Fock matrix 
 !corresponding to filled electronic states in the ao basis
 END_DOC
 integer :: i,j,j_plus
 do j =1, elec_num
  do i = 1, 2*dirac_ao_num
   dirac_mo_coef_electronic(i,j) = dirac_mo_coef(i,j + 2*small_ao_num)
  enddo
 enddo
 END_PROVIDER

 BEGIN_PROVIDER [complex*16, dirac_SCF_density_matrix_ao, (2*dirac_ao_num,2*dirac_ao_num) ]
   implicit none
   BEGIN_DOC
   ! S^{-1}.P.S^{-1}  where P = C.C^t
   END_DOC
   call zgemm('N','C',2*dirac_ao_num,2*dirac_ao_num,elec_num,(1.d0,0.d0), &
        dirac_mo_coef_electronic, size(dirac_mo_coef_electronic,1), &
        dirac_mo_coef_electronic, size(dirac_mo_coef_electronic,1), (0.d0,0.d0), &
        dirac_SCF_density_matrix_ao, size(dirac_SCF_density_matrix_ao,1))
 ! call zgemm('N','C',2*dirac_ao_num,2*dirac_ao_num,elec_num,(1.d0,0.d0), &
 !       dirac_mo_coef, size(dirac_mo_coef,1), &
 !       dirac_mo_coef, size(dirac_mo_coef,1), (0.d0,0.d0), &
 !       dirac_SCF_density_matrix_ao, size(dirac_SCF_density_matrix_ao,1))
  END_PROVIDER

