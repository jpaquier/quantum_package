 BEGIN_PROVIDER [complex*16,dirac_mo_coef_electronic, (2*dirac_ao_num,elec_num)]
 implicit none
 BEGIN_DOC
 !This is the matrix of the eigenvectors of the Fock matrix 
 ! corresponding to filled electronic states in the ao basis
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
   !Density matrix in the AOs
   ! C.C* 
   ! with C the dirac_mo_coef_electronic 
   END_DOC
   complex*16, Allocatable    :: dirac_SCF_density_matrix_ao_tmp(:,:)
   Allocate (dirac_SCF_density_matrix_ao_tmp(2*dirac_ao_num,2*dirac_ao_num))
   call zgemm('N','C',2*dirac_ao_num,2*dirac_ao_num,elec_num,(1.d0,0.d0), &
        dirac_mo_coef_electronic, size(dirac_mo_coef_electronic,1), &
        dirac_mo_coef_electronic, size(dirac_mo_coef_electronic,1), (0.d0,0.d0), &
        dirac_SCF_density_matrix_ao_tmp, size(dirac_SCF_density_matrix_ao_tmp,1))
  dirac_SCF_density_matrix_ao = (dirac_SCF_density_matrix_ao_tmp)
  deallocate(dirac_SCF_density_matrix_ao_tmp)
 END_PROVIDER

