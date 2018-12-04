program Dirac_Hartree_Fock
  implicit none
  integer :: i
  BEGIN_DOC
! TODO
  END_DOC

 print*, "2*small_ao_num =",2*small_ao_num
 do i=1,2*dirac_mo_tot_num
  print*,i,eigenvalues_dirac_fock_matrix_C_G_mo(i)
 enddo


end
