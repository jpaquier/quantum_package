program scf
  BEGIN_DOC
! Produce `Hartree_Fock` MO orbital 
! output: mo_basis.mo_tot_num mo_basis.mo_label mo_basis.ao_md5 mo_basis.mo_coef mo_basis.mo_occ
! output: hartree_fock.energy
! optional: mo_basis.mo_coef
  END_DOC
  call orthonormalize_mos
  call run
end
subroutine run

  use bitmasks
  implicit none
  BEGIN_DOC
! Run SCF calculation
  END_DOC
  double precision               :: SCF_energy_before,SCF_energy_after,diag_H_mat_elem
  double precision               :: E0
  integer                        :: i_it, i, j, k
   
  E0 = RS_KS_energy 

  mo_label = "Canonical"
  call damping_SCF
  
 write(*, '(A22,X,F16.10)') 'one_electron_energy = ',one_electron_energy
 write(*, '(A22,X,F16.10)') 'two_electron_energy = ',two_electron_energy
 write(*, '(A22,X,F16.10)') 'e_exchange_dft      = ',e_exchange_dft
 write(*, '(A22,X,F16.10)') 'e_correlation_dft   = ',e_correlation_dft
 write(*, '(A22,X,F16.10)') 'Fock_matrix_energy  = ',Fock_matrix_energy
  
end
