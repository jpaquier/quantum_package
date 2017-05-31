
subroutine i_H_j_dft_general(key_i,key_j,Nint,hij_core, hij_hartree, hij_erf)
  use bitmasks
  implicit none
  BEGIN_DOC
  ! Returns <i|H|j> where i and j are determinants
  END_DOC
  integer, intent(in)            :: Nint
  integer(bit_kind), intent(in)  :: key_i(Nint,2), key_j(Nint,2)
  double precision, intent(out)  :: hij_core, hij_hartree, hij_erf
 
  call i_H_j_erf(key_i,key_j,Nint,hij_erf)
! print*, hij_erf
  call i_H_j_monoelec_dft(key_i,key_j,Nint,hij_core, hij_hartree)
! print*, hij_core, hij_hartree
  
end



