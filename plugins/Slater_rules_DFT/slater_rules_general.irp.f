
subroutine i_H_j_dft_no_vxc(key_i,key_j,Nint,hij_core, hij_hartree, hij_erf)
  use bitmasks
  implicit none
  BEGIN_DOC
  ! Returns <i| h_{core} + W_{ee}^{lr} + v_{H}^{sr} |j> where i and j are determinants
  ! hij_core    = <i| h_{core}   |j> 
  ! hij_hartree = <i| v_{H}^{sr} |j> 
  ! hij_erf     = <i| W_{ee}^{lr}|j> 
  END_DOC
  integer, intent(in)            :: Nint
  integer(bit_kind), intent(in)  :: key_i(Nint,2), key_j(Nint,2)
  double precision, intent(out)  :: hij_core, hij_hartree, hij_erf
 
  call i_H_j_erf(key_i,key_j,Nint,hij_erf)
  call i_H_j_monoelec_dft_no_vxc(key_i,key_j,Nint,hij_core, hij_hartree)
  
end




subroutine i_H_j_dft_with_vxc(key_i,key_j,Nint,hij_core,hij_hartree,hij_vxc,hij_erf)
  use bitmasks
  implicit none
  BEGIN_DOC
  ! Returns <i| h_{core} + W_{ee}^{lr} + v_{H}^{sr} + v_{xc}^{sr}|j> where i and j are determinants
  ! hij_core    = <i| h_{core}   |j> 
  ! hij_hartree = <i| v_{H}^{sr} |j> 
  ! hij_erf     = <i| W_{ee}^{lr}|j> 
  END_DOC
  integer, intent(in)            :: Nint
  integer(bit_kind), intent(in)  :: key_i(Nint,2), key_j(Nint,2)
  double precision, intent(out)  :: hij_core, hij_hartree, hij_erf, hij_vxc
 
  call i_H_j_erf(key_i,key_j,Nint,hij_erf)
  call i_H_j_monoelec_dft_with_vxc(key_i,key_j,Nint,hij_core, hij_hartree,hij_vxc)
  
end



