

subroutine i_H_j_mono_spin_monoelec_dft(key_i,key_j,Nint,spin,hij)
  use bitmasks
  implicit none
  BEGIN_DOC
  ! Returns <i|H|j> where i and j are determinants differing by a single excitation
  END_DOC
  integer, intent(in)            :: Nint, spin
  integer(bit_kind), intent(in)  :: key_i(Nint,2), key_j(Nint,2)
  double precision, intent(out)  :: hij
  
  integer                        :: exc(0:2,2)
  double precision               :: phase

  call get_mono_excitation_spin(key_i(1,spin),key_j(1,spin),exc,phase,Nint)
  integer :: m,p
  m = exc(1,1)
  p = exc(1,2)
  hij = phase * one_e_energy_potential(m,p)
end


double precision function diag_H_mat_elem_monoelec_dft(det_in,Nint)
  use bitmasks
  implicit none
  BEGIN_DOC
  ! Computes <i|H|i>
  END_DOC
  integer,intent(in)             :: Nint
  integer(bit_kind),intent(in)   :: det_in(Nint,2)
  
  integer(bit_kind)              :: hole(Nint,2)
  integer(bit_kind)              :: particle(Nint,2)
  integer                        :: i, nexc(2), ispin
  integer                        :: occ_particle(Nint*bit_kind_size,2)
  integer                        :: occ_hole(Nint*bit_kind_size,2)
  integer(bit_kind)              :: det_tmp(Nint,2)
  integer                        :: na, nb
  
  ASSERT (Nint > 0)
  ASSERT (sum(popcnt(det_in(:,1))) == elec_alpha_num)
  ASSERT (sum(popcnt(det_in(:,2))) == elec_beta_num)
  
  diag_H_mat_elem_monoelec_dft = 0.d0
  
  !call debug_det(det_in,Nint)
  integer                        :: tmp(2)
  !DIR$ FORCEINLINE
  call bitstring_to_list_ab(det_in, occ_particle, tmp, Nint)
  do ispin = 1,2 
   do i = 1, tmp(ispin)
    diag_H_mat_elem_monoelec_dft +=  one_e_energy_potential(occ_particle(i,ispin),occ_particle(i,ispin))
   enddo
  enddo
  
end

