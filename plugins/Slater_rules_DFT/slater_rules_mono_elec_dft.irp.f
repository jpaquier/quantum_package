
subroutine i_H_j_monoelec_dft_no_vxc(key_i,key_j,Nint,hij_core, hij_hartree)
  use bitmasks
  implicit none
  BEGIN_DOC
  ! Returns hij_hartree = <i|v_{H}^{sr}|i> and hij_core = <i|h_{core}|j> where i and j are determinants differing by a single excitation
  ! and v_{H}^{sr} is the short-range part of the hartree term and h_{core} = T + v_{n-e}
  END_DOC
  integer, intent(in)            :: Nint
  integer(bit_kind), intent(in)  :: key_i(Nint,2), key_j(Nint,2)
  double precision, intent(out)  :: hij_core, hij_hartree
  
  integer                        :: exc(0:2,2,2), degree
  double precision               :: phase

  hij_core = 0.d0
  hij_hartree = 0.d0
  call get_excitation_degree(key_i,key_j,degree,Nint)
  if (degree == 1)then
   call get_mono_excitation(key_i,key_j,exc,phase,Nint)
   call decode_exc(exc,degree,h1,p1,h2,p2,s1,s2)
   integer :: h1,h2,p1,p2,s1,s2
   hij_core = phase * ( mo_kinetic_integral(h1,p1) + mo_nucl_elec_integral(h1,p1) )
   hij_hartree = phase *  short_range_Hartree_operator(h1,p1) 
  else if (degree == 0)then
   call diag_H_mat_elem_monoelec_dft_no_vxc_components(key_i,Nint, hij_core, hij_hartree)
  endif
end

subroutine i_H_j_monoelec_dft_with_vxc(key_i,key_j,Nint,hij_core, hij_hartree, hij_vxc)
  use bitmasks
  implicit none
  BEGIN_DOC
  ! Returns hij_hartree = <i|v_{H}^{sr}|i>, hij_core = <i|h_{core}|j>, hij_vxc = <i|v_{xc}^{sr}|j> where i and j are determinants differing by a single excitation
  ! v_{H}^{sr} is the short-range part of the hartree term, h_{core} = T + v_{n-e}, and v_{xc}^{sr} is the short-range exchange-correlation potential
  END_DOC
  integer, intent(in)            :: Nint
  integer(bit_kind), intent(in)  :: key_i(Nint,2), key_j(Nint,2)
  double precision, intent(out)  :: hij_core, hij_hartree,hij_vxc
  
  integer                        :: exc(0:2,2,2), degree
  double precision               :: phase

  hij_core = 0.d0
  hij_hartree = 0.d0
  call get_excitation_degree(key_i,key_j,degree,Nint)
  if (degree == 1)then
   call get_mono_excitation(key_i,key_j,exc,phase,Nint)
   call decode_exc(exc,degree,h1,p1,h2,p2,s1,s2)
   integer :: h1,h2,p1,p2,s1,s2
   hij_core    = phase * ( mo_kinetic_integral(h1,p1) + mo_nucl_elec_integral(h1,p1) )
   hij_hartree = phase *   short_range_Hartree_operator(h1,p1) 
   hij_vxc     = phase *   0.5d0 * (potential_x_alpha_mo(h1,p1,1) + potential_c_alpha_mo(h1,p1,1)                               &
                                   + potential_x_beta_mo(h1,p1,1) + potential_c_beta_mo(h1,p1,1)   )

  else if (degree == 0)then
   call diag_H_mat_elem_monoelec_dft_with_vxc_components(key_i,Nint, hij_core, hij_hartree,hij_vxc)
  endif
end



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

subroutine diag_H_mat_elem_monoelec_dft_no_vxc_components(det_in,Nint, hij_core, hij_hartree)
  use bitmasks
  implicit none
  BEGIN_DOC
  ! Computes <i|H|i>
  END_DOC
  integer,intent(in)             :: Nint
  integer(bit_kind),intent(in)   :: det_in(Nint,2)
  double precision, intent(out)  :: hij_core, hij_hartree
  
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
  
  hij_core = 0.d0
  hij_hartree = 0.d0
  
  !call debug_det(det_in,Nint)
  integer                        :: tmp(2)
  !DIR$ FORCEINLINE
  call bitstring_to_list_ab(det_in, occ_particle, tmp, Nint)
  do ispin = 1,2 
   do i = 1, tmp(ispin)
    hij_core +=  mo_kinetic_integral(occ_particle(i,ispin),occ_particle(i,ispin)) + mo_nucl_elec_integral(occ_particle(i,ispin),occ_particle(i,ispin))
    hij_hartree +=  0.5d0 * short_range_Hartree_operator(occ_particle(i,ispin),occ_particle(i,ispin))
   enddo
  enddo
  
end


subroutine diag_H_mat_elem_monoelec_dft_with_vxc_components(det_in,Nint, hij_core, hij_hartree, hij_vxc)
  use bitmasks
  implicit none
  BEGIN_DOC
  ! Computes <i|H|i>
  END_DOC
  integer,intent(in)             :: Nint
  integer(bit_kind),intent(in)   :: det_in(Nint,2)
  double precision, intent(out)  :: hij_core, hij_hartree, hij_vxc
  
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
  
  hij_core = 0.d0
  hij_hartree = 0.d0
  hij_vxc = 0.d0
  
  !call debug_det(det_in,Nint)
  integer                        :: tmp(2)
  !DIR$ FORCEINLINE
  call bitstring_to_list_ab(det_in, occ_particle, tmp, Nint)
  do ispin = 1,2 
   do i = 1, tmp(ispin)
    hij_core    +=  mo_kinetic_integral(occ_particle(i,ispin),occ_particle(i,ispin)) + mo_nucl_elec_integral(occ_particle(i,ispin),occ_particle(i,ispin))
    hij_hartree +=  0.5d0 * short_range_Hartree_operator(occ_particle(i,ispin),occ_particle(i,ispin))
    hij_vxc     +=  0.5d0 * (potential_x_alpha_mo(i,i,1)+ potential_c_alpha_mo(i,i,1)+ potential_x_beta_mo(i,i,1)+ potential_c_beta_mo(i,i,1))
   enddo
  enddo
  
end

subroutine i_H_j_mono(key_i,key_j,Nint,hij)
  use bitmasks
  implicit none
  BEGIN_DOC
  ! Returns <i|h_{core}|j> where i and j are determinants
  ! and h_{core} is the sum of the kinetic operator and the nuclear-electron attraction
  END_DOC
  integer, intent(in)            :: Nint
  integer(bit_kind), intent(in)  :: key_i(Nint,2), key_j(Nint,2)
  double precision, intent(out)  :: hij
  
  integer                        :: exc(0:2,2,2)
  integer                        :: degree
  integer                        :: m,n,p,q
  integer                        :: i,j,k
  integer                        :: occ(Nint*bit_kind_size,2)
  double precision               :: diag_H_mat_elem_mono, phase,phase_2
  integer                        :: n_occ_ab(2)
  
  ASSERT (Nint > 0)
  ASSERT (Nint == N_int)
  ASSERT (sum(popcnt(key_i(:,1))) == elec_alpha_num)
  ASSERT (sum(popcnt(key_i(:,2))) == elec_beta_num)
  ASSERT (sum(popcnt(key_j(:,1))) == elec_alpha_num)
  ASSERT (sum(popcnt(key_j(:,2))) == elec_beta_num)
  
  hij = 0.d0
  !DIR$ FORCEINLINE
  call get_excitation_degree(key_i,key_j,degree,Nint)
  integer :: spin
  select case (degree)
    case (1)
      call get_mono_excitation(key_i,key_j,exc,phase,Nint)
      if (exc(0,1,1) == 1) then
        ! Mono alpha
        m = exc(1,1,1)
        p = exc(1,2,1)
      else
        ! Mono beta
        m = exc(1,1,2)
        p = exc(1,2,2)
      endif
      hij = phase * (mo_nucl_elec_integral(p,m) + mo_kinetic_integral(p,m))
    case (0)
      hij = diag_H_mat_elem_mono(key_i,Nint)
  end select
end

double precision function diag_H_mat_elem_mono(key_i,Nint)
 implicit none
 integer(bit_kind), intent(in) :: key_i(N_int,2)
 integer, intent(in)  :: Nint
 integer :: i,j
 integer                        :: occ(Nint*bit_kind_size,2)
 integer                        :: n_occ_ab(2)
 call bitstring_to_list_ab(key_i, occ, n_occ_ab, Nint)
 diag_H_mat_elem_mono = 0.d0
 ! alpha - alpha
 do i = 1, n_occ_ab(1)
  diag_H_mat_elem_mono += (mo_nucl_elec_integral(occ(i,1),occ(i,1)) + mo_kinetic_integral(occ(i,1),occ(i,1)))
 enddo

 ! beta - beta 
 do i = 1, n_occ_ab(2)
  diag_H_mat_elem_mono += (mo_nucl_elec_integral(occ(i,2),occ(i,2)) + mo_kinetic_integral(occ(i,2),occ(i,2)))
 enddo

end


subroutine i_H_j_hartree(key_i,key_j,Nint,hij)
  use bitmasks
  implicit none
  BEGIN_DOC
  ! Returns <i|H|j> where i and j are determinants
  END_DOC
  integer, intent(in)            :: Nint
  integer(bit_kind), intent(in)  :: key_i(Nint,2), key_j(Nint,2)
  double precision, intent(out)  :: hij
  
  integer                        :: exc(0:2,2,2)
  integer                        :: degree
  integer                        :: m,n,p,q
  integer                        :: i,j,k
  integer                        :: occ(Nint*bit_kind_size,2)
  double precision               :: diag_H_mat_elem_hartree, phase,phase_2
  integer                        :: n_occ_ab(2)
  
  ASSERT (Nint > 0)
  ASSERT (Nint == N_int)
  ASSERT (sum(popcnt(key_i(:,1))) == elec_alpha_num)
  ASSERT (sum(popcnt(key_i(:,2))) == elec_beta_num)
  ASSERT (sum(popcnt(key_j(:,1))) == elec_alpha_num)
  ASSERT (sum(popcnt(key_j(:,2))) == elec_beta_num)
  
  hij = 0.d0
  !DIR$ FORCEINLINE
  call get_excitation_degree(key_i,key_j,degree,Nint)
  integer :: spin
  select case (degree)
    case (1)
      call get_mono_excitation(key_i,key_j,exc,phase,Nint)
      if (exc(0,1,1) == 1) then
        ! Mono alpha
        m = exc(1,1,1)
        p = exc(1,2,1)
      else
        ! Mono beta
        m = exc(1,1,2)
        p = exc(1,2,2)
      endif
      hij = phase * (short_range_Hartree_operator(p,m) + short_range_Hartree_operator(p,m))
    case (0)
      hij = diag_H_mat_elem_hartree(key_i,Nint)
  end select
end




double precision function diag_H_mat_elem_hartree(key_i,Nint)
 implicit none
 integer(bit_kind), intent(in) :: key_i(N_int,2)
 integer, intent(in)  :: Nint
 integer :: i,j
 integer                        :: occ(Nint*bit_kind_size,2)
 integer                        :: n_occ_ab(2)
 call bitstring_to_list_ab(key_i, occ, n_occ_ab, Nint)
 diag_H_mat_elem_hartree = 0.d0
 ! alpha - alpha
 do i = 1, n_occ_ab(1)
  diag_H_mat_elem_hartree += (short_range_Hartree_operator(occ(i,1),occ(i,1)) + short_range_Hartree_operator(occ(i,1),occ(i,1)))
 enddo

 ! beta - beta 
 do i = 1, n_occ_ab(2)
  diag_H_mat_elem_hartree += (short_range_Hartree_operator(occ(i,2),occ(i,2)) + short_range_Hartree_operator(occ(i,2),occ(i,2)))
 enddo

end



