

subroutine i_H_j_mono_spin(key_i,key_j,Nint,spin,hij)
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

  PROVIDE big_array_exchange_integrals mo_bielec_integrals_in_map

  call get_mono_excitation_spin(key_i(1,spin),key_j(1,spin),exc,phase,Nint)
  call get_mono_excitation_from_fock(key_i,key_j,exc(1,1),exc(1,2),spin,phase,hij)
end

subroutine i_H_j_double_spin(key_i,key_j,Nint,hij)
  use bitmasks
  implicit none
  BEGIN_DOC
  ! Returns <i|H|j> where i and j are determinants differing by a same-spin double excitation
  END_DOC
  integer, intent(in)            :: Nint
  integer(bit_kind), intent(in)  :: key_i(Nint), key_j(Nint)
  double precision, intent(out)  :: hij
  
  integer                        :: exc(0:2,2)
  double precision               :: phase
  double precision, external     :: get_mo_bielec_integral

  PROVIDE big_array_exchange_integrals mo_bielec_integrals_in_map

  call get_double_excitation_spin(key_i,key_j,exc,phase,Nint)
  hij = phase*(get_mo_bielec_integral(                               &
      exc(1,1),                                                      &
      exc(2,1),                                                      &
      exc(1,2),                                                      &
      exc(2,2), mo_integrals_map) -                                  &
      get_mo_bielec_integral(                                        &
      exc(1,1),                                                      &
      exc(2,1),                                                      &
      exc(2,2),                                                      &
      exc(1,2), mo_integrals_map) )
end

subroutine i_H_j_double_alpha_beta(key_i,key_j,Nint,hij)
  use bitmasks
  implicit none
  BEGIN_DOC
  ! Returns <i|H|j> where i and j are determinants differing by an opposite-spin double excitation
  END_DOC
  integer, intent(in)            :: Nint
  integer(bit_kind), intent(in)  :: key_i(Nint,2), key_j(Nint,2)
  double precision, intent(out)  :: hij
  
  integer                        :: exc(0:2,2,2)
  double precision               :: phase, phase2
  double precision, external     :: get_mo_bielec_integral

  PROVIDE big_array_exchange_integrals mo_bielec_integrals_in_map

  call get_mono_excitation_spin(key_i(1,1),key_j(1,1),exc(0,1,1),phase,Nint)
  call get_mono_excitation_spin(key_i(1,2),key_j(1,2),exc(0,1,2),phase2,Nint)
  phase = phase*phase2
  if (exc(1,1,1) == exc(1,2,2)) then
    hij = phase * big_array_exchange_integrals(exc(1,1,1),exc(1,1,2),exc(1,2,1))
  else if (exc(1,2,1) == exc(1,1,2)) then
    hij = phase * big_array_exchange_integrals(exc(1,2,1),exc(1,1,1),exc(1,2,2))
  else
    hij = phase*get_mo_bielec_integral(                              &
        exc(1,1,1),                                                  &
        exc(1,1,2),                                                  &
        exc(1,2,1),                                                  &
        exc(1,2,2) ,mo_integrals_map)
  endif
end



subroutine i_H_j_mono_spin_bielec(key_i,key_j,Nint,spin,hij)
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

  PROVIDE big_array_exchange_integrals mo_bielec_integrals_in_map

  call get_mono_excitation_spin(key_i(1,spin),key_j(1,spin),exc,phase,Nint)
  call get_mono_excitation_from_fock_bielec(key_i,key_j,exc(1,1),exc(1,2),spin,phase,hij)
end


double precision function diag_H_mat_elem_bielec(det_in,Nint)
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
  
  nexc(1) = 0
  nexc(2) = 0
  do i=1,Nint
    hole(i,1)     = xor(det_in(i,1),ref_bitmask(i,1))
    hole(i,2)     = xor(det_in(i,2),ref_bitmask(i,2))
    particle(i,1) = iand(hole(i,1),det_in(i,1))
    particle(i,2) = iand(hole(i,2),det_in(i,2))
    hole(i,1)     = iand(hole(i,1),ref_bitmask(i,1))
    hole(i,2)     = iand(hole(i,2),ref_bitmask(i,2))
    nexc(1)       = nexc(1) + popcnt(hole(i,1))
    nexc(2)       = nexc(2) + popcnt(hole(i,2))
  enddo
  
  diag_H_mat_elem_bielec = bi_elec_ref_bitmask_energy
  if (nexc(1)+nexc(2) == 0) then
    return
  endif
  
  !call debug_det(det_in,Nint)
  integer                        :: tmp(2)
  !DIR$ FORCEINLINE
  call bitstring_to_list_ab(particle, occ_particle, tmp, Nint)
  ASSERT (tmp(1) == nexc(1))
  ASSERT (tmp(2) == nexc(2))
  !DIR$ FORCEINLINE
  call bitstring_to_list_ab(hole, occ_hole, tmp, Nint)
  ASSERT (tmp(1) == nexc(1))
  ASSERT (tmp(2) == nexc(2))
  
  det_tmp = ref_bitmask
  do ispin=1,2
    na = elec_num_tab(ispin)
    nb = elec_num_tab(iand(ispin,1)+1)
    do i=1,nexc(ispin)
      !DIR$ FORCEINLINE
      call ac_operator_bielec( occ_particle(i,ispin), ispin, det_tmp, diag_H_mat_elem_bielec, Nint,na,nb)
      !DIR$ FORCEINLINE
      call a_operator_bielec ( occ_hole    (i,ispin), ispin, det_tmp, diag_H_mat_elem_bielec, Nint,na,nb)
    enddo
  enddo
end


subroutine a_operator_bielec(iorb,ispin,key,hjj,Nint,na,nb)
  use bitmasks
  implicit none
  BEGIN_DOC
  ! Needed for diag_H_mat_elem
  END_DOC
  integer, intent(in)            :: iorb, ispin, Nint
  integer, intent(inout)         :: na, nb
  integer(bit_kind), intent(inout) :: key(Nint,2)
  double precision, intent(inout) :: hjj
  
  integer                        :: occ(Nint*bit_kind_size,2)
  integer                        :: other_spin
  integer                        :: k,l,i
  integer                        :: tmp(2)
  
  ASSERT (iorb > 0)
  ASSERT (ispin > 0)
  ASSERT (ispin < 3)
  ASSERT (Nint > 0)
  
  k = ishft(iorb-1,-bit_kind_shift)+1
  ASSERT (k > 0)
  l = iorb - ishft(k-1,bit_kind_shift)-1
  key(k,ispin) = ibclr(key(k,ispin),l)
  other_spin = iand(ispin,1)+1
  
  !DIR$ FORCEINLINE
  call bitstring_to_list_ab(key, occ, tmp, Nint)
  na = na-1
  
  ! Same spin
  do i=1,na
    hjj = hjj - mo_bielec_integral_jj_anti(occ(i,ispin),iorb)
  enddo
  
  ! Opposite spin
  do i=1,nb
    hjj = hjj - mo_bielec_integral_jj(occ(i,other_spin),iorb)
  enddo
  
end


subroutine ac_operator_bielec(iorb,ispin,key,hjj,Nint,na,nb)
  use bitmasks
  implicit none
  BEGIN_DOC
  ! Needed for diag_H_mat_elem
  END_DOC
  integer, intent(in)            :: iorb, ispin, Nint
  integer, intent(inout)         :: na, nb
  integer(bit_kind), intent(inout) :: key(Nint,2)
  double precision, intent(inout) :: hjj
  
  integer                        :: occ(Nint*bit_kind_size,2)
  integer                        :: other_spin
  integer                        :: k,l,i
  
  ASSERT (iorb > 0)
  ASSERT (ispin > 0)
  ASSERT (ispin < 3)
  ASSERT (Nint > 0)
  
  integer                        :: tmp(2)
  !DIR$ FORCEINLINE
  call bitstring_to_list_ab(key, occ, tmp, Nint)
  ASSERT (tmp(1) == elec_alpha_num)
  ASSERT (tmp(2) == elec_beta_num)
  
  k = ishft(iorb-1,-bit_kind_shift)+1
  ASSERT (k > 0)
  l = iorb - ishft(k-1,bit_kind_shift)-1
  key(k,ispin) = ibset(key(k,ispin),l)
  other_spin = iand(ispin,1)+1
  
  
  ! Same spin
  do i=1,na
    hjj = hjj + mo_bielec_integral_jj_anti(occ(i,ispin),iorb)
  enddo
  
  ! Opposite spin
  do i=1,nb
    hjj = hjj + mo_bielec_integral_jj(occ(i,other_spin),iorb)
  enddo
  na = na+1
end



subroutine i_H_j_mono_spin_monoelec(key_i,key_j,Nint,spin,hij)
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
  hij = phase * mo_mono_elec_integral(m,p)
end


double precision function diag_H_mat_elem_monoelec(det_in,Nint)
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
  
  diag_H_mat_elem_monoelec = 0.d0
  
  !call debug_det(det_in,Nint)
  integer                        :: tmp(2)
  !DIR$ FORCEINLINE
  call bitstring_to_list_ab(det_in, occ_particle, tmp, Nint)
  do ispin = 1,2 
   do i = 1, tmp(ispin)
    diag_H_mat_elem_monoelec +=  mo_mono_elec_integral(occ_particle(i,ispin),occ_particle(i,ispin))
   enddo
  enddo
  
end

subroutine i_H_j_monoelec(key_i,key_j,Nint,hij)
  use bitmasks
  implicit none
  BEGIN_DOC
  ! Returns <i|H|j> where i and j are determinants
  END_DOC
  integer, intent(in)            :: Nint
  integer(bit_kind), intent(in)  :: key_i(Nint,2), key_j(Nint,2)
  double precision, intent(out)  :: hij

  integer :: degree,m,p
  double precision :: diag_H_mat_elem_monoelec,phase
  integer                        :: exc(0:2,2,2)
  call get_excitation_degree(key_i,key_j,degree,Nint)
  hij = 0.d0
  if(degree>1)then
   return
  endif
  if(degree==0)then
   hij = diag_H_mat_elem_monoelec(key_i,N_int)
  else 
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
   hij = phase * mo_mono_elec_integral(m,p)
  endif

end

subroutine i_H_j_bielec(key_i,key_j,Nint,hij)
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
  double precision               :: get_mo_bielec_integral
  integer                        :: m,n,p,q
  integer                        :: i,j,k
  integer                        :: occ(Nint*bit_kind_size,2)
  double precision               :: diag_H_mat_elem, phase,phase_2
  integer                        :: n_occ_ab(2)
  PROVIDE mo_bielec_integrals_in_map mo_integrals_map big_array_exchange_integrals bi_elec_ref_bitmask_energy
  
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
    case (2)
      call get_double_excitation(key_i,key_j,exc,phase,Nint)
      if (exc(0,1,1) == 1) then
        ! Mono alpha, mono beta
        if(exc(1,1,1) == exc(1,2,2) )then
         hij = phase * big_array_exchange_integrals(exc(1,1,1),exc(1,1,2),exc(1,2,1))
        else if (exc(1,2,1) ==exc(1,1,2))then
         hij = phase * big_array_exchange_integrals(exc(1,2,1),exc(1,1,1),exc(1,2,2))
        else
         hij = phase*get_mo_bielec_integral(                          &
             exc(1,1,1),                                              &
             exc(1,1,2),                                              &
             exc(1,2,1),                                              &
             exc(1,2,2) ,mo_integrals_map)
        endif
      else if (exc(0,1,1) == 2) then
        ! Double alpha
        hij = phase*(get_mo_bielec_integral(                         &
            exc(1,1,1),                                              &
            exc(2,1,1),                                              &
            exc(1,2,1),                                              &
            exc(2,2,1) ,mo_integrals_map) -                          &
            get_mo_bielec_integral(                                  &
            exc(1,1,1),                                              &
            exc(2,1,1),                                              &
            exc(2,2,1),                                              &
            exc(1,2,1) ,mo_integrals_map) )
      else if (exc(0,1,2) == 2) then
        ! Double beta
        hij = phase*(get_mo_bielec_integral(                         &
            exc(1,1,2),                                              &
            exc(2,1,2),                                              &
            exc(1,2,2),                                              &
            exc(2,2,2) ,mo_integrals_map) -                          &
            get_mo_bielec_integral(                                  &
            exc(1,1,2),                                              &
            exc(2,1,2),                                              &
            exc(2,2,2),                                              &
            exc(1,2,2) ,mo_integrals_map) )
      endif
    case (1)
      call get_mono_excitation(key_i,key_j,exc,phase,Nint)
      !DIR$ FORCEINLINE
      call bitstring_to_list_ab(key_i, occ, n_occ_ab, Nint)
      if (exc(0,1,1) == 1) then
        ! Mono alpha
        m = exc(1,1,1)
        p = exc(1,2,1)
        spin = 1
      else
        ! Mono beta
        m = exc(1,1,2)
        p = exc(1,2,2)
        spin = 2
      endif
      call get_mono_excitation_from_fock_bielec(key_i,key_j,p,m,spin,phase,hij)
    case (0)
      double precision :: diag_H_mat_elem_bielec
      hij = diag_H_mat_elem_bielec(key_i,Nint)
  end select
end

