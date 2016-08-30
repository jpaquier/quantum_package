subroutine apply_exc_to_psi(orb,hole_particle,spin_exc, & 
           norm_out,psi_in_out,psi_in_out_coef, ndet,dim_psi_in,dim_psi_coef,N_states_in)
  use bitmasks
 implicit none
 integer, intent(in) :: orb, hole_particle,spin_exc,N_states_in,ndet,dim_psi_in,dim_psi_coef
 double precision, intent(out)  :: norm_out(N_states_in)
 integer(bit_kind), intent(inout) :: psi_in_out(N_int,2,dim_psi_in)
 double precision,  intent(inout) :: psi_in_out_coef(dim_psi_coef,N_states_in)
  BEGIN_DOC
  ! apply a contracted excitation to psi_in_out whose coefficients 
  ! are psi_in_out_coef
  ! hole_particle =  1  ===> creation     of an electron in psi_in_out
  !               = -1  ===> annhilation  of an electron in psi_in_out
  ! orb ===> is the index of orbital where you want wether to create or 
  !          annhilate an electron
  ! spin_exc ===> is the spin of the electron (1 == alpha) (2 == beta)
  ! the wave function gets out normalized to unity
  !
  ! norm_out is the sum of the squared of the coefficients 
  ! on which the excitation has been possible
  END_DOC

  integer :: elec_num_tab_local(2)
  elec_num_tab_local = 0
  do i = 1, ndet
   if( psi_in_out_coef (i,1) .ne. 0.d0)then
    do j = 1, N_int
     elec_num_tab_local(1) += popcnt(psi_in_out(j,1,i))
     elec_num_tab_local(2) += popcnt(psi_in_out(j,2,i))
    enddo
    exit
   endif
  enddo
  integer :: i,j,accu_elec
  if(hole_particle == 1)then
   do i = 1, ndet
     call set_bit_to_integer(orb,psi_in_out(1,spin_exc,i),N_int) 
     accu_elec = 0
     do j = 1, N_int
      accu_elec += popcnt(psi_in_out(j,spin_exc,i))
     enddo 
     if(accu_elec .ne. elec_num_tab_local(spin_exc)+1)then
      do j = 1, N_int
       psi_in_out(j,1,i) = 0_bit_kind
       psi_in_out(j,2,i) = 0_bit_kind
      enddo
      do j = 1, N_states_in
       psi_in_out_coef(i,j) = 0.d0
      enddo
     endif
   enddo
  else if (hole_particle == -1)then
   do i = 1, ndet
     call clear_bit_to_integer(orb,psi_in_out(1,spin_exc,i),N_int) 
     accu_elec = 0
     do j = 1, N_int
      accu_elec += popcnt(psi_in_out(j,spin_exc,i))
     enddo 
     if(accu_elec .ne. elec_num_tab_local(spin_exc)-1)then
      do j = 1, N_int
       psi_in_out(j,1,i) = 0_bit_kind
       psi_in_out(j,2,i) = 0_bit_kind
      enddo
      do j = 1, N_states_in
       psi_in_out_coef(i,j) = 0.d0
      enddo
     endif
   enddo
  endif
  norm_out = 0.d0
  double precision :: norm_factor
  do j = 1, N_states_in
   do i = 1, ndet
    norm_out(j) += psi_in_out_coef(i,j) * psi_in_out_coef(i,j)
   enddo
   if(norm_out(j).le.1.d-10)then
    norm_factor = 0.d0
   else 
    norm_factor = 1.d0/(dsqrt(norm_out(j)))
   endif
   do i = 1, ndet
    psi_in_out_coef(i,j) = psi_in_out_coef(i,j) * norm_factor
   enddo
  enddo
end


double precision function diag_H_mat_elem_no_elec_check(det_in,Nint)
  implicit none
  BEGIN_DOC
  ! Computes <i|H|i>
  END_DOC
  integer,intent(in)             :: Nint
  integer(bit_kind),intent(in)   :: det_in(Nint,2)
  
  integer                        :: i, j, iorb, jorb 
  integer                        :: occ(Nint*bit_kind_size,2)
  integer                        ::  elec_num_tab_local(2)

  diag_H_mat_elem_no_elec_check = 0.d0
  call bitstring_to_list(det_in(1,1), occ(1,1), elec_num_tab_local(1), N_int)
  call bitstring_to_list(det_in(1,2), occ(1,2), elec_num_tab_local(2), N_int)
  ! alpha - alpha 
  do i = 1, elec_num_tab_local(1)
   iorb =  occ(i,1)
   diag_H_mat_elem_no_elec_check += mo_mono_elec_integral(iorb,iorb)
   do j = i+1, elec_num_tab_local(1)
    jorb = occ(j,1)
    diag_H_mat_elem_no_elec_check +=  mo_bielec_integral_jj_anti(jorb,iorb)
   enddo
  enddo 

  ! beta - beta   
  do i = 1, elec_num_tab_local(2)
   iorb =  occ(i,2)
   diag_H_mat_elem_no_elec_check += mo_mono_elec_integral(iorb,iorb)
   do j = i+1, elec_num_tab_local(2)
    jorb = occ(j,2)
    diag_H_mat_elem_no_elec_check +=  mo_bielec_integral_jj_anti(jorb,iorb)
   enddo
  enddo 
  

  ! alpha - beta   
  do i = 1, elec_num_tab_local(2)
   iorb =  occ(i,2)
   do j = 1, elec_num_tab_local(1)
    jorb = occ(j,1)
    diag_H_mat_elem_no_elec_check +=  mo_bielec_integral_jj(jorb,iorb)
   enddo
  enddo 
  
end

subroutine a_operator_no_check(iorb,ispin,key,hjj,Nint,na,nb)
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
  
  hjj = hjj - mo_mono_elec_integral(iorb,iorb)
  
  ! Same spin
  do i=1,na
    hjj = hjj - mo_bielec_integral_jj_anti(occ(i,ispin),iorb)
  enddo
  
  ! Opposite spin
  do i=1,nb
    hjj = hjj - mo_bielec_integral_jj(occ(i,other_spin),iorb)
  enddo
  
end


subroutine ac_operator_no_check(iorb,ispin,key,hjj,Nint,na,nb)
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
  
  k = ishft(iorb-1,-bit_kind_shift)+1
  ASSERT (k > 0)
  l = iorb - ishft(k-1,bit_kind_shift)-1
  key(k,ispin) = ibset(key(k,ispin),l)
  other_spin = iand(ispin,1)+1
  
  hjj = hjj + mo_mono_elec_integral(iorb,iorb)
  
  print*,'na.nb = ',na,nb
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


subroutine i_H_j_dyall(key_i,key_j,Nint,hij)
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
  double precision               :: get_mo_bielec_integral_schwartz
  integer                        :: m,n,p,q
  integer                        :: i,j,k
  integer                        :: occ(Nint*bit_kind_size,2)
  double precision               :: diag_H_mat_elem_no_elec_check, phase,phase_2
  integer                        :: n_occ_ab(2)
  logical                        :: has_mipi(Nint*bit_kind_size)
  double precision               :: mipi(Nint*bit_kind_size), miip(Nint*bit_kind_size)
  PROVIDE mo_bielec_integrals_in_map mo_integrals_map
  
  ASSERT (Nint > 0)
  ASSERT (Nint == N_int)
  
  hij = 0.d0
  !DIR$ FORCEINLINE
  call get_excitation_degree(key_i,key_j,degree,Nint)
  select case (degree)
    case (2)
      call get_double_excitation(key_i,key_j,exc,phase,Nint)
      if (exc(0,1,1) == 1) then
        ! Mono alpha, mono beta
        hij = phase*get_mo_bielec_integral_schwartz(                          &
            exc(1,1,1),                                              &
            exc(1,1,2),                                              &
            exc(1,2,1),                                              &
            exc(1,2,2) ,mo_integrals_map)
      else if (exc(0,1,1) == 2) then
        ! Double alpha
        hij = phase*(get_mo_bielec_integral_schwartz(                         &
            exc(1,1,1),                                              &
            exc(2,1,1),                                              &
            exc(1,2,1),                                              &
            exc(2,2,1) ,mo_integrals_map) -                          &
            get_mo_bielec_integral_schwartz(                                  &
            exc(1,1,1),                                              &
            exc(2,1,1),                                              &
            exc(2,2,1),                                              &
            exc(1,2,1) ,mo_integrals_map) )
      else if (exc(0,1,2) == 2) then
        ! Double beta
        hij = phase*(get_mo_bielec_integral_schwartz(                         &
            exc(1,1,2),                                              &
            exc(2,1,2),                                              &
            exc(1,2,2),                                              &
            exc(2,2,2) ,mo_integrals_map) -                          &
            get_mo_bielec_integral_schwartz(                                  &
            exc(1,1,2),                                              &
            exc(2,1,2),                                              &
            exc(2,2,2),                                              &
            exc(1,2,2) ,mo_integrals_map) )
      endif
    case (1)
      call get_mono_excitation(key_i,key_j,exc,phase,Nint)
      !DIR$ FORCEINLINE
      call bitstring_to_list_ab(key_i, occ, n_occ_ab, Nint)
      has_mipi = .False.
      if (exc(0,1,1) == 1) then
        ! Mono alpha
        m = exc(1,1,1)
        p = exc(1,2,1)
        do k = 1, n_occ_ab(1)
          i = occ(k,1)
          if (.not.has_mipi(i)) then
            mipi(i) = get_mo_bielec_integral_schwartz(m,i,p,i,mo_integrals_map)
            miip(i) = get_mo_bielec_integral_schwartz(m,i,i,p,mo_integrals_map)
            has_mipi(i) = .True.
          endif
        enddo
        do k = 1, n_occ_ab(2)
          i = occ(k,2)
          if (.not.has_mipi(i)) then
            mipi(i) = get_mo_bielec_integral_schwartz(m,i,p,i,mo_integrals_map)
            has_mipi(i) = .True.
          endif
        enddo
        
        do k = 1, n_occ_ab(1)
          hij = hij + mipi(occ(k,1)) - miip(occ(k,1))
        enddo
        do k = 1, n_occ_ab(2)
          hij = hij + mipi(occ(k,2))
        enddo
        
      else
        ! Mono beta
        m = exc(1,1,2)
        p = exc(1,2,2)
        do k = 1, n_occ_ab(2)
          i = occ(k,2)
          if (.not.has_mipi(i)) then
            mipi(i) = get_mo_bielec_integral_schwartz(m,i,p,i,mo_integrals_map)
            miip(i) = get_mo_bielec_integral_schwartz(m,i,i,p,mo_integrals_map)
            has_mipi(i) = .True.
          endif
        enddo
        do k = 1, n_occ_ab(1)
          i = occ(k,1)
          if (.not.has_mipi(i)) then
            mipi(i) = get_mo_bielec_integral_schwartz(m,i,p,i,mo_integrals_map)
            has_mipi(i) = .True.
          endif
        enddo
        
        do k = 1, n_occ_ab(1)
          hij = hij + mipi(occ(k,1))
        enddo
        do k = 1, n_occ_ab(2)
          hij = hij + mipi(occ(k,2)) - miip(occ(k,2))
        enddo
        
      endif
      hij = phase*(hij + mo_mono_elec_integral(m,p) + fock_operator_active_from_core_inact(m,p) )
      
    case (0)
      hij = diag_H_mat_elem_no_elec_check(key_i,Nint)
  end select
end


subroutine u0_H_dyall_u0(energies,psi_in,psi_in_coef,ndet,dim_psi_in,dim_psi_coef,N_states_in,state_target)
  use bitmasks
 implicit none
 integer, intent(in) :: N_states_in,ndet,dim_psi_in,dim_psi_coef,state_target
 integer(bit_kind), intent(in) :: psi_in(N_int,2,dim_psi_in)
 double precision,  intent(in) :: psi_in_coef(dim_psi_coef,N_states_in)
 double precision,  intent(out) :: energies(N_states_in)
 
 integer :: i,j 
 double precision :: hij,accu
 energies = 0.d0
 accu = 0.d0
 double precision, allocatable :: psi_coef_tmp(:)
 allocate(psi_coef_tmp(ndet))
 
 do i = 1, ndet
  psi_coef_tmp(i) = psi_in_coef(i,state_target)
 enddo

 double precision :: hij_bis
 do i = 1, ndet
  if(psi_coef_tmp(i)==0.d0)cycle
  do j = 1, ndet
   if(psi_coef_tmp(j)==0.d0)cycle
   call i_H_j_dyall(psi_in(1,1,i),psi_in(1,1,j),N_int,hij)
!  call i_H_j(psi_in(1,1,i),psi_in(1,1,j),N_int,hij_bis)
!  print*, hij_bis,hij
   accu += psi_coef_tmp(i) * psi_coef_tmp(j) * hij
  enddo
 enddo
 energies(state_target) = accu
 deallocate(psi_coef_tmp)
end