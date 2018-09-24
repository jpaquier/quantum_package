!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! THIS FILE CONTAINS EVERYTHING YOU NEED TO COMPUTE THE LONG RANGE PART OF THE INTERACTION
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine i_H_j_erf_mu_of_r(key_i,key_j,Nint,hij)
  use bitmasks
  implicit none
  BEGIN_DOC
  ! Returns <i|W_{ee}^{lr}|j> where i and j are determinants 
  ! and the W_{ee}^{lr} is the long range two-body interaction
  END_DOC
  integer, intent(in)            :: Nint
  integer(bit_kind), intent(in)  :: key_i(Nint,2), key_j(Nint,2)
  double precision, intent(out)  :: hij
  
  integer                        :: exc(0:2,2,2)
  integer                        :: degree
  double precision               :: get_mo_bielec_integral_erf_mu_of_r
  integer                        :: m,n,p,q
  integer                        :: i,j,k
  integer                        :: occ(Nint*bit_kind_size,2)
  double precision               :: diag_H_mat_elem_erf_mu_of_r, phase,phase_2
  integer                        :: n_occ_ab(2)
  PROVIDE mo_bielec_integrals_erf_mu_of_r_in_map mo_integrals_erf_mu_of_r_map big_array_exchange_integrals_erf_mu_of_r
  
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
         hij = phase * big_array_exchange_integrals_erf_mu_of_r(exc(1,1,1),exc(1,1,2),exc(1,2,1))
        else if (exc(1,2,1) ==exc(1,1,2))then
         hij = phase * big_array_exchange_integrals_erf_mu_of_r(exc(1,2,1),exc(1,1,1),exc(1,2,2))
        else
         hij = phase*get_mo_bielec_integral_erf_mu_of_r(                          &
             exc(1,1,1),                                              &
             exc(1,1,2),                                              &
             exc(1,2,1),                                              &
             exc(1,2,2) ,mo_integrals_erf_mu_of_r_map)
        endif
      else if (exc(0,1,1) == 2) then
        ! Double alpha
        hij = phase*(get_mo_bielec_integral_erf_mu_of_r(                         &
            exc(1,1,1),                                              &
            exc(2,1,1),                                              &
            exc(1,2,1),                                              &
            exc(2,2,1) ,mo_integrals_erf_mu_of_r_map) -                          &
            get_mo_bielec_integral_erf_mu_of_r(                                  &
            exc(1,1,1),                                              &
            exc(2,1,1),                                              &
            exc(2,2,1),                                              &
            exc(1,2,1) ,mo_integrals_erf_mu_of_r_map) )
      else if (exc(0,1,2) == 2) then
        ! Double beta
        hij = phase*(get_mo_bielec_integral_erf_mu_of_r(                         &
            exc(1,1,2),                                              &
            exc(2,1,2),                                              &
            exc(1,2,2),                                              &
            exc(2,2,2) ,mo_integrals_erf_mu_of_r_map) -                          &
            get_mo_bielec_integral_erf_mu_of_r(                                  &
            exc(1,1,2),                                              &
            exc(2,1,2),                                              &
            exc(2,2,2),                                              &
            exc(1,2,2) ,mo_integrals_erf_mu_of_r_map) )
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
        do i = 1, n_occ_ab(1)
         hij += -big_array_exchange_integrals_erf_mu_of_r(occ(i,1),m,p) + big_array_coulomb_integrals_erf_mu_of_r(occ(i,1),m,p)
        enddo
        do i = 1, n_occ_ab(2)
         hij += big_array_coulomb_integrals_erf_mu_of_r(occ(i,2),m,p)
        enddo
      else
        ! Mono beta
        m = exc(1,1,2)
        p = exc(1,2,2)
        spin = 2
        do i = 1, n_occ_ab(2)
         hij += -big_array_exchange_integrals_erf_mu_of_r(occ(i,2),m,p) + big_array_coulomb_integrals_erf_mu_of_r(occ(i,2),m,p)
        enddo
        do i = 1, n_occ_ab(1)
         hij += big_array_coulomb_integrals_erf_mu_of_r(occ(i,1),m,p)
        enddo
      endif
      hij = hij * phase
    case (0)
      hij = diag_H_mat_elem_erf_mu_of_r(key_i,Nint)
  end select
end


double precision function diag_H_mat_elem_erf_mu_of_r(key_i,Nint)
 BEGIN_DOC 
! returns <i|W_{ee}^{lr}|i> where |i> is a determinant and 
! W_{ee}^{lr} is the two body long-range interaction
 END_DOC
 implicit none
 integer(bit_kind), intent(in) :: key_i(N_int,2)
 integer, intent(in)  :: Nint
 integer :: i,j
 integer                        :: occ(Nint*bit_kind_size,2)
 integer                        :: n_occ_ab(2)
 call bitstring_to_list_ab(key_i, occ, n_occ_ab, Nint)
 diag_H_mat_elem_erf_mu_of_r = 0.d0
 ! alpha - alpha
 do i = 1, n_occ_ab(1)
  do j = i+1, n_occ_ab(1)
   diag_H_mat_elem_erf_mu_of_r += mo_bielec_integral_erf_mu_of_r_jj_anti(occ(i,1),occ(j,1))
  enddo
 enddo

 ! beta - beta 
 do i = 1, n_occ_ab(2)
  do j = i+1, n_occ_ab(2)
   diag_H_mat_elem_erf_mu_of_r += mo_bielec_integral_erf_mu_of_r_jj_anti(occ(i,2),occ(j,2))
  enddo
 enddo

 ! alpha - beta 
 do i = 1, n_occ_ab(1)
  do j = 1, n_occ_ab(2)
   diag_H_mat_elem_erf_mu_of_r += mo_bielec_integral_erf_mu_of_r_jj(occ(i,1),occ(j,2))
  enddo
 enddo
end







subroutine i_H_j_erf_mu_of_r_and_short_coulomb(key_i,key_j,Nint,hij)
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
  double precision               :: get_mo_bielec_integral_erf_mu_of_r
  integer                        :: m,n,p,q
  integer                        :: i,j,k
  integer                        :: occ(Nint*bit_kind_size,2)
  double precision               :: diag_H_mat_elem_erf_mu_of_r, phase,phase_2
  integer                        :: n_occ_ab(2)
  PROVIDE mo_bielec_integrals_erf_mu_of_r_in_map mo_integrals_erf_mu_of_r_map big_array_exchange_integrals_erf_mu_of_r
  
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
         hij = phase*get_mo_bielec_integral_erf_mu_of_r(                          &
             exc(1,1,1),                                              &
             exc(1,1,2),                                              &
             exc(1,2,1),                                              &
             exc(1,2,2) ,mo_integrals_erf_mu_of_r_map)
        endif
      else if (exc(0,1,1) == 2) then
        ! Double alpha
        hij = phase*(get_mo_bielec_integral_erf_mu_of_r(                         &
            exc(1,1,1),                                              &
            exc(2,1,1),                                              &
            exc(1,2,1),                                              &
            exc(2,2,1) ,mo_integrals_erf_mu_of_r_map) -                          &
            get_mo_bielec_integral_erf_mu_of_r(                                  &
            exc(1,1,1),                                              &
            exc(2,1,1),                                              &
            exc(2,2,1),                                              &
            exc(1,2,1) ,mo_integrals_erf_mu_of_r_map) )
      else if (exc(0,1,2) == 2) then
        ! Double beta
        hij = phase*(get_mo_bielec_integral_erf_mu_of_r(                         &
            exc(1,1,2),                                              &
            exc(2,1,2),                                              &
            exc(1,2,2),                                              &
            exc(2,2,2) ,mo_integrals_erf_mu_of_r_map) -                          &
            get_mo_bielec_integral_erf_mu_of_r(                                  &
            exc(1,1,2),                                              &
            exc(2,1,2),                                              &
            exc(2,2,2),                                              &
            exc(1,2,2) ,mo_integrals_erf_mu_of_r_map) )
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
        do i = 1, n_occ_ab(1)
         hij += -big_array_exchange_integrals_erf_mu_of_r(occ(i,1),m,p) + big_array_coulomb_integrals_erf_mu_of_r(occ(i,1),m,p)
        enddo
        do i = 1, n_occ_ab(2)
         hij += big_array_coulomb_integrals_erf_mu_of_r(occ(i,2),m,p)
        enddo
      else
        ! Mono beta
        m = exc(1,1,2)
        p = exc(1,2,2)
        spin = 2
        do i = 1, n_occ_ab(2)
         hij += -big_array_exchange_integrals_erf_mu_of_r(occ(i,2),m,p) + big_array_coulomb_integrals_erf_mu_of_r(occ(i,2),m,p)
        enddo
        do i = 1, n_occ_ab(1)
         hij += big_array_coulomb_integrals_erf_mu_of_r(occ(i,1),m,p)
        enddo
      endif
      hij = hij + mo_nucl_elec_integral(m,p) + mo_kinetic_integral(m,p) + short_range_Hartree_operator(m,p)
      hij = hij * phase
    case (0)
      hij = diag_H_mat_elem_erf_mu_of_r(key_i,Nint)
  end select
end

double precision function diag_H_mat_elem_erf_mu_of_r_and_short_coulomb(key_i,Nint)
 implicit none
 integer(bit_kind), intent(in) :: key_i(N_int,2)
 integer, intent(in)  :: Nint
 integer :: i,j
 integer                        :: occ(Nint*bit_kind_size,2)
 integer                        :: n_occ_ab(2)

 call bitstring_to_list_ab(key_i, occ, n_occ_ab, Nint)
 diag_H_mat_elem_erf_mu_of_r_and_short_coulomb = 0.d0
 ! alpha - alpha
 do i = 1, n_occ_ab(1)
  diag_H_mat_elem_erf_mu_of_r_and_short_coulomb += mo_nucl_elec_integral(occ(i,1),mo_nucl_elec_integral(i,1)) + mo_kinetic_integral(occ(i,1),mo_nucl_elec_integral(i,1)) &
                                         + short_range_Hartree_operator(occ(i,1),occ(i,1))
  do j = i+1, n_occ_ab(1)
   diag_H_mat_elem_erf_mu_of_r_and_short_coulomb += mo_bielec_integral_erf_mu_of_r_jj_anti(occ(i,1),occ(j,1))
  enddo
 enddo

 ! beta - beta 
 do i = 1, n_occ_ab(2)
  diag_H_mat_elem_erf_mu_of_r_and_short_coulomb += mo_nucl_elec_integral(occ(i,2),mo_nucl_elec_integral(i,2)) + mo_kinetic_integral(occ(i,2),mo_nucl_elec_integral(i,2)) & 
                                         + short_range_Hartree_operator(occ(i,2),occ(i,2))
  do j = i+1, n_occ_ab(2)
   diag_H_mat_elem_erf_mu_of_r_and_short_coulomb += mo_bielec_integral_erf_mu_of_r_jj_anti(occ(i,2),occ(j,2))
  enddo
 enddo

 ! alpha - beta 
 do i = 1, n_occ_ab(1)
  do j = 1, n_occ_ab(2)
   diag_H_mat_elem_erf_mu_of_r_and_short_coulomb += mo_bielec_integral_erf_mu_of_r_jj(occ(i,1),occ(j,2))
  enddo
 enddo

end


subroutine i_H_j_erf_mu_of_r_component(key_i,key_j,Nint,hij_core,hij_hartree,hij_erf_mu_of_r,hij_total)
  use bitmasks
  implicit none
  BEGIN_DOC
  ! Returns <i|H|j> where i and j are determinants
  END_DOC
  integer, intent(in)            :: Nint
  integer(bit_kind), intent(in)  :: key_i(Nint,2), key_j(Nint,2)
  double precision, intent(out)  :: hij_core
  double precision, intent(out)  :: hij_hartree
  double precision, intent(out)  :: hij_erf_mu_of_r
  double precision, intent(out)  :: hij_total
  
  integer                        :: exc(0:2,2,2)
  integer                        :: degree
  double precision               :: get_mo_bielec_integral_erf_mu_of_r
  integer                        :: m,n,p,q
  integer                        :: i,j,k
  integer                        :: occ(Nint*bit_kind_size,2)
  double precision               :: diag_H_mat_elem_erf_mu_of_r, phase,phase_2
  integer                        :: n_occ_ab(2)
  PROVIDE mo_bielec_integrals_erf_mu_of_r_in_map mo_integrals_erf_mu_of_r_map big_array_exchange_integrals_erf_mu_of_r
  
  ASSERT (Nint > 0)
  ASSERT (Nint == N_int)
  ASSERT (sum(popcnt(key_i(:,1))) == elec_alpha_num)
  ASSERT (sum(popcnt(key_i(:,2))) == elec_beta_num)
  ASSERT (sum(popcnt(key_j(:,1))) == elec_alpha_num)
  ASSERT (sum(popcnt(key_j(:,2))) == elec_beta_num)
  
  hij_core = 0.d0
  hij_hartree = 0.d0
  hij_erf_mu_of_r = 0.d0
 
  !DIR$ FORCEINLINE
  call get_excitation_degree(key_i,key_j,degree,Nint)
  integer :: spin
  select case (degree)
    case (2)
      call get_double_excitation(key_i,key_j,exc,phase,Nint)
      if (exc(0,1,1) == 1) then
        ! Mono alpha, mono beta
        if(exc(1,1,1) == exc(1,2,2) )then
         hij_erf_mu_of_r = phase * big_array_exchange_integrals(exc(1,1,1),exc(1,1,2),exc(1,2,1))
        else if (exc(1,2,1) ==exc(1,1,2))then
         hij_erf_mu_of_r = phase * big_array_exchange_integrals(exc(1,2,1),exc(1,1,1),exc(1,2,2))
        else
         hij_erf_mu_of_r = phase*get_mo_bielec_integral_erf_mu_of_r(                          &
             exc(1,1,1),                                              &
             exc(1,1,2),                                              &
             exc(1,2,1),                                              &
             exc(1,2,2) ,mo_integrals_erf_mu_of_r_map)
        endif
      else if (exc(0,1,1) == 2) then
        ! Double alpha
        hij_erf_mu_of_r = phase*(get_mo_bielec_integral_erf_mu_of_r(                         &
            exc(1,1,1),                                              &
            exc(2,1,1),                                              &
            exc(1,2,1),                                              &
            exc(2,2,1) ,mo_integrals_erf_mu_of_r_map) -                          &
            get_mo_bielec_integral_erf_mu_of_r(                                  &
            exc(1,1,1),                                              &
            exc(2,1,1),                                              &
            exc(2,2,1),                                              &
            exc(1,2,1) ,mo_integrals_erf_mu_of_r_map) )
      else if (exc(0,1,2) == 2) then
        ! Double beta
        hij_erf_mu_of_r = phase*(get_mo_bielec_integral_erf_mu_of_r(                         &
            exc(1,1,2),                                              &
            exc(2,1,2),                                              &
            exc(1,2,2),                                              &
            exc(2,2,2) ,mo_integrals_erf_mu_of_r_map) -                          &
            get_mo_bielec_integral_erf_mu_of_r(                                  &
            exc(1,1,2),                                              &
            exc(2,1,2),                                              &
            exc(2,2,2),                                              &
            exc(1,2,2) ,mo_integrals_erf_mu_of_r_map) )
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
        do i = 1, n_occ_ab(1)
         hij_erf_mu_of_r += -big_array_exchange_integrals_erf_mu_of_r(occ(i,1),m,p) + big_array_coulomb_integrals_erf_mu_of_r(occ(i,1),m,p)
        enddo
        do i = 1, n_occ_ab(2)
         hij_erf_mu_of_r += big_array_coulomb_integrals_erf_mu_of_r(occ(i,2),m,p)
        enddo
      else
        ! Mono beta
        m = exc(1,1,2)
        p = exc(1,2,2)
        spin = 2
        do i = 1, n_occ_ab(2)
         hij_erf_mu_of_r += -big_array_exchange_integrals_erf_mu_of_r(occ(i,2),m,p) + big_array_coulomb_integrals_erf_mu_of_r(occ(i,2),m,p)
        enddo
        do i = 1, n_occ_ab(1)
         hij_erf_mu_of_r += big_array_coulomb_integrals_erf_mu_of_r(occ(i,1),m,p)
        enddo
      endif
      hij_core = mo_nucl_elec_integral(m,p) + mo_kinetic_integral(m,p) 
      hij_hartree =  short_range_Hartree_operator(m,p)
      hij_total = (hij_erf_mu_of_r + hij_core + hij_hartree) * phase
    case (0)
      call diag_H_mat_elem_erf_mu_of_r_component(key_i,hij_core,hij_hartree,hij_erf_mu_of_r,hij_total,Nint)
  end select
end

subroutine diag_H_mat_elem_erf_mu_of_r_component(key_i,hij_core,hij_hartree,hij_erf_mu_of_r,hij_total,Nint)
 implicit none
 integer(bit_kind), intent(in) :: key_i(N_int,2)
 integer, intent(in)  :: Nint
 double precision, intent(out)  :: hij_core
 double precision, intent(out)  :: hij_hartree
 double precision, intent(out)  :: hij_erf_mu_of_r
 double precision, intent(out)  :: hij_total
 integer :: i,j
 integer                        :: occ(Nint*bit_kind_size,2)
 integer                        :: n_occ_ab(2)

 call bitstring_to_list_ab(key_i, occ, n_occ_ab, Nint)
 hij_core = 0.d0
 hij_hartree = 0.d0
 hij_erf_mu_of_r = 0.d0
 ! alpha - alpha
 do i = 1, n_occ_ab(1)
  hij_core += mo_nucl_elec_integral(occ(i,1),mo_nucl_elec_integral(i,1)) + mo_kinetic_integral(occ(i,1),mo_nucl_elec_integral(i,1))  
  hij_hartree += short_range_Hartree_operator(occ(i,1),occ(i,1))
  do j = i+1, n_occ_ab(1)
   hij_erf_mu_of_r += mo_bielec_integral_erf_mu_of_r_jj_anti(occ(i,1),occ(j,1))
  enddo
 enddo

 ! beta - beta 
 do i = 1, n_occ_ab(2)
  hij_core += mo_nucl_elec_integral(occ(i,2),mo_nucl_elec_integral(i,2)) + mo_kinetic_integral(occ(i,2),mo_nucl_elec_integral(i,2))   
  hij_hartree += short_range_Hartree_operator(occ(i,2),occ(i,2))
  do j = i+1, n_occ_ab(2)
   hij_erf_mu_of_r += mo_bielec_integral_erf_mu_of_r_jj_anti(occ(i,2),occ(j,2))
  enddo
 enddo

 ! alpha - beta 
 do i = 1, n_occ_ab(1)
  do j = 1, n_occ_ab(2)
   hij_erf_mu_of_r += mo_bielec_integral_erf_mu_of_r_jj(occ(i,1),occ(j,2))
  enddo
 enddo
 hij_total = hij_erf_mu_of_r + hij_hartree + hij_core

end



subroutine i_H_j_mono_spin_erf_mu_of_r(key_i,key_j,Nint,spin,hij)
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

  PROVIDE big_array_exchange_integrals_erf_mu_of_r mo_bielec_integrals_erf_mu_of_r_in_map

  call i_H_j_erf_mu_of_r(key_i,key_j,Nint,hij)
! call get_mono_excitation_spin(key_i(1,spin),key_j(1,spin),exc,phase,Nint)
! call get_mono_excitation_from_fock(key_i,key_j,exc(1,1),exc(1,2),spin,phase,hij)
end




subroutine i_H_j_double_spin_erf_mu_of_r(key_i,key_j,Nint,hij)
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
  double precision, external     :: get_mo_bielec_integral_erf_mu_of_r

  PROVIDE big_array_exchange_integrals_erf_mu_of_r mo_bielec_integrals_erf_mu_of_r_in_map

  call get_double_excitation_spin(key_i,key_j,exc,phase,Nint)
  hij = phase*(get_mo_bielec_integral_erf_mu_of_r(                               &
      exc(1,1),                                                      &
      exc(2,1),                                                      &
      exc(1,2),                                                      &
      exc(2,2), mo_integrals_erf_mu_of_r_map) -                                  &
      get_mo_bielec_integral_erf_mu_of_r(                                        &
      exc(1,1),                                                      &
      exc(2,1),                                                      &
      exc(2,2),                                                      &
      exc(1,2), mo_integrals_erf_mu_of_r_map) )
end

subroutine i_H_j_double_alpha_beta_erf_mu_of_r(key_i,key_j,Nint,hij)
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
  double precision, external     :: get_mo_bielec_integral_erf_mu_of_r

  PROVIDE big_array_exchange_integrals_erf_mu_of_r mo_bielec_integrals_erf_mu_of_r_in_map

  call get_mono_excitation_spin(key_i(1,1),key_j(1,1),exc(0,1,1),phase,Nint)
  call get_mono_excitation_spin(key_i(1,2),key_j(1,2),exc(0,1,2),phase2,Nint)
  phase = phase*phase2
  if (exc(1,1,1) == exc(1,2,2)) then
    hij = phase * big_array_exchange_integrals_erf_mu_of_r(exc(1,1,1),exc(1,1,2),exc(1,2,1))
  else if (exc(1,2,1) == exc(1,1,2)) then
    hij = phase * big_array_exchange_integrals_erf_mu_of_r(exc(1,2,1),exc(1,1,1),exc(1,2,2))
  else
    hij = phase*get_mo_bielec_integral_erf_mu_of_r(                              &
        exc(1,1,1),                                                  &
        exc(1,1,2),                                                  &
        exc(1,2,1),                                                  &
        exc(1,2,2) ,mo_integrals_erf_mu_of_r_map)
  endif
end


