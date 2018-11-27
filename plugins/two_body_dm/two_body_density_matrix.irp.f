 BEGIN_PROVIDER [double precision, two_bod_alpha_beta_mo, (mo_tot_num,mo_tot_num,mo_tot_num,mo_tot_num,N_states)]
 implicit none
 BEGIN_DOC
 !  two_bod_alpha_beta(i,j,k,l) = <Psi| a^{dagger}_{j,alpha} a^{dagger}_{l,beta} a_{k,beta} a_{i,alpha} | Psi>
 !                     1 1 2 2  = chemist notations 
 !  note that no 1/2 factor is introduced in order to take into acccount for the spin symmetry
 !  
 END_DOC
 integer :: dim1,dim2,dim3,dim4
 double precision :: cpu_0,cpu_1
 dim1 = mo_tot_num
 dim2 = mo_tot_num
 dim3 = mo_tot_num
 dim4 = mo_tot_num
 two_bod_alpha_beta_mo = 0.d0
 print*,'providing two_bod_alpha_beta ...'
 call cpu_time(cpu_0)
 call two_body_dm_nstates_openmp(two_bod_alpha_beta_mo,dim1,dim2,dim3,dim4,psi_coef,size(psi_coef,2),size(psi_coef,1))
 call cpu_time(cpu_1)
 print*,'two_bod_alpha_beta provided in',dabs(cpu_1-cpu_0)

 END_PROVIDER 


 BEGIN_PROVIDER [double precision, two_bod_alpha_beta_mo_physician, (mo_tot_num,mo_tot_num,mo_tot_num,mo_tot_num,N_states)]
 implicit none
 BEGIN_DOC
 !  two_bod_alpha_beta_mo_physician,(i,j,k,l) = <Psi| a^{dagger}_{k,alpha} a^{dagger}_{l,beta} a_{j,beta} a_{i,alpha} | Psi>
 !                                   1 2 1 2  = physicist notations 
 !  note that no 1/2 factor is introduced in order to take into acccount for the spin symmetry
 !  
 END_DOC
 integer :: i,j,k,l,istate
 double precision :: cpu_0,cpu_1
 two_bod_alpha_beta_mo_physician = 0.d0
 print*,'providing two_bod_alpha_beta_mo_physician ...'
 call cpu_time(cpu_0)
 do istate = 1, N_states 
  do i = 1, mo_tot_num
   do j = 1, mo_tot_num
    do k = 1, mo_tot_num
     do l = 1, mo_tot_num
      !                               1 2 1 2                                 1 1 2 2 
      two_bod_alpha_beta_mo_physician(l,k,i,j,istate) = two_bod_alpha_beta_mo(i,l,j,k,istate)
     enddo
    enddo
   enddo
  enddo
 enddo
 print*,'two_bod_alpha_beta_mo_physician provided in',dabs(cpu_1-cpu_0)

 END_PROVIDER 


 subroutine two_body_dm_nstates_openmp(big_array,dim1,dim2,dim3,dim4,u_0,N_st,sze)
  use bitmasks
  implicit none
  BEGIN_DOC
  ! Computes v_0 = H|u_0> and s_0 = S^2 |u_0>
  !
  ! Assumes that the determinants are in psi_det
  !
  ! istart, iend, ishift, istep are used in ZMQ parallelization.
  END_DOC
  integer, intent(in)            :: N_st,sze
  integer, intent(in) :: dim1,dim2,dim3,dim4
  double precision, intent(inout) :: big_array(dim1,dim2,dim3,dim4,N_states)
  double precision, intent(inout) :: u_0(sze,N_st)
  integer :: k
  double precision, allocatable  :: u_t(:,:)
  !DIR$ ATTRIBUTES ALIGN : $IRP_ALIGN :: u_t
  allocate(u_t(N_st,N_det))
  do k=1,N_st
    call dset_order(u_0(1,k),psi_bilinear_matrix_order,N_det)
  enddo
  call dtranspose(                                                   &
      u_0,                                                           &
      size(u_0, 1),                                                  &
      u_t,                                                           &
      size(u_t, 1),                                                  &
      N_det, N_st)

  call two_body_dm_nstates_openmp_work(big_array,dim1,dim2,dim3,dim4,u_t,N_st,sze,1,N_det,0,1)
  deallocate(u_t)

  do k=1,N_st
    call dset_order(u_0(1,k),psi_bilinear_matrix_order_reverse,N_det)
  enddo

 end


 subroutine two_body_dm_nstates_openmp_work(big_array,dim1,dim2,dim3,dim4,u_t,N_st,sze,istart,iend,ishift,istep)
  use bitmasks
  implicit none
  BEGIN_DOC
  ! Computes v_0 = H|u_0> and s_0 = S^2 |u_0>
  !
  ! Default should be 1,N_det,0,1
  END_DOC
  integer, intent(in)            :: N_st,sze,istart,iend,ishift,istep
  integer, intent(in) :: dim1,dim2,dim3,dim4
  double precision, intent(inout) :: big_array(dim1,dim2,dim3,dim4,N_states)
  double precision, intent(in)   :: u_t(N_st,N_det)

  
  PROVIDE N_int 

  select case (N_int)
    case (1)
      call two_body_dm_nstates_openmp_work_1(big_array,dim1,dim2,dim3,dim4,u_t,N_st,sze,istart,iend,ishift,istep)
    case (2)
      call two_body_dm_nstates_openmp_work_2(big_array,dim1,dim2,dim3,dim4,u_t,N_st,sze,istart,iend,ishift,istep)
    case (3)
      call two_body_dm_nstates_openmp_work_3(big_array,dim1,dim2,dim3,dim4,u_t,N_st,sze,istart,iend,ishift,istep)
    case (4)
      call two_body_dm_nstates_openmp_work_4(big_array,dim1,dim2,dim3,dim4,u_t,N_st,sze,istart,iend,ishift,istep)
    case default
      call two_body_dm_nstates_openmp_work_N_int(big_array,dim1,dim2,dim3,dim4,u_t,N_st,sze,istart,iend,ishift,istep)
  end select
 end
 BEGIN_TEMPLATE

 subroutine two_body_dm_nstates_openmp_work_$N_int(big_array,dim1,dim2,dim3,dim4,u_t,N_st,sze,istart,iend,ishift,istep)
  use bitmasks
  implicit none
  integer, intent(in)            :: N_st,sze,istart,iend,ishift,istep
  integer, intent(in) :: dim1,dim2,dim3,dim4
  double precision, intent(inout) :: big_array(dim1,dim2,dim3,dim4,N_states)
  double precision, intent(in)   :: u_t(N_st,N_det)

  double precision               :: hij, sij
  integer                        :: i,j,k,l
  integer                        :: k_a, k_b, l_a, l_b, m_a, m_b
  integer                        :: istate
  integer                        :: krow, kcol, krow_b, kcol_b
  integer                        :: lrow, lcol
  integer                        :: mrow, mcol
  integer(bit_kind)              :: spindet($N_int)
  integer(bit_kind)              :: tmp_det($N_int,2)
  integer(bit_kind)              :: tmp_det2($N_int,2)
  integer(bit_kind)              :: tmp_det3($N_int,2)
  integer(bit_kind), allocatable :: buffer(:,:)
  integer                        :: n_doubles
  integer, allocatable           :: doubles(:)
  integer, allocatable           :: singles_a(:)
  integer, allocatable           :: singles_b(:)
  integer, allocatable           :: idx(:), idx0(:)
  integer                        :: maxab, n_singles_a, n_singles_b, kcol_prev, nmax
  integer*8                      :: k8

  maxab = max(N_det_alpha_unique, N_det_beta_unique)+1
  allocate(idx0(maxab))
  
  do i=1,maxab
    idx0(i) = i
  enddo

  ! Prepare the array of all alpha single excitations
  ! -------------------------------------------------

  PROVIDE N_int nthreads_davidson
  
  ! Alpha/Beta double excitations
  ! =============================
    
  allocate( buffer($N_int,maxab),                                     &
      singles_a(maxab),                                              &
      singles_b(maxab),                                              &
      doubles(maxab),                                                &
      idx(maxab))

  kcol_prev=-1

  ASSERT (iend <= N_det)
  ASSERT (istart > 0)
  ASSERT (istep  > 0)

  do k_a=istart+ishift,iend,istep

    krow = psi_bilinear_matrix_rows(k_a)
    ASSERT (krow <= N_det_alpha_unique)

    kcol = psi_bilinear_matrix_columns(k_a)
    ASSERT (kcol <= N_det_beta_unique)

    tmp_det(1:$N_int,1) = psi_det_alpha_unique(1:$N_int, krow)
    tmp_det(1:$N_int,2) = psi_det_beta_unique (1:$N_int, kcol)

    if (kcol /= kcol_prev) then
      call get_all_spin_singles_$N_int(                              &
          psi_det_beta_unique, idx0,                                 &
          tmp_det(1,2), N_det_beta_unique,                           &
          singles_b, n_singles_b)
    endif
    kcol_prev = kcol

    ! Loop over singly excited beta columns
    ! -------------------------------------

    do i=1,n_singles_b
      lcol = singles_b(i)

      tmp_det2(1:$N_int,2) = psi_det_beta_unique(1:$N_int, lcol)

      l_a = psi_bilinear_matrix_columns_loc(lcol)
      ASSERT (l_a <= N_det)

      do j=1,psi_bilinear_matrix_columns_loc(lcol+1) - l_a
        lrow = psi_bilinear_matrix_rows(l_a)
        ASSERT (lrow <= N_det_alpha_unique)

        buffer(1:$N_int,j) = psi_det_alpha_unique(1:$N_int, lrow)

        ASSERT (l_a <= N_det)
        idx(j) = l_a
        l_a = l_a+1
      enddo
      j = j-1

      call get_all_spin_singles_$N_int(                              &
          buffer, idx, tmp_det(1,1), j,                              &
          singles_a, n_singles_a )

      ! Loop over alpha singles
      ! -----------------------

      do k = 1,n_singles_a
        l_a = singles_a(k)
        ASSERT (l_a <= N_det)

        lrow = psi_bilinear_matrix_rows(l_a)
        ASSERT (lrow <= N_det_alpha_unique)

        tmp_det2(1:$N_int,1) = psi_det_alpha_unique(1:$N_int, lrow)
        !!!!!!!!!!!!!!!!!! ALPHA BETA 
        do l= 1, N_states
         c_1(l) = u_t(l,l_a)
         c_2(l) = u_t(l,k_a)
        enddo
        call off_diagonal_double_to_two_body_ab_dm(tmp_det,tmp_det2,c_1,c_2,big_array,dim1,dim2,dim3,dim4)
      enddo

    enddo

  enddo


  do k_a=istart+ishift,iend,istep


    ! Single and double alpha excitations
    ! ===================================
    
    
    ! Initial determinant is at k_a in alpha-major representation
    ! -----------------------------------------------------------------------
    
    krow = psi_bilinear_matrix_rows(k_a)
    ASSERT (krow <= N_det_alpha_unique)

    kcol = psi_bilinear_matrix_columns(k_a)
    ASSERT (kcol <= N_det_beta_unique)
    
    tmp_det(1:$N_int,1) = psi_det_alpha_unique(1:$N_int, krow)
    tmp_det(1:$N_int,2) = psi_det_beta_unique (1:$N_int, kcol)
    
    ! Initial determinant is at k_b in beta-major representation
    ! ----------------------------------------------------------------------
    
    k_b = psi_bilinear_matrix_order_transp_reverse(k_a)

    spindet(1:$N_int) = tmp_det(1:$N_int,1)
    
    ! Loop inside the beta column to gather all the connected alphas
    lcol = psi_bilinear_matrix_columns(k_a)
    l_a = psi_bilinear_matrix_columns_loc(lcol)
    do i=1,N_det_alpha_unique
      if (l_a > N_det) exit
      lcol = psi_bilinear_matrix_columns(l_a)
      if (lcol /= kcol) exit
      lrow = psi_bilinear_matrix_rows(l_a)
      ASSERT (lrow <= N_det_alpha_unique)

      buffer(1:$N_int,i) = psi_det_alpha_unique(1:$N_int, lrow)
      idx(i) = l_a
      l_a = l_a+1
    enddo
    i = i-1
    
    call get_all_spin_singles_and_doubles_$N_int(                    &
        buffer, idx, spindet, i,                                     &
        singles_a, doubles, n_singles_a, n_doubles )

    ! Compute Hij for all alpha singles
    ! ----------------------------------

    tmp_det2(1:$N_int,2) = psi_det_beta_unique (1:$N_int, kcol)
    do i=1,n_singles_a
      l_a = singles_a(i)
      ASSERT (l_a <= N_det)

      lrow = psi_bilinear_matrix_rows(l_a)
      ASSERT (lrow <= N_det_alpha_unique)

      tmp_det2(1:$N_int,1) = psi_det_alpha_unique(1:$N_int, lrow)
      !!!! MONO SPIN 
      do l= 1, N_states
       c_1(l) = u_t(l,l_a)
       c_2(l) = u_t(l,k_a)
      enddo
      call off_diagonal_single_to_two_body_ab_dm(tmp_det, tmp_det2,c_1,c_2,big_array,dim1,dim2,dim3,dim4)

    enddo

    
   !! Compute Hij for all alpha doubles
   !! ----------------------------------
   !
   !do i=1,n_doubles
   !  l_a = doubles(i)
   !  ASSERT (l_a <= N_det)

   !  lrow = psi_bilinear_matrix_rows(l_a)
   !  ASSERT (lrow <= N_det_alpha_unique)

   !  call i_H_j_double_spin_erf( tmp_det(1,1), psi_det_alpha_unique(1, lrow), $N_int, hij)
   !  do l=1,N_st
   !    v_t(l,k_a) = v_t(l,k_a) + hij * u_t(l,l_a)
   !    ! same spin => sij = 0
   !  enddo
   !enddo
    


    ! Single and double beta excitations
    ! ==================================

    
    ! Initial determinant is at k_a in alpha-major representation
    ! -----------------------------------------------------------------------
    
    krow = psi_bilinear_matrix_rows(k_a)
    kcol = psi_bilinear_matrix_columns(k_a)
    
    tmp_det(1:$N_int,1) = psi_det_alpha_unique(1:$N_int, krow)
    tmp_det(1:$N_int,2) = psi_det_beta_unique (1:$N_int, kcol)
    
    spindet(1:$N_int) = tmp_det(1:$N_int,2)
    
    ! Initial determinant is at k_b in beta-major representation
    ! -----------------------------------------------------------------------

    k_b = psi_bilinear_matrix_order_transp_reverse(k_a) 
    
    ! Loop inside the alpha row to gather all the connected betas
    lrow = psi_bilinear_matrix_transp_rows(k_b)
    l_b = psi_bilinear_matrix_transp_rows_loc(lrow)
    do i=1,N_det_beta_unique
      if (l_b > N_det) exit
      lrow = psi_bilinear_matrix_transp_rows(l_b)
      if (lrow /= krow) exit
      lcol = psi_bilinear_matrix_transp_columns(l_b)
      ASSERT (lcol <= N_det_beta_unique)

      buffer(1:$N_int,i) = psi_det_beta_unique(1:$N_int, lcol)
      idx(i) = l_b
      l_b = l_b+1
    enddo
    i = i-1
  
    call get_all_spin_singles_and_doubles_$N_int(                    &
        buffer, idx, spindet, i,                                     &
        singles_b, doubles, n_singles_b, n_doubles )
    
    ! Compute Hij for all beta singles
    ! ----------------------------------
    
    tmp_det2(1:$N_int,1) = psi_det_alpha_unique(1:$N_int, krow)
    do i=1,n_singles_b
      l_b = singles_b(i)
      ASSERT (l_b <= N_det)

      lcol = psi_bilinear_matrix_transp_columns(l_b)
      ASSERT (lcol <= N_det_beta_unique)

      tmp_det2(1:$N_int,2) = psi_det_beta_unique (1:$N_int, lcol)
      l_a = psi_bilinear_matrix_transp_order(l_b)
      do l= 1, N_states
       c_1(l) = u_t(l,l_a)
       c_2(l) = u_t(l,k_a)
      enddo
      call off_diagonal_single_to_two_body_ab_dm(tmp_det, tmp_det2,c_1,c_2,big_array,dim1,dim2,dim3,dim4)
      ASSERT (l_a <= N_det)
    enddo
   !
   !! Compute Hij for all beta doubles
   !! ----------------------------------
   !
   !do i=1,n_doubles
   !  l_b = doubles(i)
   !  ASSERT (l_b <= N_det)

   !  lcol = psi_bilinear_matrix_transp_columns(l_b)
   !  ASSERT (lcol <= N_det_beta_unique)

   !  call i_H_j_double_spin_erf( tmp_det(1,2), psi_det_beta_unique(1, lcol), $N_int, hij)
   !  l_a = psi_bilinear_matrix_transp_order(l_b)
   !  ASSERT (l_a <= N_det)

   !  do l=1,N_st
   !    v_t(l,k_a) = v_t(l,k_a) + hij * u_t(l,l_a)
   !    ! same spin => sij = 0 
   !  enddo
   !enddo


    ! Diagonal contribution
    ! =====================

    
    ! Initial determinant is at k_a in alpha-major representation
    ! -----------------------------------------------------------------------
    
    krow = psi_bilinear_matrix_rows(k_a)
    ASSERT (krow <= N_det_alpha_unique)

    kcol = psi_bilinear_matrix_columns(k_a)
    ASSERT (kcol <= N_det_beta_unique)
    
    tmp_det(1:$N_int,1) = psi_det_alpha_unique(1:$N_int, krow)
    tmp_det(1:$N_int,2) = psi_det_beta_unique (1:$N_int, kcol)
    
    double precision, external :: diag_H_mat_elem_erf, diag_S_mat_elem
    double precision :: c_1(N_states),c_2(N_states)
    do l = 1, N_states
     c_1(l) = u_t(l,k_a)
    enddo
  
    call diagonal_contrib_to_two_body_ab_dm(tmp_det,c_1,big_array,dim1,dim2,dim3,dim4)

  end do
  deallocate(buffer, singles_a, singles_b, doubles, idx)

 end

 SUBST [ N_int ]

 1;;
 2;;
 3;;
 4;;
 N_int;;
 
 END_TEMPLATE

 subroutine diagonal_contrib_to_two_body_ab_dm(det_1,c_1,big_array,dim1,dim2,dim3,dim4)
 use bitmasks
 implicit none
 integer, intent(in) :: dim1,dim2,dim3,dim4
 double precision, intent(inout) :: big_array(dim1,dim2,dim3,dim4,N_states)
 integer(bit_kind), intent(in)  :: det_1(N_int,2)
 double precision, intent(in)   :: c_1(N_states)
 integer                        :: occ(N_int*bit_kind_size,2)
 integer                        :: n_occ_ab(2)
 integer :: i,j,h1,h2,istate
 double precision               :: c_1_bis
 call bitstring_to_list_ab(det_1, occ, n_occ_ab, N_int)
 do istate = 1, N_states
  c_1_bis = c_1(istate) * c_1(istate)
  do i = 1, n_occ_ab(1)
   h1 = occ(i,1)
   do j = 1, n_occ_ab(2)
    h2 = occ(j,2)
    big_array(h1,h1,h2,h2,istate) += c_1_bis 
   enddo 
  enddo
 enddo
 end

 subroutine diagonal_contrib_to_all_two_body_dm(det_1,c_1,big_array_ab,big_array_aa,big_array_bb,dim1,dim2,dim3,dim4)
 use bitmasks
 implicit none
 integer, intent(in) :: dim1,dim2,dim3,dim4
 double precision, intent(inout) :: big_array_ab(dim1,dim2,dim3,dim4,N_states)
 double precision, intent(inout) :: big_array_aa(dim1,dim2,dim3,dim4,N_states)
 double precision, intent(inout) :: big_array_bb(dim1,dim2,dim3,dim4,N_states)
 integer(bit_kind), intent(in)  :: det_1(N_int,2)
 double precision, intent(in)   :: c_1(N_states)
 integer                        :: occ(N_int*bit_kind_size,2)
 integer                        :: n_occ_ab(2)
 integer :: i,j,h1,h2,istate
 double precision               :: c_1_bis
 BEGIN_DOC
! no factor 1/2 have to be taken into account as the permutations are already taken into account
 END_DOC
 call bitstring_to_list_ab(det_1, occ, n_occ_ab, N_int)
 do istate = 1, N_states
  c_1_bis = c_1(istate) * c_1(istate)
  do i = 1, n_occ_ab(1)
   h1 = occ(i,1)
   do j = 1, n_occ_ab(2)
    h2 = occ(j,2)
    big_array_ab(h1,h1,h2,h2,istate) += c_1_bis 
   enddo 
   do j = 1, n_occ_ab(1)
    h2 = occ(j,1)
    big_array_aa(h1,h2,h1,h2,istate) -= c_1_bis 
    big_array_aa(h1,h1,h2,h2,istate) += c_1_bis 
   enddo
  enddo
  do i = 1, n_occ_ab(2)
   h1 = occ(i,2)
   do j = 1, n_occ_ab(2)
    h2 = occ(j,2)
    big_array_bb(h1,h1,h2,h2,istate) += c_1_bis 
    big_array_bb(h1,h2,h1,h2,istate) -= c_1_bis 
   enddo
  enddo
 enddo
 end


 subroutine off_diagonal_double_to_two_body_ab_dm(det_1,det_2,c_1,c_2,big_array,dim1,dim2,dim3,dim4)
 use bitmasks
 implicit none
 integer, intent(in) :: dim1,dim2,dim3,dim4
 double precision, intent(inout) :: big_array(dim1,dim2,dim3,dim4,N_states)
 integer(bit_kind), intent(in)  :: det_1(N_int,2),det_2(N_int,2)
 double precision, intent(in)   :: c_1(N_states),c_2(N_states)
 integer :: i,j,h1,h2,p1,p2,istate
 integer                        :: exc(0:2,2,2)
 double precision               :: phase
 call get_double_excitation(det_1,det_2,exc,phase,N_int)
 h1 = exc(1,1,1) 
 h2 = exc(1,1,2) 
 p1 = exc(1,2,1)
 p2 = exc(1,2,2)
 do istate = 1, N_states
  big_array(h1,p1,h2,p2,istate) += c_1(istate) * phase * c_2(istate)
! big_array(p1,h1,p2,h2,istate) += c_1(istate) * phase * c_2(istate)
 enddo
 end

 subroutine off_diagonal_single_to_two_body_ab_dm(det_1,det_2,c_1,c_2,big_array,dim1,dim2,dim3,dim4)
 use bitmasks
 implicit none
 integer, intent(in) :: dim1,dim2,dim3,dim4
 double precision, intent(inout) :: big_array(dim1,dim2,dim3,dim4,N_states)
 integer(bit_kind), intent(in)  :: det_1(N_int,2),det_2(N_int,2)
 double precision, intent(in)   :: c_1(N_states),c_2(N_states)
 integer                        :: occ(N_int*bit_kind_size,2)
 integer                        :: n_occ_ab(2)
 integer :: i,j,h1,h2,istate,p1
 integer                        :: exc(0:2,2,2)
 double precision               :: phase
 call bitstring_to_list_ab(det_1, occ, n_occ_ab, N_int)
 call get_mono_excitation(det_1,det_2,exc,phase,N_int)
 if (exc(0,1,1) == 1) then
  ! Mono alpha
  h1 = exc(1,1,1)
  p1 = exc(1,2,1)
  do istate = 1, N_states
   do i = 1, n_occ_ab(2)
    h2 = occ(i,2)
    big_array(h1,p1,h2,h2,istate) += 1.d0 * c_1(istate) * c_2(istate) * phase
    enddo 
  enddo
 else 
  ! Mono beta
  h1 = exc(1,1,2)
  p1 = exc(1,2,2)
  do istate = 1, N_states
   do i = 1, n_occ_ab(1)
    h2 = occ(i,1)
    big_array(h2,h2,h1,p1,istate) += 1.d0 * c_1(istate) * c_2(istate) * phase
    enddo 
  enddo
 endif
 end


 BEGIN_PROVIDER [double precision, two_bod_alpha_beta_ao, (ao_num,ao_num,ao_num,ao_num,N_states)]
 implicit none
 BEGIN_DOC
 !  two_bod_alpha_beta(i,j,k,l) = <Psi| a^{dagger}_{j,alpha} a^{dagger}_{l,beta} a_{k,beta} a_{i,alpha} | Psi>
 !  note that no 1/2 factor is introduced in order to take into acccount for the spin symmetry
 END_DOC
  two_bod_alpha_beta_ao = 0.d0
 END_PROVIDER 




