subroutine u_0_H_u_0_monoelec_dft(e_0,u_0,n,keys_tmp,Nint,N_st,sze)
  use bitmasks
  implicit none
  BEGIN_DOC
  ! Computes e_0 = <u_0|H|u_0>/<u_0|u_0>
  !
  ! n : number of determinants
  !
  END_DOC
  integer, intent(in)            :: n,Nint, N_st, sze
  double precision, intent(out)  :: e_0(N_st)
  double precision, intent(inout):: u_0(sze,N_st)
  integer(bit_kind),intent(in)   :: keys_tmp(Nint,2,n)
  
  double precision, allocatable  :: v_0(:,:), s_0(:,:)
  double precision               :: u_dot_u,u_dot_v,diag_H_mat_elem
  integer :: i,j
  allocate (v_0(sze,N_st),s_0(sze,N_st))
  call H_S2_u_0_monoelec_dft_nstates_openmp(v_0,s_0,u_0,N_st,sze)
  do i=1,N_st
    e_0(i) = u_dot_v(v_0(1,i),u_0(1,i),n)/u_dot_u(u_0(1,i),n)
  enddo
  deallocate (s_0, v_0)
end

BEGIN_PROVIDER [ double precision, psi_energy_monoelec_dft, (N_states) ]
  implicit none
  BEGIN_DOC
! Energy of the current wave function
  END_DOC
! call u_0_H_u_0_monoelec_dft(psi_energy_monoelec_dft,psi_coef,N_det,psi_det,N_int,N_states,psi_det_size)
  integer :: i,j
  double precision :: accu, hij
  do i = 1, N_det
   do j = 1, N_det
    call i_H_j_mono_monoelec_dft(psi_det(1,1,j),psi_det(1,1,i),N_int,hij)
   enddo
  enddo
END_PROVIDER



subroutine H_S2_u_0_monoelec_dft_nstates_openmp(v_0,s_0,u_0,N_st,sze)
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
  double precision, intent(inout)  :: v_0(sze,N_st), s_0(sze,N_st), u_0(sze,N_st)
  integer :: k
  double precision, allocatable  :: u_t(:,:)
  !DIR$ ATTRIBUTES ALIGN : $IRP_ALIGN :: u_t
  allocate(u_t(N_st,N_det))
  do k=1,N_st
    call dset_order(u_0(1,k),psi_bilinear_matrix_order,N_det)
  enddo
  v_0 = 0.d0
  s_0 = 0.d0
  call dtranspose(                                                   &
      u_0,                                                           &
      size(u_0, 1),                                                  &
      u_t,                                                           &
      size(u_t, 1),                                                  &
      N_det, N_st)

  call H_S2_u_0_monoelec_dft_nstates_openmp_work(v_0,s_0,u_t,N_st,sze,1,N_det,0,1)
  deallocate(u_t)

  do k=1,N_st
    call dset_order(v_0(1,k),psi_bilinear_matrix_order_reverse,N_det)
    call dset_order(s_0(1,k),psi_bilinear_matrix_order_reverse,N_det)
    call dset_order(u_0(1,k),psi_bilinear_matrix_order_reverse,N_det)
  enddo

end


subroutine H_S2_u_0_monoelec_dft_nstates_openmp_work(v_0,s_0,u_t,N_st,sze,istart,iend,ishift,istep)
  use bitmasks
  implicit none
  BEGIN_DOC
  ! Computes v_0 = H|u_0> and s_0 = S^2 |u_0>
  !
  ! Default should be 1,N_det,0,1
  END_DOC
  integer, intent(in)            :: N_st,sze,istart,iend,ishift,istep
  double precision, intent(in)   :: u_t(N_st,N_det)
  double precision, intent(out)  :: v_0(sze,N_st), s_0(sze,N_st) 

  
  PROVIDE bi_elec_ref_bitmask_energy short_range_Hartree N_int

  select case (N_int)
    case (1)
      call H_S2_u_0_monoelec_dft_nstates_openmp_work_1(v_0,s_0,u_t,N_st,sze,istart,iend,ishift,istep)
    case (2)
      call H_S2_u_0_monoelec_dft_nstates_openmp_work_2(v_0,s_0,u_t,N_st,sze,istart,iend,ishift,istep)
    case (3)
      call H_S2_u_0_monoelec_dft_nstates_openmp_work_3(v_0,s_0,u_t,N_st,sze,istart,iend,ishift,istep)
    case (4)
      call H_S2_u_0_monoelec_dft_nstates_openmp_work_4(v_0,s_0,u_t,N_st,sze,istart,iend,ishift,istep)
    case default
      call H_S2_u_0_monoelec_dft_nstates_openmp_work_N_int(v_0,s_0,u_t,N_st,sze,istart,iend,ishift,istep)
  end select
end
BEGIN_TEMPLATE

subroutine H_S2_u_0_monoelec_dft_nstates_openmp_work_$N_int(v_0,s_0,u_t,N_st,sze,istart,iend,ishift,istep)
  use bitmasks
  implicit none
  BEGIN_DOC
  ! Computes v_0 = H|u_0> and s_0 = S^2 |u_0>
  !
  ! Default should be 1,N_det,0,1
  END_DOC
  integer, intent(in)            :: N_st,sze,istart,iend,ishift,istep
  double precision, intent(in)   :: u_t(N_st,N_det)
  double precision, intent(out)  :: v_0(sze,N_st), s_0(sze,N_st) 

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
  double precision, allocatable  :: v_t(:,:), s_t(:,:)
  !DIR$ ATTRIBUTES ALIGN : $IRP_ALIGN :: v_t, s_t

  maxab = max(N_det_alpha_unique, N_det_beta_unique)+1
  allocate(idx0(maxab))
  
  do i=1,maxab
    idx0(i) = i
  enddo

  ! Prepare the array of all alpha single excitations
  ! -------------------------------------------------

  PROVIDE N_int
  !$OMP PARALLEL DEFAULT(NONE)                                       &
      !$OMP   SHARED(psi_bilinear_matrix_rows, N_det,                &
      !$OMP          psi_bilinear_matrix_columns,                    &
      !$OMP          psi_det_alpha_unique, psi_det_beta_unique,      &
      !$OMP          n_det_alpha_unique, n_det_beta_unique, N_int,   &
      !$OMP          psi_bilinear_matrix_transp_rows,                &
      !$OMP          psi_bilinear_matrix_transp_columns,             &
      !$OMP          psi_bilinear_matrix_transp_order, N_st,         &
      !$OMP          psi_bilinear_matrix_order_transp_reverse,       &
      !$OMP          psi_bilinear_matrix_columns_loc,                &
      !$OMP          istart, iend, istep,        &
      !$OMP          ishift, idx0, u_t, maxab, v_0, s_0)             &
      !$OMP   PRIVATE(krow, kcol, tmp_det, spindet, k_a, k_b, i,     &
      !$OMP          lcol, lrow, l_a, l_b, nmax,         &
      !$OMP          buffer, doubles, n_doubles, &
      !$OMP          tmp_det2, hij, sij, idx, l, kcol_prev, v_t,     &
      !$OMP          singles_a, n_singles_a, singles_b,              &
      !$OMP          n_singles_b, s_t, k8)
  
  ! Alpha/Beta double excitations
  ! =============================
    
  allocate( buffer($N_int,maxab),                                     &
      singles_a(maxab),                                              &
      singles_b(maxab),                                              &
      doubles(maxab),                                                &
      idx(maxab),                                                    &
      v_t(N_st,N_det), s_t(N_st,N_det))
  kcol_prev=-1

  v_t = 0.d0
  s_t = 0.d0


! !$OMP DO SCHEDULE(dynamic,64)
! do k_a=istart+ishift,iend,istep

!   krow = psi_bilinear_matrix_rows(k_a)
!   kcol = psi_bilinear_matrix_columns(k_a)
!   
!   tmp_det(1:$N_int,1) = psi_det_alpha_unique(1:$N_int, krow)
!   tmp_det(1:$N_int,2) = psi_det_beta_unique (1:$N_int, kcol)
!   
!   if (kcol /= kcol_prev) then
!     call get_all_spin_singles_$N_int(                              &
!         psi_det_beta_unique(1,kcol+1), idx0(kcol+1),               &
!         tmp_det(1,2), N_det_beta_unique-kcol,                      &
!         singles_b, n_singles_b)
!   endif
!   kcol_prev = kcol

!   ! Loop over singly excited beta columns > current column
!   ! ------------------------------------------------------

!   do i=1,n_singles_b
!     lcol = singles_b(i)

!     tmp_det2(1:$N_int,2) = psi_det_beta_unique(1:$N_int, lcol)

!     l_a = psi_bilinear_matrix_columns_loc(lcol)

!     nmax = psi_bilinear_matrix_columns_loc(lcol+1) - l_a
!     do j=1,nmax
!       lrow = psi_bilinear_matrix_rows(l_a)
!       buffer(1:$N_int,j) = psi_det_alpha_unique(1:$N_int, lrow)
!       idx(j) = l_a
!       l_a = l_a+1
!     enddo
!     j = j-1
!     
!     call get_all_spin_singles_$N_int(                              &
!         buffer, idx, tmp_det(1,1), j,                              &
!         singles_a, n_singles_a )

!   enddo

! enddo
! !$OMP END DO 

  !$OMP DO SCHEDULE(dynamic,64)
  do k_a=istart+ishift,iend,istep


    ! Single and double alpha excitations
    ! ===================================
    
    
    ! Initial determinant is at k_a in alpha-major representation
    ! -----------------------------------------------------------------------
    
    krow = psi_bilinear_matrix_rows(k_a)
    kcol = psi_bilinear_matrix_columns(k_a)
    
    tmp_det(1:$N_int,1) = psi_det_alpha_unique(1:$N_int, krow)
    tmp_det(1:$N_int,2) = psi_det_beta_unique (1:$N_int, kcol)
    
    ! Initial determinant is at k_b in beta-major representation
    ! ----------------------------------------------------------------------
    
    k_b = psi_bilinear_matrix_order_transp_reverse(k_a)

    spindet(1:$N_int) = tmp_det(1:$N_int,1)
    
    ! Loop inside the beta column to gather all the connected alphas
    l_a = k_a+1
    nmax = min(N_det_alpha_unique, N_det - l_a)
    do i=1,nmax
      lcol = psi_bilinear_matrix_columns(l_a)
      if (lcol /= kcol) exit
      lrow = psi_bilinear_matrix_rows(l_a)
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
      lrow = psi_bilinear_matrix_rows(l_a)
      tmp_det2(1:$N_int,1) = psi_det_alpha_unique(1:$N_int, lrow)
      call i_H_j_mono_spin_monoelec_dft( tmp_det, tmp_det2, $N_int, 1, hij)
      do l=1,N_st
        v_t(l,l_a) = v_t(l,l_a) + hij * u_t(l,k_a)
        v_t(l,k_a) = v_t(l,k_a) + hij * u_t(l,l_a)
        ! single => sij = 0 
      enddo
    enddo

    
    ! Single beta excitations
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
    l_b = k_b+1
    nmax = min(N_det_beta_unique, N_det - l_b)
    do i=1,nmax
      lrow = psi_bilinear_matrix_transp_rows(l_b)
      if (lrow /= krow) exit
      lcol = psi_bilinear_matrix_transp_columns(l_b)
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
      lcol = psi_bilinear_matrix_transp_columns(l_b)
      tmp_det2(1:$N_int,2) = psi_det_beta_unique (1:$N_int, lcol)
      call i_H_j_mono_spin_monoelec_dft( tmp_det, tmp_det2, $N_int, 2, hij)
      l_a = psi_bilinear_matrix_transp_order(l_b)
      do l=1,N_st
        v_t(l,l_a) = v_t(l,l_a) + hij * u_t(l,k_a)
        v_t(l,k_a) = v_t(l,k_a) + hij * u_t(l,l_a)
        ! single => sij = 0 
      enddo
    enddo
    
    ! Diagonal contribution
    ! =====================

    
    ! Initial determinant is at k_a in alpha-major representation
    ! -----------------------------------------------------------------------
    
    krow = psi_bilinear_matrix_rows(k_a)
    kcol = psi_bilinear_matrix_columns(k_a)
    
    tmp_det(1:$N_int,1) = psi_det_alpha_unique(1:$N_int, krow)
    tmp_det(1:$N_int,2) = psi_det_beta_unique (1:$N_int, kcol)
    
    double precision, external :: diag_H_mat_elem_monoelec_dft, diag_S_mat_elem
  
    hij = diag_H_mat_elem_monoelec_dft(tmp_det,$N_int) 
    sij = diag_S_mat_elem(tmp_det,$N_int)
    do l=1,N_st
      v_t(l,k_a) = v_t(l,k_a) + hij * u_t(l,k_a)
      s_t(l,k_a) = s_t(l,k_a) + sij * u_t(l,k_a)
    enddo

  end do
  !$OMP END DO NOWAIT
  deallocate(buffer, singles_a, singles_b, doubles, idx)

  !$OMP CRITICAL
  do l=1,N_st
    do i=1, N_det
      v_0(i,l) = v_0(i,l) + v_t(l,i)
      s_0(i,l) = s_0(i,l) + s_t(l,i)
    enddo
  enddo
  !$OMP END CRITICAL
  deallocate(v_t, s_t)

  !$OMP BARRIER
  !$OMP END PARALLEL

end

SUBST [ N_int ]

1;;
2;;
3;;
4;;
N_int;;

END_TEMPLATE


