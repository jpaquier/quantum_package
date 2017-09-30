subroutine davidson_diag_general_bis(u_in,dim_in,energies,sze,N_st,N_st_diag,Nint,iunit)
  use bitmasks
  implicit none
  BEGIN_DOC
  ! Davidson diagonalization.
  !
  ! u_in : guess coefficients on the various states. Overwritten
  !   on exit
  !
  ! dim_in : leftmost dimension of u_in
  !
  ! sze : Number of determinants
  !
  ! N_st : Number of eigenstates
  !
  ! iunit : Unit number for the I/O
  !
  ! Initial guess vectors are not necessarily orthonormal
  END_DOC
  integer, intent(in)            :: dim_in, sze, N_st, N_st_diag, Nint, iunit
  double precision, intent(inout) :: u_in(dim_in,N_st_diag)
  double precision, intent(out)  :: energies(N_st_diag)
  double precision, allocatable  :: H_jj(:)
  
  double precision               :: diag_H_mat_elem, diag_S_mat_elem
  integer                        :: i
  ASSERT (N_st > 0)
  ASSERT (sze > 0)
  ASSERT (Nint > 0)
  ASSERT (Nint == N_int)
  PROVIDE mo_bielec_integrals_in_map
  allocate(H_jj(sze) )
 print*, 'dim_in',dim_in
 print*, 'sze ',sze
  
  do i=1,sze
    H_jj(i) = diagonal_superci_matrix(i)
  enddo

  call davidson_diag_hjj_bis(u_in,H_jj,energies,dim_in,sze,N_st,N_st_diag,Nint,iunit)
  deallocate (H_jj)
end


subroutine davidson_diag_hjj_bis(u_in,H_jj,energies,dim_in,sze,N_st,N_st_diag,Nint,iunit)
  use bitmasks
  implicit none
  BEGIN_DOC
  ! Davidson diagonalization with specific diagonal elements of the H matrix
  !
  ! H_jj : specific diagonal H matrix elements to diagonalize de Davidson
  !
  ! dets_in : bitmasks corresponding to determinants
  !
  ! u_in : guess coefficients on the various states. Overwritten
  !   on exit
  !
  ! dim_in : leftmost dimension of u_in
  !
  ! sze : Number of determinants
  !
  ! N_st : Number of eigenstates
  ! 
  ! N_st_diag : Number of states in which H is diagonalized. Assumed > sze
  !
  ! iunit : Unit for the I/O
  !
  ! Initial guess vectors are not necessarily orthonormal
  END_DOC
  integer, intent(in)            :: dim_in, sze, N_st, N_st_diag, Nint
  double precision,  intent(in)  :: H_jj(sze)
  integer,  intent(in)           :: iunit
  double precision, intent(inout) :: u_in(dim_in,N_st_diag)
  double precision, intent(out)  :: energies(N_st_diag)
  
  integer                        :: iter
  integer                        :: i,j,k,l,m
  logical                        :: converged
  
  double precision               :: u_dot_v, u_dot_u
  
  integer                        :: k_pairs, kl
  
  integer                        :: iter2
  double precision, allocatable  :: W(:,:),  U(:,:), S(:,:), overlap(:,:)
  double precision, allocatable  :: y(:,:), h(:,:), lambda(:), s2(:)
  double precision, allocatable  :: c(:), s_(:,:), s_tmp(:,:)
  double precision               :: diag_h_mat_elem
  double precision, allocatable  :: residual_norm(:)
  character*(16384)              :: write_buffer
  double precision               :: to_print(3,N_st)
  double precision               :: cpu, wall
  integer                        :: shift, shift2, itermax
  double precision               :: r1, r2
  logical                        :: state_ok(N_st_diag*davidson_sze_max)
  include 'constants.include.F'
 print*, 'dim_in',dim_in
 print*, 'sze ',sze
  
  !DIR$ ATTRIBUTES ALIGN : $IRP_ALIGN :: U, W, S, y, h, lambda
  if (N_st_diag*3 > sze) then
    print *,  'error in Davidson :'
    print *,  'Increase n_det_max_jacobi to ', N_st_diag*3
    stop -1
  endif
  
  integer, external              :: align_double
  itermax = max(3,min(davidson_sze_max, sze/N_st_diag))
  
  call write_time(iunit)
  call wall_time(wall)
  call cpu_time(cpu)
  write(iunit,'(A)') ''
  write(iunit,'(A)') 'Davidson Diagonalization'
  write(iunit,'(A)') '------------------------'
  write(iunit,'(A)') ''
  call write_int(iunit,N_st,'Number of states')
  call write_int(iunit,N_st_diag,'Number of states in diagonalization')
  call write_int(iunit,sze,'Number of determinants')
  r1 = 8.d0*(3.d0*dble(sze*N_st_diag*itermax+5.d0*(N_st_diag*itermax)**2 & 
    + 4.d0*(N_st_diag*itermax)+nproc*(4.d0*N_det_alpha_unique+2.d0*N_st_diag*sze)))/(1024.d0**3)
  call write_double(iunit, r1, 'Memory(Gb)')
  write(iunit,'(A)') ''
  write_buffer = '===== '
  do i=1,N_st
    write_buffer = trim(write_buffer)//' ================ =========== ==========='
  enddo
  write(iunit,'(A)') trim(write_buffer)
  write_buffer = ' Iter'
  do i=1,N_st
    write_buffer = trim(write_buffer)//'      Energy          S^2      Residual  '
  enddo
  write(iunit,'(A)') trim(write_buffer)
  write_buffer = '===== '
  do i=1,N_st
    write_buffer = trim(write_buffer)//' ================ =========== ==========='
  enddo
  write(iunit,'(A)') trim(write_buffer)
  

  allocate(                                                          &
      ! Large
      W(sze,N_st_diag*itermax),                                    &
      U(sze,N_st_diag*itermax),                                    &
      S(sze,N_st_diag*itermax),                                    &

      ! Small
      h(N_st_diag*itermax,N_st_diag*itermax),                        &
      y(N_st_diag*itermax,N_st_diag*itermax),                        &
      s_(N_st_diag*itermax,N_st_diag*itermax),                       &
      s_tmp(N_st_diag*itermax,N_st_diag*itermax),                    &
      residual_norm(N_st_diag),                                      &
      c(N_st_diag*itermax),                                          &
      s2(N_st_diag*itermax),                                         &
      overlap(N_st_diag*itermax, N_st_diag*itermax),                 &
      lambda(N_st_diag*itermax))
  
  h = 0.d0
  U = 0.d0
  W = 0.d0
  S = 0.d0
  y = 0.d0
  s_ = 0.d0
  s_tmp = 0.d0


  ASSERT (N_st > 0)
  ASSERT (N_st_diag >= N_st)
  ASSERT (sze > 0)
  ASSERT (Nint > 0)
  ASSERT (Nint == N_int)
  
  ! Davidson iterations
  ! ===================
  
  converged = .False.
  
  do k=N_st+1,N_st_diag
    u_in(k,k) = 10.d0
    do i=1,sze
      call random_number(r1)
      call random_number(r2)
      r1 = dsqrt(-2.d0*dlog(r1))
      r2 = dtwo_pi*r2
      u_in(i,k) = r1*dcos(r2)
    enddo
  enddo
  do k=1,N_st_diag
    call normalize(u_in(1,k),sze)
  enddo

  
  do while (.not.converged)
    
    do k=1,N_st_diag
      do i=1,sze
        U(i,k) = u_in(i,k)
      enddo
    enddo
    
    do iter=1,itermax-1
      
      shift  = N_st_diag*(iter-1)
      shift2 = N_st_diag*iter
      
      call ortho_qr(U,size(U,1),sze,shift2)

      ! Compute |W_k> = \sum_i |i><i|H|u_k>
      ! -----------------------------------------
      
       
      
      call apply_H_superci_to_vector(U(1,shift+1),W(1,shift+1))
      
      ! Compute h_kl = <u_k | W_l> = <u_k| H |u_l>
      ! -------------------------------------------

      call dgemm('T','N', shift2, shift2, sze,                       &
          1.d0, U, size(U,1), W, size(W,1),                          &
          0.d0, h, size(h,1))
      
      
      ! Diagonalize h
      ! -------------

      call lapack_diag(lambda,y,h,size(h,1),shift2)
      

      ! Express eigenvectors of h in the determinant basis
      ! --------------------------------------------------
      
      call dgemm('N','N', sze, N_st_diag, shift2,                    &
          1.d0, U, size(U,1), y, size(y,1), 0.d0, U(1,shift2+1), size(U,1))
      call dgemm('N','N', sze, N_st_diag, shift2,                    &
          1.d0, W, size(W,1), y, size(y,1), 0.d0, W(1,shift2+1), size(W,1))
      call dgemm('N','N', sze, N_st_diag, shift2,                    &
          1.d0, S, size(S,1), y, size(y,1), 0.d0, S(1,shift2+1), size(S,1))

      ! Compute residual vector and davidson step
      ! -----------------------------------------
      
      do k=1,N_st_diag
        do i=1,sze
          U(i,shift2+k) =  &
            (lambda(k) * U(i,shift2+k) - W(i,shift2+k) )      &
            /max(H_jj(i) - lambda (k),1.d-2)
        enddo

        if (k <= N_st) then
          residual_norm(k) = u_dot_u(U(1,shift2+k),sze)
          to_print(1,k) = lambda(k) + nuclear_repulsion
          to_print(3,k) = residual_norm(k)
        endif
      enddo
      
      write(iunit,'(1X,I3,1X,100(1X,F16.10,1X,F11.6,1X,E11.3))')  iter, to_print(1:3,1:N_st)
      call davidson_converged(lambda,residual_norm,wall,iter,cpu,N_st,converged)
      do k=1,N_st
        if (residual_norm(k) > 1.e8) then
        print *,  ''
          stop 'Davidson failed'
        endif
      enddo
      if (converged) then
        exit
      endif
      
    enddo

    ! Re-contract to u_in
    ! -----------
    
    call dgemm('N','N', sze, N_st_diag, shift2, 1.d0,      &
        U, size(U,1), y, size(y,1), 0.d0, u_in, size(u_in,1))

  enddo

  do k=1,N_st_diag
    energies(k) = lambda(k)
  enddo
  write_buffer = '===== '
  do i=1,N_st
    write_buffer = trim(write_buffer)//' ================ =========== ==========='
  enddo
  write(iunit,'(A)') trim(write_buffer)
  write(iunit,'(A)') ''
  call write_time(iunit)

  deallocate (                                                       &
      W, residual_norm,                                              &
      U, overlap,                                                    &
      c, S,                                                          &
      h,                                                             &
      y, s_, s_tmp,                                                  &
      lambda                                                         &
      )
end

