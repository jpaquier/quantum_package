
subroutine apply_H_superci_to_vector(u0,u1)
 implicit none
 double precision, intent(in) :: u0(size_super_ci)
 double precision, intent(out) :: u1(size_super_ci)

 integer :: i,iorb,j,jorb,k,korb,l,lorb
 integer :: index_i,index_j
 double precision :: dsqrt_2
 dsqrt_2 = dsqrt(2.d0)

 u1 = 0.d0
 
 if(type_of_superci == 0)then
   do j = 1, n_virt_orb
    jorb = list_virt(j)
    do i = 1, n_core_inact_orb
     iorb = list_core_inact(i)
     index_i = index_rotation_CI(i,j) 
     ! Diagonal and Brillouin matrix elements 
     u1(1) += u0(index_i) * dsqrt_2 * Fock_matrix_alpha_beta_average_mo(iorb,jorb)
     u1(index_i) += u0(1) * dsqrt_2 * Fock_matrix_alpha_beta_average_mo(iorb,jorb)
     u1(index_i) += u0(index_i) * diagonal_superci_matrix(index_i)
     ! Interaction through the virt-virt Fock operator
     do k = j+1, n_virt_orb
      korb = list_virt(k)
      index_j = index_rotation_CI(i,k)
      u1(index_i) += u0(index_j) * Fock_matrix_alpha_beta_average_mo(jorb,korb)
      u1(index_j) += u0(index_i) * Fock_matrix_alpha_beta_average_mo(jorb,korb)
     enddo
     ! Interaction through the core-core Fock operator
     do k = i+1, n_core_inact_orb
      korb = list_core_inact(k)
      index_j = index_rotation_CI(k,j)
      u1(index_i) -= u0(index_j) *  Fock_matrix_alpha_beta_average_mo(iorb,korb)
      u1(index_j) -= u0(index_i) *  Fock_matrix_alpha_beta_average_mo(iorb,korb)
     enddo
    enddo
   enddo

 else if(type_of_superci == 1)then

  do i = 1, n_core_inact_orb
   iorb = list_core_inact(i)
   do j = 1, n_virt_orb
    jorb = list_virt(j)
    index_i = index_rotation_CI(i,j) 
    ! Diagonal and Brillouin matrix elements 
    u1(1) += u0(index_i) * dsqrt_2 * Fock_matrix_alpha_beta_average_mo(iorb,jorb)
    u1(index_i) += u0(1) * dsqrt_2 * Fock_matrix_alpha_beta_average_mo(iorb,jorb)
    u1(index_i) += u0(index_i) * diagonal_superci_matrix(index_i)
    ! Interaction through the virt-virt Fock operator
    do k = j+1, n_virt_orb
     korb = list_virt(k)
     index_j = index_rotation_CI(i,k)
     u1(index_i) += u0(index_j) * (Fock_matrix_alpha_beta_average_mo(jorb,korb) &
                                                                   - transformed_occ1_virt2_virt2(i,k,j) + 2.d0 * transformed_occ1_virt1_occ2_virt2(i,k,i,j))
     u1(index_j) += u0(index_i) * (Fock_matrix_alpha_beta_average_mo(jorb,korb) & 
                                                                   - transformed_occ1_virt2_virt2(i,k,j) + 2.d0 * transformed_occ1_virt1_occ2_virt2(i,k,i,j))
    enddo
    ! Interaction through the core-core Fock operator
    do k = i+1, n_core_inact_orb
     korb = list_core_inact(k)
     index_j = index_rotation_CI(k,j)
     u1(index_i) -= u0(index_j) *  (Fock_matrix_alpha_beta_average_mo(iorb,korb) & 
                                                                     - transformed_virt1_occ2_occ2(j,k,i) + 2.d0 * transformed_occ1_virt1_occ2_virt2(i,j,k,j))
     u1(index_j) -= u0(index_i) *  (Fock_matrix_alpha_beta_average_mo(iorb,korb) & 
                                                                     - transformed_virt1_occ2_occ2(j,k,i) + 2.d0 * transformed_occ1_virt1_occ2_virt2(i,j,k,j))
    enddo
    
   enddo
  enddo

 else if(type_of_superci == 2)then

  do i = 1, n_core_inact_orb
   iorb = list_core_inact(i)
   do j = 1, n_virt_orb
    jorb = list_virt(j)
    index_i = index_rotation_CI(i,j) 
    ! Diagonal and Brillouin matrix elements 
    u1(1) += u0(index_i) * dsqrt_2 * Fock_matrix_alpha_beta_average_mo(iorb,jorb)
    u1(index_i) += u0(1) * dsqrt_2 * Fock_matrix_alpha_beta_average_mo(iorb,jorb)
    u1(index_i) += u0(index_i) * diagonal_superci_matrix(index_i)
    ! Interaction through the virt-virt Fock operator
    do k = j+1, n_virt_orb
     korb = list_virt(k)
     index_j = index_rotation_CI(i,k)
     u1(index_i) += u0(index_j) * (Fock_matrix_alpha_beta_average_mo(jorb,korb) &
                                                                   - transformed_occ1_virt2_virt2(i,k,j) + 2.d0 * transformed_occ1_virt1_occ2_virt2(i,k,i,j))
     u1(index_j) += u0(index_i) * (Fock_matrix_alpha_beta_average_mo(jorb,korb) & 
                                                                   - transformed_occ1_virt2_virt2(i,k,j) + 2.d0 * transformed_occ1_virt1_occ2_virt2(i,k,i,j))
    enddo
    ! Interaction through the core-core Fock operator
    do k = i+1, n_core_inact_orb
     korb = list_core_inact(k)
     index_j = index_rotation_CI(k,j)
     u1(index_i) -= u0(index_j) *  (Fock_matrix_alpha_beta_average_mo(iorb,korb) & 
                                                                     - transformed_virt1_occ2_occ2(j,k,i) + 2.d0 * transformed_occ1_virt1_occ2_virt2(i,j,k,j))
     u1(index_j) -= u0(index_i) *  (Fock_matrix_alpha_beta_average_mo(iorb,korb) & 
                                                                     - transformed_virt1_occ2_occ2(j,k,i) + 2.d0 * transformed_occ1_virt1_occ2_virt2(i,j,k,j))
    enddo
    ! Hole-particle interaction 
    do l = 1, j-1
     do k = 1, i-1
      index_j = index_rotation_CI(k,l)
      u1(index_i) += u0(index_j) * ( 2.d0 * transformed_occ1_virt1_occ2_virt2(k,l,i,j) - transformed_virt1_virt1_occ2_occ2(l,j,k,i) )
     enddo
    enddo

    do l = 1, j-1
     do k = i+1, n_core_inact_orb
      index_j = index_rotation_CI(k,l)
      u1(index_i) += u0(index_j) * ( 2.d0 * transformed_occ1_virt1_occ2_virt2(k,l,i,j) - transformed_virt1_virt1_occ2_occ2(l,j,k,i) )
     enddo
    enddo

    do l = j+1, n_virt_orb
     do k = 1, i-1
      index_j = index_rotation_CI(k,l)
      u1(index_i) += u0(index_j) * ( 2.d0 * transformed_occ1_virt1_occ2_virt2(k,l,i,j) - transformed_virt1_virt1_occ2_occ2(l,j,k,i) )
      superci_matrix(index_rotation_CI(i,j),index_rotation_CI(k,l)) = 2.d0 * transformed_occ1_virt1_occ2_virt2(k,l,i,j) - transformed_virt1_virt1_occ2_occ2(l,j,k,i)
     enddo
    enddo

    do l = j+1, n_virt_orb
     do k = i+1, n_core_inact_orb
      index_j = index_rotation_CI(k,l)
      u1(index_i) += u0(index_j) * ( 2.d0 * transformed_occ1_virt1_occ2_virt2(k,l,i,j) - transformed_virt1_virt1_occ2_occ2(l,j,k,i) )
     enddo
    enddo
    
   enddo
  enddo
 endif
end 



subroutine create_guess_super_ci(u_guess,e_guess)
 implicit none
 integer :: i,j
 double precision, intent(out) :: u_guess(size_super_ci),e_guess
 double precision, allocatable :: u0(:),u1(:),u2(:)
 allocate(u0(size_super_ci),u1(size_super_ci),u2(size_super_ci))
 u0 = 0.d0
 u0(1) = 1.d0 
 call apply_H_superci_to_vector(u0,u1)
 double precision :: interaction,e1,u_dot_v,delta_e,c1
 double precision :: norm
 norm = u_dot_v(u0,u1,size_super_ci)
 u1 = u1 - norm * u0 

 norm = u_dot_v(u1,u1,size_super_ci)
 norm = 1.d0/dsqrt(norm)
 u1 = u1 * norm
 call apply_H_superci_to_vector(u1,u2)
 interaction = u_dot_v(u0,u2,size_super_ci)
 e1 = u_dot_v(u2,u1,size_super_ci)
 delta_e = e1 
 if(delta_e > 0.d0)then
  e_guess =  0.5d0 * (delta_e - dsqrt(delta_e * delta_e + 4.d0 * interaction **2 ))
 else 
  e_guess =  0.5d0 * (delta_e + dsqrt(delta_e * delta_e + 4.d0 * interaction **2 ))
 endif
 if(interaction.gt.1.d-10)then
  c1 = e_guess/interaction 
 else 
  c1 = 0.d0
 endif

 u_guess = u0 + c1 * u1
 norm = u_dot_v(u_guess,u_guess,size_super_ci)
 norm = 1.d0/dsqrt(norm)
 u_guess = u_guess * norm

 double precision :: e0
 call apply_H_superci_to_vector(u_guess,u1)
 e0 = u_dot_v(u_guess,u1,size_super_ci) 
 if(dabs(e0 - e_guess).gt.1.d-10)then
  print*, 'PB in create_guess_super_ci'
  print*, e0,e_guess
 stop
 endif
 deallocate(u1,u0)

end


subroutine davidson_diag_general(u_in,energies,dim_in,sze,N_st,N_st_diag,Nint,iunit)
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
  double precision, intent(out)  :: energies(N_st)
  double precision, allocatable  :: H_jj(:)
  
  integer                        :: i
  ASSERT (N_st > 0)
  ASSERT (sze > 0)
  ASSERT (Nint > 0)
  ASSERT (Nint == N_int)
  PROVIDE mo_bielec_integrals_in_map
  allocate(H_jj(sze))
  
  do i=1,sze
    H_jj(i) = diagonal_superci_matrix(i)
    print*, H_jj(i),u_in(i,1)
  enddo

  call davidson_diag_hjj_general(u_in,H_jj,energies,dim_in,sze,N_st,N_st_diag,Nint,iunit)
  deallocate (H_jj)
end




subroutine davidson_diag_hjj_general(u_in,H_jj,energies,dim_in,sze,N_st,N_st_diag,Nint,iunit)
  implicit none
  BEGIN_DOC
  ! Davidson diagonalization with specific diagonal elements of the H matrix
  !
  ! H_jj : specific diagonal H matrix elements to diagonalize de Davidson
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
  ! N_st_diag : Number of states in which H is diagonalized
  !
  ! iunit : Unit for the I/O
  !
  ! Initial guess vectors are not necessarily orthonormal
  END_DOC
  integer, intent(in)            :: dim_in, sze, N_st, N_st_diag, Nint
  double precision,  intent(in)  :: H_jj(sze)
  integer,  intent(in)  :: iunit
  double precision, intent(inout) :: u_in(dim_in,N_st_diag)
  double precision, intent(out)  :: energies(N_st_diag)
  
  integer                        :: iter
  integer                        :: i,j,k,l,m
  logical                        :: converged
  
  double precision, allocatable  :: overlap(:,:)
  double precision               :: u_dot_v, u_dot_u
  
  integer, allocatable           :: kl_pairs(:,:)
  integer                        :: k_pairs, kl
  
  integer                        :: iter2
  double precision, allocatable  :: W(:,:,:),  U(:,:,:), R(:,:)
  double precision, allocatable  :: y(:,:,:,:), h(:,:,:,:), lambda(:)
  double precision, allocatable  :: c(:), H_small(:,:)
  double precision, allocatable  :: residual_norm(:)
  character*(16384)              :: write_buffer
  double precision               :: to_print(2,N_st)
  double precision               :: cpu, wall
  include 'constants.include.F'
  

  !DIR$ ATTRIBUTES ALIGN : $IRP_ALIGN :: U, W, R, y, h, lambda

  PROVIDE nuclear_repulsion

  call write_time(iunit)
  call wall_time(wall)
  call cpu_time(cpu)
  write(iunit,'(A)') ''
  write(iunit,'(A)') 'Davidson Diagonalization'
  write(iunit,'(A)') '------------------------'
  write(iunit,'(A)') ''
  call write_int(iunit,N_st,'Number of states')
  call write_int(iunit,N_st_diag,'Number of states in diagonalization')
  call write_int(iunit,sze,'Number of components ')
  write(iunit,'(A)') ''
  write_buffer = '===== '
  do i=1,N_st
    write_buffer = trim(write_buffer)//' ================ ================'
  enddo
  write(iunit,'(A)') trim(write_buffer)
  write_buffer = ' Iter'
  do i=1,N_st
    write_buffer = trim(write_buffer)//'            Energy           Residual'
  enddo
  write(iunit,'(A)') trim(write_buffer)
  write_buffer = '===== '
  do i=1,N_st
    write_buffer = trim(write_buffer)//' ================ ================'
  enddo
  write(iunit,'(A)') trim(write_buffer)

  integer, external :: align_double

  allocate(                                                          &
      kl_pairs(2,N_st_diag*(N_st_diag+1)/2),                         &
      W(sze,N_st_diag,davidson_sze_max),                           &
      U(sze,N_st_diag,davidson_sze_max),                           &
      R(sze,N_st_diag),                                            &
      h(N_st_diag,davidson_sze_max,N_st_diag,davidson_sze_max),      &
      y(N_st_diag,davidson_sze_max,N_st_diag,davidson_sze_max),      &
      residual_norm(N_st_diag),                                      &
      overlap(N_st_diag,N_st_diag),                                  &
      c(N_st_diag*davidson_sze_max),                                 &
      H_small(N_st_diag,N_st_diag),                                  &
      lambda(N_st_diag*davidson_sze_max))
  
  ASSERT (N_st > 0)
  ASSERT (N_st_diag >= N_st)
  ASSERT (sze > 0)
  ASSERT (Nint > 0)
  ASSERT (Nint == N_int)
  
  ! Davidson iterations
  ! ===================
  
  converged = .False.
  
 !do k=1,N_st_diag

 ! !if (k > N_st) then
 ! !  do i=1,sze
 ! !    double precision               :: r1, r2
 ! !    call random_number(r1)
 ! !    call random_number(r2)
 ! !    u_in(i,k) = dsqrt(-2.d0*dlog(r1))*dcos(dtwo_pi*r2)
 ! !  enddo
 ! !endif
 !  
 ! !! Gram-Schmidt
 ! !! ------------
 ! !call dgemv('T',sze,k-1,1.d0,u_in,size(u_in,1),                   &
 ! !    u_in(1,k),1,0.d0,c,1)
 ! !call dgemv('N',sze,k-1,-1.d0,u_in,size(u_in,1),                  &
 ! !    c,1,1.d0,u_in(1,k),1)
 ! !call normalize(u_in(1,k),sze)
 !enddo


  
  do while (.not.converged)
    
    do k=1,N_st_diag
      do i=1,sze
        U(i,k,1) = u_in(i,k)
      enddo
    enddo

    do iter=1,davidson_sze_max-1
      
      ! Compute |W_k> = \sum_i |i><i|H|u_k>
      ! -----------------------------------------
      
      call apply_H_superci_to_vector(U(1,1,iter),W(1,1,iter))

      call dgemm('T','N', N_st_diag*iter, N_st_diag, sze,            &
          1.d0, U, size(U,1), W(1,1,iter), size(W,1),                &
          0.d0, h(1,1,1,iter), size(h,1)*size(h,2))

      ! Diagonalize h
      ! -------------
      call lapack_diag(lambda,y,h,N_st_diag*davidson_sze_max,N_st_diag*iter)
      
      ! Express eigenvectors of h in the determinant basis
      ! --------------------------------------------------
      
      do k=1,N_st_diag
        do i=1,sze
          U(i,k,iter+1) = 0.d0
          W(i,k,iter+1) = 0.d0
        enddo
      enddo
!
!
      call dgemm('N','N', sze, N_st_diag, N_st_diag*iter,            &
          1.d0, U, size(U,1), y, size(y,1)*size(y,2), 0.d0, U(1,1,iter+1), size(U,1))
      call dgemm('N','N',sze,N_st_diag,N_st_diag*iter,               &
          1.d0, W, size(W,1), y, size(y,1)*size(y,2), 0.d0, W(1,1,iter+1), size(W,1))


      ! Compute residual vector
      ! -----------------------
      
      do k=1,N_st_diag
        do i=1,sze
          R(i,k) = lambda(k) * U(i,k,iter+1) - W(i,k,iter+1)
        enddo
        if (k <= N_st) then
          residual_norm(k) = u_dot_u(R(1,k),sze)
          to_print(1,k) = lambda(k) + nuclear_repulsion
          to_print(2,k) = residual_norm(k)
        endif
      enddo
      
      write(iunit,'(1X,I3,1X,100(1X,F16.10,1X,E16.6))')  iter, to_print(:,1:N_st)
      call davidson_converged(lambda,residual_norm,wall,iter,cpu,N_st,converged)
      if (converged) then
        exit
      endif
      
      ! Davidson step
      ! -------------
      
      do k=1,N_st_diag
        do i=1,sze
          U(i,k,iter+1) = -1.d0/max(H_jj(i) - lambda(k),1.d-2) * R(i,k)
        enddo
      enddo
      
      ! Gram-Schmidt
      ! ------------
      
      do k=1,N_st_diag

        call dgemv('T',sze,N_st_diag*iter,1.d0,U,size(U,1),  &
              U(1,k,iter+1),1,0.d0,c,1)
        call dgemv('N',sze,N_st_diag*iter,-1.d0,U,size(U,1), &
              c,1,1.d0,U(1,k,iter+1),1)
!
        call dgemv('T',sze,k-1,1.d0,U(1,1,iter+1),size(U,1),   &
            U(1,k,iter+1),1,0.d0,c,1)
        call dgemv('N',sze,k-1,-1.d0,U(1,1,iter+1),size(U,1),        &
            c,1,1.d0,U(1,k,iter+1),1)

        call normalize( U(1,k,iter+1), sze )
      enddo

    enddo

    if (.not.converged) then
      iter = davidson_sze_max-1
    endif
    
    ! Re-contract to u_in
    ! -----------
    
    do k=1,N_st_diag
      energies(k) = lambda(k)
      do i=1,sze
        u_in(i,k) = 0.d0
      enddo
    enddo

    call dgemm('N','N', sze, N_st_diag, N_st_diag*iter, 1.d0,      &
        U, size(U,1), y, N_st_diag*davidson_sze_max, &
        0.d0, u_in, size(u_in,1))

  enddo

  write_buffer = '===== '
  do i=1,N_st
    write_buffer = trim(write_buffer)//' ================ ================'
  enddo
  write(iunit,'(A)') trim(write_buffer)
  write(iunit,'(A)') ''
  call write_time(iunit)

  deallocate (                                                       &
      kl_pairs,                                                      &
      W, residual_norm,                                              &
      U, overlap,                                                    &
      R, c,                                                          &
      h,                                                             &
      y,                                                             &
      lambda                                                         &
      )
end

