 BEGIN_PROVIDER [ double precision, dirac_mo_coef_re ,(twice_dirac_ao_num,twice_dirac_mo_tot_num) ]
  implicit none
  BEGIN_DOC
  ! Real part of the coefficient of the ith ao on the jth mo
  END_DOC
  logical                        :: has
  PROVIDE ezfio_filename
  if (mpi_master) then
    if (size(dirac_mo_coef_re) == 0) return
    call ezfio_has_dirac_mo_basis_dirac_mo_coef_re(has)
    if (has) then
     call ezfio_get_dirac_mo_basis_dirac_mo_coef_re(dirac_mo_coef_re)
    endif
  endif
  IRP_IF MPI
    include 'mpif.h'
    integer :: ierr
    call MPI_BCAST( dirac_mo_coef_re, (twice_dirac_ao_num)*(twice_dirac_mo_tot_num), MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
    if (ierr /= MPI_SUCCESS) then
      stop 'Unable to read dirac_mo_coef_re with MPI'
    endif
  IRP_ENDIF
  call write_time(6)
  if (mpi_master) then
    write(6, *) 'Read  dirac_mo_coef_re'
  endif
 END_PROVIDER

 BEGIN_PROVIDER [ double precision, dirac_mo_coef_im ,(twice_dirac_ao_num,twice_dirac_mo_tot_num) ]
  implicit none
  BEGIN_DOC
  ! Imaginary part of the coefficient of the ith ao on the jth mo
  END_DOC
  logical                        :: has
  PROVIDE ezfio_filename
  if (mpi_master) then
    if (size(dirac_mo_coef_im) == 0) return

    call ezfio_has_dirac_mo_basis_dirac_mo_coef_im(has)
    if (has) then
      call ezfio_get_dirac_mo_basis_dirac_mo_coef_im(dirac_mo_coef_im)
    endif
  endif
  IRP_IF MPI
    include 'mpif.h'
    integer :: ierr
    call MPI_BCAST( dirac_mo_coef_im, (twice_dirac_ao_num)*(twice_dirac_mo_tot_num),MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
    if (ierr /= MPI_SUCCESS) then
      stop 'Unable to read dirac_mo_coef_im with MPI'
    endif
  IRP_ENDIF
  call write_time(6)
  if (mpi_master) then
    write(6, *) 'Read  dirac_mo_coef_im'
  endif
 END_PROVIDER


 BEGIN_PROVIDER [ complex*16, dirac_mo_coef, (2*dirac_ao_num,2*dirac_mo_tot_num) ]
  implicit none
  BEGIN_DOC
  ! Molecular orbital coefficients on AO basis set
  ! dirac_mo_coef(i,j) = coefficient of the ith ao on the jth mo
  ! dirac_mo_label : Label characterizing the MOS (local, canonical, natural, etc)
  END_DOC
  integer                        :: i,i_minus,j,j_minus
  double precision, allocatable  :: buffer(:,:)
  logical                        :: exists_Re,exists_Im
  PROVIDE ezfio_filename 
  if (mpi_master) then
    ! Coefs
    call ezfio_has_dirac_mo_basis_dirac_mo_coef_Re(exists_Re)
    call ezfio_has_dirac_mo_basis_dirac_mo_coef_Im(exists_Im)
  endif
  IRP_IF MPI
    include 'mpif.h'
    integer :: ierr
    call MPI_BCAST(exists_Re, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
    if (ierr /= MPI_SUCCESS) then
      stop 'Unable to read dirac_mo_coef_Re with MPI'
    endif
    call MPI_BCAST(exists_Im, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
    if (ierr /= MPI_SUCCESS) then
      stop 'Unable to read dirac_mo_coef_Im with MPI'
    endif
  IRP_ENDIF
  if (exists_Re .and. exists_Im)  then
    if (mpi_master) then
      call ezfio_get_dirac_mo_basis_dirac_mo_coef_Re(dirac_mo_coef_Re)
     !write(*,*) 'Read  dirac_mo_coef_Re'
      call ezfio_get_dirac_mo_basis_dirac_mo_coef_Im(dirac_mo_coef_Im)
     !write(*,*) 'Read  dirac_mo_coef_Im'
      dirac_mo_coef = (1,0)*dirac_mo_coef_Re  + (0,1)*dirac_mo_coef_Im
    endif
    IRP_IF MPI
      call MPI_BCAST( dirac_mo_coef_Re, 2*dirac_mo_tot_num*2*dirac_ao_num, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
      if (ierr /= MPI_SUCCESS) then
        stop 'Unable to read dirac_mo_coef_Re with MPI'
      endif
      call MPI_BCAST( dirac_mo_coef_Im, 2*dirac_mo_tot_num*2*dirac_ao_num,MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
      if (ierr /= MPI_SUCCESS) then
        stop 'Unable to read dirac_mo_coef_Im with MPI'
      endif
    IRP_ENDIF
  else
    ! Orthonormalized AO basis
   dirac_mo_coef = (0.d0,0.d0)
   do j=1, 2*dirac_mo_tot_num
    if (j .le. large_mo_tot_num) then
     do i=1, 2*dirac_ao_num
      if (i .le. large_ao_num) then
       dirac_mo_coef(i,j) = (1.d0,0.d0)*large_mo_coef(i,j)
      endif
     enddo
    elseif (j .gt. large_mo_tot_num .and. j .le. 2*large_mo_tot_num) then
     j_minus = j - large_mo_tot_num
     do i=1, 2*dirac_ao_num
      if (i .gt. large_ao_num .and. i .le. 2*large_ao_num) then
       i_minus = i - large_ao_num
       dirac_mo_coef(i,j) = (1.d0,0.d0)*large_mo_coef(i_minus,j_minus)
      endif
     enddo
    elseif (j.gt. 2*large_mo_tot_num .and. j .le. (2*large_mo_tot_num+small_mo_tot_num)) then
     j_minus = j - 2*large_mo_tot_num
     do i=1, 2*dirac_ao_num
     i_minus = i - 2*large_ao_num
      if (i .gt. 2*large_ao_num .and. i .le. (2*large_ao_num+small_ao_num)) then
       dirac_mo_coef(i,j) = (1.d0,0.d0)*small_mo_coef(i_minus,j_minus)
      endif
     enddo
    elseif (j .gt. (2*large_ao_num + small_ao_num)) then
     j_minus = j - (2*large_mo_tot_num + small_mo_tot_num)
     do i=1, 2*(large_ao_num+small_ao_num)
      if (i .gt. (2*large_ao_num+small_ao_num)) then
       i_minus = i - (2*large_ao_num + small_ao_num)
       dirac_mo_coef(i,j) = (1.d0,0.d0)*small_mo_coef(i_minus,j_minus)
      endif
     enddo
    endif
   enddo
  endif
 END_PROVIDER


 BEGIN_PROVIDER [ complex*16, dirac_mo_coef_S, (2*(dirac_ao_num),2*(dirac_mo_tot_num))
 implicit none
  BEGIN_DOC
  !Molecular orbital coefficients on AO basis set diagonalizing the overlap matrix S
  !dirac_mo_coef_S(i,j) = coefficient of the ith ao on the jth mo
  END_DOC
  integer                        :: i,i_minus,j,j_minus
  PROVIDE ezfio_filename
  dirac_mo_coef_S = (0.d0,0.d0)
  do j=1, 2*dirac_mo_tot_num
   if (j .le. large_mo_tot_num) then
    do i=1, 2*dirac_ao_num
     if (i .le. large_ao_num) then
      dirac_mo_coef_S(i,j) = (1.d0,0.d0)*large_mo_coef(i,j)
     endif
    enddo
   elseif (j .gt. large_mo_tot_num .and. j .le. 2*large_mo_tot_num) then
    j_minus = j - large_mo_tot_num
    do i=1, 2*dirac_ao_num
     if (i .gt. large_ao_num .and. i .le. 2*large_ao_num) then
      i_minus = i - large_ao_num
      dirac_mo_coef_S(i,j) = (1.d0,0.d0)*large_mo_coef(i_minus,j_minus)
     endif
    enddo
   elseif (j.gt. 2*large_mo_tot_num .and. j .le. (2*large_mo_tot_num+small_mo_tot_num)) then
    j_minus = j - 2*large_mo_tot_num
    do i=1, 2*dirac_ao_num
    i_minus = i - 2*large_ao_num
     if (i .gt. 2*large_ao_num .and. i .le. (2*large_ao_num+small_ao_num)) then
      dirac_mo_coef_S(i,j) = (1.d0,0.d0)*small_mo_coef(i_minus,j_minus)
     endif
    enddo
   elseif (j .gt. (2*large_ao_num + small_ao_num)) then
    j_minus = j - (2*large_mo_tot_num + small_mo_tot_num)
    do i=1, 2*(large_ao_num+small_ao_num)
     if (i .gt. (2*large_ao_num+small_ao_num)) then
      i_minus = i - (2*large_ao_num + small_ao_num)
      dirac_mo_coef_S(i,j) = (1.d0,0.d0)*small_mo_coef(i_minus,j_minus)
     endif
    enddo
   endif
  enddo
 END_PROVIDER
 

!BEGIN_PROVIDER [ double precision, dirac_mo_coef_guess, (2*(dirac_ao_num),2*(dirac_mo_tot_num))
!implicit none
! BEGIN_DOC
! !Molecular orbital coefficients on AO basis set obtained from non-relativistic
! ! calculations and transformed as a relativistic guess
! !dirac_mo_coef(i,j) = coefficient of the ith ao on the jth mo
! END_DOC
! integer                        :: i,i_minus,i_plus,j,j_minus,j_plus
! PROVIDE ezfio_filename
! dirac_mo_coef_guess = 0.d0
! do j=1, 2*dirac_mo_tot_num
!  if (j .le. 2*small_mo_tot_num) then
!   j_plus = j + 2*large_ao_num
!   dirac_mo_coef_guess(j_plus,j) = 1.d0
!  elseif (j.gt. 2*small_mo_tot_num .and. j .le. (2*small_mo_tot_num+large_mo_tot_num)) then
!   j_minus = j - 2*small_mo_tot_num
!   do i=1, large_ao_num
!    dirac_mo_coef_guess(i,j) = mo_coef(i,j_minus)
!   enddo
!  elseif (j .gt. (2*small_mo_tot_num + large_mo_tot_num)) then
!   j_minus = j - 2*small_mo_tot_num 
!   dirac_mo_coef_guess(j_minus,j) = 1.d0
!  endif
! enddo
!END_PROVIDER

!BEGIN_PROVIDER [ complex*16, dirac_mo_coef_guess_ortho_canonical, (2*(dirac_ao_num),2*(dirac_mo_tot_num))
!implicit none
! BEGIN_DOC
! !Molecular orbital coefficients on AO basis set obtained from non-relativistic
! ! calculations and transformed as an ortho_canonical relativistic guess
! !dirac_mo_coef(i,j) = coefficient of the ith ao on the jth mo
! END_DOC
! integer                        :: i,i_minus,j,j_minus
! double precision               :: dirac_mo_coef_guess_ortho_canonical_tmp(2*(dirac_ao_num),2*(dirac_mo_tot_num))
! PROVIDE ezfio_filename
! dirac_mo_coef_guess_ortho_canonical_tmp = 0.d0
! do i = 1, 2*dirac_ao_num
!  dirac_mo_coef_guess_ortho_canonical_tmp(i,i) += 1.d0
! enddo
!!dirac_mo_ortho_canonical_num = 2*dirac_ao_num
! call ortho_canonical(dirac_mo_coef_guess,size(dirac_mo_coef_guess,1), 2*dirac_ao_num, dirac_mo_coef_guess_ortho_canonical_tmp, size(dirac_mo_coef_guess_ortho_canonical_tmp,1), 2*dirac_ao_num)
! dirac_mo_coef_guess_ortho_canonical = (1.d0,0.d0)*dirac_mo_coef_guess_ortho_canonical_tmp
!END_PROVIDER

 subroutine dirac_ao_to_mo(A_ao,LDA_ao,A_mo,LDA_mo)
  implicit none
  BEGIN_DOC
  ! Transform A from the AO basis to the MO basis
  !
  ! Ct.A_ao.C
  !
  ! For complex matrices
  END_DOC
  integer, intent(in)            :: LDA_ao,LDA_mo
  complex*16, intent(in)         :: A_ao(LDA_ao,2*dirac_ao_num)
  double precision, intent(out)  :: A_mo(LDA_mo,2*dirac_mo_tot_num)
  complex*16, allocatable        :: T(:,:)
  allocate ( T(2*(dirac_ao_num),2*(dirac_mo_tot_num)) )
  call zgemm('N','N', 2*(dirac_ao_num), 2*(dirac_mo_tot_num), 2*(dirac_ao_num), &
      (1.d0,0.d0), A_ao,LDA_ao,                                             &
      dirac_mo_coef_S, size(dirac_mo_coef_S,1),                          &
      (0.d0,0.d0), T, size(T,1))
  call zgemm('T','N', 2*(dirac_mo_tot_num), 2*(dirac_mo_tot_num), 2*(dirac_ao_num), &
      (1.d0,0.d0), dirac_mo_coef_S,size(dirac_mo_coef_S,1),                     &
      T, 2*(dirac_ao_num),                                    &
      (0.d0,0.d0), A_mo, size(A_mo,1))
  deallocate(T)
 end


 BEGIN_PROVIDER [ complex*16, dirac_S_mo_coef, (2*dirac_ao_num, 2*dirac_mo_tot_num) ]
  implicit none
  BEGIN_DOC
  ! Product S.C where S is the overlap matrix in the AO basis and C the mo_coef matrix.
  END_DOC
  call zgemm('N','N', 2*dirac_ao_num, 2*dirac_mo_tot_num, 2*dirac_ao_num,            &
     (1.d0,0.d0), dirac_ao_overlap,size(dirac_ao_overlap,1),                         &
     dirac_mo_coef, size(dirac_mo_coef,1),                                           &
     (0.d0,0.d0), dirac_S_mo_coef, size(dirac_S_mo_coef,1))
 END_PROVIDER
 

 subroutine dirac_mo_to_ao(A_mo,LDA_mo,A_ao,LDA_ao)
  implicit none
  BEGIN_DOC
  ! Transform A from the MO basis to the AO basis
  !
  ! (S.C).A_mo.(S.C)*
  !
  ! For complex matrices
  END_DOC
  integer, intent(in)            :: LDA_ao,LDA_mo
  complex*16, intent(in)         :: A_mo(LDA_mo,2*dirac_mo_tot_num)
  complex*16, intent(out)        :: A_ao(LDA_ao,2*dirac_ao_num)
  complex*16, allocatable        :: T(:,:)
  allocate ( T(2*dirac_mo_tot_num,2*dirac_ao_num) )
  call zgemm('N','T', 2*dirac_mo_tot_num, 2*dirac_ao_num, 2*dirac_mo_tot_num,        &
       (1.d0,0.d0), A_mo,size(A_mo,1),                                               &
       dirac_S_mo_coef, size(dirac_S_mo_coef,1),                                     &
       (0.d0,0.d0), T, size(T,1))
  call zgemm('N','N', 2*dirac_ao_num, 2*dirac_ao_num, 2*dirac_mo_tot_num,            &
       (1.d0,0.d0), dirac_S_mo_coef, size(dirac_S_mo_coef,1),                        &
       T, size(T,1),                                                                 &
       (0.d0,0.d0), A_ao, size(A_ao,1))
  deallocate(T)
 end


 BEGIN_PROVIDER [ complex*16, dirac_mo_overlap_bis,(2*dirac_mo_tot_num,2*dirac_mo_tot_num)]
  implicit none
  BEGIN_DOC
  !Array of the overlap matrix of the MOs
  END_DOC
  integer :: i,j,n,l
  double precision :: f
  integer :: lmax
  lmax = (2*dirac_ao_num/4) * 4
  do j=1,2*dirac_mo_tot_num
   do i= 1,2*dirac_mo_tot_num
    dirac_mo_overlap_bis(i,j) = 0.d0
    do n = 1, lmax,4
     do l = 1, 2*dirac_ao_num
      dirac_mo_overlap_bis(i,j) += dirac_mo_coef(l,i) * &
                              ( dirac_mo_coef(n  ,j) * dirac_ao_overlap(l,n  )  &
                              + dirac_mo_coef(n+1,j) * dirac_ao_overlap(l,n+1)  &
                              + dirac_mo_coef(n+2,j) * dirac_ao_overlap(l,n+2)  &
                              + dirac_mo_coef(n+3,j) * dirac_ao_overlap(l,n+3)  )
     enddo
    enddo
    do n = lmax+1, 2*dirac_ao_num
     do l = 1, 2*dirac_ao_num
      dirac_mo_overlap_bis(i,j) +=  dirac_mo_coef(n,j) * dirac_mo_coef(l,i) * dirac_ao_overlap(l,n)
     enddo
    enddo
   enddo
  enddo
 END_PROVIDER
 
 BEGIN_PROVIDER [complex*16, dirac_mo_overlap,(2*(dirac_mo_tot_num),2*(dirac_mo_tot_num))]
  implicit none
  BEGIN_DOC
  !Array of the overlap matrix of the MOs
  END_DOC
    call dirac_ao_to_mo(                                          &
        dirac_ao_overlap,                                         &
        size(dirac_ao_overlap,1),                                 &
        dirac_mo_overlap,                                         &
        size(dirac_mo_overlap,1)                                  &
        )
 END_PROVIDER
 

 BEGIN_PROVIDER[ character*(64), dirac_mo_label ]
  implicit none
  BEGIN_DOC
  ! Molecular orbital coefficients on AO basis set
  ! mo_coef(i,j) = coefficient of the ith ao on the jth mo
  ! mo_label : Label characterizing the MOS (local, canonical, natural, etc)
  END_DOC
  logical                        :: exists
  PROVIDE ezfio_filename
  if (mpi_master) then
    call ezfio_has_dirac_mo_basis_dirac_mo_label(exists)
    if (exists) then
      call ezfio_get_dirac_mo_basis_dirac_mo_label(dirac_mo_label)
      dirac_mo_label = trim(dirac_mo_label)
    else
      dirac_mo_label = 'no_label'
    endif
    write(*,*) '* dirac_mo_label          ', trim(dirac_mo_label)
  endif
 !IRP_IF MPI
 !  include 'mpif.h'
 !  integer :: ierr
 !  call MPI_BCAST( mo_label, 64, MPI_CHARACTER, 0, MPI_COMM_WORLD, ierr)
 !  if (ierr /= MPI_SUCCESS) then
 !    stop 'Unable to read mo_label with MPI'
 !  endif
 !IRP_ENDIF
 END_PROVIDER

!subroutine orthonormalize_dirac_mos
! implicit none
! BEGIN_DOC
! !Orthonormalized complex MOs
! END_DOC 
! call ortho_canonical_complex(dirac_mo_overlap,2*dirac_mo_tot_num,2*dirac_mo_tot_num,dirac_mo_coef,2*dirac_ao_num,2*dirac_ao_num)
! dirac_mo_label = 'Orthonormalized'
! SOFT_TOUCH dirac_mo_coef dirac_mo_label
!end

 

