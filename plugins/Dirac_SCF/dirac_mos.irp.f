 BEGIN_PROVIDER [ complex*16, dirac_mo_coef, (2*(dirac_ao_num),2*(dirac_mo_tot_num))
 implicit none
  BEGIN_DOC
  ! Molecular orbital coefficients on AO basis set
  ! dirac_mo_coef(i,j) = coefficient of the ith ao on the jth mo
  ! mo_label : Label characterizing the MOS (local, canonical, natural,
  ! etc)
  END_DOC
  integer                        :: i,i_minus,j,j_minus
  PROVIDE ezfio_filename
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
 END_PROVIDER


 subroutine dirac_ao_to_mo(A_ao,LDA_ao,A_mo,LDA_mo)
  implicit none
  BEGIN_DOC
  ! Transform A from the AO basis to the MO basis
  !
  ! Ct.A_ao.C
  END_DOC
  integer, intent(in)            :: LDA_ao,LDA_mo
  complex*16, intent(in)         :: A_ao(LDA_ao,2*dirac_ao_num)
  double precision, intent(out)  :: A_mo(LDA_mo,2*dirac_mo_tot_num)
  complex*16, allocatable        :: T(:,:)
  allocate ( T(2*(dirac_ao_num),2*(dirac_mo_tot_num)) )
  call zgemm('N','N', 2*(dirac_ao_num), 2*(dirac_mo_tot_num), 2*(dirac_ao_num), &
      (1.d0,0.d0), A_ao,LDA_ao,                                             &
      dirac_mo_coef, size(dirac_mo_coef,1),                          &
      (0.d0,0.d0), T, size(T,1))
  call zgemm('T','N', 2*(dirac_mo_tot_num), 2*(dirac_mo_tot_num), 2*(dirac_ao_num), &
      (1.d0,0.d0), dirac_mo_coef,size(dirac_mo_coef,1),                     &
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
  ! (S.C).A_mo.(S.C)t
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
    call ezfio_has_mo_basis_mo_label(exists)
    if (exists) then
      call ezfio_get_mo_basis_mo_label(mo_label)
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

