
 BEGIN_PROVIDER [ integer, dirac_mo_tot_num ]
  implicit none
  BEGIN_DOC
  !Concatenation of the large_mo_tot_num and small_mo_tot_num 
  END_DOC
  dirac_mo_tot_num = large_mo_tot_num + small_mo_tot_num
 END_PROVIDER

 BEGIN_PROVIDER [ integer, twice_dirac_mo_tot_num ]
  implicit none
  BEGIN_DOC
  ! twice the number of dirac MOs
  END_DOC
  logical                        :: has
  PROVIDE ezfio_filename
  if (mpi_master) then
    call ezfio_has_dirac_mo_basis_twice_dirac_mo_tot_num(has)
  endif
  IRP_IF MPI
    include 'mpif.h'
    integer                        :: ierr
    call MPI_BCAST( has, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
    if (ierr /= MPI_SUCCESS) then
      stop 'Unable to read mo_tot_num with MPI'
    endif
  IRP_ENDIF
  if (.not.has) then
    twice_dirac_mo_tot_num = 2*dirac_mo_tot_num 
    call ezfio_set_dirac_mo_basis_twice_dirac_mo_tot_num(twice_dirac_mo_tot_num)
  else
    if (mpi_master) then
      call ezfio_get_dirac_mo_basis_twice_dirac_mo_tot_num(twice_dirac_mo_tot_num)
    endif
    IRP_IF MPI
      call MPI_BCAST( mo_tot_num, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
      if (ierr /= MPI_SUCCESS) then
        stop 'Unable to read mo_tot_num with MPI'
      endif
    IRP_ENDIF
  endif
 !call write_int(6,twice_dirac_mo_tot_num,'twice_dirac_mo_tot_num')
  ASSERT (twice_dirac_mo_tot_num > 0)
 END_PROVIDER
