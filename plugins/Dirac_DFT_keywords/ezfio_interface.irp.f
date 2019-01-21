! DO NOT MODIFY BY HAND
! Created by $QP_ROOT/scripts/ezfio_interface/ei_handler.py
! from file /nas/home_lct/jpaquier/quantum_package/src/Dirac_DFT_keywords/EZFIO.cfg


BEGIN_PROVIDER [ character*(32), dirac_exchange_functional  ]
  implicit none
  BEGIN_DOC
! name of the relativistic exchange functional. Possibles choices are [ dirac_short_range_LDA_P2 | dirac_short_range_LDA_P4 | dirac_short_range_LDA_P6 ]
  END_DOC

  logical                        :: has
  PROVIDE ezfio_filename
  if (mpi_master) then
    
    call ezfio_has_dirac_dft_keywords_dirac_exchange_functional(has)
    if (has) then
      write(6,'(A)') '.. >>>>> [ IO READ: dirac_exchange_functional ] <<<<< ..'
      call ezfio_get_dirac_dft_keywords_dirac_exchange_functional(dirac_exchange_functional)
    else
      print *, 'dirac_dft_keywords/dirac_exchange_functional not found in EZFIO file'
      stop 1
    endif
  endif
  IRP_IF MPI_DEBUG
    print *,  irp_here, mpi_rank
    call MPI_BARRIER(MPI_COMM_WORLD, ierr)
  IRP_ENDIF
  IRP_IF MPI
    include 'mpif.h'
    integer :: ierr
    call MPI_BCAST( dirac_exchange_functional, 1*32, MPI_CHARACTER, 0, MPI_COMM_WORLD, ierr)
    if (ierr /= MPI_SUCCESS) then
      stop 'Unable to read dirac_exchange_functional with MPI'
    endif
  IRP_ENDIF

  call write_time(6)

END_PROVIDER

BEGIN_PROVIDER [ character*(32), dirac_correlation_functional  ]
  implicit none
  BEGIN_DOC
! name of the correlation relativistic functional
  END_DOC

  logical                        :: has
  PROVIDE ezfio_filename
  if (mpi_master) then
    
    call ezfio_has_dirac_dft_keywords_dirac_correlation_functional(has)
    if (has) then
      write(6,'(A)') '.. >>>>> [ IO READ: dirac_correlation_functional ] <<<<< ..'
      call ezfio_get_dirac_dft_keywords_dirac_correlation_functional(dirac_correlation_functional)
    else
      print *, 'dirac_dft_keywords/dirac_correlation_functional not found in EZFIO file'
      stop 1
    endif
  endif
  IRP_IF MPI_DEBUG
    print *,  irp_here, mpi_rank
    call MPI_BARRIER(MPI_COMM_WORLD, ierr)
  IRP_ENDIF
  IRP_IF MPI
    include 'mpif.h'
    integer :: ierr
    call MPI_BCAST( dirac_correlation_functional, 1*32, MPI_CHARACTER, 0, MPI_COMM_WORLD, ierr)
    if (ierr /= MPI_SUCCESS) then
      stop 'Unable to read dirac_correlation_functional with MPI'
    endif
  IRP_ENDIF

  call write_time(6)

END_PROVIDER
