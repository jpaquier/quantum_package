! DO NOT MODIFY BY HAND
! Created by $QP_ROOT/scripts/ezfio_interface/ei_handler.py
! from file /nas/home_lct/jpaquier/quantum_package/src/Dirac_Integrals_Bielec/EZFIO.cfg


BEGIN_PROVIDER [ double precision, dirac_ao_integrals_threshold  ]
  implicit none
  BEGIN_DOC
! If |<pq|rs>| < dirac_ao_integrals_threshold then <pq|rs> is zero
  END_DOC

  logical                        :: has
  PROVIDE ezfio_filename
  if (mpi_master) then
    
    call ezfio_has_dirac_integrals_bielec_dirac_threshold_ao(has)
    if (has) then
      call ezfio_get_dirac_integrals_bielec_dirac_threshold_ao(dirac_ao_integrals_threshold)
    else
      print *, 'dirac_integrals_bielec/dirac_threshold_ao not found in EZFIO file'
      stop 1
    endif
  endif
  IRP_IF MPI
    include 'mpif.h'
    integer :: ierr
    call MPI_BCAST( dirac_ao_integrals_threshold, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
    if (ierr /= MPI_SUCCESS) then
      stop 'Unable to read dirac_ao_integrals_threshold with MPI'
    endif
  IRP_ENDIF

  call write_time(6)
  if (mpi_master) then
    write(6, *) 'Read  dirac_ao_integrals_threshold'
  endif

END_PROVIDER

BEGIN_PROVIDER [ character*(32), dirac_interaction  ]
  implicit none
  BEGIN_DOC
! Type of electron-electron interaction used. Possibles choices are [ Coulomb | Coulomb_Gaunt | Coulomb_Breit]
  END_DOC

  logical                        :: has
  PROVIDE ezfio_filename
  if (mpi_master) then
    
    call ezfio_has_dirac_integrals_bielec_dirac_interaction(has)
    if (has) then
      call ezfio_get_dirac_integrals_bielec_dirac_interaction(dirac_interaction)
    else
      print *, 'dirac_integrals_bielec/dirac_interaction not found in EZFIO file'
      stop 1
    endif
  endif
  IRP_IF MPI
    include 'mpif.h'
    integer :: ierr
    call MPI_BCAST( dirac_interaction, 1*32, MPI_CHARACTER, 0, MPI_COMM_WORLD, ierr)
    if (ierr /= MPI_SUCCESS) then
      stop 'Unable to read dirac_interaction with MPI'
    endif
  IRP_ENDIF

  call write_time(6)
  if (mpi_master) then
    write(6, *) 'Read  dirac_interaction'
  endif

END_PROVIDER
