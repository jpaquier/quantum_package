! DO NOT MODIFY BY HAND
! Created by $QP_ROOT/scripts/ezfio_interface/ei_handler.py
! from file /nas/home_lct/jpaquier/programs/quantum_package_master/src/Dirac_SCF/EZFIO.cfg


BEGIN_PROVIDER [ character*(32), dirac_correlation_functional  ]
  implicit none
  BEGIN_DOC
! name of the relativistic correlation functional
  END_DOC

  logical                        :: has
  PROVIDE ezfio_filename
  if (mpi_master) then
    
    call ezfio_has_dirac_scf_dirac_correlation_functional(has)
    if (has) then
      call ezfio_get_dirac_scf_dirac_correlation_functional(dirac_correlation_functional)
    else
      print *, 'dirac_scf/dirac_correlation_functional not found in EZFIO file'
      stop 1
    endif
  endif
  IRP_IF MPI
    include 'mpif.h'
    integer :: ierr
    call MPI_BCAST( dirac_correlation_functional, 1*32, MPI_CHARACTER, 0, MPI_COMM_WORLD, ierr)
    if (ierr /= MPI_SUCCESS) then
      stop 'Unable to read dirac_correlation_functional with MPI'
    endif
  IRP_ENDIF

  call write_time(6)
  if (mpi_master) then
    write(6, *) 'Read  dirac_correlation_functional'
  endif

END_PROVIDER

BEGIN_PROVIDER [ double precision, dirac_thresh_scf  ]
  implicit none
  BEGIN_DOC
! Threshold on the convergence of the Hartree Fock energy.
  END_DOC

  logical                        :: has
  PROVIDE ezfio_filename
  if (mpi_master) then
    
    call ezfio_has_dirac_scf_dirac_thresh_scf(has)
    if (has) then
      call ezfio_get_dirac_scf_dirac_thresh_scf(dirac_thresh_scf)
    else
      print *, 'dirac_scf/dirac_thresh_scf not found in EZFIO file'
      stop 1
    endif
  endif
  IRP_IF MPI
    include 'mpif.h'
    integer :: ierr
    call MPI_BCAST( dirac_thresh_scf, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
    if (ierr /= MPI_SUCCESS) then
      stop 'Unable to read dirac_thresh_scf with MPI'
    endif
  IRP_ENDIF

  call write_time(6)
  if (mpi_master) then
    write(6, *) 'Read  dirac_thresh_scf'
  endif

END_PROVIDER

BEGIN_PROVIDER [ double precision, dirac_ao_integrals_threshold  ]
  implicit none
  BEGIN_DOC
! If |<pq|rs>| < dirac_ao_integrals_threshold then <pq|rs> is zero
  END_DOC

  logical                        :: has
  PROVIDE ezfio_filename
  if (mpi_master) then
    
    call ezfio_has_dirac_scf_dirac_threshold_ao(has)
    if (has) then
      call ezfio_get_dirac_scf_dirac_threshold_ao(dirac_ao_integrals_threshold)
    else
      print *, 'dirac_scf/dirac_threshold_ao not found in EZFIO file'
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

BEGIN_PROVIDER [ double precision, speed_of_light  ]
  implicit none
  BEGIN_DOC
! Speed of light
  END_DOC

  logical                        :: has
  PROVIDE ezfio_filename
  if (mpi_master) then
    
    call ezfio_has_dirac_scf_speed_of_light(has)
    if (has) then
      call ezfio_get_dirac_scf_speed_of_light(speed_of_light)
    else
      print *, 'dirac_scf/speed_of_light not found in EZFIO file'
      stop 1
    endif
  endif
  IRP_IF MPI
    include 'mpif.h'
    integer :: ierr
    call MPI_BCAST( speed_of_light, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
    if (ierr /= MPI_SUCCESS) then
      stop 'Unable to read speed_of_light with MPI'
    endif
  IRP_ENDIF

  call write_time(6)
  if (mpi_master) then
    write(6, *) 'Read  speed_of_light'
  endif

END_PROVIDER

BEGIN_PROVIDER [ character*(32), dirac_mo_guess_type  ]
  implicit none
  BEGIN_DOC
! Initial MO guess. Can be [ HCore ]
  END_DOC

  logical                        :: has
  PROVIDE ezfio_filename
  if (mpi_master) then
    
    call ezfio_has_dirac_scf_dirac_mo_guess_type(has)
    if (has) then
      call ezfio_get_dirac_scf_dirac_mo_guess_type(dirac_mo_guess_type)
    else
      print *, 'dirac_scf/dirac_mo_guess_type not found in EZFIO file'
      stop 1
    endif
  endif
  IRP_IF MPI
    include 'mpif.h'
    integer :: ierr
    call MPI_BCAST( dirac_mo_guess_type, 1*32, MPI_CHARACTER, 0, MPI_COMM_WORLD, ierr)
    if (ierr /= MPI_SUCCESS) then
      stop 'Unable to read dirac_mo_guess_type with MPI'
    endif
  IRP_ENDIF

  call write_time(6)
  if (mpi_master) then
    write(6, *) 'Read  dirac_mo_guess_type'
  endif

END_PROVIDER

BEGIN_PROVIDER [ logical, dirac_ao_cartesian  ]
  implicit none
  BEGIN_DOC
! If true, use dirac AOs in Cartesian coordinates (6d,10f,...)
  END_DOC

  logical                        :: has
  PROVIDE ezfio_filename
  if (mpi_master) then
    
    call ezfio_has_dirac_scf_dirac_ao_cartesian(has)
    if (has) then
      call ezfio_get_dirac_scf_dirac_ao_cartesian(dirac_ao_cartesian)
    else
      print *, 'dirac_scf/dirac_ao_cartesian not found in EZFIO file'
      stop 1
    endif
  endif
  IRP_IF MPI
    include 'mpif.h'
    integer :: ierr
    call MPI_BCAST( dirac_ao_cartesian, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
    if (ierr /= MPI_SUCCESS) then
      stop 'Unable to read dirac_ao_cartesian with MPI'
    endif
  IRP_ENDIF

  call write_time(6)
  if (mpi_master) then
    write(6, *) 'Read  dirac_ao_cartesian'
  endif

END_PROVIDER

BEGIN_PROVIDER [ character*(32), dirac_exchange_functional  ]
  implicit none
  BEGIN_DOC
! name of the relativistic exchange functional
  END_DOC

  logical                        :: has
  PROVIDE ezfio_filename
  if (mpi_master) then
    
    call ezfio_has_dirac_scf_dirac_exchange_functional(has)
    if (has) then
      call ezfio_get_dirac_scf_dirac_exchange_functional(dirac_exchange_functional)
    else
      print *, 'dirac_scf/dirac_exchange_functional not found in EZFIO file'
      stop 1
    endif
  endif
  IRP_IF MPI
    include 'mpif.h'
    integer :: ierr
    call MPI_BCAST( dirac_exchange_functional, 1*32, MPI_CHARACTER, 0, MPI_COMM_WORLD, ierr)
    if (ierr /= MPI_SUCCESS) then
      stop 'Unable to read dirac_exchange_functional with MPI'
    endif
  IRP_ENDIF

  call write_time(6)
  if (mpi_master) then
    write(6, *) 'Read  dirac_exchange_functional'
  endif

END_PROVIDER

BEGIN_PROVIDER [ character*(32), dirac_scf_algorithm  ]
  implicit none
  BEGIN_DOC
! Type of SCF algorithm used. Possible choices are [ Simple | DIIS]
  END_DOC

  logical                        :: has
  PROVIDE ezfio_filename
  if (mpi_master) then
    
    call ezfio_has_dirac_scf_dirac_scf_algorithm(has)
    if (has) then
      call ezfio_get_dirac_scf_dirac_scf_algorithm(dirac_scf_algorithm)
    else
      print *, 'dirac_scf/dirac_scf_algorithm not found in EZFIO file'
      stop 1
    endif
  endif
  IRP_IF MPI
    include 'mpif.h'
    integer :: ierr
    call MPI_BCAST( dirac_scf_algorithm, 1*32, MPI_CHARACTER, 0, MPI_COMM_WORLD, ierr)
    if (ierr /= MPI_SUCCESS) then
      stop 'Unable to read dirac_scf_algorithm with MPI'
    endif
  IRP_ENDIF

  call write_time(6)
  if (mpi_master) then
    write(6, *) 'Read  dirac_scf_algorithm'
  endif

END_PROVIDER

BEGIN_PROVIDER [ integer, dirac_n_it_scf_max  ]
  implicit none
  BEGIN_DOC
! Maximum number of SCF iterations
  END_DOC

  logical                        :: has
  PROVIDE ezfio_filename
  if (mpi_master) then
    
    call ezfio_has_dirac_scf_dirac_n_it_scf_max(has)
    if (has) then
      call ezfio_get_dirac_scf_dirac_n_it_scf_max(dirac_n_it_scf_max)
    else
      print *, 'dirac_scf/dirac_n_it_scf_max not found in EZFIO file'
      stop 1
    endif
  endif
  IRP_IF MPI
    include 'mpif.h'
    integer :: ierr
    call MPI_BCAST( dirac_n_it_scf_max, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    if (ierr /= MPI_SUCCESS) then
      stop 'Unable to read dirac_n_it_scf_max with MPI'
    endif
  IRP_ENDIF

  call write_time(6)
  if (mpi_master) then
    write(6, *) 'Read  dirac_n_it_scf_max'
  endif

END_PROVIDER

BEGIN_PROVIDER [ character*(32), dirac_range_separation  ]
  implicit none
  BEGIN_DOC
! Use of full-range interaction or only long-range interaction. Possible choices are [ Full_range | Long_range]
  END_DOC

  logical                        :: has
  PROVIDE ezfio_filename
  if (mpi_master) then
    
    call ezfio_has_dirac_scf_dirac_range_separation(has)
    if (has) then
      call ezfio_get_dirac_scf_dirac_range_separation(dirac_range_separation)
    else
      print *, 'dirac_scf/dirac_range_separation not found in EZFIO file'
      stop 1
    endif
  endif
  IRP_IF MPI
    include 'mpif.h'
    integer :: ierr
    call MPI_BCAST( dirac_range_separation, 1*32, MPI_CHARACTER, 0, MPI_COMM_WORLD, ierr)
    if (ierr /= MPI_SUCCESS) then
      stop 'Unable to read dirac_range_separation with MPI'
    endif
  IRP_ENDIF

  call write_time(6)
  if (mpi_master) then
    write(6, *) 'Read  dirac_range_separation'
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
    
    call ezfio_has_dirac_scf_dirac_interaction(has)
    if (has) then
      call ezfio_get_dirac_scf_dirac_interaction(dirac_interaction)
    else
      print *, 'dirac_scf/dirac_interaction not found in EZFIO file'
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
