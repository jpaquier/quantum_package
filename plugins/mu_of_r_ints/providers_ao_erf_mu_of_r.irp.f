
BEGIN_PROVIDER [ logical, ao_bielec_integrals_erf_mu_of_r_in_map ]
  implicit none
  use f77_zmq
  use map_module
  BEGIN_DOC
  !  Map of Atomic integrals
  !     i(r1) j(r2) 1/r12 k(r1) l(r2)
  END_DOC
  
  integer                        :: i,j,k,l
  double precision               :: erf_mu_of_r_ao,cpu_1,cpu_2, wall_1, wall_2
  double precision               :: integral, wall_0
  include 'Utils/constants.include.F'
  
  ! For integrals file
  integer(key_kind),allocatable  :: buffer_i(:)
  integer,parameter              :: size_buffer = 1024*64
  real(integral_kind),allocatable :: buffer_value(:)
  
  integer                        :: n_integrals, rc
  integer                        :: kk, m, j1, i1, lmax
  character*(64)                 :: fmt
  
  integral = erf_mu_of_r_ao(1,1,1,1)
  double precision :: integrals_matrix(ao_num,ao_num)
  call give_all_erf_mu_of_r_kl(1,1,integrals_matrix)

  
  double precision               :: map_mb
  PROVIDE read_ao_integrals_mu_of_r disk_access_ao_integrals_mu_of_r
  if (read_ao_integrals_mu_of_r) then
    print*,'Reading the AO ERF mu of r integrals'
      call map_load_from_disk(trim(ezfio_filename)//'/work/ao_ints_erf_mu_of_r',ao_integrals_erf_mu_of_r_map)
      print*, 'AO ERF mu of r integrals provided'
      ao_bielec_integrals_erf_mu_of_r_in_map = .True.
      return
  endif
  
  print*, 'Providing the AO ERF mu of r integrals'
  call wall_time(wall_0)
  call wall_time(wall_1)
  call cpu_time(cpu_1)

  integer(ZMQ_PTR) :: zmq_to_qp_run_socket, zmq_socket_pull
  call new_parallel_job(zmq_to_qp_run_socket,zmq_socket_pull,'ao_integrals_erf_mu_of_r')

  character(len=:), allocatable :: task
  allocate(character(len=ao_num*12) :: task)
  write(fmt,*) '(', ao_num, '(I5,X,I5,''|''))'
  do l=1,ao_num
    write(task,fmt) (i,l, i=1,l)
    integer, external :: add_task_to_taskserver
    if (add_task_to_taskserver(zmq_to_qp_run_socket,trim(task)) == -1) then
      stop 'Unable to add task to server'
    endif
  enddo
  deallocate(task)
  
  integer, external :: zmq_set_running
  if (zmq_set_running(zmq_to_qp_run_socket) == -1) then
    print *,  irp_here, ': Failed in zmq_set_running'
  endif

  PROVIDE nproc
  !$OMP PARALLEL DEFAULT(shared) private(i) num_threads(nproc+1)
      i = omp_get_thread_num()
      if (i==0) then
        call ao_bielec_integrals_erf_mu_of_r_in_map_collector(zmq_socket_pull)
      else
        call ao_bielec_integrals_erf_mu_of_r_in_map_slave_inproc(i)
      endif
  !$OMP END PARALLEL

  call end_parallel_job(zmq_to_qp_run_socket, zmq_socket_pull, 'ao_integrals_erf_mu_of_r')


  print*, 'Sorting the map'
  call map_sort(ao_integrals_erf_mu_of_r_map)
  call cpu_time(cpu_2)
  call wall_time(wall_2)
  integer(map_size_kind)         :: get_ao_erf_mu_of_r_map_size, ao_erf_mu_of_r_map_size
  ao_erf_mu_of_r_map_size = get_ao_erf_mu_of_r_map_size()
  
  print*, 'AO ERF mu of r integrals provided:'
  print*, ' Size of AO ERF mu of r map :         ', map_mb(ao_integrals_erf_mu_of_r_map) ,'MB'
  print*, ' Number of AO ERF mu of r integrals :', ao_erf_mu_of_r_map_size
  print*, ' cpu  time :',cpu_2 - cpu_1, 's'
  print*, ' wall time :',wall_2 - wall_1, 's  ( x ', (cpu_2-cpu_1)/(wall_2-wall_1+tiny(1.d0)), ' )'
  
  ao_bielec_integrals_erf_mu_of_r_in_map = .True.

  if (write_ao_integrals_mu_of_r) then
    call ezfio_set_work_empty(.False.)
    call map_save_to_disk(trim(ezfio_filename)//'/work/ao_ints_erf_mu_of_r',ao_integrals_erf_mu_of_r_map)
    call ezfio_set_mu_of_r_ints_disk_access_ao_integrals_mu_of_r("Read")
  endif
  
END_PROVIDER
 

BEGIN_PROVIDER [ logical, ao_bielec_integrals_erf_mu_of_r_in_map_left ]
  implicit none
  use f77_zmq
  use map_module
  BEGIN_DOC
  !  Map of Atomic integrals
  !     i(r1) j(r2) 1/r12 k(r1) l(r2)
  END_DOC
  
  integer                        :: i,j,k,l
  double precision               :: erf_mu_of_r_ao,cpu_1,cpu_2, wall_1, wall_2
  double precision               :: integral, wall_0
  include 'Utils/constants.include.F'
  
  ! For integrals file
  integer(key_kind),allocatable  :: buffer_i(:)
  integer,parameter              :: size_buffer = 1024*64
  real(integral_kind),allocatable :: buffer_value(:)
  
  integer                        :: n_integrals, rc
  integer                        :: kk, m, j1, i1, lmax
  character*(64)                 :: fmt
  
  integral = erf_mu_of_r_ao(1,1,1,1)
  double precision :: integrals_matrix(ao_num,ao_num)
  call give_all_erf_mu_of_r_kl(1,1,integrals_matrix)

  
  double precision               :: map_mb
  PROVIDE read_ao_integrals_mu_of_r disk_access_ao_integrals_mu_of_r
  if (read_ao_integrals_mu_of_r_left) then
    print*,'Reading the AO ERF mu of r integrals'
      call map_load_from_disk(trim(ezfio_filename)//'/work/ao_ints_erf_mu_of_r_left',ao_integrals_erf_mu_of_r_map_left)
      print*, 'AO ERF mu of r integrals left provided'
      ao_bielec_integrals_erf_mu_of_r_in_map_left = .True.
      return
  endif
  
  print*, 'Providing the AO ERF mu of r integrals left'
  call wall_time(wall_0)
  call wall_time(wall_1)
  call cpu_time(cpu_1)

  integer(ZMQ_PTR) :: zmq_to_qp_run_socket, zmq_socket_pull
  call new_parallel_job(zmq_to_qp_run_socket,zmq_socket_pull,'ao_integrals_erf_mu_of_r_left')

  character(len=:), allocatable :: task
  allocate(character(len=ao_num*12) :: task)
  write(fmt,*) '(', ao_num, '(I5,X,I5,''|''))'
  do l=1,ao_num
    write(task,fmt) (i,l, i=1,l)
    integer, external :: add_task_to_taskserver
    if (add_task_to_taskserver(zmq_to_qp_run_socket,trim(task)) == -1) then
      stop 'Unable to add task to server'
    endif
  enddo
  deallocate(task)
  
  integer, external :: zmq_set_running
  if (zmq_set_running(zmq_to_qp_run_socket) == -1) then
    print *,  irp_here, ': Failed in zmq_set_running'
  endif

  PROVIDE nproc
  !$OMP PARALLEL DEFAULT(shared) private(i) num_threads(nproc+1)
      i = omp_get_thread_num()
      if (i==0) then
        call ao_bielec_integrals_erf_mu_of_r_in_map_collector_left(zmq_socket_pull)
      else
        call ao_bielec_integrals_erf_mu_of_r_in_map_slave_inproc_left(i)
      endif
  !$OMP END PARALLEL

  call end_parallel_job(zmq_to_qp_run_socket, zmq_socket_pull, 'ao_integrals_erf_mu_of_r_left')


  print*, 'Sorting the map'
  call map_sort(ao_integrals_erf_mu_of_r_map)
  call cpu_time(cpu_2)
  call wall_time(wall_2)
  integer(map_size_kind)         :: get_ao_erf_mu_of_r_map_size_left, ao_erf_mu_of_r_map_size_left
  ao_erf_mu_of_r_map_size_left = get_ao_erf_mu_of_r_map_size_left()
  
  print*, 'AO ERF mu of r integrals left provided:'
  print*, ' Size of AO ERF mu of r map :         ', map_mb(ao_integrals_erf_mu_of_r_map_left) ,'MB'
  print*, ' Number of AO ERF mu of r integrals :', ao_erf_mu_of_r_map_size_left
  print*, ' cpu  time :',cpu_2 - cpu_1, 's'
  print*, ' wall time :',wall_2 - wall_1, 's  ( x ', (cpu_2-cpu_1)/(wall_2-wall_1+tiny(1.d0)), ' )'
  
  ao_bielec_integrals_erf_mu_of_r_in_map_left = .True.

  if (write_ao_integrals_mu_of_r) then
    call ezfio_set_work_empty(.False.)
    call map_save_to_disk(trim(ezfio_filename)//'/work/ao_ints_erf_mu_of_r_left',ao_integrals_erf_mu_of_r_map_left)
    call ezfio_set_mu_of_r_ints_disk_access_ao_integrals_mu_of_r_left("Read")
  endif
  
END_PROVIDER
 

