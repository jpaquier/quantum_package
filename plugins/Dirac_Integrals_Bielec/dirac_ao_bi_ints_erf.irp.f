 use map_module
 double precision function dirac_ao_bielec_integral_erf(i,j,k,l)
  implicit none
  BEGIN_DOC
  !  integral of the dirac AO basis <ik|jl> or (ij|kl)
  !     i(r1) j(r1) 1/r12 k(r2) l(r2)
  END_DOC
  integer,intent(in)             :: i,j,k,l
  integer                        :: p,q,r,s
  double precision               :: I_center(3),J_center(3),K_center(3),L_center(3)
  integer                        :: num_i,num_j,num_k,num_l,dim1,I_power(3),J_power(3),K_power(3),L_power(3)
  double precision               :: integral
  include 'Utils/constants.include.F'
  double precision               :: P_new(0:max_dim,3),P_center(3),fact_p,pp
  double precision               :: Q_new(0:max_dim,3),Q_center(3),fact_q,qq
  integer                        :: iorder_p(3), iorder_q(3)
  double precision               :: dirac_ao_bielec_integral_schwartz_accel_erf
  !First, remove from the start the integrals that have no physical meaning
  if ((i .le. large_ao_num .and. j .le. large_ao_num .and. k .le. large_ao_num .and. l .gt. large_ao_num) .or.  &
      (i .le. large_ao_num .and. j .le. large_ao_num .and. k .gt. large_ao_num .and. l .le. large_ao_num) .or.  &
      (i .le. large_ao_num .and. j .gt. large_ao_num .and. k .le. large_ao_num .and. l .le. large_ao_num) .or.  &
      (i .gt. large_ao_num .and. j .le. large_ao_num .and. k .le. large_ao_num .and. l .le. large_ao_num) .or.  &
      (i .le. large_ao_num .and. j .gt. large_ao_num .and. k .gt. large_ao_num .and. l .gt. large_ao_num) .or.  &
      (i .gt. large_ao_num .and. j .le. large_ao_num .and. k .gt. large_ao_num .and. l .gt. large_ao_num) .or.  &
      (i .gt. large_ao_num .and. j .gt. large_ao_num .and. k .le. large_ao_num .and. l .gt. large_ao_num) .or.  &
      (i .gt. large_ao_num .and. j .gt. large_ao_num .and. k .gt. large_ao_num .and. l .le. large_ao_num)) then
   dirac_ao_bielec_integral_erf = 0.d0
  else 
   if (dirac_ao_prim_num(i) * dirac_ao_prim_num(j) * dirac_ao_prim_num(k) * dirac_ao_prim_num(l) > 1024 ) then
    dirac_ao_bielec_integral_erf = dirac_ao_bielec_integral_schwartz_accel_erf(i,j,k,l)
    return
   endif
   dim1 = n_pt_max_integrals
   num_i = dirac_ao_nucl(i)
   num_j = dirac_ao_nucl(j)
   num_k = dirac_ao_nucl(k)
   num_l = dirac_ao_nucl(l)
   dirac_ao_bielec_integral_erf = 0.d0
   if (num_i /= num_j .or. num_k /= num_l .or. num_j /= num_k)then
    do p = 1, 3
     I_power(p) = dirac_ao_power(i,p)
     J_power(p) = dirac_ao_power(j,p)
     K_power(p) = dirac_ao_power(k,p)
     L_power(p) = dirac_ao_power(l,p)
     I_center(p) = nucl_coord(num_i,p)
     J_center(p) = nucl_coord(num_j,p)
     K_center(p) = nucl_coord(num_k,p)
     L_center(p) = nucl_coord(num_l,p)
    enddo
   double precision               :: coef1, coef2, coef3, coef4
   double precision               :: p_inv,q_inv
   double precision               :: general_primitive_integral_erf
   do p = 1, dirac_ao_prim_num(i)  
    coef1 = dirac_ao_coef_normalized_ordered_transp(p,i)
    do q = 1, dirac_ao_prim_num(j)
     coef2 = coef1*dirac_ao_coef_normalized_ordered_transp(q,j)
     call give_explicit_poly_and_gaussian(P_new,P_center,pp,fact_p,iorder_p,   &
          dirac_ao_expo_ordered_transp(p,i),dirac_ao_expo_ordered_transp(q,j), &
          I_power,J_power,I_center,J_center,dim1)
     p_inv = 1.d0/pp
     do r = 1, dirac_ao_prim_num(k)
      coef3 = coef2*dirac_ao_coef_normalized_ordered_transp(r,k)
      do s = 1, dirac_ao_prim_num(l)
       coef4 = coef3*dirac_ao_coef_normalized_ordered_transp(s,l)
       call give_explicit_poly_and_gaussian(Q_new,Q_center,qq,fact_q,iorder_q,   &
            dirac_ao_expo_ordered_transp(r,k),dirac_ao_expo_ordered_transp(s,l), &
            K_power,L_power,K_center,L_center,dim1)
       q_inv = 1.d0/qq
       integral = general_primitive_integral_erf(dim1,           &
                  P_new,P_center,fact_p,pp,p_inv,iorder_p,   &
                  Q_new,Q_center,fact_q,qq,q_inv,iorder_q)
       dirac_ao_bielec_integral_erf +=+  coef4 * integral
      enddo ! s
     enddo  ! r
    enddo   ! q
   enddo    ! p
   else
    do p = 1, 3
     I_power(p) = dirac_ao_power(i,p)
     J_power(p) = dirac_ao_power(j,p)
     K_power(p) = dirac_ao_power(k,p)
     L_power(p) = dirac_ao_power(l,p)
    enddo
   double  precision              :: ERI_erf
    do p = 1, dirac_ao_prim_num(i)
     coef1 = dirac_ao_coef_normalized_ordered_transp(p,i)
     do q = 1, dirac_ao_prim_num(j)
      coef2 = coef1*dirac_ao_coef_normalized_ordered_transp(q,j)
      do r = 1, dirac_ao_prim_num(k)
       coef3 = coef2*dirac_ao_coef_normalized_ordered_transp(r,k)
       do s = 1, dirac_ao_prim_num(l)
        coef4 = coef3*dirac_ao_coef_normalized_ordered_transp(s,l)
        integral = ERI_erf(                                                 &
                   dirac_ao_expo_ordered_transp(p,i),dirac_ao_expo_ordered_transp(q,j),dirac_ao_expo_ordered_transp(r,k),dirac_ao_expo_ordered_transp(s,l),&
                   I_power(1),J_power(1),K_power(1),L_power(1),         &
                   I_power(2),J_power(2),K_power(2),L_power(2),         &
                   I_power(3),J_power(3),K_power(3),L_power(3))
        dirac_ao_bielec_integral_erf += coef4 * integral
       enddo ! s
      enddo  ! r
     enddo   ! q
    enddo    ! p
   endif
  endif
 end



 double precision function dirac_ao_bielec_integral_schwartz_accel_erf(i,j,k,l)
  implicit none
  BEGIN_DOC
  !  integral of the AO basis <ik|jl> or (ij|kl)
  !     i(r1) j(r1) 1/r12 k(r2) l(r2)
  END_DOC
  integer,intent(in)             :: i,j,k,l
  integer                        :: p,q,r,s
  double precision               :: I_center(3),J_center(3),K_center(3),L_center(3)
  integer                        :: num_i,num_j,num_k,num_l,dim1,I_power(3),J_power(3),K_power(3),L_power(3)
  double precision               :: integral
  include 'Utils/constants.include.F'
  double precision               :: P_new(0:max_dim,3),P_center(3),fact_p,pp
  double precision               :: Q_new(0:max_dim,3),Q_center(3),fact_q,qq
  integer                        :: iorder_p(3), iorder_q(3)
  double precision, allocatable  :: schwartz_kl(:,:)
  double precision               :: schwartz_ij
  dim1 = n_pt_max_integrals
  num_i = dirac_ao_nucl(i)
  num_j = dirac_ao_nucl(j)
  num_k = dirac_ao_nucl(k)
  num_l = dirac_ao_nucl(l)
  dirac_ao_bielec_integral_schwartz_accel_erf = 0.d0
  double precision               :: thr
  thr = ao_integrals_threshold*ao_integrals_threshold
  allocate(schwartz_kl(0:dirac_ao_prim_num(l),0:dirac_ao_prim_num(k)))
  if (num_i /= num_j .or. num_k /= num_l .or. num_j /= num_k)then
   do p = 1, 3
    I_power(p) = dirac_ao_power(i,p)
    J_power(p) = dirac_ao_power(j,p)
    K_power(p) = dirac_ao_power(k,p)
    L_power(p) = dirac_ao_power(l,p)
    I_center(p) = nucl_coord(num_i,p)
    J_center(p) = nucl_coord(num_j,p)
    K_center(p) = nucl_coord(num_k,p)
    L_center(p) = nucl_coord(num_l,p)
   enddo
   schwartz_kl(0,0) = 0.d0
   do r = 1, dirac_ao_prim_num(k)
    coef1 = dirac_ao_coef_normalized_ordered_transp(r,k)*dirac_ao_coef_normalized_ordered_transp(r,k)
    schwartz_kl(0,r) = 0.d0
    do s = 1, dirac_ao_prim_num(l)
     coef2 = coef1 * dirac_ao_coef_normalized_ordered_transp(s,l) * dirac_ao_coef_normalized_ordered_transp(s,l)
     call give_explicit_poly_and_gaussian(Q_new,Q_center,qq,fact_q,iorder_q,&
          dirac_ao_expo_ordered_transp(r,k),dirac_ao_expo_ordered_transp(s,l),                 &
          K_power,L_power,K_center,L_center,dim1)
     q_inv = 1.d0/qq
     schwartz_kl(s,r) = general_primitive_integral_erf(dim1,          &
                        Q_new,Q_center,fact_q,qq,q_inv,iorder_q,                 &
                        Q_new,Q_center,fact_q,qq,q_inv,iorder_q)                 &
                        * coef2 
     schwartz_kl(0,r) = max(schwartz_kl(0,r),schwartz_kl(s,r))
    enddo
    schwartz_kl(0,0) = max(schwartz_kl(0,r),schwartz_kl(0,0))
   enddo
   do p = 1, dirac_ao_prim_num(i)
   double precision               :: coef1
    coef1 = dirac_ao_coef_normalized_ordered_transp(p,i)
    do q = 1, dirac_ao_prim_num(j)
     double precision               :: coef2
     coef2 = coef1*dirac_ao_coef_normalized_ordered_transp(q,j)
     double precision               :: p_inv,q_inv
     call give_explicit_poly_and_gaussian(P_new,P_center,pp,fact_p,iorder_p,&
          dirac_ao_expo_ordered_transp(p,i),dirac_ao_expo_ordered_transp(q,j),                 &
          I_power,J_power,I_center,J_center,dim1)
     p_inv = 1.d0/pp
     schwartz_ij = general_primitive_integral_erf(dim1,               &
                   P_new,P_center,fact_p,pp,p_inv,iorder_p,                 &
                   P_new,P_center,fact_p,pp,p_inv,iorder_p)                 &
                   *coef2*coef2
     if (schwartz_kl(0,0)*schwartz_ij < thr) then
      cycle
     endif
     do r = 1, dirac_ao_prim_num(k)
      if (schwartz_kl(0,r)*schwartz_ij < thr) then
       cycle
      endif
      double precision               :: coef3
      coef3 = coef2*dirac_ao_coef_normalized_ordered_transp(r,k)
      do s = 1, dirac_ao_prim_num(l)
       double precision               :: coef4
       if (schwartz_kl(s,r)*schwartz_ij < thr) then
        cycle
       endif
       coef4 = coef3*dirac_ao_coef_normalized_ordered_transp(s,l)
       double precision               :: general_primitive_integral_erf
       call give_explicit_poly_and_gaussian(Q_new,Q_center,qq,fact_q,iorder_q,&
            dirac_ao_expo_ordered_transp(r,k),dirac_ao_expo_ordered_transp(s,l),             &
            K_power,L_power,K_center,L_center,dim1)
       q_inv = 1.d0/qq
       integral = general_primitive_integral_erf(dim1,              &
                  P_new,P_center,fact_p,pp,p_inv,iorder_p,             &
                  Q_new,Q_center,fact_q,qq,q_inv,iorder_q)
       dirac_ao_bielec_integral_schwartz_accel_erf +=  coef4 * integral
      enddo ! s
     enddo  ! r
    enddo   ! q
   enddo    ! p
  else
   do p = 1, 3
    I_power(p) = dirac_ao_power(i,p)
    J_power(p) = dirac_ao_power(j,p)
    K_power(p) = dirac_ao_power(k,p)
    L_power(p) = dirac_ao_power(l,p)
   enddo
   double  precision              :: ERI_erf
   schwartz_kl(0,0) = 0.d0
   do r = 1, dirac_ao_prim_num(k)
    coef1 = dirac_ao_coef_normalized_ordered_transp(r,k)*dirac_ao_coef_normalized_ordered_transp(r,k)
    schwartz_kl(0,r) = 0.d0
    do s = 1, dirac_ao_prim_num(l)
     coef2 = coef1*dirac_ao_coef_normalized_ordered_transp(s,l)*dirac_ao_coef_normalized_ordered_transp(s,l)
     schwartz_kl(s,r) = ERI_erf(                                      &
                        dirac_ao_expo_ordered_transp(r,k),dirac_ao_expo_ordered_transp(s,l),dirac_ao_expo_ordered_transp(r,k),dirac_ao_expo_ordered_transp(s,l),&
                        K_power(1),L_power(1),K_power(1),L_power(1),             &
                        K_power(2),L_power(2),K_power(2),L_power(2),             &
                        K_power(3),L_power(3),K_power(3),L_power(3))             &
                        *coef2
     schwartz_kl(0,r) = max(schwartz_kl(0,r),schwartz_kl(s,r))
    enddo
    schwartz_kl(0,0) = max(schwartz_kl(0,r),schwartz_kl(0,0))
   enddo
   do p = 1, dirac_ao_prim_num(i)
    coef1 = dirac_ao_coef_normalized_ordered_transp(p,i)
    do q = 1, dirac_ao_prim_num(j)
     coef2 = coef1*dirac_ao_coef_normalized_ordered_transp(q,j)
     schwartz_ij = ERI_erf(                                          &
                   dirac_ao_expo_ordered_transp(p,i),dirac_ao_expo_ordered_transp(q,j),dirac_ao_expo_ordered_transp(p,i),dirac_ao_expo_ordered_transp(q,j),&
                   I_power(1),J_power(1),I_power(1),J_power(1),         &
                   I_power(2),J_power(2),I_power(2),J_power(2),         &
                   I_power(3),J_power(3),I_power(3),J_power(3))*coef2*coef2
     if (schwartz_kl(0,0)*schwartz_ij < thr) then
      cycle
     endif
     do r = 1, dirac_ao_prim_num(k)
      if (schwartz_kl(0,r)*schwartz_ij < thr) then
       cycle
      endif
      coef3 = coef2*dirac_ao_coef_normalized_ordered_transp(r,k)
      do s = 1, dirac_ao_prim_num(l)
       if (schwartz_kl(s,r)*schwartz_ij < thr) then
        cycle
       endif
       coef4 = coef3*dirac_ao_coef_normalized_ordered_transp(s,l)
       integral = ERI_erf(                                          &
                  dirac_ao_expo_ordered_transp(p,i),dirac_ao_expo_ordered_transp(q,j),dirac_ao_expo_ordered_transp(r,k),dirac_ao_expo_ordered_transp(s,l),&
                  I_power(1),J_power(1),K_power(1),L_power(1),         &
                  I_power(2),J_power(2),K_power(2),L_power(2),         &
                  I_power(3),J_power(3),K_power(3),L_power(3))
       dirac_ao_bielec_integral_schwartz_accel_erf += coef4 * integral
      enddo ! s
     enddo  ! r
    enddo   ! q
   enddo    ! p
  endif
  deallocate (schwartz_kl)
 end

 subroutine compute_dirac_ao_bielec_integrals_erf(j,k,l,sze,buffer_value)
  implicit none
  use map_module
  BEGIN_DOC
  ! Compute AO 1/r12 integrals for all i and fixed j,k,l
  END_DOC
  include 'Utils/constants.include.F'
  integer, intent(in)            :: j,k,l,sze
  real(integral_kind), intent(out) :: buffer_value(sze)
  double precision               :: dirac_ao_bielec_integral_erf
  integer                        :: i
  if (dirac_ao_overlap_abs(j,l) < thresh) then
   buffer_value = 0._integral_kind
   return
  endif
  if (dirac_ao_bielec_integral_erf_schwartz(j,l) < thresh ) then
   buffer_value = 0._integral_kind
   return
  endif
  do i = 1, dirac_ao_num
   if (dirac_ao_overlap_abs(i,k)*dirac_ao_overlap_abs(j,l) < thresh) then
    buffer_value(i) = 0._integral_kind
    cycle
   endif
   if (dirac_ao_bielec_integral_erf_schwartz(i,k)*dirac_ao_bielec_integral_erf_schwartz(j,l) < thresh ) then
    buffer_value(i) = 0._integral_kind
    cycle
   endif
 !DIR$ FORCEINLINE
   buffer_value(i) = dirac_ao_bielec_integral_erf(i,k,j,l)
  enddo
 end


 BEGIN_PROVIDER [ logical, dirac_ao_bielec_integrals_erf_in_map ]
  implicit none
  use f77_zmq
  use map_module
  BEGIN_DOC
 ! Map of Atomic integrals
 ! i(r1) j(r2) 1/r12 k(r1) l(r2)
  END_DOC
  integer                        :: i,j,k,l
  double precision               :: dirac_ao_bielec_integral_erf,cpu_1,cpu_2, wall_1, wall_2
  double precision               :: integral, wall_0
  include 'Utils/constants.include.F'
  ! For integrals file
  integer(key_kind),allocatable  :: buffer_i(:)
  integer,parameter              :: size_buffer = 1024*64
  real(integral_kind),allocatable :: buffer_value(:)
  integer                        :: n_integrals, rc
  integer                        :: kk, m, j1, i1, lmax
  character*(64)                 :: fmt
  integral = dirac_ao_bielec_integral_erf(1,1,1,1)
  double precision               :: map_mb
 !PROVIDE read_ao_integrals disk_access_ao_integrals
 !if (read_ao_integrals) then
 ! print*,'Reading the AO integrals'
 ! call map_load_from_disk(trim(ezfio_filename)//'/work/ao_ints',ao_integrals_map)
 ! print*, 'AO integrals provided'
 ! ao_bielec_integrals_in_map = .True.
 ! return
 !endif
  print*, 'Providing the Dirac AO ERF integrals'
  call wall_time(wall_0)
  call wall_time(wall_1)
  call cpu_time(cpu_1)
  integer(ZMQ_PTR) :: zmq_to_qp_run_socket, zmq_socket_pull
  call new_parallel_job(zmq_to_qp_run_socket,zmq_socket_pull,'dirac_ao_integrals_erf')
  character(len=:), allocatable :: task
  allocate(character(len=dirac_ao_num*12) :: task)
  write(fmt,*) '(', dirac_ao_num, '(I5,X,I5,''|''))'
  do l=1,dirac_ao_num
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
        call dirac_ao_bielec_integrals_erf_in_map_collector(zmq_socket_pull)
      else
        call dirac_ao_bielec_integrals_erf_in_map_slave_inproc(i)
      endif
 !$OMP END PARALLEL
  call end_parallel_job(zmq_to_qp_run_socket, zmq_socket_pull, 'dirac_ao_integrals_erf')
  print*, 'Sorting the map'
  call map_sort(dirac_ao_integrals_erf_map)
  call cpu_time(cpu_2)
  call wall_time(wall_2)
  integer(map_size_kind)         :: get_dirac_ao_erf_map_size, dirac_ao_erf_map_size
  dirac_ao_erf_map_size = get_dirac_ao_erf_map_size()
  print*, 'DIRAC AO ERF integrals provided:'
  print*, ' Size of DIRAC AO ERF map :         ', map_mb(dirac_ao_integrals_erf_map) ,'MB'
  print*, ' Number of DIRAC AO integrals  :', dirac_ao_erf_map_size
  print*, ' cpu  time :',cpu_2 - cpu_1, 's'
  print*, ' wall time :',wall_2 - wall_1, 's  ( x ', (cpu_2-cpu_1)/(wall_2-wall_1+tiny(1.d0)), ' )'
  dirac_ao_bielec_integrals_erf_in_map = .True.
 !if (write_dirac_ao_integrals.and.mpi_master) then
 !  call ezfio_set_work_empty(.False.)
 !  call map_save_to_disk(trim(ezfio_filename)//'/work/ao_ints',ao_integrals_map)
 !  call ezfio_set_integrals_bielec_disk_access_ao_integrals("Read")
 !endif
 END_PROVIDER

 BEGIN_PROVIDER [ double precision, dirac_ao_bielec_integral_erf_schwartz,(dirac_ao_num,dirac_ao_num)  ]
  implicit none
  BEGIN_DOC
  !  Needed to compute Schwartz inequalities
  END_DOC
  integer                        :: i,k
  double precision               :: dirac_ao_bielec_integral_erf,cpu_1,cpu_2, wall_1, wall_2
  dirac_ao_bielec_integral_erf_schwartz(1,1) = dirac_ao_bielec_integral_erf(1,1,1,1)
 !$OMP PARALLEL DO PRIVATE(i,k)                                     &
 !$OMP DEFAULT(NONE)                                            &
 !$OMP SHARED (dirac_ao_num,dirac_ao_bielec_integral_erf_schwartz)              &
 !$OMP SCHEDULE(dynamic)
  do i=1,dirac_ao_num
    do k=1,i
      dirac_ao_bielec_integral_erf_schwartz(i,k) = dsqrt(dirac_ao_bielec_integral_erf(i,k,i,k))
      dirac_ao_bielec_integral_erf_schwartz(k,i) = dirac_ao_bielec_integral_erf_schwartz(i,k)
    enddo
  enddo
    !$OMP END PARALLEL DO
 END_PROVIDER

 subroutine dirac_ao_bielec_integrals_erf_in_map_collector(zmq_socket_pull)
  use map_module
  use f77_zmq
  implicit none
  BEGIN_DOC
 !Collects results from the AO integral calculation
  END_DOC
  integer(ZMQ_PTR), intent(in)   :: zmq_socket_pull
  integer                        :: j,l,n_integrals
  integer                        :: rc
  real(integral_kind), allocatable :: buffer_value(:)
  integer(key_kind), allocatable :: buffer_i(:)
  integer(ZMQ_PTR),external      :: new_zmq_to_qp_run_socket
  integer(ZMQ_PTR)               :: zmq_to_qp_run_socket
  integer(ZMQ_PTR), external     :: new_zmq_pull_socket
  integer*8                      :: control, accu, sze
  integer                        :: task_id, more
  zmq_to_qp_run_socket = new_zmq_to_qp_run_socket()
  sze = dirac_ao_num*dirac_ao_num
  allocate ( buffer_i(sze), buffer_value(sze) )
  accu = 0_8
  more = 1
  do while (more == 1)
   rc = f77_zmq_recv( zmq_socket_pull, n_integrals, 4, 0)
   if (rc == -1) then
    n_integrals = 0
    return
   endif
   if (rc /= 4) then
    print *, irp_here,  ': f77_zmq_recv( zmq_socket_pull, n_integrals, 4, 0)'
    stop 'error'
   endif
   if (n_integrals >= 0) then
    if (n_integrals > sze) then
     deallocate (buffer_value, buffer_i)
     sze = n_integrals
     allocate (buffer_value(sze), buffer_i(sze))
    endif
    rc = f77_zmq_recv( zmq_socket_pull, buffer_i, key_kind*n_integrals, 0)
    if (rc /= key_kind*n_integrals) then
     print *,  rc, key_kind, n_integrals
     print *, irp_here,  ': f77_zmq_recv( zmq_socket_pull, buffer_i, key_kind*n_integrals, 0)'
     stop 'error'
    endif
    rc = f77_zmq_recv( zmq_socket_pull, buffer_value, integral_kind*n_integrals, 0)
    if (rc /= integral_kind*n_integrals) then
     print *, irp_here,  ': f77_zmq_recv( zmq_socket_pull, buffer_value, integral_kind*n_integrals, 0)'
     stop 'error'
    endif
    rc = f77_zmq_recv( zmq_socket_pull, task_id, 4, 0)
 IRP_IF ZMQ_PUSH
 IRP_ELSE
    rc = f77_zmq_send( zmq_socket_pull, 0, 4, 0)
    if (rc /= 4) then
     print *,  irp_here, ' : f77_zmq_send (zmq_socket_pull,...'
     stop 'error'
    endif
 IRP_ENDIF
    call insert_into_dirac_ao_integrals_erf_map(n_integrals,buffer_i,buffer_value)
    accu += n_integrals
    if (task_id /= 0) then
     integer, external :: zmq_delete_task
     if (zmq_delete_task(zmq_to_qp_run_socket,zmq_socket_pull,task_id,more) == -1) then
      stop 'Unable to delete task'
     endif
    endif
   endif
  enddo
  deallocate( buffer_i, buffer_value )
  integer (map_size_kind) :: get_dirac_ao_erf_map_size 
  control = get_dirac_ao_erf_map_size(dirac_ao_integrals_erf_map)
  if (control /= accu) then
   print *, ''
   print *, irp_here
   print *, 'Control : ', control
   print *, 'Accu    : ', accu
   print *, 'Some integrals were lost during the parallel computation.'
   print *, 'Try to reduce the number of threads.'
   stop
  endif
  call end_zmq_to_qp_run_socket(zmq_to_qp_run_socket)
 end

 subroutine insert_into_dirac_ao_integrals_erf_map(n_integrals,buffer_i, buffer_values)
  use map_module
  implicit none
  BEGIN_DOC
  ! Create new entry into AO map
  END_DOC
  integer, intent(in)                :: n_integrals
  integer(key_kind), intent(inout)   :: buffer_i(n_integrals)
  real(integral_kind), intent(inout) :: buffer_values(n_integrals)
  call map_append(dirac_ao_integrals_erf_map, buffer_i, buffer_values, n_integrals)
 end

 subroutine dirac_ao_bielec_integrals_erf_in_map_slave_inproc(i)
  implicit none
  integer, intent(in)            :: i
  BEGIN_DOC
 ! Computes a buffer of integrals. i is the ID of the current thread.
  END_DOC
  call dirac_ao_bielec_integrals_erf_in_map_slave(1,i)
 end

 subroutine dirac_ao_bielec_integrals_erf_in_map_slave(thread,iproc)
  use map_module
  use f77_zmq
  implicit none
  BEGIN_DOC
 ! Computes a buffer of integrals
  END_DOC
  integer, intent(in)            :: thread, iproc
  integer                        :: j,l,n_integrals
  integer                        :: rc
  real(integral_kind), allocatable :: buffer_value(:)
  integer(key_kind), allocatable :: buffer_i(:)
  integer                        :: worker_id, task_id
  character*(512)                :: task
  integer(ZMQ_PTR),external      :: new_zmq_to_qp_run_socket
  integer(ZMQ_PTR)               :: zmq_to_qp_run_socket
  integer(ZMQ_PTR), external     :: new_zmq_push_socket
  integer(ZMQ_PTR)               :: zmq_socket_push
  character*(64)                 :: state
  zmq_to_qp_run_socket = new_zmq_to_qp_run_socket()
  integer, external :: connect_to_taskserver
  if (connect_to_taskserver(zmq_to_qp_run_socket,worker_id,thread) == -1) then
   call end_zmq_to_qp_run_socket(zmq_to_qp_run_socket)
   return
  endif
  zmq_socket_push      = new_zmq_push_socket(thread)
  allocate ( buffer_i(dirac_ao_num*dirac_ao_num), buffer_value(dirac_ao_num*dirac_ao_num) )
  do
   integer, external :: get_task_from_taskserver
   if (get_task_from_taskserver(zmq_to_qp_run_socket,worker_id, task_id, task) == -1) then
    exit
   endif
   if (task_id == 0) exit
   read(task,*) j, l
   integer, external :: task_done_to_taskserver
   call compute_dirac_ao_integrals_erf_jl(j,l,n_integrals,buffer_i,buffer_value)
   if (task_done_to_taskserver(zmq_to_qp_run_socket,worker_id,task_id) == -1)  then
    stop 'Unable to send task_done'
   endif
   call push_integrals(zmq_socket_push, n_integrals, buffer_i, buffer_value, task_id)
  enddo
  integer, external :: disconnect_from_taskserver
  if (disconnect_from_taskserver(zmq_to_qp_run_socket,worker_id) == -1) then
   continue
  endif
  deallocate( buffer_i, buffer_value )
  call end_zmq_to_qp_run_socket(zmq_to_qp_run_socket)
  call end_zmq_push_socket(zmq_socket_push,thread)
 end

 subroutine compute_dirac_ao_integrals_erf_jl(j,l,n_integrals,buffer_i,buffer_value)
  implicit none
  use map_module
  BEGIN_DOC
  !  Parallel client for AO integrals
  END_DOC
  integer, intent(in)                   :: j,l
  integer,intent(out)                   :: n_integrals
  integer(key_kind),intent(out)         :: buffer_i(dirac_ao_num*dirac_ao_num)
  real(integral_kind),intent(out)       :: buffer_value(dirac_ao_num*dirac_ao_num)
  integer                               :: i,k
  double precision                      :: dirac_ao_bielec_integral_erf,cpu_1,cpu_2, wall_0,wall_1, wall_2
  double precision                      :: integral
  double precision                      :: thr
  integer                               :: kk, m, j1, i1
  thr = dirac_ao_integrals_threshold
  n_integrals = 0
  j1 = j+ishft(l*l-l,-1)
  do k = 1, dirac_ao_num           ! r1
    i1 = ishft(k*k-k,-1)
    if (i1 > j1) then
      exit
    endif
    do i = 1, k
      i1 += 1
      if (i1 > j1) then
        exit
      endif
      if (dirac_ao_overlap_abs(i,k)*dirac_ao_overlap_abs(j,l) < thr) then
        cycle
      endif
      if (dirac_ao_bielec_integral_erf_schwartz(i,k)*dirac_ao_bielec_integral_erf_schwartz(j,l) < thr ) then
        cycle
      endif
      !DIR$ FORCEINLINE
      integral = dirac_ao_bielec_integral_erf(i,k,j,l)  ! i,k : r1    j,l : r2
      if (abs(integral) < thr) then
        cycle
      endif
      n_integrals += 1
      !DIR$ FORCEINLINE
      call bielec_integrals_index(i,j,k,l,buffer_i(n_integrals))
      buffer_value(n_integrals) = integral
    enddo
  enddo
 end
                              

 
 
 ! AO Map
 ! ======

 BEGIN_PROVIDER [ type(map_type), dirac_ao_integrals_erf_map ]
  implicit none
  BEGIN_DOC
  ! AO integrals
  END_DOC
  integer(key_kind)              :: key_max
  integer(map_size_kind)         :: sze
  call bielec_integrals_index(dirac_ao_num,dirac_ao_num,dirac_ao_num,dirac_ao_num,key_max)
  sze = key_max
  call map_init(dirac_ao_integrals_erf_map,sze)
  print*,  'DIRAC AO map initialized : ', sze
 END_PROVIDER


 BEGIN_PROVIDER [ integer, dirac_ao_integrals_erf_cache_min ]
 &BEGIN_PROVIDER [ integer, dirac_ao_integrals_erf_cache_max ]
  implicit none
  BEGIN_DOC
  ! Min and max values of the AOs for which the integrals are in the cache
  END_DOC
  dirac_ao_integrals_erf_cache_min = max(1,dirac_ao_num - 63)
  dirac_ao_integrals_erf_cache_max = dirac_ao_num
 END_PROVIDER

 BEGIN_PROVIDER [ double precision, dirac_ao_integrals_erf_cache, (0:64*64*64*64) ]
  implicit none
  BEGIN_DOC
  ! Cache of AO integrals for fast access
  END_DOC
  PROVIDE dirac_ao_bielec_integrals_erf_in_map
  integer                        :: i,j,k,l,ii
  integer(key_kind)              :: idx
  real(integral_kind)            :: integral
 !$OMP PARALLEL DO PRIVATE (i,j,k,l,idx,ii,integral)
  do l=dirac_ao_integrals_erf_cache_min,dirac_ao_integrals_erf_cache_max
   do k=dirac_ao_integrals_erf_cache_min,dirac_ao_integrals_erf_cache_max
    do j=dirac_ao_integrals_erf_cache_min,dirac_ao_integrals_erf_cache_max
     do i=dirac_ao_integrals_erf_cache_min,dirac_ao_integrals_erf_cache_max
 !DIR$ FORCEINLINE
      call bielec_integrals_index(i,j,k,l,idx)
 !DIR$ FORCEINLINE
      call map_get(dirac_ao_integrals_erf_map,idx,integral)
      ii = l-dirac_ao_integrals_erf_cache_min
      ii = ior( ishft(ii,6), k-dirac_ao_integrals_erf_cache_min)
      ii = ior( ishft(ii,6), j-dirac_ao_integrals_erf_cache_min)
      ii = ior( ishft(ii,6), i-dirac_ao_integrals_erf_cache_min)
      dirac_ao_integrals_erf_cache(ii) = integral
     enddo
    enddo
   enddo
  enddo
 !$OMP END PARALLEL DO
 END_PROVIDER

double precision function get_dirac_ao_bielec_integral_erf(i,j,k,l,map) result(result)
  use map_module
  implicit none
  BEGIN_DOC
  ! Gets one AO bi-electronic integral from the AO map
  ! i,k :: r1, j,l :: r2
  END_DOC
  integer, intent(in)            :: i,j,k,l
  integer(key_kind)              :: idx
  type(map_type), intent(inout)  :: map
  integer                        :: ii
  real(integral_kind)            :: tmp
  PROVIDE dirac_ao_bielec_integrals_erf_in_map dirac_ao_integrals_erf_cache dirac_ao_integrals_erf_cache_min
  !DIR$ FORCEINLINE
  if (dirac_ao_overlap_abs(i,k)*dirac_ao_overlap_abs(j,l) < dirac_ao_integrals_threshold ) then
    tmp = 0.d0
  else if (dirac_ao_bielec_integral_erf_schwartz(i,k)*dirac_ao_bielec_integral_erf_schwartz(j,l) < dirac_ao_integrals_threshold) then
    tmp = 0.d0
  else
    ii = l-dirac_ao_integrals_erf_cache_min
    ii = ior(ii, k-dirac_ao_integrals_erf_cache_min)
    ii = ior(ii, j-dirac_ao_integrals_erf_cache_min)
    ii = ior(ii, i-dirac_ao_integrals_erf_cache_min)
    if (iand(ii, -64) /= 0) then
      !DIR$ FORCEINLINE
      call bielec_integrals_index(i,j,k,l,idx)
      !DIR$ FORCEINLINE
      call map_get(map,idx,tmp)
    else
      ii = l-dirac_ao_integrals_erf_cache_min
      ii = ior( ishft(ii,6), k-dirac_ao_integrals_erf_cache_min)
      ii = ior( ishft(ii,6), j-dirac_ao_integrals_erf_cache_min)
      ii = ior( ishft(ii,6), i-dirac_ao_integrals_erf_cache_min)
      tmp = dirac_ao_integrals_erf_cache(ii)
    endif
  endif
  result = tmp
end


 subroutine get_dirac_ao_bielec_integrals_erf(j,k,l,sze,out_val)
  use map_module
  BEGIN_DOC
  ! Gets multiple AO bi-electronic integral from the AO map .
  ! All i are retrieved for j,k,l fixed.
  ! j,l :: r1  ;;; l :: r2
  END_DOC
  implicit none
  integer, intent(in)            :: j,k,l, sze
  real(integral_kind), intent(out) :: out_val(sze)
  integer                        :: i
  integer(key_kind)              :: hash
  double precision               :: thresh
  PROVIDE dirac_ao_bielec_integrals_erf_in_map dirac_ao_integrals_erf_map
  thresh = dirac_ao_integrals_threshold
  if (dirac_ao_overlap_abs(j,l) < thresh) then
   out_val = 0.d0
   return
  endif
  double precision :: get_dirac_ao_bielec_integral_erf
  do i=1,sze
   out_val(i) = get_dirac_ao_bielec_integral_erf(i,j,k,l,dirac_ao_integrals_erf_map)
  enddo
 end

 subroutine get_dirac_ao_bielec_integrals_erf_non_zero(j,k,l,sze,out_val,out_val_index,non_zero_int)
  use map_module
  implicit none
  BEGIN_DOC
  ! Gets multiple AO bi-electronic integral from the AO map .
  ! All non-zero i are retrieved for j,k,l fixed.
  END_DOC
  integer, intent(in)            :: j,k,l, sze
  real(integral_kind), intent(out) :: out_val(sze)
  integer, intent(out)           :: out_val_index(sze),non_zero_int
  integer                        :: i
  integer(key_kind)              :: hash
  double precision               :: thresh,tmp
  PROVIDE dirac_ao_bielec_integrals_erf_in_map
  thresh = dirac_ao_integrals_threshold
  non_zero_int = 0
  if (dirac_ao_overlap_abs(j,l) < thresh) then
   out_val = 0.d0
   return
  endif
  non_zero_int = 0
  do i=1,sze
   integer, external :: dirac_ao_l4
   double precision, external :: dirac_ao_bielec_integral
 !DIR$ FORCEINLINE
   if (dirac_ao_bielec_integral_erf_schwartz(i,k)*dirac_ao_bielec_integral_erf_schwartz(j,l) < thresh) then
    cycle
   endif
   call bielec_integrals_index(i,j,k,l,hash)
   call map_get(dirac_ao_integrals_erf_map, hash,tmp)
   if (dabs(tmp) < thresh ) cycle
   non_zero_int = non_zero_int+1
   out_val_index(non_zero_int) = i
   out_val(non_zero_int) = tmp
  enddo
 end


 function get_dirac_ao_erf_map_size()
  implicit none
  integer (map_size_kind) :: get_dirac_ao_erf_map_size
  BEGIN_DOC
  ! Returns the number of elements in the AO map
  END_DOC
  get_dirac_ao_erf_map_size = dirac_ao_integrals_erf_map % n_elements
 end
 
 subroutine clear_dirac_ao_erf_map
  implicit none
  BEGIN_DOC
  ! Frees the memory of the AO map
  END_DOC
  call map_deinit(dirac_ao_integrals_erf_map)
  FREE dirac_ao_integrals_erf_map
 end

