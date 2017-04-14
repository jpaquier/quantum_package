program selection_slave
  implicit none
  BEGIN_DOC
! Helper program to compute the PT2 in distributed mode.
  END_DOC

  read_wf = .False.
  SOFT_TOUCH read_wf
  call provide_everything
  call switch_qp_run_to_master
  call run_wf
end

subroutine provide_everything
  PROVIDE H_apply_buffer_allocated mo_bielec_integrals_in_map psi_det_generators psi_coef_generators psi_det_sorted_bit psi_selectors n_det_generators n_states generators_bitmask zmq_context
   PROVIDE pt2_e0_denominator mo_tot_num N_int fragment_count
end

subroutine run_wf
  use f77_zmq
  implicit none

  integer(ZMQ_PTR), external :: new_zmq_to_qp_run_socket
  integer(ZMQ_PTR) :: zmq_to_qp_run_socket
  double precision :: energy(N_states)
  character*(64) :: states(1)
  integer :: rc, i
  
  call provide_everything
  
  zmq_context = f77_zmq_ctx_new ()
  states(1) = 'selection'

  zmq_to_qp_run_socket = new_zmq_to_qp_run_socket()

  do

    call wait_for_states(states,zmq_state,1)

    if(trim(zmq_state) == 'Stopped') then

      exit

    else if (trim(zmq_state) == 'selection') then

      ! Selection
      ! ---------

      print *,  'Selection'
      call zmq_get_psi(zmq_to_qp_run_socket,1,energy,N_states)
  
      !$OMP PARALLEL PRIVATE(i)
      i = omp_get_thread_num()
      call selection_slave_tcp(i, energy)
      !$OMP END PARALLEL
      print *,  'Selection done'

    endif

  end do
end


subroutine selection_slave_tcp(i,energy)
  implicit none
  double precision, intent(in) :: energy(N_states)
  integer, intent(in)            :: i

  call run_selection_slave(0,i,energy)
end
