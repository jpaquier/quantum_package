
  use bitmasks
BEGIN_PROVIDER [integer(bit_kind), psi_active, (N_int,2,psi_det_size)]
 BEGIN_DOC
! active part of psi
 END_DOC
 implicit none
  use bitmasks
 integer :: i,j,k,l
 provide cas_bitmask
 do i = 1, N_det
  do j = 1, N_int
   psi_active(j,1,i) = iand(psi_det(j,1,i),cas_bitmask(j,1,1))
   psi_active(j,2,i) = iand(psi_det(j,2,i),cas_bitmask(j,1,1))
  enddo
 enddo
END_PROVIDER


subroutine give_holes_and_particles_in_active_space(det_1,det_2,n_holes_spin,n_particles_spin,n_holes,n_particles,& 
                                                    holes_active_list,particles_active_list)
 implicit none
 use bitmasks
 integer(bit_kind),intent(in)  :: det_1(N_int,2)
 integer(bit_kind),intent(in ) :: det_2(N_int,2)
 integer, intent(out) :: n_holes_spin(2),n_particles_spin(2)
 integer, intent(out) :: n_holes,n_particles
 integer, intent(out) :: holes_active_list(2 * n_act_orb,2)
 integer, intent(out) :: particles_active_list(2 * n_act_orb,2)
 integer :: i
 integer(bit_kind) :: holes(N_int,2)
 integer(bit_kind) :: particles(N_int,2)
 integer(bit_kind) :: det_tmp_2(N_int,2),det_tmp_1(N_int,2)
 BEGIN_DOC
! returns the holes and particles operators WITHIN THE ACTIVE SPACE 
! that connect det_1 and det_2. By definition, the holes/particles 
! are such that one starts from det_1 and goes to det_2
!        
! n_holes is the total number of holes
! n_particles is the total number of particles
! n_holes_spin is the number of number of holes per spin (1=alpha, 2=beta)
! n_particles_spin is the number of number of particles per spin (1=alpha, 2=beta)
! holes_active_list is the index of the holes per spin, that ranges from 1 to n_act_orb
! particles_active_list is the index of the particles per spin, that ranges from 1 to n_act_orb
 END_DOC
 
 call give_active_part_determinant(det_1,det_tmp_1)
 call give_active_part_determinant(det_2,det_tmp_2)
 do i = 1, N_int
  holes(i,1) = iand(det_tmp_1(i,1),xor(det_tmp_1(i,1),det_tmp_2(i,1)))
  holes(i,2) = iand(det_tmp_1(i,2),xor(det_tmp_1(i,2),det_tmp_2(i,2)))
  particles(i,1) = iand(det_tmp_2(i,1),xor(det_tmp_1(i,1),det_tmp_2(i,1)))
  particles(i,2) = iand(det_tmp_2(i,2),xor(det_tmp_1(i,2),det_tmp_2(i,2)))
 enddo

 integer :: holes_list(N_int*bit_kind_size,2)
 holes_list = 0
 call bitstring_to_list(holes(1,1), holes_list(1,1), n_holes_spin(1), N_int)
 call bitstring_to_list(holes(1,2), holes_list(1,2), n_holes_spin(2), N_int)

 n_holes = 0
 do i = 1, n_holes_spin(1)
  n_holes +=1
  holes_active_list(i,1) = list_act_reverse(holes_list(i,1))
 enddo
 do i = 1, n_holes_spin(2)
  n_holes +=1
  holes_active_list(i,2) = list_act_reverse(holes_list(i,2))
 enddo


 integer :: particles_list(N_int*bit_kind_size,2)
 particles_list = 0
 call bitstring_to_list(particles(1,1), particles_list(1,1), n_particles_spin(1), N_int)
 call bitstring_to_list(particles(1,2), particles_list(1,2), n_particles_spin(2), N_int)
 n_particles = 0
 do i = 1, n_particles_spin(1)
  n_particles += 1
  particles_active_list(i,1) = list_act_reverse(particles_list(i,1))
 enddo
 do i = 1, n_particles_spin(2)
  n_particles += 1
  particles_active_list(i,2) = list_act_reverse(particles_list(i,2))
 enddo

end

subroutine give_holes_in_inactive_space(det_1,n_holes_spin,n_holes,holes_list)
 BEGIN_DOC
! returns the holes operators WITHIN THE INACTIVE SPACE 
! that has lead to det_1.  
!        
! n_holes is the total number of holes
! n_holes_spin is the number of number of holes per spin (1=alpha, 2=beta)
! holes_inactive_list is the index of the holes per spin, that ranges from 1 to mo_tot_num
 END_DOC
 implicit none
 use bitmasks
 integer(bit_kind),intent(in)  :: det_1(N_int,2)
 integer, intent(out) :: n_holes_spin(2)
 integer, intent(out) :: n_holes
 integer, intent(out) :: holes_list(N_int*bit_kind_size,2)
 integer :: i
 integer(bit_kind) :: holes(N_int,2)
 integer(bit_kind) :: det_tmp_1(N_int,2)
 
 call give_core_inactive_part_determinant(det_1,det_tmp_1)
 
 do i = 1, N_int
  holes(i,1) = iand(reunion_of_core_inact_bitmask(i,1),xor(det_tmp_1(i,1),reunion_of_core_inact_bitmask(i,1)))
  holes(i,2) = iand(reunion_of_core_inact_bitmask(i,2),xor(det_tmp_1(i,2),reunion_of_core_inact_bitmask(i,2)))
 enddo
 holes_list = 0
 call bitstring_to_list(holes(1,1), holes_list(1,1), n_holes_spin(1), N_int)
 call bitstring_to_list(holes(1,2), holes_list(1,2), n_holes_spin(2), N_int)
 n_holes =  n_holes_spin(1) +  n_holes_spin(2)

end

subroutine give_particles_in_virt_space(det_1,n_particles_spin,n_particles,particles_list)
 BEGIN_DOC
! returns the holes operators WITHIN THE VIRTUAL SPACE 
! that has lead to det_1.  
!        
! n_particles is the total number of particles
! n_particles_spin is the number of number of particles per spin (1=alpha, 2=beta)
! particles_inactive_list is the index of the particles per spin, that ranges from 1 to mo_tot_num
 END_DOC
 implicit none
 use bitmasks
 integer(bit_kind),intent(in)  :: det_1(N_int,2)
 integer, intent(out) :: n_particles_spin(2)
 integer, intent(out) :: n_particles
 integer, intent(out) :: particles_list(N_int*bit_kind_size,2)
 integer :: i
 integer(bit_kind) :: det_tmp_1(N_int,2)
 integer(bit_kind) :: particles(N_int,2)

 call give_virt_part_determinant(det_1,det_tmp_1)

 do i = 1, N_int
  particles(i,1) = iand(virt_bitmask(i,1),det_tmp_1(i,1))
  particles(i,2) = iand(virt_bitmask(i,2),det_tmp_1(i,2))
 enddo

 particles_list = 0
 call bitstring_to_list(particles(1,1), particles_list(1,1), n_particles_spin(1), N_int)
 call bitstring_to_list(particles(1,2), particles_list(1,2), n_particles_spin(2), N_int)
 n_particles = n_particles_spin(1) + n_particles_spin(2)
 

end

subroutine get_delta_e_dyall(det_1,det_2,delta_e_final)
 implicit none
  use bitmasks
 double precision, intent(out) :: delta_e_final
 integer(bit_kind), intent(in) :: det_1(N_int,2),det_2(N_int,2)
 integer :: i,j,k,l
 
 integer :: n_holes_spin(2)
 integer :: n_holes
 integer :: holes_list(N_int*bit_kind_size,2)


 double precision :: delta_e_inactive
 integer :: i_hole_inact

 call give_holes_in_inactive_space(det_2,n_holes_spin,n_holes,holes_list)
 delta_e_inactive = 0.d0
 do i = 1, n_holes_spin(1)
  i_hole_inact = holes_list(i,1)
  delta_e_inactive += fock_core_inactive_total_spin_trace(i_hole_inact)
 enddo

 do i = 1, n_holes_spin(2)
  i_hole_inact = holes_list(i,2)
  delta_e_inactive += fock_core_inactive_total_spin_trace(i_hole_inact)
 enddo

 double precision :: delta_e_virt
 integer :: i_part_virt
 integer :: n_particles_spin(2)
 integer :: n_particles
 integer :: particles_list(N_int*bit_kind_size,2)
 
 call give_particles_in_virt_space(det_2,n_particles_spin,n_particles,particles_list)
 delta_e_virt = 0.d0
 do i = 1, n_particles_spin(1)
  i_part_virt = particles_list(i,1)
  delta_e_virt += fock_virt_total_spin_trace(i_part_virt)
 enddo

 do i = 1, n_particles_spin(2)
  i_part_virt = particles_list(i,2)
  delta_e_virt += fock_virt_total_spin_trace(i_part_virt)
 enddo


 integer :: n_holes_spin_act(2),n_particles_spin_act(2)
 integer :: n_holes_act,n_particles_act
 integer :: holes_active_list(2*n_act_orb,2)
 integer :: holes_active_list_spin_traced(4*n_act_orb)
 integer :: particles_active_list(2*n_act_orb,2)
 integer :: particles_active_list_spin_traced(4*n_act_orb)
 double precision :: delta_e_act
 delta_e_act = 0.d0
 call give_holes_and_particles_in_active_space(det_1,det_2,n_holes_spin_act,n_particles_spin_act, &
                                               n_holes_act,n_particles_act,holes_active_list,particles_active_list)
 integer :: icount,icountbis
 integer :: hole_list_practical(2,elec_num_tab(1)+elec_num_tab(2)), particle_list_practical(2,elec_num_tab(1)+elec_num_tab(2))
 icount = 0
 icountbis = 0
 do i = 1, n_holes_spin_act(1)
  icount += 1
  icountbis += 1
  hole_list_practical(1,icountbis) = 1
  hole_list_practical(2,icountbis) =  holes_active_list(i,1)
  holes_active_list_spin_traced(icount) = holes_active_list(i,1)
 enddo
 do i = 1, n_holes_spin_act(2)
  icount += 1
  icountbis += 1
  hole_list_practical(1,icountbis) = 2
  hole_list_practical(2,icountbis) =  holes_active_list(i,2)
  holes_active_list_spin_traced(icount) = holes_active_list(i,2)
 enddo
 if(icount .ne. n_holes_act) then
  print*,''
  print*, icount, n_holes_act
  print * , 'pb in holes_active_list_spin_traced !!' 
  stop
 endif

 icount = 0
 icountbis = 0
 do i = 1, n_particles_spin_act(1)
  icount += 1
  icountbis += 1
  particle_list_practical(1,icountbis) = 1
  particle_list_practical(2,icountbis) =  particles_active_list(i,1)
  particles_active_list_spin_traced(icount) = particles_active_list(i,1)
 enddo
 do i = 1, n_particles_spin_act(2)
  icount += 1
  icountbis += 1
  particle_list_practical(1,icountbis) = 2
  particle_list_practical(2,icountbis) =  particles_active_list(i,2)
  particles_active_list_spin_traced(icount) = particles_active_list(i,2)
 enddo
 if(icount .ne. n_particles_act) then
  print*, icount, n_particles_act
  print * , 'pb in particles_active_list_spin_traced !!' 
  stop
 endif


 integer :: i_hole_act, j_hole_act, k_hole_act
 integer :: i_particle_act, j_particle_act, k_particle_act
 

 integer :: ispin,jspin,kspin
 if      (n_holes_act == 0 .and. n_particles_act == 1) then
! i_particle_act =  particles_active_list_spin_traced(1)
! delta_e_act += one_creat_spin_trace(i_particle_act )
  ispin = particle_list_practical(1,1)
  i_particle_act =  particle_list_practical(2,1)
  delta_e_act += one_creat(i_particle_act,ispin)

 else if (n_holes_act == 1 .and. n_particles_act == 0) then
! i_hole_act =  holes_active_list_spin_traced(1)
! delta_e_act += one_anhil_spin_trace(i_hole_act )
  ispin = hole_list_practical(1,1)
  i_hole_act =  hole_list_practical(2,1)
  delta_e_act += one_anhil(i_hole_act , ispin)

 else if (n_holes_act == 1 .and. n_particles_act == 1) then
! i_hole_act =  holes_active_list_spin_traced(1)
! i_particle_act =  particles_active_list_spin_traced(1)
! delta_e_act += one_anhil_one_creat_spin_trace(i_hole_act,i_particle_act)
  ! first hole
  ispin = hole_list_practical(1,1)
  i_hole_act =  hole_list_practical(2,1)
  ! first particle
  jspin = particle_list_practical(1,1)
  i_particle_act =  particle_list_practical(2,1)
  delta_e_act += one_anhil_one_creat(i_particle_act,i_hole_act,jspin,ispin)

 else if (n_holes_act == 2 .and. n_particles_act == 0) then
! i_hole_act =  holes_active_list_spin_traced(1)
! j_hole_act =  holes_active_list_spin_traced(1)
! delta_e_act += two_anhil_spin_trace(i_hole_act,j_hole_act)
  ispin = hole_list_practical(1,1)
  i_hole_act =  hole_list_practical(2,1)
  jspin = hole_list_practical(1,2)
  j_hole_act =  hole_list_practical(2,2)
  delta_e_act += two_anhil(i_hole_act,j_hole_act,ispin,jspin)

 else if (n_holes_act == 0 .and. n_particles_act == 2) then
! i_particle_act =  particles_active_list_spin_traced(1)
! j_particle_act =  particles_active_list_spin_traced(2)
! delta_e_act += two_creat_spin_trace(i_particle_act,j_particle_act)
  ispin = particle_list_practical(1,1)
  i_particle_act =  particle_list_practical(2,1)
  jspin = particle_list_practical(1,2)
  j_particle_act =  particle_list_practical(2,2)
  delta_e_act += two_creat(i_particle_act,j_particle_act,ispin,jspin)

 else if (n_holes_act == 2 .and. n_particles_act == 1) then
! i_hole_act =  holes_active_list_spin_traced(1)
! j_hole_act =  holes_active_list_spin_traced(2)
! i_particle_act =  particles_active_list_spin_traced(1)
! print*, 'i_hole_act,j_hole_act,i_particle_act'
! print*, i_hole_act,j_hole_act,i_particle_act
! print*, two_anhil_one_creat_spin_trace(i_hole_act,j_hole_act,i_particle_act)
! delta_e_act += two_anhil_one_creat_spin_trace(i_hole_act,j_hole_act,i_particle_act)

  ! first hole
  ispin = hole_list_practical(1,1)
  i_hole_act =  hole_list_practical(2,1)
  ! second hole
  jspin = hole_list_practical(1,2)
  j_hole_act =  hole_list_practical(2,2)
  ! first particle
  kspin = particle_list_practical(1,1)
  i_particle_act =  particle_list_practical(2,1)
  delta_e_act += two_anhil_one_creat(i_particle_act,i_hole_act,j_hole_act,kspin,ispin,jspin)

 else if (n_holes_act == 1 .and. n_particles_act == 2) then
! i_hole_act =  holes_active_list_spin_traced(1)
! i_particle_act =  particles_active_list_spin_traced(1)
! j_particle_act =  particles_active_list_spin_traced(2)
! delta_e_act += two_creat_one_anhil_spin_trace(i_hole_act,i_particle_act,j_particle_act)

  ! first hole
  ispin = hole_list_practical(1,1)
  i_hole_act =  hole_list_practical(2,1)
  ! first particle
  jspin = particle_list_practical(1,1)
  i_particle_act =  particle_list_practical(2,1)
  ! second particle
  kspin = particle_list_practical(1,2)
  j_particle_act =  particle_list_practical(2,2)

  delta_e_act += two_creat_one_anhil(i_particle_act,j_particle_act,i_hole_act,jspin,kspin,ispin)

 else if (n_holes_act == 3 .and. n_particles_act == 0) then
! i_hole_act =  holes_active_list_spin_traced(1)
! j_hole_act =  holes_active_list_spin_traced(2)
! k_hole_act =  holes_active_list_spin_traced(3)
! delta_e_act += three_anhil_spin_trace(i_hole_act,j_hole_act,k_hole_act)

  ! first hole
  ispin = hole_list_practical(1,1)
  i_hole_act =  hole_list_practical(2,1)
  ! second hole
  jspin = hole_list_practical(1,2)
  j_hole_act =  hole_list_practical(2,2)
  ! third hole
  kspin = hole_list_practical(1,3)
  k_hole_act =  hole_list_practical(2,3)
  delta_e_act += three_anhil(i_hole_act,j_hole_act,k_hole_act,ispin,jspin,kspin)

 else if (n_holes_act == 0 .and. n_particles_act == 3) then
! i_particle_act =  particles_active_list_spin_traced(1)
! j_particle_act =  particles_active_list_spin_traced(2)
! k_particle_act =  particles_active_list_spin_traced(3)
! delta_e_act += three_creat_spin_trace(i_particle_act,j_particle_act,k_particle_act)
  ! first particle
  ispin = particle_list_practical(1,1)
  i_particle_act =  particle_list_practical(2,1)
  ! second particle
  jspin = particle_list_practical(1,2)
  j_particle_act =  particle_list_practical(2,2)
  ! second particle
  kspin = particle_list_practical(1,3)
  k_particle_act =  particle_list_practical(2,3)
 
  delta_e_act += three_creat(i_particle_act,j_particle_act,k_particle_act,ispin,jspin,kspin)

 else if (n_holes_act .ge. 2 .and. n_particles_act .ge.2) then

  delta_e_act = -10000000.d0

 endif

!print*, 'one_anhil_spin_trace'
!print*,  one_anhil_spin_trace(1), one_anhil_spin_trace(2) 


 delta_e_final = delta_e_act  + delta_e_inactive - delta_e_virt
!if(delta_e_final .le. -100d0.or.delta_e_final > 0.d0 .or. delta_e_final == 0.d0)then
!if(delta_e_final == 0.d0)then
 if(.False.)then
 call debug_det(det_1,N_int)
 call debug_det(det_2,N_int)
 print*, 'n_holes_act,n_particles_act'
 print*,  n_holes_act,n_particles_act 
 print*, 'delta_e_act,delta_e_inactive,delta_e_vir'
 print*,  delta_e_act,delta_e_inactive,delta_e_virt 
 delta_e_final = -1000.d0
!stop
  
 endif

end

subroutine get_delta_e_dyall_verbose(det_1,det_2,delta_e_final)
 implicit none
  use bitmasks
 double precision, intent(out) :: delta_e_final
 integer(bit_kind), intent(in) :: det_1(N_int,2),det_2(N_int,2)
 integer :: i,j,k,l
 
 integer :: n_holes_spin(2)
 integer :: n_holes
 integer :: holes_list(N_int*bit_kind_size,2)


 double precision :: delta_e_inactive
 integer :: i_hole_inact

 call give_holes_in_inactive_space(det_2,n_holes_spin,n_holes,holes_list)
 delta_e_inactive = 0.d0
 do i = 1, n_holes_spin(1)
  i_hole_inact = holes_list(i,1)
  delta_e_inactive += fock_core_inactive_total_spin_trace(i_hole_inact)
 enddo

 do i = 1, n_holes_spin(2)
  i_hole_inact = holes_list(i,2)
  delta_e_inactive += fock_core_inactive_total_spin_trace(i_hole_inact)
 enddo

 double precision :: delta_e_virt
 integer :: i_part_virt
 integer :: n_particles_spin(2)
 integer :: n_particles
 integer :: particles_list(N_int*bit_kind_size,2)
 
 call give_particles_in_virt_space(det_2,n_particles_spin,n_particles,particles_list)
 delta_e_virt = 0.d0
 do i = 1, n_particles_spin(1)
  i_part_virt = particles_list(i,1)
  delta_e_virt += fock_virt_total_spin_trace(i_part_virt)
 enddo

 do i = 1, n_particles_spin(2)
  i_part_virt = particles_list(i,2)
  delta_e_virt += fock_virt_total_spin_trace(i_part_virt)
 enddo


 integer :: n_holes_spin_act(2),n_particles_spin_act(2)
 integer :: n_holes_act,n_particles_act
 integer :: holes_active_list(2*n_act_orb,2)
 integer :: holes_active_list_spin_traced(4*n_act_orb)
 integer :: particles_active_list(2*n_act_orb,2)
 integer :: particles_active_list_spin_traced(4*n_act_orb)
 double precision :: delta_e_act
 delta_e_act = 0.d0
 call give_holes_and_particles_in_active_space(det_1,det_2,n_holes_spin_act,n_particles_spin_act, &
                                               n_holes_act,n_particles_act,holes_active_list,particles_active_list)
 integer :: icount,icountbis
 integer :: hole_list_practical(2,elec_num_tab(1)+elec_num_tab(2)), particle_list_practical(2,elec_num_tab(1)+elec_num_tab(2))
 icount = 0
 icountbis = 0
 do i = 1, n_holes_spin_act(1)
  icount += 1
  icountbis += 1
  hole_list_practical(1,icountbis) = 1
  hole_list_practical(2,icountbis) =  holes_active_list(i,1)
  holes_active_list_spin_traced(icount) = holes_active_list(i,1)
 enddo
 do i = 1, n_holes_spin_act(2)
  icount += 1
  icountbis += 1
  hole_list_practical(1,icountbis) = 2
  hole_list_practical(2,icountbis) =  holes_active_list(i,2)
  holes_active_list_spin_traced(icount) = holes_active_list(i,2)
 enddo
 if(icount .ne. n_holes_act) then
  print*,''
  print*, icount, n_holes_act
  print * , 'pb in holes_active_list_spin_traced !!' 
  stop
 endif

 icount = 0
 icountbis = 0
 do i = 1, n_particles_spin_act(1)
  icount += 1
  icountbis += 1
  particle_list_practical(1,icountbis) = 1
  particle_list_practical(2,icountbis) =  particles_active_list(i,1)
  particles_active_list_spin_traced(icount) = particles_active_list(i,1)
 enddo
 do i = 1, n_particles_spin_act(2)
  icount += 1
  icountbis += 1
  particle_list_practical(1,icountbis) = 2
  particle_list_practical(2,icountbis) =  particles_active_list(i,2)
  particles_active_list_spin_traced(icount) = particles_active_list(i,2)
 enddo
 if(icount .ne. n_particles_act) then
  print*, icount, n_particles_act
  print * , 'pb in particles_active_list_spin_traced !!' 
  stop
 endif


 integer :: i_hole_act, j_hole_act, k_hole_act
 integer :: i_particle_act, j_particle_act, k_particle_act
 

 integer :: ispin,jspin,kspin
 if      (n_holes_act == 0 .and. n_particles_act == 1) then
! i_particle_act =  particles_active_list_spin_traced(1)
! delta_e_act += one_creat_spin_trace(i_particle_act )
  ispin = particle_list_practical(1,1)
  i_particle_act =  particle_list_practical(2,1)
  delta_e_act += one_creat(i_particle_act,ispin)

 else if (n_holes_act == 1 .and. n_particles_act == 0) then
! i_hole_act =  holes_active_list_spin_traced(1)
! delta_e_act += one_anhil_spin_trace(i_hole_act )
  ispin = hole_list_practical(1,1)
  i_hole_act =  hole_list_practical(2,1)
  delta_e_act += one_anhil(i_hole_act , ispin)

 else if (n_holes_act == 1 .and. n_particles_act == 1) then
! i_hole_act =  holes_active_list_spin_traced(1)
! i_particle_act =  particles_active_list_spin_traced(1)
! delta_e_act += one_anhil_one_creat_spin_trace(i_hole_act,i_particle_act)
  ! first hole
  ispin = hole_list_practical(1,1)
  i_hole_act =  hole_list_practical(2,1)
  ! first particle
  jspin = particle_list_practical(1,1)
  i_particle_act =  particle_list_practical(2,1)
  delta_e_act += one_anhil_one_creat(i_particle_act,i_hole_act,jspin,ispin)

 else if (n_holes_act == 2 .and. n_particles_act == 0) then
! i_hole_act =  holes_active_list_spin_traced(1)
! j_hole_act =  holes_active_list_spin_traced(1)
! delta_e_act += two_anhil_spin_trace(i_hole_act,j_hole_act)
  ispin = hole_list_practical(1,1)
  i_hole_act =  hole_list_practical(2,1)
  jspin = hole_list_practical(1,2)
  j_hole_act =  hole_list_practical(2,2)
  delta_e_act += two_anhil(i_hole_act,j_hole_act,ispin,jspin)

 else if (n_holes_act == 0 .and. n_particles_act == 2) then
! i_particle_act =  particles_active_list_spin_traced(1)
! j_particle_act =  particles_active_list_spin_traced(2)
! delta_e_act += two_creat_spin_trace(i_particle_act,j_particle_act)
  ispin = particle_list_practical(1,1)
  i_particle_act =  particle_list_practical(2,1)
  jspin = particle_list_practical(1,2)
  j_particle_act =  particle_list_practical(2,2)
  delta_e_act += two_creat(i_particle_act,j_particle_act,ispin,jspin)

 else if (n_holes_act == 2 .and. n_particles_act == 1) then
! i_hole_act =  holes_active_list_spin_traced(1)
! j_hole_act =  holes_active_list_spin_traced(2)
! i_particle_act =  particles_active_list_spin_traced(1)
! print*, 'i_hole_act,j_hole_act,i_particle_act'
! print*, i_hole_act,j_hole_act,i_particle_act
! print*, two_anhil_one_creat_spin_trace(i_hole_act,j_hole_act,i_particle_act)
! delta_e_act += two_anhil_one_creat_spin_trace(i_hole_act,j_hole_act,i_particle_act)

  ! first hole
  ispin = hole_list_practical(1,1)
  i_hole_act =  hole_list_practical(2,1)
  ! second hole
  jspin = hole_list_practical(1,2)
  j_hole_act =  hole_list_practical(2,2)
  ! first particle
  kspin = particle_list_practical(1,1)
  i_particle_act =  particle_list_practical(2,1)
  delta_e_act += two_anhil_one_creat(i_particle_act,i_hole_act,j_hole_act,kspin,ispin,jspin)

 else if (n_holes_act == 1 .and. n_particles_act == 2) then
! i_hole_act =  holes_active_list_spin_traced(1)
! i_particle_act =  particles_active_list_spin_traced(1)
! j_particle_act =  particles_active_list_spin_traced(2)
! delta_e_act += two_creat_one_anhil_spin_trace(i_hole_act,i_particle_act,j_particle_act)

  ! first hole
  ispin = hole_list_practical(1,1)
  i_hole_act =  hole_list_practical(2,1)
  ! first particle
  jspin = particle_list_practical(1,1)
  i_particle_act =  particle_list_practical(2,1)
  ! second particle
  kspin = particle_list_practical(1,2)
  j_particle_act =  particle_list_practical(2,2)

  delta_e_act += two_creat_one_anhil(i_particle_act,j_particle_act,i_hole_act,jspin,kspin,ispin)

 else if (n_holes_act == 3 .and. n_particles_act == 0) then
! i_hole_act =  holes_active_list_spin_traced(1)
! j_hole_act =  holes_active_list_spin_traced(2)
! k_hole_act =  holes_active_list_spin_traced(3)
! delta_e_act += three_anhil_spin_trace(i_hole_act,j_hole_act,k_hole_act)

  ! first hole
  ispin = hole_list_practical(1,1)
  i_hole_act =  hole_list_practical(2,1)
  ! second hole
  jspin = hole_list_practical(1,2)
  j_hole_act =  hole_list_practical(2,2)
  ! third hole
  kspin = hole_list_practical(1,3)
  k_hole_act =  hole_list_practical(2,3)
  delta_e_act += three_anhil(i_hole_act,j_hole_act,k_hole_act,ispin,jspin,kspin)

 else if (n_holes_act == 0 .and. n_particles_act == 3) then
! i_particle_act =  particles_active_list_spin_traced(1)
! j_particle_act =  particles_active_list_spin_traced(2)
! k_particle_act =  particles_active_list_spin_traced(3)
! delta_e_act += three_creat_spin_trace(i_particle_act,j_particle_act,k_particle_act)
  ! first particle
  ispin = particle_list_practical(1,1)
  i_particle_act =  particle_list_practical(2,1)
  ! second particle
  jspin = particle_list_practical(1,2)
  j_particle_act =  particle_list_practical(2,2)
  ! second particle
  kspin = particle_list_practical(1,3)
  k_particle_act =  particle_list_practical(2,3)
 
  delta_e_act += three_creat(i_particle_act,j_particle_act,k_particle_act,ispin,jspin,kspin)

 endif

!print*, 'one_anhil_spin_trace'
!print*,  one_anhil_spin_trace(1), one_anhil_spin_trace(2) 


 delta_e_final = delta_e_act  + delta_e_inactive - delta_e_virt
!if(delta_e_final .le. -100d0.or.delta_e_final > 0.d0 .or. delta_e_final == 0.d0)then
!if(delta_e_final == 0.d0)then
 call debug_det(det_1,N_int)
 call debug_det(det_2,N_int)
 print*, 'n_holes_act,n_particles_act'
 print*,  n_holes_act,n_particles_act 
 print*, 'delta_e_act,delta_e_inactive,delta_e_vir'
 print*,  delta_e_act,delta_e_inactive,delta_e_virt 
 delta_e_final = -1000.d0

end