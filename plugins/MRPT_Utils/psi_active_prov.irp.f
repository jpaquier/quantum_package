
  use bitmasks

subroutine give_holes_and_particles_inactive_active_virtual_space(det_ref,det_pert,list_holes_inactive,list_holes_active,list_particle_active,list_particle_virt,n_holes_inactive,n_holes_active,n_particle_active,n_particle_virt)
 implicit none
 use bitmasks
 integer(bit_kind),intent(in)  :: det_ref(N_int,2)
 integer(bit_kind),intent(in ) :: det_pert(N_int,2)
 integer, intent(out) :: n_holes_inactive(2)
 integer, intent(out) :: n_holes_active(2)
 integer, intent(out) :: n_particle_active(2)
 integer, intent(out) :: n_particle_virt(2)
 integer, intent(out) :: list_holes_inactive(N_int*bit_kind_size,2)
 integer, intent(out) :: list_holes_active(N_int*bit_kind_size,2)
 integer, intent(out) :: list_particle_active(N_int*bit_kind_size,2)
 integer, intent(out) :: list_particle_virt(N_int*bit_kind_size,2)
 integer :: list_tmp(N_int*bit_kind_size,2)

 integer :: exc(0:2,2,2)
 integer :: degree
 double precision :: phase
 integer  :: h1,h2,p1,p2,s1,s2
 integer(bit_kind) :: key_holes(N_int,2), key_part(N_int,2)
 integer(bit_kind) :: key_holes_inact(N_int,2)
 integer(bit_kind) :: key_holes_act(N_int,2), key_part_act(N_int,2)
 integer(bit_kind) :: key_part_virt(N_int,2)
 list_holes_inactive = 0
 list_holes_active = 0
 list_particle_active = 0
 list_particle_virt = 0
 
 call get_excitation(det_ref,det_pert,exc,degree,phase,N_int)
 call decode_exc(exc,degree,h1,p1,h2,p2,s1,s2)
 key_holes = 0_bit_kind
 key_part  = 0_bit_kind
 if(degree==2)then
  if(s1==s2)then
   if(s1==1)then
    call set_bit_to_integer(h1,key_holes(1,1),N_int)
    call set_bit_to_integer(h2,key_holes(1,1),N_int)
    call set_bit_to_integer(p1,key_part(1,1),N_int)
    call set_bit_to_integer(p2,key_part(1,1),N_int)
   else
    call set_bit_to_integer(h1,key_holes(1,2),N_int)
    call set_bit_to_integer(h2,key_holes(1,2),N_int)
    call set_bit_to_integer(p1,key_part(1,2),N_int)
    call set_bit_to_integer(p2,key_part(1,2),N_int)
   endif
  else 
    call set_bit_to_integer(h1,key_holes(1,1),N_int)
    call set_bit_to_integer(p1,key_part(1,1),N_int)
    call set_bit_to_integer(h2,key_holes(1,2),N_int)
    call set_bit_to_integer(p2,key_part(1,2),N_int)
  endif
 else 
  if (s1==1)then
   call set_bit_to_integer(h1,key_holes(1,1),N_int)
   call set_bit_to_integer(p1,key_part(1,1),N_int)
  else
   call set_bit_to_integer(h1,key_holes(1,2),N_int)
   call set_bit_to_integer(p1,key_part(1,2),N_int)
  endif
 endif
 integer :: k
 do k = 1, N_int
  key_holes_inact(k,1) = iand(key_holes(k,1),inact_bitmask(k,1))
  key_holes_inact(k,2) = iand(key_holes(k,2),inact_bitmask(k,2))
  key_holes_act(k,1) = iand(key_holes(k,1),act_bitmask(k,1))
  key_holes_act(k,2) = iand(key_holes(k,2),act_bitmask(k,2))
  key_part_act(k,1) = iand(key_part(k,1),act_bitmask(k,1))
  key_part_act(k,2) = iand(key_part(k,2),act_bitmask(k,2))
  key_part_virt(k,1) = iand(key_part(k,1),virt_bitmask(k,1))
  key_part_virt(k,2) = iand(key_part(k,2),virt_bitmask(k,2))
 enddo
 call bitstring_to_list(key_holes_inact(1,1), list_holes_inactive(1,1), n_holes_inactive(1), N_int)
 call bitstring_to_list(key_holes_inact(1,2), list_holes_inactive(1,2), n_holes_inactive(2), N_int)
 
 call bitstring_to_list(key_holes_act(1,1), list_tmp(1,1), n_holes_active(1), N_int)
 do k = 1, n_holes_active(1)
  list_holes_active(k,1) = list_act_reverse(list_tmp(k,1))
 enddo
 call bitstring_to_list(key_holes_act(1,2), list_tmp(1,2), n_holes_active(2), N_int)
 do k = 1, n_holes_active(2)
  list_holes_active(k,2) = list_act_reverse(list_tmp(k,2))
 enddo

 call bitstring_to_list(key_part_act(1,1), list_tmp(1,1), n_particle_active(1), N_int)
 do k = 1, n_particle_active(1)
  list_particle_active(k,1) = list_act_reverse(list_tmp(k,1))
 enddo
 call bitstring_to_list(key_part_act(1,2), list_tmp(1,2), n_particle_active(2), N_int)
 do k = 1, n_particle_active(2)
  list_particle_active(k,2) = list_act_reverse(list_tmp(k,2))
 enddo

 call bitstring_to_list(key_part_virt(1,1), list_particle_virt(1,1), n_particle_virt(1), N_int)
 call bitstring_to_list(key_part_virt(1,2), list_particle_virt(1,2), n_particle_virt(2), N_int)
end


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
 BEGIN_DOC
 ! routine that returns the delta_e with the Moller Plesset and Dyall operators
 !
 ! with det_1 being a determinant from the cas, and det_2 being a perturber
 !
 ! Delta_e(det_1,det_2) = sum (hole) epsilon(hole) + sum(part) espilon(part) + delta_e(act)
 !
 ! where hole is necessary in the inactive, part necessary in the virtuals
 ! 
 ! and delta_e(act) is obtained from the contracted application of the excitation 
 !
 ! operator in the active space that lead from det_1 to det_2
 END_DOC
 implicit none
 use bitmasks
 double precision, intent(out) :: delta_e_final(N_states)
 integer(bit_kind), intent(in) :: det_1(N_int,2),det_2(N_int,2)
 integer :: i,j,k,l
 integer :: i_state
 
 integer :: n_holes_inactive(2)
 integer :: n_holes_active(2)
 integer :: n_particle_active(2)
 integer :: n_particle_virt(2)
 integer :: list_holes_inactive(N_int*bit_kind_size,2)
 integer :: list_holes_active(N_int*bit_kind_size,2)
 integer :: list_particle_active(N_int*bit_kind_size,2)
 integer :: list_particle_virt(N_int*bit_kind_size,2)


 double precision :: delta_e_inactive(N_states)
 integer :: i_hole_inact

!print*,'det print'
!call debug_det(det_2,N_int)


 call get_excitation_degree(det_1,det_2,degree,N_int)
 if(degree>2)then
  delta_e_final = -1.d+10
  return
 endif

 
 call give_holes_and_particles_inactive_active_virtual_space(det_1,det_2,list_holes_inactive,list_holes_active,list_particle_active,list_particle_virt,n_holes_inactive,n_holes_active,n_particle_active,n_particle_virt)

 delta_e_inactive = 0.d0
 do i = 1, n_holes_inactive(1)
  i_hole_inact = list_holes_inactive(i,1)
  do i_state = 1, N_states
   delta_e_inactive += fock_core_inactive_total_spin_trace(i_hole_inact,i_state)
  enddo
 enddo

 do i = 1, n_holes_inactive(2)
  i_hole_inact = list_holes_inactive(i,2)
  do i_state = 1, N_states
   delta_e_inactive(i_state) += fock_core_inactive_total_spin_trace(i_hole_inact,i_state)
  enddo
 enddo

 double precision :: delta_e_virt(N_states)
 integer :: i_part_virt
 
 delta_e_virt = 0.d0
 do i = 1, n_particle_virt(1)
  i_part_virt = list_particle_virt(i,1)
  do i_state = 1, N_states
   delta_e_virt += fock_virt_total_spin_trace(i_part_virt,i_state)
  enddo
 enddo

 do i = 1, n_particle_virt(2)
  i_part_virt = list_particle_virt(i,2)
  do i_state = 1, N_states
   delta_e_virt += fock_virt_total_spin_trace(i_part_virt,i_state)
  enddo
 enddo

 integer :: give_spin_exc_1, spin_exc_2(2),hp_2(2), spin_exc_3(3),hp_3(3)
 double precision :: delta_e_act(N_states)

 delta_e_act = 0.d0
 integer :: n_holes_act,ispin,jspin,kspin
 integer :: i_particle_act,i_hole_act
 integer :: j_particle_act,j_hole_act
 integer :: k_particle_act,k_hole_act
 integer :: n_particles_act 
 n_holes_act = n_holes_active(1) + n_holes_active(2)
 n_particles_act = n_particle_active(1) + n_particle_active(2)
 if      (n_holes_act == 0 .and. n_particles_act == 1) then
  ispin = give_spin_exc_1(n_particle_active)
  i_particle_act =  list_particle_active(1,ispin)
  do i_state = 1, N_states
   delta_e_act(i_state) += one_creat(i_particle_act,ispin,i_state)
  enddo

 else if (n_holes_act == 1 .and. n_particles_act == 0) then
  ispin = give_spin_exc_1(n_holes_active)
  i_hole_act =  list_holes_active(1,ispin)
  do i_state = 1, N_states
   delta_e_act(i_state) += one_anhil(i_hole_act , ispin,i_state)
  enddo

 else if (n_holes_act == 1 .and. n_particles_act == 1) then
  ! first hole
  ispin = give_spin_exc_1(n_holes_active)
  i_hole_act =  list_holes_active(1,ispin)
  ! first particle
  jspin = give_spin_exc_1(n_particle_active)
  i_particle_act =  list_particle_active(1,jspin)
! print*, 'ispin,jspin',ispin,jspin
! call debug_det(det_1,N_int)
! call debug_det(det_2,N_int)
  do i_state = 1, N_states
   delta_e_act(i_state) += one_anhil_one_creat(i_particle_act,i_hole_act,jspin,ispin,i_state)
  enddo

 else if (n_holes_act == 2 .and. n_particles_act == 0) then
  call give_spin_exc_2(n_holes_active,list_holes_active,spin_exc_2,hp_2) 
  ispin = spin_exc_2(1)
  i_hole_act =  hp_2(1)
  jspin = spin_exc_2(2)
  j_hole_act =  hp_2(2)
  do i_state = 1, N_states
   delta_e_act(i_state) += two_anhil(i_hole_act,j_hole_act,ispin,jspin,i_state)
  enddo

 else if (n_holes_act == 0 .and. n_particles_act == 2) then
  call give_spin_exc_2(n_particle_active,list_particle_active,spin_exc_2,hp_2) 
  ! first particle 
  ispin = spin_exc_2(1)
  i_particle_act = hp_2(1)
  ! second particle
  jspin = spin_exc_2(2)
  j_particle_act = hp_2(2)
  do i_state = 1, N_states
   delta_e_act(i_state) += two_creat(i_particle_act,j_particle_act,ispin,jspin,i_state)
  enddo

 else if (n_holes_act == 2 .and. n_particles_act == 1) then
  call give_spin_exc_2(n_holes_active,list_holes_active,spin_exc_2,hp_2) 
  ! first hole
  ispin = spin_exc_2(1)
  i_hole_act =  hp_2(1)
  ! second hole
  jspin = spin_exc_2(2)
  j_hole_act =  hp_2(2)
  ! first particle
  kspin = give_spin_exc_1(n_particle_active)
  i_particle_act =  list_particle_active(1,kspin)
  do i_state = 1, N_states
   delta_e_act(i_state) += two_anhil_one_creat(i_particle_act,i_hole_act,j_hole_act,kspin,ispin,jspin,i_state)
  enddo

 else if (n_holes_act == 1 .and. n_particles_act == 2) then
  ! first hole
  ispin = give_spin_exc_1(n_holes_active)
  i_hole_act =  list_holes_active(1,ispin)
  call give_spin_exc_2(n_particle_active,list_particle_active,spin_exc_2,hp_2) 
  ! first particle 
  jspin = spin_exc_2(1)
  i_particle_act = hp_2(1)
  ! second particle
  kspin = spin_exc_2(2)
  j_particle_act = hp_2(2)

  do i_state = 1, N_states
   delta_e_act(i_state) += two_creat_one_anhil(i_particle_act,j_particle_act,i_hole_act,jspin,kspin,ispin,i_state)
  enddo

 else if (n_holes_act == 3 .and. n_particles_act == 0) then
  call give_spin_exc_3(n_particle_active,list_particle_active,spin_exc_3,hp_3)
  ! first hole
  ispin = spin_exc_3(1)
  i_hole_act =  hp_3(1)
  ! second hole
  jspin = spin_exc_3(2)
  j_hole_act = hp_3(2) 
  ! third hole
  kspin = spin_exc_3(3)
  k_hole_act = hp_3(3) 
  do i_state = 1, N_states
   delta_e_act(i_state) += three_anhil(i_hole_act,j_hole_act,k_hole_act,ispin,jspin,kspin,i_state)
  enddo

 else if (n_holes_act == 0 .and. n_particles_act == 3) then
  call give_spin_exc_3(n_holes_active,list_holes_active,spin_exc_3,hp_3)
  ! first particle
  ispin = spin_exc_3(1)
  i_particle_act = hp_3(1) 
  ! second particle
  jspin = spin_exc_3(2)
  j_particle_act =  hp_3(2) 
  ! second particle
  kspin = spin_exc_3(3)
  k_particle_act =  hp_3(3) 
  do i_state = 1, N_states
   delta_e_act(i_state) += three_creat(i_particle_act,j_particle_act,k_particle_act,ispin,jspin,kspin,i_state)
  enddo
 
 else if (n_holes_act .eq. 0 .and. n_particles_act .eq.0)then
  integer :: degree
  integer(bit_kind) :: det_1_active(N_int,2)
  integer          :: h1,h2,p1,p2,s1,s2
  integer          :: exc(0:2,2,2)
  integer          :: i_hole, i_part
  double precision :: phase
  call get_excitation_degree(det_1,det_2,degree,N_int)
  if(degree == 1)then
   call get_excitation(det_1,det_2,exc,degree,phase,N_int)
   call decode_exc(exc,degree,h1,p1,h2,p2,s1,s2)
!  print*, 'h1 print',h1
!  print*, 'p1 print',p1
   i_hole =  list_inact_reverse(h1)
   i_part =  list_virt_reverse(p1)
   do i_state = 1, N_states
     delta_e_act(i_state) += one_anhil_one_creat_inact_virt(i_hole,i_part,i_state) 
   enddo
  endif
 else if (n_holes_act .ge. 2 .and. n_particles_act .ge.2) then
  delta_e_act = -1.d+20
 endif

!print*, 'one_anhil_spin_trace'
!print*,  one_anhil_spin_trace(1), one_anhil_spin_trace(2) 


 do i_state = 1, n_states
  delta_e_final(i_state) = delta_e_act(i_state)  + delta_e_inactive(i_state) - delta_e_virt(i_state)
 enddo
 logical :: test,test2
 test = (n_holes_act == 1 .and. n_particles_act == 1)
 !print*,n_holes_active
 !print*,n_particle_active
 ispin = give_spin_exc_1(n_holes_active) 
 jspin = give_spin_exc_1(n_particle_active)
 test2 = (ispin == jspin)
!i_particle_act =  list_particle_active(1,jspin)
!i_hole_act =  list_holes_active(1,ispin)
!print*,test,test2
!if(.not.test)then
!  do i_state = 1, n_states
!   delta_e_final(i_state) = 1.d+20
!  enddo
!endif

!if(test)then
! if(test2.eqv..False.)then
!  do i_state = 1, n_states
!   delta_e_final(i_state) = 1.d+20
!  enddo
! else 
! endif
!endif

 !if(n_holes_inactive(2)==1)then
 ! do i_state = 1, n_states
 !  delta_e_final(i_state) = 1.d+20
 ! enddo
 !endif
 !do i_state = 1, N_states
 ! if(delta_e_final(i_state).lt.1.d+19)then
 !  print*,'i_hole_act,i_particle_act',list_act(i_hole_act),list_act(i_particle_act),ispin,jspin
 ! endif
 !enddo


!write(*,'(100(f16.10,X))'), delta_e_final(1) , delta_e_act(1)  , delta_e_inactive(1) , delta_e_virt(1)

end


integer function give_spin_exc_1(array)
 implicit none
 use bitmasks
 integer, intent(in) :: array(2)
 if (array(1).ne.0)then
  give_spin_exc_1 = 1
 else 
  give_spin_exc_1 = 2
 endif
end

subroutine give_spin_exc_2(array_n,array_hp,spin_exc_2,hp_2)
 implicit none
 use bitmasks
 integer, intent(in) :: array_hp(N_int*bit_kind_size,2)
 integer, intent(in) :: array_n(2)
 integer, intent(out) :: spin_exc_2(2),hp_2(2)
 if (array_n(1).eq.2)then
  spin_exc_2(1)= 1
  hp_2(1) = array_hp(1,1)
  spin_exc_2(2)= 1
  hp_2(2) = array_hp(2,1)
 else if(array_n(2) ==2)then
  spin_exc_2(1)= 2
  hp_2(1) = array_hp(1,2)
  spin_exc_2(2)= 2
  hp_2(2) = array_hp(2,2)
 else 
  spin_exc_2(1)= 1
  hp_2(1) = array_hp(1,1)
  spin_exc_2(2)= 2
  hp_2(2) = array_hp(1,2)
 endif
end

subroutine give_spin_exc_3(array_n,array_hp,spin_exc_3,hp_3)
 implicit none
 integer, intent(in) :: array_hp(N_int*bit_kind_size,2)
 integer, intent(in) :: array_n(2)
 integer, intent(out) :: spin_exc_3(3),hp_3(3)
 if (array_n(1).eq.3)then
  spin_exc_3(1)= 1
  hp_3(1) = array_hp(1,1)
  spin_exc_3(2)= 1
  hp_3(2) = array_hp(2,1)
  spin_exc_3(3)= 1
  hp_3(3) = array_hp(3,1)
 else if (array_n(2).eq.3)then
  spin_exc_3(1)= 2
  hp_3(1) = array_hp(1,2)
  spin_exc_3(2)= 2
  hp_3(2) = array_hp(2,2)
  spin_exc_3(3)= 2
  hp_3(3) = array_hp(3,2)
 else if (array_n(1).eq.1)then
  spin_exc_3(1)= 1
  hp_3(1) = array_hp(1,1)
  spin_exc_3(2)= 2
  hp_3(2) = array_hp(2,2)
  spin_exc_3(3)= 2
  hp_3(3) = array_hp(3,2)
 else if (array_n(1).eq.2)then
  spin_exc_3(1)= 1
  hp_3(1) = array_hp(1,1)
  spin_exc_3(2)= 1
  hp_3(2) = array_hp(2,1)
  spin_exc_3(3)= 2
  hp_3(3) = array_hp(3,2)
 endif


end

subroutine get_delta_e_dyall_old(det_1,det_2,delta_e_final)
 BEGIN_DOC
 ! routine that returns the delta_e with the Moller Plesset and Dyall operators
 !
 ! with det_1 being a determinant from the cas, and det_2 being a perturber
 !
 ! Delta_e(det_1,det_2) = sum (hole) epsilon(hole) + sum(part) espilon(part) + delta_e(act)
 !
 ! where hole is necessary in the inactive, part necessary in the virtuals
 ! 
 ! and delta_e(act) is obtained from the contracted application of the excitation 
 !
 ! operator in the active space that lead from det_1 to det_2
 END_DOC
 implicit none
  use bitmasks
 double precision, intent(out) :: delta_e_final(N_states)
 integer(bit_kind), intent(in) :: det_1(N_int,2),det_2(N_int,2)
 integer :: i,j,k,l
 integer :: i_state
 
 integer :: n_holes_spin(2)
 integer :: n_holes
 integer :: holes_list(N_int*bit_kind_size,2)


 double precision :: delta_e_inactive(N_states)
 integer :: i_hole_inact

 call get_excitation_degree(det_1,det_2,degree,N_int)
 if(degree>2)then
  delta_e_final = -1.d+10
  return
 endif

 call give_holes_in_inactive_space(det_2,n_holes_spin,n_holes,holes_list)
 delta_e_inactive = 0.d0
 do i = 1, n_holes_spin(1)
  i_hole_inact = holes_list(i,1)
  do i_state = 1, N_states
   delta_e_inactive += fock_core_inactive_total_spin_trace(i_hole_inact,i_state)
  enddo
 enddo

 do i = 1, n_holes_spin(2)
  i_hole_inact = holes_list(i,2)
  do i_state = 1, N_states
   delta_e_inactive(i_state) += fock_core_inactive_total_spin_trace(i_hole_inact,i_state)
  enddo
 enddo

 double precision :: delta_e_virt(N_states)
 integer :: i_part_virt
 integer :: n_particles_spin(2)
 integer :: n_particles
 integer :: particles_list(N_int*bit_kind_size,2)
 
 call give_particles_in_virt_space(det_2,n_particles_spin,n_particles,particles_list)
 delta_e_virt = 0.d0
 do i = 1, n_particles_spin(1)
  i_part_virt = particles_list(i,1)
  do i_state = 1, N_states
   delta_e_virt += fock_virt_total_spin_trace(i_part_virt,i_state)
  enddo
 enddo

 do i = 1, n_particles_spin(2)
  i_part_virt = particles_list(i,2)
  do i_state = 1, N_states
   delta_e_virt += fock_virt_total_spin_trace(i_part_virt,i_state)
  enddo
 enddo


 integer :: n_holes_spin_act(2),n_particles_spin_act(2)
 integer :: n_holes_act,n_particles_act
 integer :: holes_active_list(2*n_act_orb,2)
 integer :: holes_active_list_spin_traced(4*n_act_orb)
 integer :: particles_active_list(2*n_act_orb,2)
 integer :: particles_active_list_spin_traced(4*n_act_orb)
 double precision :: delta_e_act(N_states)
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
  ispin = particle_list_practical(1,1)
  i_particle_act =  particle_list_practical(2,1)
   do i_state = 1, N_states
    delta_e_act(i_state) += one_creat(i_particle_act,ispin,i_state)
   enddo

 else if (n_holes_act == 1 .and. n_particles_act == 0) then
  ispin = hole_list_practical(1,1)
  i_hole_act =  hole_list_practical(2,1)
  print*, 'i_hole_act print = ',i_hole_act
   do i_state = 1, N_states
    delta_e_act(i_state) += one_anhil(i_hole_act , ispin,i_state)
   enddo

 else if (n_holes_act == 1 .and. n_particles_act == 1) then
  ! first hole
  ispin = hole_list_practical(1,1)
  i_hole_act =  hole_list_practical(2,1)
  ! first particle
  jspin = particle_list_practical(1,1)
  i_particle_act =  particle_list_practical(2,1)
  do i_state = 1, N_states
   delta_e_act(i_state) += one_anhil_one_creat(i_particle_act,i_hole_act,jspin,ispin,i_state)
  enddo

 else if (n_holes_act == 2 .and. n_particles_act == 0) then
  ispin = hole_list_practical(1,1)
  i_hole_act =  hole_list_practical(2,1)
  jspin = hole_list_practical(1,2)
  j_hole_act =  hole_list_practical(2,2)
  do i_state = 1, N_states
   delta_e_act(i_state) += two_anhil(i_hole_act,j_hole_act,ispin,jspin,i_state)
  enddo

 else if (n_holes_act == 0 .and. n_particles_act == 2) then
  ispin = particle_list_practical(1,1)
  i_particle_act =  particle_list_practical(2,1)
  jspin = particle_list_practical(1,2)
  j_particle_act =  particle_list_practical(2,2)
  do i_state = 1, N_states
   delta_e_act(i_state) += two_creat(i_particle_act,j_particle_act,ispin,jspin,i_state)
  enddo

 else if (n_holes_act == 2 .and. n_particles_act == 1) then
  ! first hole
  ispin = hole_list_practical(1,1)
  i_hole_act =  hole_list_practical(2,1)
  ! second hole
  jspin = hole_list_practical(1,2)
  j_hole_act =  hole_list_practical(2,2)
  ! first particle
  kspin = particle_list_practical(1,1)
  i_particle_act =  particle_list_practical(2,1)
  do i_state = 1, N_states
   delta_e_act(i_state) += two_anhil_one_creat(i_particle_act,i_hole_act,j_hole_act,kspin,ispin,jspin,i_state)
  enddo

 else if (n_holes_act == 1 .and. n_particles_act == 2) then
  ! first hole
  ispin = hole_list_practical(1,1)
  i_hole_act =  hole_list_practical(2,1)
  ! first particle
  jspin = particle_list_practical(1,1)
  i_particle_act =  particle_list_practical(2,1)
  ! second particle
  kspin = particle_list_practical(1,2)
  j_particle_act =  particle_list_practical(2,2)

  do i_state = 1, N_states
   delta_e_act(i_state) += two_creat_one_anhil(i_particle_act,j_particle_act,i_hole_act,jspin,kspin,ispin,i_state)
  enddo

 else if (n_holes_act == 3 .and. n_particles_act == 0) then
  ! first hole
  ispin = hole_list_practical(1,1)
  i_hole_act =  hole_list_practical(2,1)
  ! second hole
  jspin = hole_list_practical(1,2)
  j_hole_act =  hole_list_practical(2,2)
  ! third hole
  kspin = hole_list_practical(1,3)
  k_hole_act =  hole_list_practical(2,3)
  do i_state = 1, N_states
   delta_e_act(i_state) += three_anhil(i_hole_act,j_hole_act,k_hole_act,ispin,jspin,kspin,i_state)
  enddo

 else if (n_holes_act == 0 .and. n_particles_act == 3) then
  ! first particle
  ispin = particle_list_practical(1,1)
  i_particle_act =  particle_list_practical(2,1)
  ! second particle
  jspin = particle_list_practical(1,2)
  j_particle_act =  particle_list_practical(2,2)
  ! second particle
  kspin = particle_list_practical(1,3)
  k_particle_act =  particle_list_practical(2,3)
  do i_state = 1, N_states
   delta_e_act(i_state) += three_creat(i_particle_act,j_particle_act,k_particle_act,ispin,jspin,kspin,i_state)
  enddo
 
 else if (n_holes_act .eq. 0 .and. n_particles_act .eq.0)then
  integer :: degree
  integer(bit_kind) :: det_1_active(N_int,2)
  integer          :: h1,h2,p1,p2,s1,s2
  integer          :: exc(0:2,2,2)
  integer          :: i_hole, i_part
  double precision :: phase
  call get_excitation_degree(det_1,det_2,degree,N_int)
  if(degree == 1)then
   call get_excitation(det_1,det_2,exc,degree,phase,N_int)
   call decode_exc(exc,degree,h1,p1,h2,p2,s1,s2)
   i_hole =  list_inact_reverse(h1)
   i_part =  list_virt_reverse(p1)
   do i_state = 1, N_states
     delta_e_act(i_state) += one_anhil_one_creat_inact_virt(i_hole,i_part,i_state) 
   enddo
  endif
 else if (n_holes_act .ge. 2 .and. n_particles_act .ge.2) then
  delta_e_act = -10000000.d0
 endif

!print*, 'one_anhil_spin_trace'
!print*,  one_anhil_spin_trace(1), one_anhil_spin_trace(2) 


 do i_state = 1, n_states
  delta_e_final(i_state) = delta_e_act(i_state)  + delta_e_inactive(i_state) - delta_e_virt(i_state)
 enddo
!write(*,'(100(f16.10,X))'), delta_e_final(1) , delta_e_act(1)  , delta_e_inactive(1) , delta_e_virt(1)

end



subroutine get_delta_e_dyall_general_mp(det_1,det_2,delta_e_final)
 BEGIN_DOC
 ! routine that returns the delta_e with the Moller Plesset and Dyall operators
 !
 ! with det_1 being a determinant from the cas, and det_2 being a perturber
 !
 ! Delta_e(det_1,det_2) = sum (hole) epsilon(hole) + sum(part) espilon(part) + delta_e(act)
 !
 ! where hole is necessary in the inactive, part necessary in the virtuals
 ! 
 ! and delta_e(act) is obtained as the sum of energies of excitations a la MP
 !
 END_DOC
 implicit none
  use bitmasks
 double precision, intent(out) :: delta_e_final(N_states)
 integer(bit_kind), intent(in) :: det_1(N_int,2),det_2(N_int,2)
 integer :: i,j,k,l
 integer :: i_state
 
 integer :: n_holes_spin(2)
 integer :: n_holes
 integer :: holes_list(N_int*bit_kind_size,2)


 double precision :: delta_e_inactive(N_states)
 integer :: i_hole_inact


 call give_holes_in_inactive_space(det_2,n_holes_spin,n_holes,holes_list)
 delta_e_inactive = 0.d0
 do i = 1, n_holes_spin(1)
  i_hole_inact = holes_list(i,1)
  do i_state = 1, N_states
   delta_e_inactive += fock_core_inactive_total_spin_trace(i_hole_inact,i_state)
  enddo
 enddo

 do i = 1, n_holes_spin(2)
  i_hole_inact = holes_list(i,2)
  do i_state = 1, N_states
   delta_e_inactive(i_state) += fock_core_inactive_total_spin_trace(i_hole_inact,i_state)
  enddo
 enddo

 double precision :: delta_e_virt(N_states)
 integer :: i_part_virt
 integer :: n_particles_spin(2)
 integer :: n_particles
 integer :: particles_list(N_int*bit_kind_size,2)
 
 call give_particles_in_virt_space(det_2,n_particles_spin,n_particles,particles_list)
 delta_e_virt = 0.d0
 do i = 1, n_particles_spin(1)
  i_part_virt = particles_list(i,1)
  do i_state = 1, N_states
   delta_e_virt += fock_virt_total_spin_trace(i_part_virt,i_state)
  enddo
 enddo

 do i = 1, n_particles_spin(2)
  i_part_virt = particles_list(i,2)
  do i_state = 1, N_states
   delta_e_virt += fock_virt_total_spin_trace(i_part_virt,i_state)
  enddo
 enddo


 integer :: n_holes_spin_act(2),n_particles_spin_act(2)
 integer :: n_holes_act,n_particles_act
 integer :: holes_active_list(2*n_act_orb,2)
 integer :: holes_active_list_spin_traced(4*n_act_orb)
 integer :: particles_active_list(2*n_act_orb,2)
 integer :: particles_active_list_spin_traced(4*n_act_orb)
 double precision :: delta_e_act(N_states)
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
  hole_list_practical(1,icountbis) = 1  ! spin 
  hole_list_practical(2,icountbis) =  holes_active_list(i,1) ! index of active orb
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

 do i = 1, n_holes_act 
  ispin = hole_list_practical(1,i)
  i_hole_act =  hole_list_practical(2,i)
  do i_state = 1, N_states
   delta_e_act(i_state) += one_anhil(i_hole_act , ispin,i_state)
  enddo
 enddo

 do i = 1, n_particles_act
  ispin = particle_list_practical(1,i)
  i_particle_act =  particle_list_practical(2,i)
  do i_state = 1, N_states
   delta_e_act(i_state) += one_creat(i_particle_act, ispin,i_state)
  enddo
 enddo

 do i_state = 1, n_states
  delta_e_final(i_state) = delta_e_act(i_state)  + delta_e_inactive(i_state) - delta_e_virt(i_state)
 enddo

end

