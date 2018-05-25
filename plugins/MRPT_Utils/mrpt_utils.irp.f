
 BEGIN_PROVIDER [ double precision, delta_ij_mrpt, (N_det_ref,N_det_ref,N_states) ]
&BEGIN_PROVIDER [ double precision, second_order_pt_new, (N_states) ]
&BEGIN_PROVIDER [ double precision, second_order_pt_new_1h, (N_states) ]
&BEGIN_PROVIDER [ double precision, second_order_pt_new_1p, (N_states) ]
&BEGIN_PROVIDER [ double precision, second_order_pt_new_1h1p, (N_states) ]
&BEGIN_PROVIDER [ double precision, second_order_pt_new_2h, (N_states) ]
&BEGIN_PROVIDER [ double precision, second_order_pt_new_2p, (N_states) ]
&BEGIN_PROVIDER [ double precision, second_order_pt_new_1h2p, (N_states) ]
&BEGIN_PROVIDER [ double precision, second_order_pt_new_2h1p, (N_states) ]
&BEGIN_PROVIDER [ double precision, second_order_pt_new_2h2p, (N_states) ]
  use bitmasks
 implicit none
 BEGIN_DOC
 ! Dressing matrix in N_det_ref basis
 END_DOC
 integer :: i,j,m
 integer :: idet,jdet,hole,part
 integer :: i_state,ispin
 double precision :: accu(N_states)
 double precision, allocatable :: delta_ij_tmp(:,:,:)
 integer           :: degree(N_det_ref)
 integer           :: idx(0:N_det_ref)
 integer :: exc(0:2,2,2)
 double precision :: phase
 integer :: occ(N_int*bit_kind_size,2)
 integer  :: n_elec_tmp(2)
 integer :: iorb_a


 delta_ij_mrpt = 0.d0

 allocate (delta_ij_tmp(N_det_ref,N_det_ref,N_states))


 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 ! 1h 
 delta_ij_tmp = 0.d0
 accu = 0.d0
 if(orbital_ordered)then
  call set_bitmask_for_1h
 else 
  call set_bitmask_particl_as_input(reunion_of_bitmask)
  call set_bitmask_hole_as_input(reunion_of_bitmask)
 endif
 call H_apply_mrpt_1h(delta_ij_tmp,N_det_ref)
 do i_state = 1, N_states
 do i = 1, N_det_ref
  do j = 1, N_det_ref
   accu(i_state) += delta_ij_tmp(j,i,i_state) * psi_coef(i,i_state) * psi_coef(j,i_state)
   delta_ij_mrpt(j,i,i_state) += delta_ij_tmp(j,i,i_state)
  enddo
 enddo
   second_order_pt_new_1h(i_state) = accu(i_state) 
 enddo
 print*, '1h   = ',accu
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 ! 1p 
 delta_ij_tmp = 0.d0
 if(orbital_ordered)then
  call set_bitmask_for_1p
 else 
  call set_bitmask_particl_as_input(reunion_of_bitmask)
  call set_bitmask_hole_as_input(reunion_of_bitmask)
 endif
 call H_apply_mrpt_1p(delta_ij_tmp,N_det_ref)
 accu = 0.d0
 do i_state = 1, N_states
 do i = 1, N_det_ref
  do j = 1, N_det_ref
   accu(i_state) += delta_ij_tmp(j,i,i_state) * psi_coef(i,i_state) * psi_coef(j,i_state)
   delta_ij_mrpt(j,i,i_state) += delta_ij_tmp(j,i,i_state)
  enddo
 enddo
 second_order_pt_new_1p(i_state) = accu(i_state) 
 enddo
 print*, '1p   = ',accu
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 ! 1h1p 
 delta_ij_tmp = 0.d0
 if(orbital_ordered)then
  call set_bitmask_for_1h1p
 else 
  call set_bitmask_particl_as_input(reunion_of_bitmask)
  call set_bitmask_hole_as_input(reunion_of_bitmask)
 endif
 call H_apply_mrpt_1h1p(delta_ij_tmp,N_det_ref)
 double precision :: e_corr_from_1h1p_singles(N_states)
 accu = 0.d0
 do i_state = 1, N_states
 do i = 1, N_det_ref
  do j = 1, N_det_ref
   accu(i_state) += delta_ij_tmp(j,i,i_state) * psi_coef(i,i_state) * psi_coef(j,i_state)
   delta_ij_mrpt(j,i,i_state) += delta_ij_tmp(j,i,i_state)
  enddo
 enddo
 second_order_pt_new_1h1p(i_state) = accu(i_state) 
 enddo
 print*, '1h1p = ',accu
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 ! 1h1p third order
 if(do_third_order_1h1p)then
  delta_ij_tmp = 0.d0
  call give_1h1p_sec_order_singles_contrib(delta_ij_tmp)
  accu = 0.d0
  do i_state = 1, N_states
  do i = 1, N_det_ref
   do j = 1, N_det_ref
    accu(i_state) += delta_ij_tmp(j,i,i_state) * psi_coef(i,i_state) * psi_coef(j,i_state)
    delta_ij_mrpt(j,i,i_state) += delta_ij_tmp(j,i,i_state)
   enddo
  enddo
  second_order_pt_new_1h1p(i_state) = accu(i_state) 
  enddo
  print*, '1h1p(3)',accu
 endif
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 ! 2h   
 delta_ij_tmp = 0.d0
 if(orbital_ordered)then
  call set_bitmask_for_2h
 else 
  call set_bitmask_particl_as_input(reunion_of_bitmask)
  call set_bitmask_hole_as_input(reunion_of_bitmask)
 endif
 call H_apply_mrpt_2h(delta_ij_tmp,N_det_ref)
 accu = 0.d0
 do i_state = 1, N_states
 do i = 1, N_det_ref
  do j = 1, N_det_ref
   accu(i_state) += delta_ij_tmp(j,i,i_state) * psi_coef(i,i_state) * psi_coef(j,i_state)
   delta_ij_mrpt(j,i,i_state) += delta_ij_tmp(j,i,i_state)
  enddo
 enddo
 second_order_pt_new_2h(i_state) = accu(i_state) 
 enddo
 print*, '2h   = ',accu
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 ! 2p   
 delta_ij_tmp = 0.d0
 if(orbital_ordered)then
  call set_bitmask_for_2p
 else 
  call set_bitmask_particl_as_input(reunion_of_bitmask)
  call set_bitmask_hole_as_input(reunion_of_bitmask)
 endif
 call H_apply_mrpt_2p(delta_ij_tmp,N_det_ref)
 accu = 0.d0
 do i_state = 1, N_states
 do i = 1, N_det_ref
  do j = 1, N_det_ref
   accu(i_state) += delta_ij_tmp(j,i,i_state) * psi_coef(i,i_state) * psi_coef(j,i_state)
   delta_ij_mrpt(j,i,i_state) += delta_ij_tmp(j,i,i_state)
  enddo
 enddo
 second_order_pt_new_2p(i_state) = accu(i_state) 
 enddo
 print*, '2p   = ',accu
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 ! 1h2p   
!!!QQcall H_apply_mrpt_1h2p(delta_ij_tmp,N_det_ref)  very old way to do it
!!!QQcall give_1h2p_contrib(delta_ij_tmp)   old way to do it 

 accu = 0.d0
 do i_state = 1, N_states
  do idet = 1, N_det_ref
   !!! Diagonal element 
   call bitstring_to_list_ab(psi_active(1,1,idet), occ, n_elec_tmp, N_int)
   do ispin = 1, 2
    do i = 1, n_elec_tmp(ispin)
     iorb_a = list_act_reverse(occ(i,ispin))
     delta_ij_mrpt(idet,idet,i_state) += effective_fock_operator_1h2p(iorb_a,iorb_a,ispin,i_state)
     accu(i_state) += effective_fock_operator_1h2p(iorb_a,iorb_a,ispin,i_state)* psi_coef(idet,i_state)**2
    enddo
   enddo
   !!! Extra diagonal elements 
   call get_excitation_degree_vector_mono(psi_det,psi_det(1,1,idet),degree,N_int,N_det_ref,idx)
   do jdet = 1, idx(0)
    if(idx(jdet)==idet)cycle
    call get_mono_excitation(psi_det(1,1,idet),psi_det(1,1,idx(jdet)),exc,phase,N_int)
    if (exc(0,1,1) == 1) then
       ! Mono alpha
       hole = list_act_reverse(exc(1,1,1))   !!!  a_a
       part = list_act_reverse(exc(1,2,1))   !!!  a^{\dagger}_{b}
       ispin =  1
    else
       ! Mono beta
       hole = list_act_reverse(exc(1,1,2))   !!!  a_a
       part = list_act_reverse(exc(1,2,2))   !!!  a^{\dagger}_{b}
       ispin =  2
    endif
    delta_ij_mrpt(idet,idx(jdet),i_state) += effective_fock_operator_1h2p(hole,part,ispin,i_state) * phase
    accu(i_state) += effective_fock_operator_1h2p(hole,part,ispin,i_state) * phase * psi_coef(idet,i_state) * psi_coef(idx(jdet),i_state)
   enddo
  enddo
  print*, '1h2p =',accu(i_state)
 enddo

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!QQcall H_apply_mrpt_2h1p(delta_ij_tmp,N_det_ref)  very old way to do it
!!!QQcall give_2h1p_contrib(delta_ij_tmp)   old way to do it 


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 integer(bit_kind) :: key_tmp(N_int,2)
 accu = 0.d0
 do i_state = 1, N_states
  do idet = 1, N_det_ref
   !!! Diagonal element 
   do i = 1, N_int
    key_tmp(i,1) = xor(psi_active(i,1,idet),act_bitmask(i,1))
    key_tmp(i,2) = xor(psi_active(i,2,idet),act_bitmask(i,2))
   enddo
   call bitstring_to_list_ab(key_tmp, occ, n_elec_tmp, N_int)
   do ispin = 1, 2
    do i = 1, n_elec_tmp(ispin)
     iorb_a = list_act_reverse(occ(i,ispin))
     delta_ij_mrpt(idet,idet,i_state) += effective_fock_operator_2h1p(iorb_a,iorb_a,ispin,i_state)
     accu(i_state) += effective_fock_operator_2h1p(iorb_a,iorb_a,ispin,i_state)* psi_coef(idet,i_state)**2
    enddo
   enddo
   !!! Extra diagonal elements 
   call get_excitation_degree_vector_mono(psi_det,psi_det(1,1,idet),degree,N_int,N_det_ref,idx)
   do jdet = 1, idx(0)
    if(idx(jdet)==idet)cycle
    call get_mono_excitation(psi_det(1,1,idet),psi_det(1,1,idx(jdet)),exc,phase,N_int)
    if (exc(0,1,1) == 1) then
       ! Mono alpha
       hole = list_act_reverse(exc(1,1,1))   !!!  a_a
       part = list_act_reverse(exc(1,2,1))   !!!  a^{\dagger}_{b}
       ispin =  1
    else
       ! Mono beta
       hole = list_act_reverse(exc(1,1,2))   !!!  a_a
       part = list_act_reverse(exc(1,2,2))   !!!  a^{\dagger}_{b}
       ispin =  2
    endif
    delta_ij_mrpt(idet,idx(jdet),i_state) += - effective_fock_operator_2h1p(hole,part,ispin,i_state) * phase !! the "-" sign comes from the a_a a^{\dagger}_b
    accu(i_state) += - effective_fock_operator_2h1p(hole,part,ispin,i_state) * phase * psi_coef(idet,i_state) * psi_coef(idx(jdet),i_state)!! the "-" sign comes from the a_a a^{\dagger}_b
   enddo
  enddo
 enddo
 print*, '2h1p = ',accu


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 double precision :: contrib_2h2p(N_states)
 call give_2h2p(contrib_2h2p)
 do i_state = 1, N_states
 do i = 1, N_det_ref
   delta_ij_mrpt(i,i,i_state) += contrib_2h2p(i_state)
 enddo
 second_order_pt_new_2h2p(i_state) = contrib_2h2p(i_state) 
 enddo
 print*, '2h2p = ',contrib_2h2p(1) 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 ! total  
 
 print*, 'Hamiltonian dressing'
 accu = 0.d0
 do i_state = 1, N_states
 do i = 1, N_det_ref
  write(*,'(1000(F16.10,x))')delta_ij_mrpt(i,:,:)
  do j = i_state, N_det_ref
   accu(i_state) += delta_ij_mrpt(j,i,i_state) * psi_coef(i,i_state) * psi_coef(j,i_state)
  enddo
 enddo
 second_order_pt_new(i_state) = accu(i_state) 
 print*, 'total= ',accu(i_state)
 enddo




END_PROVIDER

 BEGIN_PROVIDER [double precision, Hmatrix_dressed_pt2_new, (N_det_ref,N_det_ref,N_states)]
 implicit none
 integer :: i,j,i_state
 do i_state = 1, N_states
  do i = 1,N_det_ref
   do j = 1,N_det_ref
    Hmatrix_dressed_pt2_new(j,i,i_state) = H_matrix_all_dets(j,i) + delta_ij_mrpt(j,i,i_state)
   enddo
  enddo
 enddo
 print*, 'Total dressed matrix '
 do i = 1, N_det_ref
  write(*,'(1000(F16.10,x))')Hmatrix_dressed_pt2_new(i,:,:)
 enddo
 END_PROVIDER 

 

 BEGIN_PROVIDER [double precision, Hmatrix_dressed_pt2_new_symmetrized, (N_det_ref,N_det_ref,N_states)]
 implicit none
 integer :: i,j,i_state
 do i_state = 1, N_states
  do i = 1,N_det_ref
   do j = i,N_det_ref
    Hmatrix_dressed_pt2_new_symmetrized(j,i,i_state) =  H_matrix_all_dets(j,i) & 
                                            + 0.5d0 * ( delta_ij_mrpt(j,i,i_state) + delta_ij_mrpt(i,j,i_state) )
    Hmatrix_dressed_pt2_new_symmetrized(i,j,i_state) =  Hmatrix_dressed_pt2_new_symmetrized(j,i,i_state) 
   enddo
  enddo
 enddo
 END_PROVIDER 

  BEGIN_PROVIDER [ double precision, CI_electronic_dressed_pt2_new_energy, (N_states) ]
 &BEGIN_PROVIDER [ double precision, CI_dressed_pt2_new_eigenvectors, (N_det_ref,N_states) ]
 &BEGIN_PROVIDER [ double precision, CI_dressed_pt2_new_eigenvectors_s2, (N_states) ]
  BEGIN_DOC
  ! Eigenvectors/values of the CI matrix
  END_DOC
  implicit none
  double precision :: ovrlp,u_dot_v
  integer :: i_good_state
  integer, allocatable  :: index_good_state_array(:)
  logical, allocatable  :: good_state_array(:)
  double precision, allocatable :: s2_values_tmp(:)
  integer :: i_other_state
  double precision, allocatable  :: eigenvectors(:,:), eigenvalues(:)
  integer :: i_state
  double precision :: s2,e_0
  integer                        :: i,j,k
  double precision, allocatable :: s2_eigvalues(:)
  double precision, allocatable :: e_array(:)
  integer, allocatable :: iorder(:)
  
  ! Guess values for the "N_states" states of the CI_dressed_pt2_new_eigenvectors 
  do j=1,min(N_states,N_det_ref)
    do i=1,N_det_ref
      CI_dressed_pt2_new_eigenvectors(i,j) = psi_coef(i,j)
    enddo
  enddo

  do j=N_det_ref+1,N_states
    do i=1,N_det_ref
      CI_dressed_pt2_new_eigenvectors(i,j) = 0.d0
    enddo
  enddo
  
  if (diag_algorithm == "Davidson") then
    
    print*, 'Davidson not yet implemented for the dressing ... '
    stop
    
   else if (diag_algorithm == "Lapack") then

     allocate (eigenvectors(size(H_matrix_all_dets,1),N_det_ref))
     allocate (eigenvalues(N_det_ref))
     call lapack_diag(eigenvalues,eigenvectors,                      &
         Hmatrix_dressed_pt2_new_symmetrized(1,1,1),size(H_matrix_all_dets,1),N_det_ref)
     CI_electronic_dressed_pt2_new_energy(:) = 0.d0
     if (s2_eig) then
       i_state = 0
       allocate (s2_eigvalues(N_det_ref))
       allocate(index_good_state_array(N_det_ref),good_state_array(N_det_ref))
       good_state_array = .False.
       call u_0_S2_u_0(s2_eigvalues,eigenvectors,N_det_ref,psi_ref,N_int,&
         N_det_ref,size(eigenvectors,1))
       do j=1,N_det_ref
         ! Select at least n_states states with S^2 values closed to "expected_s2"
         if(dabs(s2_eigvalues(j)-expected_s2).le.0.5d0)then
           i_state +=1
           index_good_state_array(i_state) = j
           good_state_array(j) = .True.
         endif
         if(i_state.eq.N_states) then
           exit
         endif
       enddo
       if(i_state .ne.0)then
         ! Fill the first "i_state" states that have a correct S^2 value
         do j = 1, i_state
           do i=1,N_det_ref
             CI_dressed_pt2_new_eigenvectors(i,j) = eigenvectors(i,index_good_state_array(j))
           enddo
           CI_electronic_dressed_pt2_new_energy(j) = eigenvalues(index_good_state_array(j))
           CI_dressed_pt2_new_eigenvectors_s2(j) = s2_eigvalues(index_good_state_array(j))
         enddo
         i_other_state = 0
         do j = 1, N_det_ref
           if(good_state_array(j))cycle
           i_other_state +=1
           if(i_state+i_other_state.gt.n_states)then
             exit
           endif
           do i=1,N_det_ref
             CI_dressed_pt2_new_eigenvectors(i,i_state+i_other_state) = eigenvectors(i,j)
           enddo
           CI_electronic_dressed_pt2_new_energy(i_state+i_other_state) = eigenvalues(j)
           CI_dressed_pt2_new_eigenvectors_s2(i_state+i_other_state) = s2_eigvalues(i_state+i_other_state)
         enddo

       else
         print*,''
         print*,'!!!!!!!!   WARNING  !!!!!!!!!'
         print*,'  Within the ',N_det_ref,'determinants selected'
         print*,'  and the ',N_states,'states requested'
         print*,'  We did not find any state with S^2 values close to ',expected_s2
         print*,'  We will then set the first N_states eigenvectors of the H matrix'
         print*,'  as the CI_dressed_pt2_new_eigenvectors'
         print*,'  You should consider more states and maybe ask for s2_eig to be .True. or just enlarge the CI space'
         print*,''
         do j=1,min(N_states,N_det_ref)
           do i=1,N_det_ref
             CI_dressed_pt2_new_eigenvectors(i,j) = eigenvectors(i,j)
           enddo
           CI_electronic_dressed_pt2_new_energy(j) = eigenvalues(j)
           CI_dressed_pt2_new_eigenvectors_s2(j) = s2_eigvalues(j)
         enddo
       endif
       deallocate(index_good_state_array,good_state_array)
       deallocate(s2_eigvalues)
     else
       call u_0_S2_u_0(CI_dressed_pt2_new_eigenvectors_s2,eigenvectors,N_det_ref,psi_ref,N_int,&
          min(N_det_ref,N_states),size(eigenvectors,1))
       ! Select the "N_states" states of lowest energy
       do j=1,min(N_det_ref,N_states)
         do i=1,N_det_ref
           CI_dressed_pt2_new_eigenvectors(i,j) = eigenvectors(i,j)
         enddo
         CI_electronic_dressed_pt2_new_energy(j) = eigenvalues(j)
       enddo
     endif
     deallocate(eigenvectors,eigenvalues)
   endif

 
END_PROVIDER
 

BEGIN_PROVIDER [ double precision, CI_dressed_pt2_new_energy, (N_states) ]
  implicit none
  BEGIN_DOC
  ! N_states lowest eigenvalues of the CI matrix
  END_DOC
  
  integer                        :: j
  character*(8)                  :: st
  call write_time(6)
  do j=1,N_states
    CI_dressed_pt2_new_energy(j) = CI_electronic_dressed_pt2_new_energy(j) + nuclear_repulsion
    write(st,'(I4)') j
    call write_double(6,CI_dressed_pt2_new_energy(j),'Energy of state '//trim(st))
    call write_double(6,CI_dressed_pt2_new_eigenvectors_s2(j),'S^2 of state '//trim(st))
  enddo

END_PROVIDER
