!BEGIN_PROVIDER [ integer, number_of_small_component_expo_per_shell_per_atom,(0:7,nucl_num)] 
!&BEGIN_PROVIDER [ integer, nmax_of_small_component_expo] 
!&BEGIN_PROVIDER [ integer, number_of_small_component_ao_per_shell_per_atom,(0:7,nucl_num)]
!&BEGIN_PROVIDER [ integer, number_of_small_component_ao_per_atom,(nucl_num)]
!&BEGIN_PROVIDER [ integer, small_ao_num]
!implicit none
!BEGIN_DOC
!!Number_of_small_component_expo_per_shell_per_atom = number of ao exponents in
!! the shell "i" of the nucleus "j" for the small component basis in the 
!! unrestricted kinetic balance scheme
!END_DOC                                                                                                 
!integer :: i,i_count,l_type,j,index_ao_previous,index_ao,k,l,l_count                                                      
!integer :: imax
!index_ao_previous = 0 
!index_ao = 0
!imax = 0
!number_of_small_component_expo_per_shell_per_atom(:,:) = 0
!number_of_small_component_ao_per_shell_per_atom(:,:) = 0
!number_of_small_component_ao_per_atom(:) = 0
!small_ao_num = 0
!do i = 1, nucl_num
! do l_type = 0,7
!  if(Nucl_num_l_type_Aos(l_type,i) == 0)cycle                                                           
!  do j = 1, Nucl_num_l_type_Aos(l_type,i)
!   !"physical" index of all the AO corresponding to the shell type l_type                              
!   ! attached to atom i                                                                                 
!   index_ao_previous = index_ao
!   index_ao = Nucl_list_l_type_Aos(j,l_type,i)                                                          
!   if (l_type == 0) then
!    number_of_small_component_expo_per_shell_per_atom(l_type+1,i) += 1
!   else
!    if (j == 1) then
!     number_of_small_component_expo_per_shell_per_atom(l_type-1,i) += 1
!     number_of_small_component_expo_per_shell_per_atom(l_type+1,i) += 1
!    elseif ( j /= 1 .and. ao_expo_ordered_transp(1,index_ao) /= ao_expo_ordered_transp(1,index_ao_previous)) then
!     number_of_small_component_expo_per_shell_per_atom(l_type-1,i) += 1
!     number_of_small_component_expo_per_shell_per_atom(l_type+1,i) += 1
!    endif                                  
!   endif
!  enddo
! enddo 
! !Second loop to count the number of AOs of the small component basis
! do l_type = 0,7  
!  if(imax .lt. number_of_small_component_expo_per_shell_per_atom(l_type,i)) then
!   imax =  number_of_small_component_expo_per_shell_per_atom(l_type,i)
!  endif 
!  number_of_small_component_ao_per_shell_per_atom(l_type,i) += number_of_small_component_expo_per_shell_per_atom(l_type,i)*(l_type+1)*(l_type+2)/2 
!  number_of_small_component_ao_per_atom(i) += number_of_small_component_ao_per_shell_per_atom(l_type,i) 
! enddo
! small_ao_num += number_of_small_component_ao_per_atom(i)
!enddo
!nmax_of_small_component_expo = imax
!END_PROVIDER


!BEGIN_PROVIDER [ double precision, small_component_expo_per_shell_per_atom,(nmax_of_small_component_expo,0:7,nucl_num)]
!implicit none
!BEGIN_DOC
!!Arrays containing the values of the exponents for the small component in the
!! unrestricted kinetic balance scheme
!END_DOC                                                                       
!integer :: i,l_type,j,index_ao_previous,index_ao,k,l,l_count
!integer :: expo_count(0:7) 
!integer :: iorder(nmax_of_small_component_expo)
!double precision :: small_component_expo_per_shell(nmax_of_small_component_expo,0:7)
!double precision :: small_component_expo(nmax_of_small_component_expo)
!small_component_expo_per_shell_per_atom(:,:,:) = 0
!small_component_expo_per_shell(:,:) = 0
!index_ao_previous = 0
!index_ao = 0
!do i = 1, nucl_num
! iorder(:) = 0
! expo_count(:) = 0
! do l_type = 0,7
!  if(Nucl_num_l_type_Aos(l_type,i) == 0)cycle
!  do j = 1, Nucl_num_l_type_Aos(l_type,i)
!   !"physical" index of all the AO corresponding to the shell type l_type
!   ! attached to atom i
!   index_ao_previous = index_ao 
!   index_ao = Nucl_list_l_type_Aos(j,l_type,i)
!   if (l_type == 0) then
!    expo_count(l_type+1) += 1
!    small_component_expo_per_shell(expo_count(l_type+1),l_type+1) = ao_expo_ordered_transp(1,index_ao)
!   else
!    if (j == 1) then
!     expo_count(l_type-1) += 1                              
!     small_component_expo_per_shell(expo_count(l_type-1),l_type-1) = ao_expo_ordered_transp(1,index_ao) 
!     expo_count(l_type+1) += 1
!     small_component_expo_per_shell(expo_count(l_type+1),l_type+1) = ao_expo_ordered_transp(1,index_ao) 
!    elseif ( j /= 1 .and. ao_expo_ordered_transp(1,index_ao) /=ao_expo_ordered_transp(1,index_ao_previous)) then
!     expo_count(l_type-1) += 1
!     small_component_expo_per_shell(expo_count(l_type-1),l_type-1) = ao_expo_ordered_transp(1,index_ao)
!     expo_count(l_type+1) += 1
!     small_component_expo_per_shell(expo_count(l_type+1),l_type+1) = ao_expo_ordered_transp(1,index_ao) 
!    endif
!   endif
!  enddo
! enddo 
! !Second loop to sort the AO exponents in the right order
! do l_type = 0,7
!  small_component_expo(:)=0
!  do k = 1, expo_count(l_type)
!   iorder(k) = k
!   small_component_expo(k) = small_component_expo_per_shell(k,l_type)
!  enddo
!   call dsort(small_component_expo,iorder,expo_count(l_type))
!  do l =1, expo_count(l_type)
!   small_component_expo_per_shell_per_atom(l,l_type,i) = small_component_expo(expo_count(l_type)+1-l)
!  enddo 
! enddo 
!enddo
!END_PROVIDER

!BEGIN_PROVIDER [ integer, small_ao_prim_num, (small_ao_num) ]
!implicit none
!BEGIN_DOC
!!number of small component primitives
!END_DOC
!integer :: i
!do i = 1, small_ao_num
! small_ao_prim_num(i) = 1
!enddo
!END_PROVIDER

!BEGIN_PROVIDER [ integer, small_ao_prim_num_max ]
!implicit none
!BEGIN_DOC
!!max number of primitives of the small component
!END_DOC
!small_ao_prim_num_max = maxval(small_ao_prim_num)
!END_PROVIDER



!BEGIN_PROVIDER [ integer, small_ao_l, (small_ao_num) ]
!&BEGIN_PROVIDER [ integer, small_ao_l_max  ]
!&BEGIN_PROVIDER [ integer, small_ao_nucl, (small_ao_num) ]  
!&BEGIN_PROVIDER [ double precision, small_ao_expo, (small_ao_num) ]
!&BEGIN_PROVIDER [ double precision, small_ao_expo_ordered_transp, (small_ao_prim_num_max, small_ao_num) ]
!implicit none
!BEGIN_DOC
!!small_ao_l = l value of : (a+b+c) in x^a y^b z^c for the AOs of 
!! the small component basis
!!small_ao_l_max is the max value of small_ao_l
!!small_ao_nucl is the nucleus on which the small component AO is located
!!small_ao_expo is the exponent of the small component AO
!!small_ao_expo_ordered_transp is the transposed ordered small_ao_expo
!END_DOC
!integer :: i,l,l_type,j,j_count,k,k_count
!j_count = 0
!k_count = 0
!do i = 1,nucl_num
! do l_type = 0,7
!  do j = 1, number_of_small_component_ao_per_shell_per_atom(l_type,i)  
!   j_count += 1
!   small_ao_l(j_count) = l_type
!   small_ao_nucl(j_count) = i 
!  enddo
!  do k = 1, number_of_small_component_expo_per_shell_per_atom(l_type,i)
!   do l =1, (l_type+1)*(l_type+2)/2
!    k_count += 1
!    small_ao_expo(k_count) = small_component_expo_per_shell_per_atom(k,l_type,i)
!   enddo
!  enddo
! enddo
!enddo
!do i = 1, small_ao_num
! small_ao_expo_ordered_transp(1,i) = small_ao_expo (i)
!enddo
!small_ao_l_max = maxval(small_ao_l)
!END_PROVIDER


!BEGIN_PROVIDER [ integer, small_ao_power, (small_ao_num,3)]
!implicit none
!BEGIN_DOC
!!small_ao_power is the AO power for x,y,z for the AOs of 
!! the small component ao basis
!END_DOC
!integer :: i,j,j_count,k,k_count,k_j,l,l_type,l_count
!small_ao_power(:,:) = 0
!k_count = 0
!k_j = 0
!do i = 1,nucl_num
! do l_type = 0,7
!  do k = 1, number_of_small_component_expo_per_shell_per_atom(l_type,i)
!   do j =1, l_type + 1
!    if (j == 1) then
!     k_count += 1
!     k_j = k_count
!     small_ao_power(k_count,1) = l_type
!    else
!     do l = 1, j
!      k_count += 1
!      l_count += 1  
!      small_ao_power(k_count,1) = small_ao_power(k_j,1) - j +1
!      small_ao_power(k_count,2) = small_ao_power(k_j,2) + j -l 
!      small_ao_power(k_count,3) = small_ao_power(k_j,3) +l -1 
!     enddo
!    endif
!   enddo
!  enddo
! enddo 
!enddo  
!END_PROVIDER
!
! 
!BEGIN_PROVIDER [ double precision, small_ao_coef, (small_ao_num,small_ao_prim_num_max)]
!&BEGIN_PROVIDER [ double precision, small_ao_coef_normalized_ordered_transp, (small_ao_prim_num_max,small_ao_num)]
!BEGIN_DOC
! !Normalized and ordered coefficient of the small component AOs in the
! ! unrestricted kinetic balance scheme
! END_DOC
! double precision               :: small_norm,small_overlap_x,small_overlap_y,small_overlap_z,C_A(3), c
! integer                        :: l, powA(3), nz
! integer                        :: i,j,k
! do i = 1, small_ao_num
!  small_ao_coef(i,1) = 1.d0
! enddo
! nz=100
! C_A(1) = 0.d0
! C_A(2) = 0.d0
! C_A(3) = 0.d0
! small_ao_coef_normalized_ordered_transp = 0.d0
! do i=1, small_ao_num
!  powA(1) = small_ao_power(i,1)
!  powA(2) = small_ao_power(i,2)
!  powA(3) = small_ao_power(i,3) 
!  do j = 1, small_ao_prim_num(i)
!   call overlap_gaussian_xyz(C_A,C_A,small_ao_expo_ordered_transp(j,i),small_ao_expo_ordered_transp(j,i),powA,powA,small_overlap_x,small_overlap_y,small_overlap_z,small_norm,nz)
!   small_ao_coef_normalized_ordered_transp(j,i) = small_ao_coef(i,j)/sqrt(small_norm)
!  enddo
! enddo
!END_PROVIDER
!

!BEGIN_PROVIDER [ double precision, small_ao_ortho_canonical_coef, (small_ao_num,small_ao_num)]
!&BEGIN_PROVIDER [ integer, small_ao_ortho_canonical_num ]
! implicit none
! BEGIN_DOC
!!Matrix of the coefficients of the small component MOs generated by the 
!! orthonormalization by the S^{-1/2} canonical transformation of the small component AOs. 
!!small_ao_ortho_canonical_coef(i,j) = coefficient of the ith small component AO
!! on the jth small_ao_ortho_canonical orbital
! END_DOC
! integer :: i
! small_ao_ortho_canonical_coef = 0.d0
! do i=1, small_ao_num
!   small_ao_ortho_canonical_coef(i,i) = 1.d0
! enddo
! small_ao_ortho_canonical_num = small_ao_num
! call ortho_canonical(small_ao_overlap,size(small_ao_overlap,1), small_ao_num,small_ao_ortho_canonical_coef, size(small_ao_ortho_canonical_coef,1), small_ao_ortho_canonical_num)
!END_PROVIDER
! 
!BEGIN_PROVIDER [double precision, small_ao_ortho_canonical_overlap, (small_ao_ortho_canonical_num, small_ao_ortho_canonical_num)]
! implicit none
! BEGIN_DOC
!!Overlap matrix of the small_ao_ortho_canonical.
!!Expected to be the Identity
! END_DOC
! integer                        :: i,j,k,l
! double precision               :: c
! do j=1, small_ao_ortho_canonical_num
!   do i=1, small_ao_ortho_canonical_num
!     small_ao_ortho_canonical_overlap(i,j) = 0.d0
!   enddo
! enddo
! do j=1, small_ao_ortho_canonical_num
!   do k=1, small_ao_num
!     c = 0.d0
!     do l=1, small_ao_num
!       c +=  small_ao_ortho_canonical_coef(l,j) * small_ao_overlap(l,k)
!     enddo
!     do i=1, small_ao_ortho_canonical_num
!       small_ao_ortho_canonical_overlap(i,j) += small_ao_ortho_canonical_coef(k,i) * c
!     enddo
!   enddo
! enddo
!END_PROVIDER

!BEGIN_PROVIDER [ integer, small_mo_tot_num ]
! implicit none
! BEGIN_DOC
! !Number of small component MOs
! END_DOC
! small_mo_tot_num = small_ao_ortho_canonical_num
! ASSERT (small_mo_tot_num > 0)
!END_PROVIDER


!BEGIN_PROVIDER [ double precision, small_mo_coef, (small_ao_num,small_mo_tot_num) ]
! implicit none
! BEGIN_DOC
! !Molecular orbital coefficients on small component AO basis set
! !small_mo_coef(i,j) = coefficient of the ith small component AO on the jth small component MO
! END_DOC
! integer                        :: i, j
! double precision, allocatable  :: buffer(:,:)
! logical                        :: exists
! do i=1,small_mo_tot_num
!   do j=1,small_ao_num
!     small_mo_coef(j,i) = small_ao_ortho_canonical_coef(j,i)
!   enddo
! enddo
!END_PROVIDER


!BEGIN_PROVIDER [double precision, small_ao_ortho_lowdin_coef, (small_ao_num,small_ao_num)]
! implicit none
! BEGIN_DOC
! !Matrix of the coefficients of the mos generated by the Lowdin orthonormalization of the small component AOs
! !ao_ortho_lowdin_coef(i,j) = coefficient of the ith AO on the jth ao_ortho_lowdin orbital
! END_DOC
! integer                        :: i,j,k,l
! double precision               :: accu
! double precision, allocatable  :: tmp_matrix(:,:)
! allocate (tmp_matrix(small_ao_num,small_ao_num))
! tmp_matrix(:,:) = 0.d0
! do j=1, small_ao_num
!   tmp_matrix(j,j) = 1.d0
! enddo
! call ortho_lowdin(small_ao_overlap,small_ao_num,small_ao_num,tmp_matrix,small_ao_num,small_ao_num)
! do i=1, small_ao_num
!   do j=1, small_ao_num
!     small_ao_ortho_lowdin_coef(j,i) = tmp_matrix(i,j)
!   enddo
! enddo
! deallocate(tmp_matrix)
!END_PROVIDER

!BEGIN_PROVIDER [double precision, small_ao_ortho_lowdin_overlap, (small_ao_num,small_ao_num)]
! implicit none
! BEGIN_DOC
! !Overlap matrix of the ao_ortho_lowdin
! !Supposed to be the Identity
! END_DOC
! integer                        :: i,j,k,l
! double precision               :: c
! do j=1, small_ao_num
!  do i=1,small_ao_num
!   small_ao_ortho_lowdin_overlap(i,j) = (0.d0,0.d0)
!  enddo
! enddo
! do k=1, small_ao_num
!  do j=1, small_ao_num
!   c = (0.d0,0.d0)
!   do l=1, small_ao_num
!    c +=  small_ao_ortho_lowdin_coef(j,l) * small_ao_overlap(k,l)
!   enddo
!   do i=1, small_ao_num
!    small_ao_ortho_lowdin_overlap(i,j) += small_ao_ortho_lowdin_coef(i,k) * c
!   enddo
!  enddo
! enddo
!END_PROVIDER
