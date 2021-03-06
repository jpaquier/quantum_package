!BEGIN_PROVIDER [ integer, large_ao_num ]
! implicit none
! BEGIN_DOC
! !Number of large component AO
! END_DOC  
!large_ao_num = ao_num
!END_PROVIDER

!BEGIN_PROVIDER [ integer, large_ao_prim_num, (large_ao_num) ]
!implicit none
! BEGIN_DOC
! !Number of large component primitives
! END_DOC  
! large_ao_prim_num = 1
!END_PROVIDER

!BEGIN_PROVIDER [ integer, large_ao_prim_num_max ]
!implicit none
! BEGIN_DOC
! !max number of primitives of the large component
! END_DOC  
! large_ao_prim_num_max = maxval(large_ao_prim_num)
!END_PROVIDER

!BEGIN_PROVIDER [ integer, large_ao_power, (large_ao_num,3) ]
!&BEGIN_PROVIDER [ integer, large_ao_l, (large_ao_num) ]
!&BEGIN_PROVIDER [ integer, large_ao_nucl, (large_ao_num) ]
!&BEGIN_PROVIDER [ double precision, large_ao_expo_ordered_transp, (large_ao_prim_num_max,large_ao_num) ]
!implicit none
!BEGIN_DOC
!!Nucleus on which the AOs are centered
!!large_ao_power is the AO power for x,y,z for the AOs of 
!! the large component ao basis
!!large_ao_l = l value of : (a+b+c) in x^a y^b z^c for the AOs of 
!! the large component basis
!!large_ao_nucl is the nucleus on which the large component AO is located
!!large_ao_expo_ordered_transp is the transposed ordered large_ao_expo (which
!! are not defined given that they are the non-relativistic  ao_expo ) 
!END_DOC  
! large_ao_nucl = ao_nucl
! large_ao_power = ao_power
! large_ao_l = ao_l
! large_ao_expo_ordered_transp = ao_expo_ordered_transp
!END_PROVIDER


!BEGIN_PROVIDER [ double precision, large_ao_coef, (large_ao_num,large_ao_prim_num_max)]
!&BEGIN_PROVIDER [ double precision, large_ao_coef_normalized_ordered_transp, (large_ao_prim_num_max,large_ao_num)]
!BEGIN_DOC
! !Normalized and ordered coefficient of the large component AOs 
! END_DOC
! double precision               :: large_norm,large_overlap_x,large_overlap_y,large_overlap_z,C_A(3), c
! integer                        :: l, powA(3), nz
! integer                        :: i,j,k
! do i = 1, large_ao_num
!  large_ao_coef(i,1) = 1
! enddo
! nz=100
! C_A(1) = 0.d0
! C_A(2) = 0.d0
! C_A(3) = 0.d0
! large_ao_coef_normalized_ordered_transp = 0.d0
! do i=1, large_ao_num
!  powA(1) = large_ao_power(i,1)
!  powA(2) = large_ao_power(i,2)
!  powA(3) = large_ao_power(i,3)
!  do j = 1, large_ao_prim_num(i)
!   call overlap_gaussian_xyz(C_A,C_A,large_ao_expo_ordered_transp(j,i),large_ao_expo_ordered_transp(j,i),powA,powA,large_overlap_x,large_overlap_y,large_overlap_z,large_norm,nz)
!   large_ao_coef_normalized_ordered_transp(j,i) = large_ao_coef(i,j)/sqrt(large_norm)
!  enddo
! enddo
!END_PROVIDER


!BEGIN_PROVIDER [ double precision, large_ao_ortho_canonical_coef, (large_ao_num,large_ao_num)]
!&BEGIN_PROVIDER [ integer, large_ao_ortho_canonical_num ]
! implicit none
! BEGIN_DOC
!!Matrix of the coefficients of the large component MOs generated by the 
!! orthonormalization by the S^{-1/2} canonical transformation of the large
!! component AOs. 
!!large_ao_ortho_canonical_coef(i,j) = coefficient of the ith ao on the jth
!! large_ao_ortho_canonical orbital
! END_DOC
! integer :: i
! large_ao_ortho_canonical_coef = 0.d0
! do i=1, large_ao_num
!   large_ao_ortho_canonical_coef(i,i) = 1.d0
! enddo
! large_ao_ortho_canonical_num = large_ao_num
! call ortho_canonical(large_ao_overlap,size(large_ao_overlap,1), large_ao_num,large_ao_ortho_canonical_coef, size(large_ao_ortho_canonical_coef,1), large_ao_ortho_canonical_num)
!END_PROVIDER
! 
!BEGIN_PROVIDER [double precision, large_ao_ortho_canonical_overlap, (large_ao_ortho_canonical_num, large_ao_ortho_canonical_num)]
! implicit none
! BEGIN_DOC
!!Overlap matrix of the large_ao_ortho_canonical.
!!Expected to be the Identity
! END_DOC
! integer                        :: i,j,k,l
! double precision               :: c
! do j=1, large_ao_ortho_canonical_num
!   do i=1, large_ao_ortho_canonical_num
!     large_ao_ortho_canonical_overlap(i,j) = 0.d0
!   enddo
! enddo
! do j=1, large_ao_ortho_canonical_num
!   do k=1, large_ao_num
!     c = 0.d0
!     do l=1, large_ao_num
!       c +=  large_ao_ortho_canonical_coef(l,j) * large_ao_overlap(l,k)
!     enddo
!     do i=1, large_ao_ortho_canonical_num
!       large_ao_ortho_canonical_overlap(i,j) += large_ao_ortho_canonical_coef(k,i) * c
!     enddo
!   enddo
! enddo
!END_PROVIDER

!BEGIN_PROVIDER [ integer, large_mo_tot_num ]
! implicit none
! BEGIN_DOC
! !Number of large component MOs
! END_DOC
! large_mo_tot_num = large_ao_ortho_canonical_num
! ASSERT (large_mo_tot_num > 0)
!END_PROVIDER


!BEGIN_PROVIDER [ double precision, large_mo_coef, (large_ao_num,large_mo_tot_num) ]
! implicit none
! BEGIN_DOC
! !Molecular orbital coefficients on large component AO basis set
! ! large_mo_coef(i,j) = coefficient of the ith large component ao on the jth
! ! large component mo
! END_DOC
! integer                        :: i, j
! double precision, allocatable  :: buffer(:,:)
! logical                        :: exists
! do i=1,large_mo_tot_num
!   do j=1,large_ao_num
!     large_mo_coef(j,i) = large_ao_ortho_canonical_coef(j,i)
!   enddo
! enddo
!END_PROVIDER



!BEGIN_PROVIDER [double precision, large_ao_ortho_lowdin_coef, (large_ao_num,large_ao_num)]
! implicit none
! BEGIN_DOC
! !Matrix of the coefficients of the MOs generated by the Lowdin orthonormalization of the large component AOs
! ! large_ao_ortho_lowdin_coef(i,j) = coefficient of the ith ao on the jth large_ao_ortho_lowdin orbital
! END_DOC
! integer                        :: i,j,k,l
! double precision               :: accu
! double precision, allocatable  :: tmp_matrix(:,:)
! allocate (tmp_matrix(large_ao_num,large_ao_num))
! tmp_matrix(:,:) = 0.d0
! do j=1, large_ao_num
!   tmp_matrix(j,j) = 1.d0
! enddo
! call ortho_lowdin(large_ao_overlap,large_ao_num,large_ao_num,tmp_matrix,large_ao_num,large_ao_num)
! do i=1, large_ao_num
!   do j=1, large_ao_num
!     large_ao_ortho_lowdin_coef(j,i) = tmp_matrix(i,j)
!   enddo
! enddo
! deallocate(tmp_matrix)
!END_PROVIDER

!BEGIN_PROVIDER [double precision, large_ao_ortho_lowdin_overlap, (large_ao_num,large_ao_num)]
! implicit none
! BEGIN_DOC
! !Overlap matrix of the large_ao_ortho_lowdin
! ! supposed to be the Identity
! END_DOC
! integer                        :: i,j,k,l
! double precision               :: c
! do j=1, large_ao_num
!   do i=1, large_ao_num
!     large_ao_ortho_lowdin_overlap(i,j) = 0.d0
!   enddo
! enddo
! do k=1, large_ao_num
!   do j=1, large_ao_num
!     c = 0.d0
!     do l=1, large_ao_num
!       c +=  large_ao_ortho_lowdin_coef(j,l) * large_ao_overlap(k,l)
!     enddo
!     do i=1, large_ao_num
!       large_ao_ortho_lowdin_overlap(i,j) += large_ao_ortho_lowdin_coef(i,k) * c
!     enddo
!   enddo
! enddo
!END_PROVIDER


