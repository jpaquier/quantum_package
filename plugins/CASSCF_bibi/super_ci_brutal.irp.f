
 use bitmasks

 BEGIN_PROVIDER [integer(bit_kind), psi_superci, (N_int,2,N_det,2,size_super_ci)]
&BEGIN_PROVIDER [integer(bit_kind), psi_superci_size, (2,size_super_ci)]
&BEGIN_PROVIDER [double precision, contraction_coef_psi_superci, (N_det,2,size_super_ci,N_states)]
 implicit none
 BEGIN_DOC
! psi_superci(:,:,i,1,j) = alpha component of the excitation 
! psi_superci(:,:,i,2,j) = beta  component of the excitation 
 END_DOC


 use bitmasks
 integer  :: i,j,k,l
 integer :: i_State
 integer(bit_kind) :: det_tmp(N_int,2)
 i = 1
 do j = 1, N_det
  do k = 1, N_int
   psi_superci(k,1,j,1,i) = psi_det(k,1,j)
   psi_superci(k,2,j,1,i) = psi_det(k,2,j)
  enddo
  do i_state = 1, N_states
   contraction_coef_psi_superci(j,1,i,i_state) = psi_coef(j,i_state)
  enddo
 enddo
 
 integer :: i_core,i_hole
 integer :: j_virt,j_part
 integer :: i_ok
 double precision :: norm(N_states),phase
 integer :: Ndet
 integer         :: exc(0:2,2,2)
 integer         :: degree
 do i = 2, size_super_ci
  i_core = index_rotation_CI_reverse(i,1)
  i_hole = list_core_inact(i_core)
  j_virt = index_rotation_CI_reverse(i,2)
  j_part = list_virt(j_virt)
  norm = 0.d0

  psi_superci_size(1,i) = 0
  psi_superci_size(2,i) = 0
  do j = 1, N_det
   do k = 1, N_int
    det_tmp(k,1) = psi_det(k,1,j)
    det_tmp(k,2) = psi_det(k,2,j)
   enddo
   call do_mono_excitation(det_tmp,i_hole,j_part,1,i_ok)
   if(i_ok ==1)then
    do k = 1, N_int
     psi_superci(k,1,j,1,i) = det_tmp(k,1)
     psi_superci(k,2,j,1,i) = det_tmp(k,2)
    enddo
    psi_superci_size(1,i) += 1
    call get_excitation(det_tmp,psi_det(1,1,j),exc,degree,phase,N_int)
    do i_State = 1, N_states
     contraction_coef_psi_superci(j,1,i,i_state) = psi_coef(j,i_state) * phase
    enddo
   else 
    do k = 1, N_int
     psi_superci(k,1,j,1,i) = 0_bit_kind
     psi_superci(k,2,j,1,i) = 0_bit_kind
    enddo
    do i_state = 1, N_states
     contraction_coef_psi_superci(j,1,i,i_state) = 0.d0
    enddo
   endif
   do i_state = 1, N_states 
    norm(i_state) += contraction_coef_psi_superci(j,1,i,i_state) * contraction_coef_psi_superci(j,1,i,i_state)
   enddo
  enddo
  do j = 1, N_det
   do k = 1, N_int
    det_tmp(k,1) = psi_det(k,1,j)
    det_tmp(k,2) = psi_det(k,2,j)
   enddo
   call do_mono_excitation(det_tmp,i_hole,j_part,2,i_ok)
   if(i_ok ==1)then
    do k = 1, N_int
     psi_superci(k,1,j,2,i) = det_tmp(k,1)
     psi_superci(k,2,j,2,i) = det_tmp(k,2)
    enddo
    psi_superci_size(2,i) += 1
    call get_excitation(det_tmp,psi_det(1,1,j),exc,degree,phase,N_int)
    do i_state = 1, N_states
     contraction_coef_psi_superci(j,2,i,i_state) = psi_coef(j,i_state) * phase 
    enddo
   else 
    do k = 1, N_int
     psi_superci(k,1,j,2,i) = 0_bit_kind
     psi_superci(k,2,j,2,i) = 0_bit_kind
    enddo
    do i_state = 1, N_states
     contraction_coef_psi_superci(j,2,i,i_state) = 0.d0
    enddo
   endif
   do i_state = 1, N_states
    norm(i_state) += contraction_coef_psi_superci(j,2,i,i_state) * contraction_coef_psi_superci(j,2,i,i_state)
   enddo
  enddo
  norm = 1.d0/dsqrt(norm)
  do i_state = 1, N_states 
   do j = 1, N_det
    contraction_coef_psi_superci(j,1,i,i_state) *= norm(i_state)
    contraction_coef_psi_superci(j,2,i,i_state) *= norm(i_state)
   enddo
  enddo

 enddo


END_PROVIDER 


 BEGIN_PROVIDER [double precision, diagonal_matrix_superci_brutal, (size_super_ci,N_states)]
&BEGIN_PROVIDER [double precision, brillouin_matrix_superci_brutal, (size_super_ci,N_states)]
&BEGIN_PROVIDER [double precision,total_matrix_superci_brutal, (size_super_ci,size_super_ci,N_states)]
 use bitmasks
 implicit none
 integer :: i,j,k,l,idet
 double precision :: hij,  accu
 integer(bit_kind), allocatable :: psi_1(:,:,:)
 double precision, allocatable :: psi_1_coef(:)
 integer(bit_kind), allocatable :: psi_2(:,:,:)
 double precision, allocatable :: psi_2_coef(:)
 integer :: i_state
 i = 1 
 do i_state = 1, N_states
  diagonal_matrix_superci_brutal(1,i_state) = 0.d0
  total_matrix_superci_brutal(1,1,i_state) = 0.d0 
  do i = 2, size_super_ci
   allocate(psi_1(N_int,2,psi_superci_size(1,i)+psi_superci_size(2,i)),psi_1_coef(psi_superci_size(1,i)+psi_superci_size(2,i)))
   idet = 0
   do j = 1, N_det
    if(psi_superci(1,1,j,1,i) == 0_bit_kind)cycle
     idet+=1
     do k = 1, N_int
      psi_1(k,1,idet) = psi_superci(k,1,j,1,i)
      psi_1(k,2,idet) = psi_superci(k,2,j,1,i)
     enddo
!!   call debug_Det(psi_1(1,1,idet),N_int)
     psi_1_coef(idet) = contraction_coef_psi_superci(j,1,i,i_state)
   enddo
   do j = 1, N_det
    if(psi_superci(1,1,j,2,i) == 0_bit_kind)cycle
     idet+=1
     do k = 1, N_int
      psi_1(k,1,idet) = psi_superci(k,1,j,2,i)
      psi_1(k,2,idet) = psi_superci(k,2,j,2,i)
     enddo
     psi_1_coef(idet) = contraction_coef_psi_superci(j,2,i,i_state)
   enddo
 
   ! < Psi_1 | H | Psi_1> 
   accu = 0.d0
   do k = 1, idet  
    do j = 1, idet
     call i_H_j(psi_1(1,1,k),psi_1(1,1,j),N_int,hij) 
     accu += psi_1_coef(k) * psi_1_coef(j) * hij
    enddo
   enddo
   diagonal_matrix_superci_brutal(i,i_State) = accu - psi_energy(i_state)
   total_matrix_superci_brutal(i,i,i_state) = accu - psi_energy(i_state)
   ! < Psi_0 | H | Psi_1> 
   accu = 0.d0
   do k = 1, N_det
    do j = 1, idet
     call i_H_j(psi_det(1,1,k),psi_1(1,1,j),N_int,hij) 
     accu += psi_1_coef(j) * psi_coef(k,i_state) * hij
    enddo
   enddo
   brillouin_matrix_superci_brutal(i,i_State) = accu 
   total_matrix_superci_brutal(i,1,i_state) = accu 
   total_matrix_superci_brutal(1,i,i_state) = accu 
 
   integer :: jdet
   do l = i+1, size_super_ci
     allocate(psi_2(N_int,2,psi_superci_size(1,l)+psi_superci_size(2,l)),psi_2_coef(psi_superci_size(1,l)+psi_superci_size(2,l)))
     jdet = 0
     do j = 1, N_det
      if(psi_superci(1,1,j,1,l) == 0_bit_kind)cycle
       jdet+=1
       do k = 1, N_int
        psi_2(k,1,jdet) = psi_superci(k,1,j,1,l)
        psi_2(k,2,jdet) = psi_superci(k,2,j,1,l)
       enddo
       psi_2_coef(jdet) = contraction_coef_psi_superci(j,1,l,i_state)
     enddo
     do j = 1, N_det
      if(psi_superci(1,1,j,2,i) == 0_bit_kind)cycle
       jdet+=1
       do k = 1, N_int
        psi_2(k,1,jdet) = psi_superci(k,1,j,2,l)
        psi_2(k,2,jdet) = psi_superci(k,2,j,2,l)
       enddo
       psi_2_coef(jdet) = contraction_coef_psi_superci(j,2,l,i_state)
     enddo
     ! < Psi_1 | H | Psi_1> 
     accu = 0.d0
     do k = 1, idet  
      do j = 1, jdet
       call i_H_j(psi_1(1,1,k),psi_2(1,1,j),N_int,hij) 
       accu += psi_1_coef(k) * psi_2_coef(j) * hij
      enddo
     enddo
     total_matrix_superci_brutal(i,l,i_state) = accu
     total_matrix_superci_brutal(l,i,i_state) = accu
     deallocate(psi_2,psi_2_coef)
     
   enddo
 
   deallocate(psi_1,psi_1_coef)
  enddo
 enddo

END_PROVIDER 




 BEGIN_PROVIDER [double precision, eigenvectors_sci_brutal, (size_super_ci,size_super_ci)]
&BEGIN_PROVIDER [double precision, eigenvalues_sci_brutal, (size_super_ci)]
 implicit none
 integer :: i_state

 double precision, allocatable :: matrix_tmp(:,:),eigenvectors(:,:),eigenvalues(:)
 allocate(matrix_tmp(size_super_ci,size_super_ci),eigenvectors(size_super_ci,size_super_ci),eigenvalues(size_super_ci))
 do i_state = 1, N_states
  matrix_tmp(:,:) = total_matrix_superci_brutal(:,:,i_state)
  call lapack_diag(eigenvalues,eigenvectors,matrix_tmp,size_super_ci,size_super_ci)
  eigenvectors_sci_brutal(:,i_state) = eigenvectors(:,1)
  eigenvalues_sci_brutal(i_state) = eigenvalues(1)
 enddo


END_PROVIDER 
