
 use bitmasks

 BEGIN_PROVIDER [integer(bit_kind), psi_superci, (N_int,2,N_det,2,size_super_ci)]
&BEGIN_PROVIDER [integer(bit_kind), psi_superci_size, (2,size_super_ci)]
&BEGIN_PROVIDER [double precision, contraction_coef_psi_superci, (N_det,2,size_super_ci)]
 implicit none
 BEGIN_DOC
! psi_superci(:,:,i,1,j) = alpha component of the excitation 
! psi_superci(:,:,i,2,j) = beta  component of the excitation 
 END_DOC


 use bitmasks
 integer  :: i,j,k,l
 integer(bit_kind) :: det_tmp(N_int,2)
 i = 1
 do j = 1, N_det
  do k = 1, N_int
   psi_superci(k,1,j,1,i) = psi_det(k,1,j)
   psi_superci(k,2,j,1,i) = psi_det(k,2,j)
  enddo
  contraction_coef_psi_superci(j,1,i) = psi_coef(j,1)
 enddo
 
 integer :: i_core,i_hole
 integer :: j_virt,j_part
 integer :: i_ok
 double precision :: norm,phase
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
    contraction_coef_psi_superci(j,1,i) = psi_coef(j,1) * phase
   else 
    do k = 1, N_int
     psi_superci(k,1,j,1,i) = 0_bit_kind
     psi_superci(k,2,j,1,i) = 0_bit_kind
    enddo
    contraction_coef_psi_superci(j,1,i) = 0.d0
   endif
   norm += contraction_coef_psi_superci(j,1,i) * contraction_coef_psi_superci(j,1,i)
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
    contraction_coef_psi_superci(j,2,i) = psi_coef(j,1) * phase 
   else 
    do k = 1, N_int
     psi_superci(k,1,j,2,i) = 0_bit_kind
     psi_superci(k,2,j,2,i) = 0_bit_kind
    enddo
    contraction_coef_psi_superci(j,2,i) = 0.d0
   endif
   norm += contraction_coef_psi_superci(j,1,i) * contraction_coef_psi_superci(j,1,i)
  enddo
  norm = 1.d0/dsqrt(norm)
  do j = 1, N_det
   contraction_coef_psi_superci(j,1,i) *= norm
   contraction_coef_psi_superci(j,2,i) *= norm
  enddo

 enddo


END_PROVIDER 


 BEGIN_PROVIDER [double precision, diagonal_matrix_superci_brutal, (size_super_ci)]
&BEGIN_PROVIDER [double precision,brillouin_matrix_superci_brutal, (size_super_ci)]
&BEGIN_PROVIDER [double precision,total_matrix_superci_brutal, (size_super_ci,size_super_ci)]
 use bitmasks
 implicit none
 integer :: i,j,k,l,idet
 double precision :: hij,  accu
 integer(bit_kind), allocatable :: psi_1(:,:,:)
 double precision, allocatable :: psi_1_coef(:)
 integer(bit_kind), allocatable :: psi_2(:,:,:)
 double precision, allocatable :: psi_2_coef(:)
 i = 1 
 diagonal_matrix_superci_brutal(1) = 0.d0
 total_matrix_superci_brutal(1,1) = 0.d0 
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
!   call debug_Det(psi_1(1,1,idet),N_int)
    psi_1_coef(idet) = contraction_coef_psi_superci(j,1,i)
  enddo
  do j = 1, N_det
   if(psi_superci(1,1,j,2,i) == 0_bit_kind)cycle
    idet+=1
    do k = 1, N_int
     psi_1(k,1,idet) = psi_superci(k,1,j,2,i)
     psi_1(k,2,idet) = psi_superci(k,2,j,2,i)
    enddo
    psi_1_coef(idet) = contraction_coef_psi_superci(j,2,i)
  enddo

  ! < Psi_1 | H | Psi_1> 
  accu = 0.d0
  do k = 1, idet  
   do j = 1, idet
    call i_H_j(psi_1(1,1,k),psi_1(1,1,j),N_int,hij) 
    accu += psi_1_coef(k) * psi_1_coef(j) * hij
   enddo
  enddo
  diagonal_matrix_superci_brutal(i) = accu - psi_energy(1)
  total_matrix_superci_brutal(i,i) = accu - psi_energy(1)
  ! < Psi_0 | H | Psi_1> 
  accu = 0.d0
  do k = 1, N_det
   do j = 1, idet
    call i_H_j(psi_det(1,1,k),psi_1(1,1,j),N_int,hij) 
    accu += psi_1_coef(j) * psi_coef(k,1) * hij
   enddo
  enddo
  brillouin_matrix_superci_brutal(i) = accu 
  total_matrix_superci_brutal(i,1) = accu 
  total_matrix_superci_brutal(1,i) = accu 

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
      psi_2_coef(jdet) = contraction_coef_psi_superci(j,1,l)
    enddo
    do j = 1, N_det
     if(psi_superci(1,1,j,2,i) == 0_bit_kind)cycle
      jdet+=1
      do k = 1, N_int
       psi_2(k,1,jdet) = psi_superci(k,1,j,2,l)
       psi_2(k,2,jdet) = psi_superci(k,2,j,2,l)
      enddo
      psi_2_coef(jdet) = contraction_coef_psi_superci(j,2,l)
    enddo
    ! < Psi_1 | H | Psi_1> 
    accu = 0.d0
    do k = 1, idet  
     do j = 1, jdet
      call i_H_j(psi_1(1,1,k),psi_2(1,1,j),N_int,hij) 
      accu += psi_1_coef(k) * psi_2_coef(j) * hij
     enddo
    enddo
    total_matrix_superci_brutal(i,l) = accu
    total_matrix_superci_brutal(l,i) = accu
    deallocate(psi_2,psi_2_coef)
    
  enddo

  deallocate(psi_1,psi_1_coef)
 enddo

END_PROVIDER 




 BEGIN_PROVIDER [double precision, eigenvectors_sci_brutal, (size_super_ci,size_super_ci)]
&BEGIN_PROVIDER [double precision, eigenvalues_sci_brutal, (size_super_ci)]
 implicit none
 integer :: i,j

 call lapack_diag(eigenvalues_sci_brutal,eigenvectors_sci_brutal,total_matrix_superci_brutal,size_super_ci,size_super_ci)


END_PROVIDER 
