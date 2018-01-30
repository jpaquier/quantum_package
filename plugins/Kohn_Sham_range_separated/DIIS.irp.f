BEGIN_PROVIDER [ double precision, threshold_DIIS_nonzero ]
 implicit none
 BEGIN_DOC
 ! If threshold_DIIS is zero, choose sqrt(thresh_scf)
 END_DOC
 if (threshold_DIIS == 0.d0) then
   threshold_DIIS_nonzero = dsqrt(thresh_scf)
 else
   threshold_DIIS_nonzero = threshold_DIIS
 endif
 ASSERT (threshold_DIIS_nonzero >= 0.d0)

END_PROVIDER

BEGIN_PROVIDER [double precision, FPS_SPF_Matrix_ao, (ao_num, ao_num)]
  implicit none
  BEGIN_DOC
  !   Commutator FPS - SPF
  END_DOC
  double precision, allocatable  :: scratch(:,:)
  allocate(                                                          &
      scratch(ao_num, ao_num)                                  &
      )
  
  ! Compute FP
  
  call dgemm('N','N',ao_num,ao_num,ao_num,                           &
      1.d0,                                                          &
      Fock_Matrix_ao,Size(Fock_Matrix_ao,1),                         &
      RS_KS_Density_Matrix_ao,Size(RS_KS_Density_Matrix_ao,1),             &
      0.d0,                                                          &
      scratch,Size(scratch,1))
  
  ! Compute FPS
  
  call dgemm('N','N',ao_num,ao_num,ao_num,                           &
      1.d0,                                                          &
      scratch,Size(scratch,1),                                       &
      ao_Overlap,Size(ao_Overlap,1),                                 &
      0.d0,                                                          &
      FPS_SPF_Matrix_ao,Size(FPS_SPF_Matrix_ao,1))
  
  ! Compute SP
  
  call dgemm('N','N',ao_num,ao_num,ao_num,                           &
      1.d0,                                                          &
      ao_Overlap,Size(ao_Overlap,1),                                 &
      RS_KS_Density_Matrix_ao,Size(RS_KS_Density_Matrix_ao,1),             &
      0.d0,                                                          &
      scratch,Size(scratch,1))
  
  ! Compute FPS - SPF
  
  call dgemm('N','N',ao_num,ao_num,ao_num,                           &
      -1.d0,                                                         &
      scratch,Size(scratch,1),                                       &
      Fock_Matrix_ao,Size(Fock_Matrix_ao,1),                         &
      1.d0,                                                          &
      FPS_SPF_Matrix_ao,Size(FPS_SPF_Matrix_ao,1))

END_PROVIDER

bEGIN_PROVIDER [double precision, FPS_SPF_Matrix_mo, (mo_tot_num, mo_tot_num)]
  implicit none
  begin_doc 
!   Commutator FPS - SPF in mo basis
  end_doc
  call ao_to_mo(FPS_SPF_Matrix_ao, size(FPS_SPF_Matrix_ao,1), &
     FPS_SPF_Matrix_mo, size(FPS_SPF_Matrix_mo,1))
END_PROVIDER


 BEGIN_PROVIDER [ double precision, eigenvalues_Fock_matrix_ao, (ao_num) ]
&BEGIN_PROVIDER [ double precision, eigenvectors_Fock_matrix_ao, (ao_num,ao_num) ]

   BEGIN_DOC
   ! Eigenvalues and eigenvectors of the Fock matrix over the ao basis
   END_DOC

   implicit none
   
   double precision, allocatable  :: scratch(:,:),work(:),Xt(:,:)
   integer                        :: lwork,info
   integer                        :: i,j
   
   lwork = 3*ao_num - 1
   allocate(                                                         &
       scratch(ao_num,ao_num),                                 &
       work(lwork),                                                  &
       Xt(ao_num,ao_num)                                             &
       )
 
! Calculate Xt

  do i=1,ao_num
    do j=1,ao_num
      Xt(i,j) = S_half_inv(j,i)
    enddo
  enddo

! Calculate Fock matrix in orthogonal basis: F' = Xt.F.X

  call dgemm('N','N',ao_num,ao_num,ao_num,     &
       1.d0,                                   &
       Fock_matrix_ao,size(Fock_matrix_ao,1),  &
       S_half_inv,size(S_half_inv,1),        &
       0.d0,                                   &
       eigenvectors_Fock_matrix_ao,size(eigenvectors_Fock_matrix_ao,1))       

  call dgemm('N','N',ao_num,ao_num,ao_num,                              &
       1.d0,                                                            &
       Xt,size(Xt,1),                                                   &
       eigenvectors_Fock_matrix_ao,size(eigenvectors_Fock_matrix_ao,1), &
       0.d0,                                                            &
       scratch,size(scratch,1))
     
! Diagonalize F' to obtain eigenvectors in orthogonal basis C' and eigenvalues
  
   call dsyev('V','U',ao_num,       &
        scratch,size(scratch,1),    &
        eigenvalues_Fock_matrix_ao, &
        work,lwork,info)
 
   if(info /= 0) then
     print *,  irp_here//' failed : ', info
     stop 1
   endif

! Back-transform eigenvectors: C =X.C'

  call dgemm('N','N',ao_num,ao_num,ao_num,     &
       1.d0,                                   &
       S_half_inv,size(S_half_inv,1),        &
       scratch,size(scratch,1),                &
       0.d0,                                   &
       eigenvectors_Fock_matrix_ao,size(eigenvectors_Fock_matrix_ao,1))       
   
END_PROVIDER

