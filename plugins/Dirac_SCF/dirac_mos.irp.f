 BEGIN_PROVIDER [ complex*16,dirac_mo_coef,(2*(ao_num+small_ao_num),2*(mo_tot_num+small_mo_tot_num))
 implicit none
  BEGIN_DOC
  ! Molecular orbital coefficients on AO basis set
  ! dirac_mo_coef(i,j) = coefficient of the ith ao on the jth mo
  ! mo_label : Label characterizing the MOS (local, canonical, natural,
  ! etc)
  END_DOC
  integer                        :: i, j, k, l
  double precision, allocatable  :: buffer(:,:)
  logical                        :: exists
  PROVIDE ezfio_filename
  do i=1, 2*(mo_tot_num+small_mo_tot_num)
   if (i .le. mo_tot_num) then
    l = i - 0
    do j=1, 2*(ao_num+small_ao_num)
     if (j .le. ao_num) then
      k = j - 0
      dirac_mo_coef(j,i) = (1.d0,0.d0)*mo_coef(k,l)
     elseif (j .gt. mo_tot_num) then
      dirac_mo_coef(j,i) = (1.d0,0.d0)*0.d0
     endif
    enddo
   elseif (i.gt. mo_tot_num .and. i .le. 2*mo_tot_num) then
    l = i - mo_tot_num
    do j=1, 2*(ao_num+small_ao_num)
     if (j .le. ao_num) then
      dirac_mo_coef(j,i) = (1.d0,0.d0)*0.d0
     elseif (j .gt. ao_num .and. j .le. 2*ao_num) then
      k = j - ao_num
      dirac_mo_coef(j,i) = (1.d0,0.d0)*mo_coef(k,l)
     elseif (j .gt. 2*ao_num) then
      dirac_mo_coef(j,i) = (1.d0,0.d0)*0.d0
     endif
    enddo
   elseif (i.gt. 2*mo_tot_num .and. i .le. (2*mo_tot_num+small_mo_tot_num)) then
    l = i - 2*mo_tot_num
    do j=1, 2*(ao_num+small_ao_num)
     if (j .le. 2*ao_num) then
      dirac_mo_coef(j,i) = (1.d0,0.d0)*0.d0
     elseif (j .gt. 2*ao_num .and. j .le. (2*ao_num+small_ao_num)) then
      k = j - 2*ao_num
      dirac_mo_coef(j,i) = (1.d0,0.d0)*small_mo_coef(k,l)
     elseif (j .gt. (2*ao_num+small_ao_num)) then
      dirac_mo_coef(j,i) = (1.d0,0.d0)*0.d0
     endif
    enddo
   else
    l = i - (2*mo_tot_num+small_mo_tot_num)
    do j=1, 2*(ao_num+small_ao_num)
     if (j .le. (2*ao_num+small_ao_num)) then
      dirac_mo_coef(j,i) = (1.d0,0.d0)*0.d0
     else
      k = j - (2*ao_num+small_ao_num)
      dirac_mo_coef(j,i) = (1.d0,0.d0)*small_mo_coef(k,l)
     endif
    enddo
   endif
  enddo
 END_PROVIDER


 subroutine dirac_ao_to_mo(A_ao,LDA_ao,A_mo,LDA_mo)
  implicit none
  BEGIN_DOC
  ! Transform A from the AO basis to the MO basis
  !
  ! Ct.A_ao.C
  END_DOC
  integer, intent(in)            :: LDA_ao,LDA_mo
  complex*16, intent(in)         :: A_ao(LDA_ao,2*(ao_num+small_ao_num))
  double precision, intent(out)  :: A_mo(LDA_mo,2*(mo_tot_num+small_mo_tot_num))
  complex*16, allocatable        :: T(:,:)
  allocate ( T(2*(ao_num+small_ao_num),2*(mo_tot_num+small_mo_tot_num)) )
  !DIR$ ATTRIBUTES ALIGN : $IRP_ALIGN :: T
  call zgemm('N','N', 2*(ao_num+small_ao_num), 2*(mo_tot_num+small_mo_tot_num), 2*(ao_num+small_ao_num), &
      (1.d0,0.d0), A_ao,LDA_ao,                                             &
      dirac_mo_coef, size(dirac_mo_coef,1),                          &
      (0.d0,0.d0), T, size(T,1))
  call zgemm('T','N', 2*(mo_tot_num+small_mo_tot_num), 2*(mo_tot_num+small_mo_tot_num), 2*(ao_num+small_ao_num), &
      (1.d0,0.d0), dirac_mo_coef,size(dirac_mo_coef,1),                     &
      T, 2*(ao_num+small_ao_num),                                    &
      (0.d0,0.d0), A_mo, size(A_mo,1))
  deallocate(T)
 end


 BEGIN_PROVIDER [double precision,dirac_fock_matrix_eigenvalues,(2*(mo_tot_num+small_mo_tot_num))]
 &BEGIN_PROVIDER [complex*16, dirac_fock_matrix_eigenvectors,(2*(mo_tot_num+small_mo_tot_num),2*(mo_tot_num+small_mo_tot_num))]
 implicit none
 integer :: n,nmax
 double precision :: eigenvalues( 2*(mo_tot_num+small_mo_tot_num))
 complex*16       :: eigenvectors(2*(mo_tot_num+small_mo_tot_num),2*(mo_tot_num+small_mo_tot_num))
 n = 2*(mo_tot_num+small_mo_tot_num)
 nmax = n
 call lapack_diag_complex(eigenvalues,eigenvectors,dirac_mo_mono_elec_integral,nmax,n)
 dirac_fock_matrix_eigenvalues = eigenvalues
 dirac_fock_matrix_eigenvectors = eigenvectors
 END_PROVIDER
 
