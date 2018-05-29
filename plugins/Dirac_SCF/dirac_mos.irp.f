 BEGIN_PROVIDER [ complex*16, dirac_mo_coef, (2*(dirac_ao_num),2*(dirac_mo_tot_num))
 implicit none
  BEGIN_DOC
  ! Molecular orbital coefficients on AO basis set
  ! dirac_mo_coef(i,j) = coefficient of the ith ao on the jth mo
  ! mo_label : Label characterizing the MOS (local, canonical, natural,
  ! etc)
  END_DOC
  integer                        :: i,i_minus,j,j_minus
  PROVIDE ezfio_filename
  dirac_mo_coef = (0.d0,0.d0)
  do j=1, 2*dirac_mo_tot_num
   if (j .le. large_mo_tot_num) then
    do i=1, 2*dirac_ao_num
     if (i .le. large_ao_num) then
      dirac_mo_coef(i,j) = (1.d0,0.d0)*large_mo_coef(i,j)
     endif
    enddo
   elseif (j .gt. large_mo_tot_num .and. j .le. 2*large_mo_tot_num) then
    j_minus = j - large_mo_tot_num
    do i=1, 2*dirac_ao_num
     if (i .gt. large_ao_num .and. i .le. 2*large_ao_num) then
      i_minus = i - large_ao_num
      dirac_mo_coef(i,j) = (1.d0,0.d0)*large_mo_coef(i_minus,j_minus)
     endif
    enddo
   elseif (j.gt. 2*large_mo_tot_num .and. j .le. (2*large_mo_tot_num+small_mo_tot_num)) then
    j_minus = j - 2*large_mo_tot_num
    do i=1, 2*dirac_ao_num
    i_minus = i - 2*large_ao_num
     if (i .gt. 2*large_ao_num .and. i .le. (2*large_ao_num+small_ao_num)) then
      dirac_mo_coef(i,j) = (1.d0,0.d0)*small_mo_coef(i_minus,j_minus)
     endif
    enddo
   elseif (j .gt. (2*large_ao_num + small_ao_num)) then
    j_minus = j - (2*large_mo_tot_num + small_mo_tot_num)
    do i=1, 2*(large_ao_num+small_ao_num)
     if (i .gt. (2*large_ao_num+small_ao_num)) then
      i_minus = i - (2*large_ao_num + small_ao_num)
      dirac_mo_coef(i,j) = (1.d0,0.d0)*small_mo_coef(i_minus,j_minus)
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
  allocate ( T(2*(dirac_ao_num),2*(dirac_mo_tot_num)) )
  call zgemm('N','N', 2*(dirac_ao_num), 2*(dirac_mo_tot_num), 2*(dirac_ao_num), &
      (1.d0,0.d0), A_ao,LDA_ao,                                             &
      dirac_mo_coef, size(dirac_mo_coef,1),                          &
      (0.d0,0.d0), T, size(T,1))
  call zgemm('T','N', 2*(dirac_mo_tot_num), 2*(dirac_mo_tot_num), 2*(dirac_ao_num), &
      (1.d0,0.d0), dirac_mo_coef,size(dirac_mo_coef,1),                     &
      T, 2*(dirac_ao_num),                                    &
      (0.d0,0.d0), A_mo, size(A_mo,1))
  deallocate(T)
 end


 BEGIN_PROVIDER [double precision,dirac_fock_matrix_eigenvalues,(2*(dirac_mo_tot_num))]
 &BEGIN_PROVIDER [complex*16, dirac_fock_matrix_eigenvectors,(2*(dirac_mo_tot_num),2*(dirac_mo_tot_num))]
 implicit none
 integer :: n,nmax
 double precision :: eigenvalues( 2*(dirac_mo_tot_num))
 complex*16       :: eigenvectors(2*(dirac_mo_tot_num),2*(dirac_mo_tot_num))
 n = 2*(dirac_mo_tot_num)
 nmax = n
 call lapack_diag_complex(eigenvalues,eigenvectors,dirac_mo_mono_elec_integral,nmax,n)
 dirac_fock_matrix_eigenvalues = eigenvalues
 dirac_fock_matrix_eigenvectors = eigenvectors
 END_PROVIDER
 
