 BEGIN_PROVIDER [ complex*16, dirac_ao_ortho_canonical_coef, (2*(dirac_ao_num),2*(dirac_mo_tot_num))
 implicit none
  BEGIN_DOC
  ! Molecular orbital coefficients on AO basis set
  ! dirac_ao_ortho_canonical_coef(i,j) = coefficient of the ith ao on the jth mo
  ! ao_ortho_canonical_label : Label characterizing the MOS (local, canonical, natural,
  ! etc)
  END_DOC
  integer                        :: i,i_minus,j,j_minus
  PROVIDE ezfio_filename
  dirac_ao_ortho_canonical_coef = (0.d0,0.d0)
  do j=1, 2*dirac_mo_tot_num
   if (j .le. large_mo_tot_num) then
    do i=1, 2*dirac_ao_num
     if (i .le. large_ao_num) then
      dirac_ao_ortho_canonical_coef(i,j) = (1.d0,0.d0)*large_ao_ortho_canonical_coef(i,j)
     endif
    enddo
   elseif (j .gt. large_mo_tot_num .and. j .le. 2*large_mo_tot_num) then
    j_minus = j - large_mo_tot_num
    do i=1, 2*dirac_ao_num
     if (i .gt. large_ao_num .and. i .le. 2*large_ao_num) then
      i_minus = i - large_ao_num
      dirac_ao_ortho_canonical_coef(i,j) = (1.d0,0.d0)*large_ao_ortho_canonical_coef(i_minus,j_minus)
     endif
    enddo
   elseif (j.gt. 2*large_mo_tot_num .and. j .le. (2*large_mo_tot_num+small_mo_tot_num)) then
    j_minus = j - 2*large_mo_tot_num
    do i=1, 2*dirac_ao_num
    i_minus = i - 2*large_ao_num
     if (i .gt. 2*large_ao_num .and. i .le. (2*large_ao_num+small_ao_num)) then
      dirac_ao_ortho_canonical_coef(i,j) = (1.d0,0.d0)*small_ao_ortho_canonical_coef(i_minus,j_minus)
     endif
    enddo
   elseif (j .gt. (2*large_ao_num + small_ao_num)) then
    j_minus = j - (2*large_mo_tot_num + small_mo_tot_num)
    do i=1, 2*(large_ao_num+small_ao_num)
     if (i .gt. (2*large_ao_num+small_ao_num)) then
      i_minus = i - (2*large_ao_num + small_ao_num)
      dirac_ao_ortho_canonical_coef(i,j) = (1.d0,0.d0)*small_ao_ortho_canonical_coef(i_minus,j_minus)
     endif
    enddo
   endif
  enddo
 END_PROVIDER

 BEGIN_PROVIDER [ complex*16, dirac_ao_ortho_lowdin_coef, (2*(dirac_ao_num),2*(dirac_mo_tot_num))
 implicit none
  BEGIN_DOC
  ! Molecular orbital coefficients on AO basis set
  ! dirac_ao_ortho_lowdin_coef(i,j) = coefficient of the ith ao on the jth mo
  ! ao_ortho_lowdin_label : Label characterizing the MOS (local, canonical, natural,
  ! etc)
  END_DOC
  integer                        :: i,i_minus,j,j_minus
  PROVIDE ezfio_filename
  dirac_ao_ortho_lowdin_coef = (0.d0,0.d0)
  do j=1, 2*dirac_mo_tot_num
   if (j .le. large_mo_tot_num) then
    do i=1, 2*dirac_ao_num
     if (i .le. large_ao_num) then
      dirac_ao_ortho_lowdin_coef(i,j) = (1.d0,0.d0)*large_ao_ortho_lowdin_coef(i,j)
     endif
    enddo
   elseif (j .gt. large_mo_tot_num .and. j .le. 2*large_mo_tot_num) then
    j_minus = j - large_mo_tot_num
    do i=1, 2*dirac_ao_num
     if (i .gt. large_ao_num .and. i .le. 2*large_ao_num) then
      i_minus = i - large_ao_num
      dirac_ao_ortho_lowdin_coef(i,j) = (1.d0,0.d0)*large_ao_ortho_lowdin_coef(i_minus,j_minus)
     endif
    enddo
   elseif (j.gt. 2*large_mo_tot_num .and. j .le. (2*large_mo_tot_num+small_mo_tot_num)) then
    j_minus = j - 2*large_mo_tot_num
    do i=1, 2*dirac_ao_num
    i_minus = i - 2*large_ao_num
     if (i .gt. 2*large_ao_num .and. i .le. (2*large_ao_num+small_ao_num)) then
      dirac_ao_ortho_lowdin_coef(i,j) = (1.d0,0.d0)*small_ao_ortho_lowdin_coef(i_minus,j_minus)
     endif
    enddo
   elseif (j .gt. (2*large_ao_num + small_ao_num)) then
    j_minus = j - (2*large_mo_tot_num + small_mo_tot_num)
    do i=1, 2*(large_ao_num+small_ao_num)
     if (i .gt. (2*large_ao_num+small_ao_num)) then
      i_minus = i - (2*large_ao_num + small_ao_num)
      dirac_ao_ortho_lowdin_coef(i,j) = (1.d0,0.d0)*small_ao_ortho_lowdin_coef(i_minus,j_minus)
     endif
    enddo
   endif
  enddo
 END_PROVIDER


!subroutine dirac_huckel_guess
! implicit none
! BEGIN_DOC
! ! Build the dirac MOs using the extended Huckel model (This guess may have
! ! absolutely no meaning nor usage)
! END_DOC
! integer                        :: i,j,j_plus
! double precision               :: accu
! double precision               :: c
! character*(64)                 :: label
! complex*16, allocatable        :: A(:,:)
! label = "Guess"
! c = 0.5d0 * 1.75d0
! allocate (A(2*dirac_ao_num, 2*dirac_ao_num))
! A = (0.d0,0.d0)
! do j=1,2*dirac_ao_num
!  do i=1,2*dirac_ao_num
!   A(i,j) = c * dirac_ao_overlap(i,j) * (dirac_ao_mono_elec_integral_diag(i) + dirac_ao_mono_elec_integral_diag(j))
!  enddo
!  A(j,j) = dirac_ao_mono_elec_integral_diag(j) + dirac_ao_bi_elec_integral(j,j)
! enddo
! dirac_Fock_matrix_ao(1:2*dirac_ao_num,1:2*dirac_ao_num) = A(1:2*dirac_ao_num,1:2*dirac_ao_num)
! TOUCH dirac_Fock_matrix_ao
! do j =1, elec_num
!  j_plus = j+2*small_ao_num
!  do i = 1, 2*dirac_ao_num
!   dirac_mo_coef (i,j) = eigenvectors_dirac_fock_matrix_mo (i,j_plus)
!  enddo
! enddo
! SOFT_TOUCH dirac_mo_coef
!!call save_mos
! deallocate(A)
!end
