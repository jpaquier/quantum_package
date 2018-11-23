
 BEGIN_PROVIDER [ complex*16, dirac_mo_coef_S, (2*(dirac_ao_num),2*(dirac_mo_tot_num))
 implicit none
  BEGIN_DOC
  !Molecular orbital coefficients on AO basis set diagonalizing the overlap matrix S
  !dirac_mo_coef_S(i,j) = coefficient of the ith ao on the jth mo
  END_DOC
  integer                        :: i,i_minus,j,j_minus
  PROVIDE ezfio_filename
  dirac_mo_coef_S = (0.d0,0.d0)
  do j=1, 2*dirac_mo_tot_num
   if (j .le. large_mo_tot_num) then
    do i=1, 2*dirac_ao_num
     if (i .le. large_ao_num) then
      dirac_mo_coef_S(i,j) = (1.d0,0.d0)*large_mo_coef(i,j)
     endif
    enddo
   elseif (j .gt. large_mo_tot_num .and. j .le. 2*large_mo_tot_num) then
    j_minus = j - large_mo_tot_num
    do i=1, 2*dirac_ao_num
     if (i .gt. large_ao_num .and. i .le. 2*large_ao_num) then
      i_minus = i - large_ao_num
      dirac_mo_coef_S(i,j) = (1.d0,0.d0)*large_mo_coef(i_minus,j_minus)
     endif
    enddo
   elseif (j.gt. 2*large_mo_tot_num .and. j .le. (2*large_mo_tot_num+small_mo_tot_num)) then
    j_minus = j - 2*large_mo_tot_num
    do i=1, 2*dirac_ao_num
    i_minus = i - 2*large_ao_num
     if (i .gt. 2*large_ao_num .and. i .le. (2*large_ao_num+small_ao_num)) then
      dirac_mo_coef_S(i,j) = (1.d0,0.d0)*small_mo_coef(i_minus,j_minus)
     endif
    enddo
   elseif (j .gt. (2*large_ao_num + small_ao_num)) then
    j_minus = j - (2*large_mo_tot_num + small_mo_tot_num)
    do i=1, 2*(large_ao_num+small_ao_num)
     if (i .gt. (2*large_ao_num+small_ao_num)) then
      i_minus = i - (2*large_ao_num + small_ao_num)
      dirac_mo_coef_S(i,j) = (1.d0,0.d0)*small_mo_coef(i_minus,j_minus)
     endif
    enddo
   endif
  enddo
 END_PROVIDER
 
