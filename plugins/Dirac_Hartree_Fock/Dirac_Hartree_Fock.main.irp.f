program Dirac_Hartree_Fock
  implicit none

  integer :: i,j
  double precision :: pi
  complex*16 :: ortho(2*dirac_ao_num)
 pi = 4d0*atan(1d0)

 print*, "small_ao_num =",small_ao_num
 print*, "large_ao_num =",large_ao_num

!do i=1,2*dirac_mo_tot_num
! print*,i,eigenvalues_dirac_fock_matrix_C_G_mo(i)
!enddo



 !! Print the positive energy molecular orbitals
!do j= 2*small_ao_num+1, 2*small_ao_num + elec_num
! do i= 1, 2*dirac_ao_num
!  if (i .le. large_ao_num) then
!   print*, i,j, dirac_mo_coef(i,j), dirac_ao_nucl(i),dirac_ao_power(i,1), dirac_ao_power(i,2), dirac_ao_power(i,3) 
!  elseif (i .gt. large_ao_num .and. i .le. 2*large_ao_num) then
!   print*, i,j, dirac_mo_coef(i,j), dirac_ao_nucl(i-large_ao_num),dirac_ao_power(i-large_ao_num,1), dirac_ao_power(i-large_ao_num,2), dirac_ao_power(i-large_ao_num,3)
!  elseif (i .gt. 2*large_ao_num .and. i .le. 2*large_ao_num + small_ao_num) then
!   print*, i,j, dirac_mo_coef(i,j), dirac_ao_nucl(i-large_ao_num),dirac_ao_power(i-large_ao_num,1),dirac_ao_power(i-large_ao_num,2), dirac_ao_power(i-large_ao_num,3)
!  elseif (i .gt. large_ao_num + small_ao_num) then
!   print*, i,j, dirac_mo_coef(i,j), dirac_ao_nucl(i-large_ao_num-small_ao_num),dirac_ao_power(i-large_ao_num-small_ao_num,1), dirac_ao_power(i-large_ao_num-small_ao_num,2), dirac_ao_power(i-large_ao_num-small_ao_num,3)
!  endif
! enddo
! print*,'************'
!enddo

 print*,'************************************************************************'
 do j =1,2*dirac_ao_num
  do i = 1,2*dirac_ao_num 
   print*, i, j, dirac_SCF_density_matrix_ao(i,j)
  enddo
  print*, '***************'
 enddo
 

!!! Test to check if the Kramer pair symetries are indeed present 
!!!(Does not work for the orbitals corresponding to P-1 or P+1, and most likely not for D-1 and D+1) 
!do j= 2*small_ao_num+1, 2*small_ao_num + 2*large_ao_num
! if ( modulo(j,2) == 1) then 
! do i= 1, large_ao_num
!  if (Abs(dirac_mo_coef(i,j)) .lt. 1.e-9) then
!   print*,i,j, Abs(dirac_mo_coef(i,j))
!  else
!   print*,i,j, Abs(dirac_mo_coef(i,j))/Abs(dirac_mo_coef(i+large_ao_num,j+1)), atan(aimag(dirac_mo_coef(i,j))/real(dirac_mo_coef(i,j))) -atan(-aimag(dirac_mo_coef(i+large_ao_num,j+1))/real(dirac_mo_coef(i+large_ao_num,j+1))), dirac_ao_power(i,1), dirac_ao_power(i,2), dirac_ao_power(i,3)
!  endif 
! enddo
! print*,'*********************'
! do i = large_ao_num+1, 2*large_ao_num
!  if (Abs(dirac_mo_coef(i,j)) .lt. 1.e-9) then
!   print*,i,j, Abs(dirac_mo_coef(i,j))
!  else
!   print*,i,j, Abs(dirac_mo_coef(i,j))/Abs(dirac_mo_coef(i-large_ao_num,j+1)), atan(aimag(dirac_mo_coef(i,j))/real(dirac_mo_coef(i,j))) -atan(-aimag(dirac_mo_coef(i-large_ao_num,j+1))/real(dirac_mo_coef(i-large_ao_num,j+1))), dirac_ao_power(i-large_ao_num,1), dirac_ao_power(i-large_ao_num,2), dirac_ao_power(i-large_ao_num,3)  
!  endif 
! enddo
! print*,'*********************'
! do i = 2*large_ao_num+1, 2*large_ao_num + small_ao_num
!  if (Abs(dirac_mo_coef(i,j)) .lt. 1.e-9) then
!   print*,i,j, Abs(dirac_mo_coef(i,j))
!  else
!   print*,i,j, Abs(dirac_mo_coef(i,j))/Abs(dirac_mo_coef(i+small_ao_num,j+1)),atan(aimag(dirac_mo_coef(i,j))/real(dirac_mo_coef(i,j))) -atan(-aimag(dirac_mo_coef(i+small_ao_num,j+1))/real(dirac_mo_coef(i+small_ao_num,j+1))), dirac_ao_power(i-large_ao_num,1), dirac_ao_power(i-large_ao_num,2), dirac_ao_power(i-large_ao_num,3)
!  endif
! enddo
! print*,'*********************'
! do i = 2*large_ao_num + small_ao_num + 1, 2*large_ao_num + 2*small_ao_num
!  if (Abs(dirac_mo_coef(i,j)) .lt. 1.e-9) then
!   print*,i,j, Abs(dirac_mo_coef(i,j))
!  else
!   print*,i,j, Abs(dirac_mo_coef(i,j))/Abs(dirac_mo_coef(i-small_ao_num,j+1)),atan(aimag(dirac_mo_coef(i,j))/real(dirac_mo_coef(i,j))) -atan(-aimag(dirac_mo_coef(i-small_ao_num,j+1))/real(dirac_mo_coef(i-small_ao_num,j+1))), dirac_ao_power(i-dirac_ao_num,1), dirac_ao_power(i-dirac_ao_num,2), dirac_ao_power(i-dirac_ao_num,3)
!  endif
! enddo
!print*,'*****************************************************'
! endif 
!enddo


!!! Due to P1 = (Px + i Py)/Sqrt[2] and P-1 = (Px - i Py)/Sqrt[2]
!!! I couldn't make this test work, most likely the 4 degenerate P orbitals are combinations and cannot therefore be straightforwardly compared.
!do j= 2*small_ao_num+7, 2*small_ao_num +7
! do i= 1, large_ao_num
!  if (Abs(dirac_mo_coef(i,j)) .lt. 1.e-9) then
!   print*,i,j, Abs(dirac_mo_coef(i,j))
!  else
!   print*,i,j, Abs(dirac_mo_coef(i,j))/Abs(dirac_mo_coef(i+large_ao_num,j+3)), atan(aimag(dirac_mo_coef(i,j))/real(dirac_mo_coef(i,j))) -atan(-aimag(dirac_mo_coef(i+large_ao_num,j+3))/real(dirac_mo_coef(i+large_ao_num,j+3))), dirac_ao_power(i,1), dirac_ao_power(i,2), dirac_ao_power(i,3)
!  endif 
! enddo
!print*,'*****************************************************'
!enddo
end
