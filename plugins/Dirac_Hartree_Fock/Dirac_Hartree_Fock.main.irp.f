program Dirac_Hartree_Fock
  implicit none
  integer :: i,j
  double precision :: pi
  complex*16 :: ortho(2*dirac_ao_num)

 print*, "small_ao_num =",small_ao_num
 print*, "large_ao_num =",large_ao_num

!do i = 1,large_ao_num
! print*,i, large_ao_expo_ordered_transp(1,i)
!enddo
!do i = 1,small_ao_num
! print*,i, small_ao_expo_ordered_transp(1,i)
!enddo
!print*,'***'
!do i = 1,dirac_ao_num
! print*,i, dirac_ao_expo_ordered_transp(1,i), dirac_ao_power(i,1), dirac_ao_power(i,2), dirac_ao_power(i,3)
!enddo
 do i=1,2*dirac_mo_tot_num
  print*,i,eigenvalues_dirac_fock_matrix_C_G_mo(i)
 enddo
 do j= 2*small_ao_num+1, 2*small_ao_num + 2*large_ao_num
  do i= 1, large_ao_num
   print*, i,j, dirac_mo_coef(i,j)
  enddo
  print*,'************'
 enddo

 pi = 4d0*atan(1d0)
 do j= 2*small_ao_num+1, 2*small_ao_num + 2*large_ao_num
  if ( modulo(j,2) == 1) then 
  do i= 1, large_ao_num
   if (Abs(dirac_mo_coef(i,j)) .lt. 10**(-5)) then
    print*,i,j, 0,0
   else
    print*,i,j, Abs(dirac_mo_coef(i,j))/Abs(dirac_mo_coef(i+large_ao_num,j+1)), atan(aimag(dirac_mo_coef(i,j))/real(dirac_mo_coef(i,j))) -atan(-aimag(dirac_mo_coef(i+large_ao_num,j+1))/real(dirac_mo_coef(i+large_ao_num,j+1)))
   endif 
  enddo
 print*,'*********************'
  do i = large_ao_num+1, 2*large_ao_num
   print*,i,j, Abs(dirac_mo_coef(i,j))/Abs(dirac_mo_coef(i-large_ao_num,j+1)), atan(aimag(dirac_mo_coef(i,j))/real(dirac_mo_coef(i,j))) -atan(-aimag(dirac_mo_coef(i-large_ao_num,j+1))/real(dirac_mo_coef(i-large_ao_num,j+1)))  
  enddo
 print*,'*******************************************'
  endif 
 enddo

! ortho = (0.d0,0.d0)
! do j= 2*small_ao_num, 2*small_ao_num+3
!  do i= 1,2*dirac_ao_num
!  !print*, i, j, dirac_mo_coef_guess(i,j)
!  !print*,i,j,dirac_mo_coef(i,j)
!  !print*,i,j,dirac_mo_coef_S(i,j) 
!  !print*,i,j,dirac_mo_overlap(i,j) 
!  !print*,i,j,dirac_mo_overlap_bis(i,j) 
!  !print*,i,j,dirac_SCF_density_matrix_ao(i+large_ao_num,j)
!  !print*,i,j,dirac_SCF_density_matrix_ao(j,i+large_ao_num)
!  !print*,i,j,dirac_ao_bi_elec_integral(i,j)
!  !print*,i,j,dirac_ao_bi_elec_integral_qq(i,j)
!  !print*,i,j,dirac_ao_mono_elec_integral (i,j)
!  !print*,i,j,dirac_Fock_matrix_ao (i,j)
!   ortho(j) += Conjg(dirac_mo_coef(i,j))*dirac_mo_coef(i,j)
!  enddo
! print*,'********************'     
!  print*,j,ortho(j)
! enddo
end
