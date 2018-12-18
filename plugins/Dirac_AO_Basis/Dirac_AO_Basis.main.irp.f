program Dirac_AO_Basis
  implicit none
  integer :: i,j
  double precision :: r(3)
  complex*16 :: dirac_aos_array(2*dirac_ao_num)
  integer :: dirac_power_ao(3)
  r= 0.d0
 do i=1, dirac_ao_num
 !print*, dirac_ao_nucl(i)
 enddo
 
 print*,large_ao_num,small_ao_num
 
 do j=1,nucl_num
  do i =1, Nucl_N_dirac_aos(j)
   print*,i,j,Nucl_dirac_Aos_transposed(i,j),dirac_ao_coef_normalized_ordered_transp_per_nucl(1,i,j),dirac_ao_expo_ordered_transp_per_nucl(1,i,j), dirac_ao_power_ordered_transp_per_nucl(1,i,j),dirac_ao_power_ordered_transp_per_nucl(2,i,j),dirac_ao_power_ordered_transp_per_nucl(3,i,j)
  enddo
 !print*,'******************'
 enddo

 call give_all_dirac_aos_at_r(r,dirac_aos_array)
!do i = 1,dirac_ao_num
! print*,i,dirac_aos_array_tmp(i)
!enddo
 print*,'**********************************************'
 do i = 1,2*dirac_ao_num
  print*,i,dirac_aos_array(i)
 enddo

end


