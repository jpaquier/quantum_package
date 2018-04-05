double precision function two_dm_in_r(r1,r2)
 implicit none
 double precision, intent(in) :: r1(3),r2(3)
 double precision :: mos_array_r1(mo_tot_num)
 double precision :: mos_array_r2(mo_tot_num)
 integer :: i,j,k,l
 call give_all_mos_at_r(r1,mos_array_r1) 
 call give_all_mos_at_r(r2,mos_array_r2) 
 two_dm_in_r = 0.d0
 do i = 1, mo_tot_num
  do j = 1, mo_tot_num
   do k = 1, mo_tot_num
    do l = 1, mo_tot_num
     two_dm_in_r += two_bod_alpha_beta_mo_transposed(l,k,j,i,1) * mos_array_r1(i) * mos_array_r1(l)  * mos_array_r2(k) * mos_array_r2(j)
    enddo
   enddo
  enddo
 enddo
end
