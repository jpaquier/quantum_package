 subroutine give_all_dirac_aos_at_r(r,aos_array)
 implicit none
 BEGIN_dOC
 ! input : r == r(1) = x and so on
 ! aos_array(i) = aos(i) evaluated in r
 END_DOC
 double precision, intent(in) :: r(3)
 double precision, intent(out) :: aos_array(dirac_ao_num)
 integer :: power_ao(3) 
 integer :: i,j,k,l,m
 double precision :: dx,dy,dz,r2
 double precision ::      dx2,dy2,dz2
 double precision :: center_ao(3)
 double precision :: beta
  do i = 1, nucl_num
   center_ao(1:3) = nucl_coord(i,1:3)
   dx = (r(1) - center_ao(1)) 
   dy = (r(2) - center_ao(2)) 
   dz = (r(3) - center_ao(3)) 
   r2 = dx*dx + dy*dy + dz*dz
   do j = 1,Nucl_N_dirac_Aos(i) 
    k = Nucl_dirac_Aos_transposed(j,i) ! index of the dirac ao in the ordered format 
    aos_array(k) = 0.d0
    power_ao(1:3)= dirac_ao_power_ordered_transp_per_nucl(1:3,j,i)
    double precision :: power
    dx2 = power(power_ao(1),dx)
    dy2 = power(power_ao(2),dy)
    dz2 = power(power_ao(3),dz)
    do l = 1,dirac_ao_prim_num(k) 
     beta = dirac_ao_expo_ordered_transp_per_nucl(l,j,i)
     aos_array(k)+= dirac_ao_coef_normalized_ordered_transp_per_nucl(l,j,i) * dexp(-beta*r2) 
    enddo
    aos_array(k) = aos_array(k) * dx2 * dy2 * dz2
   enddo
  enddo
 end

