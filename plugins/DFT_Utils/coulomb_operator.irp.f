subroutine coulomb_operator_in_real_space(r1,r2,coulomb)
use map_module
 implicit none
 double precision, intent(in) :: r1(3), r2(3)
 double precision, intent(out):: coulomb
 integer :: i,j,k,l 
 double precision :: integrals_array(mo_tot_num,mo_tot_num)
 double precision :: mos_array_r1(mo_tot_num)
 double precision :: mos_array_r2(mo_tot_num)
 double precision :: a1,a2,a3
 double precision :: c1,c2,c3
 double precision :: threshold
 double precision :: integral,get_mo_bielec_integral
 threshold = 1.d-12

 call give_all_mos_at_r(r1,mos_array_r1) 
 call give_all_mos_at_r(r2,mos_array_r2) 

!call give_all_mos_at_r_old(r1,mos_array_r1) 
!call give_all_mos_at_r_old(r2,mos_array_r2) 

 coulomb = 0.d0
 do i = 1, mo_tot_num
  a1 = mos_array_r1(i) 
  c1 = dabs(a1)
  if(c1.le.threshold)cycle
  do j = 1, mo_tot_num 
   a2 = a1 * mos_array_r2(j)
   c2 = dabs(a2)
   coulomb += coulomb_integrals_contracted_with_overlap(j,i) * a2 
   enddo
  enddo
 !coulomb = coulomb * 0.5d0 
end

BEGIN_PROVIDER [double precision, coulomb_integrals_contracted_with_overlap, (mo_tot_num,mo_tot_num)]
 implicit none
 BEGIN_DOC
! coulomb_integrals_contracted_with_overlap(i,j) = \sum_{kl} <ij|kl> * (\int_r1 phi_k(r1)) * (\int_r2 phi_l(r2))
 END_DOC
 integer :: i,j,k,l
 double precision :: integrals_array(mo_tot_num,mo_tot_num)
 coulomb_integrals_contracted_with_overlap = 0.d0
 do i = 1, mo_tot_num
  do j = 1, mo_tot_num 
   call get_mo_bielec_integrals_ij(i,j,mo_tot_num,integrals_array,mo_integrals_map) 
    do k = 1, mo_tot_num
     do l = 1, mo_tot_num
      coulomb_integrals_contracted_with_overlap(j,i)+=  integrals_array(k,l) * integrals_mo_over_r(k) * integrals_mo_over_r(l)
     enddo
    enddo
   enddo
  enddo


END_PROVIDER 


subroutine erf_coulomb_operator_in_real_space(r1,r2,coulomb)
use map_module
 implicit none
 double precision, intent(in) :: r1(3), r2(3)
 double precision, intent(out):: coulomb
 integer :: i,j,k,l 
 double precision :: integrals_array(mo_tot_num,mo_tot_num)
 double precision :: mos_array_r1(mo_tot_num)
 double precision :: mos_array_r2(mo_tot_num)
 double precision :: a1,a2,a3
 double precision :: c1,c2,c3
 double precision :: threshold
 double precision :: integral,get_mo_bielec_integral_erf
 threshold = 1.d-12

 call give_all_mos_at_r(r1,mos_array_r1) 
 call give_all_mos_at_r(r2,mos_array_r2) 

!call give_all_mos_at_r_old(r1,mos_array_r1) 
!call give_all_mos_at_r_old(r2,mos_array_r2) 

 coulomb = 0.d0
 do i = 1, mo_tot_num
  a1 = mos_array_r1(i) 
  c1 = dabs(a1)
  if(c1.le.threshold)cycle
  do j = 1, mo_tot_num 
   a2 = a1 * mos_array_r2(j)
   c2 = dabs(a2)
   if(c2.le.threshold)cycle
    do k = 1, mo_tot_num
     a3 = a2 * mos_array_r1(k)
     c3 = dabs(a3)
     if(c3.le.threshold)cycle
     do l = 1, mo_tot_num
      integral = get_mo_bielec_integral_erf(i,j,k,l,mo_integrals_erf_map)
      coulomb += a3 * mos_array_r2(l) * integral
!     coulomb += mos_array_r1(i) * mos_array_r1(k) * mos_array_r2(j) * mos_array_r2(l) * integrals_array(k,l)
     enddo
    enddo
   enddo
  enddo
 coulomb = coulomb * 0.5d0 

end

subroutine nuclear_coulomb_operator_in_real_space(r1,coulomb)
 implicit none
 double precision, intent(in) :: r1(3)
 double precision, intent(out):: coulomb
 integer :: i,j 
 double precision :: mos_array_r1(mo_tot_num)

 call give_all_mos_at_r(r1,mos_array_r1) 
 coulomb = 0.d0
 do i = 1, mo_tot_num
  do j = 1, mo_tot_num 
   coulomb += mo_nucl_elec_integral(j,i) * mos_array_r1(j) * integrals_mo_over_r(i)
  enddo
 enddo

end


subroutine expectation_value_in_real_space_old(r1,r2,coulomb,two_dm)
use map_module
 implicit none
 double precision, intent(in) :: r1(3), r2(3)
 double precision, intent(out):: coulomb, two_dm
 integer :: i,j,k,l,m,n  
 double precision :: integrals_array(mo_tot_num,mo_tot_num)
 double precision :: mos_array_r1(mo_tot_num)
 double precision :: mos_array_r2(mo_tot_num)
 double precision :: a1,a2,a3,a4
 double precision :: c1,c2,c3,c4
 double precision :: threshold
 double precision :: integral,get_mo_bielec_integral
 double precision :: two_dm_orb
 threshold = 0.d-10

 call give_all_mos_at_r(r1,mos_array_r1) 
 call give_all_mos_at_r(r2,mos_array_r2) 

 two_dm = 0.d0
 coulomb = 0.d0
 do i = 1, mo_tot_num
  a1 = mos_array_r1(i) 
  c1 = dabs(a1)
  if(c1.le.threshold)cycle
  do j = 1, mo_tot_num 
   a2 = a1 * mos_array_r2(j)
   c2 = dabs(a2)
   if(c2.le.threshold)cycle
   call get_mo_bielec_integrals_ij(i,j,mo_tot_num,integrals_array,mo_integrals_map) 
   do m = 1, mo_tot_num
    a3 = a2 * mos_array_r1(m)
    c3 = dabs(a3)
    if(c3.le.threshold)cycle
    do n = 1, mo_tot_num
     a4 = a3 * mos_array_r2(n)
     c4 = dabs(a4)
!    two_dm += two_bod_alpha_beta_mo_transposed(n,m,j,i,1) * mos_array_r1(i) * mos_array_r1(n) * mos_array_r2(j) * mos_array_r2(m)
     if(c4.le.threshold)cycle
       do l = 1, mo_tot_num
        do k = 1, mo_tot_num
         coulomb += a4 * integrals_array(k,l) * two_bod_alpha_beta_mo_transposed(k,l,n,m,1)
        !                               <ij|kl>    * gamma(kl,mn)                              
        enddo
       enddo
    enddo
   enddo
  enddo
 enddo
 do i = 1, mo_tot_num
  do j = 1, mo_tot_num 
   do k = 1, mo_tot_num
    do l = 1, mo_tot_num
     two_dm += two_bod_alpha_beta_mo_transposed(l,k,j,i,1) * mos_array_r1(i) * mos_array_r1(l)  * mos_array_r2(k) * mos_array_r2(j)  
    enddo
   enddo
  enddo
 enddo
!coulomb = coulomb * 0.5d0 
!print*,coulomb,two_dm
 coulomb = coulomb / two_dm

end

subroutine expectation_value_in_real_space(r1,r2,coulomb)
use map_module
 implicit none
 double precision, intent(in) :: r1(3), r2(3)
 double precision, intent(out):: coulomb
 integer :: i,j,k,l,m,n  
 double precision :: integrals_array(mo_tot_num,mo_tot_num)
 double precision :: mos_array_r1(mo_tot_num)
 double precision :: mos_array_r2(mo_tot_num)
 double precision :: a1,a2,a3,a4
 double precision :: c1,c2,c3,c4
 double precision :: threshold
 double precision :: integral,get_mo_bielec_integral
 double precision :: two_dm
 threshold = 1.d-12

 call give_all_mos_at_r(r1,mos_array_r1) 
 call give_all_mos_at_r(r2,mos_array_r2) 

 two_dm = 0.d0
 coulomb = 0.d0
 do i = 1, mo_tot_num
  a1 = mos_array_r1(i) 
  c1 = dabs(a1)
  if(c1.le.threshold)cycle
  do j = 1, mo_tot_num 
   a2 = a1 * mos_array_r2(j)
   c2 = dabs(a2)
   if(c2.le.threshold)cycle
   do m = 1, mo_tot_num
    a3 = a2 * mos_array_r1(m)
    c3 = dabs(a3)
    if(c3.le.threshold)cycle
    do n = 1, mo_tot_num
     a4 = a3 * mos_array_r2(n)
     c4 = dabs(a4)
     two_dm += two_bod_alpha_beta_mo_transposed(n,m,j,i,1) * a4 
     coulomb += a4 * two_bod_alpha_beta_mo_contracted(n,m,j,i,1)
    !if(dabs(a4 * two_bod_alpha_beta_mo_contracted(n,m,j,i,1)).gt.1.d-10)then
    ! print*,'#######'
    ! print*,n,m,j,i
    ! print*,coulomb
    ! print*,a4,two_bod_alpha_beta_mo_contracted(n,m,j,i,1),a4 * two_bod_alpha_beta_mo_contracted(n,m,j,i,1)
    ! print*,two_dm,coulomb/two_dm
    !endif
    enddo
   enddo
  enddo
 enddo
!if(two_dm.lt.0.d0)then
! print*,two_dm 
! print*,r1
! print*,r2
! pause
!endif

!if(coulomb.lt.0.d0)then
! print*,'coulomb'
! print*,coulomb
! print*,r1
! print*,r2
! pause
!endif
!coulomb = coulomb * 0.5d0 
!if(dabs(two_dm).lt.1.d-6.and.dabs(coulomb).lt.1.d-6)then
! coulomb = 0.d0 
!else
  coulomb = coulomb / two_dm
!endif

end
