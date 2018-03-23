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
 double precision :: integral,get_mo_bielec_integral,mo_bielec_integral
 double precision :: two_dm_orb
 threshold = 0.d-10

 call give_all_mos_at_r(r1,mos_array_r1) 
 call give_all_mos_at_r(r2,mos_array_r2) 

 two_dm = 0.d0
 coulomb = 0.d0
 do i = 1, mo_tot_num
  a1 = mos_array_r1(i)  ! electron 1 
  c1 = dabs(a1)
  if(c1.le.threshold)cycle
  do j = 1, mo_tot_num 
   a2 = a1 * mos_array_r2(j)  ! electron 2
   c2 = dabs(a2)
   if(c2.le.threshold)cycle
!  call get_mo_bielec_integrals_ij(i,j,mo_tot_num,integrals_array,mo_integrals_map) 
   do m = 1, mo_tot_num
    a3 = a2 * mos_array_r1(m)  ! electron 1 
    c3 = dabs(a3)
    if(c3.le.threshold)cycle
    do n = 1, mo_tot_num
     a4 = a3 * mos_array_r2(n)   ! electron 2
     c4 = dabs(a4)
!    two_dm += two_bod_alpha_beta_mo_transposed(n,m,j,i,1) * mos_array_r1(i) * mos_array_r1(n) * mos_array_r2(j) * mos_array_r2(m)
     if(c4.le.threshold)cycle
       do l = 1, mo_tot_num
        do k = 1, mo_tot_num
        integral = mo_bielec_integral(i,j,k,l)
!        coulomb += a4 * integrals_array(k,l) * two_bod_alpha_beta_mo_transposed(k,l,n,m,1)
         coulomb += integral * two_bod_alpha_beta_mo_transposed(k,l,n,m,1) * mos_array_r2(n) * mos_array_r2(j) * mos_array_r1(i) * mos_array_r1(m)
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
!coulomb = coulomb / two_dm
!coulomb = exp(log(coulomb) - log(two_dm))

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
     two_dm += two_bod_alpha_beta_mo_transposed(n,m,j,i,1) * mos_array_r1(i) * mos_array_r1(n) * mos_array_r2(j) * mos_array_r2(m)
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


subroutine numerical_delta_function_coulomb(r1,r2,integral)
 implicit none
 BEGIN_DOC
! compute numerically the following integral: 
! int_r \sum_m \phi_m(r)\phi_m(r2) * 1/|r - r1|
 END_DOC
 double precision,intent(in) :: r1(3),r2(3)
 double precision,intent(out):: integral
 double precision :: r(3),mos_array_r2(mo_tot_num),mos_array_r(mo_tot_num),distance
 integer :: j,k,l ,m
 integral = 0.d0
 call give_all_mos_at_r(r2,mos_array_r2) 
 do j = 1, nucl_num
  do k = 1, n_points_radial_grid  -1
   do l = 1, n_points_integration_angular 
    r(1) = grid_points_per_atom(1,l,k,j)
    r(2) = grid_points_per_atom(2,l,k,j)
    r(3) = grid_points_per_atom(3,l,k,j)
    distance  = (r1(1) - r(1))*(r1(1) - r(1))
    distance += (r1(2) - r(2))*(r1(2) - r(2))
    distance += (r1(3) - r(3))*(r1(3) - r(3))
    distance = dsqrt(distance)
    if(dabs(distance).lt.1.d-10)cycle
    call give_all_mos_at_r(r,mos_array_r) 
    do m = 1, mo_tot_num
     integral += mos_array_r(m) * mos_array_r2(m) * 1.d0/(distance) * final_weight_functions_at_grid_points(l,k,j)
    enddo
   enddo
  enddo
 enddo
end

subroutine local_r12_operator_on_hf(r1,r2,integral_psi)
 implicit none
 BEGIN_DOC
! computes the following ANALYTICAL integral
! \sum_{k,l} phi_l(r1) phi_k(r2) int_{r,r'} phi_k(x) phi_l(x') 1s(x) 1s(x') * 1 /|x-x'|
 END_DOC
 double precision, intent(in) :: r1(3), r2(3)
 double precision, intent(out):: integral_psi
 integer :: k,l
 double precision :: mo_bielec_integral
 double precision :: mos_array_r1(mo_tot_num)
 double precision :: mos_array_r2(mo_tot_num)
 double precision :: get_mo_bielec_integral
 call give_all_mos_at_r(r1,mos_array_r1) 
 call give_all_mos_at_r(r2,mos_array_r2) 

 integral_psi = 0.d0
 do k = 1, mo_tot_num
  do l = 1, mo_tot_num
   integral_psi += get_mo_bielec_integral(1,1,k,l,mo_integrals_map) * mos_array_r2(l) * mos_array_r1(k)
  enddo
 enddo
 integral_psi = integral_psi / (mos_array_r1(1) * mos_array_r2(1))

end


subroutine local_erf_r12_operator_on_hf(r1,r2,integral_psi)
 implicit none
 BEGIN_DOC
! computes the following ANALYTICAL integral
! \sum_{k,l} phi_l(r1) phi_k(r2) int_{r,r'} phi_k(r) phi_l(r') 1s(r) 1s(r') * erf(mu_erf |r-r'|) /|r-r'|
 END_DOC
 double precision, intent(in) :: r1(3), r2(3)
 double precision, intent(out):: integral_psi
 integer :: k,l
 double precision :: mo_bielec_integral
 double precision :: mos_array_r1(mo_tot_num)
 double precision :: mos_array_r2(mo_tot_num)
 double precision :: get_mo_bielec_integral_erf
 call give_all_mos_at_r(r1,mos_array_r1) 
 call give_all_mos_at_r(r2,mos_array_r2) 

 integral_psi = 0.d0
 do k = 1, mo_tot_num
  do l = 1, mo_tot_num
   integral_psi += integrals_for_hf_potential(l,k) * mos_array_r2(l) * mos_array_r1(k)
  enddo
 enddo
 integral_psi = integral_psi / (mos_array_r1(1) * mos_array_r2(1))

end

BEGIN_PROVIDER [double precision, integrals_for_hf_potential, (mo_tot_num,mo_tot_num)]
 implicit none
 integer :: k,l
 double precision :: get_mo_bielec_integral_erf
 do k = 1, mo_tot_num
  do l = 1, mo_tot_num
   integrals_for_hf_potential(l,k) = get_mo_bielec_integral_erf(1,1,k,l,mo_integrals_erf_map) 
  enddo
 enddo


END_PROVIDER 

double precision function mu_coulomb(y,x)
 implicit none
 BEGIN_DOC
! finds the correct mu that such that 
! y = erf(mu_coulomb * x)/x 
 ENd_DOC
 double precision, intent(in) :: y,x
 double precision :: y_tmp,x_tmp,thr,inverse_erf
 thr = 1.d-10
 y_tmp = y * x
 x_tmp =inverse_erf(y_tmp,thr)
 if(x.ne.0.d0)then
  mu_coulomb = x_tmp / x
 else
  mu_coulomb = x_tmp / 1.d-15
 endif
end
