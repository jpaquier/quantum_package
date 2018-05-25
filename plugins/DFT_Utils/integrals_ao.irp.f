BEGIN_PROVIDER [double precision, integrals_ao_over_r, (ao_num)]
 implicit none
 BEGIN_DOC
!  integral over r of the AOs; note that it is not the integral of the AO squared !!!!
 END_DOC
 integer :: i,j,k
 integer :: power_A(3)
 double precision :: F_integral,integral(3),integral_tot,alpha
 integrals_ao_over_r = 0.d0 
 do i = 1, ao_num
  power_A(1:3)= ao_power(i,1:3)
  do j = 1, ao_prim_num(i)  
   alpha = ao_expo_ordered_transp(j,i)
   do k = 1, 3
    integral(k) = F_integral(power_A(k),alpha)
   enddo
   integral_tot = integral(1) * integral(2) * integral(3)
   integrals_ao_over_r(i) +=  ao_coef_normalized_ordered_transp(j,i) * integral_tot
  enddo
 enddo

END_PROVIDER 

BEGIN_PROVIDER [double precision, integrals_mo_over_r, (mo_tot_num)]
 implicit none
 BEGIN_DOC
!  integral over r of the MOs; note that it is not the integral of the MO squared !!!!
 END_DOC
 integrals_mo_over_r = 0.d0
 integer :: i,j
 do i = 1, mo_tot_num
  do j = 1, ao_num
   integrals_mo_over_r(i) += integrals_ao_over_r(j) * mo_coef(j,i)
  enddo
 enddo



END_PROVIDER 
