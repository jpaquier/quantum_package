 BEGIN_PROVIDER [ double precision, small_ao_nucl_elec_integral,(small_ao_num,small_ao_num)]
   BEGIN_DOC
   ! interaction nuclear electron for the AOs within the small component basis
   ! set in unrestricted kinetic balance
   END_DOC
   implicit none
   double precision               :: alpha, beta, gama, delta,c,Z
   integer                        :: num_A,num_B
   double precision               :: A_center(3),B_center(3),C_center(3)
   integer                        :: power_A(3),power_B(3)
   integer                        :: i,j,k,l,n_pt_in,m
   double precision               :: NAI_pol_mult 
! overlap_x,overlap_y,overlap_z,overlap,dx
!  if (read_ao_one_integrals) then
!   call read_one_e_integrals('ao_ne_integral', ao_nucl_elec_integral, size(ao_nucl_elec_integral,1), size(ao_nucl_elec_integral,2))
!    print *,  'AO N-e integrals read from disk'
!  else
   print*,' computing the small_ao_nucl_elec_integral'
     small_ao_nucl_elec_integral = 0.d0
     !        _
     ! /|  / |_)
     !  | /  | \
     !
   n_pt_in = n_pt_max_integrals
   do j = 1, small_ao_num
    num_A = small_ao_nucl(j)
    power_A(1:3)= small_ao_power(j,1:3)
    A_center(1:3) = nucl_coord(num_A,1:3)
    alpha = small_ao_expo(j) 
    do i = 1, small_ao_num
     num_B = small_ao_nucl(i)
     power_B(1:3)= small_ao_power(i,1:3)
     B_center(1:3) = nucl_coord(num_B,1:3)
     beta = small_ao_expo(i)
     c = 0.d0
     do k = 1, nucl_num
      Z = nucl_charge(k)
      C_center(1:3) = nucl_coord(k,1:3)
      c -= Z*NAI_pol_mult(A_center,B_center,power_A,power_B,alpha,beta,C_center,n_pt_in)
     enddo
    small_ao_nucl_elec_integral(i,j) += c*small_ao_coef_normalized(j)*small_ao_coef_normalized(i)
    enddo
   enddo
!  endif
!  if (write_ao_one_integrals) then
!   call write_one_e_integrals('ao_ne_integral', ao_nucl_elec_integral, size(ao_nucl_elec_integral,1), size(ao_nucl_elec_integral,2))
!    print *,  'AO N-e integrals written to disk'
!  endif
 END_PROVIDER


 BEGIN_PROVIDER [ double precision, small_ao_nucl_elec_integral_per_atom,(small_ao_num,small_ao_num,nucl_num)]
   BEGIN_DOC
   ! interaction nuclear electron for the AOs within the small component basis
   ! set in unrestricted kinetic balance
   END_DOC
   implicit none
   double precision               :: alpha, beta, gama, delta,c,Z
   integer                        :: num_A,num_B
   double precision               :: A_center(3),B_center(3),C_center(3)
   integer                        :: power_A(3),power_B(3)
   integer                        :: i,j,k,l,n_pt_in,m
   double precision               :: NAI_pol_mult 
! overlap_x,overlap_y,overlap_z,overlap,dx
!  if (read_ao_one_integrals) then
!   call read_one_e_integrals('ao_ne_integral', ao_nucl_elec_integral, size(ao_nucl_elec_integral,1), size(ao_nucl_elec_integral,2))
!    print *,  'AO N-e integrals read from disk'
!  else
   print*,' computing the small_ao_nucl_elec_integral_per_atom'
     small_ao_nucl_elec_integral_per_atom = 0.d0
     !        _
     ! /|  / |_)
     !  | /  | \
     !
   n_pt_in = n_pt_max_integrals
   do j = 1, small_ao_num
    num_A = small_ao_nucl(j)
    power_A(1:3)= small_ao_power(j,1:3)
    A_center(1:3) = nucl_coord(num_A,1:3)
    alpha = small_ao_expo(j) 
    do i = 1, small_ao_num
     num_B = small_ao_nucl(i)
     power_B(1:3)= small_ao_power(i,1:3)
     B_center(1:3) = nucl_coord(num_B,1:3)
     beta = small_ao_expo(i)
     c = 0.d0
     do k = 1, nucl_num
      c = 0d0
      Z = nucl_charge(k)
      C_center(1:3) = nucl_coord(k,1:3)
      c -= Z*NAI_pol_mult(A_center,B_center,power_A,power_B,alpha,beta,C_center,n_pt_in)
      small_ao_nucl_elec_integral_per_atom(i,j,k) += c*small_ao_coef_normalized(j)*small_ao_coef_normalized(i)
     enddo
    enddo
   enddo
!  endif
!  if (write_ao_one_integrals) then
!   call write_one_e_integrals('ao_ne_integral', ao_nucl_elec_integral, size(ao_nucl_elec_integral,1), size(ao_nucl_elec_integral,2))
!    print *,  'AO N-e integrals written to disk'
!  endif
 END_PROVIDER

!double precision function NAI_pol_mult(A_center,B_center,power_A,power_B,alpha,beta,C_center,n_pt_in)
!! function that calculate the folowing integral :
!!       int{dr} of (x-A_x)^ax (x-B_X)^bx exp(-alpha (x-A_x)^2 - beta (x-B_x)^2 )
!!       1/(r-R_c)
!implicit none
!integer, intent(in) :: n_pt_in
!double precision,intent(in) :: C_center(3),A_center(3),B_center(3),alpha,beta
!integer :: power_A(3),power_B(3)
!integer :: i,j,k,l,n_pt
!double precision :: P_center(3)
!double precision :: d(0:n_pt_in),pouet,coeff,rho,dist,const,pouet_2,p,p_inv,factor
!double precision :: I_n_special_exact,integrate_bourrin,I_n_bibi
!double precision ::  V_e_n,const_factor,dist_integral,tmp
!double precision :: accu,epsilo,rint
!integer :: n_pt_out,lmax
!include 'Utils/constants.include.F'
! if ( (A_center(1)/=B_center(1)).or. &
!      (A_center(2)/=B_center(2)).or. &
!      (A_center(3)/=B_center(3)).or. &
!      (A_center(1)/=C_center(1)).or. &
!      (A_center(2)/=C_center(2)).or. &
!      (A_center(3)/=C_center(3))) then
!      continue
! else
!       NAI_pol_mult = V_e_n(power_A(1),power_A(2),power_A(3),power_B(1),power_B(2),power_B(3),alpha,beta)
!      return
! endif
! p = alpha + beta
!! print*, "a"
! p_inv = 1.d0/p
! rho = alpha * beta * p_inv
! dist = 0.d0
! dist_integral = 0.d0
! do i = 1, 3
!  P_center(i) = (alpha * A_center(i) + beta * B_center(i)) * p_inv
!  dist += (A_center(i) - B_center(i))*(A_center(i) - B_center(i))
!  dist_integral += (P_center(i) - C_center(i))*(P_center(i) - C_center(i))
! enddo
! const_factor = dist*rho
! const = p * dist_integral
! if(const_factor > 80.d0)then
!  NAI_pol_mult = 0.d0
!  return
! endif
! p = alpha + beta
!! print*, "a"
! p_inv = 1.d0/p
! rho = alpha * beta * p_inv
! dist = 0.d0
! dist_integral = 0.d0
! do i = 1, 3
!  P_center(i) = (alpha * A_center(i) + beta * B_center(i)) * p_inv
!  dist += (A_center(i) - B_center(i))*(A_center(i) - B_center(i))
!  dist_integral += (P_center(i) - C_center(i))*(P_center(i) - C_center(i))
! enddo
! const_factor = dist*rho
! const = p * dist_integral
! if(const_factor > 80.d0)then
!  NAI_pol_mult = 0.d0
!  return
! endif
! factor = dexp(-const_factor)
! coeff = dtwo_pi * factor * p_inv
! lmax = 20

!!  print*, "b"
! do i = 0, n_pt_in
!   d(i) = 0.d0
! enddo
! n_pt =  2 * ( (power_A(1) + power_B(1)) +(power_A(2) + power_B(2)) +(power_A(3) + power_B(3)) )
! if (n_pt == 0) then
!  epsilo = 1.d0
!  pouet = rint(0,const)
!  NAI_pol_mult = coeff * pouet
!  return
! endif

! call give_polynom_mult_center_mono_elec(A_center,B_center,alpha,beta,power_A,power_B,C_center,n_pt_in,d,n_pt_out)

! if(n_pt_out<0)then
!  NAI_pol_mult = 0.d0
!  return
! endif
! accu = 0.d0

!! 1/r1 standard attraction integral
! epsilo = 1.d0
!! sum of integrals of type : int {t,[0,1]}  exp-(rho.(P-Q)^2 * t^2) * t^i
! do i =0 ,n_pt_out,2
!  accu +=  d(i) * rint(i/2,const)
! enddo
! NAI_pol_mult = accu * coeff

!end



 double precision function G_V_e_n(a_x,a_y,a_z,b_x,b_y,b_z,alpha,beta)
 implicit none
 !!! primitive nuclear attraction between the two primitives centered on the same atom ::
 !!!!         primitive_1 = x**(a_x) y**(a_y) z**(a_z) exp(-alpha * r**2)
 !!!!         primitive_2 = x**(b_x) y**(b_y) z**(b_z) exp(- beta * r**2)
 !!! for gaussian charge density of the nucleus
 integer :: a_x,a_y,a_z,b_x,b_y,b_z
 double precision :: alpha ,beta
 double precision :: G_V_r, V_phi, V_theta
 if(iand((a_x+b_x),1)==1.or.iand(a_y+b_y,1)==1.or.iand((a_z+b_z),1)==1)then
  G_V_e_n = 0.d0
 else
  G_V_e_n =   G_V_r(a_x+b_x+a_y+b_y+a_z+b_z+1,alpha+beta)    &
 &          * V_phi(a_x+b_x,a_y+b_y)                         &
 &          * V_theta(a_z+b_z,a_x+b_x+a_y+b_y+1)
 endif
 end
 
 double precision function G_V_r(n,alpha)
 !!!! calculate the radial part of the nuclear attraction integral which is the following integral :
 !!     integral on "r" with boundaries ( 0 ; + infinity) of [ erf(G*r)  r**n  exp(-alpha * r**2) ]
 !!! CAUTION  :: this function requires the constant sqpi = dsqrt(pi)
 implicit none
 double precision ::  alpha , fact 
 integer :: n
 include 'Utils/constants.include.F'
 if(iand(n,1)==1)then
  G_V_r = 0
 !V_r = 0.5d0 * fact(ishft(n,-1)) / (alpha ** (ishft(n,-1) + 1)) 
 else
  G_V_r = 0
 !V_r = sqpi * fact(n) / fact(ishft(n,-1)) * (0.5d0/sqrt(alpha)) ** (n+1)
 endif
 end



