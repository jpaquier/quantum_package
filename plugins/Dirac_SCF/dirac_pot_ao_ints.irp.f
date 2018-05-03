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

 BEGIN_PROVIDER[double precision, teeest, (small_ao_num,small_ao_num)]
  BEGIN_DOC
 !This is a teeest, because iand, ior, ishift are a mystery in my mind 
  END_DOC
   implicit none
   double precision               :: alpha , beta 
   integer                        :: i,j,k
   do i = 1, small_ao_num
    do j = 1, small_ao_num
     teeest = iand(i+j,1)
    enddo
   enddo 
 END_PROVIDER

!double precision function G_V_e_n(a_x,a_y,a_z,b_x,b_y,b_z,alpha,beta)
!implicit none
!!!! primitive nuclear attraction between the two primitives centered on the same atom ::
!!!!!         primitive_1 = x**(a_x) y**(a_y) z**(a_z) exp(-alpha * r**2)
!!!!!         primitive_2 = x**(b_x) y**(b_y) z**(b_z) exp(- beta * r**2)
!!!! for gaussian charge density of the nucleus
!integer :: a_x,a_y,a_z,b_x,b_y,b_z
!double precision :: alpha ,beta
!double precision :: V_r, V_phi, V_theta
!if(iand((a_x+b_x),1)==1.or.iand(a_y+b_y,1)==1.or.iand((a_z+b_z),1)==1)then
! G_V_e_n = 0.d0
!else
! G_V_e_n =   G_V_r(a_x+b_x+a_y+b_y+a_z+b_z+1,alpha+beta)    &
!&          * V_phi(a_x+b_x,a_y+b_y)                         &
!&          * V_theta(a_z+b_z,a_x+b_x+a_y+b_y+1)
!endif
!end
!
!double precision function G_V_r(n,alpha)
!!!!! calculate the radial part of the nuclear attraction integral which is the following integral :
!!!     integral on "r" with boundaries ( 0 ; + infinity) of [ erf(G*r)  r**n  exp(-alpha * r**2) ]
!!!! CAUTION  :: this function requires the constant sqpi = dsqrt(pi)
!implicit none
!double precision ::  alpha , fact 
!integer :: n
!include 'Utils/constants.include.F'
!if(iand(n,1)==1)then
! V_r = 0.5d0 * fact(ishft(n,-1)) / (alpha ** (ishft(n,-1) + 1)) 
!else
! V_r = sqpi * fact(n) / fact(ishft(n,-1)) * (0.5d0/sqrt(alpha)) ** (n+1)
!endif
!end


!double precision function V_phi(n,m)
!implicit none
!!!!! calculate the angular "phi" part of the nuclear attraction integral wich is the following integral :
!!!     integral on "phi" with boundaries ( 0 ; 2 pi) of [ cos(phi) **n  sin(phi) **m ]
!integer :: n,m, i
!double precision ::  prod, Wallis
! prod = 1.d0
! do i = 0,ishft(n,-1)-1
!  prod = prod/ (1.d0 + dfloat(m+1)/dfloat(n-i-i-1))
! enddo
! V_phi = 4.d0 * prod * Wallis(m)
!end

!double precision function V_theta(n,m)
!implicit none
!!!!! calculate the angular "theta" part of the nuclear attraction integral wich is the following integral :
!!!     integral on "theta" with boundaries ( 0 ; pi) of [ cos(theta) **n  sin(theta) **m ]
!integer :: n,m,i
!double precision ::  Wallis, prod
!include 'Utils/constants.include.F'
! V_theta = 0.d0
! prod = 1.d0
! do i = 0,ishft(n,-1)-1
!  prod = prod / (1.d0 + dfloat(m+1)/dfloat(n-i-i-1))
! enddo
! V_theta = (prod+prod) * Wallis(m)
!end

!double precision function Wallis(n)
!!!!! calculate the Wallis integral :
!!!     integral on "theta" with boundaries ( 0 ; pi/2) of [ cos(theta) **n ]
!implicit none
!double precision :: fact
!integer :: n,p
!include 'Utils/constants.include.F'
!if(iand(n,1)==0)then
! Wallis = fact(ishft(n,-1))
! Wallis = pi  * fact(n) / (dble(ibset(0_8,n)) * (Wallis+Wallis)*Wallis) 
!else
! p = ishft(n,-1)
! Wallis = fact(p)
! Wallis = dble(ibset(0_8,p+p)) * Wallis*Wallis / fact(p+p+1)
!endif
!end

