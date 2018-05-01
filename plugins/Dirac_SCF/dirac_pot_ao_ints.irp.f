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


 BEGIN_PROVIDER [ double precision, small_ao_nucl_elec_integral_per_atom,(small_ao_num,small_ao_num)]
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
