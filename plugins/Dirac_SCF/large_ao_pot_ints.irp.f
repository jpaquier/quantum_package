!BEGIN_PROVIDER [ double precision, large_ao_nucl_elec_integral,(large_ao_num,large_ao_num)]
!  BEGIN_DOC
!  !Interaction nucleus-electron for the AOs within the large component 
!  ! basis set in the unrestricted kinetic balance scheme
!  END_DOC
!  implicit none
!  double precision               :: alpha, beta, gama, delta,c,Z,coef
!  integer                        :: num_A,num_B
!  double precision               :: A_center(3),B_center(3),C_center(3)
!  integer                        :: power_A(3),power_B(3)
!  integer                        :: i,j,k,l,n_pt_in,m
!  double precision               :: NAI_pol_mult 
!    large_ao_nucl_elec_integral = 0.d0
!    !        _
!    ! /|  / |_)
!    !  | /  | \
!    !
!  n_pt_in = n_pt_max_integrals
!  do j = 1, large_ao_num
!   num_A = large_ao_nucl(j)
!   power_A(1:3)= large_ao_power(j,1:3)
!   A_center(1:3) = nucl_coord(num_A,1:3)
!   do i = 1, large_ao_num
!    num_B = large_ao_nucl(i)
!    power_B(1:3)= large_ao_power(i,1:3)
!    B_center(1:3) = nucl_coord(num_B,1:3)
!    do l = 1,large_ao_prim_num(j)
!     alpha  = large_ao_expo_ordered_transp(l,j)
!     do m = 1,large_ao_prim_num(i)
!      beta  = large_ao_expo_ordered_transp(m,i)
!      c = 0.d0
!      do k = 1, nucl_num
!       Z = nucl_charge(k)
!       C_center(1:3) = nucl_coord(k,1:3)
!       c -= Z*NAI_pol_mult(A_center,B_center,power_A,power_B,alpha,beta,C_center,n_pt_in)
!      enddo
!      coef = large_ao_coef_normalized_ordered_transp(l,j)*large_ao_coef_normalized_ordered_transp(m,i)
!      large_ao_nucl_elec_integral(i,j) += c*coef
!     enddo
!    enddo
!   enddo
!  enddo
!END_PROVIDER


!BEGIN_PROVIDER [ double precision, large_ao_nucl_elec_integral_per_atom, (large_ao_num,large_ao_num,nucl_num)]
!  BEGIN_DOC
!  !Interaction nucleus-electron per atom for the AOs within the large component
!  ! basis set in unrestricted kinetic balance
!  END_DOC
!  implicit none
!  double precision               :: alpha, beta, delta,c,Z,coef
!  integer                        :: num_A,num_B
!  double precision               :: A_center(3),B_center(3),C_center(3)
!  integer                        :: power_A(3),power_B(3)
!  integer                        :: i,j,k,l,n_pt_in,m
!  double precision               :: NAI_pol_mult 
!  print*,' computing the large_ao_nucl_elec_integral_per_atom'
!  large_ao_nucl_elec_integral_per_atom = 0.d0
!    !        _
!    ! /|  / |_)
!    !  | /  | \
!    !
!  n_pt_in = n_pt_max_integrals
!  do j = 1, large_ao_num
!   num_A = large_ao_nucl(j)
!   power_A(1:3)= large_ao_power(j,1:3)
!   A_center(1:3) = nucl_coord(num_A,1:3)
!   do k = 1, nucl_num
!    Z = nucl_charge(k)
!    C_center(1:3) = nucl_coord(k,1:3)
!    do i = 1, large_ao_num
!     num_B = large_ao_nucl(i)
!     power_B(1:3)= large_ao_power(i,1:3)
!     B_center(1:3) = nucl_coord(num_B,1:3)
!    c = 0.d0
!    coef = 0.d0
!     do l = 1,large_ao_prim_num(j)
!      alpha   = large_ao_expo_ordered_transp(l,j)
!      do m = 1,large_ao_prim_num(i)
!       beta   = large_ao_expo_ordered_transp(m,i)
!       c -= Z*NAI_pol_mult(A_center,B_center,power_A,power_B,alpha,beta,C_center,n_pt_in)
!       coef = large_ao_coef_normalized_ordered_transp(l,j)*large_ao_coef_normalized_ordered_transp(m,i)
!       large_ao_nucl_elec_integral_per_atom(i,j,k) += c*coef
!      enddo
!     enddo
!    enddo
!   enddo
!  enddo
!END_PROVIDER


