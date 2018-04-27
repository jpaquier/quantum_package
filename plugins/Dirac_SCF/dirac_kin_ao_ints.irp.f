 BEGIN_PROVIDER [ double precision, small_ao_deriv_1_x, (small_ao_num,ao_num)]
 &BEGIN_PROVIDER [ double precision, small_ao_deriv_1_y, (small_ao_num,ao_num)]
 &BEGIN_PROVIDER [ double precision, small_ao_deriv_1_z, (small_ao_num,ao_num)]
 BEGIN_DOC
 ! array of the integrals of AO_i * d/dx  AO_j
 ! array of the integrals of AO_i * d/dy  AO_j
 ! array of the integrals of AO_i * d/dz  AO_j
 ! for AO_i within the small component basis and AO_j within the large component basis
 END_DOC
 implicit none
  integer :: i,j,n,l
  double precision :: f, tmp
  integer :: dim1
  double precision :: overlap, overlap_x, overlap_y, overlap_z
  double precision :: alpha , beta 
  double precision :: A_center(3), B_center(3)
  integer :: power_A(3), power_B(3)
  double precision :: lower_exp_val, dx, c,accu_x,accu_y,accu_z
  integer :: i_component
  dim1=500
  lower_exp_val = 40.d0
  small_ao_deriv_1_x= 0.d0
  small_ao_deriv_1_y= 0.d0
  small_ao_deriv_1_z= 0.d0
  do j=1,ao_num
   A_center(1) = nucl_coord( ao_nucl(j), 1 )
   A_center(2) = nucl_coord( ao_nucl(j), 2 )
   A_center(3) = nucl_coord( ao_nucl(j), 3 )
   power_A(1)  = ao_power( j, 1 )
   power_A(2)  = ao_power( j, 2 )
   power_A(3)  = ao_power( j, 3 )
   alpha   = ao_expo_ordered_transp(1,j)
   do i= 1,small_ao_num
    B_center(1) = nucl_coord( small_ao_nucl(i), 1 )
    B_center(2) = nucl_coord( small_ao_nucl(i), 2 )
    B_center(3) = nucl_coord( small_ao_nucl(i), 3 )
    power_B(1)  = small_ao_power( i, 1 )
    power_B(2)  = small_ao_power( i, 2 )
    power_B(3)  = small_ao_power( i, 3 )
    accu_x = 0.d0
    accu_y = 0.d0
    accu_z = 0.d0
    beta  = small_ao_expo(i)
    call overlap_gaussian_xyz(A_center,B_center,alpha,beta,power_A,power_B,overlap_x,overlap_y,overlap_z,overlap,dim1)
    c = small_ao_coef_normalized(i) * ao_coef_normalized_ordered_transp(1,j)
    i_component = 1
    call overlap_bourrin_deriv_x(i_component,A_center,B_center,alpha,beta,power_A,power_B,dx,lower_exp_val,tmp,dim1)
    accu_x += c*(tmp*overlap_y*overlap_z)
    i_component = 2
    call overlap_bourrin_deriv_x(i_component,A_center,B_center,alpha,beta,power_A,power_B,dx,lower_exp_val,tmp,dim1)
    accu_y += c*(tmp*overlap_x*overlap_z)
    i_component = 3
    call overlap_bourrin_deriv_x(i_component,A_center,B_center,alpha,beta,power_A,power_B,dx,lower_exp_val,tmp,dim1)
    accu_z += c*(tmp*overlap_y*overlap_x)
    small_ao_deriv_1_x(i,j) = accu_x
    small_ao_deriv_1_y(i,j) = accu_y
    small_ao_deriv_1_z(i,j) = accu_z
   enddo
  enddo
 END_PROVIDER


 BEGIN_PROVIDER [double precision, dirac_ao_kinetic_integral, (small_ao_num,ao_num)]
 BEGIN_DOC
 ! array of the priminitve basis kinetic integrals
 !  \langle \chi_i |\hat{T}| \chi_j \rangle
 END_DOC
 integer :: i,j
!if (read_ao_one_integrals) then
! call read_one_e_integrals('ao_kinetic_integral', ao_kinetic_integral, size(ao_kinetic_integral,1), size(ao_kinetic_integral,2))
!   print *,  'AO kinetic integrals read from disk'
! else
!   print *,  'Computing AO kinetic integrals '
! do j = 1, small_ao_num
!  do i = 1, ao_num
!   ao_kinetic_integral(i,j) = -0.5d0 * (ao_deriv2_x(i,j) + ao_deriv2_y(i,j) + ao_deriv2_z(i,j) )
!  enddo
! enddo
!endif
! if (write_ao_one_integrals) then
!  call write_one_e_integrals('ao_kinetic_integral', ao_kinetic_integral, size(ao_kinetic_integral,1), size(ao_kinetic_integral,2))
!  print *,  'AO kinetic integrals written to disk'
! endif
 END_PROVIDER
