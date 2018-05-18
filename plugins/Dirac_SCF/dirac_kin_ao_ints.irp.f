 BEGIN_PROVIDER [ double precision, small_ao_deriv1_x, (small_ao_num,ao_num)]
 &BEGIN_PROVIDER [ double precision, small_ao_deriv1_y, (small_ao_num,ao_num)]
 &BEGIN_PROVIDER [ double precision, small_ao_deriv1_z, (small_ao_num,ao_num)]
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
  double precision :: overlap_x0, overlap_y0, overlap_z0
  double precision :: alpha , beta, c 
  double precision :: A_center(3), B_center(3)
  integer :: power_A(3), power_B(3)
  double precision :: deriv_tmp1, deriv_tmp2
  dim1=100
  small_ao_deriv1_x= 0.d0
  small_ao_deriv1_y= 0.d0
  small_ao_deriv1_z= 0.d0
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
    beta  = small_ao_expo(i)
    call overlap_gaussian_xyz(A_center,B_center,alpha,beta,power_A,power_B,overlap_x0,overlap_y0,overlap_z0,overlap,dim1)
    c = small_ao_coef_normalized(i) * ao_coef_normalized_ordered_transp(1,j)

    power_A(1) = power_A(1)-1
    if (power_A(1)>-1) then
     call overlap_gaussian_xyz(A_center,B_center,alpha,beta,power_A,power_B,deriv_tmp1,overlap_y,overlap_z,overlap,dim1)
    else
     deriv_tmp1 = 0.d0
    endif
    power_A(1) = power_A(1)+2 
    call overlap_gaussian_xyz(A_center,B_center,alpha,beta,power_A,power_B,deriv_tmp2,overlap_y,overlap_z,overlap,dim1)
    power_A(1) = power_A(1)-1   
    small_ao_deriv1_x(i,j) += c*(Power_A(1)*deriv_tmp1 - 2.d0*alpha*deriv_tmp2) * overlap_y0 * overlap_z0
 
    power_A(2) = power_A(2)-1
    if (power_A(2)>-1) then
     call overlap_gaussian_xyz(A_center,B_center,alpha,beta,power_A,power_B,overlap_x,deriv_tmp1,overlap_z,overlap,dim1)
    else
     deriv_tmp1 = 0.d0
    endif
    power_A(2) = power_A(2)+2 
    call overlap_gaussian_xyz(A_center,B_center,alpha,beta,power_A,power_B,overlap_x,deriv_tmp2,overlap_z,overlap,dim1)
    power_A(2) = power_A(2)-1   
    small_ao_deriv1_y(i,j) += c*(Power_A(2)*deriv_tmp1 - 2.d0*alpha*deriv_tmp2) * overlap_x0 * overlap_z0

    power_A(3) = power_A(3)-1
    if (power_A(3)>-1) then
     call overlap_gaussian_xyz(A_center,B_center,alpha,beta,power_A,power_B,overlap_x,overlap_y,deriv_tmp1,overlap,dim1)
    else
     deriv_tmp1 = 0.d0
    endif
    power_A(3) = power_A(3)+2 
    call overlap_gaussian_xyz(A_center,B_center,alpha,beta,power_A,power_B,overlap_x,overlap_y,deriv_tmp2,overlap,dim1)
    power_A(3) = power_A(3)-1   
    small_ao_deriv1_z(i,j) += c*(Power_A(3)*deriv_tmp1 - 2.d0*alpha*deriv_tmp2) * overlap_x0 * overlap_y0
   enddo
  enddo
 END_PROVIDER


 BEGIN_PROVIDER [complex*16, dirac_ao_kinetic_integral_z, (small_ao_num,ao_num)]
 &BEGIN_PROVIDER [complex*16, dirac_ao_kinetic_integral_plus, (small_ao_num,ao_num)]
 &BEGIN_PROVIDER [complex*16, dirac_ao_kinetic_integral_minus, (small_ao_num,ao_num)]
 BEGIN_DOC
 ! array of the priminitve basis kinetic integrals
 !  \langle \chi_i |c{\alpha}.\hat_{p}| \chi_j \rangle
 ! for AO_i within the small component basis and AO_j within the large component basis
 END_DOC
 integer :: i,j
 double precision :: c
 c = 137.0359895
!if (read_ao_one_integrals) then
! call read_one_e_integrals('ao_kinetic_integral', ao_kinetic_integral, size(ao_kinetic_integral,1), size(ao_kinetic_integral,2))
!   print *,  'AO kinetic integrals read from disk'
! else
 print *,  'Computing AO kinetic integrals (SL bloc)'
  do j = 1, ao_num
   do i = 1, small_ao_num
    dirac_ao_kinetic_integral_z(i,j) = - c * (0.d0,1.d0) * small_ao_deriv1_z(i,j)
    dirac_ao_kinetic_integral_plus(i,j) = - c * (0.d0,1.d0) * (small_ao_deriv1_x(i,j) + (0.d0,1.d0)*small_ao_deriv1_y(i,j))
    dirac_ao_kinetic_integral_minus(i,j) = - c * (0.d0,1.d0) * (small_ao_deriv1_x(i,j) - (0.d0,1.d0)*small_ao_deriv1_y(i,j))
    enddo
  enddo
!endif
! if (write_ao_one_integrals) then
!  call write_one_e_integrals('ao_kinetic_integral', ao_kinetic_integral, size(ao_kinetic_integral,1), size(ao_kinetic_integral,2))
!  print *,  'AO kinetic integrals written to disk'
! endif
 END_PROVIDER
