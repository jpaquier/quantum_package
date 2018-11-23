!BEGIN_PROVIDER [ double precision, large_ao_overlap,(large_ao_num,large_ao_num) ]
!&BEGIN_PROVIDER [ double precision, large_ao_overlap_x,(large_ao_num,large_ao_num) ]
!&BEGIN_PROVIDER [ double precision, large_ao_overlap_y,(large_ao_num,large_ao_num) ]
!&BEGIN_PROVIDER [ double precision, large_ao_overlap_z,(large_ao_num,large_ao_num) ]
! implicit none
! BEGIN_DOC
!!Overlap between atomic basis functions of the large component:
!!:math:`\int \chi_i(r) \chi_j(r) dr)`
!!with correct overlap_x*overlap_y*overlap_z=overlap
! END_DOC
! integer :: i,j,n,l
! double precision :: f
! integer :: dim1
! double precision :: large_overlap, large_overlap_x, large_overlap_y, large_overlap_z
! double precision :: alpha, beta, c
! double precision :: A_center(3), B_center(3)
! integer :: power_A(3), power_B(3)
! dim1=500
! c = 1
! do j=1, large_ao_num
!  A_center(1) = nucl_coord( large_ao_nucl(j), 1 )
!  A_center(2) = nucl_coord( large_ao_nucl(j), 2 )
!  A_center(3) = nucl_coord( large_ao_nucl(j), 3 )
!  power_A(1)  = large_ao_power( j, 1 )
!  power_A(2)  = large_ao_power( j, 2 )
!  power_A(3)  = large_ao_power( j, 3 )
!  do i= 1, large_ao_num
!   large_ao_overlap(i,j)= 0.d0
!   large_ao_overlap_x(i,j)= 0.d0
!   large_ao_overlap_y(i,j)= 0.d0
!   large_ao_overlap_z(i,j)= 0.d0
!   B_center(1) = nucl_coord( large_ao_nucl(i), 1 )
!   B_center(2) = nucl_coord( large_ao_nucl(i), 2 )
!   B_center(3) = nucl_coord( large_ao_nucl(i), 3 )
!   power_B(1)  = large_ao_power( i, 1 )
!   power_B(2)  = large_ao_power( i, 2 )
!   power_B(3)  = large_ao_power( i, 3 )
!   do n = 1,large_ao_prim_num(j)
!    alpha  = large_ao_expo_ordered_transp(n,j)
!    do l = 1, large_ao_prim_num(i)
!     beta  = large_ao_expo_ordered_transp(l,i)
!     call overlap_gaussian_xyz(A_center,B_center,alpha,beta,power_A,power_B,large_overlap_x,large_overlap_y,large_overlap_z,large_overlap,dim1)
!     c = large_ao_coef_normalized_ordered_transp(n,j) * large_ao_coef_normalized_ordered_transp(l,i)
!     large_ao_overlap(i,j) += c*large_overlap
!     large_ao_overlap_x(i,j) += (c**(1/3.))*large_overlap_x
!     large_ao_overlap_y(i,j) += (c**(1/3.))*large_overlap_y
!     large_ao_overlap_z(i,j) += (c**(1/3.))*large_overlap_z
!    enddo
!   enddo
!  enddo
! enddo
!END_PROVIDER

!BEGIN_PROVIDER [ double precision, large_ao_overlap_abs,(large_ao_num,large_ao_num) ]
! implicit none
! BEGIN_DOC  
!!Overlap between absolute value of large component atomic basis functions:
!!:math:`\int |\chi_i(r)| |\chi_j(r)| dr)`
! END_DOC
! integer :: i,j,n,l
! double precision :: f
! integer :: dim1
! double precision :: large_overlap, large_overlap_x, large_overlap_y, large_overlap_z
! double precision :: alpha , beta , c 
! double precision :: A_center(3), B_center(3)
! integer :: power_A(3), power_B(3)
! double precision :: lower_exp_val, dx
! dim1=500
! lower_exp_val = 40.d0
! c = 1
! do j=1, large_ao_num
!  A_center(1) = nucl_coord( large_ao_nucl(j), 1 )
!  A_center(2) = nucl_coord( large_ao_nucl(j), 2 )
!  A_center(3) = nucl_coord( large_ao_nucl(j), 3 )
!  power_A(1)  = large_ao_power( j, 1 )
!  power_A(2)  = large_ao_power( j, 2 )
!  power_A(3)  = large_ao_power( j, 3 )
!  do i= 1, large_ao_num
!   B_center(1) = nucl_coord( large_ao_nucl(i), 1 )
!   B_center(2) = nucl_coord( large_ao_nucl(i), 2 )
!   B_center(3) = nucl_coord( large_ao_nucl(i), 3 )
!   power_B(1)  = large_ao_power( i, 1 )
!   power_B(2)  = large_ao_power( i, 2 )
!   power_B(3)  = large_ao_power( i, 3 )
!   do n = 1,large_ao_prim_num(j)
!    alpha  = large_ao_expo_ordered_transp(n,j)
!    do l = 1, large_ao_prim_num(i)
!     beta  = large_ao_expo_ordered_transp(l,i)
!     call overlap_x_abs(A_center(1),B_center(1),alpha,beta,power_A(1),power_B(1),large_overlap_x,lower_exp_val,dx,dim1)
!     call overlap_x_abs(A_center(2),B_center(2),alpha,beta,power_A(2),power_B(2),large_overlap_y,lower_exp_val,dx,dim1)
!     call overlap_x_abs(A_center(3),B_center(3),alpha,beta,power_A(3),power_B(3),large_overlap_z,lower_exp_val,dx,dim1)
!     c = large_ao_coef_normalized_ordered_transp(n,j) * large_ao_coef_normalized_ordered_transp(l,i)
!     large_ao_overlap_abs(i,j) += abs(c) * large_overlap_x * large_overlap_y * large_overlap_z
!    enddo
!   enddo
!  enddo
! enddo 
!END_PROVIDER

