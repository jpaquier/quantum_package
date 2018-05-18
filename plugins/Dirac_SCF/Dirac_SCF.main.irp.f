program Dirac_SCF
  implicit none
  BEGIN_DOC
! Print the small component basis, after unrestricted kinetic balance
  END_DOC
  integer :: i,j,k,l_type,n
  double precision :: fact
  include 'Utils/constants.include.F'

!print*,''
!do j = 1,small_ao_num 
! print*,'***************************************************'
! print*,'nucleus =',small_ao_nucl(j), 'nuclear_coordinate =', nucl_coord( small_ao_nucl(j), 1 ), nucl_coord( small_ao_nucl(j), 2 ), nucl_coord(small_ao_nucl(j), 3 )
! print*,'small_ao_expo =', small_ao_expo(j)
! print*,'small_ao_power =', small_ao_power(j,1), small_ao_power(j,2), small_ao_power(j,3)
! print*,'small_ao_coef_normalized =', small_ao_coef_normalized(j)
! print*,'***************************************************'
! do i= 1,ao_num
!  print*,'nucleus =',ao_nucl(i), 'nuclear_coordinate =', nucl_coord( ao_nucl(i), 1 ), nucl_coord( ao_nucl(i), 2 ), nucl_coord(ao_nucl(i), 3 )
!  print*,'ao_expo_ordered_transp =', ao_expo_ordered_transp(1,i)
!  print*,'ao_power =', ao_power(i,1), ao_power(i,2), ao_power(i,3) 
!  print*,'ao_coef_normalized_ordered_transp =', ao_coef_normalized_ordered_transp(1,i)
!  print*,'small_ao_deriv_1_x =', small_ao_deriv_1_x(j,i)
!  print*,'small_ao_deriv_1_y =', small_ao_deriv_1_y(j,i) 
!  print*,'small_ao_deriv_1_z =', small_ao_deriv_1_z(j,i) 
!  print*,'dirac_ao_kinetic_integral_z =', dirac_ao_kinetic_integral_plus(j,i)
!  print*,'dirac_ao_kinetic_integral_z =', dirac_ao_kinetic_integral_minus(j,i)
!  print*,'dirac_ao_kinetic_integral_z =', dirac_ao_kinetic_integral_z(j,i)
!  print*,'******************'
! enddo
!enddo

!print*,'ao_num =', ao_num
!do i = 1, ao_num
!print*,'*********************************************'
! print*,'nucleus =', ao_nucl(i) 
!!print*,'nuclear_coordinate =', nucl_coord( ao_nucl(i), 1 ), nucl_coord( ao_nucl(i), 2 ), nucl_coord(ao_nucl(i), 3 )
! print*,'ao_expo_ordered_transp =', ao_expo_ordered_transp(1,i)
! print*,'ao_power =', ao_power(i,1), ao_power(i,2), ao_power(i,3) 
!!print*,'ao_coef_normalized_ordered_transp =', ao_coef_normalized_ordered_transp(1,i)
! print*,'******************'
!  do j= 1, ao_num
!    print*,'nucleus =', ao_nucl(j) 
!   !print*,'nuclear_coordinate =', nucl_coord( ao_nucl(j), 1 ), nucl_coord( ao_nucl(j), 2 ), nucl_coord(ao_nucl(j), 3 )
!    print*,'ao_expo_ordered_transp =', ao_expo_ordered_transp(1,j)
!    print*,'ao_power =', ao_power(j,1), ao_power(j,2), ao_power(j,3)  
!   !print*,'ao_coef_normalized_ordered_transp =', ao_coef_normalized_ordered_transp(1,j)
!    print*,'ao_overlap =', ao_overlap(i,j)
!    print*,'ao_nucl_elec_integral =', ao_nucl_elec_integral(i,j)
!    print*,'******************'
! enddo
!enddo


 print*,'small_ao_num =', small_ao_num
 do i = 1,small_ao_num 
 print*,'*********************************************'
  print*,'nucleus =',small_ao_nucl(i) 
 !print*,'nuclear_coordinate =', nucl_coord( small_ao_nucl(i), 1 ), nucl_coord( small_ao_nucl(i), 2 ), nucl_coord(small_ao_nucl(i), 3 )
  print*, 'small_ao_expo =', small_ao_expo(i)
  print*,'small_ao_power =', small_ao_power(i,1), small_ao_power(i,2), small_ao_power(i,3)
 !print*,'small_ao_coef_normalized =', small_ao_coef_normalized(i) 
  print*,'******************'
  print*, i, dble(i-1) 
 ! do j= 1,small_ao_num
 ! print*,'nucleus =',small_ao_nucl(j) 
 !!print*,'nuclear_coordinate =', nucl_coord( small_ao_nucl(j), 1 ), nucl_coord( small_ao_nucl(j), 2 ), nucl_coord(small_ao_nucl(j), 3 )
 ! print*, 'small_ao_expo =', small_ao_expo(j)
 ! print*,'small_ao_power =', small_ao_power(j,1), small_ao_power(j,2), small_ao_power(j,3)
 !!print*,'small_ao_coef_normalized =', small_ao_coef_normalized(j)
 !!print*,'small_ao_overlap_x =', small_ao_overlap_x(i,j)
 !!print*,'small_ao_overlap_y =', small_ao_overlap_y(i,j)
 !!print*,'small_ao_overlap_z =', small_ao_overlap_z(i,j)
 ! print*,'small_ao_overlap =', small_ao_overlap(i,j)
 ! print*,'small_ao_mass_energy =', small_ao_mass_energy(i,j)
 ! print*,'small_ao_nucl_elec_integral =', small_ao_nucl_elec_integral(i,j)
 ! do k = 1, nucl_num 
 !  print*,'small_ao_nucl_elec_integral_per_atom =', small_ao_nucl_elec_integral_per_atom(i,j,k) 
 ! enddo 
 ! print*, iand((small_ao_power(i,1)+small_ao_power(j,1)),1)
 ! print*, iand((small_ao_power(i,2)+small_ao_power(j,2)),1)
 ! print*, iand((small_ao_power(i,3)+small_ao_power(j,3)),1)
 !print*,'******************'
 !enddo
 enddo




!do i = 1, nucl_num
! print*,''
! print*,'**************************************************'
! print*,'nucleus = ',i
! print*,''
! print*, "n_max_of_small_component_expo =", nmax_of_small_component_expo
! print*, "number_of_small_component_ao =", number_of_small_component_ao_per_atom(i) 
! print*,'***********************'
! do l_type = 0,7
!  print*,"l_type = ", l_type
!  print*,"number_of_small_component_expo =", number_of_small_component_expo_per_shell_per_atom(l_type,i) 
!  do k = 1, number_of_small_component_expo_per_shell_per_atom(l_type,i)
!   print*, k,small_component_expo_per_shell_per_atom(k,l_type,i)
!  enddo
!  print*, 'number_of_small_component_ao =', number_of_small_component_ao_per_shell_per_atom(l_type,i)
!  print*,'****************'  
! enddo
!enddo



end
