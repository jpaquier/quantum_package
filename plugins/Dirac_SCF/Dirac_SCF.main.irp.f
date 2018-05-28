program Dirac_SCF
  implicit none
  BEGIN_DOC
! Print the small component basis, after unrestricted kinetic balance
  END_DOC
  integer :: i,j,k,l,l_type,n,g,h,df_L
  double precision :: integral, ao_bielec_integral, dirac_ao_bielec_integral
  include 'Utils/constants.include.F'


!print*,'***************'
!integral =0
!do i = 1, 2*ao_num
! print*,i,d_L(i)
!!integral = dirac_ao_bielec_integral(1,1,1,d_L(i))
!!print*, integral
!enddo

 print*,"*******************************"
 print*,'mo_tot_num =', mo_tot_num
 print*,'small_mo_tot_num =', small_mo_tot_num
 print*,'dirac_mo_tot_num =', dirac_mo_tot_num
 do j = 1,2*(dirac_mo_tot_num)
  print*,j, dirac_fock_matrix_eigenvalues(j)
 enddo
!print*,'**************************************************'
!print*, dirac_ao_num
!print*, dirac_mo_tot_num
!do i = 1,ao_num
! print*,'**************************************************'
! do j = 1,ao_num
! !print*, i, j, dirac_ao_bi_elec_integral(i,j)
! !print*,i, j, dirac_SCF_density_matrix_ao(i,j)
! !print*, i, j, dirac_ao_mono_elec_integral(i,j)
!  print*,i,j, dirac_mo_coef(i,j) 
! !print*,i, j, dirac_mo_mono_elec_integral(i,j)
! enddo
!enddo


!print*,''
!do j = 1,small_ao_num 
! print*,'***************************************************'
! print*,'nucleus =',small_ao_nucl(j), 'nuclear_coordinate =', nucl_coord( small_ao_nucl(j), 1 ), nucl_coord( small_ao_nucl(j), 2 ), nucl_coord(small_ao_nucl(j), 3 )
! print*,'small_ao_expo =', small_ao_expo(j)
! print*,'small_ao_power =', small_ao_power(j,1), small_ao_power(j,2), small_ao_power(j,3)
!!print*,'small_ao_coef_normalized =', small_ao_coef_normalized(j)
! print*,'***************************************************'
! do i = 1,ao_num
!  print*,'nucleus =',ao_nucl(i), 'nuclear_coordinate =', nucl_coord( ao_nucl(i), 1 ), nucl_coord( ao_nucl(i), 2 ), nucl_coord(ao_nucl(i), 3 )
!  print*,'ao_expo_ordered_transp =', ao_expo_ordered_transp(1,i)
!  print*,'ao_power =', ao_power(i,1), ao_power(i,2), ao_power(i,3) 
!! print*,'ao_coef_normalized_ordered_transp =', ao_coef_normalized_ordered_transp(1,i)
!! print*,'small_ao_deriv1_x =', small_ao_deriv1_x(j,i)
!! print*,'small_ao_deriv1_y =', small_ao_deriv1_y(j,i)
!! print*,'small_ao_deriv1_z =', small_ao_deriv1_z(j,i) 
!print*,'dirac_ao_kinetic_integral_plus =', dirac_ao_kinetic_integral_plus(j,i),Conjg(dirac_ao_kinetic_integral_plus(j,i))
!  print*,'dirac_ao_kinetic_integral_minus =', dirac_ao_kinetic_integral_minus(j,i), Conjg(dirac_ao_kinetic_integral_minus(j,i))
!  print*,'dirac_ao_kinetic_integral_z =', dirac_ao_kinetic_integral_z(j,i), Conjg(dirac_ao_kinetic_integral_z(j,i))
!  print*,'******************'
! enddo
!enddo



!print*,'small_ao_num =', small_ao_num
!do j = 1,small_ao_num 
!print*,'*********************************************'
!!print*,'nucleus =',small_ao_nucl(i) 
!!print*,'nuclear_coordinate =', nucl_coord( small_ao_nucl(i), 1 ), nucl_coord( small_ao_nucl(i), 2 ), nucl_coord(small_ao_nucl(i), 3 )
!!print*, 'small_ao_expo =', small_ao_expo(i)
!!print*,'small_ao_power =', small_ao_power(i,1), small_ao_power(i,2), small_ao_power(i,3)
!!print*,'small_ao_coef_normalized =', small_ao_coef_normalized(i) 
! print*,'******************'
!  do i= 1,small_ao_num
! !print*,'nucleus =',small_ao_nucl(j) 
! !print*,'nuclear_coordinate =', nucl_coord( small_ao_nucl(j), 1 ), nucl_coord( small_ao_nucl(j), 2 ), nucl_coord(small_ao_nucl(j), 3 )
! !print*, 'small_ao_expo =', small_ao_expo(j)
! !print*,'small_ao_power =', small_ao_power(j,1), small_ao_power(j,2), small_ao_power(j,3)
! !print*,'small_ao_coef_normalized =', small_ao_coef_normalized(j)
! !print*,'small_ao_overlap_x =', small_ao_overlap_x(i,j)
! !print*,'small_ao_overlap_y =', small_ao_overlap_y(i,j)
! !print*,'small_ao_overlap_z =', small_ao_overlap_z(i,j)
! !print*,'small_ao_overlap =', small_ao_overlap(i,j)
! !print*,'small_ao_mass_energy =', small_ao_mass_energy(i,j)
! !print*,'small_ao_nucl_elec_integral =', small_ao_nucl_elec_integral(i,j)
! !do k = 1, nucl_num 
! ! print*,'small_ao_nucl_elec_integral_per_atom =', small_ao_nucl_elec_integral_per_atom(i,j,k) 
! !enddo 
!  print*,i,j,small_ao_overlap_abs(i,j)
! !print*,i, j, dirac_SCF_density_matrix_ao(i,j)
! !print*,i,j,dirac_ao_bi_elec_integral_S_beta_S_beta(i,j)
!  print*,'******************'
! enddo
!enddo


!print*,'ao_num =', ao_num
!print*,'******************'
!print*,'*********************************************'
!print*, 'dirac_SCF_density_matrix_ao ==== '
!print*,''
!do j = 1,ao_num
! write(*,'(100((F10.4,"+i",F10.4)),X)')dirac_SCF_density_matrix_ao(j,1:ao_num)
! write(*,'(100((F10.4,"+i"),X))')dirac_SCF_density_matrix_ao(j,1:ao_num)
!!print*,'nucleus =',ao_nucl(j)
!!print*,'nuclear_coordinate =', nucl_coord(ao_nucl(j), 1 ), nucl_coord( ao_nucl(j), 2 ), nucl_coord(ao_nucl(j), 3 )
!!print*,'ao_expo_ordered_transp =', ao_expo_ordered_transp(1,i)
!!print*,'ao_power =', ao_power(i,1), ao_power(i,2), ao_power(i,3)  
! do i= 1,ao_num
! !print*,'nucleus =',ao_nucl(j)
! !print*,'nuclear_coordinate =', nucl_coord(ao_nucl(j), 1 ), nucl_coord(ao_nucl(j), 2 ), nucl_coord(ao_nucl(j), 3 )
! !print*,'ao_expo_ordered_transp =', ao_expo_ordered_transp(1,j)
! !print*,'ao_power =', ao_power(j,1), ao_power(j,2), ao_power(j,3)
! !do k = 1,nucl_num 
! ! print*, 'ao_nucl_elec_integral_per_atom =', ao_nucl_elec_integral_per_atom(i,j,k)
! !enddo
! !print*,'ao_ortho_canonical_coef =', ao_ortho_canonical_coef(j,i)
! !print*,'ao_ortho_canonical_overlap =', ao_ortho_canonical_overlap(j,i)
! !print*,'mo_mono_elec_integral =', mo_mono_elec_integral(i,j)
! !print*,'ao_mono_elec_integral =', ao_mono_elec_integral(i,j)
!   
! !integral = dirac_ao_bielec_integral(1,i,1,j)
! !print*, integral
! !integral = ao_bielec_integral(1,i,1,j)
! !print*, integral
!
! !print*,i,j,dirac_ao_bi_elec_integral_L_alpha_L_alpha(i,j)
! !print*,i,j,ao_bi_elec_integral_alpha(i,j)
!   
! !print*,i, j, dirac_SCF_density_matrix_ao(i,j)
! !print*,i, j, SCF_density_matrix_ao_alpha(i,j)
!  
!  
! !print*,i,j, dirac_mo_coef(i,j)
! !print*,i,j, mo_coef(i,j)   
!!print*,'******************'
! enddo
!enddo

!print*,'ao_num =', ao_num
!print*,'******************'
!print*,'*********************************************'
!print*, 'SCF_density_matrix_ao_alpha ==== '
!print*,''
!do j = 1,ao_num
! write(*,'(100((F10.4)),X)')SCF_density_matrix_ao_alpha(j,1:ao_num)
!enddo

!print*,'*********************'
!print*,'*********************'
!print*,'*********************'
!print*,'mo_coef ===='
!do j = 1,ao_num
!  write(*,'(100((F10.4)),X)')mo_coef(j,1:ao_num)
!enddo
!print*,'                     '
!print*,'                     '
!print*,'*********************'
!print*,'*********************'
!print*,'*********************'
!print*,'dirac_mo_coef ===='
!do j = 1,ao_num
!  write(*,'(100((F10.4)),X)')dirac_mo_coef(j,1:ao_num)
!enddo

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
 
!double precision ::  get_dirac_ao_bielec_integral,dirac_map,dirac_map_bis, tmp
!print*,''
!print*,''
!print*,''
!print*,'integrals '
!print*,''
!print*,''
!print*,''
!tmp = 0.d0
!double precision :: threshold 
!threshold = 1.d-10
!do i = 1, dirac_ao_num
! do j = 1, dirac_ao_num
!  do k = 1, dirac_ao_num
!   do l = 1, dirac_ao_num
!   dirac_map = get_dirac_ao_bielec_integral(i,k,j,l,dirac_ao_integrals_map)
!   dirac_map_bis = dirac_ao_bielec_integral(i,j,k,l)
!   if(dabs(dirac_map-dirac_map_bis).gt.threshold)then
!    print*,i,j,k,l
!    print*,dirac_map,dirac_map_bis
!    pause
!   endif
!   tmp += dabs(dirac_map-dirac_map_bis)
!   enddo
!  enddo
! enddo
!enddo
!print*,tmp, tmp/dble(dirac_ao_num**4)


end
