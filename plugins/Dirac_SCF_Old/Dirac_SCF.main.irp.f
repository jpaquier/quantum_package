program Dirac_SCF
  implicit none
  BEGIN_DOC
  !Print whatever tests one may deem useful or remotely interesting 
  END_DOC
  integer :: i,j,k,l,l_type,n,g,h,df_L
  complex*16       :: ortho(2*dirac_ao_num),testm(2,2), eigenvectors(2,2)
  include 'Utils/constants.include.F'

 do j = 1, 2*dirac_ao_num
! print*,'********************************'
 !print*,'integral =',  ao_bielec_integral(1,1,1,1)
 !print*,'dirac_integral =',  dirac_ao_bielec_integral(i,i,i,i)
 !print*,'dirac_ao_bi_elec_integral =', dirac_ao_bi_elec_integral(i,i)
 !print*,'integral_erf =',  ao_bielec_integral_erf(1,1,1,1)
 !print*,'dirac_integral_erf =',  dirac_ao_bielec_integral_erf(i,i,i,i)
 !print*,'dirac_ao_bi_elec_integral =', dirac_ao_bi_elec_integral(i,i)
 !print*,'dirac_ao_bi_elec_erf_integral =', dirac_ao_bi_elec_erf_integral(i,i)
 !print*,'dirac_ao_bi_elec_Gaunt_integral =', dirac_ao_bi_elec_Gaunt_integral(i,i)
 !print*,'dirac_ao_bi_elec_Gaunt_erf_integral =', dirac_ao_bi_elec_Gaunt_erf_integral(i,i)  
 !print*,'**************************'
  do i = 1, 2*dirac_ao_num
!  print*, i, j, dirac_mo_coef(i,j)
  enddo
 enddo
 
 do i =1, 2*dirac_ao_num
 !print*, i, eigenvalues_dirac_mono_elec_mo(i)
 enddo
 
!print*,'************'
!print*, large_ao_num,small_ao_num,dirac_ao_num
!do j = 1, 2*dirac_ao_num
! print*,'************************************************************'
! print*,'nucleus =',dirac_ao_nucl(d_L(j))
!!print*,'nuclear_coordinate =', nucl_coord(dirac_ao_nucl(j), 1 ), nucl_coord( !dirac_ao_nucl(j), 2 ), nucl_coord(dirac_ao_nucl(j), 3 )
! print*,'dirac_ao_expo_ordered_transp =', dirac_ao_expo_ordered_transp(1,d_L(j))
! print*,'dirac_ao_power =', dirac_ao_power(d_L(j),1), dirac_ao_power(d_L(j),2), dirac_ao_power(d_L(j),3)
! print*,'************************************'
! do i = 1, 2*dirac_ao_num
! print*,'nucleus =',dirac_ao_nucl(d_L(i))
!!print*,'nuclear_coordinate =', nucl_coord(dirac_ao_nucl(i), 1 ), nucl_coord(  !dirac_ao_nucl(i), 2 ), nucl_coord(dirac_ao_nucl(i), 3 )
! print*,'dirac_ao_expo_ordered_transp =', dirac_ao_expo_ordered_transp(1,d_L(i))
! print*,'dirac_ao_power =', dirac_ao_power(d_L(i),1), dirac_ao_power(d_L(i),2),dirac_ao_power(d_L(i),3)
! print*,i,j,dirac_ao_mono_elec_nucl_integral(i,j)
! print*,'****************'
!enddo
!enddo

!print*,'***************'
!integral =0
!do i = 1, dirac_ao_num
! print*,'nucleus =',dirac_ao_nucl(i) 
!!print*,'nuclear_coordinate =', nucl_coord(dirac_ao_nucl(i), 1 ), nucl_coord( dirac_ao_nucl(i), 2 ), nucl_coord(dirac_ao_nucl(i), 3 )
! print*,'dirac_ao_expo_ordered_transp =', dirac_ao_expo_ordered_transp(1,i)
! print*,'dirac_ao_power =', dirac_ao_power(i,1), dirac_ao_power(i,2), dirac_ao_power(i,3) 
! print*,'dirac_ao_coef_normalized_ordered_transp =',dirac_ao_coef_normalized_ordered_transp(1,i)
! integral =  dirac_ao_bielec_integral(1,2,1,i)
! print*, i,integral
! print*,'*************************'
!enddo


!do j = 1, 2*dirac_ao_num
!!print*,j,eigenvalues_dirac_fock_matrix_mo(j)
! do i = 1, 2*dirac_ao_num
!  print*,i,j,eigenvectors_dirac_fock_matrix_mo(i,j)
! enddo
! print*,'***********************'
!enddo

!print*,large_ao_num,small_ao_num
!do j = 1, 2*dirac_ao_num
! do i = 1, 2*dirac_ao_num
!  print*,i,j,dirac_SCF_density_matrix_ao(i,j)
! enddo
! print*,'*************'
!enddo
 
!print*,large_ao_num,small_ao_num
!do j = 1, elec_num
! do i = 1, 2*dirac_ao_num
!  print*,i,j,dirac_mo_coef_electronic(i,j)
! enddo
! print*,'*************'
!enddo

!print*,large_ao_num,small_ao_num
!do j = 1, large_ao_num
! do i = 1, large_ao_num
!  print*,i,j,large_ao_ortho_lowdin_overlap(i,j)
! !print*,i,j,ao_mono_elec_integral(i,j)
! !print*,i,j,ao_bi_elec_integral(i,j)
! enddo
! print*,'**************'
!enddo

 
!print*,large_ao_num,small_ao_num
!do j = 1, 2*dirac_ao_num
! do i = 1, 2*dirac_ao_num
! !print*,i,j,mo_overlap(i,j)
!  print*,i,j,dirac_mo_overlap(i,j)
! !print*,i,j,large_ao_ortho_canonical_overlap(i,j)
! !print*,i,j,large_ao_ortho_lowdin_overlap(i,j)
! !print*,i,j,dirac_ao_mono_elec_integral(i,j)
! !print*,i,j,dirac_mo_mono_elec_integral(i,j)
! !print*,i,j,dirac_Fock_matrix_ao(i,j)
! !print*,i,j,dirac_Fock_matrix_mo(i,j)  
! !print*,i,j,dirac_S_mo_coef(i,j)
!  print*,'******' 
! enddo
! print*,'******************'
!enddo

!print*,large_ao_num,small_ao_num
!do j = 1, 2*dirac_ao_num
! print*,j,eigenvalues_dirac_fock_matrix_mo(j)
!!do i = 1, 2*dirac_ao_num
!! print*,i,j,eigenvectors_dirac_fock_matrix_mo(i,j)
!!!print*,i,j,dirac_ao_ortho_lowdin_coef(i,j)
!!!print*,i,j,dirac_mo_coef(i,j)
!!!print*,'*****'
!!!print*,i,j,dirac_ao_mono_elec_integral(i,j)
!!!print*,i,j,dirac_ao_mono_elec_nucl_integral(i,j)
!!enddo
!!print*,'*************'
!enddo

!print*,  HF_one_electron_energy 
!print*,  dirac_HF_one_electron_energy
!print*,  dirac_HF_one_electron_energy_complex
!print*,'***********' 
!print*,  HF_two_electron_energy
!print*,  dirac_HF_two_electron_energy
!print*,  dirac_HF_two_electron_energy_complex
!print*,'***********'
!print*, HF_energy
!print*, dirac_HF_energy 

!print*,'large_ao_num =',large_ao_num
!do j = 1, large_ao_num
! do i = 1, large_ao_num
!  print*,i,j,ao_bi_elec_integral_alpha(i,j)
!  print*,i,j,dirac_ao_bi_elec_integral(i,j)
! !print*,i,j,dirac_ao_bi_elec_integral_L_alpha_L_alpha(i,j)
!  print*,'******'
! !print*,i,j,dirac_ao_bi_elec_integral(i+large_ao_num,j+large_ao_num)   
! !print*,i,j,dirac_ao_bi_elec_integral_L_beta_L_beta(i,j)  
! !print*,'**********' 
! enddo
! print*,'*************************************'
!enddo

!print*,'small_ao_num =',small_ao_num
!do j = 1, small_ao_num
! do i = 1, small_ao_num
!  print*,i,j,dirac_ao_bi_elec_integral(i+2*large_ao_num,j+2*large_ao_num)
!! print*,i,j,dirac_ao_bi_elec_integral_S_alpha_S_alpha(i,j)
!! print*,'******'
!! print*,i,j,dirac_ao_bi_elec_integral(i+(2*large_ao_num+small_ao_num),j+(2*large_ao_num+small_ao_num))   
!! print*,i,j,dirac_ao_bi_elec_integral_S_beta_S_beta(i,j)  
!! print*,'********************'
! enddo
! print*,'******************************************'
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
