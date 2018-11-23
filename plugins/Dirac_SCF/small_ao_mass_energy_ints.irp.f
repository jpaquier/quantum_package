!BEGIN_PROVIDER [double precision, small_ao_mass_energy, (small_ao_num, small_ao_num)]
!implicit none
! BEGIN_DOC
! !Mass energy times the overlap matrix between two AOs within the small
! ! component basis, for a -mc**2 shifted mass energy
! END_DOC
! integer :: i,j
! double precision :: m,c
! m = 1.d0
!!c = 137.0359998d0 
! c = speed_of_light
! do i = 1, small_ao_num
!  do j = 1, small_ao_num
!   small_ao_mass_energy(i,j) = -2 * m * (c**2.) * small_ao_overlap(i,j)
!  enddo
! enddo  
!END_PROVIDER
