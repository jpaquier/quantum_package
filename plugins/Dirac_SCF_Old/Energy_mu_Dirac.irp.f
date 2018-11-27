!program energy_mu_dirac
! BEGIN_DOC
! ! produce the dirac energy
! END_DOC
!  call run_energy_mu 
!end

!subroutine run_energy_mu
! BEGIN_DOC
! ! Gives the energy for a given value of mu_erf
! END_DOC
! use bitmasks
! implicit none 
! integer :: i,length
!!Choose Interaction
! if (dirac_interaction == "Coulomb") then
!  print*,'**********'
!  print*,'Long-range Coulomb interaction'
!  print*, 'mu_erf =',mu_erf
!  print*, 'dirac_SCF_Coulomb_erf_energy =',dirac_SCF_Coulomb_erf_energy
!  open (10, file='Energy_DHF_LRC.dat',position ='append') 
!  write(10,*) mu_erf,dirac_SCF_Coulomb_erf_energy
! elseif (dirac_interaction == "Coulomb_Gaunt") then
!  print*,'**********'
!  print*,'Long-range Coulomb-Gaunt interaction'
!  print*, 'mu_erf =',mu_erf
!  print*, 'dirac_SCF_Coulomb_Gaunt_erf_energy =',dirac_SCF_Coulomb_Gaunt_erf_energy
!  open (10, file='Energy_DHF_LRCG.dat', position='append')
!  write(10,*) mu_erf,dirac_SCF_Coulomb_Gaunt_erf_energy
! else
!  print *,  'Unrecognized dirac_interaction : '//dirac_interaction
!  stop 1
! endif
!!mu_erf+=0.05
! call ezfio_set_dft_keywords_mu_erf(mu_erf)
!end

