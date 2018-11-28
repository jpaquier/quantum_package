 program dirac_hartree_mu
  BEGIN_DOC
  ! produce the dirac energy
  END_DOC
   call run_hartree_mu 
 end

 subroutine run_hartree_mu
  BEGIN_DOC
  ! Gives the energy for a given value of mu_erf
  END_DOC
  use bitmasks
  implicit none 
  integer :: i,length
 !Choose Interaction
  if (dirac_interaction == "Coulomb") then
   print*,'**********'
   print*,'Long-range Coulomb interaction'   
   print*, 'mu_erf =',mu_erf
   print*, 'dirac_HF_two_electron_C_Hartree_energy=', dirac_HF_two_electron_C_Hartree_energy
   open (10, file='Hartree_DHF_LRC.dat',position ='append') 
   write(10,*) mu_erf, dirac_HF_two_electron_C_Hartree_energy
  elseif (dirac_interaction == "Coulomb_Gaunt") then
   print*,'**********'
   print*,'Long-range Coulomb-Gaunt interaction'
   print*, 'mu_erf =',mu_erf
   print*, 'dirac_HF_two_electron_C_G_Hartree_energy=', dirac_HF_two_electron_C_G_Hartree_energy
   open (10, file='Hartree_DHF_LRCG.dat', position='append')
   write(10,*) mu_erf, dirac_HF_two_electron_C_G_Hartree_energy
  else
   print *,  'Unrecognized dirac_interaction : '//dirac_interaction
   stop 1
  endif
  mu_erf+=0.05d0
  call ezfio_set_dft_keywords_mu_erf(mu_erf)
 end

