 program dirac_energy_mono
  BEGIN_DOC
  ! produce the dirac energy
  END_DOC
   call run_energy_mono
 end

 subroutine run_energy_mono
  BEGIN_DOC
  ! Gives the detailed monoelectronic energy
  END_DOC
  use bitmasks
  implicit none 
  integer :: i,length
   print*,'**********'
   print*, 'dirac_HF_one_electron_energy=', dirac_HF_one_electron_energy
   print*, 'dirac_HF_one_electron_mass_energy=', dirac_HF_one_electron_mass_energy
   print*, 'dirac_HF_one_electron_kinetic_energy=', dirac_HF_one_electron_kinetic_energy
   print*, 'dirac_HF_one_electron_nucl_energy=', dirac_HF_one_electron_nucl_energy
   open (10, file='Energy_DHF_mono.dat',position ='append') 
   write(10,*)  dirac_HF_one_electron_energy
   open (11, file='Energy_DHF_mono_mass.dat',position ='append')
   write(11,*)  dirac_HF_one_electron_mass_energy
   open (12, file='Energy_DHF_mono_kin.dat',position ='append')
   write(12,*)  dirac_HF_one_electron_kinetic_energy
   open (13, file='Energy_DHF_mono_nucl.dat',position ='append')
   write(13,*)  dirac_HF_one_electron_nucl_energy
 end

