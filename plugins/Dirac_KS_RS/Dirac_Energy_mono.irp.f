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
   open (10, file='Energy_DHF_LRC.dat',position ='append') 
   write(10,*) 'dirac_HF_one_electron_energy=', dirac_HF_one_electron_energy
   write(10,*) 'dirac_HF_one_electron_mass_energy=', dirac_HF_one_electron_mass_energy
   write(10,*) 'dirac_HF_one_electron_kinetic_energy=', dirac_HF_one_electron_kinetic_energy
   write(10,*) 'dirac_HF_one_electron_nucl_energy=', dirac_HF_one_electron_nucl_energy
 end

