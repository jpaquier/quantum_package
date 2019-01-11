 BEGIN_PROVIDER [ complex*16, dirac_HF_one_electron_energy_complex]
 &BEGIN_PROVIDER [ double precision, dirac_HF_one_electron_energy]
  implicit none
  BEGIN_DOC
  !One-electron energy of the Nucleus-Electron interaction
  !The energy is supposed to be a real, thus we check for its complex part to be
  ! a VERY small artifact and take only its real part
  END_DOC
  integer :: i,j
  dirac_HF_one_electron_energy_complex = (0.d0,0.d0)
  do j=1, 2*dirac_ao_num
   do i=1, 2*dirac_ao_num
    dirac_HF_one_electron_energy_complex += dirac_ao_mono_elec_integral(i,j)* dirac_SCF_density_matrix_ao(j,i) 
   enddo
  enddo
  dirac_HF_one_electron_energy = real(dirac_HF_one_electron_energy_complex)
  if (aimag(dirac_HF_one_electron_energy_complex) .gt. 1.d-10 ) then
  print*, 'Warning! The energy is not real'
  print*, 'dirac_HF_one_electron_energy_complex =',dirac_HF_one_electron_energy_complex
  STOP
  endif
 END_PROVIDER

