 program dirac_scf
  BEGIN_DOC
  ! Produces `Dirac_Hartree_Fock` MO orbital 
  END_DOC
  call create_dirac_guess
  call run_dirac
  call print_dirac_energies
  call print_dirac_mo_coef
 end

 subroutine create_dirac_guess
  implicit none
  BEGIN_DOC
  ! Create a MO guess if no MOs are present in the EZFIO directory
  END_DOC
 logical                        :: exists_Re, exists_Im
  PROVIDE ezfio_filename
  call ezfio_has_dirac_mo_basis_dirac_mo_coef_Re(exists_Re)
  call ezfio_has_dirac_mo_basis_dirac_mo_coef_Im(exists_Im)
  if (exists_Re .and. exists_Im) then
   dirac_mo_coef = (1,0)*dirac_mo_coef_Re + (0,1)*dirac_mo_coef_Im
   print*, 'dirac_mo_coef from previous calculation'
  else 
   if (dirac_mo_guess_type == "HCore") then
    dirac_mo_coef = eigenvectors_dirac_mono_elec_ao
    TOUCH dirac_mo_coef       
    dirac_mo_label = 'Guess'
    SOFT_TOUCH dirac_mo_coef dirac_mo_label
    print*, 'HCore Guess'
   else
    print *,  'Unrecognized dirac_guess_type : '// dirac_mo_guess_type
    stop 1
   endif
  endif              
 end

 subroutine run_dirac
  BEGIN_DOC
  ! Run SCF_dirac calculation
  END_DOC
  use bitmasks
  implicit none
  call damping_Dirac_SCF  
 end

 subroutine print_dirac_energies
  BEGIN_DOC
  ! Print dirac bi-electronic energies 
  END_DOC
  use bitmasks
  implicit none
  PROVIDE ezfio_filename
 if (dirac_interaction == "Coulomb") then
  dirac_C_Hartree_Energy = dirac_HF_two_electron_C_Hartree_energy 
  call ezfio_set_dirac_scf_utils_dirac_c_hartree_energy(dirac_C_Hartree_Energy)
  dirac_C_Exchange_Energy = dirac_HF_two_electron_C_Exchange_energy
  call ezfio_set_dirac_scf_utils_dirac_c_exchange_energy(dirac_C_Exchange_Energy)
 elseif (dirac_interaction == "Coulomb_Gaunt") then
  dirac_C_Hartree_Energy = dirac_HF_two_electron_C_Hartree_energy
  call ezfio_set_dirac_scf_utils_dirac_c_hartree_energy(dirac_C_Hartree_Energy)
  dirac_C_Exchange_Energy = dirac_HF_two_electron_C_Exchange_energy
  call ezfio_set_dirac_scf_utils_dirac_c_exchange_energy(dirac_C_Exchange_Energy) 
  dirac_G_Hartree_Energy = dirac_HF_two_electron_G_Hartree_energy
  call ezfio_set_dirac_scf_utils_dirac_g_hartree_energy(dirac_G_Hartree_Energy)
  dirac_G_Exchange_Energy = dirac_HF_two_electron_G_Exchange_energy
  call ezfio_set_dirac_scf_utils_dirac_g_exchange_energy(dirac_G_Exchange_Energy) 
  dirac_C_G_Hartree_Energy = dirac_HF_two_electron_C_G_Hartree_energy
  call ezfio_set_dirac_scf_utils_dirac_c_g_hartree_energy(dirac_C_G_Hartree_Energy)
  dirac_C_G_Exchange_Energy = dirac_HF_two_electron_C_G_Exchange_energy
  call ezfio_set_dirac_scf_utils_dirac_c_g_exchange_energy(dirac_C_G_Exchange_Energy) 
 else
  print *,  'Unrecognized dirac_interaction : '//dirac_interaction
  stop 1
 endif
 end

 subroutine print_dirac_mo_coef
  BEGIN_DOC 
  ! Print dirac_mo_coef in dirac_mo_coef_Re and dirac_mo_coef_Im
  END_DOC
  use bitmasks
  implicit none
  integer :: i,j
  complex*16 :: dirac_mo_coef_temp
  PROVIDE ezfio_filename
  dirac_mo_coef_temp = dirac_mo_coef(1,2*small_ao_num+1)
  do j= 1, 2*dirac_ao_num
   do i= 1, 2*dirac_ao_num
    dirac_mo_coef(i,j) =dirac_mo_coef(i,j)*(Abs(dirac_mo_coef_temp)/dirac_mo_coef_temp)
   enddo
  enddo
  dirac_mo_coef_Re = Real(dirac_mo_coef)
  dirac_mo_coef_Im = Aimag(dirac_mo_coef)
  call ezfio_set_dirac_mo_basis_dirac_mo_coef_Re(dirac_mo_coef_Re)
  call ezfio_set_dirac_mo_basis_dirac_mo_coef_Im(dirac_mo_coef_Im)
 end

