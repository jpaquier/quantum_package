 program dirac_scf
  BEGIN_DOC
  ! Produce `Dirac_Hartree_Fock` MO orbital 
  ! output: none for the moment
  END_DOC
  call create_dirac_guess
  call run_dirac
  call print_dirac_mo_coef
 end

 subroutine create_dirac_guess
  implicit none
  BEGIN_DOC
  ! Create a MO guess if no MOs are present in the EZFIO directory
  END_DOC
 logical                        :: exists_Re, exists_Im
  PROVIDE ezfio_filename
  call ezfio_has_dirac_scf_dirac_mo_coef_Re(exists_Re)
  call ezfio_has_dirac_scf_dirac_mo_coef_Im(exists_Im)
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
    print *,  'Unrecognized MO guess type : '//mo_guess_type
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

 subroutine print_dirac_mo_coef
  BEGIN_DOC 
  ! Print dirac_mo_coef in dirac_mo_coef_Re and dirac_mo_coef_Im
  END_DOC
  use bitmasks
  implicit none
  PROVIDE ezfio_filename
  dirac_mo_coef_Re = Real(dirac_mo_coef)
  dirac_mo_coef_Im = Aimag(dirac_mo_coef)
  call ezfio_set_dirac_scf_dirac_mo_coef_Re(dirac_mo_coef_Re)
  call ezfio_set_dirac_scf_dirac_mo_coef_Im(dirac_mo_coef_Im)
 end

