 program scf_dirac
  BEGIN_DOC
  ! Produce `Hartree_Fock` MO orbital 
  ! output: mo_basis.mo_tot_num mo_basis.mo_label mo_basis.ao_md5 mo_basis.mo_coef mo_basis.mo_occ
  ! output: hartree_fock.energy
  ! optional: mo_basis.mo_coef
  END_DOC
  call create_dirac_guess
  call orthonormalize_dirac_mos
  call run_dirac
 end


 subroutine create_dirac_guess
  implicit none
  BEGIN_DOC
  ! Create a MO guess if no MOs are present in the EZFIO directory
  END_DOC
  logical                        :: exists
  PROVIDE ezfio_filename
 !call ezfio_has_mo_basis_mo_coef(exists)
 !if (.not.exists) then
    if (mo_guess_type == "HCore") then
      dirac_mo_coef = dirac_ao_ortho_canonical_coef
      TOUCH dirac_mo_coef
      dirac_mo_label = 'Guess'
      call dirac_mo_as_eigvectors_of_dirac_mo_matrix(dirac_mo_mono_elec_integral,size(dirac_mo_mono_elec_integral,1),size(dirac_mo_mono_elec_integral,2),dirac_mo_label,1,.false.)
      SOFT_TOUCH dirac_mo_coef dirac_mo_label
    else if (mo_guess_type == "Huckel") then
      call dirac_huckel_guess
    else
      print *,  'Unrecognized MO guess type : '//mo_guess_type
      stop 1
    endif
 !endif
 end

 subroutine run_dirac
  BEGIN_DOC
  ! Run SCF_dirac calculation
  END_DOC
  use bitmasks
  implicit none
  dirac_mo_label = "Canonical"
  soft_touch dirac_mo_label
 !Choose SCF algorithm
  call damping_dirac_SCF  
 !call Roothaan_Hall_SCF
 end


