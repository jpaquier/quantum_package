program scf
  BEGIN_DOC
! Produce `Kohn_Sham` MO orbital 
! output: mo_basis.mo_tot_num mo_basis.mo_label mo_basis.ao_md5 mo_basis.mo_coef mo_basis.mo_occ
! output: kohn_sham.energy
! optional: mo_basis.mo_coef
  END_DOC
  read_wf = .False.
  touch read_wf
  call check_coherence_functional
  call create_guess
  call orthonormalize_mos
  call run
end

subroutine check_coherence_functional
 implicit none
 integer :: ifound_x,ifound_c
 ifound_x = index(exchange_functional,"short_range")
 ifound_c = index(correlation_functional,"short_range")
 if(ifound_x .ne.0 .or. ifound_c .ne. 0)then
  print*,'YOU ARE USING THE REGULAR KS PROGRAM BUT YOUR INPUT KEYWORD FOR '
  print*,'exchange_functional is ',exchange_functional
  print*,'correlation_functional is ',correlation_functional
  print*,'CHANGE THE exchange_functional and correlation_functional keywords to regular functionals'
  print*,'or switch to the RS_KS_SCF program that uses range separated functionals'
  stop
 endif

end

subroutine create_guess
  implicit none
  BEGIN_DOC
!   Create a MO guess if no MOs are present in the EZFIO directory
  END_DOC
  logical                        :: exists
  PROVIDE ezfio_filename
  call ezfio_has_mo_basis_mo_coef(exists)
  if (.not.exists) then
    if (mo_guess_type == "HCore") then
      mo_coef = ao_ortho_lowdin_coef
      TOUCH mo_coef
      mo_label = 'Guess'
      call mo_as_eigvectors_of_mo_matrix(mo_mono_elec_integral,size(mo_mono_elec_integral,1),size(mo_mono_elec_integral,2),mo_label,.false.)
      SOFT_TOUCH mo_coef mo_label
    else if (mo_guess_type == "Huckel") then
      call huckel_guess
    else
      print *,  'Unrecognized MO guess type : '//mo_guess_type
      stop 1
    endif
  endif
end

subroutine run

  BEGIN_DOC
!   Run SCF calculation
  END_DOC

  use bitmasks
  implicit none

  double precision               :: EHF
   
  EHF = KS_energy 

  mo_label = "Canonical"

! Choose SCF algorithm

!    call damping_SCF   ! Deprecated routine
  call Roothaan_Hall_SCF
  
end


