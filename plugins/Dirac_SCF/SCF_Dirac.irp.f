 program scf_dirac
  BEGIN_DOC
  ! Produce `Dirac-Hartree_Fock` MO orbital 
  ! output: none for the moment
  END_DOC
  call create_dirac_guess
 !call test_dirac
  call run_dirac
 !call test_dirac_end
 end

 subroutine create_dirac_guess
  implicit none
  BEGIN_DOC
  ! Create a MO guess if no MOs are present in the EZFIO directory
  END_DOC
 logical                        :: exists
  PROVIDE ezfio_filename
  call ezfio_has_mo_basis_mo_coef(exists)
 !if (exists) then
 ! dirac_mo_coef = dirac_mo_coef_guess
 ! print*, 'Non-relativistic Guess'
 !elseif (.not.exists) then
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
 !endif              
 end

 subroutine test_dirac
  implicit none
  integer   :: i,j
  complex*16 :: ortho(2*dirac_ao_num)
  ortho = (0.d0,0.d0)
  do j= 1, 2*dirac_ao_num
   do i=1, 2*dirac_ao_num
   !print*, i, j, dirac_mo_coef_guess(i,j)
   !print*,i,j,dirac_mo_coef(i,j)
   !print*,i,j,dirac_mo_coef_S(i,j) 
   !print*,i,j,dirac_mo_overlap(i,j) 
   !print*,i,j,dirac_mo_overlap_bis(i,j) 
   !print*,i,j,dirac_SCF_density_matrix_ao(i+large_ao_num,j)
   !print*,i,j,dirac_SCF_density_matrix_ao(j,i+large_ao_num)
   !print*,i,j,dirac_ao_bi_elec_integral(i,j)
   !print*,i,j,dirac_ao_bi_elec_integral_qq(i,j)
   !print*,i,j,dirac_ao_mono_elec_integral (i,j)
   !print*,i,j,dirac_Fock_matrix_ao (i,j)
    ortho(j) += Conjg(dirac_mo_coef(i,2))*dirac_mo_coef(i,j)
   enddo
  !print*,'********************'     
  !print*,j,ortho(j)
  enddo
  print*, nuclear_repulsion
  print*, dirac_HF_one_electron_energy
  print*, dirac_HF_two_electron_energy
  print*,'***********'
  print*, dirac_SCF_energy 
 !do i = 1, 2*dirac_ao_num
 ! print*,i,eigenvalues_dirac_fock_matrix_mo(i)
 !!print*,i,eigenvalues_dirac_mono_elec_mo(i)
 !!print*,'*********'
 !enddo
 end



 subroutine run_dirac
  BEGIN_DOC
  ! Run SCF_dirac calculation
  END_DOC
  use bitmasks
  implicit none
 !dirac_mo_label = "Canonical"
 !soft_touch dirac_mo_label
 !Choose SCF algorithm
  if (dirac_interaction == "Coulomb") then
   call damping_Dirac_SCF  
  !call Roothaan_Hall_Dirac_SCF
  elseif (dirac_interaction == "Coulomb_Gaunt") then
   call damping_Dirac_Gaunt_SCF
  else
   print *,  'Unrecognized dirac_interaction : '//dirac_interaction
   stop 1
  endif
 end

 subroutine test_dirac_end
  implicit none
  integer   :: i,j
  complex*16 :: ortho(2*dirac_ao_num)
  ortho = (0.d0,0.d0)
  do j= 1, 2*dirac_ao_num
   do i=1, 2*dirac_ao_num
   !print*, i, j, dirac_mo_coef_guess(i,j)
   !print*,i,j,dirac_mo_coef(i,j)
   !print*,i,j,dirac_mo_coef_S(i,j) 
   !print*,i,j,dirac_mo_overlap(i,j) 
   !print*,i,j,dirac_mo_overlap_bis(i,j) 
   !print*,i,j,dirac_SCF_density_matrix_ao(i+large_ao_num,j)
   !print*,i,j,dirac_SCF_density_matrix_ao(j,i+large_ao_num)
   !print*,i,j,dirac_ao_bi_elec_integral(i,j)
   !print*,i,j,dirac_ao_bi_elec_integral_qq(i,j)
   !print*,i,j,dirac_ao_mono_elec_integral (i,j)
   !print*,i,j,dirac_Fock_matrix_ao (i,j)
    ortho(j) += Conjg(dirac_mo_coef(i,2))*dirac_mo_coef(i,j)
   enddo
  !print*,'********************'     
 ! print*,j,ortho(j)
  enddo
  print*, nuclear_repulsion
  print*, dirac_HF_one_electron_energy
  print*, dirac_HF_two_electron_energy
  print*,'***********'
  print*, dirac_SCF_energy
  do i = 1, 2*dirac_ao_num
   print*,i,eigenvalues_dirac_fock_matrix_mo(i)
  !print*,i,eigenvalues_dirac_mono_elec_mo(i)
  !print*,'*********'
  enddo
 end

