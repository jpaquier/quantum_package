!program scf_dirac
! BEGIN_DOC
! ! Produce `Dirac-Hartree_Fock` MO orbital 
! ! output: none for the moment
! END_DOC
! call check_range_separation
! call print_dirac_mo_coef
!!call test_dirac_end
!end

!subroutine check_range_separation
! implicit none
! BEGIN_DOC
! ! Checks for the use of fullrange or long-range Hartree-Fock
! END_DOC
! integer :: ifound_x,ifound_c
! if (dirac_range_separation == "Full_range") then
!  print *, 'DHF calculation'
!  call create_dirac_guess
!! call test_dirac
!  call run_dirac
! elseif (dirac_range_separation == "Long_range") then
!  print *, '          '
!  print *, 'Range-separated RDFT calculation'
!  print *, '**************************'
!  print *, 'mu_erf = ',mu_erf
!  print *, '**************************'
!  if (dirac_exchange_functional == "None") then
!   ifound_x = 1
!  else
!   ifound_x = index(dirac_exchange_functional,"short_range")
!  endif
!  if (dirac_correlation_functional == "None") then
!   ifound_c = 1
!  else
!   ifound_c = index(dirac_correlation_functional,"short_range")
!  endif
!! print*,ifound_x,ifound_c
!  if(ifound_x ==0 .or. ifound_c == 0)then
!   print *,'YOU ARE USING THE RELATIVISTIC RANGE SEPARATED KS PROGRAM BUT YOUR INPUT KEYWORD FOR '
!   print *,'dirac_exchange_functional is ', dirac_exchange_functional
!   print *,'dirac_correlation_functional is ', dirac_correlation_functional
!   print *,'CHANGE THE dirac_exchange_functional AND dirac_correlation_functional KEYWORDS TO RELATIVISTIC RANGE SEPARATED FUNCTIONALS'
!   stop
!  else
!   print *,'dirac_exchange_functional =',dirac_exchange_functional
!   print *,'dirac_correlation_functional =',dirac_correlation_functional
!   call create_dirac_guess
!   call test_dirac
!   call run_dirac_erf
!  endif
! else
!  print *,  'Unrecognized dirac_range_separation : '//dirac_range_separation
!  stop 1
! endif 
!end
!

!subroutine create_dirac_guess
! implicit none
! BEGIN_DOC
! ! Create a MO guess if no MOs are present in the EZFIO directory
! END_DOC
!logical                        :: exists_Re, exists_Im
! PROVIDE ezfio_filename
! call ezfio_has_dirac_scf_dirac_mo_coef_Re(exists_Re)
! call ezfio_has_dirac_scf_dirac_mo_coef_Im(exists_Im)
! if (exists_Re .and. exists_Im) then
!  dirac_mo_coef = (1,0)*dirac_mo_coef_Re + (0,1)*dirac_mo_coef_Im
!  print*, 'dirac_mo_coef from previous calculation'
! else 
!  if (dirac_mo_guess_type == "HCore") then
!   dirac_mo_coef = eigenvectors_dirac_mono_elec_ao
!   TOUCH dirac_mo_coef       
!   dirac_mo_label = 'Guess'
!   SOFT_TOUCH dirac_mo_coef dirac_mo_label
!   print*, 'HCore Guess'
!  else
!   print *,  'Unrecognized MO guess type : '//mo_guess_type
!   stop 1
!  endif
! endif              
!end

!subroutine test_dirac
! implicit none
! integer   :: i,j
! complex*16 :: ortho(2*dirac_ao_num)
! ortho = (0.d0,0.d0)
! do j= 2*small_ao_num+2, 2*small_ao_num+2
!  do i= 1, 2*dirac_ao_num
!  !print*, i, j, dirac_mo_coef_guess(i,j)
!  !print*,i,j,dirac_mo_coef(i,j)
!  !print*,i,j,dirac_mo_coef_S(i,j) 
!  !print*,i,j,dirac_mo_overlap(i,j) 
!  !print*,i,j,dirac_mo_overlap_bis(i,j) 
!  !print*,i,j,dirac_SCF_density_matrix_ao(i+large_ao_num,j)
!  !print*,i,j,dirac_SCF_density_matrix_ao(j,i+large_ao_num)
!  !print*,i,j,dirac_ao_bi_elec_integral(i,j)
!  !print*,i,j,dirac_ao_bi_elec_integral_qq(i,j)
!  !print*,i,j,dirac_ao_mono_elec_integral (i,j)
!  !print*,i,j,dirac_Fock_matrix_ao (i,j)
!  !ortho(j) += Conjg(dirac_mo_coef(i,2))*dirac_mo_coef(i,j)
!  enddo
!  print*,'********************'     
! !print*,j,ortho(j)
! enddo
!!print*, nuclear_repulsion
!!print*, dirac_HF_one_electron_energy
!!print*, dirac_HF_two_electron_energy
! print*,'***********'
!!print*, dirac_SCF_energy 
!!do i = 1, 2*dirac_ao_num
!! print*,i,eigenvalues_dirac_fock_matrix_mo(i)
!!!print*,i,eigenvalues_dirac_mono_elec_mo(i)
!!!print*,'*********'
!!enddo
!end



!subroutine run_dirac
! BEGIN_DOC
! ! Run SCF_dirac calculation
! END_DOC
! use bitmasks
! implicit none
!!dirac_mo_label = "Canonical"
!!soft_touch dirac_mo_label
!!Choose SCF algorithm
!!if (dirac_interaction == "Coulomb") then
!  call damping_Dirac_SCF  
! !call Roothaan_Hall_Dirac_SCF
!!elseif (dirac_interaction == "Coulomb_Gaunt") then
! !call damping_Dirac_Gaunt_SCF
!!else
!! print *,  'Unrecognized dirac_interaction : '//dirac_interaction
!! stop 1
!!endif
!end

!subroutine run_dirac_erf
! BEGIN_DOC
! ! Run SCF_dirac calculation
! END_DOC
! use bitmasks
! implicit none
!!dirac_mo_label = "Canonical"
!!soft_touch dirac_mo_label
!!Choose SCF algorithm
! if (dirac_interaction == "Coulomb") then
!! call damping_Dirac_erf_SCF  
! elseif (dirac_interaction == "Coulomb_Gaunt") then
!! call damping_Dirac_Gaunt_erf_SCF
! else
!  print *,  'Unrecognized dirac_interaction : '//dirac_interaction
!  stop 1
! endif
!end
!
!subroutine print_dirac_mo_coef
! BEGIN_DOC 
! ! Print dirac_mo_coef in dirac_mo_coef_Re and dirac_mo_coef_Im
! END_DOC
! use bitmasks
! implicit none
! PROVIDE ezfio_filename
! dirac_mo_coef_Re = Real(dirac_mo_coef)
! dirac_mo_coef_Im = Aimag(dirac_mo_coef)
! call ezfio_set_dirac_scf_dirac_mo_coef_Re(dirac_mo_coef_Re)
! call ezfio_set_dirac_scf_dirac_mo_coef_Im(dirac_mo_coef_Im)
!end

!subroutine test_dirac_end
! implicit none
! integer   :: i,j
! complex*16 :: ortho(2*dirac_ao_num)
! ortho = (0.d0,0.d0)
! do j= 2*small_ao_num+2, 2*small_ao_num+2
!  do i= 1,2*dirac_ao_num
!  !print*, i, j, dirac_mo_coef_guess(i,j)
!   print*,i,j,dirac_mo_coef(i,j)
!  !print*,i,j,dirac_mo_coef_S(i,j) 
!  !print*,i,j,dirac_mo_overlap(i,j) 
!  !print*,i,j,dirac_mo_overlap_bis(i,j) 
!  !print*,i,j,dirac_SCF_density_matrix_ao(i+large_ao_num,j)
!  !print*,i,j,dirac_SCF_density_matrix_ao(j,i+large_ao_num)
!  !print*,i,j,dirac_ao_bi_elec_integral(i,j)
!  !print*,i,j,dirac_ao_bi_elec_integral_qq(i,j)
!  !print*,i,j,dirac_ao_mono_elec_integral (i,j)
!  !print*,i,j,dirac_Fock_matrix_ao (i,j)
!  !ortho(j) += Conjg(dirac_mo_coef(i,2))*dirac_mo_coef(i,j)
!  enddo
! print*,'********************'     
!! print*,j,ortho(j)
! enddo
!!print*, nuclear_repulsion
!!print*, dirac_HF_one_electron_energy
!!print*, dirac_HF_two_electron_energy
!!print*,'***********'
!!print*, dirac_SCF_energy
! do i = 1, 2*dirac_ao_num
! !print*,i,eigenvalues_dirac_fock_matrix_mo(i)
! !print*,i,eigenvalues_dirac_mono_elec_mo(i)
! !print*,'*********'
! enddo
!end

