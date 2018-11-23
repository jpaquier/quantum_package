!subroutine damping_dirac_SCF_old
! implicit none
! double precision               :: E, E_min, delta_E, E_half, E_new, a, b, lambda,delta_D
! complex*16, allocatable        :: D(:,:)
! complex*16, allocatable        :: D_new(:,:), F_new(:,:)
! complex*16, allocatable        :: delta(:,:)
! complex*16                     :: delta_D_complex
! integer                        :: i,j,k
! logical                        :: saving
! character                      :: save_char
! allocate(D( 2*dirac_ao_num, 2*dirac_ao_num ),                           &
!          F_new( 2*dirac_ao_num, 2*dirac_ao_num ),                       &
!          D_new( 2*dirac_ao_num, 2*dirac_ao_num ),                       &
!          delta( 2*dirac_ao_num, 2*dirac_ao_num ))
! do j=1,2*dirac_ao_num
!  do i=1,2*dirac_ao_num
!   D(i,j) = dirac_SCF_density_matrix_ao(i,j)
!  enddo
! enddo
! call write_time(6)
! write(6,'(A4,1X,A16, 1X, A16, 1X, A16, 1X, A4 )')  &
!   '====','================','================','================', '===='
! write(6,'(A4,1X,A16, 1X, A16, 1X, A16, 1X, A4 )')  &
!   '  N ', 'Energy  ', 'Energy diff  ', 'Density diff  ', 'Save'
! write(6,'(A4,1X,A16, 1X, A16, 1X, A16, 1X, A4 )')  &
!   '====','================','================','================', '===='
! E = dirac_SCF_energy + 1.d0
! E_min = dirac_SCF_energy
! delta_D_complex = (0.d0,0.d0)
! delta_D = 0.d0
! do k=1,n_it_scf_max
!  delta_E = dirac_SCF_energy - E
!  E = dirac_SCF_energy
!  if (dabs(delta_E) < thresh_scf)  then
!   exit
!  endif
!  saving = E < E_min
!  if (saving) then
! ! call save_mos
!   save_char = 'X'
!   E_min = E
!  else
!    save_char = ' '
!  endif
!  write(6,'(I4, 1X, F16.8, 1X, F16.8, 1X, F16.8, 3X, A4 )')  &
!   k, dirac_SCF_energy, delta_E, delta_D, save_char
!  D = dirac_SCF_density_matrix_ao
!  dirac_mo_coef = eigenvectors_dirac_fock_matrix_ao
!  TOUCH dirac_mo_coef
!  D_new = dirac_SCF_density_matrix_ao
!  F_new = dirac_Fock_matrix_ao
!  E_new = dirac_SCF_energy
!  delta = D_new - D
!  lambda = .5d0
!  E_half = 0.d0
!  do while (E_half > E)
!   dirac_SCF_density_matrix_ao = D + lambda * delta
!   TOUCH dirac_SCF_density_matrix_ao
!   dirac_mo_coef = eigenvectors_dirac_fock_matrix_ao
!   TOUCH dirac_mo_coef
!   E_half = dirac_SCF_energy
!   if ((E_half > E).and.(E_new < E)) then
!    lambda = 1.d0
!    exit
!   else if ((E_half > E).and.(lambda > 5.d-4)) then
!    lambda = 0.5d0 * lambda
!    E_new = E_half
!   else
!    exit
!   endif
!  enddo
!  a = (E_new + E - 2.d0*E_half)*2.d0
!  b = -E_new - 3.d0*E + 4.d0*E_half
!  lambda = -lambda*b/(a+1.d-16)
!  D = (1.d0-lambda) * D + lambda * D_new
!  delta_E = dirac_SCF_energy - E
!  do j=1,2*dirac_ao_num
!   do i=1,2*dirac_ao_num
!    delta_D_complex += D(i,j) - dirac_SCF_density_matrix_ao(i,j)
!   enddo
!  enddo 
!  delta_D_complex = delta_D_complex/dble(2*dirac_ao_num)
!  delta_D = real(delta_D_complex)
!  if (aimag(delta_D_complex) .gt. 1.d-10) then  
!   print*, 'Warning! The delta_D is complex'
!   print*, 'delta_D_complex =', delta_D_complex
!   STOP
!  endif
!  dirac_SCF_density_matrix_ao = D
!  TOUCH dirac_SCF_density_matrix_ao
!  dirac_mo_coef = eigenvectors_dirac_fock_matrix_ao
!  TOUCH dirac_mo_coef
! enddo
! write(6,'(A4,1X,A16, 1X, A16, 1X, A16, 1X, A4 )')  '====','================','================','================', '===='
! write(6,*)
! call write_double(6, dirac_SCF_energy, 'Coulomb Hartree-Fock energy')
!!call ezfio_set_hartree_fock_energy(E_min)
! call write_time(6)
! deallocate(D,F_new,D_new,delta)
!end
