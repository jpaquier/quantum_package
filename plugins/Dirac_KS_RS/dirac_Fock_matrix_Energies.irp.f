 BEGIN_PROVIDER [ complex*16, dirac_ao_bi_elec_C_Hartree_integral, (2*dirac_ao_num, 2*dirac_ao_num) ]
  use map_module
  implicit none
  BEGIN_DOC
  !Array of the bi-electronic Fock matrix for the Coulomb interaction
  ! in Dirac AO basis set
  !Take care, the density matrix index
  ! have been correctly inverse, unlike
  ! in the non-relativistic code
  END_DOC
  PROVIDE dirac_ao_bielec_integrals_erf_in_map
  integer                        :: i,j,k,l,k1 
  integer                        :: r,s,p,q
  double precision               :: dirac_ao_bielec_integral, local_threshold
  double precision               :: integral
  complex*16, allocatable        :: dirac_ao_bi_elec_C_Hartree_integral_tmp(:,:)
  integer(omp_lock_kind)         :: lck(dirac_ao_num)
  integer*8                      :: i8
  integer                        :: ii(8), jj(8), kk(8), ll(8), k2
  integer(cache_map_size_kind)   :: n_elements_max, n_elements
  integer(key_kind), allocatable :: keys(:)
  double precision, allocatable  :: values(:)
  dirac_ao_bi_elec_C_Hartree_integral = (0.d0,0.d0)
!!$OMP PARALLEL DEFAULT(NONE)                                      &
!!$OMP PRIVATE(i,j,l,k1,k,integral,ii,jj,kk,ll,i8,keys,values,n_elements_max, &
!!$OMP  n_elements,dirac_ao_bi_elec_C_Hartree_integral_tmp)&
!!$OMP SHARED(dirac_ao_num,dirac_SCF_density_matrix_ao,&
!!$OMP  dirac_ao_integrals_map, dirac_ao_bi_elec_C_Hartree_integral) 
  call get_cache_map_n_elements_max(dirac_ao_integrals_erf_map,n_elements_max)
  allocate(keys(n_elements_max), values(n_elements_max))
  allocate(dirac_ao_bi_elec_C_Hartree_integral_tmp(2*dirac_ao_num,2*dirac_ao_num))
  dirac_ao_bi_elec_C_Hartree_integral_tmp = (0.d0,0.d0)
!!$OMP DO SCHEDULE(dynamic,64)
!!DIR$ NOVECTOR
  do i8=0_8,dirac_ao_integrals_erf_map%map_size
   n_elements = n_elements_max
   call get_cache_map(dirac_ao_integrals_erf_map,i8,keys,values,n_elements)
   do k1=1,n_elements
    call bielec_integrals_index_reverse(kk,ii,ll,jj,keys(k1))
    do k2=1,8
     if (kk(k2)==0) then
      cycle
     endif
     i = ii(k2) ! electron 1
     j = jj(k2) ! electron 1
     k = kk(k2) ! electron 2
     l = ll(k2) ! electron 2
     ! values(k1) = (ij|kl) <=> <ik|jl>
     integral = values (k1)
     if ((i .le. large_ao_num .and. j .le. large_ao_num .and. k .le. large_ao_num .and. l .le. large_ao_num) .or.  &
         (i .gt. large_ao_num .and. j .gt. large_ao_num .and. k .gt. large_ao_num .and. l .gt. large_ao_num)) then
      !L_alpha L_alpha .or. S_alpha S_alpha
      dirac_ao_bi_elec_C_Hartree_integral_tmp(dirac_inverse_list(i,1),dirac_inverse_list(j,1)) += (dirac_SCF_density_matrix_ao(dirac_inverse_list(l,1),dirac_inverse_list(k,1))) * integral 
      dirac_ao_bi_elec_C_Hartree_integral_tmp(dirac_inverse_list(i,1),dirac_inverse_list(j,1)) += (dirac_SCF_density_matrix_ao(dirac_inverse_list(l,2),dirac_inverse_list(k,2))) * integral
      !L_beta L_beta .or. S_beta S_beta
      dirac_ao_bi_elec_C_Hartree_integral_tmp(dirac_inverse_list(i,2),dirac_inverse_list(j,2)) += (dirac_SCF_density_matrix_ao(dirac_inverse_list(l,1),dirac_inverse_list(k,1))) * integral
      dirac_ao_bi_elec_C_Hartree_integral_tmp(dirac_inverse_list(i,2),dirac_inverse_list(j,2)) += (dirac_SCF_density_matrix_ao(dirac_inverse_list(l,2),dirac_inverse_list(k,2))) * integral
     elseif((i .le. large_ao_num .and. j .le. large_ao_num .and. k .gt. large_ao_num .and. l .gt. large_ao_num) .or. &
            (i .gt. large_ao_num .and. j .gt. large_ao_num .and. k .le. large_ao_num .and. l .le. large_ao_num))then
       !L_alpha L_alpha .or. S_alpha S_alpha
      dirac_ao_bi_elec_C_Hartree_integral_tmp(dirac_inverse_list(i,1),dirac_inverse_list(j,1)) += (dirac_SCF_density_matrix_ao(dirac_inverse_list(l,1),dirac_inverse_list(k,1))) * integral
      dirac_ao_bi_elec_C_Hartree_integral_tmp(dirac_inverse_list(i,1),dirac_inverse_list(j,1)) += (dirac_SCF_density_matrix_ao(dirac_inverse_list(l,2),dirac_inverse_list(k,2))) * integral
      !L_beta L_beta .or. S_beta S_beta
      dirac_ao_bi_elec_C_Hartree_integral_tmp(dirac_inverse_list(i,2),dirac_inverse_list(j,2)) += (dirac_SCF_density_matrix_ao(dirac_inverse_list(l,1),dirac_inverse_list(k,1))) * integral
      dirac_ao_bi_elec_C_Hartree_integral_tmp(dirac_inverse_list(i,2),dirac_inverse_list(j,2)) += (dirac_SCF_density_matrix_ao(dirac_inverse_list(l,2),dirac_inverse_list(k,2))) * integral
     endif
    enddo
   enddo
  enddo
!!$OMP END DO NOWAIT
!!$OMP CRITICAL
  dirac_ao_bi_elec_C_Hartree_integral += dirac_ao_bi_elec_C_Hartree_integral_tmp
!!$OMP END CRITICAL
  deallocate(keys,values,dirac_ao_bi_elec_C_Hartree_integral_tmp)
!!$OMP END PARALLEL
 END_PROVIDER


 BEGIN_PROVIDER [ complex*16, dirac_HF_two_electron_C_Hartree_energy_complex]
 &BEGIN_PROVIDER [ double precision, dirac_HF_two_electron_C_Hartree_energy] 
  implicit none
  BEGIN_DOC
  !Two-electrons energy of the Coulomb ee interaction
  !The energy is supposed to be a real, thus we check for its complex part to be
  ! a VERY small artifact and take only its real part
  END_DOC
  integer :: i,j
  dirac_HF_two_electron_C_Hartree_energy_complex = (0.d0,0.d0)
  do j=1, 2*dirac_ao_num
   do i=1, 2*dirac_ao_num
    dirac_HF_two_electron_C_Hartree_energy_complex += 0.5d0* (dirac_ao_bi_elec_C_Hartree_integral(i,j)) * dirac_SCF_density_matrix_ao(j,i)
   enddo
  enddo
  dirac_HF_two_electron_C_Hartree_energy = real(dirac_HF_two_electron_C_Hartree_energy_complex)
  if (aimag(dirac_HF_two_electron_C_Hartree_energy_complex) .gt. 1.d-10) then
  print*, 'Warning! The energy is not real'
  print*, 'dirac_HF_two_electron_C_Hartree_energy_complex =', dirac_HF_two_electron_C_Hartree_energy_complex
  STOP
  endif
 END_PROVIDER


 BEGIN_PROVIDER [ complex*16, dirac_ao_bi_elec_C_Exchange_integral, (2*dirac_ao_num, 2*dirac_ao_num) ]
  use map_module
  implicit none
  BEGIN_DOC
  !Array of the bi-electronic Fock matrix for the Coulomb interaction
  ! in Dirac AO basis set
  !Take care, the density matrix index
  ! have been correctly inverse, unlike
  ! in the non-relativistic code
  END_DOC
  PROVIDE dirac_ao_bielec_integrals_erf_in_map
  integer                        :: i,j,k,l,k1 
  integer                        :: r,s,p,q
  double precision               :: dirac_ao_bielec_integral, local_threshold
  double precision               :: integral
  complex*16, allocatable        :: dirac_ao_bi_elec_C_Exchange_integral_tmp(:,:)
  integer(omp_lock_kind)         :: lck(dirac_ao_num)
  integer*8                      :: i8
  integer                        :: ii(8), jj(8), kk(8), ll(8), k2
  integer(cache_map_size_kind)   :: n_elements_max, n_elements
  integer(key_kind), allocatable :: keys(:)
  double precision, allocatable  :: values(:)
  dirac_ao_bi_elec_C_Exchange_integral = (0.d0,0.d0)
!!$OMP PARALLEL DEFAULT(NONE)                                      &
!!$OMP PRIVATE(i,j,l,k1,k,integral,ii,jj,kk,ll,i8,keys,values,n_elements_max, &
!!$OMP  n_elements,dirac_ao_bi_elec_C_Exchange_integral_tmp)&
!!$OMP SHARED(dirac_ao_num,dirac_SCF_density_matrix_ao,&
!!$OMP  dirac_ao_integrals_map, dirac_ao_bi_elec_C_Exchange_integral) 
  call get_cache_map_n_elements_max(dirac_ao_integrals_erf_map,n_elements_max)
  allocate(keys(n_elements_max), values(n_elements_max))
  allocate(dirac_ao_bi_elec_C_Exchange_integral_tmp(2*dirac_ao_num,2*dirac_ao_num))
  dirac_ao_bi_elec_C_Exchange_integral_tmp = (0.d0,0.d0)
!!$OMP DO SCHEDULE(dynamic,64)
!!DIR$ NOVECTOR
  do i8=0_8,dirac_ao_integrals_erf_map%map_size
   n_elements = n_elements_max
   call get_cache_map(dirac_ao_integrals_erf_map,i8,keys,values,n_elements)
   do k1=1,n_elements
    call bielec_integrals_index_reverse(kk,ii,ll,jj,keys(k1))
    do k2=1,8
     if (kk(k2)==0) then
      cycle
     endif
     i = ii(k2) ! electron 1
     j = jj(k2) ! electron 1
     k = kk(k2) ! electron 2
     l = ll(k2) ! electron 2
     ! values(k1) = (ij|kl) <=> <ik|jl>
     integral = values (k1)
     if ((i .le. large_ao_num .and. j .le. large_ao_num .and. k .le. large_ao_num .and. l .le. large_ao_num) .or.  &
         (i .gt. large_ao_num .and. j .gt. large_ao_num .and. k .gt. large_ao_num .and. l .gt. large_ao_num)) then
      !L_alpha L_alpha .or. S_alpha S_alpha
      dirac_ao_bi_elec_C_Exchange_integral_tmp(dirac_inverse_list(i,1),dirac_inverse_list(l,1)) -= (dirac_SCF_density_matrix_ao(dirac_inverse_list(j,1),dirac_inverse_list(k,1))) * integral
      !L_beta L_beta .or. S_beta S_beta
      dirac_ao_bi_elec_C_Exchange_integral_tmp(dirac_inverse_list(i,2),dirac_inverse_list(l,2)) -= (dirac_SCF_density_matrix_ao(dirac_inverse_list(j,2),dirac_inverse_list(k,2))) * integral   
      !L_beta L_alpha .or. S_beta S_alpha
      dirac_ao_bi_elec_C_Exchange_integral_tmp(dirac_inverse_list(i,2),dirac_inverse_list(l,1)) -= (dirac_SCF_density_matrix_ao(dirac_inverse_list(j,2),dirac_inverse_list(k,1))) * integral
      !L_alpha L_beta .or. S_alpha S_beta
      dirac_ao_bi_elec_C_Exchange_integral_tmp(dirac_inverse_list(i,1),dirac_inverse_list(l,2)) -= (dirac_SCF_density_matrix_ao(dirac_inverse_list(j,1),dirac_inverse_list(k,2))) * integral
     elseif((i .le. large_ao_num .and. j .le. large_ao_num .and. k .gt. large_ao_num .and. l .gt. large_ao_num) .or. &
            (i .gt. large_ao_num .and. j .gt. large_ao_num .and. k .le. large_ao_num .and. l .le. large_ao_num))then
      !L_alpha S_alpha .or S_alpha L_alpha
      dirac_ao_bi_elec_C_Exchange_integral_tmp(dirac_inverse_list(i,1),dirac_inverse_list(l,1)) -= (dirac_SCF_density_matrix_ao(dirac_inverse_list(j,1),dirac_inverse_list(k,1))) * integral
      !L_beta S_beta .or. S_beta L_beta
      dirac_ao_bi_elec_C_Exchange_integral_tmp(dirac_inverse_list(i,2),dirac_inverse_list(l,2)) -= (dirac_SCF_density_matrix_ao(dirac_inverse_list(j,2),dirac_inverse_list(k,2))) * integral
      !L_beta S_alpha .or. S_beta L_alpha
      dirac_ao_bi_elec_C_Exchange_integral_tmp(dirac_inverse_list(i,2),dirac_inverse_list(l,1)) -= (dirac_SCF_density_matrix_ao(dirac_inverse_list(j,2),dirac_inverse_list(k,1))) * integral
      !L_alpha S_beta .or. S_alpha L_beta
      dirac_ao_bi_elec_C_Exchange_integral_tmp(dirac_inverse_list(i,1),dirac_inverse_list(l,2)) -= (dirac_SCF_density_matrix_ao(dirac_inverse_list(j,1),dirac_inverse_list(k,2))) * integral
     endif
    enddo
   enddo
  enddo
!!$OMP END DO NOWAIT
!!$OMP CRITICAL
  dirac_ao_bi_elec_C_Exchange_integral += dirac_ao_bi_elec_C_Exchange_integral_tmp
!!$OMP END CRITICAL
  deallocate(keys,values,dirac_ao_bi_elec_C_Exchange_integral_tmp)
!!$OMP END PARALLEL
 END_PROVIDER

 BEGIN_PROVIDER [ complex*16, dirac_HF_two_electron_C_Exchange_energy_complex]
 &BEGIN_PROVIDER [ double precision, dirac_HF_two_electron_C_Exchange_energy] 
  implicit none
  BEGIN_DOC
  !Two-electrons energy of the Coulomb ee interaction
  !The energy is supposed to be a real, thus we check for its complex part to be
  ! a VERY small artifact and take only its real part
  END_DOC
  integer :: i,j
  dirac_HF_two_electron_C_Exchange_energy_complex = (0.d0,0.d0)
  do j=1, 2*dirac_ao_num
   do i=1, 2*dirac_ao_num
    dirac_HF_two_electron_C_Exchange_energy_complex += 0.5d0* (dirac_ao_bi_elec_C_Exchange_integral(i,j)) * dirac_SCF_density_matrix_ao(j,i)
   enddo
  enddo
  dirac_HF_two_electron_C_Exchange_energy = real(dirac_HF_two_electron_C_Exchange_energy_complex)
  if (aimag(dirac_HF_two_electron_C_Exchange_energy_complex) .gt. 1.d-10) then
  print*, 'Warning! The energy is not real'
  print*, 'dirac_HF_two_electron_C_Exchange_energy_complex =', dirac_HF_two_electron_C_Exchange_energy_complex
  STOP
  endif
 END_PROVIDER


 BEGIN_PROVIDER [ complex*16, dirac_ao_bi_elec_G_Hartree_integral, (2*dirac_ao_num, 2*dirac_ao_num) ]
  use map_module
  implicit none
  BEGIN_DOC
  !Array of the bi-electronic Fock matrix for the Coulomb-Gaunt interaction,
  ! in Dirac AO basis set.
  !Take care, the density matrix index
  ! have been correctly inverse, unlike
  ! in the non-relativistic code
  END_DOC
  PROVIDE dirac_ao_bielec_integrals_erf_in_map
  integer                        :: i,j,k,l,k1 
  integer                        :: r,s,p,q
  double precision               :: dirac_ao_bielec_integral, local_threshold
  double precision               :: integral
  complex*16, allocatable        :: dirac_ao_bi_elec_G_Hartree_integral_tmp(:,:)
  integer(omp_lock_kind)         :: lck(dirac_ao_num)
  integer*8                      :: i8
  integer                        :: ii(8), jj(8), kk(8), ll(8), k2
  integer(cache_map_size_kind)   :: n_elements_max, n_elements
  integer(key_kind), allocatable :: keys(:)
  double precision, allocatable  :: values(:)
  dirac_ao_bi_elec_G_Hartree_integral = (0.d0,0.d0)
!!$OMP PARALLEL DEFAULT(NONE)                                      &
!!$OMP PRIVATE(i,j,l,k1,k,integral,ii,jj,kk,ll,i8,keys,values,n_elements_max, &
!!$OMP  n_elements,dirac_ao_bi_elec_G_Hartree_integral_tmp)&
!!$OMP SHARED(dirac_ao_num,dirac_SCF_density_matrix_ao,&
!!$OMP  dirac_ao_integrals_map, dirac_ao_bi_elec_G_Hartree_integral) 
  call get_cache_map_n_elements_max(dirac_ao_integrals_erf_map,n_elements_max)
  allocate(keys(n_elements_max), values(n_elements_max))
  allocate(dirac_ao_bi_elec_G_Hartree_integral_tmp(2*dirac_ao_num,2*dirac_ao_num))
  dirac_ao_bi_elec_G_Hartree_integral_tmp = (0.d0,0.d0)
!!$OMP DO SCHEDULE(dynamic,64)
!!DIR$ NOVECTOR
  do i8=0_8,dirac_ao_integrals_erf_map%map_size
   n_elements = n_elements_max
   call get_cache_map(dirac_ao_integrals_erf_map,i8,keys,values,n_elements)
   do k1=1,n_elements
    call bielec_integrals_index_reverse(kk,ii,ll,jj,keys(k1))
    do k2=1,8
     if (kk(k2)==0) then
      cycle
     endif
     i = ii(k2) ! electron 1
     j = jj(k2) ! electron 1
     k = kk(k2) ! electron 2
     l = ll(k2) ! electron 2
     ! values(k1) = (ij|kl) <=> <ik|jl>
     integral = values (k1)
 !!! Gaunt part of the interaction !!!
     if ((i .le. large_ao_num .and. j .gt. large_ao_num .and. k .gt. large_ao_num .and. l .le. large_ao_num)  .or.  &
         (i .gt. large_ao_num .and. j .le. large_ao_num .and. k .le. large_ao_num .and. l .gt. large_ao_num)) then  
      !L_alpha S_alpha .or S_alpha L_alpha
      dirac_ao_bi_elec_G_Hartree_integral_tmp(dirac_inverse_list(i,1),dirac_inverse_list(j,1)) -= (dirac_SCF_density_matrix_ao(dirac_inverse_list(l,1),dirac_inverse_list(k,1))) * integral
      dirac_ao_bi_elec_G_Hartree_integral_tmp(dirac_inverse_list(i,1),dirac_inverse_list(j,1)) += (dirac_SCF_density_matrix_ao(dirac_inverse_list(l,2),dirac_inverse_list(k,2))) * integral
      !L_beta S_beta .or S_beta L_beta
      dirac_ao_bi_elec_G_Hartree_integral_tmp(dirac_inverse_list(i,2),dirac_inverse_list(j,2)) -= (dirac_SCF_density_matrix_ao(dirac_inverse_list(l,2),dirac_inverse_list(k,2))) * integral
      dirac_ao_bi_elec_G_Hartree_integral_tmp(dirac_inverse_list(i,2),dirac_inverse_list(j,2)) += (dirac_SCF_density_matrix_ao(dirac_inverse_list(l,1),dirac_inverse_list(k,1))) * integral
      !L_beta S_alpha .or S_beta L_alpha
      dirac_ao_bi_elec_G_Hartree_integral_tmp(dirac_inverse_list(i,2),dirac_inverse_list(j,1)) -= 2*(dirac_SCF_density_matrix_ao(dirac_inverse_list(l,2),dirac_inverse_list(k,1))) * integral
      !L_alpha S_beta .or S_alpha L_beta
      dirac_ao_bi_elec_G_Hartree_integral_tmp(dirac_inverse_list(i,1),dirac_inverse_list(j,2)) -= 2*(dirac_SCF_density_matrix_ao(dirac_inverse_list(l,1),dirac_inverse_list(k,2))) * integral
     elseif ((i .le. large_ao_num .and. j .gt. large_ao_num .and. k .le. large_ao_num .and. l .gt. large_ao_num)  .or.  &
             (i .gt. large_ao_num .and. j .le. large_ao_num .and. k .gt. large_ao_num .and. l .le. large_ao_num)) then
      !L_alpha S_alpha .or S_alpha L_alpha
      dirac_ao_bi_elec_G_Hartree_integral_tmp(dirac_inverse_list(i,1),dirac_inverse_list(j,1)) -= (dirac_SCF_density_matrix_ao(dirac_inverse_list(l,1),dirac_inverse_list(k,1))) * integral
      dirac_ao_bi_elec_G_Hartree_integral_tmp(dirac_inverse_list(i,1),dirac_inverse_list(j,1)) += (dirac_SCF_density_matrix_ao(dirac_inverse_list(l,2),dirac_inverse_list(k,2))) * integral
      !L_beta S_beta .or S_beta L_beta
      dirac_ao_bi_elec_G_Hartree_integral_tmp(dirac_inverse_list(i,2),dirac_inverse_list(j,2)) -= (dirac_SCF_density_matrix_ao(dirac_inverse_list(l,2),dirac_inverse_list(k,2))) * integral
      dirac_ao_bi_elec_G_Hartree_integral_tmp(dirac_inverse_list(i,2),dirac_inverse_list(j,2)) += (dirac_SCF_density_matrix_ao(dirac_inverse_list(l,1),dirac_inverse_list(k,1))) * integral
      !L_beta S_alpha .or S_beta L_alpha
      dirac_ao_bi_elec_G_Hartree_integral_tmp(dirac_inverse_list(i,2),dirac_inverse_list(j,1)) -= 2*(dirac_SCF_density_matrix_ao(dirac_inverse_list(l,2),dirac_inverse_list(k,1))) * integral
      !L_alpha S_beta .or S_alpha L_beta
      dirac_ao_bi_elec_G_Hartree_integral_tmp(dirac_inverse_list(i,1),dirac_inverse_list(j,2)) -= 2*(dirac_SCF_density_matrix_ao(dirac_inverse_list(l,1),dirac_inverse_list(k,2))) * integral
     endif
    enddo
   enddo
  enddo
!!$OMP END DO NOWAIT
!!$OMP CRITICAL
  dirac_ao_bi_elec_G_Hartree_integral += dirac_ao_bi_elec_G_Hartree_integral_tmp
!!$OMP END CRITICAL
  deallocate(keys,values,dirac_ao_bi_elec_G_Hartree_integral_tmp)
!!$OMP END PARALLEL
 END_PROVIDER

 BEGIN_PROVIDER [ complex*16, dirac_HF_two_electron_G_Hartree_energy_complex]
 &BEGIN_PROVIDER [ double precision, dirac_HF_two_electron_G_Hartree_energy] 
  implicit none
  BEGIN_DOC
  !Two-electrons energy of the Coulomb_Gaunt ee interaction
  !The energy is supposed to be a real, thus we check for its complex part to be
  ! a VERY small artifact and take only its real part
  END_DOC
  integer :: i,j
  dirac_HF_two_electron_G_Hartree_energy_complex = (0.d0,0.d0)
  do j=1, 2*dirac_ao_num
   do i=1, 2*dirac_ao_num
    dirac_HF_two_electron_G_Hartree_energy_complex += 0.5d0* (dirac_ao_bi_elec_G_Hartree_integral(i,j)) * dirac_SCF_density_matrix_ao(j,i)
   enddo
  enddo
  dirac_HF_two_electron_G_Hartree_energy = real(dirac_HF_two_electron_G_Hartree_energy_complex)
  if (aimag(dirac_HF_two_electron_G_Hartree_energy_complex) .gt. 1.d-10) then
  print*, 'Warning! The energy is not real'
  print*, 'dirac_HF_two_electron_G_Hartree_energy_complex =', dirac_HF_two_electron_G_Hartree_energy_complex
  STOP
  endif
 END_PROVIDER


 BEGIN_PROVIDER [ complex*16, dirac_ao_bi_elec_G_Exchange_integral, (2*dirac_ao_num, 2*dirac_ao_num) ]
  use map_module
  implicit none
  BEGIN_DOC
  !Array of the bi-electronic Fock matrix for the Coulomb-Gaunt interaction,
  ! in Dirac AO basis set.
  !Take care, the density matrix index
  ! have been correctly inverse, unlike
  ! in the non-relativistic code
  END_DOC
  PROVIDE dirac_ao_bielec_integrals_erf_in_map
  integer                        :: i,j,k,l,k1 
  integer                        :: r,s,p,q
  double precision               :: dirac_ao_bielec_integral, local_threshold
  double precision               :: integral
  complex*16, allocatable        :: dirac_ao_bi_elec_G_Exchange_integral_tmp(:,:)
  integer(omp_lock_kind)         :: lck(dirac_ao_num)
  integer*8                      :: i8
  integer                        :: ii(8), jj(8), kk(8), ll(8), k2
  integer(cache_map_size_kind)   :: n_elements_max, n_elements
  integer(key_kind), allocatable :: keys(:)
  double precision, allocatable  :: values(:)
  dirac_ao_bi_elec_G_Exchange_integral = (0.d0,0.d0)
!!$OMP PARALLEL DEFAULT(NONE)                                      &
!!$OMP PRIVATE(i,j,l,k1,k,integral,ii,jj,kk,ll,i8,keys,values,n_elements_max, &
!!$OMP  n_elements,dirac_ao_bi_elec_G_Exchange_integral_tmp)&
!!$OMP SHARED(dirac_ao_num,dirac_SCF_density_matrix_ao,&
!!$OMP  dirac_ao_integrals_map, dirac_ao_bi_elec_G_Exchange_integral) 
  call get_cache_map_n_elements_max(dirac_ao_integrals_erf_map,n_elements_max)
  allocate(keys(n_elements_max), values(n_elements_max))
  allocate(dirac_ao_bi_elec_G_Exchange_integral_tmp(2*dirac_ao_num,2*dirac_ao_num))
  dirac_ao_bi_elec_G_Exchange_integral_tmp = (0.d0,0.d0)
!!$OMP DO SCHEDULE(dynamic,64)
!!DIR$ NOVECTOR
  do i8=0_8,dirac_ao_integrals_erf_map%map_size
   n_elements = n_elements_max
   call get_cache_map(dirac_ao_integrals_erf_map,i8,keys,values,n_elements)
   do k1=1,n_elements
    call bielec_integrals_index_reverse(kk,ii,ll,jj,keys(k1))
    do k2=1,8
     if (kk(k2)==0) then
      cycle
     endif
     i = ii(k2) ! electron 1
     j = jj(k2) ! electron 1
     k = kk(k2) ! electron 2
     l = ll(k2) ! electron 2
     ! values(k1) = (ij|kl) <=> <ik|jl>
     integral = values (k1)
 !!! Gaunt part of the interaction !!!
     if ((i .le. large_ao_num .and. j .gt. large_ao_num .and. k .gt. large_ao_num .and. l .le. large_ao_num)  .or.  &
         (i .gt. large_ao_num .and. j .le. large_ao_num .and. k .le. large_ao_num .and. l .gt. large_ao_num)) then  
      !L_alpha L_alpha .or. S_alpha S_alpha
      dirac_ao_bi_elec_G_Exchange_integral_tmp(dirac_inverse_list(i,1),dirac_inverse_list(l,1)) += 2*(dirac_SCF_density_matrix_ao(dirac_inverse_list(j,2),dirac_inverse_list(k,2))) * integral 
      dirac_ao_bi_elec_G_Exchange_integral_tmp(dirac_inverse_list(i,1),dirac_inverse_list(l,1)) += (dirac_SCF_density_matrix_ao(dirac_inverse_list(j,1),dirac_inverse_list(k,1))) * integral
      !L_beta L_beta .or. S_beta S_beta
      dirac_ao_bi_elec_G_Exchange_integral_tmp(dirac_inverse_list(i,2),dirac_inverse_list(l,2)) += 2*(dirac_SCF_density_matrix_ao(dirac_inverse_list(j,1),dirac_inverse_list(k,1))) * integral
      dirac_ao_bi_elec_G_Exchange_integral_tmp(dirac_inverse_list(i,2),dirac_inverse_list(l,2)) += (dirac_SCF_density_matrix_ao(dirac_inverse_list(j,2),dirac_inverse_list(k,2))) * integral
      !L_beta L_alpha .or. S_beta S_alpha
      dirac_ao_bi_elec_G_Exchange_integral_tmp(dirac_inverse_list(i,2),dirac_inverse_list(l,1)) -= (dirac_SCF_density_matrix_ao(dirac_inverse_list(j,2),dirac_inverse_list(k,1))) * integral
      !L_alpha L_beta .or. S_alpha S_beta
      dirac_ao_bi_elec_G_Exchange_integral_tmp(dirac_inverse_list(i,1),dirac_inverse_list(l,2)) -= (dirac_SCF_density_matrix_ao(dirac_inverse_list(j,1),dirac_inverse_list(k,2))) * integral
     elseif ((i .le. large_ao_num .and. j .gt. large_ao_num .and. k .le. large_ao_num .and. l .gt. large_ao_num)  .or.  &
             (i .gt. large_ao_num .and. j .le. large_ao_num .and. k .gt. large_ao_num .and. l .le. large_ao_num)) then
      !L_alpha S_alpha .or S_alpha L_alpha
      dirac_ao_bi_elec_G_Exchange_integral_tmp(dirac_inverse_list(i,1),dirac_inverse_list(l,1)) += 2*(dirac_SCF_density_matrix_ao(dirac_inverse_list(j,2),dirac_inverse_list(k,2))) * integral
      dirac_ao_bi_elec_G_Exchange_integral_tmp(dirac_inverse_list(i,1),dirac_inverse_list(l,1)) += (dirac_SCF_density_matrix_ao(dirac_inverse_list(j,1),dirac_inverse_list(k,1))) * integral
      !L_beta S_beta .or S_beta L_beta
      dirac_ao_bi_elec_G_Exchange_integral_tmp(dirac_inverse_list(i,2),dirac_inverse_list(l,2)) += 2*(dirac_SCF_density_matrix_ao(dirac_inverse_list(j,1),dirac_inverse_list(k,1))) * integral
      dirac_ao_bi_elec_G_Exchange_integral_tmp(dirac_inverse_list(i,2),dirac_inverse_list(l,2)) += (dirac_SCF_density_matrix_ao(dirac_inverse_list(j,2),dirac_inverse_list(k,2))) * integral
      !L_beta S_alpha .or S_beta L_alpha
      dirac_ao_bi_elec_G_Exchange_integral_tmp(dirac_inverse_list(i,2),dirac_inverse_list(l,1)) -= (dirac_SCF_density_matrix_ao(dirac_inverse_list(j,2),dirac_inverse_list(k,1))) * integral
      !L_alpha S_beta .or S_alpha L_beta
      dirac_ao_bi_elec_G_Exchange_integral_tmp(dirac_inverse_list(i,1),dirac_inverse_list(l,2)) -= (dirac_SCF_density_matrix_ao(dirac_inverse_list(j,1),dirac_inverse_list(k,2))) * integral
     endif
    enddo
   enddo
  enddo
!!$OMP END DO NOWAIT
!!$OMP CRITICAL
  dirac_ao_bi_elec_G_Exchange_integral += dirac_ao_bi_elec_G_Exchange_integral_tmp
!!$OMP END CRITICAL
  deallocate(keys,values,dirac_ao_bi_elec_G_Exchange_integral_tmp)
!!$OMP END PARALLEL
 END_PROVIDER

 BEGIN_PROVIDER [ complex*16, dirac_HF_two_electron_G_Exchange_energy_complex]
 &BEGIN_PROVIDER [ double precision, dirac_HF_two_electron_G_Exchange_energy] 
  implicit none
  BEGIN_DOC
  !Two-electrons energy of the Coulomb_Gaunt ee interaction
  !The energy is supposed to be a real, thus we check for its complex part to be
  ! a VERY small artifact and take only its real part
  END_DOC
  integer :: i,j
  dirac_HF_two_electron_G_Exchange_energy_complex = (0.d0,0.d0)
  do j=1, 2*dirac_ao_num
   do i=1, 2*dirac_ao_num
    dirac_HF_two_electron_G_Exchange_energy_complex += 0.5d0* (dirac_ao_bi_elec_G_Exchange_integral(i,j)) * dirac_SCF_density_matrix_ao(j,i)
   enddo
  enddo
  dirac_HF_two_electron_G_Exchange_energy = real(dirac_HF_two_electron_G_Exchange_energy_complex)
  if (aimag(dirac_HF_two_electron_G_Exchange_energy_complex) .gt. 1.d-10) then
  print*, 'Warning! The energy is not real'
  print*, 'dirac_HF_two_electron_G_Exchange_energy_complex =', dirac_HF_two_electron_G_Exchange_energy_complex
  STOP
  endif
 END_PROVIDER
 
 BEGIN_PROVIDER [ complex*16, dirac_ao_bi_elec_C_G_Hartree_integral, (2*dirac_ao_num, 2*dirac_ao_num) ]
  use map_module
  implicit none
  BEGIN_DOC
  !Array of the bi-electronic Fock matrix for the Coulomb-Gaunt interaction,
  ! in Dirac AO basis set.
  !Take care, the density matrix index
  ! have been correctly inverse, unlike
  ! in the non-relativistic code
  END_DOC
  PROVIDE dirac_ao_bielec_integrals_erf_in_map
  integer                        :: i,j,k,l,k1 
  integer                        :: r,s,p,q
  double precision               :: dirac_ao_bielec_integral, local_threshold
  double precision               :: integral
  complex*16, allocatable        :: dirac_ao_bi_elec_C_G_Hartree_integral_tmp(:,:)
  integer(omp_lock_kind)         :: lck(dirac_ao_num)
  integer*8                      :: i8
  integer                        :: ii(8), jj(8), kk(8), ll(8), k2
  integer(cache_map_size_kind)   :: n_elements_max, n_elements
  integer(key_kind), allocatable :: keys(:)
  double precision, allocatable  :: values(:)
  dirac_ao_bi_elec_C_G_Hartree_integral = (0.d0,0.d0)
!!$OMP PARALLEL DEFAULT(NONE)                                      &
!!$OMP PRIVATE(i,j,l,k1,k,integral,ii,jj,kk,ll,i8,keys,values,n_elements_max, &
!!$OMP  n_elements,dirac_ao_bi_elec_C_G_Hartree_integral_tmp)&
!!$OMP SHARED(dirac_ao_num,dirac_SCF_density_matrix_ao,&
!!$OMP  dirac_ao_integrals_map, dirac_ao_bi_elec_C_G_Hartree_integral) 
  call get_cache_map_n_elements_max(dirac_ao_integrals_erf_map,n_elements_max)
  allocate(keys(n_elements_max), values(n_elements_max))
  allocate(dirac_ao_bi_elec_C_G_Hartree_integral_tmp(2*dirac_ao_num,2*dirac_ao_num))
  dirac_ao_bi_elec_C_G_Hartree_integral_tmp = (0.d0,0.d0)
!!$OMP DO SCHEDULE(dynamic,64)
!!DIR$ NOVECTOR
  do i8=0_8,dirac_ao_integrals_erf_map%map_size
   n_elements = n_elements_max
   call get_cache_map(dirac_ao_integrals_erf_map,i8,keys,values,n_elements)
   do k1=1,n_elements
    call bielec_integrals_index_reverse(kk,ii,ll,jj,keys(k1))
    do k2=1,8
     if (kk(k2)==0) then
      cycle
     endif
     i = ii(k2) ! electron 1
     j = jj(k2) ! electron 1
     k = kk(k2) ! electron 2
     l = ll(k2) ! electron 2
     ! values(k1) = (ij|kl) <=> <ik|jl>
     integral = values (k1)
 !!! Coulomb part of the interaction !!!
     if ((i .le. large_ao_num .and. j .le. large_ao_num .and. k .le. large_ao_num .and. l .le. large_ao_num) .or.  &
         (i .gt. large_ao_num .and. j .gt. large_ao_num .and. k .gt. large_ao_num .and. l .gt. large_ao_num)) then
      !L_alpha L_alpha .or. S_alpha S_alpha
      dirac_ao_bi_elec_C_G_Hartree_integral_tmp(dirac_inverse_list(i,1),dirac_inverse_list(j,1)) += (dirac_SCF_density_matrix_ao(dirac_inverse_list(l,1),dirac_inverse_list(k,1))) * integral 
      dirac_ao_bi_elec_C_G_Hartree_integral_tmp(dirac_inverse_list(i,1),dirac_inverse_list(j,1)) += (dirac_SCF_density_matrix_ao(dirac_inverse_list(l,2),dirac_inverse_list(k,2))) * integral
      !L_beta L_beta .or. S_beta S_beta
      dirac_ao_bi_elec_C_G_Hartree_integral_tmp(dirac_inverse_list(i,2),dirac_inverse_list(j,2)) += (dirac_SCF_density_matrix_ao(dirac_inverse_list(l,1),dirac_inverse_list(k,1))) * integral
      dirac_ao_bi_elec_C_G_Hartree_integral_tmp(dirac_inverse_list(i,2),dirac_inverse_list(j,2)) += (dirac_SCF_density_matrix_ao(dirac_inverse_list(l,2),dirac_inverse_list(k,2))) * integral
     elseif((i .le. large_ao_num .and. j .le. large_ao_num .and. k .gt. large_ao_num .and. l .gt. large_ao_num) .or. &
            (i .gt. large_ao_num .and. j .gt. large_ao_num .and. k .le. large_ao_num .and. l .le. large_ao_num))then
       !L_alpha L_alpha .or. S_alpha S_alpha
      dirac_ao_bi_elec_C_G_Hartree_integral_tmp(dirac_inverse_list(i,1),dirac_inverse_list(j,1)) += (dirac_SCF_density_matrix_ao(dirac_inverse_list(l,1),dirac_inverse_list(k,1))) * integral
      dirac_ao_bi_elec_C_G_Hartree_integral_tmp(dirac_inverse_list(i,1),dirac_inverse_list(j,1)) += (dirac_SCF_density_matrix_ao(dirac_inverse_list(l,2),dirac_inverse_list(k,2))) * integral
      !L_beta L_beta .or. S_beta S_beta
      dirac_ao_bi_elec_C_G_Hartree_integral_tmp(dirac_inverse_list(i,2),dirac_inverse_list(j,2)) += (dirac_SCF_density_matrix_ao(dirac_inverse_list(l,1),dirac_inverse_list(k,1))) * integral
      dirac_ao_bi_elec_C_G_Hartree_integral_tmp(dirac_inverse_list(i,2),dirac_inverse_list(j,2)) += (dirac_SCF_density_matrix_ao(dirac_inverse_list(l,2),dirac_inverse_list(k,2))) * integral
     endif
 !!! Gaunt part of the interaction !!!
     if ((i .le. large_ao_num .and. j .gt. large_ao_num .and. k .gt. large_ao_num .and. l .le. large_ao_num)  .or.  &
         (i .gt. large_ao_num .and. j .le. large_ao_num .and. k .le. large_ao_num .and. l .gt. large_ao_num)) then  
      !L_alpha S_alpha .or S_alpha L_alpha
      dirac_ao_bi_elec_C_G_Hartree_integral_tmp(dirac_inverse_list(i,1),dirac_inverse_list(j,1)) -= (dirac_SCF_density_matrix_ao(dirac_inverse_list(l,1),dirac_inverse_list(k,1))) * integral
      dirac_ao_bi_elec_C_G_Hartree_integral_tmp(dirac_inverse_list(i,1),dirac_inverse_list(j,1)) += (dirac_SCF_density_matrix_ao(dirac_inverse_list(l,2),dirac_inverse_list(k,2))) * integral
      !L_beta S_beta .or S_beta L_beta
      dirac_ao_bi_elec_C_G_Hartree_integral_tmp(dirac_inverse_list(i,2),dirac_inverse_list(j,2)) -= (dirac_SCF_density_matrix_ao(dirac_inverse_list(l,2),dirac_inverse_list(k,2))) * integral
      dirac_ao_bi_elec_C_G_Hartree_integral_tmp(dirac_inverse_list(i,2),dirac_inverse_list(j,2)) += (dirac_SCF_density_matrix_ao(dirac_inverse_list(l,1),dirac_inverse_list(k,1))) * integral
      !L_beta S_alpha .or S_beta L_alpha
      dirac_ao_bi_elec_C_G_Hartree_integral_tmp(dirac_inverse_list(i,2),dirac_inverse_list(j,1)) -= 2*(dirac_SCF_density_matrix_ao(dirac_inverse_list(l,2),dirac_inverse_list(k,1))) * integral
      !L_alpha S_beta .or S_alpha L_beta
      dirac_ao_bi_elec_C_G_Hartree_integral_tmp(dirac_inverse_list(i,1),dirac_inverse_list(j,2)) -= 2*(dirac_SCF_density_matrix_ao(dirac_inverse_list(l,1),dirac_inverse_list(k,2))) * integral
     elseif ((i .le. large_ao_num .and. j .gt. large_ao_num .and. k .le. large_ao_num .and. l .gt. large_ao_num)  .or.  &
             (i .gt. large_ao_num .and. j .le. large_ao_num .and. k .gt. large_ao_num .and. l .le. large_ao_num)) then
      !L_alpha S_alpha .or S_alpha L_alpha
      dirac_ao_bi_elec_C_G_Hartree_integral_tmp(dirac_inverse_list(i,1),dirac_inverse_list(j,1)) -= (dirac_SCF_density_matrix_ao(dirac_inverse_list(l,1),dirac_inverse_list(k,1))) * integral
      dirac_ao_bi_elec_C_G_Hartree_integral_tmp(dirac_inverse_list(i,1),dirac_inverse_list(j,1)) += (dirac_SCF_density_matrix_ao(dirac_inverse_list(l,2),dirac_inverse_list(k,2))) * integral
      !L_beta S_beta .or S_beta L_beta
      dirac_ao_bi_elec_C_G_Hartree_integral_tmp(dirac_inverse_list(i,2),dirac_inverse_list(j,2)) -= (dirac_SCF_density_matrix_ao(dirac_inverse_list(l,2),dirac_inverse_list(k,2))) * integral
      dirac_ao_bi_elec_C_G_Hartree_integral_tmp(dirac_inverse_list(i,2),dirac_inverse_list(j,2)) += (dirac_SCF_density_matrix_ao(dirac_inverse_list(l,1),dirac_inverse_list(k,1))) * integral
      !L_beta S_alpha .or S_beta L_alpha
      dirac_ao_bi_elec_C_G_Hartree_integral_tmp(dirac_inverse_list(i,2),dirac_inverse_list(j,1)) -= 2*(dirac_SCF_density_matrix_ao(dirac_inverse_list(l,2),dirac_inverse_list(k,1))) * integral
      !L_alpha S_beta .or S_alpha L_beta
      dirac_ao_bi_elec_C_G_Hartree_integral_tmp(dirac_inverse_list(i,1),dirac_inverse_list(j,2)) -= 2*(dirac_SCF_density_matrix_ao(dirac_inverse_list(l,1),dirac_inverse_list(k,2))) * integral
     endif
    enddo
   enddo
  enddo
!!$OMP END DO NOWAIT
!!$OMP CRITICAL
  dirac_ao_bi_elec_C_G_Hartree_integral += dirac_ao_bi_elec_C_G_Hartree_integral_tmp
!!$OMP END CRITICAL
  deallocate(keys,values,dirac_ao_bi_elec_C_G_Hartree_integral_tmp)
!!$OMP END PARALLEL
 END_PROVIDER

 BEGIN_PROVIDER [ complex*16, dirac_HF_two_electron_C_G_Hartree_energy_complex]
 &BEGIN_PROVIDER [ double precision, dirac_HF_two_electron_C_G_Hartree_energy] 
  implicit none
  BEGIN_DOC
  !Two-electrons energy of the Coulomb_Gaunt ee interaction
  !The energy is supposed to be a real, thus we check for its complex part to be
  ! a VERY small artifact and take only its real part
  END_DOC
  integer :: i,j
  dirac_HF_two_electron_C_G_Hartree_energy_complex = (0.d0,0.d0)
  do j=1, 2*dirac_ao_num
   do i=1, 2*dirac_ao_num
    dirac_HF_two_electron_C_G_Hartree_energy_complex += 0.5d0* (dirac_ao_bi_elec_C_G_Hartree_integral(i,j)) * dirac_SCF_density_matrix_ao(j,i)
   enddo
  enddo
  dirac_HF_two_electron_C_G_Hartree_energy = real(dirac_HF_two_electron_C_G_Hartree_energy_complex)
  if (aimag(dirac_HF_two_electron_C_G_Hartree_energy_complex) .gt. 1.d-10) then
  print*, 'Warning! The energy is not real'
  print*, 'dirac_HF_two_electron_C_G_Hartree_energy_complex =', dirac_HF_two_electron_C_G_Hartree_energy_complex
  STOP
  endif
 END_PROVIDER
 
 BEGIN_PROVIDER [ complex*16, dirac_ao_bi_elec_C_G_Exchange_integral, (2*dirac_ao_num, 2*dirac_ao_num) ]
  use map_module
  implicit none
  BEGIN_DOC
  !Array of the bi-electronic Fock matrix for the Coulomb-Gaunt interaction,
  ! in Dirac AO basis set.
  !Take care, the density matrix index
  ! have been correctly inverse, unlike
  ! in the non-relativistic code
  END_DOC
  PROVIDE dirac_ao_bielec_integrals_erf_in_map
  integer                        :: i,j,k,l,k1 
  integer                        :: r,s,p,q
  double precision               :: dirac_ao_bielec_integral, local_threshold
  double precision               :: integral
  complex*16, allocatable        :: dirac_ao_bi_elec_C_G_Exchange_integral_tmp(:,:)
  integer(omp_lock_kind)         :: lck(dirac_ao_num)
  integer*8                      :: i8
  integer                        :: ii(8), jj(8), kk(8), ll(8), k2
  integer(cache_map_size_kind)   :: n_elements_max, n_elements
  integer(key_kind), allocatable :: keys(:)
  double precision, allocatable  :: values(:)
  dirac_ao_bi_elec_C_G_Exchange_integral = (0.d0,0.d0)
!!$OMP PARALLEL DEFAULT(NONE)                                      &
!!$OMP PRIVATE(i,j,l,k1,k,integral,ii,jj,kk,ll,i8,keys,values,n_elements_max, &
!!$OMP  n_elements,dirac_ao_bi_elec_C_G_Exchange_integral_tmp)&
!!$OMP SHARED(dirac_ao_num,dirac_SCF_density_matrix_ao,&
!!$OMP  dirac_ao_integrals_map, dirac_ao_bi_elec_C_G_Exchange_integral) 
  call get_cache_map_n_elements_max(dirac_ao_integrals_erf_map,n_elements_max)
  allocate(keys(n_elements_max), values(n_elements_max))
  allocate(dirac_ao_bi_elec_C_G_Exchange_integral_tmp(2*dirac_ao_num,2*dirac_ao_num))
  dirac_ao_bi_elec_C_G_Exchange_integral_tmp = (0.d0,0.d0)
!!$OMP DO SCHEDULE(dynamic,64)
!!DIR$ NOVECTOR
  do i8=0_8,dirac_ao_integrals_erf_map%map_size
   n_elements = n_elements_max
   call get_cache_map(dirac_ao_integrals_erf_map,i8,keys,values,n_elements)
   do k1=1,n_elements
    call bielec_integrals_index_reverse(kk,ii,ll,jj,keys(k1))
    do k2=1,8
     if (kk(k2)==0) then
      cycle
     endif
     i = ii(k2) ! electron 1
     j = jj(k2) ! electron 1
     k = kk(k2) ! electron 2
     l = ll(k2) ! electron 2
     ! values(k1) = (ij|kl) <=> <ik|jl>
     integral = values (k1)
 !!! Coulomb part of the interaction !!!
     if ((i .le. large_ao_num .and. j .le. large_ao_num .and. k .le. large_ao_num .and. l .le. large_ao_num) .or.  &
         (i .gt. large_ao_num .and. j .gt. large_ao_num .and. k .gt. large_ao_num .and. l .gt. large_ao_num)) then
      !L_alpha L_alpha .or. S_alpha S_alpha
      dirac_ao_bi_elec_C_G_Exchange_integral_tmp(dirac_inverse_list(i,1),dirac_inverse_list(l,1)) -= (dirac_SCF_density_matrix_ao(dirac_inverse_list(j,1),dirac_inverse_list(k,1))) * integral
      !L_beta L_beta .or. S_beta S_beta
      dirac_ao_bi_elec_C_G_Exchange_integral_tmp(dirac_inverse_list(i,2),dirac_inverse_list(l,2)) -= (dirac_SCF_density_matrix_ao(dirac_inverse_list(j,2),dirac_inverse_list(k,2))) * integral   
      !L_beta L_alpha .or. S_beta S_alpha
      dirac_ao_bi_elec_C_G_Exchange_integral_tmp(dirac_inverse_list(i,2),dirac_inverse_list(l,1)) -= (dirac_SCF_density_matrix_ao(dirac_inverse_list(j,2),dirac_inverse_list(k,1))) * integral
      !L_alpha L_beta .or. S_alpha S_beta
      dirac_ao_bi_elec_C_G_Exchange_integral_tmp(dirac_inverse_list(i,1),dirac_inverse_list(l,2)) -= (dirac_SCF_density_matrix_ao(dirac_inverse_list(j,1),dirac_inverse_list(k,2))) * integral
     elseif((i .le. large_ao_num .and. j .le. large_ao_num .and. k .gt. large_ao_num .and. l .gt. large_ao_num) .or. &
            (i .gt. large_ao_num .and. j .gt. large_ao_num .and. k .le. large_ao_num .and. l .le. large_ao_num))then
      !L_alpha S_alpha .or S_alpha L_alpha
      dirac_ao_bi_elec_C_G_Exchange_integral_tmp(dirac_inverse_list(i,1),dirac_inverse_list(l,1)) -= (dirac_SCF_density_matrix_ao(dirac_inverse_list(j,1),dirac_inverse_list(k,1))) * integral
      !L_beta S_beta .or. S_beta L_beta
      dirac_ao_bi_elec_C_G_Exchange_integral_tmp(dirac_inverse_list(i,2),dirac_inverse_list(l,2)) -= (dirac_SCF_density_matrix_ao(dirac_inverse_list(j,2),dirac_inverse_list(k,2))) * integral
      !L_beta S_alpha .or. S_beta L_alpha
      dirac_ao_bi_elec_C_G_Exchange_integral_tmp(dirac_inverse_list(i,2),dirac_inverse_list(l,1)) -= (dirac_SCF_density_matrix_ao(dirac_inverse_list(j,2),dirac_inverse_list(k,1))) * integral
      !L_alpha S_beta .or. S_alpha L_beta
      dirac_ao_bi_elec_C_G_Exchange_integral_tmp(dirac_inverse_list(i,1),dirac_inverse_list(l,2)) -= (dirac_SCF_density_matrix_ao(dirac_inverse_list(j,1),dirac_inverse_list(k,2))) * integral
     endif
 !!! Gaunt part of the interaction !!!
     if ((i .le. large_ao_num .and. j .gt. large_ao_num .and. k .gt. large_ao_num .and. l .le. large_ao_num)  .or.  &
         (i .gt. large_ao_num .and. j .le. large_ao_num .and. k .le. large_ao_num .and. l .gt. large_ao_num)) then  
      !L_alpha L_alpha .or. S_alpha S_alpha
      dirac_ao_bi_elec_C_G_Exchange_integral_tmp(dirac_inverse_list(i,1),dirac_inverse_list(l,1)) += 2*(dirac_SCF_density_matrix_ao(dirac_inverse_list(j,2),dirac_inverse_list(k,2))) * integral 
      dirac_ao_bi_elec_C_G_Exchange_integral_tmp(dirac_inverse_list(i,1),dirac_inverse_list(l,1)) += (dirac_SCF_density_matrix_ao(dirac_inverse_list(j,1),dirac_inverse_list(k,1))) * integral
      !L_beta L_beta .or. S_beta S_beta
      dirac_ao_bi_elec_C_G_Exchange_integral_tmp(dirac_inverse_list(i,2),dirac_inverse_list(l,2)) += 2*(dirac_SCF_density_matrix_ao(dirac_inverse_list(j,1),dirac_inverse_list(k,1))) * integral
      dirac_ao_bi_elec_C_G_Exchange_integral_tmp(dirac_inverse_list(i,2),dirac_inverse_list(l,2)) += (dirac_SCF_density_matrix_ao(dirac_inverse_list(j,2),dirac_inverse_list(k,2))) * integral
      !L_beta L_alpha .or. S_beta S_alpha
      dirac_ao_bi_elec_C_G_Exchange_integral_tmp(dirac_inverse_list(i,2),dirac_inverse_list(l,1)) -= (dirac_SCF_density_matrix_ao(dirac_inverse_list(j,2),dirac_inverse_list(k,1))) * integral
      !L_alpha L_beta .or. S_alpha S_beta
      dirac_ao_bi_elec_C_G_Exchange_integral_tmp(dirac_inverse_list(i,1),dirac_inverse_list(l,2)) -= (dirac_SCF_density_matrix_ao(dirac_inverse_list(j,1),dirac_inverse_list(k,2))) * integral
     elseif ((i .le. large_ao_num .and. j .gt. large_ao_num .and. k .le. large_ao_num .and. l .gt. large_ao_num)  .or.  &
             (i .gt. large_ao_num .and. j .le. large_ao_num .and. k .gt. large_ao_num .and. l .le. large_ao_num)) then
      !L_alpha S_alpha .or S_alpha L_alpha
      dirac_ao_bi_elec_C_G_Exchange_integral_tmp(dirac_inverse_list(i,1),dirac_inverse_list(l,1)) += 2*(dirac_SCF_density_matrix_ao(dirac_inverse_list(j,2),dirac_inverse_list(k,2))) * integral
      dirac_ao_bi_elec_C_G_Exchange_integral_tmp(dirac_inverse_list(i,1),dirac_inverse_list(l,1)) += (dirac_SCF_density_matrix_ao(dirac_inverse_list(j,1),dirac_inverse_list(k,1))) * integral
      !L_beta S_beta .or S_beta L_beta
      dirac_ao_bi_elec_C_G_Exchange_integral_tmp(dirac_inverse_list(i,2),dirac_inverse_list(l,2)) += 2*(dirac_SCF_density_matrix_ao(dirac_inverse_list(j,1),dirac_inverse_list(k,1))) * integral
      dirac_ao_bi_elec_C_G_Exchange_integral_tmp(dirac_inverse_list(i,2),dirac_inverse_list(l,2)) += (dirac_SCF_density_matrix_ao(dirac_inverse_list(j,2),dirac_inverse_list(k,2))) * integral
      !L_beta S_alpha .or S_beta L_alpha
      dirac_ao_bi_elec_C_G_Exchange_integral_tmp(dirac_inverse_list(i,2),dirac_inverse_list(l,1)) -= (dirac_SCF_density_matrix_ao(dirac_inverse_list(j,2),dirac_inverse_list(k,1))) * integral
      !L_alpha S_beta .or S_alpha L_beta
      dirac_ao_bi_elec_C_G_Exchange_integral_tmp(dirac_inverse_list(i,1),dirac_inverse_list(l,2)) -= (dirac_SCF_density_matrix_ao(dirac_inverse_list(j,1),dirac_inverse_list(k,2))) * integral
     endif
    enddo
   enddo
  enddo
!!$OMP END DO NOWAIT
!!$OMP CRITICAL
  dirac_ao_bi_elec_C_G_Exchange_integral += dirac_ao_bi_elec_C_G_Exchange_integral_tmp
!!$OMP END CRITICAL
  deallocate(keys,values,dirac_ao_bi_elec_C_G_Exchange_integral_tmp)
!!$OMP END PARALLEL
 END_PROVIDER

 BEGIN_PROVIDER [ complex*16, dirac_HF_two_electron_C_G_Exchange_energy_complex]
 &BEGIN_PROVIDER [ double precision, dirac_HF_two_electron_C_G_Exchange_energy] 
  implicit none
  BEGIN_DOC
  !Two-electrons energy of the Coulomb_Gaunt ee interaction
  !The energy is supposed to be a real, thus we check for its complex part to be
  ! a VERY small artifact and take only its real part
  END_DOC
  integer :: i,j
  dirac_HF_two_electron_C_G_Exchange_energy_complex = (0.d0,0.d0)
  do j=1, 2*dirac_ao_num
   do i=1, 2*dirac_ao_num
    dirac_HF_two_electron_C_G_Exchange_energy_complex += 0.5d0* (dirac_ao_bi_elec_C_G_Exchange_integral(i,j)) * dirac_SCF_density_matrix_ao(j,i)
   enddo
  enddo
  dirac_HF_two_electron_C_G_Exchange_energy = real(dirac_HF_two_electron_C_G_Exchange_energy_complex)
  if (aimag(dirac_HF_two_electron_C_G_Exchange_energy_complex) .gt. 1.d-10) then
  print*, 'Warning! The energy is not real'
  print*, 'dirac_HF_two_electron_C_G_Exchange_energy_complex =', dirac_HF_two_electron_C_G_Exchange_energy_complex
  STOP
  endif
 END_PROVIDER
