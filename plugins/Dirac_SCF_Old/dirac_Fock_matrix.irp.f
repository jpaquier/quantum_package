!BEGIN_PROVIDER [ complex*16, dirac_ao_mono_elec_integral,(2*dirac_ao_num,2*dirac_ao_num)]
!&BEGIN_PROVIDER [ complex*16, dirac_ao_mono_elec_integral_diag, (2*dirac_ao_num) ] 
!implicit none
! integer          :: i,j
! BEGIN_DOC
! !Array of the mono electronic hamiltonian on the dirac AO basis
! ! in the 4x4 component formalism with cartesian basis and 
! ! the unrestricted kinetic-balance scheme  
! END_DOC
! print*,'Computing the mono-electronic Fock matrix'
! dirac_ao_mono_elec_integral = (0.d0,0.d0)
! do j = 1, 2*(dirac_ao_num)
!  do i = 1, 2*(dirac_ao_num)
!   dirac_ao_mono_elec_integral(i,j) += (dirac_ao_mono_elec_nucl_integral(i,j) + dirac_ao_mono_elec_mass_integral(i,j) + dirac_ao_mono_elec_kinetic_integral(i,j) )
!  enddo
! enddo
! do j = 1, 2*dirac_ao_num
!  dirac_ao_mono_elec_integral_diag(j) = dirac_ao_mono_elec_integral(j,j)
! enddo
!END_PROVIDER



!BEGIN_PROVIDER [ complex*16, dirac_ao_bi_elec_integral, (2*dirac_ao_num, 2*dirac_ao_num) ]
! use map_module
! implicit none
! BEGIN_DOC
! !Array of the bi-electronic Fock matrix for the Coulomb interaction
! ! in Dirac AO basis set
! !Take care, the density matrix index
! ! have been correctly inverse, unlike
! ! in the non-relativistic code
! END_DOC
! PROVIDE dirac_ao_bielec_integrals_in_map
! integer                        :: i,j,k,l,k1 
! integer                        :: r,s,p,q
! double precision               :: dirac_ao_bielec_integral, local_threshold
! double precision               :: integral
! complex*16, allocatable        :: dirac_ao_bi_elec_integral_tmp(:,:)
! integer(omp_lock_kind)         :: lck(dirac_ao_num)
! integer*8                      :: i8
! integer                        :: ii(8), jj(8), kk(8), ll(8), k2
! integer(cache_map_size_kind)   :: n_elements_max, n_elements
! integer(key_kind), allocatable :: keys(:)
! double precision, allocatable  :: values(:)
! dirac_ao_bi_elec_integral = (0.d0,0.d0)
!!$OMP PARALLEL DEFAULT(NONE)                                      &
!!$OMP PRIVATE(i,j,l,k1,k,integral,ii,jj,kk,ll,i8,keys,values,n_elements_max, &
!!$OMP  n_elements,dirac_ao_bi_elec_integral_tmp)&
!!$OMP SHARED(dirac_ao_num,dirac_SCF_density_matrix_ao,&
!!$OMP  dirac_ao_integrals_map, dirac_ao_bi_elec_integral) 
! call get_cache_map_n_elements_max(dirac_ao_integrals_map,n_elements_max)
! allocate(keys(n_elements_max), values(n_elements_max))
! allocate(dirac_ao_bi_elec_integral_tmp(2*dirac_ao_num,2*dirac_ao_num))
! dirac_ao_bi_elec_integral_tmp = (0.d0,0.d0)
!!$OMP DO SCHEDULE(dynamic,64)
!!DIR$ NOVECTOR
! do i8=0_8,dirac_ao_integrals_map%map_size
!  n_elements = n_elements_max
!  call get_cache_map(dirac_ao_integrals_map,i8,keys,values,n_elements)
!  do k1=1,n_elements
!   call bielec_integrals_index_reverse(kk,ii,ll,jj,keys(k1))
!   do k2=1,8
!    if (kk(k2)==0) then
!     cycle
!    endif
!    i = ii(k2) ! electron 1
!    j = jj(k2) ! electron 1
!    k = kk(k2) ! electron 2
!    l = ll(k2) ! electron 2
!    ! values(k1) = (ij|kl) <=> <ik|jl>
!    integral = values (k1)
!    if ((i .le. large_ao_num .and. j .le. large_ao_num .and. k .le. large_ao_num .and. l .le. large_ao_num) .or.  &
!        (i .gt. large_ao_num .and. j .gt. large_ao_num .and. k .gt. large_ao_num .and. l .gt. large_ao_num)) then
!     !L_alpha L_alpha .or. S_alpha S_alpha
!     dirac_ao_bi_elec_integral_tmp(dirac_inverse_list(i,1),dirac_inverse_list(j,1)) += (dirac_SCF_density_matrix_ao(dirac_inverse_list(l,1),dirac_inverse_list(k,1))) * integral 
!     dirac_ao_bi_elec_integral_tmp(dirac_inverse_list(i,1),dirac_inverse_list(j,1)) += (dirac_SCF_density_matrix_ao(dirac_inverse_list(l,2),dirac_inverse_list(k,2))) * integral
!     dirac_ao_bi_elec_integral_tmp(dirac_inverse_list(i,1),dirac_inverse_list(l,1)) -= (dirac_SCF_density_matrix_ao(dirac_inverse_list(j,1),dirac_inverse_list(k,1))) * integral
!     !L_beta L_beta .or. S_beta S_beta
!     dirac_ao_bi_elec_integral_tmp(dirac_inverse_list(i,2),dirac_inverse_list(j,2)) += (dirac_SCF_density_matrix_ao(dirac_inverse_list(l,1),dirac_inverse_list(k,1))) * integral
!     dirac_ao_bi_elec_integral_tmp(dirac_inverse_list(i,2),dirac_inverse_list(j,2)) += (dirac_SCF_density_matrix_ao(dirac_inverse_list(l,2),dirac_inverse_list(k,2))) * integral
!     dirac_ao_bi_elec_integral_tmp(dirac_inverse_list(i,2),dirac_inverse_list(l,2)) -= (dirac_SCF_density_matrix_ao(dirac_inverse_list(j,2),dirac_inverse_list(k,2))) * integral   
!     !L_beta L_alpha .or. S_beta S_alpha
!     dirac_ao_bi_elec_integral_tmp(dirac_inverse_list(i,2),dirac_inverse_list(l,1)) -= (dirac_SCF_density_matrix_ao(dirac_inverse_list(j,2),dirac_inverse_list(k,1))) * integral
!     !L_alpha L_beta .or. S_alpha S_beta
!     dirac_ao_bi_elec_integral_tmp(dirac_inverse_list(i,1),dirac_inverse_list(l,2)) -= (dirac_SCF_density_matrix_ao(dirac_inverse_list(j,1),dirac_inverse_list(k,2))) * integral
!    elseif((i .le. large_ao_num .and. j .le. large_ao_num .and. k .gt. large_ao_num .and. l .gt. large_ao_num) .or. &
!           (i .gt. large_ao_num .and. j .gt. large_ao_num .and. k .le. large_ao_num .and. l .le. large_ao_num))then
!      !L_alpha L_alpha .or. S_alpha S_alpha
!     dirac_ao_bi_elec_integral_tmp(dirac_inverse_list(i,1),dirac_inverse_list(j,1)) += (dirac_SCF_density_matrix_ao(dirac_inverse_list(l,1),dirac_inverse_list(k,1))) * integral
!     dirac_ao_bi_elec_integral_tmp(dirac_inverse_list(i,1),dirac_inverse_list(j,1)) += (dirac_SCF_density_matrix_ao(dirac_inverse_list(l,2),dirac_inverse_list(k,2))) * integral
!     !L_beta L_beta .or. S_beta S_beta
!     dirac_ao_bi_elec_integral_tmp(dirac_inverse_list(i,2),dirac_inverse_list(j,2)) += (dirac_SCF_density_matrix_ao(dirac_inverse_list(l,1),dirac_inverse_list(k,1))) * integral
!     dirac_ao_bi_elec_integral_tmp(dirac_inverse_list(i,2),dirac_inverse_list(j,2)) += (dirac_SCF_density_matrix_ao(dirac_inverse_list(l,2),dirac_inverse_list(k,2))) * integral
!     !L_alpha S_alpha .or S_alpha L_alpha
!     dirac_ao_bi_elec_integral_tmp(dirac_inverse_list(i,1),dirac_inverse_list(l,1)) -= (dirac_SCF_density_matrix_ao(dirac_inverse_list(j,1),dirac_inverse_list(k,1))) * integral
!     !L_beta S_beta .or. S_beta L_beta
!     dirac_ao_bi_elec_integral_tmp(dirac_inverse_list(i,2),dirac_inverse_list(l,2)) -= (dirac_SCF_density_matrix_ao(dirac_inverse_list(j,2),dirac_inverse_list(k,2))) * integral
!     !L_beta S_alpha .or. S_beta L_alpha
!     dirac_ao_bi_elec_integral_tmp(dirac_inverse_list(i,2),dirac_inverse_list(l,1)) -= (dirac_SCF_density_matrix_ao(dirac_inverse_list(j,2),dirac_inverse_list(k,1))) * integral
!     !L_alpha S_beta .or. S_alpha L_beta
!     dirac_ao_bi_elec_integral_tmp(dirac_inverse_list(i,1),dirac_inverse_list(l,2)) -= (dirac_SCF_density_matrix_ao(dirac_inverse_list(j,1),dirac_inverse_list(k,2))) * integral
!    endif
!   enddo
!  enddo
! enddo
!!$OMP END DO NOWAIT
!!$OMP CRITICAL
! dirac_ao_bi_elec_integral += dirac_ao_bi_elec_integral_tmp
!!$OMP END CRITICAL
! deallocate(keys,values,dirac_ao_bi_elec_integral_tmp)
!!$OMP END PARALLEL
!END_PROVIDER

!BEGIN_PROVIDER [ complex*16, dirac_ao_bi_elec_erf_integral, (2*dirac_ao_num, 2*dirac_ao_num) ]
! use map_module
! implicit none
! BEGIN_DOC
! !Array of the bi-electronic Fock matrix for the long-range Coulomb interaction
! ! in Dirac AO basis set
! !Take care, the density matrix index
! ! have been correctly inverse, unlike
! ! in the non-relativistic code
! END_DOC
! PROVIDE dirac_ao_bielec_integrals_erf_in_map
! integer                        :: i,j,k,l,k1 
! integer                        :: r,s,p,q
! double precision               :: dirac_ao_bielec_integral_erf, local_threshold
! double precision               :: integral
! complex*16, allocatable        :: dirac_ao_bi_elec_erf_integral_tmp(:,:)
! integer(omp_lock_kind)         :: lck(dirac_ao_num)
! integer*8                      :: i8
! integer                        :: ii(8), jj(8), kk(8), ll(8), k2
! integer(cache_map_size_kind)   :: n_elements_max, n_elements
! integer(key_kind), allocatable :: keys(:)
! double precision, allocatable  :: values(:)
! dirac_ao_bi_elec_erf_integral = (0.d0,0.d0)
!!$OMP PARALLEL DEFAULT(NONE)                                      &
!!$OMP PRIVATE(i,j,l,k1,k,integral,ii,jj,kk,ll,i8,keys,values,n_elements_max, &
!!$OMP  n_elements,dirac_ao_bi_elec_erf_integral_tmp)&
!!$OMP SHARED(dirac_ao_num,dirac_SCF_density_matrix_ao,&
!!$OMP  dirac_ao_integrals_erf_map, dirac_ao_bi_elec_erf_integral) 
! call get_cache_map_n_elements_max(dirac_ao_integrals_erf_map,n_elements_max)
! allocate(keys(n_elements_max), values(n_elements_max))
! allocate(dirac_ao_bi_elec_erf_integral_tmp(2*dirac_ao_num,2*dirac_ao_num))
! dirac_ao_bi_elec_erf_integral_tmp = (0.d0,0.d0)
!!$OMP DO SCHEDULE(dynamic,64)
!!DIR$ NOVECTOR
! do i8=0_8,dirac_ao_integrals_erf_map%map_size
!  n_elements = n_elements_max
!  call get_cache_map(dirac_ao_integrals_erf_map,i8,keys,values,n_elements)
!  do k1=1,n_elements
!   call bielec_integrals_index_reverse(kk,ii,ll,jj,keys(k1))
!   do k2=1,8
!    if (kk(k2)==0) then
!     cycle
!    endif
!    i = ii(k2) ! electron 1
!    j = jj(k2) ! electron 1
!    k = kk(k2) ! electron 2
!    l = ll(k2) ! electron 2
!    ! values(k1) = (ij|kl) <=> <ik|jl>
!    integral = values (k1)
!    if ((i .le. large_ao_num .and. j .le. large_ao_num .and. k .le. large_ao_num .and. l .le. large_ao_num) .or.  &
!        (i .gt. large_ao_num .and. j .gt. large_ao_num .and. k .gt. large_ao_num .and. l .gt. large_ao_num)) then
!     !L_alpha L_alpha .or. S_alpha S_alpha
!     dirac_ao_bi_elec_erf_integral_tmp(dirac_inverse_list(i,1),dirac_inverse_list(j,1)) += (dirac_SCF_density_matrix_ao(dirac_inverse_list(l,1),dirac_inverse_list(k,1))) * integral 
!     dirac_ao_bi_elec_erf_integral_tmp(dirac_inverse_list(i,1),dirac_inverse_list(j,1)) += (dirac_SCF_density_matrix_ao(dirac_inverse_list(l,2),dirac_inverse_list(k,2))) * integral
!     dirac_ao_bi_elec_erf_integral_tmp(dirac_inverse_list(i,1),dirac_inverse_list(l,1)) -= (dirac_SCF_density_matrix_ao(dirac_inverse_list(j,1),dirac_inverse_list(k,1))) * integral
!     !L_beta L_beta .or. S_beta S_beta
!     dirac_ao_bi_elec_erf_integral_tmp(dirac_inverse_list(i,2),dirac_inverse_list(j,2)) += (dirac_SCF_density_matrix_ao(dirac_inverse_list(l,1),dirac_inverse_list(k,1))) * integral
!     dirac_ao_bi_elec_erf_integral_tmp(dirac_inverse_list(i,2),dirac_inverse_list(j,2)) += (dirac_SCF_density_matrix_ao(dirac_inverse_list(l,2),dirac_inverse_list(k,2))) * integral
!     dirac_ao_bi_elec_erf_integral_tmp(dirac_inverse_list(i,2),dirac_inverse_list(l,2)) -= (dirac_SCF_density_matrix_ao(dirac_inverse_list(j,2),dirac_inverse_list(k,2))) * integral   
!     !L_beta L_alpha .or. S_beta S_alpha
!     dirac_ao_bi_elec_erf_integral_tmp(dirac_inverse_list(i,2),dirac_inverse_list(l,1)) -= (dirac_SCF_density_matrix_ao(dirac_inverse_list(j,2),dirac_inverse_list(k,1))) * integral
!     !L_alpha L_beta .or. S_alpha S_beta
!     dirac_ao_bi_elec_erf_integral_tmp(dirac_inverse_list(i,1),dirac_inverse_list(l,2)) -= (dirac_SCF_density_matrix_ao(dirac_inverse_list(j,1),dirac_inverse_list(k,2))) * integral
!    elseif((i .le. large_ao_num .and. j .le. large_ao_num .and. k .gt. large_ao_num .and. l .gt. large_ao_num) .or. &
!           (i .gt. large_ao_num .and. j .gt. large_ao_num .and. k .le. large_ao_num .and. l .le. large_ao_num))then
!      !L_alpha L_alpha .or. S_alpha S_alpha
!     dirac_ao_bi_elec_erf_integral_tmp(dirac_inverse_list(i,1),dirac_inverse_list(j,1)) += (dirac_SCF_density_matrix_ao(dirac_inverse_list(l,1),dirac_inverse_list(k,1))) * integral
!     dirac_ao_bi_elec_erf_integral_tmp(dirac_inverse_list(i,1),dirac_inverse_list(j,1)) += (dirac_SCF_density_matrix_ao(dirac_inverse_list(l,2),dirac_inverse_list(k,2))) * integral
!     !L_beta L_beta .or. S_beta S_beta
!     dirac_ao_bi_elec_erf_integral_tmp(dirac_inverse_list(i,2),dirac_inverse_list(j,2)) += (dirac_SCF_density_matrix_ao(dirac_inverse_list(l,1),dirac_inverse_list(k,1))) * integral
!     dirac_ao_bi_elec_erf_integral_tmp(dirac_inverse_list(i,2),dirac_inverse_list(j,2)) += (dirac_SCF_density_matrix_ao(dirac_inverse_list(l,2),dirac_inverse_list(k,2))) * integral
!     !L_alpha S_alpha .or S_alpha L_alpha
!     dirac_ao_bi_elec_erf_integral_tmp(dirac_inverse_list(i,1),dirac_inverse_list(l,1)) -= (dirac_SCF_density_matrix_ao(dirac_inverse_list(j,1),dirac_inverse_list(k,1))) * integral
!     !L_beta S_beta .or. S_beta L_beta
!     dirac_ao_bi_elec_erf_integral_tmp(dirac_inverse_list(i,2),dirac_inverse_list(l,2)) -= (dirac_SCF_density_matrix_ao(dirac_inverse_list(j,2),dirac_inverse_list(k,2))) * integral
!     !L_beta S_alpha .or. S_beta L_alpha
!     dirac_ao_bi_elec_erf_integral_tmp(dirac_inverse_list(i,2),dirac_inverse_list(l,1)) -= (dirac_SCF_density_matrix_ao(dirac_inverse_list(j,2),dirac_inverse_list(k,1))) * integral
!     !L_alpha S_beta .or. S_alpha L_beta
!     dirac_ao_bi_elec_erf_integral_tmp(dirac_inverse_list(i,1),dirac_inverse_list(l,2)) -= (dirac_SCF_density_matrix_ao(dirac_inverse_list(j,1),dirac_inverse_list(k,2))) * integral
!    endif
!   enddo
!  enddo
! enddo
!!$OMP END DO NOWAIT
!!$OMP CRITICAL
! dirac_ao_bi_elec_erf_integral += dirac_ao_bi_elec_erf_integral_tmp
!!$OMP END CRITICAL
! deallocate(keys,values,dirac_ao_bi_elec_erf_integral_tmp)
!!$OMP END PARALLEL
!END_PROVIDER


!BEGIN_PROVIDER [ complex*16, dirac_ao_bi_elec_Gaunt_integral, (2*dirac_ao_num, 2*dirac_ao_num) ]
! use map_module
! implicit none
! BEGIN_DOC
! !Array of the bi-electronic Fock matrix for the Coulomb-Gaunt interaction,
! ! in Dirac AO basis set.
! !Take care, the density matrix index
! ! have been correctly inverse, unlike
! ! in the non-relativistic code
! END_DOC
! PROVIDE dirac_ao_bielec_integrals_in_map
! integer                        :: i,j,k,l,k1 
! integer                        :: r,s,p,q
! double precision               :: dirac_ao_bielec_integral, local_threshold
! double precision               :: integral
! complex*16, allocatable        :: dirac_ao_bi_elec_Gaunt_integral_tmp(:,:)
! integer(omp_lock_kind)         :: lck(dirac_ao_num)
! integer*8                      :: i8
! integer                        :: ii(8), jj(8), kk(8), ll(8), k2
! integer(cache_map_size_kind)   :: n_elements_max, n_elements
! integer(key_kind), allocatable :: keys(:)
! double precision, allocatable  :: values(:)
! dirac_ao_bi_elec_Gaunt_integral = (0.d0,0.d0)
!!$OMP PARALLEL DEFAULT(NONE)                                      &
!!$OMP PRIVATE(i,j,l,k1,k,integral,ii,jj,kk,ll,i8,keys,values,n_elements_max, &
!!$OMP  n_elements,dirac_ao_bi_elec_Gaunt_integral_tmp)&
!!$OMP SHARED(dirac_ao_num,dirac_SCF_density_matrix_ao,&
!!$OMP  dirac_ao_integrals_map, dirac_ao_bi_elec_Gaunt_integral) 
! call get_cache_map_n_elements_max(dirac_ao_integrals_map,n_elements_max)
! allocate(keys(n_elements_max), values(n_elements_max))
! allocate(dirac_ao_bi_elec_Gaunt_integral_tmp(2*dirac_ao_num,2*dirac_ao_num))
! dirac_ao_bi_elec_Gaunt_integral_tmp = (0.d0,0.d0)
!!$OMP DO SCHEDULE(dynamic,64)
!!DIR$ NOVECTOR
! do i8=0_8,dirac_ao_integrals_map%map_size
!  n_elements = n_elements_max
!  call get_cache_map(dirac_ao_integrals_map,i8,keys,values,n_elements)
!  do k1=1,n_elements
!   call bielec_integrals_index_reverse(kk,ii,ll,jj,keys(k1))
!   do k2=1,8
!    if (kk(k2)==0) then
!     cycle
!    endif
!    i = ii(k2) ! electron 1
!    j = jj(k2) ! electron 1
!    k = kk(k2) ! electron 2
!    l = ll(k2) ! electron 2
!    ! values(k1) = (ij|kl) <=> <ik|jl>
!    integral = values (k1)
!!!! Coulomb part of the interaction !!!
!    if ((i .le. large_ao_num .and. j .le. large_ao_num .and. k .le. large_ao_num .and. l .le. large_ao_num) .or.  &
!        (i .gt. large_ao_num .and. j .gt. large_ao_num .and. k .gt. large_ao_num .and. l .gt. large_ao_num)) then
!     !L_alpha L_alpha .or. S_alpha S_alpha
!     dirac_ao_bi_elec_Gaunt_integral_tmp(dirac_inverse_list(i,1),dirac_inverse_list(j,1)) += (dirac_SCF_density_matrix_ao(dirac_inverse_list(l,1),dirac_inverse_list(k,1))) * integral 
!     dirac_ao_bi_elec_Gaunt_integral_tmp(dirac_inverse_list(i,1),dirac_inverse_list(j,1)) += (dirac_SCF_density_matrix_ao(dirac_inverse_list(l,2),dirac_inverse_list(k,2))) * integral
!     dirac_ao_bi_elec_Gaunt_integral_tmp(dirac_inverse_list(i,1),dirac_inverse_list(l,1)) -= (dirac_SCF_density_matrix_ao(dirac_inverse_list(j,1),dirac_inverse_list(k,1))) * integral
!     !L_beta L_beta .or. S_beta S_beta
!     dirac_ao_bi_elec_Gaunt_integral_tmp(dirac_inverse_list(i,2),dirac_inverse_list(j,2)) += (dirac_SCF_density_matrix_ao(dirac_inverse_list(l,1),dirac_inverse_list(k,1))) * integral
!     dirac_ao_bi_elec_Gaunt_integral_tmp(dirac_inverse_list(i,2),dirac_inverse_list(j,2)) += (dirac_SCF_density_matrix_ao(dirac_inverse_list(l,2),dirac_inverse_list(k,2))) * integral
!     dirac_ao_bi_elec_Gaunt_integral_tmp(dirac_inverse_list(i,2),dirac_inverse_list(l,2)) -= (dirac_SCF_density_matrix_ao(dirac_inverse_list(j,2),dirac_inverse_list(k,2))) * integral   
!     !L_beta L_alpha .or. S_beta S_alpha
!     dirac_ao_bi_elec_Gaunt_integral_tmp(dirac_inverse_list(i,2),dirac_inverse_list(l,1)) -= (dirac_SCF_density_matrix_ao(dirac_inverse_list(j,2),dirac_inverse_list(k,1))) * integral
!     !L_alpha L_beta .or. S_alpha S_beta
!     dirac_ao_bi_elec_Gaunt_integral_tmp(dirac_inverse_list(i,1),dirac_inverse_list(l,2)) -= (dirac_SCF_density_matrix_ao(dirac_inverse_list(j,1),dirac_inverse_list(k,2))) * integral
!    elseif((i .le. large_ao_num .and. j .le. large_ao_num .and. k .gt. large_ao_num .and. l .gt. large_ao_num) .or. &
!           (i .gt. large_ao_num .and. j .gt. large_ao_num .and. k .le. large_ao_num .and. l .le. large_ao_num))then
!      !L_alpha L_alpha .or. S_alpha S_alpha
!     dirac_ao_bi_elec_Gaunt_integral_tmp(dirac_inverse_list(i,1),dirac_inverse_list(j,1)) += (dirac_SCF_density_matrix_ao(dirac_inverse_list(l,1),dirac_inverse_list(k,1))) * integral
!     dirac_ao_bi_elec_Gaunt_integral_tmp(dirac_inverse_list(i,1),dirac_inverse_list(j,1)) += (dirac_SCF_density_matrix_ao(dirac_inverse_list(l,2),dirac_inverse_list(k,2))) * integral
!     !L_beta L_beta .or. S_beta S_beta
!     dirac_ao_bi_elec_Gaunt_integral_tmp(dirac_inverse_list(i,2),dirac_inverse_list(j,2)) += (dirac_SCF_density_matrix_ao(dirac_inverse_list(l,1),dirac_inverse_list(k,1))) * integral
!     dirac_ao_bi_elec_Gaunt_integral_tmp(dirac_inverse_list(i,2),dirac_inverse_list(j,2)) += (dirac_SCF_density_matrix_ao(dirac_inverse_list(l,2),dirac_inverse_list(k,2))) * integral
!     !L_alpha S_alpha .or S_alpha L_alpha
!     dirac_ao_bi_elec_Gaunt_integral_tmp(dirac_inverse_list(i,1),dirac_inverse_list(l,1)) -= (dirac_SCF_density_matrix_ao(dirac_inverse_list(j,1),dirac_inverse_list(k,1))) * integral
!     !L_beta S_beta .or. S_beta L_beta
!     dirac_ao_bi_elec_Gaunt_integral_tmp(dirac_inverse_list(i,2),dirac_inverse_list(l,2)) -= (dirac_SCF_density_matrix_ao(dirac_inverse_list(j,2),dirac_inverse_list(k,2))) * integral
!     !L_beta S_alpha .or. S_beta L_alpha
!     dirac_ao_bi_elec_Gaunt_integral_tmp(dirac_inverse_list(i,2),dirac_inverse_list(l,1)) -= (dirac_SCF_density_matrix_ao(dirac_inverse_list(j,2),dirac_inverse_list(k,1))) * integral
!     !L_alpha S_beta .or. S_alpha L_beta
!     dirac_ao_bi_elec_Gaunt_integral_tmp(dirac_inverse_list(i,1),dirac_inverse_list(l,2)) -= (dirac_SCF_density_matrix_ao(dirac_inverse_list(j,1),dirac_inverse_list(k,2))) * integral
!    endif
!!!! Gaunt part of the interaction !!!
!    if ((i .le. large_ao_num .and. j .gt. large_ao_num .and. k .gt. large_ao_num .and. l .le. large_ao_num)  .or.  &
!        (i .gt. large_ao_num .and. j .le. large_ao_num .and. k .le. large_ao_num .and. l .gt. large_ao_num)) then  
!     !L_alpha L_alpha .or. S_alpha S_alpha
!     dirac_ao_bi_elec_Gaunt_integral_tmp(dirac_inverse_list(i,1),dirac_inverse_list(l,1)) += 2*(dirac_SCF_density_matrix_ao(dirac_inverse_list(j,2),dirac_inverse_list(k,2))) * integral 
!     dirac_ao_bi_elec_Gaunt_integral_tmp(dirac_inverse_list(i,1),dirac_inverse_list(l,1)) += (dirac_SCF_density_matrix_ao(dirac_inverse_list(j,1),dirac_inverse_list(k,1))) * integral
!     !L_beta L_beta .or. S_beta S_beta
!     dirac_ao_bi_elec_Gaunt_integral_tmp(dirac_inverse_list(i,2),dirac_inverse_list(l,2)) += 2*(dirac_SCF_density_matrix_ao(dirac_inverse_list(j,1),dirac_inverse_list(k,1))) * integral
!     dirac_ao_bi_elec_Gaunt_integral_tmp(dirac_inverse_list(i,2),dirac_inverse_list(l,2)) += (dirac_SCF_density_matrix_ao(dirac_inverse_list(j,2),dirac_inverse_list(k,2))) * integral
!     !L_beta L_alpha .or. S_beta S_alpha
!     dirac_ao_bi_elec_Gaunt_integral_tmp(dirac_inverse_list(i,2),dirac_inverse_list(l,1)) -= (dirac_SCF_density_matrix_ao(dirac_inverse_list(j,2),dirac_inverse_list(k,1))) * integral
!     !L_alpha L_beta .or. S_alpha S_beta
!     dirac_ao_bi_elec_Gaunt_integral_tmp(dirac_inverse_list(i,1),dirac_inverse_list(l,2)) -= (dirac_SCF_density_matrix_ao(dirac_inverse_list(j,1),dirac_inverse_list(k,2))) * integral
!     !L_alpha S_alpha .or S_alpha L_alpha
!     dirac_ao_bi_elec_Gaunt_integral_tmp(dirac_inverse_list(i,1),dirac_inverse_list(j,1)) -= (dirac_SCF_density_matrix_ao(dirac_inverse_list(l,1),dirac_inverse_list(k,1))) * integral
!     dirac_ao_bi_elec_Gaunt_integral_tmp(dirac_inverse_list(i,1),dirac_inverse_list(j,1)) += (dirac_SCF_density_matrix_ao(dirac_inverse_list(l,2),dirac_inverse_list(k,2))) * integral
!     !L_beta S_beta .or S_beta L_beta
!     dirac_ao_bi_elec_Gaunt_integral_tmp(dirac_inverse_list(i,2),dirac_inverse_list(j,2)) -= (dirac_SCF_density_matrix_ao(dirac_inverse_list(l,2),dirac_inverse_list(k,2))) * integral
!     dirac_ao_bi_elec_Gaunt_integral_tmp(dirac_inverse_list(i,2),dirac_inverse_list(j,2)) += (dirac_SCF_density_matrix_ao(dirac_inverse_list(l,1),dirac_inverse_list(k,1))) * integral
!     !L_beta S_alpha .or S_beta L_alpha
!     dirac_ao_bi_elec_Gaunt_integral_tmp(dirac_inverse_list(i,2),dirac_inverse_list(j,1)) -= 2*(dirac_SCF_density_matrix_ao(dirac_inverse_list(l,2),dirac_inverse_list(k,1))) * integral
!     !L_alpha S_beta .or S_alpha L_beta
!     dirac_ao_bi_elec_Gaunt_integral_tmp(dirac_inverse_list(i,1),dirac_inverse_list(j,2)) -= 2*(dirac_SCF_density_matrix_ao(dirac_inverse_list(l,1),dirac_inverse_list(k,2))) * integral
!    elseif ((i .le. large_ao_num .and. j .gt. large_ao_num .and. k .le. large_ao_num .and. l .gt. large_ao_num)  .or.  &
!            (i .gt. large_ao_num .and. j .le. large_ao_num .and. k .gt. large_ao_num .and. l .le. large_ao_num)) then
!     !L_alpha S_alpha .or S_alpha L_alpha
!     dirac_ao_bi_elec_Gaunt_integral_tmp(dirac_inverse_list(i,1),dirac_inverse_list(j,1)) -= (dirac_SCF_density_matrix_ao(dirac_inverse_list(l,1),dirac_inverse_list(k,1))) * integral
!     dirac_ao_bi_elec_Gaunt_integral_tmp(dirac_inverse_list(i,1),dirac_inverse_list(j,1)) += (dirac_SCF_density_matrix_ao(dirac_inverse_list(l,2),dirac_inverse_list(k,2))) * integral
!     dirac_ao_bi_elec_Gaunt_integral_tmp(dirac_inverse_list(i,1),dirac_inverse_list(l,1)) += 2*(dirac_SCF_density_matrix_ao(dirac_inverse_list(j,2),dirac_inverse_list(k,2))) * integral
!     dirac_ao_bi_elec_Gaunt_integral_tmp(dirac_inverse_list(i,1),dirac_inverse_list(l,1)) += (dirac_SCF_density_matrix_ao(dirac_inverse_list(j,1),dirac_inverse_list(k,1))) * integral
!     !L_beta S_beta .or S_beta L_beta
!     dirac_ao_bi_elec_Gaunt_integral_tmp(dirac_inverse_list(i,2),dirac_inverse_list(j,2)) -= (dirac_SCF_density_matrix_ao(dirac_inverse_list(l,2),dirac_inverse_list(k,2))) * integral
!     dirac_ao_bi_elec_Gaunt_integral_tmp(dirac_inverse_list(i,2),dirac_inverse_list(j,2)) += (dirac_SCF_density_matrix_ao(dirac_inverse_list(l,1),dirac_inverse_list(k,1))) * integral
!     dirac_ao_bi_elec_Gaunt_integral_tmp(dirac_inverse_list(i,2),dirac_inverse_list(l,2)) += 2*(dirac_SCF_density_matrix_ao(dirac_inverse_list(j,1),dirac_inverse_list(k,1))) * integral
!     dirac_ao_bi_elec_Gaunt_integral_tmp(dirac_inverse_list(i,2),dirac_inverse_list(l,2)) += (dirac_SCF_density_matrix_ao(dirac_inverse_list(j,2),dirac_inverse_list(k,2))) * integral
!     !L_beta S_alpha .or S_beta L_alpha
!     dirac_ao_bi_elec_Gaunt_integral_tmp(dirac_inverse_list(i,2),dirac_inverse_list(j,1)) -= 2*(dirac_SCF_density_matrix_ao(dirac_inverse_list(l,2),dirac_inverse_list(k,1))) * integral
!     dirac_ao_bi_elec_Gaunt_integral_tmp(dirac_inverse_list(i,2),dirac_inverse_list(l,1)) -= (dirac_SCF_density_matrix_ao(dirac_inverse_list(j,2),dirac_inverse_list(k,1))) * integral
!     !L_alpha S_beta .or S_alpha L_beta
!     dirac_ao_bi_elec_Gaunt_integral_tmp(dirac_inverse_list(i,1),dirac_inverse_list(j,2)) -= 2*(dirac_SCF_density_matrix_ao(dirac_inverse_list(l,1),dirac_inverse_list(k,2))) * integral
!     dirac_ao_bi_elec_Gaunt_integral_tmp(dirac_inverse_list(i,1),dirac_inverse_list(l,2)) -= (dirac_SCF_density_matrix_ao(dirac_inverse_list(j,1),dirac_inverse_list(k,2))) * integral
!    endif
!   enddo
!  enddo
! enddo
!!$OMP END DO NOWAIT
!!$OMP CRITICAL
! dirac_ao_bi_elec_Gaunt_integral += dirac_ao_bi_elec_Gaunt_integral_tmp
!!$OMP END CRITICAL
! deallocate(keys,values,dirac_ao_bi_elec_Gaunt_integral_tmp)
!!$OMP END PARALLEL
!END_PROVIDER

!BEGIN_PROVIDER [ complex*16, dirac_ao_bi_elec_Gaunt_erf_integral, (2*dirac_ao_num, 2*dirac_ao_num) ]
! use map_module
! implicit none
! BEGIN_DOC
! !Array of the bi-electronic Fock matrix for the long-range Coulomb-Gaunt interaction,
! ! in Dirac AO basis set.
! !Take care, the density matrix index
! ! have been correctly inverse, unlike
! ! in the non-relativistic code
! END_DOC
! PROVIDE dirac_ao_bielec_integrals_erf_in_map
! integer                        :: i,j,k,l,k1 
! integer                        :: r,s,p,q
! double precision               :: dirac_ao_bielec_integral_erf, local_threshold
! double precision               :: integral
! complex*16, allocatable        :: dirac_ao_bi_elec_Gaunt_erf_integral_tmp(:,:)
! integer(omp_lock_kind)         :: lck(dirac_ao_num)
! integer*8                      :: i8
! integer                        :: ii(8), jj(8), kk(8), ll(8), k2
! integer(cache_map_size_kind)   :: n_elements_max, n_elements
! integer(key_kind), allocatable :: keys(:)
! double precision, allocatable  :: values(:)
! dirac_ao_bi_elec_Gaunt_erf_integral = (0.d0,0.d0)
!!$OMP PARALLEL DEFAULT(NONE)                                      &
!!$OMP PRIVATE(i,j,l,k1,k,integral,ii,jj,kk,ll,i8,keys,values,n_elements_max, &
!!$OMP  n_elements,dirac_ao_bi_elec_Gaunt_erf_integral_tmp)&
!!$OMP SHARED(dirac_ao_num,dirac_SCF_density_matrix_ao,&
!!$OMP  dirac_ao_integrals_map, dirac_ao_bi_elec_Gaunt_erf_integral) 
! call get_cache_map_n_elements_max(dirac_ao_integrals_erf_map,n_elements_max)
! allocate(keys(n_elements_max), values(n_elements_max))
! allocate(dirac_ao_bi_elec_Gaunt_erf_integral_tmp(2*dirac_ao_num,2*dirac_ao_num))
! dirac_ao_bi_elec_Gaunt_erf_integral_tmp = (0.d0,0.d0)
!!$OMP DO SCHEDULE(dynamic,64)
!!DIR$ NOVECTOR
! do i8=0_8,dirac_ao_integrals_erf_map%map_size
!  n_elements = n_elements_max
!  call get_cache_map(dirac_ao_integrals_erf_map,i8,keys,values,n_elements)
!  do k1=1,n_elements
!   call bielec_integrals_index_reverse(kk,ii,ll,jj,keys(k1))
!   do k2=1,8
!    if (kk(k2)==0) then
!     cycle
!    endif
!    i = ii(k2) ! electron 1
!    j = jj(k2) ! electron 1
!    k = kk(k2) ! electron 2
!    l = ll(k2) ! electron 2
!    ! values(k1) = (ij|kl) <=> <ik|jl>
!    integral = values (k1)
!!!! Coulomb part of the interaction !!!
!    if ((i .le. large_ao_num .and. j .le. large_ao_num .and. k .le. large_ao_num .and. l .le. large_ao_num) .or.  &
!        (i .gt. large_ao_num .and. j .gt. large_ao_num .and. k .gt. large_ao_num .and. l .gt. large_ao_num)) then
!     !L_alpha L_alpha .or. S_alpha S_alpha
!     dirac_ao_bi_elec_Gaunt_erf_integral_tmp(dirac_inverse_list(i,1),dirac_inverse_list(j,1)) += (dirac_SCF_density_matrix_ao(dirac_inverse_list(l,1),dirac_inverse_list(k,1))) * integral 
!     dirac_ao_bi_elec_Gaunt_erf_integral_tmp(dirac_inverse_list(i,1),dirac_inverse_list(j,1)) += (dirac_SCF_density_matrix_ao(dirac_inverse_list(l,2),dirac_inverse_list(k,2))) * integral
!     dirac_ao_bi_elec_Gaunt_erf_integral_tmp(dirac_inverse_list(i,1),dirac_inverse_list(l,1)) -= (dirac_SCF_density_matrix_ao(dirac_inverse_list(j,1),dirac_inverse_list(k,1))) * integral
!     !L_beta L_beta .or. S_beta S_beta
!     dirac_ao_bi_elec_Gaunt_erf_integral_tmp(dirac_inverse_list(i,2),dirac_inverse_list(j,2)) += (dirac_SCF_density_matrix_ao(dirac_inverse_list(l,1),dirac_inverse_list(k,1))) * integral
!     dirac_ao_bi_elec_Gaunt_erf_integral_tmp(dirac_inverse_list(i,2),dirac_inverse_list(j,2)) += (dirac_SCF_density_matrix_ao(dirac_inverse_list(l,2),dirac_inverse_list(k,2))) * integral
!     dirac_ao_bi_elec_Gaunt_erf_integral_tmp(dirac_inverse_list(i,2),dirac_inverse_list(l,2)) -= (dirac_SCF_density_matrix_ao(dirac_inverse_list(j,2),dirac_inverse_list(k,2))) * integral   
!     !L_beta L_alpha .or. S_beta S_alpha
!     dirac_ao_bi_elec_Gaunt_erf_integral_tmp(dirac_inverse_list(i,2),dirac_inverse_list(l,1)) -= (dirac_SCF_density_matrix_ao(dirac_inverse_list(j,2),dirac_inverse_list(k,1))) * integral
!     !L_alpha L_beta .or. S_alpha S_beta
!     dirac_ao_bi_elec_Gaunt_erf_integral_tmp(dirac_inverse_list(i,1),dirac_inverse_list(l,2)) -= (dirac_SCF_density_matrix_ao(dirac_inverse_list(j,1),dirac_inverse_list(k,2))) * integral
!    elseif((i .le. large_ao_num .and. j .le. large_ao_num .and. k .gt. large_ao_num .and. l .gt. large_ao_num) .or. &
!           (i .gt. large_ao_num .and. j .gt. large_ao_num .and. k .le. large_ao_num .and. l .le. large_ao_num))then
!      !L_alpha L_alpha .or. S_alpha S_alpha
!     dirac_ao_bi_elec_Gaunt_erf_integral_tmp(dirac_inverse_list(i,1),dirac_inverse_list(j,1)) += (dirac_SCF_density_matrix_ao(dirac_inverse_list(l,1),dirac_inverse_list(k,1))) * integral
!     dirac_ao_bi_elec_Gaunt_erf_integral_tmp(dirac_inverse_list(i,1),dirac_inverse_list(j,1)) += (dirac_SCF_density_matrix_ao(dirac_inverse_list(l,2),dirac_inverse_list(k,2))) * integral
!     !L_beta L_beta .or. S_beta S_beta
!     dirac_ao_bi_elec_Gaunt_erf_integral_tmp(dirac_inverse_list(i,2),dirac_inverse_list(j,2)) += (dirac_SCF_density_matrix_ao(dirac_inverse_list(l,1),dirac_inverse_list(k,1))) * integral
!     dirac_ao_bi_elec_Gaunt_erf_integral_tmp(dirac_inverse_list(i,2),dirac_inverse_list(j,2)) += (dirac_SCF_density_matrix_ao(dirac_inverse_list(l,2),dirac_inverse_list(k,2))) * integral
!     !L_alpha S_alpha .or S_alpha L_alpha
!     dirac_ao_bi_elec_Gaunt_erf_integral_tmp(dirac_inverse_list(i,1),dirac_inverse_list(l,1)) -= (dirac_SCF_density_matrix_ao(dirac_inverse_list(j,1),dirac_inverse_list(k,1))) * integral
!     !L_beta S_beta .or. S_beta L_beta
!     dirac_ao_bi_elec_Gaunt_erf_integral_tmp(dirac_inverse_list(i,2),dirac_inverse_list(l,2)) -= (dirac_SCF_density_matrix_ao(dirac_inverse_list(j,2),dirac_inverse_list(k,2))) * integral
!     !L_beta S_alpha .or. S_beta L_alpha
!     dirac_ao_bi_elec_Gaunt_erf_integral_tmp(dirac_inverse_list(i,2),dirac_inverse_list(l,1)) -= (dirac_SCF_density_matrix_ao(dirac_inverse_list(j,2),dirac_inverse_list(k,1))) * integral
!     !L_alpha S_beta .or. S_alpha L_beta
!     dirac_ao_bi_elec_Gaunt_erf_integral_tmp(dirac_inverse_list(i,1),dirac_inverse_list(l,2)) -= (dirac_SCF_density_matrix_ao(dirac_inverse_list(j,1),dirac_inverse_list(k,2))) * integral
!    endif
!!!! Gaunt part of the interaction !!!
!    if ((i .le. large_ao_num .and. j .gt. large_ao_num .and. k .gt. large_ao_num .and. l .le. large_ao_num)  .or.  &
!        (i .gt. large_ao_num .and. j .le. large_ao_num .and. k .le. large_ao_num .and. l .gt. large_ao_num)) then  
!     !L_alpha L_alpha .or. S_alpha S_alpha
!     dirac_ao_bi_elec_Gaunt_erf_integral_tmp(dirac_inverse_list(i,1),dirac_inverse_list(l,1)) += 2*(dirac_SCF_density_matrix_ao(dirac_inverse_list(j,2),dirac_inverse_list(k,2))) * integral 
!     dirac_ao_bi_elec_Gaunt_erf_integral_tmp(dirac_inverse_list(i,1),dirac_inverse_list(l,1)) += (dirac_SCF_density_matrix_ao(dirac_inverse_list(j,1),dirac_inverse_list(k,1))) * integral
!     !L_beta L_beta .or. S_beta S_beta
!     dirac_ao_bi_elec_Gaunt_erf_integral_tmp(dirac_inverse_list(i,2),dirac_inverse_list(l,2)) += 2*(dirac_SCF_density_matrix_ao(dirac_inverse_list(j,1),dirac_inverse_list(k,1))) * integral
!     dirac_ao_bi_elec_Gaunt_erf_integral_tmp(dirac_inverse_list(i,2),dirac_inverse_list(l,2)) += (dirac_SCF_density_matrix_ao(dirac_inverse_list(j,2),dirac_inverse_list(k,2))) * integral
!     !L_beta L_alpha .or. S_beta S_alpha
!     dirac_ao_bi_elec_Gaunt_erf_integral_tmp(dirac_inverse_list(i,2),dirac_inverse_list(l,1)) -= (dirac_SCF_density_matrix_ao(dirac_inverse_list(j,2),dirac_inverse_list(k,1))) * integral
!     !L_alpha L_beta .or. S_alpha S_beta
!     dirac_ao_bi_elec_Gaunt_erf_integral_tmp(dirac_inverse_list(i,1),dirac_inverse_list(l,2)) -= (dirac_SCF_density_matrix_ao(dirac_inverse_list(j,1),dirac_inverse_list(k,2))) * integral
!     !L_alpha S_alpha .or S_alpha L_alpha
!     dirac_ao_bi_elec_Gaunt_erf_integral_tmp(dirac_inverse_list(i,1),dirac_inverse_list(j,1)) -= (dirac_SCF_density_matrix_ao(dirac_inverse_list(l,1),dirac_inverse_list(k,1))) * integral
!     dirac_ao_bi_elec_Gaunt_erf_integral_tmp(dirac_inverse_list(i,1),dirac_inverse_list(j,1)) += (dirac_SCF_density_matrix_ao(dirac_inverse_list(l,2),dirac_inverse_list(k,2))) * integral
!     !L_beta S_beta .or S_beta L_beta
!     dirac_ao_bi_elec_Gaunt_erf_integral_tmp(dirac_inverse_list(i,2),dirac_inverse_list(j,2)) -= (dirac_SCF_density_matrix_ao(dirac_inverse_list(l,2),dirac_inverse_list(k,2))) * integral
!     dirac_ao_bi_elec_Gaunt_erf_integral_tmp(dirac_inverse_list(i,2),dirac_inverse_list(j,2)) += (dirac_SCF_density_matrix_ao(dirac_inverse_list(l,1),dirac_inverse_list(k,1))) * integral
!     !L_beta S_alpha .or S_beta L_alpha
!     dirac_ao_bi_elec_Gaunt_erf_integral_tmp(dirac_inverse_list(i,2),dirac_inverse_list(j,1)) -= 2*(dirac_SCF_density_matrix_ao(dirac_inverse_list(l,2),dirac_inverse_list(k,1))) * integral
!     !L_alpha S_beta .or S_alpha L_beta
!     dirac_ao_bi_elec_Gaunt_erf_integral_tmp(dirac_inverse_list(i,1),dirac_inverse_list(j,2)) -= 2*(dirac_SCF_density_matrix_ao(dirac_inverse_list(l,1),dirac_inverse_list(k,2))) * integral
!    elseif ((i .le. large_ao_num .and. j .gt. large_ao_num .and. k .le. large_ao_num .and. l .gt. large_ao_num)  .or.  &
!            (i .gt. large_ao_num .and. j .le. large_ao_num .and. k .gt. large_ao_num .and. l .le. large_ao_num)) then
!     !L_alpha S_alpha .or S_alpha L_alpha
!     dirac_ao_bi_elec_Gaunt_erf_integral_tmp(dirac_inverse_list(i,1),dirac_inverse_list(j,1)) -= (dirac_SCF_density_matrix_ao(dirac_inverse_list(l,1),dirac_inverse_list(k,1))) * integral
!     dirac_ao_bi_elec_Gaunt_erf_integral_tmp(dirac_inverse_list(i,1),dirac_inverse_list(j,1)) += (dirac_SCF_density_matrix_ao(dirac_inverse_list(l,2),dirac_inverse_list(k,2))) * integral
!     dirac_ao_bi_elec_Gaunt_erf_integral_tmp(dirac_inverse_list(i,1),dirac_inverse_list(l,1)) += 2*(dirac_SCF_density_matrix_ao(dirac_inverse_list(j,2),dirac_inverse_list(k,2))) * integral
!     dirac_ao_bi_elec_Gaunt_erf_integral_tmp(dirac_inverse_list(i,1),dirac_inverse_list(l,1)) += (dirac_SCF_density_matrix_ao(dirac_inverse_list(j,1),dirac_inverse_list(k,1))) * integral
!     !L_beta S_beta .or S_beta L_beta
!     dirac_ao_bi_elec_Gaunt_erf_integral_tmp(dirac_inverse_list(i,2),dirac_inverse_list(j,2)) -= (dirac_SCF_density_matrix_ao(dirac_inverse_list(l,2),dirac_inverse_list(k,2))) * integral
!     dirac_ao_bi_elec_Gaunt_erf_integral_tmp(dirac_inverse_list(i,2),dirac_inverse_list(j,2)) += (dirac_SCF_density_matrix_ao(dirac_inverse_list(l,1),dirac_inverse_list(k,1))) * integral
!     dirac_ao_bi_elec_Gaunt_erf_integral_tmp(dirac_inverse_list(i,2),dirac_inverse_list(l,2)) += 2*(dirac_SCF_density_matrix_ao(dirac_inverse_list(j,1),dirac_inverse_list(k,1))) * integral
!     dirac_ao_bi_elec_Gaunt_erf_integral_tmp(dirac_inverse_list(i,2),dirac_inverse_list(l,2)) += (dirac_SCF_density_matrix_ao(dirac_inverse_list(j,2),dirac_inverse_list(k,2))) * integral
!     !L_beta S_alpha .or S_beta L_alpha
!     dirac_ao_bi_elec_Gaunt_erf_integral_tmp(dirac_inverse_list(i,2),dirac_inverse_list(j,1)) -= 2*(dirac_SCF_density_matrix_ao(dirac_inverse_list(l,2),dirac_inverse_list(k,1))) * integral
!     dirac_ao_bi_elec_Gaunt_erf_integral_tmp(dirac_inverse_list(i,2),dirac_inverse_list(l,1)) -= (dirac_SCF_density_matrix_ao(dirac_inverse_list(j,2),dirac_inverse_list(k,1))) * integral
!     !L_alpha S_beta .or S_alpha L_beta
!     dirac_ao_bi_elec_Gaunt_erf_integral_tmp(dirac_inverse_list(i,1),dirac_inverse_list(j,2)) -= 2*(dirac_SCF_density_matrix_ao(dirac_inverse_list(l,1),dirac_inverse_list(k,2))) * integral
!     dirac_ao_bi_elec_Gaunt_erf_integral_tmp(dirac_inverse_list(i,1),dirac_inverse_list(l,2)) -= (dirac_SCF_density_matrix_ao(dirac_inverse_list(j,1),dirac_inverse_list(k,2))) * integral
!    endif
!   enddo
!  enddo
! enddo
!!$OMP END DO NOWAIT
!!$OMP CRITICAL
! dirac_ao_bi_elec_Gaunt_erf_integral += dirac_ao_bi_elec_Gaunt_erf_integral_tmp
!!$OMP END CRITICAL
! deallocate(keys,values,dirac_ao_bi_elec_Gaunt_erf_integral_tmp)
!!$OMP END PARALLEL
!END_PROVIDER



!BEGIN_PROVIDER [ complex*16, dirac_Fock_matrix_Coulomb_ao, (2*dirac_ao_num,2*dirac_ao_num) ]
! implicit none
! BEGIN_DOC
! !Dirac Fock matrix in AO basis set for the Coulomb bi-eletronic interaction
! END_DOC
! integer                        :: i,j
! do j=1,2*dirac_ao_num
!  do i=1,2*dirac_ao_num
!   dirac_Fock_matrix_Coulomb_ao(i,j) = dirac_ao_mono_elec_integral(i,j) + dirac_ao_bi_elec_integral(i,j)
!  enddo
! enddo
!END_PROVIDER

!BEGIN_PROVIDER [ complex*16, dirac_Fock_matrix_Coulomb_erf_ao, (2*dirac_ao_num,2*dirac_ao_num) ]
! implicit none
! BEGIN_DOC
! !Dirac Fock matrix in AO basis set for the long-range Coulomb bi-eletronic interaction
! END_DOC
! integer                        :: i,j
! do j=1,2*dirac_ao_num
!  do i=1,2*dirac_ao_num
!   dirac_Fock_matrix_Coulomb_erf_ao(i,j) = dirac_ao_mono_elec_integral(i,j) + dirac_ao_bi_elec_erf_integral(i,j)
!  enddo
! enddo
!END_PROVIDER

!
!BEGIN_PROVIDER [ complex*16, dirac_Fock_matrix_Coulomb_Gaunt_ao, (2*dirac_ao_num,2*dirac_ao_num) ]
! implicit none
! BEGIN_DOC
! !Dirac Fock matrix in AO basis set for the Coulomb-Gaunt bi-electronic interaction
! END_DOC
! integer                        :: i,j
! do j=1,2*dirac_ao_num
!  do i=1,2*dirac_ao_num
!   dirac_Fock_matrix_Coulomb_Gaunt_ao(i,j) = dirac_ao_mono_elec_integral(i,j) + dirac_ao_bi_elec_Gaunt_integral(i,j)
!  enddo
! enddo
!END_PROVIDER

!BEGIN_PROVIDER [ complex*16, dirac_Fock_matrix_Coulomb_Gaunt_erf_ao, (2*dirac_ao_num,2*dirac_ao_num) ]
! implicit none
! BEGIN_DOC
! !Dirac Fock matrix in AO basis set for the long-range Coulomb-Gaunt bi-electronic interaction
! END_DOC
! integer                        :: i,j
! do j=1,2*dirac_ao_num
!  do i=1,2*dirac_ao_num
!   dirac_Fock_matrix_Coulomb_Gaunt_erf_ao(i,j) = dirac_ao_mono_elec_integral(i,j) + dirac_ao_bi_elec_Gaunt_erf_integral(i,j)
!  enddo
! enddo
!END_PROVIDER



!BEGIN_PROVIDER [ double precision, dirac_extra_energy_contrib_from_density]
!implicit none
! BEGIN_DOC
! !Contribution from electronic density
! END_DOC
! dirac_extra_energy_contrib_from_density = 0.d0
!END_PROVIDER



!BEGIN_PROVIDER [ complex*16, dirac_HF_one_electron_energy_complex]
!&BEGIN_PROVIDER [ double precision, dirac_HF_one_electron_energy]
! implicit none
! BEGIN_DOC
! !One-electron energy of the Nucleus-Electron interaction
! !The energy is supposed to be a real, thus we check for its complex part to be
! ! a VERY small artifact and take only its real part
! END_DOC
! integer :: i,j
! dirac_HF_one_electron_energy_complex = (0.d0,0.d0)
! do j=1, 2*dirac_ao_num
!  do i=1, 2*dirac_ao_num
!   dirac_HF_one_electron_energy_complex += dirac_ao_mono_elec_integral(i,j)* dirac_SCF_density_matrix_ao(j,i) 
!  enddo
! enddo
! dirac_HF_one_electron_energy = real(dirac_HF_one_electron_energy_complex)
! if (aimag(dirac_HF_one_electron_energy_complex) .gt. 1.d-10 ) then
! print*, 'Warning! The energy is not real'
! print*, 'dirac_HF_one_electron_energy_complex =',dirac_HF_one_electron_energy_complex
! STOP
! endif
!END_PROVIDER


!
!BEGIN_PROVIDER [ complex*16, dirac_HF_two_electron_Coulomb_energy_complex]
!&BEGIN_PROVIDER [ double precision, dirac_HF_two_electron_Coulomb_energy] 
! implicit none
! BEGIN_DOC
! !Two-electrons energy of the Coulomb ee interaction
! !The energy is supposed to be a real, thus we check for its complex part to be
! ! a VERY small artifact and take only its real part
! END_DOC
! integer :: i,j
! dirac_HF_two_electron_Coulomb_energy_complex = (0.d0,0.d0)
! do j=1, 2*dirac_ao_num
!  do i=1, 2*dirac_ao_num
!   dirac_HF_two_electron_Coulomb_energy_complex += 0.5d0* (dirac_ao_bi_elec_integral(i,j)) * dirac_SCF_density_matrix_ao(j,i)
!  enddo
! enddo
! dirac_HF_two_electron_Coulomb_energy = real(dirac_HF_two_electron_Coulomb_energy_complex)
! if (aimag(dirac_HF_two_electron_Coulomb_energy_complex) .gt. 1.d-10) then
! print*, 'Warning! The energy is not real'
! print*, 'dirac_HF_two_electron_Coulomb_energy_complex =', dirac_HF_two_electron_Coulomb_energy_complex
! STOP
! endif
!END_PROVIDER

!BEGIN_PROVIDER [ complex*16, dirac_HF_two_electron_Coulomb_erf_energy_complex]
!&BEGIN_PROVIDER [ double precision, dirac_HF_two_electron_Coulomb_erf_energy] 
! implicit none
! BEGIN_DOC
! !Two-electrons energy of the Coulomb ee long-range interaction
! !The energy is supposed to be a real, thus we check for its complex part to be
! ! a VERY small artifact and take only its real part
! END_DOC
! integer :: i,j
! dirac_HF_two_electron_Coulomb_erf_energy_complex = (0.d0,0.d0)
! do j=1, 2*dirac_ao_num
!  do i=1, 2*dirac_ao_num
!   dirac_HF_two_electron_Coulomb_erf_energy_complex += 0.5d0* (dirac_ao_bi_elec_erf_integral(i,j)) * dirac_SCF_density_matrix_ao(j,i)
!  enddo
! enddo
! dirac_HF_two_electron_Coulomb_erf_energy = real(dirac_HF_two_electron_Coulomb_erf_energy_complex)
! if (aimag(dirac_HF_two_electron_Coulomb_erf_energy_complex) .gt. 1.d-10) then
! print*, 'Warning! The energy is not real'
! print*, 'dirac_HF_two_electron_Coulomb_erf_energy_complex =', dirac_HF_two_electron_Coulomb_erf_energy_complex
! STOP
! endif
!END_PROVIDER


!BEGIN_PROVIDER [ complex*16, dirac_HF_two_electron_Coulomb_Gaunt_energy_complex]
!&BEGIN_PROVIDER [ double precision, dirac_HF_two_electron_Coulomb_Gaunt_energy] 
! implicit none
! BEGIN_DOC
! !Two-electrons energy of the Coulomb_Gaunt ee interaction
! !The energy is supposed to be a real, thus we check for its complex part to be
! ! a VERY small artifact and take only its real part
! END_DOC
! integer :: i,j
! dirac_HF_two_electron_Coulomb_Gaunt_energy_complex = (0.d0,0.d0)
! do j=1, 2*dirac_ao_num
!  do i=1, 2*dirac_ao_num
!   dirac_HF_two_electron_Coulomb_Gaunt_energy_complex += 0.5d0* (dirac_ao_bi_elec_Gaunt_integral(i,j)) * dirac_SCF_density_matrix_ao(j,i)
!  enddo
! enddo
! dirac_HF_two_electron_Coulomb_Gaunt_energy = real(dirac_HF_two_electron_Coulomb_Gaunt_energy_complex)
! if (aimag(dirac_HF_two_electron_Coulomb_Gaunt_energy_complex) .gt. 1.d-10) then
! print*, 'Warning! The energy is not real'
! print*, 'dirac_HF_two_electron_Coulomb_Gaunt_energy_complex =', dirac_HF_two_electron_Coulomb_Gaunt_energy_complex
! STOP
! endif
!END_PROVIDER

!BEGIN_PROVIDER [ complex*16, dirac_HF_two_electron_Coulomb_Gaunt_erf_energy_complex]
!&BEGIN_PROVIDER [ double precision, dirac_HF_two_electron_Coulomb_Gaunt_erf_energy] 
! implicit none
! BEGIN_DOC
! !Two-electrons energy of the long-range Coulomb_Gaunt ee interaction
! !The energy is supposed to be a real, thus we check for its complex part to be
! ! a VERY small artifact and take only its real part
! END_DOC
! integer :: i,j
! dirac_HF_two_electron_Coulomb_Gaunt_erf_energy_complex = (0.d0,0.d0)
! do j=1, 2*dirac_ao_num
!  do i=1, 2*dirac_ao_num
!   dirac_HF_two_electron_Coulomb_Gaunt_erf_energy_complex += 0.5d0* (dirac_ao_bi_elec_Gaunt_erf_integral(i,j)) * dirac_SCF_density_matrix_ao(j,i)
!  enddo
! enddo
! dirac_HF_two_electron_Coulomb_Gaunt_erf_energy = real(dirac_HF_two_electron_Coulomb_Gaunt_erf_energy_complex)
! if (aimag(dirac_HF_two_electron_Coulomb_Gaunt_erf_energy_complex) .gt. 1.d-10) then
! print*, 'Warning! The energy is not real'
! print*, 'dirac_HF_two_electron_Coulomb_Gaunt_erf_energy_complex =', dirac_HF_two_electron_Coulomb_Gaunt_erf_energy_complex
! STOP
! endif
!END_PROVIDER



!BEGIN_PROVIDER [ double precision, dirac_HF_Coulomb_energy]
! implicit none
! BEGIN_DOC
! !Dirac-Hartree-Fock energy for the Coulomb ee interaction
! END_DOC
! dirac_HF_Coulomb_energy = nuclear_repulsion + dirac_HF_two_electron_Coulomb_energy + dirac_HF_one_electron_energy
!END_PROVIDER

!BEGIN_PROVIDER [ double precision, dirac_HF_Coulomb_erf_energy]
! implicit none
! BEGIN_DOC
! !Dirac-Hartree-Fock energy for the long-range Coulomb ee interaction
! END_DOC
! dirac_HF_Coulomb_erf_energy = nuclear_repulsion + dirac_HF_two_electron_Coulomb_erf_energy + dirac_HF_one_electron_energy
!END_PROVIDER


!BEGIN_PROVIDER [ double precision, dirac_HF_Coulomb_Gaunt_energy]
! implicit none
! BEGIN_DOC
! !Dirac-Hartree-Fock energy for the Coulomb_Gaunt ee interaction
! END_DOC
! dirac_HF_Coulomb_Gaunt_energy = nuclear_repulsion + dirac_HF_two_electron_Coulomb_Gaunt_energy + dirac_HF_one_electron_energy
!END_PROVIDER

!BEGIN_PROVIDER [ double precision, dirac_HF_Coulomb_Gaunt_erf_energy]
! implicit none
! BEGIN_DOC
! !Dirac-Hartree-Fock energy for the long-range Coulomb_Gaunt ee interaction
! END_DOC
! dirac_HF_Coulomb_Gaunt_erf_energy = nuclear_repulsion + dirac_HF_two_electron_Coulomb_Gaunt_erf_energy + dirac_HF_one_electron_energy
!END_PROVIDER



!BEGIN_PROVIDER [ double precision, dirac_SCF_Coulomb_energy]
! implicit none 
! BEGIN_DOC
! !Dirac_SCF energy for a Coulomb ee interaction
! END_DOC 
! dirac_SCF_Coulomb_energy = dirac_HF_Coulomb_energy + dirac_extra_energy_contrib_from_density
!END_PROVIDER

!BEGIN_PROVIDER [ double precision, dirac_SCF_Coulomb_erf_energy]
! implicit none 
! BEGIN_DOC
! !Dirac_SCF energy for a long-range Coulomb ee interaction
! END_DOC 
! dirac_SCF_Coulomb_erf_energy = dirac_HF_Coulomb_erf_energy + dirac_extra_energy_contrib_from_density
!END_PROVIDER


!BEGIN_PROVIDER [ double precision, dirac_SCF_Coulomb_Gaunt_energy]
! implicit none
! BEGIN_DOC
! !Dirac_SCF energy for a Coulomb_Gaunt ee interaction
! END_DOC 
! dirac_SCF_Coulomb_Gaunt_energy = dirac_HF_Coulomb_Gaunt_energy + dirac_extra_energy_contrib_from_density
!END_PROVIDER

!BEGIN_PROVIDER [ double precision, dirac_SCF_Coulomb_Gaunt_erf_energy]
! implicit none
! BEGIN_DOC
! !Dirac_SCF energy for a long-range Coulomb_Gaunt ee interaction
! END_DOC 
! dirac_SCF_Coulomb_Gaunt_erf_energy = dirac_HF_Coulomb_Gaunt_erf_energy + dirac_extra_energy_contrib_from_density
!END_PROVIDER


!BEGIN_PROVIDER [complex*16, dirac_Fock_matrix_Coulomb_mo,(2*dirac_mo_tot_num,2*dirac_mo_tot_num)]
! implicit none
! BEGIN_DOC
! ! Fock matrix in the MO basis for a Coulomb ee interaction
! END_DOC
!   call dirac_ao_to_mo(                                              &
!       dirac_Fock_matrix_Coulomb_ao,                                         &
!       size(dirac_Fock_matrix_Coulomb_ao,1),                                 &
!       dirac_Fock_matrix_Coulomb_mo,                                         &
!       size(dirac_Fock_matrix_Coulomb_mo,1)                                  &
!       )
!END_PROVIDER

!BEGIN_PROVIDER [complex*16, dirac_Fock_matrix_Coulomb_erf_mo,(2*dirac_mo_tot_num,2*dirac_mo_tot_num)]
! implicit none
! BEGIN_DOC
! ! Fock matrix in the MO basis for a long-range Coulomb ee interaction
! END_DOC
!   call dirac_ao_to_mo(                                              &
!       dirac_Fock_matrix_Coulomb_erf_ao,                                         &
!       size(dirac_Fock_matrix_Coulomb_ao,1),                                 &
!       dirac_Fock_matrix_Coulomb_erf_mo,                                         &
!       size(dirac_Fock_matrix_Coulomb_mo,1)                                  &
!       )
!END_PROVIDER
!

!BEGIN_PROVIDER [complex*16, dirac_Fock_matrix_Coulomb_Gaunt_mo,(2*dirac_mo_tot_num,2*dirac_mo_tot_num)]
! implicit none
! BEGIN_DOC
! ! Fock matrix in the MO basis for a Coulomb_Gaunt ee interaction
! END_DOC
!   call dirac_ao_to_mo(                                              &
!       dirac_Fock_matrix_Coulomb_Gaunt_ao,                                         &
!       size(dirac_Fock_matrix_Coulomb_Gaunt_ao,1),                                 &
!       dirac_Fock_matrix_Coulomb_Gaunt_mo,                                         &
!       size(dirac_Fock_matrix_Coulomb_Gaunt_mo,1)                                  &
!       )
!END_PROVIDER

!BEGIN_PROVIDER [complex*16, dirac_Fock_matrix_Coulomb_Gaunt_erf_mo,(2*dirac_mo_tot_num,2*dirac_mo_tot_num)]
! implicit none
! BEGIN_DOC
! ! Fock matrix in the MO basis for a long-range Coulomb_Gaunt ee interaction
! END_DOC
!   call dirac_ao_to_mo(                                              &
!       dirac_Fock_matrix_Coulomb_Gaunt_erf_ao,                                         &
!       size(dirac_Fock_matrix_Coulomb_Gaunt_ao,1),                                 &
!       dirac_Fock_matrix_Coulomb_Gaunt_erf_mo,                                         &
!       size(dirac_Fock_matrix_Coulomb_Gaunt_mo,1)                                  &
!       )
!END_PROVIDER


!
!BEGIN_PROVIDER [complex*16, dirac_Fock_matrix_diag_Coulomb_mo,(2*dirac_mo_tot_num)]
! implicit none
! BEGIN_DOC
! !Diagonal of the Fock matrix in the MO basis for a Coulomb ee interaction
! END_DOC
! integer :: i
! do i = 1, 2*dirac_mo_tot_num
!  dirac_Fock_matrix_diag_Coulomb_mo(i) = dirac_Fock_matrix_Coulomb_mo (i,i)
! enddo
!END_PROVIDER

!BEGIN_PROVIDER [complex*16, dirac_Fock_matrix_diag_Coulomb_erf_mo,(2*dirac_mo_tot_num)]
! implicit none
! BEGIN_DOC
! !Diagonal of the Fock matrix in the MO basis for a long-range Coulomb ee interaction
! END_DOC
! integer :: i
! do i = 1, 2*dirac_mo_tot_num
!  dirac_Fock_matrix_diag_Coulomb_erf_mo(i) = dirac_Fock_matrix_Coulomb_erf_mo (i,i)
! enddo
!END_PROVIDER


!BEGIN_PROVIDER [complex*16, dirac_Fock_matrix_diag_Coulomb_Gaunt_mo,(2*dirac_mo_tot_num)]
! implicit none
! BEGIN_DOC
! !Diagonal of the Fock matrix in the MO basis for a Coulomt_Gaunt ee
! !interaction
! END_DOC
! integer :: i
! do i = 1, 2*dirac_mo_tot_num
!  dirac_Fock_matrix_diag_Coulomb_Gaunt_mo(i) = dirac_Fock_matrix_Coulomb_Gaunt_mo (i,i)
! enddo
!END_PROVIDER

!BEGIN_PROVIDER [complex*16, dirac_Fock_matrix_diag_Coulomb_Gaunt_erf_mo,(2*dirac_mo_tot_num)]
! implicit none
! BEGIN_DOC
! !Diagonal of the Fock matrix in the MO basis for a long-range Coulomt_Gaunt ee
! !interaction
! END_DOC
! integer :: i
! do i = 1, 2*dirac_mo_tot_num
!  dirac_Fock_matrix_diag_Coulomb_Gaunt_erf_mo(i) = dirac_Fock_matrix_Coulomb_Gaunt_erf_mo (i,i)
! enddo
!END_PROVIDER



!BEGIN_PROVIDER [double precision,eigenvalues_dirac_fock_matrix_Coulomb_mo,(2*(dirac_mo_tot_num))]
!&BEGIN_PROVIDER [complex*16, eigenvectors_dirac_fock_matrix_Coulomb_mo,(2*(dirac_mo_tot_num),2*(dirac_mo_tot_num))]
! implicit none
! BEGIN_DOC
! !The eigenvalues and eigenvectors in the MO basis for a Coulomb ee
! ! interaction
! END_DOC
! integer :: n,nmax
! double precision :: eigenvalues( 2*(dirac_mo_tot_num))
! complex*16       :: eigenvectors(2*(dirac_mo_tot_num),2*(dirac_mo_tot_num))
! n = 2*(dirac_mo_tot_num)
! nmax = n
! call lapack_diag_complex(eigenvalues,eigenvectors,dirac_Fock_matrix_Coulomb_mo,nmax,n)
! eigenvalues_dirac_fock_matrix_Coulomb_mo = eigenvalues
! eigenvectors_dirac_fock_matrix_Coulomb_mo = eigenvectors
!END_PROVIDER

!BEGIN_PROVIDER [double precision,eigenvalues_dirac_fock_matrix_Coulomb_erf_mo,(2*(dirac_mo_tot_num))]
!&BEGIN_PROVIDER [complex*16, eigenvectors_dirac_fock_matrix_Coulomb_erf_mo,(2*(dirac_mo_tot_num),2*(dirac_mo_tot_num))]
! implicit none
! BEGIN_DOC
! !The eigenvalues and eigenvectors in the MO basis for a long-range Coulomb ee
! ! interaction
! END_DOC
! integer :: n,nmax
! double precision :: eigenvalues( 2*(dirac_mo_tot_num))
! complex*16       :: eigenvectors(2*(dirac_mo_tot_num),2*(dirac_mo_tot_num))
! n = 2*(dirac_mo_tot_num)
! nmax = n
! call lapack_diag_complex(eigenvalues,eigenvectors,dirac_Fock_matrix_Coulomb_erf_mo,nmax,n)
! eigenvalues_dirac_fock_matrix_Coulomb_erf_mo = eigenvalues
! eigenvectors_dirac_fock_matrix_Coulomb_erf_mo = eigenvectors
!END_PROVIDER


!BEGIN_PROVIDER [double precision,eigenvalues_dirac_fock_matrix_Coulomb_Gaunt_mo,(2*(dirac_mo_tot_num))]
!&BEGIN_PROVIDER [complex*16, eigenvectors_dirac_fock_matrix_Coulomb_Gaunt_mo,(2*(dirac_mo_tot_num),2*(dirac_mo_tot_num))]
!BEGIN_DOC
!!The eigenvalues and eigenvectors in the MO basis for a Coulomb_Gaunt ee
!! interaction
!END_DOC
! implicit none
! integer :: n,nmax
! double precision :: eigenvalues( 2*(dirac_mo_tot_num))
! complex*16       :: eigenvectors(2*(dirac_mo_tot_num),2*(dirac_mo_tot_num))
! n = 2*(dirac_mo_tot_num)
! nmax = n
! call lapack_diag_complex(eigenvalues,eigenvectors,dirac_Fock_matrix_Coulomb_Gaunt_mo,nmax,n)
! eigenvalues_dirac_fock_matrix_Coulomb_Gaunt_mo = eigenvalues
! eigenvectors_dirac_fock_matrix_Coulomb_Gaunt_mo = eigenvectors
!END_PROVIDER

!BEGIN_PROVIDER [double precision,eigenvalues_dirac_fock_matrix_Coulomb_Gaunt_erf_mo,(2*(dirac_mo_tot_num))]
!&BEGIN_PROVIDER [complex*16, eigenvectors_dirac_fock_matrix_Coulomb_Gaunt_erf_mo,(2*(dirac_mo_tot_num),2*(dirac_mo_tot_num))]
!BEGIN_DOC
!!The eigenvalues and eigenvectors in the MO basis for a long-range Coulomb_Gaunt ee
!! interaction
!END_DOC
! implicit none
! integer :: n,nmax
! double precision :: eigenvalues( 2*(dirac_mo_tot_num))
! complex*16       :: eigenvectors(2*(dirac_mo_tot_num),2*(dirac_mo_tot_num))
! n = 2*(dirac_mo_tot_num)
! nmax = n
! call lapack_diag_complex(eigenvalues,eigenvectors,dirac_Fock_matrix_Coulomb_Gaunt_erf_mo,nmax,n)
! eigenvalues_dirac_fock_matrix_Coulomb_Gaunt_erf_mo = eigenvalues
! eigenvectors_dirac_fock_matrix_Coulomb_Gaunt_erf_mo = eigenvectors
!END_PROVIDER



!BEGIN_PROVIDER [complex*16, eigenvectors_dirac_Fock_matrix_Coulomb_ao, (2*(dirac_mo_tot_num),2*(dirac_mo_tot_num))]
!implicit none
!BEGIN_DOC
!!The eigenvectors in the AO basis, which does not diagonalize S, 
!! for a Coulomb ee interaction
!END_DOC
!integer :: n,nmax
! call zgemm('N','N', 2*(dirac_ao_num), 2*(dirac_mo_tot_num), 2*(dirac_ao_num),              &
!     (1.d0,0.d0), dirac_mo_coef_S,size(dirac_mo_coef_S,1),                                      &
!     eigenvectors_dirac_Fock_matrix_Coulomb_mo, size(eigenvectors_dirac_Fock_matrix_Coulomb_mo,1),              &
!     (0.d0,0.d0), eigenvectors_dirac_Fock_matrix_Coulomb_ao, size(eigenvectors_dirac_Fock_matrix_Coulomb_ao,1)) 
!END_PROVIDER

!BEGIN_PROVIDER [complex*16, eigenvectors_dirac_Fock_matrix_Coulomb_erf_ao, (2*(dirac_mo_tot_num),2*(dirac_mo_tot_num))]
!implicit none
!BEGIN_DOC
!!The eigenvectors in the AO basis, which does not diagonalize S, 
!! for a long-range Coulomb ee interaction
!END_DOC
!integer :: n,nmax
! call zgemm('N','N', 2*(dirac_ao_num), 2*(dirac_mo_tot_num), 2*(dirac_ao_num),              &
!     (1.d0,0.d0), dirac_mo_coef_S,size(dirac_mo_coef_S,1),                                      &
!     eigenvectors_dirac_Fock_matrix_Coulomb_erf_mo, size(eigenvectors_dirac_Fock_matrix_Coulomb_erf_mo,1),              &
!     (0.d0,0.d0), eigenvectors_dirac_Fock_matrix_Coulomb_erf_ao, size(eigenvectors_dirac_Fock_matrix_Coulomb_erf_ao,1)) 
!END_PROVIDER


!BEGIN_PROVIDER [complex*16, eigenvectors_dirac_Fock_matrix_Coulomb_Gaunt_ao, (2*(dirac_mo_tot_num),2*(dirac_mo_tot_num))]
!implicit none
!BEGIN_DOC
!!The eigenvectors in the AO basis, which does not diagonalize S,
!! for a Coulomb_Gaunt ee interaction
!END_DOC
!integer :: n,nmax
! call zgemm('N','N', 2*(dirac_ao_num), 2*(dirac_mo_tot_num), 2*(dirac_ao_num),              &
!     (1.d0,0.d0), dirac_mo_coef_S,size(dirac_mo_coef_S,1),                                      &
!     eigenvectors_dirac_Fock_matrix_Coulomb_Gaunt_mo, size(eigenvectors_dirac_Fock_matrix_Coulomb_Gaunt_mo,1),              &
!     (0.d0,0.d0), eigenvectors_dirac_Fock_matrix_Coulomb_Gaunt_ao, size(eigenvectors_dirac_Fock_matrix_Coulomb_Gaunt_ao,1)) 
!END_PROVIDER

!BEGIN_PROVIDER [complex*16, eigenvectors_dirac_Fock_matrix_Coulomb_Gaunt_erf_ao, (2*(dirac_mo_tot_num),2*(dirac_mo_tot_num))]
!implicit none
!BEGIN_DOC
!!The eigenvectors in the AO basis, which does not diagonalize S,
!! for a long-range Coulomb_Gaunt ee interaction
!END_DOC
!integer :: n,nmax
! call zgemm('N','N', 2*(dirac_ao_num), 2*(dirac_mo_tot_num), 2*(dirac_ao_num),              &
!     (1.d0,0.d0), dirac_mo_coef_S,size(dirac_mo_coef_S,1),                                      &
!     eigenvectors_dirac_Fock_matrix_Coulomb_Gaunt_erf_mo, size(eigenvectors_dirac_Fock_matrix_Coulomb_Gaunt_erf_mo,1),              &
!     (0.d0,0.d0), eigenvectors_dirac_Fock_matrix_Coulomb_Gaunt_erf_ao, size(eigenvectors_dirac_Fock_matrix_Coulomb_Gaunt_erf_ao,1)) 
!END_PROVIDER
