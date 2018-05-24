 BEGIN_PROVIDER [ complex*16, dirac_ao_bi_elec_integral, (2*dirac_ao_num, 2*dirac_ao_num) ]
  use map_module
  implicit none
  BEGIN_DOC
  ! Fock matrix in Dirac AO basis set
  END_DOC
  integer                        :: i,j,k,l,k1,r,s
  integer                        :: i0,j0,k0,l0
  integer                        :: p,q
  complex*16               :: integral, c0
  double precision               :: dirac_ao_bielec_integral, local_threshold
  double precision, allocatable  :: dirac_ao_bi_elec_integral_tmp(:,:)
  dirac_ao_bi_elec_integral = (0.d0,0.d0)
  if (do_direct_integrals) then
 !$OMP PARALLEL DEFAULT(NONE)                                      &
 !$OMP PRIVATE(i,j,l,k1,k,integral,ii,jj,kk,ll,i8,keys,values,p,q,r,s,i0,j0,k0,l0, &
 !$OMP dirac_ao_bi_elec_integral_tmp, c0, &
 !$OMP local_threshold)&
 !$OMP SHARED(dirac_ao_num,dirac_SCF_density_matrix_ao,&
 !$OMP dirac_ao_integrals_map,dirac_ao_integrals_threshold, dirac_ao_bielec_integral_schwartz, &
 !$OMP dirac_ao_overlap_abs, dirac_ao_bi_elec_integral)
   allocate(keys(1), values(1))
   allocate(dirac_ao_bi_elec_integral_tmp(dirac_ao_num,dirac_ao_num))
   dirac_ao_bi_elec_integral_tmp = 0.d0
   q = dirac_ao_num*dirac_ao_num*dirac_ao_num*dirac_ao_num
 !$OMP DO SCHEDULE(dynamic)
   do p=1_8,q
    call bielec_integrals_index_reverse(kk,ii,ll,jj,p)
    if ( (kk(1)>dirac_ao_num) .or. &
       (ii(1)>dirac_ao_num) .or. &
       (jj(1)>dirac_ao_num) .or. &
       (ll(1)>dirac_ao_num) ) then
     cycle
    endif
    k = kk(1)
    i = ii(1)
    l = ll(1)
    j = jj(1)
    if (dirac_ao_overlap_abs(k,l)*dirac_ao_overlap_abs(i,j) < dirac_ao_integrals_threshold) then
     cycle
    endif
    local_threshold = dirac_ao_bielec_integral_schwartz(k,l)*dirac_ao_bielec_integral_schwartz(i,j)
    if (local_threshold < dirac_ao_integrals_threshold) then
     cycle
    endif
    i0 = i
    j0 = j
    k0 = k
    l0 = l
    values(1) = 0.d0
    local_threshold = dirac_ao_integrals_threshold/local_threshold
    do k2=1,8
     if (kk(k2)==0) then
      cycle
     endif
     i = ii(k2)
     j = jj(k2)
     k = kk(k2)
     l = ll(k2)
     c0 = dirac_SCF_density_matrix_ao(k,l)
     if ( zabs(c0)  < local_threshold) then
      cycle
     endif
     if (values(1) == 0.d0) then
      values(1) = dirac_ao_bielec_integral(k0,l0,i0,j0)
     endif
     integral = c0 * values(1)
     dirac_ao_bi_elec_integral_tmp(i,j) += integral
    enddo
   enddo
 !$OMP END DO NOWAIT
 !$OMP CRITICAL
  dirac_ao_bi_elec_integral +=(1.d0,0.d0)* dirac_ao_bi_elec_integral_tmp
 !$OMP END CRITICAL
  deallocate(keys,values,dirac_ao_bi_elec_integral_tmp)
 !$OMP END PARALLEL
  else
   PROVIDE dirac_ao_bielec_integrals_in_map
   integer(omp_lock_kind) :: lck(dirac_ao_num)
   integer*8                      :: i8
   integer                        :: ii(8), jj(8), kk(8), ll(8), k2
   integer(cache_map_size_kind)   :: n_elements_max, n_elements
   integer(key_kind), allocatable :: keys(:)
   double precision, allocatable  :: values(:)
 !$OMP PARALLEL DEFAULT(NONE)                                      &
 !$OMP PRIVATE(i,j,l,k1,k,integral,ii,jj,kk,ll,i8,keys,values,n_elements_max, &
 !$OMP  n_elements,dirac_ao_bi_elec_integral_tmp)&
 !$OMP SHARED(dirac_ao_num,dirac_SCF_density_matrix_ao,&
 !$OMP  dirac_ao_integrals_map, dirac_ao_bi_elec_integral) 
   call get_cache_map_n_elements_max(dirac_ao_integrals_map,n_elements_max)
   allocate(keys(n_elements_max), values(n_elements_max))
   allocate(dirac_ao_bi_elec_integral_tmp(dirac_ao_num,dirac_ao_num))
   dirac_ao_bi_elec_integral_tmp = 0.d0
 !$OMP DO SCHEDULE(dynamic,64)
 !DIR$ NOVECTOR
   do i8=0_8,dirac_ao_integrals_map%map_size
    n_elements = n_elements_max
    call get_cache_map(dirac_ao_integrals_map,i8,keys,values,n_elements)
    do k1=1,n_elements
     call bielec_integrals_index_reverse(kk,ii,ll,jj,keys(k1))
     do k2=1,8
      if (kk(k2)==0) then
       cycle
      endif
      i = ii(k2)
      j = jj(k2)
      k = kk(k2)
      l = ll(k2)
      integral = (dirac_SCF_density_matrix_ao(k,l)) * values(k1)
      dirac_ao_bi_elec_integral_tmp(i,j) += integral
      integral = values(k1)
      dirac_ao_bi_elec_integral_tmp(l,j) -= dirac_SCF_density_matrix_ao(k,i) * integral
     enddo
    enddo
   enddo
 !$OMP END DO NOWAIT
 !$OMP CRITICAL
   dirac_ao_bi_elec_integral += (1.d0,1.d0)* dirac_ao_bi_elec_integral_tmp
 !$OMP END CRITICAL
   deallocate(keys,values,dirac_ao_bi_elec_integral_tmp)
 !$OMP END PARALLEL
  endif
 END_PROVIDER



 
!BEGIN_PROVIDER [ complex*16, dirac_ao_bi_elec_integralnaive, (2*dirac_ao_num, 2*dirac_ao_num) ]
!implicit none
!BEGIN_DOC
!! Bi-electronic Fock matrix in dirac AO basis set
!END_DOC
!integer                        :: i,j,k,l,k1,r,s
!integer                        :: i0,j0,k0,l0
!integer                      :: p,q
!complex*16               :: D
!double precision   :: dirac_ao_bielec_integral
!dirac_ao_bi_elec_integralnaive = (0.d0,0.d0)
!do i = 1, 2*dirac_ao_num
! do j = 1, 2*dirac_ao_num
!  if (i .le. ao_num .and. j .le. ao_num) then
!!1 L_alpha L_alpha bloc
!   do k = 1, 2*dirac_ao_num 
!    do l = 1, 2*dirac_ao_num 
!     D = dirac_SCF_density_matrix_ao(k,l)
!     if (k .le. ao_num .and. l .le. ao_num) then
!      dirac_ao_bi_elec_integralnaive(i,j) += D *(dirac_ao_bielec_integral(i,j,k,l) - dirac_ao_bielec_integral(i,l,k,j))
!     elseif (k .gt. ao_num .and. k .le. (2*ao_num) .and. l .gt. ao_num .and. l .le. (2*ao_num)) then
!      dirac_ao_bi_elec_integralnaive(i,j) += D *dirac_ao_bielec_integral(i,j,k,l) 
!     elseif (k .gt. (2*ao_num) .and. k .le. (2*ao_num+small_ao_num) .and. l .gt. (2*ao_num) .and. l .le. (2*ao_num+small_ao_num)) then
!      dirac_ao_bi_elec_integralnaive(i,j) += D *dirac_ao_bielec_integral(i,j,k,l) 
!     elseif (k .gt. (2*ao_num + small_ao_num) .and. l .gt. (2*ao_num + small_ao_num)) then
!      dirac_ao_bi_elec_integralnaive(i,j) += D *dirac_ao_bielec_integral(i,j,k,l)
!     endif  
!    enddo
!   enddo
!  elseif(i .gt. ao_num .and. i .le. (2*ao_num) .and. j .gt. ao_num .and. j .le. (2*ao_num)) then
!!2 L_beta L_beta 
!    do k = 1, 2*dirac_ao_num
!    do l = 1, 2*dirac_ao_num
!     D = dirac_SCF_density_matrix_ao(k,l)
!     if (k .le. ao_num .and. l .le. ao_num) then
!      dirac_ao_bi_elec_integralnaive(i,j) += D *(dirac_ao_bielec_integral(i,j,k,l))
!     elseif (k .gt. ao_num .and. k .le. (2*ao_num) .and. l.gt. ao_num .and. l.le. (2*ao_num)) then
!      dirac_ao_bi_elec_integralnaive(i,j) += D *(dirac_ao_bielec_integral(i,j,k,l) -dirac_ao_bielec_integral(i,l,k,j)) 
!     elseif (k .gt. (2*ao_num) .and. k .le. (2*ao_num+small_ao_num) .and. l.gt. (2*ao_num) .and. l .le. (2*ao_num+small_ao_num)) then
!      dirac_ao_bi_elec_integralnaive(i,j) += D *dirac_ao_bielec_integral(i,j,k,l) 
!     elseif (k .gt. (2*ao_num+small_ao_num) .and. l.gt.(2*ao_num+small_ao_num)) then
!      dirac_ao_bi_elec_integralnaive(i,j) += D *dirac_ao_bielec_integral(i,j,k,l)
!     endif
!    enddo
!   enddo 
!  elseif(i .le. ao_num .and. j .gt. ao_num .and. j .le. 2*ao_num) then
!!3 L_alpha L_beta
!   do k = 1, 2*dirac_ao_num
!    do l = 1, 2*dirac_ao_num
!     D = dirac_SCF_density_matrix_ao(k,l)
!     if (k .gt. ao_num .and. k .le. 2*ao_num .and. l .le. ao_num) then
!      dirac_ao_bi_elec_integralnaive(i,j) += D *( -dirac_ao_bielec_integral(i,l,k,j)) 
!     endif
!    enddo
!   enddo 
!  elseif(i .gt. ao_num .and. i .le. 2*ao_num .and. j .le. ao_num) then
!!4 L_beta L_alpha 
!   do k = 1, 2*dirac_ao_num
!    do l = 1, 2*dirac_ao_num
!     D = dirac_SCF_density_matrix_ao(k,l)
!     if (k .le. ao_num .and. l .gt. ao_num .and. l .le. 2*ao_num) then
!      dirac_ao_bi_elec_integralnaive(i,j) += D *( -dirac_ao_bielec_integral(i,l,k,j)) 
!     endif
!    enddo
!   enddo 
!  elseif (i .gt. 2*ao_num .and. i .le. (2*ao_num+small_ao_num) .and. j .gt. 2*ao_num .and. j .le. (2*ao_num+small_ao_num)) then
!!5 S_alpha S_alpha 
!   do k = 1, 2*dirac_ao_num 
!    do l = 1, 2*dirac_ao_num 
!     D = dirac_SCF_density_matrix_ao(k,l)
!     if (k .le. ao_num .and. l .le. ao_num) then
!      dirac_ao_bi_elec_integralnaive(i,j) += D * (dirac_ao_bielec_integral(i,j,k,l))
!     elseif (k .gt. ao_num .and. k .le. 2*ao_num .and. l.gt. ao_num .and. l .le. 2*ao_num) then
!      dirac_ao_bi_elec_integralnaive(i,j) += D * (dirac_ao_bielec_integral(i,j,k,l)) 
!     elseif (k .gt. 2*ao_num .and. k .le. (2*ao_num+small_ao_num) .and. l.gt. 2*ao_num .and. l .le. (2*ao_num+small_ao_num)) then
!      dirac_ao_bi_elec_integralnaive(i,j) += D * (dirac_ao_bielec_integral(i,j,k,l) - dirac_ao_bielec_integral(i,l,k,j))
!     elseif (k .gt. (2*ao_num+small_ao_num) .and. l.gt. (2*ao_num+small_ao_num)) then
!      dirac_ao_bi_elec_integralnaive(i,j) += D * (dirac_ao_bielec_integral(i,j,k,l))
!     endif  
!    enddo
!   enddo
!  elseif(i .gt. (2*ao_num+small_ao_num) .and. j .gt. (2*ao_num+small_ao_num)) then
!!6 S_beta S_beta 
!    do k = 1, 2*dirac_ao_num
!    do l = 1, 2*dirac_ao_num
!     D = dirac_SCF_density_matrix_ao(k,l)
!     if (k .le. ao_num .and. l .le. ao_num) then
!      dirac_ao_bi_elec_integralnaive(i,j) += D *(dirac_ao_bielec_integral(i,j,k,l))
!     elseif (k .gt. ao_num .and. k .le. 2*ao_num .and. l.gt. ao_num .and. l.le. 2*ao_num) then
!      dirac_ao_bi_elec_integralnaive(i,j) += D *(dirac_ao_bielec_integral(i,j,k,l)) 
!     elseif (k .gt. 2*ao_num .and. k .le. (2*ao_num+small_ao_num) .and. l.gt. 2*ao_num .and. l .le. (2*ao_num+small_ao_num)) then
!      dirac_ao_bi_elec_integralnaive(i,j) += D *dirac_ao_bielec_integral(i,j,k,l) 
!     elseif (k .gt. (2*ao_num+small_ao_num) .and. l.gt.(2*ao_num+small_ao_num)) then
!      dirac_ao_bi_elec_integralnaive(i,j) += D *(dirac_ao_bielec_integral(i,j,k,l) - dirac_ao_bielec_integral(i,l,k,j))
!     endif
!    enddo
!   enddo 
!  elseif (i .gt. 2*ao_num .and. i .le. (2*ao_num+small_ao_num) .and. j .gt. (2*ao_num+small_ao_num)) then
!!7 S_alpha S_beta
!   do k = 1, 2*dirac_ao_num
!    do l = 1, 2*dirac_ao_num
!     D = dirac_SCF_density_matrix_ao(k,l)
!     if (k .gt. (2*ao_num+small_ao_num) .and. l .gt. 2*ao_num .and. l .le. (2*ao_num+small_ao_num)) then
!      dirac_ao_bi_elec_integralnaive(i,j) += D *( -dirac_ao_bielec_integral(i,l,k,j)) 
!     endif
!    enddo
!   enddo 
!  elseif (i .gt. (2*ao_num+small_ao_num) .and. j .gt. 2*ao_num .and. j .le. (2*ao_num+small_ao_num)) then
!!8 S_beta S_alpha
!   do k = 1, 2*dirac_ao_num
!    do l = 1, 2*dirac_ao_num
!     D = dirac_SCF_density_matrix_ao(k,l)
!     if (k .gt. 2*ao_num .and. j .le. (2*ao_num+small_ao_num) .and. l .gt. (2*ao_num+small_ao_num)) then
!      dirac_ao_bi_elec_integralnaive(i,j) += D *( -dirac_ao_bielec_integral(i,l,k,j)) 
!     endif
!    enddo
!   enddo 
!  elseif (i .le. ao_num .and. j .gt. 2*ao_num .and. j .le. (2*ao_num+small_ao_num)) then
!!9 L_alpha S_alpha
!   do k = 1, 2*dirac_ao_num
!    do l = 1, 2*dirac_ao_num
!     D = dirac_SCF_density_matrix_ao(k,l)
!     if (k .gt. 2*ao_num .and. k .le. (2*ao_num+small_ao_num) .and. l .le. ao_num) then
!      dirac_ao_bi_elec_integralnaive(i,j) += D *( -dirac_ao_bielec_integral(i,l,k,j)) 
!     endif
!    enddo
!   enddo 
!  elseif (i .gt. 2*ao_num .and. i .le. (2*ao_num+small_ao_num) .and. j .le. ao_num) then
!!10 S_alpha L_alpha
!   do k = 1, 2*dirac_ao_num
!    do l = 1, 2*dirac_ao_num
!     D = dirac_SCF_density_matrix_ao(k,l)
!     if (k .le. ao_num .and. l .gt. 2*ao_num .and. l .le. (2*ao_num+small_ao_num)) then
!      dirac_ao_bi_elec_integralnaive(i,j) += D *(-dirac_ao_bielec_integral(i,l,k,j))
!     endif
!    enddo
!   enddo
!  elseif (i .gt. ao_num .and. i .le. 2*ao_num .and. j .gt. (2*ao_num+small_ao_num)) then
!!11 L_beta S_beta
!   do k = 1, 2*dirac_ao_num
!    do l = 1, 2*dirac_ao_num
!     D = dirac_SCF_density_matrix_ao(k,l)
!     if (k .gt. (2*ao_num+small_ao_num) .and. l .gt. ao_num .and. l .le. 2*ao_num) then
!      dirac_ao_bi_elec_integralnaive(i,j) += D *( -dirac_ao_bielec_integral(i,l,k,j))
!     endif
!    enddo
!   enddo
!  elseif (i .gt. (2*ao_num+small_ao_num) .and. j .gt. ao_num .and. j .le. 2*ao_num) then
!!12 S_beta L_beta
!   do k = 1, 2*dirac_ao_num
!    do l = 1, 2*dirac_ao_num
!     D = dirac_SCF_density_matrix_ao(k,l)
!     if (k .gt. ao_num .and. k .le. 2*ao_num .and. l .gt. (2*ao_num+small_ao_num)) then
!      dirac_ao_bi_elec_integralnaive(i,j) += D *(-dirac_ao_bielec_integral(i,l,k,j))
!     endif
!    enddo
!   enddo
!  elseif (i .le. ao_num .and. j .gt. (2*ao_num+small_ao_num)) then
!!13 L_alpha S_beta
!   do k = 1, 2*dirac_ao_num
!    do l = 1, 2*dirac_ao_num
!     D = dirac_SCF_density_matrix_ao(k,l)
!     if (k .gt. (2*ao_num+small_ao_num) .and. l .le. ao_num) then
!      dirac_ao_bi_elec_integralnaive(i,j) += D *( -dirac_ao_bielec_integral(i,l,k,j)) 
!     endif
!    enddo
!   enddo 
!  elseif (i .gt. (2*ao_num) .and. i .le. (2*ao_num+small_ao_num) .and. j .gt. ao_num .and. j .le. (2*ao_num)) then
!!14 S_alpha L_beta
!   do k = 1, 2*dirac_ao_num
!    do l = 1, 2*dirac_ao_num
!     D = dirac_SCF_density_matrix_ao(k,l)
!     if (k .gt. ao_num .and. k .le. (2*ao_num) .and. l .gt. (2*ao_num) .and. l .le. (2*ao_num+small_ao_num)) then 
!      dirac_ao_bi_elec_integralnaive(i,j) += D *( -dirac_ao_bielec_integral(i,l,k,j)) 
!     endif
!    enddo
!   enddo 
!  elseif (i .gt. ao_num .and. i .le. 2*ao_num .and. j .gt. 2*ao_num .and. j .le. (2*ao_num+small_ao_num)) then
!!15 L_beta S_alpha
!   do k = 1, 2*dirac_ao_num
!    do l = 1, 2*dirac_ao_num
!     D = dirac_SCF_density_matrix_ao(k,l)
!     if (k .gt. 2*ao_num .and. k .le. (2*ao_num+small_ao_num) .and. l .gt. ao_num .and. l .le. 2*ao_num) then
!      dirac_ao_bi_elec_integralnaive(i,j) += D *(-dirac_ao_bielec_integral(i,l,k,j))
!     endif
!    enddo
!   enddo
! elseif (i .gt. (2*ao_num+small_ao_num) .and. j .le. ao_num) then
!!16 S_beta L_alpha
!   do k = 1, 2*dirac_ao_num
!    do l = 1, 2*dirac_ao_num
!     D = dirac_SCF_density_matrix_ao(k,l)
!     if (k .le. ao_num .and. l .gt. (2*ao_num+small_ao_num)) then
!      dirac_ao_bi_elec_integralnaive(i,j) += D *(-dirac_ao_bielec_integral(i,l,k,j))
!     endif
!    enddo
!   enddo 
!  endif
! enddo
!enddo  
!END_PROVIDER


