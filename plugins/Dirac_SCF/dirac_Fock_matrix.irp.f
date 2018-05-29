 BEGIN_PROVIDER [ complex*16, dirac_ao_mono_elec_integral,(2*dirac_ao_num,2*dirac_ao_num)]
  implicit none
  integer          :: i,j
  BEGIN_DOC
 ! array of the mono electronic hamiltonian on the AOs basis
 ! in the 4x4 component formalism with cartesian basis and 
 ! the unrestricted kinetic-balance scheme  
  END_DOC
  do i = 1, 2*(dirac_ao_num)
   do j = 1, 2*(dirac_ao_num)
    dirac_ao_mono_elec_integral(i,j) += (dirac_ao_mono_elec_nucl_integral(i,j) + dirac_ao_mono_elec_mass_integral(i,j) + dirac_ao_mono_elec_kinetic_integral(i,j) )
   enddo
  enddo
 END_PROVIDER


 BEGIN_PROVIDER [ complex*16, dirac_ao_bi_elec_integral, (2*dirac_ao_num, 2*dirac_ao_num) ]
  use map_module
  implicit none
  BEGIN_DOC
  ! Fock matrix in Dirac AO basis set
  END_DOC
  PROVIDE dirac_ao_bielec_integrals_in_map
  integer                        :: i,j,k,l,k1 
  integer                        :: r,s,p,q
  double precision               :: dirac_ao_bielec_integral, local_threshold
  complex*16                     :: integral
  complex*16, allocatable        :: dirac_ao_bi_elec_integral_tmp(:,:)
  integer(omp_lock_kind)         :: lck(dirac_ao_num)
  integer*8                      :: i8
  integer                        :: ii(8), jj(8), kk(8), ll(8), k2
  integer(cache_map_size_kind)   :: n_elements_max, n_elements
  integer(key_kind), allocatable :: keys(:)
  double precision, allocatable  :: values(:)
  dirac_ao_bi_elec_integral = (0.d0,0.d0)
 !$OMP PARALLEL DEFAULT(NONE)                                      &
 !$OMP PRIVATE(i,j,l,k1,k,integral,ii,jj,kk,ll,i8,keys,values,n_elements_max, &
 !$OMP  n_elements,dirac_ao_bi_elec_integral_tmp)&
 !$OMP SHARED(dirac_ao_num,dirac_SCF_density_matrix_ao,&
 !$OMP  dirac_ao_integrals_map, dirac_ao_bi_elec_integral) 
  call get_cache_map_n_elements_max(dirac_ao_integrals_map,n_elements_max)
  allocate(keys(n_elements_max), values(n_elements_max))
  allocate(dirac_ao_bi_elec_integral_tmp(2*dirac_ao_num,2*dirac_ao_num))
  dirac_ao_bi_elec_integral_tmp = (0.d0,0.d0)
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
     i = ii(k2) ! electron 1
     j = jj(k2) ! electron 1
     k = kk(k2) ! electron 2
     l = ll(k2) ! electron 2
     ! values(k1) = (ij|kl) <=> <ik|jl>
     integral = (dirac_SCF_density_matrix_ao(k,l)) * values(k1)
     dirac_ao_bi_elec_integral_tmp(i,j) += integral
     integral = values(k1)
     dirac_ao_bi_elec_integral_tmp(l,j) -= dirac_SCF_density_matrix_ao(k,i) * integral
    enddo
   enddo
  enddo
 !$OMP END DO NOWAIT
 !$OMP CRITICAL
  dirac_ao_bi_elec_integral += dirac_ao_bi_elec_integral_tmp
 !$OMP END CRITICAL
  deallocate(keys,values,dirac_ao_bi_elec_integral_tmp)
 !$OMP END PARALLEL
 END_PROVIDER



 ! 1 
 BEGIN_PROVIDER [ complex*16, dirac_ao_bi_elec_integral_L_alpha_L_alpha, (ao_num, ao_num) ]
  implicit none
  BEGIN_DOC
  ! L_alpha L_alpha bloc of the bi-electronic Fock matrix in dirac AO basis set
   END_DOC
  integer                        :: i,j,k,l 
  complex*16                     :: D
  double precision               :: dirac_ao_bielec_integral   
  dirac_ao_bi_elec_integral_L_alpha_L_alpha = (0.d0,0.d0)
  do i = 1, ao_num
   do j = 1, ao_num
    do k = 1, 2*dirac_ao_num
     do l = 1, 2*dirac_ao_num
      D = dirac_SCF_density_matrix_ao(k,l)     
      if (k .le. ao_num .and. l .le. ao_num) then
       dirac_ao_bi_elec_integral_L_alpha_L_alpha(i,j) += D*(dirac_ao_bielec_integral(i,j,d_L(k),d_L(l)) - dirac_ao_bielec_integral(i,d_L(l),d_L(k),j))
      elseif ((k .gt. ao_num .and. k .le. (2*ao_num) .and. l.gt. ao_num .and. l.le. (2*ao_num)) .or. &
              (k .gt. (2*ao_num) .and. k .le. (2*ao_num+small_ao_num) .and. l .gt. (2*ao_num) .and. l .le. (2*ao_num+small_ao_num)) .or. &
              (k .gt. (2*ao_num+small_ao_num) .and. l.gt.(2*ao_num+small_ao_num))) then
       dirac_ao_bi_elec_integral_L_alpha_L_alpha(i,j) += D*dirac_ao_bielec_integral(i,j,d_L(k),d_L(l))
      endif
     enddo
    enddo
   enddo
  enddo
 END_PROVIDER 

 !2
 BEGIN_PROVIDER [ complex*16, dirac_ao_bi_elec_integral_L_beta_L_beta, (ao_num, ao_num) ]
  implicit none
  BEGIN_DOC
  ! L_beta L_beta bloc of the bi-electronic Fock matrix in dirac AO basis set
   END_DOC
  integer                        :: i,j,k,l 
  complex*16                     :: D
  double precision               :: dirac_ao_bielec_integral   
  dirac_ao_bi_elec_integral_L_beta_L_beta = (0.d0,0.d0)
  do i = 1, ao_num
   do j = 1, ao_num
    do k = 1, 2*dirac_ao_num
     do l = 1, 2*dirac_ao_num
      D = dirac_SCF_density_matrix_ao(k,l)
      if (k .gt. ao_num .and. k .le. (2*ao_num) .and. l.gt. ao_num .and. l.le. (2*ao_num)) then
       dirac_ao_bi_elec_integral_L_beta_L_beta(i,j) += D*(dirac_ao_bielec_integral(i,j,d_L(k),d_L(l)) - dirac_ao_bielec_integral(i,d_L(l),d_L(k),j))
      elseif ((k .le. ao_num .and. l .le. ao_num) .or. &
              (k > (2*ao_num) .and. k <= (2*ao_num+small_ao_num) .and. l> (2*ao_num) .and. l <= (2*ao_num+small_ao_num)) .or. &
              (k .gt. (2*ao_num+small_ao_num) .and. l.gt.(2*ao_num+small_ao_num))) then
       dirac_ao_bi_elec_integral_L_beta_L_beta(i,j) += D*dirac_ao_bielec_integral(i,j,d_L(k),d_L(l))
      endif
     enddo
    enddo
   enddo
  enddo
 END_PROVIDER 

 !3
 BEGIN_PROVIDER [ complex*16, dirac_ao_bi_elec_integral_L_beta_L_alpha, (ao_num, ao_num) ]
  implicit none
  BEGIN_DOC
  ! L_beta L_alpha bloc of the bi-electronic Fock matrix in dirac AO basis set
   END_DOC
  integer                        :: i,j,k,l 
  complex*16                     :: D
  double precision               :: dirac_ao_bielec_integral   
  dirac_ao_bi_elec_integral_L_beta_L_alpha = (0.d0,0.d0)
  do i = 1, ao_num
   do j = 1, ao_num
    do k = 1, 2*dirac_ao_num
     do l = 1, 2*dirac_ao_num
      if (k .le. ao_num .and. l .gt. ao_num .and. l .le. 2*ao_num) then
       D = dirac_SCF_density_matrix_ao(k,l)
       dirac_ao_bi_elec_integral_L_beta_L_alpha(i,j) += D*(- dirac_ao_bielec_integral(i,d_L(l),d_L(k),j))
      endif
     enddo
    enddo
   enddo
  enddo
 END_PROVIDER 

 !4
 BEGIN_PROVIDER [ complex*16, dirac_ao_bi_elec_integral_L_alpha_L_beta, (ao_num, ao_num) ]
  implicit none
  BEGIN_DOC
  ! L_beta L_alpha bloc of the bi-electronic Fock matrix in dirac AO basis set
   END_DOC
  integer                        :: i,j,k,l 
  complex*16                     :: D
  double precision               :: dirac_ao_bielec_integral   
  dirac_ao_bi_elec_integral_L_alpha_L_beta = (0.d0,0.d0)
  do i = 1, ao_num
   do j = 1, ao_num
    dirac_ao_bi_elec_integral_L_alpha_L_beta(i,j) = Conjg(dirac_ao_bi_elec_integral_L_beta_L_alpha(j,i))
   !do k = 1, 2*dirac_ao_num
   ! do l = 1, 2*dirac_ao_num
   !  D = dirac_SCF_density_matrix_ao(k,l)
   !   if (k .gt. ao_num .and. k .le. 2*ao_num .and. l .le. ao_num) then
   !    dirac_ao_bi_elec_integral_L_beta_L_alpha(i,j) += D*(- dirac_ao_bielec_integral(i,d_L(l),d_L(k),j))
   !  endif
   ! enddo
   !enddo
   enddo
  enddo
 END_PROVIDER 

 !5
 BEGIN_PROVIDER [ complex*16, dirac_ao_bi_elec_integral_S_alpha_S_alpha, (small_ao_num, small_ao_num) ]
  implicit none
  BEGIN_DOC
  ! S_alpha S_alpha bloc of the bi-electronic Fock matrix in dirac AO basis set
   END_DOC
  integer                        :: i,j,k,l 
  complex*16                     :: D
  double precision               :: dirac_ao_bielec_integral   
  dirac_ao_bi_elec_integral_S_alpha_S_alpha = (0.d0,0.d0)
  do i = 1, small_ao_num
   do j = 1, small_ao_num
    do k = 1, 2*dirac_ao_num
     do l = 1, 2*dirac_ao_num
      D = dirac_SCF_density_matrix_ao(k,l)
      if (k .gt. (2*ao_num) .and. k .le. (2*ao_num+small_ao_num) .and. l .gt. (2*ao_num) .and. l .le. (2*ao_num+small_ao_num)) then
       dirac_ao_bi_elec_integral_S_alpha_S_alpha(i,j) += D*(dirac_ao_bielec_integral(i,j,d_L(k),d_L(l)) - dirac_ao_bielec_integral(i,d_L(l),d_L(k),j))
      elseif ((k .le. ao_num .and. l .le. ao_num) .or. &
              (k .gt. ao_num .and. k .le. (2*ao_num) .and. l.gt. ao_num .and. l.le. (2*ao_num)) .or. &
              (k .gt. (2*ao_num+small_ao_num) .and. l.gt.(2*ao_num+small_ao_num))) then
       dirac_ao_bi_elec_integral_S_alpha_S_alpha(i,j) += D*dirac_ao_bielec_integral(i,j,d_L(k),d_L(l))
      endif
     enddo
    enddo
   enddo
  enddo
 END_PROVIDER 

 !6
 BEGIN_PROVIDER [ complex*16, dirac_ao_bi_elec_integral_S_beta_S_beta, (small_ao_num, small_ao_num) ]
  implicit none
  BEGIN_DOC
  ! S_beta S_beta bloc of the bi-electronic Fock matrix in dirac AO basis set
   END_DOC
  integer                        :: i,j,k,l 
  complex*16                     :: D
  double precision               :: dirac_ao_bielec_integral   
  dirac_ao_bi_elec_integral_S_beta_S_beta = (0.d0,0.d0)
  do i = 1, small_ao_num
   do j = 1, small_ao_num
    do k = 1, 2*dirac_ao_num
     do l = 1, 2*dirac_ao_num
      D = dirac_SCF_density_matrix_ao(k,l)
      if (k .gt. (2*ao_num+small_ao_num) .and. l .gt. (2*ao_num+small_ao_num)) then
       dirac_ao_bi_elec_integral_S_beta_S_beta(i,j) += D*(dirac_ao_bielec_integral(i,j,d_L(k),d_L(l)) - dirac_ao_bielec_integral(i,d_L(l),d_L(k),j))
      elseif ((k .le. ao_num .and. l .le. ao_num) .or. &
              (k .gt. ao_num .and. k .le. (2*ao_num) .and. l.gt. ao_num .and. l.le. (2*ao_num)) .or. &
              (k .gt. (2*ao_num) .and. k .le. (2*ao_num+small_ao_num) .and. l .gt. (2*ao_num) .and. l .le. (2*ao_num+small_ao_num))) then
       dirac_ao_bi_elec_integral_S_beta_S_beta(i,j) += D*dirac_ao_bielec_integral(i,j,d_L(k),d_L(l))
      endif
     enddo
    enddo
   enddo
  enddo
 END_PROVIDER 

 !7
 BEGIN_PROVIDER [ complex*16, dirac_ao_bi_elec_integral_S_beta_S_alpha, (small_ao_num, small_ao_num) ]
  implicit none
  BEGIN_DOC
  ! S_beta S_alpha bloc of the bi-electronic Fock matrix in dirac AO basis set
   END_DOC
  integer                        :: i,j,k,l 
  complex*16                     :: D
  double precision               :: dirac_ao_bielec_integral   
  dirac_ao_bi_elec_integral_S_beta_S_alpha = (0.d0,0.d0)
  do i = 1, small_ao_num
   do j = 1, small_ao_num
    do k = 1, 2*dirac_ao_num
     do l = 1, 2*dirac_ao_num
      if (k .gt. (2*ao_num) .and. k .le. (2*ao_num+small_ao_num) .and. l .gt. (2*ao_num+small_ao_num)) then
       D = dirac_SCF_density_matrix_ao(k,l)
       dirac_ao_bi_elec_integral_S_beta_S_alpha(i,j) += D*(- dirac_ao_bielec_integral(i,d_L(l),d_L(k),j))
      endif
     enddo
    enddo
   enddo
  enddo
 END_PROVIDER 

 !8
 BEGIN_PROVIDER [ complex*16, dirac_ao_bi_elec_integral_S_alpha_S_beta, (small_ao_num, small_ao_num) ]
  implicit none
  BEGIN_DOC
  ! S_alpha S_beta bloc of the bi-electronic Fock matrix in dirac AO basis set
   END_DOC
  integer                        :: i,j,k,l 
  complex*16                     :: D
  double precision               :: dirac_ao_bielec_integral   
  dirac_ao_bi_elec_integral_S_alpha_S_beta = (0.d0,0.d0)
  do i = 1, ao_num
   do j = 1, ao_num
    dirac_ao_bi_elec_integral_S_alpha_S_beta(i,j) = Conjg(dirac_ao_bi_elec_integral_S_beta_S_alpha(j,i))
   !do k = 1, 2*dirac_ao_num
   ! do l = 1, 2*dirac_ao_num
   !  D = dirac_SCF_density_matrix_ao(k,l)
   !   if (k .gt. (2*ao_num+small_ao_num) .and. l .gt. (2*ao_num) .and. l .le. (2*ao_num+small_ao_num)) then
   !    dirac_ao_bi_elec_integral_S_beta_S_alpha(i,j) += D*(- dirac_ao_bielec_integral(i,d_L(l),d_L(k),j))
   !  endif
   ! enddo
   !enddo
   enddo
  enddo
 END_PROVIDER 

 !9
 BEGIN_PROVIDER [ complex*16, dirac_ao_bi_elec_integral_S_alpha_L_alpha, (small_ao_num, ao_num) ]
  implicit none
  BEGIN_DOC
  ! S_alpha L_alpha bloc of the bi-electronic Fock matrix in dirac AO basis set
   END_DOC
  integer                        :: i,j,k,l 
  complex*16                     :: D
  double precision               :: dirac_ao_bielec_integral   
  dirac_ao_bi_elec_integral_S_alpha_L_alpha = (0.d0,0.d0)
  do i = 1, small_ao_num
   do j = 1, ao_num
    do k = 1, 2*dirac_ao_num
     do l = 1, 2*dirac_ao_num
      if (k.le. ao_num .and. l .gt. 2*ao_num .and. l .le. (2*ao_num+small_ao_num)) then
       D = dirac_SCF_density_matrix_ao(k,l)
       dirac_ao_bi_elec_integral_S_alpha_L_alpha(i,j) += D*(- dirac_ao_bielec_integral(i,d_L(l),d_L(k),j))
      endif
     enddo
    enddo
   enddo
  enddo
 END_PROVIDER

 !10
 BEGIN_PROVIDER [ complex*16, dirac_ao_bi_elec_integral_L_alpha_S_alpha, (ao_num, ao_num) ]
  implicit none
  BEGIN_DOC
  ! L_alpha S_alpha bloc of the bi-electronic Fock matrix in dirac AO basis set
   END_DOC
  integer                        :: i,j,k,l 
  complex*16                     :: D
  double precision               :: dirac_ao_bielec_integral   
  dirac_ao_bi_elec_integral_L_alpha_S_alpha = (0.d0,0.d0)
  do i = 1, ao_num
   do j = 1, small_ao_num
    dirac_ao_bi_elec_integral_L_alpha_S_alpha(i,j) = Conjg(dirac_ao_bi_elec_integral_S_alpha_L_alpha(j,i))
   !do k = 1, 2*dirac_ao_num
   ! do l = 1, 2*dirac_ao_num
   !  if (k .gt. 2*ao_num .and. k .le. (2*ao_num+small_ao_num) .and. l .le. ao_num) then
   !   D = dirac_SCF_density_matrix_ao(k,l)
   !   dirac_ao_bi_elec_integral_S_alpha_L_alpha(i,j) += D*(- dirac_ao_bielec_integral(i,d_L(l),d_L(k),j))
   !  endif
   ! enddo
   !enddo
   enddo
  enddo
 END_PROVIDER
 
 !11
 BEGIN_PROVIDER [ complex*16, dirac_ao_bi_elec_integral_S_beta_L_beta, (small_ao_num, ao_num) ]
  implicit none
  BEGIN_DOC
  ! S_beta L_beta bloc of the bi-electronic Fock matrix in dirac AO basis set
   END_DOC
  integer                        :: i,j,k,l 
  complex*16                     :: D
  double precision               :: dirac_ao_bielec_integral   
  dirac_ao_bi_elec_integral_S_beta_L_beta = (0.d0,0.d0)
  do i = 1, small_ao_num
   do j = 1, ao_num
    do k = 1, 2*dirac_ao_num
     do l = 1, 2*dirac_ao_num
      if (k .gt. ao_num .and. k .le. 2*ao_num .and. l .gt. (2*ao_num+small_ao_num)) then
       D = dirac_SCF_density_matrix_ao(k,l)
       dirac_ao_bi_elec_integral_S_alpha_L_alpha(i,j) += D*(- dirac_ao_bielec_integral(i,d_L(l),d_L(k),j))
      endif
     enddo
    enddo
   enddo
  enddo
 END_PROVIDER

 !12
 BEGIN_PROVIDER [ complex*16, dirac_ao_bi_elec_integral_L_beta_S_beta, (ao_num, small_ao_num) ]
  implicit none
  BEGIN_DOC
  ! L_beta S_beta bloc of the bi-electronic Fock matrix in dirac AO basis set
   END_DOC
  integer                        :: i,j,k,l 
  complex*16                     :: D
  double precision               :: dirac_ao_bielec_integral   
  dirac_ao_bi_elec_integral_L_beta_S_beta = (0.d0,0.d0)
  do i = 1, ao_num
   do j = 1, small_ao_num
    dirac_ao_bi_elec_integral_L_beta_S_beta(i,j) = Conjg(dirac_ao_bi_elec_integral_S_beta_L_beta(j,i))
   !do k = 1, 2*dirac_ao_num
   ! do l = 1, 2*dirac_ao_num
   !  if (k .gt. (2*ao_num+small_ao_num) .and. l .gt. ao_num .and. l .le. 2*ao_num) then
   !   D = dirac_SCF_density_matrix_ao(k,l)
   !   dirac_ao_bi_elec_integral_S_alpha_L_alpha(i,j) += D*(- dirac_ao_bielec_integral(i,d_L(l),d_L(k),j))
   !  endif
   ! enddo
   !enddo
   enddo
  enddo
 END_PROVIDER

 !13
 BEGIN_PROVIDER [ complex*16, dirac_ao_bi_elec_integral_S_alpha_L_beta, (small_ao_num, ao_num) ]
  implicit none
  BEGIN_DOC
  ! S_alpha L_beta bloc of the bi-electronic Fock matrix in dirac AO basis set
   END_DOC
  integer                        :: i,j,k,l 
  complex*16                     :: D
  double precision               :: dirac_ao_bielec_integral   
  dirac_ao_bi_elec_integral_S_alpha_L_beta = (0.d0,0.d0)
  do i = 1, small_ao_num
   do j = 1, ao_num
    do k = 1, 2*dirac_ao_num
     do l = 1, 2*dirac_ao_num
      if (k .gt. ao_num .and. k .le. 2*ao_num .and. l .gt. (2*ao_num) .and. l .le. (2*ao_num+small_ao_num)) then
       D = dirac_SCF_density_matrix_ao(k,l)
       dirac_ao_bi_elec_integral_S_alpha_L_alpha(i,j) += D*(- dirac_ao_bielec_integral(i,d_L(l),d_L(k),j))
      endif
     enddo
    enddo
   enddo
  enddo
 END_PROVIDER

 !14
 BEGIN_PROVIDER [ complex*16, dirac_ao_bi_elec_integral_L_beta_S_alpha, (ao_num, small_ao_num) ]
  implicit none
  BEGIN_DOC
  ! L_beta S_alpha bloc of the bi-electronic Fock matrix in dirac AO basis set
   END_DOC
  integer                        :: i,j,k,l 
  complex*16                     :: D
  double precision               :: dirac_ao_bielec_integral   
  dirac_ao_bi_elec_integral_L_beta_S_alpha = (0.d0,0.d0)
  do i = 1, ao_num
   do j = 1, small_ao_num
    dirac_ao_bi_elec_integral_L_beta_S_alpha(i,j) = Conjg(dirac_ao_bi_elec_integral_S_alpha_L_beta(j,i))
   !do k = 1, 2*dirac_ao_num
   ! do l = 1, 2*dirac_ao_num
   !  if (k .gt. (2*ao_num) .and. k .le. (2*ao_num+small_ao_num) .and. l .gt. ao_num .and. l .le. 2*ao_num) then
   !   D = dirac_SCF_density_matrix_ao(k,l)
   !   dirac_ao_bi_elec_integral_S_alpha_L_alpha(i,j) += D*(- dirac_ao_bielec_integral(i,d_L(l),d_L(k),j))
   !  endif
   ! enddo
   !enddo
   enddo
  enddo
 END_PROVIDER

 !15
 BEGIN_PROVIDER [ complex*16, dirac_ao_bi_elec_integral_S_beta_L_alpha, (small_ao_num, ao_num) ]
  implicit none
  BEGIN_DOC
  ! S_beta L_alpha bloc of the bi-electronic Fock matrix in dirac AO basis set
   END_DOC
  integer                        :: i,j,k,l 
  complex*16                     :: D
  double precision               :: dirac_ao_bielec_integral   
  dirac_ao_bi_elec_integral_S_beta_L_alpha = (0.d0,0.d0)
  do i = 1, small_ao_num
   do j = 1, ao_num
    do k = 1, 2*dirac_ao_num
     do l = 1, 2*dirac_ao_num
      if (k .le. ao_num .and. l .gt. (2*ao_num+small_ao_num)) then
       D = dirac_SCF_density_matrix_ao(k,l)
       dirac_ao_bi_elec_integral_S_alpha_L_alpha(i,j) += D*(- dirac_ao_bielec_integral(i,d_L(l),d_L(k),j))
      endif
     enddo
    enddo
   enddo
  enddo
 END_PROVIDER

 !16
 BEGIN_PROVIDER [ complex*16, dirac_ao_bi_elec_integral_L_alpha_S_beta, (ao_num, small_ao_num) ]
  implicit none
  BEGIN_DOC
  ! L_alpha S_beta bloc of the bi-electronic Fock matrix in dirac AO basis set
   END_DOC
  integer                        :: i,j,k,l 
  complex*16                     :: D
  double precision               :: dirac_ao_bielec_integral   
  dirac_ao_bi_elec_integral_L_alpha_S_beta = (0.d0,0.d0)
  do i = 1, ao_num
   do j = 1, small_ao_num
    dirac_ao_bi_elec_integral_L_alpha_S_beta(i,j) = Conjg(dirac_ao_bi_elec_integral_S_beta_L_alpha(j,i))
   !do k = 1, 2*dirac_ao_num
   ! do l = 1, 2*dirac_ao_num
   !  if (k .gt. (2*ao_num+small_ao_num) .and. l .le. ao_num) then
   !   D = dirac_SCF_density_matrix_ao(k,l)
   !   dirac_ao_bi_elec_integral_S_alpha_L_alpha(i,j) += D*(- dirac_ao_bielec_integral(i,d_L(l),d_L(k),j))
   !  endif
   ! enddo
   !enddo
   enddo
  enddo
 END_PROVIDER



