 double precision function dirac_ao_bielec_integral(i,j,k,l)
  implicit none
  BEGIN_DOC
  !  integral of the dirac AO basis <ik|jl> or (ij|kl)
  !     i(r1) j(r1) 1/r12 k(r2) l(r2)
  END_DOC
  integer,intent(in)             :: i,j,k,l
  integer                        :: p,q,r,s
  double precision               :: I_center(3),J_center(3),K_center(3),L_center(3)
  integer                        :: num_i,num_j,num_k,num_l,dim1,I_power(3),J_power(3),K_power(3),L_power(3)
  double precision               :: integral
  include 'Utils/constants.include.F'
  double precision               :: P_new(0:max_dim,3),P_center(3),fact_p,pp
  double precision               :: Q_new(0:max_dim,3),Q_center(3),fact_q,qq
  integer                        :: iorder_p(3), iorder_q(3)
  double precision               :: ao_bielec_integral_schwartz_accel
  dim1 = n_pt_max_integrals
  num_i = dirac_ao_nucl(i)
  num_j = dirac_ao_nucl(j)
  num_k = dirac_ao_nucl(k)
  num_l = dirac_ao_nucl(l)
  dirac_ao_bielec_integral = 0.d0
  if (num_i /= num_j .or. num_k /= num_l .or. num_j /= num_k)then
   do p = 1, 3
    I_power(p) = dirac_ao_power(i,p)
    J_power(p) = dirac_ao_power(j,p)
    K_power(p) = dirac_ao_power(k,p)
    L_power(p) = dirac_ao_power(l,p)
    I_center(p) = nucl_coord(num_i,p)
    J_center(p) = nucl_coord(num_j,p)
    K_center(p) = nucl_coord(num_k,p)
    L_center(p) = nucl_coord(num_l,p)
   enddo
  double precision               :: coef1, coef2, coef3, coef4
  double precision               :: p_inv,q_inv
  double precision               :: general_primitive_integral
   coef1 = dirac_ao_coef_normalized(i)
   coef2 = coef1*dirac_ao_coef_normalized(j)
   call give_explicit_poly_and_gaussian(P_new,P_center,pp,fact_p,iorder_p,&
        dirac_ao_expo(i),dirac_ao_expo(j),                 &
        I_power,J_power,I_center,J_center,dim1)
   p_inv = 1.d0/pp
   coef3 = coef2*dirac_ao_coef_normalized(k)
   coef4 = coef3*dirac_ao_coef_normalized(l)
   call give_explicit_poly_and_gaussian(Q_new,Q_center,qq,fact_q,iorder_q,&
        dirac_ao_expo(k),dirac_ao_expo(l),             &
        K_power,L_power,K_center,L_center,dim1)
   q_inv = 1.d0/qq
   integral = general_primitive_integral(dim1,              &
              P_new,P_center,fact_p,pp,p_inv,iorder_p,             &
              Q_new,Q_center,fact_q,qq,q_inv,iorder_q)
   dirac_ao_bielec_integral +=+  coef4 * integral
  else
   do p = 1, 3
    I_power(p) = dirac_ao_power(i,p)
    J_power(p) = dirac_ao_power(j,p)
    K_power(p) = dirac_ao_power(k,p)
    L_power(p) = dirac_ao_power(l,p)
   enddo
  double  precision              :: ERI
   coef1 = dirac_ao_coef_normalized(i)
   coef2 = coef1*dirac_ao_coef_normalized(j)
   coef3 = coef2*dirac_ao_coef_normalized(k)
   coef4 = coef3*dirac_ao_coef_normalized(l)
   integral = ERI(                                          &
              dirac_ao_expo(i),dirac_ao_expo(j),dirac_ao_expo(k),dirac_ao_expo(l),&
              I_power(1),J_power(1),K_power(1),L_power(1),         &
              I_power(2),J_power(2),K_power(2),L_power(2),         &
              I_power(3),J_power(3),K_power(3),L_power(3))
   dirac_ao_bielec_integral += coef4 * integral
  endif
 end

 BEGIN_PROVIDER [ double precision, dirac_ao_bielec_integral_schwartz,(2*dirac_ao_num,2*dirac_ao_num)  ]
  implicit none
  BEGIN_DOC
  !  Needed to compute Schwartz inequalities
  END_DOC
  integer                        :: i,k
  double precision               :: dirac_ao_bielec_integral,cpu_1,cpu_2, wall_1, wall_2
  dirac_ao_bielec_integral_schwartz(1,1) = dirac_ao_bielec_integral(1,1,1,1)
  do i=1,2*dirac_ao_num
    do k=1,i
      dirac_ao_bielec_integral_schwartz(i,k) = dsqrt(dirac_ao_bielec_integral(i,k,i,k))
      dirac_ao_bielec_integral_schwartz(k,i) = dirac_ao_bielec_integral_schwartz(i,k)
    enddo
  enddo
 END_PROVIDER


 BEGIN_PROVIDER [ complex*16, dirac_ao_bi_elec_integral, (2*dirac_ao_num, 2*dirac_ao_num) ]
 implicit none
 BEGIN_DOC
 ! Bi-electronic Fock matrix in dirac AO basis set
 END_DOC
 integer                        :: i,j,k,l,k1,r,s
 integer                        :: i0,j0,k0,l0
 integer                      :: p,q
 complex*16               :: D
 complex*16               :: dirac_ao_bielec_integral, local_threshold
 complex*16, allocatable  :: dirac_ao_bi_elec_integral_tmp(:,:)
 dirac_ao_bi_elec_integral = (0.d0,0.d0)
 do i = 1, 2*dirac_ao_num
  do j = 1, 2*dirac_ao_num
   if (i .le. ao_num .and. j .le. ao_num) then
 !1 L_alpha L_alpha bloc
    do k = 1, 2*dirac_ao_num 
     do l = 1, 2*dirac_ao_num 
      D = dirac_SCF_density_matrix_ao(k,l)
      if (k .le. ao_num .and. l .le. ao_num) then
       dirac_ao_bi_elec_integral(i,j) += D * (dirac_ao_bielec_integral(i,j,k,l) - dirac_ao_bielec_integral(i,l,k,j))
      elseif (k .gt. ao_num .and. k .le. 2*ao_num .and. l.gt. ao_num .and. l .le. 2*ao_num) then
       dirac_ao_bi_elec_integral(i,j) += D *dirac_ao_bielec_integral(i,j,k,l) 
      elseif (k .gt. 2*ao_num .and. k .le. (2*ao_num+small_ao_num) .and. l.gt. 2*ao_num .and. l .le. (2*ao_num+small_ao_num)) then
       dirac_ao_bi_elec_integral(i,j) += D *dirac_ao_bielec_integral(i,j,k,l) 
      elseif (k .gt. (2*ao_num+small_ao_num) .and. l.gt. (2*ao_num+small_ao_num))
       dirac_ao_bi_elec_integral(i,j) += D *dirac_ao_bielec_integral(i,j,k,l)
      endif  
     enddo
    enddo
   elseif((i .gt. ao_num .and. i .le. 2*ao_num .and. j .gt. ao_num .and. j .le. 2*ao_num) then
 !2 L_beta L_beta 
     do k = 1, 2*dirac_ao_num
     do l = 1, 2*dirac_ao_num
      D = dirac_SCF_density_matrix_ao(k,l)
      if (k .le. ao_num .and. l .le. ao_num) then
       dirac_ao_bi_elec_integral(i,j) += D *(dirac_ao_bielec_integral(i,j,k,l))
      elseif (k .gt. ao_num .and. k .le. 2*ao_num .and. l.gt. ao_num .and. l.le. 2*ao_num) then
       dirac_ao_bi_elec_integral(i,j) += D *(dirac_ao_bielec_integral(i,j,k,l) -dirac_ao_bielec_integral(i,l,k,j)) 
      elseif (k .gt. 2*ao_num .and. k .le. (2*ao_num+small_ao_num) .and. l.gt. 2*ao_num .and. l .le. (2*ao_num+small_ao_num)) then
       dirac_ao_bi_elec_integral(i,j) += D *dirac_ao_bielec_integral(i,j,k,l) 
      elseif (k .gt. (2*ao_num+small_ao_num) .and. l.gt.(2*ao_num+small_ao_num))
       dirac_ao_bi_elec_integral(i,j) += D *dirac_ao_bielec_integral(i,j,k,l)
      endif
     enddo
    enddo 
   elseif((i .le. ao_num .and. j .gt. ao_num .and. j .le. 2*ao_num) then
 !3 L_alpha L_beta
    do k = 1, 2*dirac_ao_num
     do l = 1, 2*dirac_ao_num
      D = dirac_SCF_density_matrix_ao(k,l)
      if (k .gt. ao_num .and. k .le. 2*ao_num .and. l .le. ao_num) then
       dirac_ao_bi_elec_integral(i,j) += D *( -dirac_ao_bielec_integral(i,l,k,j)) 
      endif
     enddo
    enddo 
   elseif((i .gt. ao_num .and. i .le. 2*ao_num .and. j .le. ao_num) then
 !4 L_beta L_alpha 
    do k = 1, 2*dirac_ao_num
     do l = 1, 2*dirac_ao_num
      D = dirac_SCF_density_matrix_ao(k,l)
      if (k .le. ao_num .and. l .gt. ao_num .and. l .le. 2*ao_num) then
       dirac_ao_bi_elec_integral(i,j) += D *( -dirac_ao_bielec_integral(i,l,k,j)) 
      endif
     enddo
    enddo 
   elseif (i .gt. 2*ao_num .and. i .le. (2*ao_num+small_ao_num) .and. j .gt. 2*ao_num .and. j .le. (2*ao_num+small_ao_num)) then
 !5 S_alpha S_alpha 
    do k = 1, 2*dirac_ao_num 
     do l = 1, 2*dirac_ao_num 
      D = dirac_SCF_density_matrix_ao(k,l)
      if (k .le. ao_num .and. l .le. ao_num) then
       dirac_ao_bi_elec_integral(i,j) += D * (dirac_ao_bielec_integral(i,j,k,l))
      elseif (k .gt. ao_num .and. k .le. 2*ao_num .and. l.gt. ao_num .and. l .le. 2*ao_num) then
       dirac_ao_bi_elec_integral(i,j) += D * (dirac_ao_bielec_integral(i,j,k,l)) 
      elseif (k .gt. 2*ao_num .and. k .le. (2*ao_num+small_ao_num) .and. l.gt. 2*ao_num .and. l .le. (2*ao_num+small_ao_num)) then
       dirac_ao_bi_elec_integral(i,j) += D * (dirac_ao_bielec_integral(i,j,k,l) - dirac_ao_bielec_integral(i,l,k,j))
      elseif (k .gt. (2*ao_num+small_ao_num) .and. l.gt. (2*ao_num+small_ao_num))
       dirac_ao_bi_elec_integral(i,j) += D * (dirac_ao_bielec_integral(i,j,k,l))
      endif  
     enddo
    enddo
   elseif((i .gt. (2*ao_num+small_ao_num) .and. j .gt. (2*ao_num+small_ao_num)) then
 !6 S_beta S_beta 
     do k = 1, 2*dirac_ao_num
     do l = 1, 2*dirac_ao_num
      D = dirac_SCF_density_matrix_ao(k,l)
      if (k .le. ao_num .and. l .le. ao_num) then
       dirac_ao_bi_elec_integral(i,j) += D *(dirac_ao_bielec_integral(i,j,k,l))
      elseif (k .gt. ao_num .and. k .le. 2*ao_num .and. l.gt. ao_num .and. l.le. 2*ao_num) then
       dirac_ao_bi_elec_integral(i,j) += D *(dirac_ao_bielec_integral(i,j,k,l)) 
      elseif (k .gt. 2*ao_num .and. k .le. (2*ao_num+small_ao_num) .and. l.gt. 2*ao_num .and. l .le. (2*ao_num+small_ao_num)) then
       dirac_ao_bi_elec_integral(i,j) += D *dirac_ao_bielec_integral(i,j,k,l) 
      elseif (k .gt. (2*ao_num+small_ao_num) .and. l.gt.(2*ao_num+small_ao_num))
       dirac_ao_bi_elec_integral(i,j) += D *(dirac_ao_bielec_integral(i,j,k,l) - dirac_ao_bielec_integral(i,l,k,j))
      endif
     enddo
    enddo 
   elseif (i .gt. 2*ao_num .and. i .le. (2*ao_num+small_ao_num) .and. j .gt. (2*ao_num+small_ao_num)) then
 !7 S_alpha S_beta
    do k = 1, 2*dirac_ao_num
     do l = 1, 2*dirac_ao_num
      D = dirac_SCF_density_matrix_ao(k,l)
      if (k .gt. (2*ao_num+small_ao_num) .and. l .gt. 2*ao_num .and. l .le. (2*ao_num+small_ao_num)) then
       dirac_ao_bi_elec_integral(i,j) += D *( -dirac_ao_bielec_integral(i,l,k,j)) 
      endif
     enddo
    enddo 
   elseif (i .gt. (2*ao_num+small_ao_num) .and. j .gt. 2*ao_num .and. j .le. (2*ao_num+small_ao_num)) then
 !8 S_beta S_alpha
    do k = 1, 2*dirac_ao_num
     do l = 1, 2*dirac_ao_num
      D = dirac_SCF_density_matrix_ao(k,l)
      if (k .gt. 2*ao_num .and. j .le. (2*ao_num+small_ao_num) .and. l .gt. (2*ao_num+small_ao_num)) then
       dirac_ao_bi_elec_integral(i,j) += D *( -dirac_ao_bielec_integral(i,l,k,j)) 
      endif
     enddo
    enddo 
   elseif (i .le. ao_num .and. j .gt. 2*ao_num .and. j .le. (2*ao_num+small_ao_num)) then
 !9 L_alpha S_alpha
    do k = 1, 2*dirac_ao_num
     do l = 1, 2*dirac_ao_num
      D = dirac_SCF_density_matrix_ao(k,l)
      if (k .gt. 2*ao_num .and. k .le. (2*ao_num+small_ao_num) .and. l .le. ao_num) then
       dirac_ao_bi_elec_integral(i,j) += D *( -dirac_ao_bielec_integral(i,l,k,j)) 
      endif
     enddo
    enddo 
   elseif (i .gt. 2*ao_num .and. i .le. (2*ao_num+small_ao_num) .and. j .le. ao_num) then
 !10 S_alpha L_alpha
    do k = 1, 2*dirac_ao_num
     do l = 1, 2*dirac_ao_num
      D = dirac_SCF_density_matrix_ao(k,l)
      if (k .le. ao_num .and. l .gt. 2*ao_num .and. l .le. (2*ao_num+small_ao_num)) then
       dirac_ao_bi_elec_integral(i,j) += D *(-dirac_ao_bielec_integral(i,l,k,j))
      endif
     enddo
    enddo
   elseif (i .gt. ao_num .and. i .le. 2*ao_num .and. j .gt. (2*ao_num+small_ao_num)) then
 !11 L_beta S_beta
    do k = 1, 2*dirac_ao_num
     do l = 1, 2*dirac_ao_num
      D = dirac_SCF_density_matrix_ao(k,l)
      if (k .gt. (2*ao_num+small_ao_num) .and. l .gt. ao_num and l .le. 2*ao_num) then
       dirac_ao_bi_elec_integral(i,j) += D *( -dirac_ao_bielec_integral(i,l,k,j))
      endif
     enddo
    enddo
   elseif (i .gt. (2*ao_num+small_ao_num) .and. j .gt. ao_num .and. j .le. 2*ao_num) then
 !12 S_beta L_beta
    do k = 1, 2*dirac_ao_num
     do l = 1, 2*dirac_ao_num
      D = dirac_SCF_density_matrix_ao(k,l)
      if (k .gt. ao_num .and. k .le. 2*ao_num .and. l .gt. (2*ao_num+small_ao_num)) then
       dirac_ao_bi_elec_integral(i,j) += D *(-dirac_ao_bielec_integral(i,l,k,j))
      endif
     enddo
    enddo
   elseif (i .le. ao_num .and. j .gt. (2*ao_num+small_ao_num)) then
 !13 L_alpha S_beta
    do k = 1, 2*dirac_ao_num
     do l = 1, 2*dirac_ao_num
      D = dirac_SCF_density_matrix_ao(k,l)
      if (k .gt. (2*ao_num+small_ao_num) .and. l .le. ao_num) then
       dirac_ao_bi_elec_integral(i,j) += D *( -dirac_ao_bielec_integral(i,l,k,j)) 
      endif
     enddo
    enddo 
   elseif (i .gt. (2*ao_num) .and. i .le. (2*ao_num+small_no_num) .and. j .gt. ao_num .and. j .le. 2*ao_num) then
 !14 S_alpha L_beta
    do k = 1, 2*dirac_ao_num
     do l = 1, 2*dirac_ao_num
      D = dirac_SCF_density_matrix_ao(k,l)
      if (k .gt. ao_num .and. k .le. 2*ao_num .and. l .gt. 2*ao_num l .le. (2*ao_num+small_ao_num)) then 
       dirac_ao_bi_elec_integral(i,j) += D *( -dirac_ao_bielec_integral(i,l,k,j)) 
      endif
     enddo
    enddo 
   elseif (i .gt. ao_num .and. i .le. 2*ao_num .and. j .gt. 2*ao_num .and. j .le. (2*ao_num+small_ao_num)) then
 !15 L_beta S_alpha
    do k = 1, 2*dirac_ao_num
     do l = 1, 2*dirac_ao_num
      D = dirac_SCF_density_matrix_ao(k,l)
      if (k .gt. 2*ao_num .and. k .le. (2*ao_num+small_ao_num) .and. l .gt. ao_num l .le. 2*ao_num) then
       dirac_ao_bi_elec_integral(i,j) += D *(-dirac_ao_bielec_integral(i,l,k,j))
      endif
     enddo
    enddo
  elseif (i .gt. (2*ao_num+small_ao_num) .and. j .le. ao_num) then
 !16 S_beta L_alpha
    do k = 1, 2*dirac_ao_num
     do l = 1, 2*dirac_ao_num
      D = dirac_SCF_density_matrix_ao(k,l)
      if (k .le. ao_num .and. l .gt. (2*ao_num+small_ao_num)) then
       dirac_ao_bi_elec_integral(i,j) += D *(-dirac_ao_bielec_integral(i,l,k,j))
      endif
     enddo
    enddo 
   endif
  enddo
 enddo  
 END_PROVIDER

!BEGIN_PROVIDER [ complex*16, dirac_ao_bi_elec_integral, (2*dirac_ao_num, 2*dirac_ao_num) ]
!use map_module
!implicit none
!BEGIN_DOC
!! Alpha Fock matrix in AO basis set
!END_DOC
!integer                        :: i,j,k,l,k1,r,s
!integer                        :: i0,j0,k0,l0
!integer                      :: p,q
!complex*16               :: dirac_integral, c0, c1, c2
!complex*16               :: dirac_ao_bielec_integral, local_threshold
!complex*16, allocatable  :: dirac_ao_bi_elec_integral_tmp(:,:)

!dirac_ao_bi_elec_integral = (0.d0,0.d0)
!if (do_direct_integrals) then
!  allocate(keys(1), values(1))
!  allocate(dirac_ao_bi_elec_integral(2*dirac_ao_num,2*dirac_ao_num))
!  dirac_ao_bi_elec_integral = (0.d0,0.d0)
!  q = 2*dirac_ao_num * 2*dirac_ao_num * 2*dirac_ao_num * 2*dirac_ao_num
!  do p=1_8,q
!          call bielec_integrals_index_reverse(kk,ii,ll,jj,p)
!          if ( (kk(1)>ao_num)∨ &
!               (ii(1)>ao_num)∨ &
!               (jj(1)>ao_num)∨ &
!               (ll(1)>ao_num) ) then
!               cycle
!          endif
!          k = kk(1)
!          i = ii(1)
!          l = ll(1)
!          j = jj(1)

!          if (ao_overlap_abs(k,l)×ao_overlap_abs(i,j)  &
!             < ao_integrals_threshold) then
!            cycle
!          endif
!          local_threshold = ao_bielec_integral_schwartz(k,l)×ao_bielec_integral_schwartz(i,j)
!          if (local_threshold < ao_integrals_threshold) then
!            cycle
!          endif
!          i0 = i
!          j0 = j
!          k0 = k
!          l0 = l
!          values(1) = 0.d0
!          local_threshold = ao_integrals_threshold/local_threshold
!          do k2=1,8
!            if (kk(k2)≡0) then
!              cycle
!            endif
!            i = ii(k2)
!            j = jj(k2)
!            k = kk(k2)
!            l = ll(k2)
!            c0 = SCF_density_matrix_ao_alpha(k,l)+SCF_density_matrix_ao_beta(k,l)
!            c1 = SCF_density_matrix_ao_alpha(k,i)
!            c2 = SCF_density_matrix_ao_beta(k,i)
!            if ( dabs(c0)+dabs(c1)+dabs(c2) < local_threshold) then
!              cycle
!            endif
!            if (values(1) == 0.d0) then
!              values(1) = ao_bielec_integral(k0,l0,i0,j0)
!            endif
!            integral = c0 * values(1)
!            ao_bi_elec_integral_alpha_tmp(i,j) += integral
!            ao_bi_elec_integral_beta_tmp (i,j) += integral
!            integral = values(1)
!            ao_bi_elec_integral_alpha_tmp(i,j) += integral
!            ao_bi_elec_integral_beta_tmp (i,j) += integral
!            integral = values(1)
!            ao_bi_elec_integral_alpha_tmp(l,j) -= c1 * integral
!            ao_bi_elec_integral_beta_tmp (l,j) -= c2 * integral
!          enddo
!  enddo
!  !$OMP END DO NOWAIT
!  !$OMP CRITICAL
!  ao_bi_elec_integral_alpha += ao_bi_elec_integral_alpha_tmp
!  !$OMP END CRITICAL
!  !$OMP CRITICAL
!  ao_bi_elec_integral_beta  += ao_bi_elec_integral_beta_tmp
!  !$OMP END CRITICAL
!  deallocate(keys,values,ao_bi_elec_integral_alpha_tmp,ao_bi_elec_integral_beta_tmp)
!  !$OMP END PARALLEL
!else
!  PROVIDE ao_bielec_integrals_in_map

!  integer(omp_lock_kind) :: lck(ao_num)
!  integer×8                      :: i8
!  integer                        :: ii(8), jj(8), kk(8), ll(8), k2
!  integer(cache_map_size_kind)   :: n_elements_max, n_elements
!  integer(key_kind), allocatable :: keys(:)
!  double precision, allocatable  :: values(:)
!  !$OMP PARALLEL DEFAULT(NONE)                                      &
!      !$OMP PRIVATE(i,j,l,k1,k,integral,ii,jj,kk,ll,i8,keys,values,n_elements_max, &
!      !$OMP  n_elements,ao_bi_elec_integral_alpha_tmp,ao_bi_elec_integral_beta_tmp)&
!      !$OMP SHARED(ao_num,SCF_density_matrix_ao_alpha,SCF_density_matrix_ao_beta,&
!      !$OMP  ao_integrals_map, ao_bi_elec_integral_alpha, ao_bi_elec_integral_beta) 

!  call get_cache_map_n_elements_max(ao_integrals_map,n_elements_max)
!  allocate(keys(n_elements_max), values(n_elements_max))
!  allocate(ao_bi_elec_integral_alpha_tmp(ao_num,ao_num), &
!           ao_bi_elec_integral_beta_tmp(ao_num,ao_num))
!  ao_bi_elec_integral_alpha_tmp = 0.d0
!  ao_bi_elec_integral_beta_tmp  = 0.d0

!  !$OMP DO SCHEDULE(dynamic,64)
!  !DIR$ NOVECTOR
!  do i8=0_8,ao_integrals_map%map_size
!    n_elements = n_elements_max
!    call get_cache_map(ao_integrals_map,i8,keys,values,n_elements)
!    do k1=1,n_elements
!      call bielec_integrals_index_reverse(kk,ii,ll,jj,keys(k1))

!      do k2=1,8
!        if (kk(k2)≡0) then
!          cycle
!        endif
!        i = ii(k2)
!        j = jj(k2)
!        k = kk(k2)
!        l = ll(k2)
!        integral = (SCF_density_matrix_ao_alpha(k,l)+SCF_density_matrix_ao_beta(k,l)) * values(k1)
!        ao_bi_elec_integral_alpha_tmp(i,j) += integral
!        ao_bi_elec_integral_beta_tmp (i,j) += integral
!        integral = values(k1)
!        ao_bi_elec_integral_alpha_tmp(l,j) -= SCF_density_matrix_ao_alpha(k,i) * integral
!        ao_bi_elec_integral_beta_tmp (l,j) -= SCF_density_matrix_ao_beta (k,i) * integral
!      enddo
!    enddo
!  enddo
!  !$OMP END DO NOWAIT
!  !$OMP CRITICAL
!  ao_bi_elec_integral_alpha += ao_bi_elec_integral_alpha_tmp
!  !$OMP END CRITICAL
!  !$OMP CRITICAL
!  ao_bi_elec_integral_beta  += ao_bi_elec_integral_beta_tmp
!  !$OMP END CRITICAL
!  deallocate(keys,values,ao_bi_elec_integral_alpha_tmp,ao_bi_elec_integral_beta_tmp)
!  !$OMP END PARALLEL
!endif
!END_PROVIDER
 

 BEGIN_PROVIDER [ integer, dirac_ao_num ]
 &BEGIN_PROVIDER [ integer, dirac_ao_nucl, (ao_num + small_ao_num) ]
 &BEGIN_PROVIDER [ double precision, dirac_ao_coef_normalized, (ao_num + small_ao_num) ]
 &BEGIN_PROVIDER [ double precision, dirac_ao_expo, (ao_num + small_ao_num) ]
 &BEGIN_PROVIDER [ integer, dirac_ao_power, (ao_num+small_ao_num,3) ]
  implicit none
  BEGIN_DOC
  ! Concatenation of the large and small components orbital properties
  ! in general arrays, for use in the bi-electronic integrals
  END_DOC
  integer                        :: i,j,k
  dirac_ao_num = (ao_num + small_ao_num)
  do i = 1, dirac_ao_num
   if (i .le. ao_num) then
    dirac_ao_nucl(i) = ao_nucl(i)   
    dirac_ao_coef_normalized(i) = ao_coef_normalized_ordered_transp(1,i)
    dirac_ao_expo(i) = ao_expo_ordered_transp(1,i)
    do k = 1, 3
     dirac_ao_power(i,k) = ao_power(i,k)
    enddo
   else 
    j = i - ao_num
    dirac_ao_nucl(i) = small_ao_nucl(j)
    dirac_ao_coef_normalized(i) = small_ao_coef_normalized(j)
    dirac_ao_expo(i) = small_ao_expo(j)
    do k = 1, 3
     dirac_ao_power(i,k) = small_ao_power(j,k)
    enddo
   endif
  enddo
 END_PROVIDER

 BEGIN_PROVIDER [double precision, dirac_ao_overlap_abs, (dirac_ao_num,dirac_ao_num) ]
 implicit none
  BEGIN_DOC
  ! Concatenation of the large and small component
  ! overlap_abs 
  END_DOC
  integer                        :: i,j,k,l
  do i = 1, dirac_ao_num
   if (i .le. ao_num) then
    do j = 1, dirac_ao_num
     if (j .le. ao_num) then
      dirac_ao_overlap_abs(i,j) = ao_overlap_abs(i,j)
     else 
      dirac_ao_overlap_abs(i,j) = 0.d0
     endif
    enddo
   else
    k = i - ao_num
    do j = 1, dirac_ao_num
     if (j .le. ao_num) then
     dirac_ao_overlap_abs(i,j) = 0.d0
     else
     l = j - ao_num
     dirac_ao_overlap_abs(i,j) = small_ao_overlap_abs(k,l)
     endif
    enddo 
   endif
  enddo
 END_PROVIDER
