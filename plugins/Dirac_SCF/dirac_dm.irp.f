 BEGIN_PROVIDER [ double precision, one_body_dm_mo_dirac, (2*dirac_mo_tot_num,2*dirac_mo_tot_num,N_states) ]
   implicit none
   BEGIN_DOC
   ! Alpha plus beta one-body density matrix for each state
   END_DOC
   integer                        :: j,k,l,m,k_a,k_b
   integer                        :: occ(N_int*bit_kind_size,2)
   double precision               :: ck, cl, ckl
   double precision               :: phase
   integer                        :: h1,h2,p1,p2,s1,s2, degree
   integer(bit_kind)              :: tmp_det(N_int,2), tmp_det2(N_int)
   integer                        :: exc(0:2,2),n_occ(2)
   double precision, allocatable  :: tmp_a(:,:,:), tmp_b(:,:,:)
   integer                        :: krow, kcol, lrow, lcol
 ! PROVIDE psi_det
 !one_body_dm_mo_alpha = 0.d0
 !one_body_dm_mo_beta  = 0.d0
 !!$OMP PARALLEL DEFAULT(NONE)                                         &
 !  !$OMP PRIVATE(j,k,k_a,k_b,l,m,occ,ck, cl, ckl,phase,h1,h2,p1,p2,s1,s2, degree,exc, &
 !  !$OMP  tmp_a, tmp_b, n_occ, krow, kcol, lrow, lcol, tmp_det, tmp_det2)&
 !  !$OMP SHARED(psi_det,psi_coef,N_int,N_states,elec_alpha_num,&
 !  !$OMP  elec_beta_num,one_body_dm_mo_alpha,one_body_dm_mo_beta,N_det,&
 !  !$OMP  mo_tot_num,psi_bilinear_matrix_rows,psi_bilinear_matrix_columns, &
 !  !$OMP  psi_bilinear_matrix_transp_rows, psi_bilinear_matrix_transp_columns, &
 !  !$OMP  psi_bilinear_matrix_order_reverse, psi_det_alpha_unique, psi_det_beta_unique, &
 !  !$OMP  psi_bilinear_matrix_values, psi_bilinear_matrix_transp_values, &
 !  !$OMP  N_det_alpha_unique,N_det_beta_unique,irp_here)
 !allocate(tmp_a(mo_tot_num,mo_tot_num,N_states), tmp_b(mo_tot_num,mo_tot_num,N_states) )
 !tmp_a = 0.d0
 !!$OMP DO SCHEDULE(dynamic,64)
 !do k_a=1,N_det
 !  krow = psi_bilinear_matrix_rows(k_a) 
 !  ASSERT (krow <= N_det_alpha_unique)

 !  kcol = psi_bilinear_matrix_columns(k_a) 
 !  ASSERT (kcol <= N_det_beta_unique)

 !  tmp_det(1:N_int,1) = psi_det_alpha_unique(1:N_int,krow)
 !  tmp_det(1:N_int,2) = psi_det_beta_unique (1:N_int,kcol)

 !  ! Diagonal part
 !  ! -------------

 !  call bitstring_to_list_ab(tmp_det, occ, n_occ, N_int)
 !  do m=1,N_states
 !    ck = psi_bilinear_matrix_values(k_a,m)*psi_bilinear_matrix_values(k_a,m)
 !    do l=1,elec_alpha_num
 !      j = occ(l,1)
 !      tmp_a(j,j,m) += ck
 !    enddo
 !  enddo

 !  if (k_a == N_det) cycle
 !  l = k_a+1
 !  lrow = psi_bilinear_matrix_rows(l) 
 !  lcol = psi_bilinear_matrix_columns(l) 
 !  ! Fix beta determinant, loop over alphas
 !  do while ( lcol == kcol )
 !    tmp_det2(:) = psi_det_alpha_unique(:, lrow)
 !    call get_excitation_degree_spin(tmp_det(1,1),tmp_det2,degree,N_int)
 !    if (degree == 1) then
 !      exc = 0
 !      call get_mono_excitation_spin(tmp_det(1,1),tmp_det2,exc,phase,N_int)
 !      call decode_exc_spin(exc,h1,p1,h2,p2)
 !      do m=1,N_states
 !        ckl = psi_bilinear_matrix_values(k_a,m)*psi_bilinear_matrix_values(l,m) * phase
 !        tmp_a(h1,p1,m) += ckl
 !        tmp_a(p1,h1,m) += ckl
 !      enddo
 !    endif
 !    l = l+1
 !    if (l>N_det) exit
 !    lrow = psi_bilinear_matrix_rows(l) 
 !    lcol = psi_bilinear_matrix_columns(l) 
 !  enddo

 !enddo
 !!$OMP END DO NOWAIT

 !!$OMP CRITICAL
 !one_body_dm_mo_alpha(:,:,:) = one_body_dm_mo_alpha(:,:,:) + tmp_a(:,:,:)
 !!$OMP END CRITICAL
 !deallocate(tmp_a)

 !tmp_b = 0.d0
 !!$OMP DO SCHEDULE(dynamic,64)
 !do k_b=1,N_det
 !  krow = psi_bilinear_matrix_transp_rows(k_b) 
 !  ASSERT (krow <= N_det_alpha_unique)

 !  kcol = psi_bilinear_matrix_transp_columns(k_b) 
 !  ASSERT (kcol <= N_det_beta_unique)

 !  tmp_det(1:N_int,1) = psi_det_alpha_unique(1:N_int,krow)
 !  tmp_det(1:N_int,2) = psi_det_beta_unique (1:N_int,kcol)

 !  ! Diagonal part
 !  ! -------------

 !  call bitstring_to_list_ab(tmp_det, occ, n_occ, N_int)
 !  do m=1,N_states
 !    ck = psi_bilinear_matrix_transp_values(k_b,m)*psi_bilinear_matrix_transp_values(k_b,m)
 !    do l=1,elec_beta_num
 !      j = occ(l,2)
 !      tmp_b(j,j,m) += ck
 !    enddo
 !  enddo

 !  if (k_b == N_det) cycle
 !  l = k_b+1
 !  lrow = psi_bilinear_matrix_transp_rows(l) 
 !  lcol = psi_bilinear_matrix_transp_columns(l) 
 !  ! Fix beta determinant, loop over alphas
 !  do while ( lrow == krow )
 !    tmp_det2(:) = psi_det_beta_unique(:, lcol)
 !    call get_excitation_degree_spin(tmp_det(1,2),tmp_det2,degree,N_int)
 !    if (degree == 1) then
 !      exc = 0
 !      call get_mono_excitation_spin(tmp_det(1,2),tmp_det2,exc,phase,N_int)
 !      call decode_exc_spin(exc,h1,p1,h2,p2)
 !      do m=1,N_states
 !        ckl = psi_bilinear_matrix_transp_values(k_b,m)*psi_bilinear_matrix_transp_values(l,m) * phase
 !        tmp_b(h1,p1,m) += ckl
 !        tmp_b(p1,h1,m) += ckl
 !      enddo
 !    endif
 !    l = l+1
 !    if (l>N_det) exit
 !    lrow = psi_bilinear_matrix_transp_rows(l) 
 !    lcol = psi_bilinear_matrix_transp_columns(l) 
 !  enddo

 !enddo
 !!$OMP END DO NOWAIT
 !!$OMP CRITICAL
 !one_body_dm_mo_beta(:,:,:)  = one_body_dm_mo_beta(:,:,:)  + tmp_b(:,:,:)
 !!$OMP END CRITICAL

 !deallocate(tmp_b)
 !!$OMP END PARALLEL

END_PROVIDER


 BEGIN_PROVIDER [double precision, one_body_dm_dirac_mo_for_dft, (2*dirac_mo_tot_num,2*dirac_mo_tot_num, N_states)]
 implicit none
 BEGIN_DOC
 ! density used for all DFT calculations based on the density 
 END_DOC
 !if(read_density_from_input)then
 ! one_body_dm_alpha_mo_for_dft = data_one_body_alpha_dm_mo
 !else 
   one_body_dm_dirac_mo_for_dft = one_body_dm_mo_dirac
 !endif
 END_PROVIDER


 
 BEGIN_PROVIDER [ double precision, one_body_dm_dirac_ao_for_dft, (2*dirac_ao_num,2*dirac_ao_num,N_states) ]
 BEGIN_DOC
 ! one body density matrix on the AO basis based on one_body_dm_alpha_mo_for_dft
 ! Not sure if it can be kept as a double precision
 END_DOC
  implicit none
  integer :: i,j,k,l,istate
  double precision :: mo_dirac
  one_body_dm_dirac_ao_for_dft = 0.d0
  do k = 1, 2*dirac_ao_num
   do l = 1, 2*dirac_ao_num
    do i = 1, 2*dirac_mo_tot_num
     do j = 1, 2*dirac_mo_tot_num
      do istate = 1, N_states
       mo_dirac = one_body_dm_dirac_mo_for_dft(j,i,istate)
       one_body_dm_dirac_ao_for_dft(l,k,istate) += Conjg(dirac_mo_coef(k,i)) * dirac_mo_coef(l,j) *  mo_dirac
       !Not sure about the conjugate but most likely it is needed
      enddo
     enddo
    enddo
   enddo
  enddo
 END_PROVIDER
