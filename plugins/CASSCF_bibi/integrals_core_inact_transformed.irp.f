
  use map_module
 BEGIN_PROVIDER [real(integral_kind), transformed_occ1_virt2_virt2, (n_core_inact_orb,n_virt_orb,n_virt_orb)]
  use map_module
 implicit none
 BEGIN_DOC
! transformed_occ1_virt2_virt2(k,j,i) = (k_core k_core | j_virt i_virt )
! BE CAREFULL :::: TO GET BACK THE INDEX OF ORBITALS YOU NEED TO USE 
! i_real = list_core_inact(i_core)
! j_real = list_core_inact(j_virt)
 END_DOC
 transformed_occ1_virt2_virt2 = 0.d0

 double precision :: cpu0,cpu1
 integer :: i,j,k,l,iorb,jorb,korb,lorb,m,n,p,q
 double precision :: c
 real(integral_kind), allocatable :: bielec_tmp_0(:,:),matrix_tmp_1(:,:)
 real(integral_kind) :: integral,ao_bielec_integral
 real(integral_kind) :: thr
 

 provide ao_bielec_integral_schwartz ao_bielec_integrals_in_map ao_overlap_abs
 provide semi_transformed_virt_virt 
 call cpu_time(cpu0) 
 thr = 0.001d0 * dsqrt(ao_integrals_threshold)
 thr = 0.d0

  !$OMP PARALLEL DEFAULT(NONE)             &
  !$OMP PRIVATE(i,j,k,m,bielec_tmp_0,matrix_tmp_1) &
  !$OMP SHARED(n_core_inact_orb,ao_num,semi_transformed_virt_virt,n_virt_orb,mo_coef_core_inact_transp,transformed_occ1_virt2_virt2)
 
  allocate(bielec_tmp_0(ao_num,ao_num),matrix_tmp_1(n_core_inact_orb,ao_num))
  !$OMP DO SCHEDULE(guided)
 do i = 1, n_virt_orb
  do j = 1, n_virt_orb
   bielec_tmp_0(1:ao_num,1:ao_num) = semi_transformed_virt_virt(j,i,1:ao_num,1:ao_num)
   matrix_tmp_1 = 0.d0
   call dgemm('N','N',n_core_inact_orb,ao_num,ao_num,1.d0,mo_coef_core_inact_transp,n_core_inact_orb,bielec_tmp_0,ao_num,0.d0,matrix_tmp_1,n_core_inact_orb)
   do m = 1, ao_num
    do k = 1, n_core_inact_orb
     transformed_occ1_virt2_virt2(k,j,i) += matrix_tmp_1(k,m) * mo_coef_core_inact_transp(k,m) 
    enddo
   enddo
  enddo
 enddo
  !$OMP END PARALLEL
 
 call cpu_time(cpu1) 
 print*, 'Time to transformed_occ1_virt2_virt2 =',cpu1-cpu0

 END_PROVIDER 


 BEGIN_PROVIDER [real(integral_kind), transformed_virt1_occ2_occ2, (n_virt_orb,n_core_inact_orb,n_core_inact_orb)]
  use map_module
 implicit none
 BEGIN_DOC
! transformed_virt1_occ2_occ2(k,j,i) = (k_virt k_virt | j_core i_core )
! BE CAREFULL :::: TO GET BACK THE INDEX OF ORBITALS YOU NEED TO USE 
! i_real = list_core_inact(i_core)
! j_real = list_core_inact(j_virt)
 END_DOC
 transformed_virt1_occ2_occ2 = 0.d0

 double precision :: cpu0,cpu1
 integer :: i,j,k,l,iorb,jorb,korb,lorb,m,n,p,q
 double precision :: c
 real(integral_kind), allocatable :: bielec_tmp_0(:,:),matrix_tmp_1(:,:)
 real(integral_kind) :: integral,ao_bielec_integral
 real(integral_kind) :: thr
 

 provide ao_bielec_integral_schwartz ao_bielec_integrals_in_map ao_overlap_abs semi_transformed_occ_occ
 call cpu_time(cpu0) 
 thr = 0.001d0 * dsqrt(ao_integrals_threshold)

  !$OMP PARALLEL DEFAULT(NONE)             &
  !$OMP PRIVATE(i,j,k,m,bielec_tmp_0,matrix_tmp_1) &
  !$OMP SHARED(n_core_inact_orb,ao_num,semi_transformed_occ_occ,n_virt_orb,mo_coef_virt_transp,transformed_virt1_occ2_occ2)
 
  allocate(bielec_tmp_0(ao_num,ao_num),matrix_tmp_1(n_virt_orb,ao_num))
  !$OMP DO SCHEDULE(guided)
 do i = 1, n_core_inact_orb
  do j = 1, n_core_inact_orb
   bielec_tmp_0(1:ao_num,1:ao_num) = semi_transformed_occ_occ(j,i,1:ao_num,1:ao_num)
   call dgemm('N','N',n_virt_orb,ao_num,ao_num,1.d0,mo_coef_virt_transp,n_virt_orb,bielec_tmp_0,ao_num,0.d0,matrix_tmp_1,n_virt_orb)
   do m = 1, ao_num
    do k = 1, n_virt_orb
     transformed_virt1_occ2_occ2(k,j,i) += matrix_tmp_1(k,m) * mo_coef_virt_transp(k,m) 
    enddo
   enddo
  enddo
 enddo
  !$OMP END PARALLEL
 
 call cpu_time(cpu1) 
 print*, 'Time to transformed_virt1_occ2_occ2 =',cpu1-cpu0

 END_PROVIDER 






 BEGIN_PROVIDER [real(integral_kind), transformed_occ1_occ2_virt2, (n_core_inact_orb,n_core_inact_orb,n_virt_orb)]
  use map_module
 implicit none
 BEGIN_DOC
! transformed_occ_virt_occ_from_semi_trans(i,j,k) = (i_core i_core | j_core k_virt)
! BE CAREFULL :::: TO GET BACK THE INDEX OF ORBITALS YOU NEED TO USE 
! i_real = list_core_inact(i_core)
! j_real = list_core_inact(j_virt)
 END_DOC
 transformed_occ1_occ2_virt2 = 0.d0

 double precision :: cpu0,cpu1
 integer :: i,j,k,l,iorb,jorb,korb,lorb,m,n,p,q
 double precision :: c
 real(integral_kind), allocatable :: bielec_tmp_0(:,:),matrix_tmp_1(:,:),matrix_final(:,:),matrix_tmp_2(:,:)
 real(integral_kind) :: integral,ao_bielec_integral
 real(integral_kind) :: thr
 

 provide ao_bielec_integral_schwartz ao_bielec_integrals_in_map ao_overlap_abs semi_transformed_occ_virt
 call cpu_time(cpu0) 
 thr = 0.001d0 * dsqrt(ao_integrals_threshold)
  !$OMP PARALLEL DEFAULT(NONE)             &
  !$OMP PRIVATE(i,j,k,m,bielec_tmp_0,matrix_tmp_1) &
  !$OMP SHARED(n_core_inact_orb,ao_num,semi_transformed_occ_virt,n_virt_orb,mo_coef_core_inact_transp,transformed_occ1_occ2_virt2)
 
  allocate(bielec_tmp_0(ao_num,ao_num),matrix_tmp_1(n_core_inact_orb,ao_num))
  !$OMP DO SCHEDULE(guided)
 do i = 1, n_virt_orb
  do j = 1, n_core_inact_orb
   bielec_tmp_0(1:ao_num,1:ao_num) = semi_transformed_occ_virt(j,i,1:ao_num,1:ao_num)
   call dgemm('N','N',n_core_inact_orb,ao_num,ao_num,1.d0,mo_coef_core_inact_transp,n_core_inact_orb,bielec_tmp_0,ao_num,0.d0,matrix_tmp_1,n_core_inact_orb)
   do m = 1, ao_num
    do k = 1, n_core_inact_orb
     transformed_occ1_occ2_virt2(k,j,i) += matrix_tmp_1(k,m) * mo_coef_core_inact_transp(k,m) 
    enddo
   enddo
  enddo
 enddo
  !$OMP END PARALLEL
 
 call cpu_time(cpu1) 
 print*, 'Time to transformed_occ1_occ2_virt2 =',cpu1-cpu0

 END_PROVIDER 




 BEGIN_PROVIDER [real(integral_kind), transformed_occ1_virt1_occ2_virt2, (n_core_inact_orb,n_virt_orb,n_core_inact_orb,n_virt_orb)]
  use map_module
 implicit none
 BEGIN_DOC
! transformed_occ1_virt1_occ2_virt2(i,j,k,l) = (i_core j_virt | k_core j_virt)
! BE CAREFULL :::: TO GET BACK THE INDEX OF ORBITALS YOU NEED TO USE 
! i_real = list_core_inact(i_core)
! j_real = list_core_inact(j_virt)
 END_DOC
 transformed_occ1_virt1_occ2_virt2 = 0.d0

 double precision :: cpu0,cpu1
 integer :: i,j,k,l,iorb,jorb,korb,lorb,m,n,p,q
 double precision :: c
 real(integral_kind), allocatable :: bielec_tmp_0(:,:),matrix_tmp_1(:,:),matrix_final(:,:),matrix_tmp_2(:,:)
 real(integral_kind) :: integral,ao_bielec_integral
 real(integral_kind) :: thr
 

 provide ao_bielec_integral_schwartz ao_bielec_integrals_in_map ao_overlap_abs semi_transformed_occ_virt
 call cpu_time(cpu0) 
 thr = 0.001d0 * dsqrt(ao_integrals_threshold)
  !$OMP PARALLEL DEFAULT(NONE)             &
  !$OMP PRIVATE(i,j,bielec_tmp_0,matrix_tmp_1,matrix_final) &
  !$OMP SHARED(n_virt_orb,n_core_inact_orb,ao_num,semi_transformed_occ_virt,mo_coef_core_inact_transp,mo_coef_virt,transformed_occ1_virt1_occ2_virt2)
 
  allocate(bielec_tmp_0(ao_num,ao_num),matrix_tmp_1(n_core_inact_orb,ao_num),matrix_final(n_core_inact_orb,n_virt_orb))
  !$OMP DO SCHEDULE(guided)
 do i = 1, n_virt_orb
  do j = 1, n_core_inact_orb
   bielec_tmp_0(1:ao_num,1:ao_num) = semi_transformed_occ_virt(j,i,1:ao_num,1:ao_num)
!  bielec_tmp_0(1:ao_num,1:ao_num) = semi_transformed_occ_virt(1:ao_num,1:ao_num,j,i)
   call dgemm('N','N',n_core_inact_orb,ao_num,ao_num,1.d0,mo_coef_core_inact_transp,n_core_inact_orb,bielec_tmp_0,ao_num,0.d0,matrix_tmp_1,n_core_inact_orb)
   call dgemm('N','N',n_core_inact_orb,n_virt_orb,ao_num,1.d0,matrix_tmp_1,n_core_inact_orb,mo_coef_virt,ao_num,0.d0,matrix_final,n_core_inact_orb)
   transformed_occ1_virt1_occ2_virt2(1:n_core_inact_orb,1:n_virt_orb,j,i) = matrix_final(1:n_core_inact_orb,1:n_virt_orb) 
  enddo
 enddo
  !$OMP END PARALLEL
 
 call cpu_time(cpu1) 
 print*, 'Time to transformed_occ1_virt1_occ2_virt2 =',cpu1-cpu0
!free semi_transformed_occ_virt

 END_PROVIDER 



 BEGIN_PROVIDER [real(integral_kind), transformed_occ1_virt1_occ2_virt2, (n_core_inact_orb,n_virt_orb,n_core_inact_orb,n_virt_orb)]
  use map_module
 implicit none
 BEGIN_DOC
! transformed_occ1_virt1_occ2_virt2(i,j,k,l) = (i_core j_virt | k_core j_virt)
! BE CAREFULL :::: TO GET BACK THE INDEX OF ORBITALS YOU NEED TO USE 
! i_real = list_core_inact(i_core)
! j_real = list_core_inact(j_virt)
 END_DOC
 transformed_occ1_virt1_occ2_virt2 = 0.d0

 double precision :: cpu0,cpu1
 integer :: i,j,k,l,iorb,jorb,korb,lorb,m,n,p,q
 double precision :: c
 real(integral_kind), allocatable :: bielec_tmp_0(:,:),matrix_tmp_1(:,:),matrix_final(:,:),matrix_tmp_2(:,:)
 real(integral_kind) :: integral,ao_bielec_integral
 real(integral_kind) :: thr
 

 provide ao_bielec_integral_schwartz ao_bielec_integrals_in_map ao_overlap_abs semi_transformed_occ_virt
 call cpu_time(cpu0) 
 thr = 0.001d0 * dsqrt(ao_integrals_threshold)
  !$OMP PARALLEL DEFAULT(NONE)             &
  !$OMP PRIVATE(i,j,bielec_tmp_0,matrix_tmp_1,matrix_final) &
  !$OMP SHARED(n_virt_orb,n_core_inact_orb,ao_num,semi_transformed_occ_virt,mo_coef_core_inact_transp,mo_coef_virt,transformed_occ1_virt1_occ2_virt2)
 
  allocate(bielec_tmp_0(ao_num,ao_num),matrix_tmp_1(n_core_inact_orb,ao_num),matrix_final(n_core_inact_orb,n_virt_orb))
  !$OMP DO SCHEDULE(guided)
 do i = 1, n_virt_orb
  do j = 1, n_core_inact_orb
   bielec_tmp_0(1:ao_num,1:ao_num) = semi_transformed_occ_virt(j,i,1:ao_num,1:ao_num)
!  bielec_tmp_0(1:ao_num,1:ao_num) = semi_transformed_occ_virt(1:ao_num,1:ao_num,j,i)
   call dgemm('N','N',n_core_inact_orb,ao_num,ao_num,1.d0,mo_coef_core_inact_transp,n_core_inact_orb,bielec_tmp_0,ao_num,0.d0,matrix_tmp_1,n_core_inact_orb)
   call dgemm('N','N',n_core_inact_orb,n_virt_orb,ao_num,1.d0,matrix_tmp_1,n_core_inact_orb,mo_coef_virt,ao_num,0.d0,matrix_final,n_core_inact_orb)
   transformed_occ1_virt1_occ2_virt2(1:n_core_inact_orb,1:n_virt_orb,j,i) = matrix_final(1:n_core_inact_orb,1:n_virt_orb) 
  enddo
 enddo
  !$OMP END PARALLEL
 
 call cpu_time(cpu1) 
 print*, 'Time to transformed_occ1_virt1_occ2_virt2 =',cpu1-cpu0
!free semi_transformed_occ_virt

 END_PROVIDER 



 BEGIN_PROVIDER [real(integral_kind), transformed_virt1_virt1_occ2_occ2, (n_virt_orb,n_virt_orb,n_core_inact_orb,n_core_inact_orb)]
  use map_module
 implicit none
 BEGIN_DOC
! transformed_virt1_virt1_occ2_occ2(i,j,k,l) = (i_core j_core | k_virt j_virt)
! BE CAREFULL :::: TO GET BACK THE INDEX OF ORBITALS YOU NEED TO USE 
! i_real = list_core_inact(i_core)
! j_real = list_core_inact(j_virt)
 END_DOC
 transformed_virt1_virt1_occ2_occ2 = 0.d0

 double precision :: cpu0,cpu1
 integer :: i,j,k,l,iorb,jorb,korb,lorb,m,n,p,q
 double precision :: c
 real(integral_kind), allocatable :: bielec_tmp_0(:,:),matrix_tmp_1(:,:),matrix_final(:,:),matrix_tmp_2(:,:)
 real(integral_kind) :: integral,ao_bielec_integral
 real(integral_kind) :: thr
 

 provide ao_bielec_integral_schwartz ao_bielec_integrals_in_map ao_overlap_abs semi_transformed_occ_occ
 call cpu_time(cpu0) 
 thr = 0.001d0 * dsqrt(ao_integrals_threshold)
  !$OMP PARALLEL DEFAULT(NONE)             &
  !$OMP PRIVATE(i,j,bielec_tmp_0,matrix_tmp_1,matrix_final) &
  !$OMP SHARED(n_virt_orb,n_core_inact_orb,ao_num,semi_transformed_occ_occ,mo_coef_virt_transp,mo_coef_virt,transformed_virt1_virt1_occ2_occ2)
 
  allocate(bielec_tmp_0(ao_num,ao_num),matrix_tmp_1(n_virt_orb,ao_num),matrix_final(n_virt_orb,n_virt_orb))
  !$OMP DO SCHEDULE(guided)
 do i = 1, n_core_inact_orb
  do j = 1, n_core_inact_orb
   bielec_tmp_0(1:ao_num,1:ao_num) = semi_transformed_occ_occ(j,i,1:ao_num,1:ao_num)
   call dgemm('N','N',n_virt_orb,ao_num,ao_num,1.d0,mo_coef_virt_transp,n_virt_orb,bielec_tmp_0,ao_num,0.d0,matrix_tmp_1,n_virt_orb)
   call dgemm('N','N',n_virt_orb,n_virt_orb,ao_num,1.d0,matrix_tmp_1,n_virt_orb,mo_coef_virt,ao_num,0.d0,matrix_final,n_virt_orb)
   transformed_virt1_virt1_occ2_occ2(1:n_virt_orb,1:n_virt_orb,j,i) = matrix_final(1:n_virt_orb,1:n_virt_orb) 
  enddo
 enddo
  !$OMP END PARALLEL
 
 call cpu_time(cpu1) 
 print*, 'Time to transformed_virt1_virt1_occ2_occ2 =',cpu1-cpu0
!free semi_transformed_occ_virt

 END_PROVIDER 



 BEGIN_PROVIDER [real(integral_kind), transformed_act1_act1_act2_act2, (n_act_orb,n_act_orb,n_act_orb,n_act_orb)]
  use map_module
 implicit none
 BEGIN_DOC
! transformed_act1_act1_act2_act2(i,j,k,l) = (i_core j_core | k_act j_act)
! BE CAREFULL :::: TO GET BACK THE INDEX OF ORBITALS YOU NEED TO USE 
! i_real = list_act(i_core)
! j_real = list_act(j_act)
 END_DOC
 transformed_act1_act1_act2_act2 = 0.d0

 double precision :: cpu0,cpu1
 integer :: i,j,k,l,iorb,jorb,korb,lorb,m,n,p,q
 double precision :: c
 real(integral_kind), allocatable :: bielec_tmp_0(:,:),matrix_tmp_1(:,:),matrix_final(:,:),matrix_tmp_2(:,:)
 real(integral_kind) :: integral,ao_bielec_integral
 real(integral_kind) :: thr
 

 provide ao_bielec_integral_schwartz ao_bielec_integrals_in_map ao_overlap_abs semi_transformed_act_act
 call cpu_time(cpu0) 
 thr = 0.001d0 * dsqrt(ao_integrals_threshold)
  !$OMP PARALLEL DEFAULT(NONE)             &
  !$OMP PRIVATE(i,j,bielec_tmp_0,matrix_tmp_1,matrix_final) &
  !$OMP SHARED(n_act_orb,ao_num,semi_transformed_act_act,mo_coef_act_transp,mo_coef_act,transformed_act1_act1_act2_act2)
 
  allocate(bielec_tmp_0(ao_num,ao_num),matrix_tmp_1(n_act_orb,ao_num),matrix_final(n_act_orb,n_act_orb))
  !$OMP DO SCHEDULE(guided)
 do i = 1, n_act_orb
  do j = 1, n_act_orb
   bielec_tmp_0(1:ao_num,1:ao_num) = semi_transformed_act_act(j,i,1:ao_num,1:ao_num)
   call dgemm('N','N',n_act_orb,ao_num,ao_num,1.d0,mo_coef_act_transp,n_act_orb,bielec_tmp_0,ao_num,0.d0,matrix_tmp_1,n_act_orb)
   call dgemm('N','N',n_act_orb,n_act_orb,ao_num,1.d0,matrix_tmp_1,n_act_orb,mo_coef_act,ao_num,0.d0,matrix_final,n_act_orb)
   transformed_act1_act1_act2_act2(1:n_act_orb,1:n_act_orb,j,i) = matrix_final(1:n_act_orb,1:n_act_orb) 
  enddo
 enddo
  !$OMP END PARALLEL
 
 call cpu_time(cpu1) 
 print*, 'Time to transformed_act1_act1_act2_act2 =',cpu1-cpu0
!free semi_transformed_act_act

 END_PROVIDER 



