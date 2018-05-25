
  use map_module
 BEGIN_PROVIDER [real(integral_kind), semi_transformed_occ_virt, (n_core_inact_orb,n_virt_orb,ao_num,ao_num)]
  use map_module
 implicit none
 semi_transformed_occ_virt = 0.d0


 integer :: i,j,k,l,iorb,jorb,korb,lorb,m,n,p,q
 double precision :: c
 real(integral_kind), allocatable :: bielec_tmp_0(:,:),matrix_tmp_1(:,:),matrix_final(:,:)
 real(integral_kind) :: integral,ao_bielec_integral,thr

 double precision :: cpu0,cpu1
 
 provide ao_bielec_integral_schwartz ao_bielec_integrals_in_map ao_overlap_abs mo_coef_core_inact_transp mo_coef_virt
 thr = 0.001d0 * dsqrt(ao_integrals_threshold)
 call cpu_time(cpu0) 

 !$OMP PARALLEL DEFAULT(NONE)             &
 !$OMP PRIVATE(i,j,p,q,bielec_tmp_0,matrix_tmp_1,matrix_final) &
 !$OMP SHARED(n_virt_orb,n_core_inact_orb,ao_num,semi_transformed_occ_virt,mo_coef_core_inact_transp,mo_coef_virt,thr, &
 !$OMP        ao_overlap_abs,ao_bielec_integral_schwartz)

 allocate(bielec_tmp_0(ao_num,ao_num),matrix_tmp_1(n_core_inact_orb,ao_num),matrix_final(n_core_inact_orb,n_virt_orb))

  !$OMP DO SCHEDULE(guided)
 do p = 1,ao_num
   do q = 1,ao_num
    if(ao_overlap_abs(p,q).le.thr)cycle
    if(ao_bielec_integral_schwartz(p,q).lt.thr) cycle
    matrix_tmp_1 = 0.d0
    matrix_final = 0.d0
    do m = 1, ao_num
     call get_ao_bielec_integrals(p,m,q,ao_num,bielec_tmp_0(1,m)) ! k,l :: r1, m :: r2
    enddo
    call dgemm('N','N',n_core_inact_orb,ao_num,ao_num,1.d0,mo_coef_core_inact_transp,n_core_inact_orb,bielec_tmp_0,ao_num,0.d0,matrix_tmp_1,n_core_inact_orb)
    call dgemm('N','N',n_core_inact_orb,n_virt_orb,ao_num,1.d0,matrix_tmp_1,n_core_inact_orb,mo_coef_virt,ao_num,0.d0,matrix_final,n_core_inact_orb)

    semi_transformed_occ_virt(1:n_core_inact_orb,1:n_virt_orb,p,q) = matrix_final(1:n_core_inact_orb,1:n_virt_orb)
   enddo
  enddo
  !$OMP END PARALLEL

 call cpu_time(cpu1) 
 print*, 'Time to do the semi transformation core - virt : ',cpu1-cpu0
 
END_PROVIDER 


 BEGIN_PROVIDER [real(integral_kind), semi_transformed_virt_virt, (n_virt_orb,n_virt_orb,ao_num,ao_num)]
  use map_module
 implicit none
 semi_transformed_virt_virt = 0.d0


 integer :: i,j,k,l,iorb,jorb,korb,lorb,m,n,p,q
 double precision :: c
 real(integral_kind), allocatable :: bielec_tmp_0(:,:),matrix_tmp_1(:,:),matrix_final(:,:)
 real(integral_kind) :: integral,ao_bielec_integral,thr

 double precision :: cpu0,cpu1
 
 provide ao_bielec_integral_schwartz ao_bielec_integrals_in_map ao_overlap_abs mo_coef_virt  mo_coef_virt_transp
 thr = 0.001d0 * dsqrt(ao_integrals_threshold)
 call cpu_time(cpu0) 

 !$OMP PARALLEL DEFAULT(NONE)             &
 !$OMP PRIVATE(i,j,p,q,bielec_tmp_0,matrix_tmp_1,matrix_final) &
 !$OMP SHARED(n_virt_orb,ao_num,semi_transformed_virt_virt,mo_coef_virt_transp,mo_coef_virt,thr, &
 !$OMP        ao_overlap_abs,ao_bielec_integral_schwartz)

 allocate(bielec_tmp_0(ao_num,ao_num),matrix_tmp_1(n_virt_orb,ao_num),matrix_final(n_virt_orb,n_virt_orb))

  !$OMP DO SCHEDULE(guided)
 do p = 1,ao_num
   do q = 1,ao_num
    if(ao_overlap_abs(p,q).le.thr)cycle
    if(ao_bielec_integral_schwartz(p,q).lt.thr) cycle
    matrix_tmp_1 = 0.d0
    matrix_final = 0.d0
    do m = 1, ao_num
     call get_ao_bielec_integrals(p,m,q,ao_num,bielec_tmp_0(1,m)) ! k,l :: r1, m :: r2
    enddo
    call dgemm('N','N',n_virt_orb,ao_num,ao_num,1.d0,mo_coef_virt_transp,n_virt_orb,bielec_tmp_0,ao_num,0.d0,matrix_tmp_1,n_virt_orb)
    call dgemm('N','N',n_virt_orb,n_virt_orb,ao_num,1.d0,matrix_tmp_1,n_virt_orb,mo_coef_virt,ao_num,0.d0,matrix_final,n_virt_orb)

    semi_transformed_virt_virt(1:n_virt_orb,1:n_virt_orb,p,q) = matrix_final(1:n_virt_orb,1:n_virt_orb)
   enddo
  enddo
  !$OMP END PARALLEL

 call cpu_time(cpu1) 
 print*, 'Time to do the semi transformation virt - virt : ',cpu1-cpu0
 
END_PROVIDER 




 BEGIN_PROVIDER [real(integral_kind), semi_transformed_occ_occ, (n_core_inact_orb,n_core_inact_orb,ao_num,ao_num)]
  use map_module
 implicit none
 semi_transformed_occ_occ  = 0.d0


 integer :: i,j,k,l,iorb,jorb,korb,lorb,m,n,p,q
 double precision :: c
 real(integral_kind), allocatable :: bielec_tmp_0(:,:),matrix_tmp_1(:,:),matrix_final(:,:)
 real(integral_kind) :: integral,ao_bielec_integral,thr

 double precision :: cpu0,cpu1
 
 provide ao_bielec_integral_schwartz ao_bielec_integrals_in_map ao_overlap_abs mo_coef_core_inact
 thr = 0.001d0 * dsqrt(ao_integrals_threshold)
 call cpu_time(cpu0) 

 !$OMP PARALLEL DEFAULT(NONE)             &
 !$OMP PRIVATE(i,j,p,q,bielec_tmp_0,matrix_tmp_1,matrix_final) &
 !$OMP SHARED(n_core_inact_orb,ao_num,semi_transformed_occ_occ,mo_coef_core_inact_transp,mo_coef_core_inact,thr, &
 !$OMP        ao_overlap_abs,ao_bielec_integral_schwartz)

 allocate(bielec_tmp_0(ao_num,ao_num),matrix_tmp_1(n_core_inact_orb,ao_num),matrix_final(n_core_inact_orb,n_core_inact_orb))

  !$OMP DO SCHEDULE(guided)
 do p = 1,ao_num
   do q = 1,ao_num
    if(ao_overlap_abs(p,q).le.thr)cycle
    if(ao_bielec_integral_schwartz(p,q).lt.thr) cycle
    matrix_tmp_1 = 0.d0
    matrix_final = 0.d0
    do m = 1, ao_num
     call get_ao_bielec_integrals(p,m,q,ao_num,bielec_tmp_0(1,m)) ! k,l :: r1, m :: r2
    enddo
    call dgemm('N','N',n_core_inact_orb,ao_num,ao_num,1.d0,mo_coef_core_inact_transp,n_core_inact_orb,bielec_tmp_0,ao_num,0.d0,matrix_tmp_1,n_core_inact_orb)
    call dgemm('N','N',n_core_inact_orb,n_core_inact_orb,ao_num,1.d0,matrix_tmp_1,n_core_inact_orb,mo_coef_core_inact,ao_num,0.d0,matrix_final,n_core_inact_orb)

    semi_transformed_occ_occ(1:n_core_inact_orb,1:n_core_inact_orb,p,q) = matrix_final(1:n_core_inact_orb,1:n_core_inact_orb)
   enddo
  enddo
  !$OMP END PARALLEL

 call cpu_time(cpu1) 
 print*, 'Time to do the semi transformation core_inact - core_inact : ',cpu1-cpu0
 
END_PROVIDER 

 BEGIN_PROVIDER [real(integral_kind), semi_transformed_act_act, (n_act_orb,n_act_orb,ao_num,ao_num)]
  use map_module
 implicit none
 semi_transformed_act_act = 0.d0


 integer :: i,j,k,l,iorb,jorb,korb,lorb,m,n,p,q
 double precision :: c
 real(integral_kind), allocatable :: bielec_tmp_0(:,:),matrix_tmp_1(:,:),matrix_final(:,:)
 real(integral_kind) :: integral,ao_bielec_integral,thr

 double precision :: cpu0,cpu1
 
 provide ao_bielec_integral_schwartz ao_bielec_integrals_in_map ao_overlap_abs mo_coef_act_transp mo_coef_act 
 thr = 0.001d0 * dsqrt(ao_integrals_threshold)
 call cpu_time(cpu0) 

 !$OMP PARALLEL DEFAULT(NONE)             &
 !$OMP PRIVATE(i,j,p,q,bielec_tmp_0,matrix_tmp_1,matrix_final) &
 !$OMP SHARED(n_act_orb,ao_num,semi_transformed_act_act,mo_coef_act_transp,mo_coef_act,thr, &
 !$OMP        ao_overlap_abs,ao_bielec_integral_schwartz)

 allocate(bielec_tmp_0(ao_num,ao_num),matrix_tmp_1(n_act_orb,ao_num),matrix_final(n_act_orb,n_act_orb))

  !$OMP DO SCHEDULE(guided)
 do p = 1,ao_num
   do q = 1,ao_num
    if(ao_overlap_abs(p,q).le.thr)cycle
    if(ao_bielec_integral_schwartz(p,q).lt.thr) cycle
    matrix_tmp_1 = 0.d0
    matrix_final = 0.d0
    do m = 1, ao_num
     call get_ao_bielec_integrals(p,m,q,ao_num,bielec_tmp_0(1,m)) ! k,l :: r1, m :: r2
    enddo
    call dgemm('N','N',n_act_orb,ao_num,ao_num,1.d0,mo_coef_act_transp,n_act_orb,bielec_tmp_0,ao_num,0.d0,matrix_tmp_1,n_act_orb)
    call dgemm('N','N',n_act_orb,n_act_orb,ao_num,1.d0,matrix_tmp_1,n_act_orb,mo_coef_act,ao_num,0.d0,matrix_final,n_act_orb)

    semi_transformed_act_act(1:n_act_orb,1:n_act_orb,p,q) = matrix_final(1:n_act_orb,1:n_act_orb)
   enddo
  enddo
  !$OMP END PARALLEL

 call cpu_time(cpu1) 
 print*, 'Time to do the semi transformation  act - act  : ',cpu1-cpu0
 
END_PROVIDER 








 BEGIN_PROVIDER [real(integral_kind), transformed_occ_virt_old, (n_core_inact_orb,n_virt_orb,n_core_inact_orb,n_virt_orb)]
  use map_module
 implicit none
 transformed_occ_virt_old = 0.d0

 double precision :: cpu0,cpu1
 integer :: i,j,k,l,iorb,jorb,korb,lorb,m,n,p,q
 double precision :: c
 real(integral_kind), allocatable :: bielec_tmp_0(:,:),matrix_tmp_1(:,:),matrix_final(:,:),matrix_tmp_2(:,:)
 real(integral_kind) :: integral,ao_bielec_integral
 real(integral_kind), allocatable :: transformed_occ_virt_old_tmp(:,:,:,:)
 real(integral_kind) :: thr
 
 provide ao_bielec_integral_schwartz ao_bielec_integrals_in_map ao_overlap_abs

 call cpu_time(cpu0) 
  !$OMP PARALLEL DEFAULT(NONE)             &
  !$OMP PRIVATE(p,q,matrix_tmp_1,matrix_final,m,bielec_tmp_0,matrix_tmp_2,i,j,transformed_occ_virt_old_tmp,thr) & 
  !$OMP SHARED(ao_num,ao_overlap_abs,ao_integrals_threshold,ao_bielec_integral_schwartz,n_core_inact_orb,mo_coef_core_inact_transp,mo_coef_virt, &
  !$OMP         mo_coef_virt_transp,n_virt_orb,transformed_occ_virt_old)
 
  thr = 0.001d0 * dsqrt(ao_integrals_threshold)
  allocate(bielec_tmp_0(ao_num,ao_num),matrix_tmp_1(n_core_inact_orb,ao_num),matrix_final(n_core_inact_orb,n_virt_orb))
  allocate(matrix_tmp_2(n_core_inact_orb,n_virt_orb),transformed_occ_virt_old_tmp(n_core_inact_orb,n_virt_orb,n_core_inact_orb,n_virt_orb))
  !$OMP DO SCHEDULE(guided)
 do p = 1,ao_num
   do q = 1,ao_num
    if(ao_overlap_abs(p,q).le.thr)cycle
    if(ao_bielec_integral_schwartz(p,q).lt.thr) cycle
    matrix_tmp_1 = 0.d0
    matrix_final = 0.d0
    do m = 1, ao_num
     call get_ao_bielec_integrals(p,m,q,ao_num,bielec_tmp_0(1,m)) ! k,l :: r1, m :: r2
    enddo
    call dgemm('N','N',n_core_inact_orb,ao_num,ao_num,1.d0,mo_coef_core_inact_transp,n_core_inact_orb,bielec_tmp_0,ao_num,0.d0,matrix_tmp_1,n_core_inact_orb)
    call dgemm('N','N',n_core_inact_orb,n_virt_orb,ao_num,1.d0,matrix_tmp_1,n_core_inact_orb,mo_coef_virt,ao_num,0.d0,matrix_final,n_core_inact_orb)

    matrix_tmp_2 = 0.d0
    call dger(n_core_inact_orb,n_virt_orb,1.d0,mo_coef_core_inact_transp(1,p),1,mo_coef_virt_transp(1,q),1,matrix_tmp_2,n_core_inact_orb)

    do j =1, n_virt_orb
     do i = 1, n_core_inact_orb
       if(dabs(matrix_final(i,j)).lt.thr)cycle
       transformed_occ_virt_old_tmp(1:n_core_inact_orb,1:n_virt_orb,i,j) += matrix_tmp_2(1:n_core_inact_orb,1:n_virt_orb) * matrix_final(i,j)
     enddo
    enddo
   enddo
  enddo
  !$OMP END DO NOWAIT
  !$OMP CRITICAL
  transformed_occ_virt_old = transformed_occ_virt_old + transformed_occ_virt_old_tmp
  !$OMP END CRITICAL
  deallocate(transformed_occ_virt_old_tmp)
  !$OMP END PARALLEL
 
 call cpu_time(cpu1) 
 print*, 'Time to transform occ-virt bielec =',cpu1-cpu0
 
END_PROVIDER 

 BEGIN_PROVIDER [double precision, mo_coef_core_inact, (ao_num, n_core_inact_orb)]
&BEGIN_PROVIDER [double precision, mo_coef_core_inact_transp, (n_core_inact_orb,ao_num)]
 implicit none
 integer :: i,j,k,iorb,jorb
 do i = 1, n_core_inact_orb
  iorb = list_core_inact(i) 
  do k = 1, ao_num
   mo_coef_core_inact(k,i) = mo_coef(k,iorb)
   mo_coef_core_inact_transp(i,k) = mo_coef(k,iorb)
  enddo
 enddo

END_PROVIDER 

 BEGIN_PROVIDER [double precision, mo_coef_virt, (ao_num, n_virt_orb)]
&BEGIN_PROVIDER [double precision, mo_coef_virt_transp, (n_virt_orb,ao_num)]
 implicit none
 integer :: i,j,k,iorb,jorb
 do i = 1, n_virt_orb
  iorb = list_virt(i) 
  do k = 1, ao_num
   mo_coef_virt(k,i) = mo_coef(k,iorb)
   mo_coef_virt_transp(i,k) = mo_coef(k,iorb)
  enddo
 enddo

END_PROVIDER 


 BEGIN_PROVIDER [double precision, mo_coef_act, (ao_num, n_act_orb)]
&BEGIN_PROVIDER [double precision, mo_coef_act_transp, (n_act_orb,ao_num)]
 implicit none
 integer :: i,j,k,iorb,jorb
 do i = 1, n_act_orb
  iorb = list_act(i) 
  do k = 1, ao_num
   mo_coef_act(k,i) = mo_coef(k,iorb)
   mo_coef_act_transp(i,k) = mo_coef(k,iorb)
  enddo
 enddo

END_PROVIDER 


subroutine get_all_core_inact_virt_integrals(iorb,vorb,matrix_integrals)
  use map_module
 implicit none
 integer, intent(in) :: iorb,vorb
 real(integral_kind), intent(out) :: matrix_integrals(n_core_inact_orb,n_virt_orb)
 real(integral_kind) :: integrals_matrix(ao_num,ao_num)
 real(integral_kind) :: matrix_tmp_1(n_core_inact_orb,ao_num)
 integer :: i,j
 integer :: m,n

 do m = 1, ao_num
  do n = 1, ao_num
   integrals_matrix(m,n) = semi_transformed_occ_virt(m,n,iorb,vorb)
  enddo
 enddo
 matrix_integrals = 0.d0
 
 call dgemm('N','N',n_core_inact_orb,ao_num,ao_num,1.d0,mo_coef_core_inact_transp,n_core_inact_orb,integrals_matrix,ao_num,0.d0,matrix_tmp_1,n_core_inact_orb)
 call dgemm('N','N',n_core_inact_orb,n_virt_orb,ao_num,1.d0,matrix_tmp_1,n_core_inact_orb,mo_coef_virt,ao_num,0.d0,matrix_integrals,n_core_inact_orb)

end
