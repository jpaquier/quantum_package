subroutine mo_bielec_integrals_erf_mu_of_r_index(i,j,k,l,i1)
  use map_module
  implicit none
  BEGIN_DOC
  ! Computes an unique index for i,j,k,l integrals
  END_DOC
  integer, intent(in)            :: i,j,k,l
  integer(key_kind), intent(out) :: i1
  integer(key_kind)              :: p,q,r,s,i2
  p = min(i,k)
  r = max(i,k)
  p = p+ishft(r*r-r,-1)
  q = min(j,l)
  s = max(j,l)
  q = q+ishft(s*s-s,-1)
  i1 = min(p,q)
  i2 = max(p,q)
  i1 = i1+ishft(i2*i2-i2,-1)
end


BEGIN_PROVIDER [ logical, mo_bielec_integrals_erf_mu_of_r_in_map ]
  use map_module
  implicit none
  integer(bit_kind)              :: mask_ijkl(N_int,4)
  integer(bit_kind)              :: mask_ijk(N_int,3)
  
  BEGIN_DOC
  ! If True, the map of MO bielectronic integrals is provided
  END_DOC
  
  real                           :: map_mb

  mo_bielec_integrals_erf_mu_of_r_in_map = .True.
  if (read_mo_integrals_mu_of_r) then
    print*,'Reading the MO integrals_erf_mu_of_r'
    call map_load_from_disk(trim(ezfio_filename)//'/work/mo_ints_erf_mu_of_r',mo_integrals_erf_mu_of_r_map)
    print*, 'MO integrals_erf_mu_of_r provided'
    return
  else
    PROVIDE ao_bielec_integrals_erf_mu_of_r_in_map
  endif
  
  !  call four_index_transform_block(ao_integrals_erf_mu_of_r_map,mo_integrals_erf_mu_of_r_map, &
  !      mo_coef, size(mo_coef,1),                                      &
  !      1, 1, 1, 1, ao_num, ao_num, ao_num, ao_num,                    &
  !      1, 1, 1, 1, mo_num, mo_num, mo_num, mo_num)
   call add_integrals_to_map_erf_mu_of_r(full_ijkl_bitmask_4)
    integer*8                      :: get_mo_erf_mu_of_r_map_size, mo_erf_mu_of_r_map_size
    mo_erf_mu_of_r_map_size = get_mo_erf_mu_of_r_map_size()
    
! print*,'Molecular integrals ERF provided:'
! print*,' Size of MO ERF map           ', map_mb(mo_integrals_erf_mu_of_r_map) ,'MB'
! print*,' Number of MO ERF integrals: ',  mo_erf_mu_of_r_map_size
  if (write_mo_integrals_mu_of_r) then
    call ezfio_set_work_empty(.False.)
    call map_save_to_disk(trim(ezfio_filename)//'/work/mo_ints_erf_mu_of_r',mo_integrals_erf_mu_of_r_map)
    call ezfio_set_mu_of_r_ints_disk_mo_integrals_mu_of_r("Read")
  endif
  
END_PROVIDER


 BEGIN_PROVIDER [ double precision, mo_bielec_integral_erf_mu_of_r_jj_from_ao, (mo_tot_num,mo_tot_num) ]
&BEGIN_PROVIDER [ double precision, mo_bielec_integral_erf_mu_of_r_jj_exchange_from_ao, (mo_tot_num,mo_tot_num) ]
&BEGIN_PROVIDER [ double precision, mo_bielec_integral_erf_mu_of_r_jj_anti_from_ao, (mo_tot_num,mo_tot_num) ]
  BEGIN_DOC
  ! mo_bielec_integral_jj_from_ao(i,j) = J_ij
  ! mo_bielec_integral_jj_exchange_from_ao(i,j) = J_ij
  ! mo_bielec_integral_jj_anti_from_ao(i,j) = J_ij - K_ij
  END_DOC
  implicit none
  integer                        :: i,j,p,q,r,s
  double precision               :: c
  real(integral_kind)            :: integral
  integer                        :: n, pp
  real(integral_kind), allocatable :: int_value(:)
  integer, allocatable           :: int_idx(:)
  
  double precision, allocatable  :: iqrs(:,:), iqsr(:,:), iqis(:), iqri(:)
  
  if (.not.do_direct_integrals) then
    PROVIDE ao_bielec_integrals_erf_mu_of_r_in_map mo_coef
  endif
  
  mo_bielec_integral_erf_mu_of_r_jj_from_ao = 0.d0
  mo_bielec_integral_erf_mu_of_r_jj_exchange_from_ao = 0.d0
  
  !DIR$ ATTRIBUTES ALIGN : $IRP_ALIGN :: iqrs, iqsr
  
  
  !$OMP PARALLEL DEFAULT(NONE)                                       &
      !$OMP PRIVATE (i,j,p,q,r,s,integral,c,n,pp,int_value,int_idx,  &
      !$OMP  iqrs, iqsr,iqri,iqis)                                   &
      !$OMP SHARED(mo_tot_num,mo_coef_transp,ao_num,&
      !$OMP  ao_integrals_threshold,do_direct_integrals)             &
      !$OMP REDUCTION(+:mo_bielec_integral_erf_mu_of_r_jj_from_ao,mo_bielec_integral_erf_mu_of_r_jj_exchange_from_ao)
  
  allocate( int_value(ao_num), int_idx(ao_num),                      &
      iqrs(mo_tot_num,ao_num), iqis(mo_tot_num), iqri(mo_tot_num),&
      iqsr(mo_tot_num,ao_num) )
  
  !$OMP DO SCHEDULE (guided)
  do s=1,ao_num
    do q=1,ao_num
      
      do j=1,ao_num
        !DIR$ VECTOR ALIGNED
        do i=1,mo_tot_num
          iqrs(i,j) = 0.d0
          iqsr(i,j) = 0.d0
        enddo
      enddo
      
        do r=1,ao_num
          call get_ao_bielec_integrals_erf_mu_of_r_non_zero(q,r,s,ao_num,int_value,int_idx,n)
          do pp=1,n
            p = int_idx(pp)
            integral = int_value(pp)
            if (abs(integral) > ao_integrals_threshold) then
              !DIR$ VECTOR ALIGNED
              do i=1,mo_tot_num
                iqrs(i,r) += mo_coef_transp(i,p) * integral
              enddo
            endif
          enddo
          call get_ao_bielec_integrals_erf_mu_of_r_non_zero(q,s,r,ao_num,int_value,int_idx,n)
          do pp=1,n
            p = int_idx(pp)
            integral = int_value(pp)
            if (abs(integral) > ao_integrals_threshold) then
              !DIR$ VECTOR ALIGNED
              do i=1,mo_tot_num
                iqsr(i,r) += mo_coef_transp(i,p) * integral
              enddo
            endif
          enddo
        enddo
!     endif
      iqis = 0.d0
      iqri = 0.d0
      do r=1,ao_num
        !DIR$ VECTOR ALIGNED
        do i=1,mo_tot_num
          iqis(i) += mo_coef_transp(i,r) * iqrs(i,r)
          iqri(i) += mo_coef_transp(i,r) * iqsr(i,r)
        enddo
      enddo
      do i=1,mo_tot_num
        !DIR$ VECTOR ALIGNED
        do j=1,mo_tot_num
          c = mo_coef_transp(j,q)*mo_coef_transp(j,s)
          mo_bielec_integral_erf_mu_of_r_jj_from_ao(j,i) += c * iqis(i)
          mo_bielec_integral_erf_mu_of_r_jj_exchange_from_ao(j,i) += c * iqri(i)
        enddo
      enddo
      
    enddo
  enddo
  !$OMP END DO NOWAIT
  deallocate(iqrs,iqsr,int_value,int_idx)
  !$OMP END PARALLEL
  
  mo_bielec_integral_erf_mu_of_r_jj_anti_from_ao = mo_bielec_integral_erf_mu_of_r_jj_from_ao - mo_bielec_integral_erf_mu_of_r_jj_exchange_from_ao
  
  
! end
END_PROVIDER 


 BEGIN_PROVIDER [ double precision, mo_bielec_integral_erf_mu_of_r_jj, (mo_tot_num,mo_tot_num) ]
&BEGIN_PROVIDER [ double precision, mo_bielec_integral_erf_mu_of_r_jj_exchange, (mo_tot_num,mo_tot_num) ]
&BEGIN_PROVIDER [ double precision, mo_bielec_integral_erf_mu_of_r_jj_anti, (mo_tot_num,mo_tot_num) ]
  implicit none
  BEGIN_DOC
  ! mo_bielec_integral_jj(i,j) = J_ij
  ! mo_bielec_integral_jj_exchange(i,j) = K_ij
  ! mo_bielec_integral_jj_anti(i,j) = J_ij - K_ij
  END_DOC
  
  integer                        :: i,j
  double precision               :: get_mo_bielec_integral_erf_mu_of_r
  
  PROVIDE mo_bielec_integrals_erf_mu_of_r_in_map
  mo_bielec_integral_erf_mu_of_r_jj = 0.d0
  mo_bielec_integral_erf_mu_of_r_jj_exchange = 0.d0
  
  do j=1,mo_tot_num
    do i=1,mo_tot_num
      mo_bielec_integral_erf_mu_of_r_jj(i,j) = get_mo_bielec_integral_erf_mu_of_r(i,j,i,j,mo_integrals_erf_mu_of_r_map)
      mo_bielec_integral_erf_mu_of_r_jj_exchange(i,j) = get_mo_bielec_integral_erf_mu_of_r(i,j,j,i,mo_integrals_erf_mu_of_r_map)
      mo_bielec_integral_erf_mu_of_r_jj_anti(i,j) = mo_bielec_integral_erf_mu_of_r_jj(i,j) - mo_bielec_integral_erf_mu_of_r_jj_exchange(i,j)
    enddo
  enddo
  
END_PROVIDER


subroutine clear_mo_erf_mu_of_r_map
  implicit none
  BEGIN_DOC
  ! Frees the memory of the MO map
  END_DOC
  call map_deinit(mo_integrals_erf_mu_of_r_map)
  FREE mo_integrals_erf_mu_of_r_map mo_bielec_integral_erf_mu_of_r_jj mo_bielec_integral_erf_mu_of_r_jj_anti
  FREE mo_bielec_integral_erf_mu_of_r_jj_exchange mo_bielec_integrals_erf_mu_of_r_in_map
  
  
end

subroutine provide_all_mo_integrals_erf_mu_of_r
  implicit none
  provide mo_integrals_erf_mu_of_r_map mo_bielec_integral_erf_mu_of_r_jj mo_bielec_integral_erf_mu_of_r_jj_anti
  provide mo_bielec_integral_erf_mu_of_r_jj_exchange mo_bielec_integrals_erf_mu_of_r_in_map
  
end


subroutine add_integrals_to_map_erf_mu_of_r(mask_ijkl)
  use bitmasks
  implicit none
  
  BEGIN_DOC
  ! Adds integrals to tha MO map according to some bitmask
  END_DOC
  
  integer(bit_kind), intent(in)  :: mask_ijkl(N_int,4)
  
  integer                        :: i,j,k,l
  integer                        :: i0,j0,k0,l0
  double precision               :: c, cpu_1, cpu_2, wall_1, wall_2, wall_0
  
  integer, allocatable           :: list_ijkl(:,:)
  integer                        :: n_i, n_j, n_k, n_l
  integer, allocatable           :: bielec_tmp_0_idx(:)
  real(integral_kind), allocatable :: bielec_tmp_0(:,:)
  double precision, allocatable  :: bielec_tmp_1(:)
  double precision, allocatable  :: bielec_tmp_2(:,:)
  double precision, allocatable  :: bielec_tmp_3(:,:,:)
  !DIR$ ATTRIBUTES ALIGN : 64    :: bielec_tmp_1, bielec_tmp_2, bielec_tmp_3
  
  integer                        :: n_integrals
  integer                        :: size_buffer
  integer(key_kind),allocatable  :: buffer_i(:)
  real(integral_kind),allocatable :: buffer_value(:)
  double precision               :: map_mb
  
  integer                        :: i1,j1,k1,l1, ii1, kmax, thread_num
  integer                        :: i2,i3,i4
  double precision,parameter     :: thr_coef = 1.d-10
  
  PROVIDE ao_bielec_integrals_in_map  mo_coef
  
  !Get list of MOs for i,j,k and l
  !-------------------------------
  
  allocate(list_ijkl(mo_tot_num,4))
  call bitstring_to_list( mask_ijkl(1,1), list_ijkl(1,1), n_i, N_int )
  call bitstring_to_list( mask_ijkl(1,2), list_ijkl(1,2), n_j, N_int )
  call bitstring_to_list( mask_ijkl(1,3), list_ijkl(1,3), n_k, N_int )
  call bitstring_to_list( mask_ijkl(1,4), list_ijkl(1,4), n_l, N_int )
  character*(2048)               :: output(1)
  print*, 'i'
  call bitstring_to_str( output(1), mask_ijkl(1,1), N_int )
  print *,  trim(output(1))
  j = 0
  do i = 1, N_int
    j += popcnt(mask_ijkl(i,1))
  enddo
  if(j==0)then
    return
  endif
  
  print*, 'j'
  call bitstring_to_str( output(1), mask_ijkl(1,2), N_int )
  print *,  trim(output(1))
  j = 0
  do i = 1, N_int
    j += popcnt(mask_ijkl(i,2))
  enddo
  if(j==0)then
    return
  endif
  
  print*, 'k'
  call bitstring_to_str( output(1), mask_ijkl(1,3), N_int )
  print *,  trim(output(1))
  j = 0
  do i = 1, N_int
    j += popcnt(mask_ijkl(i,3))
  enddo
  if(j==0)then
    return
  endif
  
  print*, 'l'
  call bitstring_to_str( output(1), mask_ijkl(1,4), N_int )
  print *,  trim(output(1))
  j = 0
  do i = 1, N_int
    j += popcnt(mask_ijkl(i,4))
  enddo
  if(j==0)then
    return
  endif
  
  size_buffer = min(ao_num*ao_num*ao_num,16000000)
  print*, 'Providing the ERF mu of r molecular integrals '
  print*, 'Buffers : ', 8.*(mo_tot_num*(n_j)*(n_k+1) + mo_tot_num+&
      ao_num+ao_num*ao_num+ size_buffer*3)/(1024*1024), 'MB / core'
  
  call wall_time(wall_1)
  call cpu_time(cpu_1)
  double precision               :: accu_bis
  accu_bis = 0.d0
  
  !$OMP PARALLEL PRIVATE(l1,k1,j1,i1,i2,i3,i4,i,j,k,l,c, ii1,kmax,   &
      !$OMP  bielec_tmp_0_idx, bielec_tmp_0, bielec_tmp_1,bielec_tmp_2,bielec_tmp_3,&
      !$OMP  buffer_i,buffer_value,n_integrals,wall_2,i0,j0,k0,l0,   &
      !$OMP  wall_0,thread_num,accu_bis)                             &
      !$OMP  DEFAULT(NONE)                                           &
      !$OMP  SHARED(size_buffer,ao_num,mo_tot_num,n_i,n_j,n_k,n_l,   &
      !$OMP  mo_coef_transp,                                         &
      !$OMP  mo_coef_transp_is_built, list_ijkl,                     &
      !$OMP  mo_coef_is_built, wall_1,                               &
      !$OMP  mo_coef,mo_integrals_threshold,mo_integrals_erf_mu_of_r_map)
  n_integrals = 0
  wall_0 = wall_1
  allocate(bielec_tmp_3(mo_tot_num, n_j, n_k),                 &
      bielec_tmp_1(mo_tot_num),                                &
      bielec_tmp_0(ao_num,ao_num),                                   &
      bielec_tmp_0_idx(ao_num),                                      &
      bielec_tmp_2(mo_tot_num, n_j),                           &
      buffer_i(size_buffer),                                         &
      buffer_value(size_buffer) )
  
  thread_num = 0
  !$  thread_num = omp_get_thread_num()
  !$OMP DO SCHEDULE(guided)
  do l1 = 1,ao_num
    bielec_tmp_3 = 0.d0
    do k1 = 1,ao_num
      bielec_tmp_2 = 0.d0
      do j1 = 1,ao_num
        call get_ao_bielec_integrals_erf_mu_of_r(j1,k1,l1,ao_num,bielec_tmp_0(1,j1)) ! all integrals for a given l1, k1
        ! call compute_ao_bielec_integrals(j1,k1,l1,ao_num,bielec_tmp_0(1,j1))
      enddo
      do j1 = 1,ao_num
        kmax = 0
        do i1 = 1,ao_num
          c = bielec_tmp_0(i1,j1)
          if (c == 0.d0) then
            cycle
          endif
          kmax += 1
          bielec_tmp_0(kmax,j1) = c
          bielec_tmp_0_idx(kmax) = i1
        enddo
        
        if (kmax==0) then
          cycle
        endif
        
        bielec_tmp_1 = 0.d0
        ii1=1
        ! sum_m c_m^i (m)
        do ii1 = 1,kmax-4,4
          i1 = bielec_tmp_0_idx(ii1)
          i2 = bielec_tmp_0_idx(ii1+1)
          i3 = bielec_tmp_0_idx(ii1+2)
          i4 = bielec_tmp_0_idx(ii1+3)
          do i = list_ijkl(1,1), list_ijkl(n_i,1)
            bielec_tmp_1(i)  =  bielec_tmp_1(i) +                    &
                mo_coef_transp(i,i1) * bielec_tmp_0(ii1,j1) +        &
                mo_coef_transp(i,i2) * bielec_tmp_0(ii1+1,j1) +      &
                mo_coef_transp(i,i3) * bielec_tmp_0(ii1+2,j1) +      &
                mo_coef_transp(i,i4) * bielec_tmp_0(ii1+3,j1)
          enddo ! i
        enddo  ! ii1
        
        i2 = ii1
        do ii1 = i2,kmax
          i1 = bielec_tmp_0_idx(ii1)
          do i = list_ijkl(1,1), list_ijkl(n_i,1)
            bielec_tmp_1(i) = bielec_tmp_1(i) + mo_coef_transp(i,i1) * bielec_tmp_0(ii1,j1)
          enddo ! i
        enddo  ! ii1
        c = 0.d0
        
        do i = list_ijkl(1,1), list_ijkl(n_i,1)
          c = max(c,abs(bielec_tmp_1(i)))
          if (c>mo_integrals_threshold) exit
        enddo
        if ( c < mo_integrals_threshold ) then
          cycle
        endif
        
        do j0 = 1, n_j
          j = list_ijkl(j0,2)
          c = mo_coef_transp(j,j1)
          if (abs(c) < thr_coef) then
            cycle
          endif
          do i = list_ijkl(1,1), list_ijkl(n_i,1)
            bielec_tmp_2(i,j0)  = bielec_tmp_2(i,j0) + c * bielec_tmp_1(i)
          enddo ! i
        enddo  ! j
      enddo !j1
      if ( maxval(abs(bielec_tmp_2)) < mo_integrals_threshold ) then
        cycle
      endif
      
      
      do k0 = 1, n_k
        k = list_ijkl(k0,3)
        c = mo_coef_transp(k,k1)
        if (abs(c) < thr_coef) then
          cycle
        endif
        
        do j0 = 1, n_j
          j = list_ijkl(j0,2)
          do i = list_ijkl(1,1), k
            bielec_tmp_3(i,j0,k0) = bielec_tmp_3(i,j0,k0) + c* bielec_tmp_2(i,j0)
          enddo!i
        enddo !j
        
      enddo  !k
    enddo   !k1
    
    
    
    do l0 = 1,n_l
      l = list_ijkl(l0,4)
      c = mo_coef_transp(l,l1)
      if (abs(c) < thr_coef) then
        cycle
      endif
      j1 = ishft((l*l-l),-1)
      do j0 = 1, n_j
        j = list_ijkl(j0,2)
        if (j > l)  then
          exit
        endif
        j1 += 1
        do k0 = 1, n_k
          k = list_ijkl(k0,3)
          i1 = ishft((k*k-k),-1)
          if (i1<=j1) then
            continue
          else
            exit
          endif
          bielec_tmp_1 = 0.d0
          do i0 = 1, n_i
            i = list_ijkl(i0,1)
            if (i>k) then
              exit
            endif
            bielec_tmp_1(i) = c*bielec_tmp_3(i,j0,k0)
            !           i1+=1
          enddo
          
          do i0 = 1, n_i
            i = list_ijkl(i0,1)
            if(i> min(k,j1-i1+list_ijkl(1,1)-1))then
              exit
            endif
            if (abs(bielec_tmp_1(i)) < mo_integrals_threshold) then
              cycle
            endif
            n_integrals += 1
            buffer_value(n_integrals) = bielec_tmp_1(i)
            !DIR$ FORCEINLINE
            call mo_bielec_integrals_index(i,j,k,l,buffer_i(n_integrals))
            if (n_integrals == size_buffer) then
              call insert_into_mo_integrals_erf_mu_of_r_map(n_integrals,buffer_i,buffer_value,&
                  real(mo_integrals_threshold,integral_kind))
              n_integrals = 0
            endif
          enddo
        enddo
      enddo
    enddo
    
    call wall_time(wall_2)
    if (thread_num == 0) then
      if (wall_2 - wall_0 > 1.d0) then
        wall_0 = wall_2
        print*, 100.*float(l1)/float(ao_num), '% in ',               &
            wall_2-wall_1, 's', map_mb(mo_integrals_erf_mu_of_r_map) ,'MB'
      endif
    endif
  enddo
  !$OMP END DO NOWAIT
  deallocate (bielec_tmp_1,bielec_tmp_2,bielec_tmp_3)
  
  integer                        :: index_needed
  
  call insert_into_mo_integrals_erf_mu_of_r_map(n_integrals,buffer_i,buffer_value,&
      real(mo_integrals_threshold,integral_kind))
  deallocate(buffer_i, buffer_value)
  !$OMP END PARALLEL
  call map_merge(mo_integrals_erf_mu_of_r_map)
  
  call wall_time(wall_2)
  call cpu_time(cpu_2)
  integer*8                      :: get_mo_erf_mu_of_r_map_size, mo_map_size
  mo_map_size = get_mo_erf_mu_of_r_map_size()
  
  deallocate(list_ijkl)
  
  
  print*,'Molecular ERF mu of r integrals provided:'
  print*,' Size of MO ERF mu of r map           ', map_mb(mo_integrals_erf_mu_of_r_map) ,'MB'
  print*,' Number of MO ERF mu of r integrals: ',  mo_map_size
  print*,' cpu  time :',cpu_2 - cpu_1, 's'
  print*,' wall time :',wall_2 - wall_1, 's  ( x ', (cpu_2-cpu_1)/(wall_2-wall_1), ')'
  
end


