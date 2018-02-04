subroutine mo_bielec_integrals_erf_index(i,j,k,l,i1)
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


BEGIN_PROVIDER [ logical, mo_bielec_integrals_erf_in_map ]
  use map_module
  implicit none
  integer(bit_kind)              :: mask_ijkl(N_int,4)
  integer(bit_kind)              :: mask_ijk(N_int,3)
  
  BEGIN_DOC
  ! If True, the map of MO bielectronic integrals is provided
  END_DOC
  
  real                           :: map_mb

  mo_bielec_integrals_erf_in_map = .True.
  if (read_mo_integrals_erf) then
    print*,'Reading the MO integrals_erf'
    call map_load_from_disk(trim(ezfio_filename)//'/work/mo_ints_erf',mo_integrals_erf_map)
    print*, 'MO integrals_erf provided'
    return
  else
    PROVIDE ao_bielec_integrals_erf_in_map
  endif
  
     call four_index_transform_block(ao_integrals_erf_map,mo_integrals_erf_map, &
         mo_coef, size(mo_coef,1),                                      &
         1, 1, 1, 1, ao_num, ao_num, ao_num, ao_num,                    &
         1, 1, 1, 1, mo_num, mo_num, mo_num, mo_num)
    integer*8                      :: get_mo_erf_map_size, mo_erf_map_size
    mo_erf_map_size = get_mo_erf_map_size()
    
  print*,'Molecular integrals ERF provided:'
  print*,' Size of MO ERF map           ', map_mb(mo_integrals_erf_map) ,'MB'
  print*,' Number of MO ERF integrals: ',  mo_erf_map_size
  if (write_mo_integrals_erf) then
    call ezfio_set_work_empty(.False.)
    call map_save_to_disk(trim(ezfio_filename)//'/work/mo_ints_erf',mo_integrals_erf_map)
    call ezfio_set_integrals_erf_disk_access_mo_integrals_erf("Read")
  endif
  
END_PROVIDER


 BEGIN_PROVIDER [ double precision, mo_bielec_integral_erf_jj_from_ao, (mo_tot_num,mo_tot_num) ]
&BEGIN_PROVIDER [ double precision, mo_bielec_integral_erf_jj_exchange_from_ao, (mo_tot_num,mo_tot_num) ]
&BEGIN_PROVIDER [ double precision, mo_bielec_integral_erf_jj_anti_from_ao, (mo_tot_num,mo_tot_num) ]
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
    PROVIDE ao_bielec_integrals_erf_in_map mo_coef
  endif
  
  mo_bielec_integral_erf_jj_from_ao = 0.d0
  mo_bielec_integral_erf_jj_exchange_from_ao = 0.d0
  
  !DIR$ ATTRIBUTES ALIGN : $IRP_ALIGN :: iqrs, iqsr
  
  
  !$OMP PARALLEL DEFAULT(NONE)                                       &
      !$OMP PRIVATE (i,j,p,q,r,s,integral,c,n,pp,int_value,int_idx,  &
      !$OMP  iqrs, iqsr,iqri,iqis)                                   &
      !$OMP SHARED(mo_tot_num,mo_coef_transp,ao_num,&
      !$OMP  ao_integrals_threshold,do_direct_integrals)             &
      !$OMP REDUCTION(+:mo_bielec_integral_erf_jj_from_ao,mo_bielec_integral_erf_jj_exchange_from_ao)
  
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
      
      if (do_direct_integrals) then
        double precision               :: ao_bielec_integral_erf
        do r=1,ao_num
          call compute_ao_bielec_integrals_erf(q,r,s,ao_num,int_value)
          do p=1,ao_num
            integral = int_value(p)
            if (abs(integral) > ao_integrals_threshold) then
              !DIR$ VECTOR ALIGNED
              do i=1,mo_tot_num
                iqrs(i,r) += mo_coef_transp(i,p) * integral
              enddo
            endif
          enddo
          call compute_ao_bielec_integrals_erf(q,s,r,ao_num,int_value)
          do p=1,ao_num
            integral = int_value(p)
            if (abs(integral) > ao_integrals_threshold) then
              !DIR$ VECTOR ALIGNED
              do i=1,mo_tot_num
                iqsr(i,r) += mo_coef_transp(i,p) * integral
              enddo
            endif
          enddo
        enddo
        
      else
        
        do r=1,ao_num
          call get_ao_bielec_integrals_erf_non_zero(q,r,s,ao_num,int_value,int_idx,n)
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
          call get_ao_bielec_integrals_erf_non_zero(q,s,r,ao_num,int_value,int_idx,n)
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
      endif
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
          mo_bielec_integral_erf_jj_from_ao(j,i) += c * iqis(i)
          mo_bielec_integral_erf_jj_exchange_from_ao(j,i) += c * iqri(i)
        enddo
      enddo
      
    enddo
  enddo
  !$OMP END DO NOWAIT
  deallocate(iqrs,iqsr,int_value,int_idx)
  !$OMP END PARALLEL
  
  mo_bielec_integral_erf_jj_anti_from_ao = mo_bielec_integral_erf_jj_from_ao - mo_bielec_integral_erf_jj_exchange_from_ao
  
  
! end
END_PROVIDER 


 BEGIN_PROVIDER [ double precision, mo_bielec_integral_erf_jj, (mo_tot_num,mo_tot_num) ]
&BEGIN_PROVIDER [ double precision, mo_bielec_integral_erf_jj_exchange, (mo_tot_num,mo_tot_num) ]
&BEGIN_PROVIDER [ double precision, mo_bielec_integral_erf_jj_anti, (mo_tot_num,mo_tot_num) ]
  implicit none
  BEGIN_DOC
  ! mo_bielec_integral_jj(i,j) = J_ij
  ! mo_bielec_integral_jj_exchange(i,j) = K_ij
  ! mo_bielec_integral_jj_anti(i,j) = J_ij - K_ij
  END_DOC
  
  integer                        :: i,j
  double precision               :: get_mo_bielec_integral_erf
  
  PROVIDE mo_bielec_integrals_erf_in_map
  mo_bielec_integral_erf_jj = 0.d0
  mo_bielec_integral_erf_jj_exchange = 0.d0
  
  do j=1,mo_tot_num
    do i=1,mo_tot_num
      mo_bielec_integral_erf_jj(i,j) = get_mo_bielec_integral_erf(i,j,i,j,mo_integrals_erf_map)
      mo_bielec_integral_erf_jj_exchange(i,j) = get_mo_bielec_integral_erf(i,j,j,i,mo_integrals_erf_map)
      mo_bielec_integral_erf_jj_anti(i,j) = mo_bielec_integral_erf_jj(i,j) - mo_bielec_integral_erf_jj_exchange(i,j)
    enddo
  enddo
  
END_PROVIDER


subroutine clear_mo_erf_map
  implicit none
  BEGIN_DOC
  ! Frees the memory of the MO map
  END_DOC
  call map_deinit(mo_integrals_erf_map)
  FREE mo_integrals_erf_map mo_bielec_integral_erf_jj mo_bielec_integral_erf_jj_anti
  FREE mo_bielec_integral_Erf_jj_exchange mo_bielec_integrals_erf_in_map
  
  
end

subroutine provide_all_mo_integrals_erf
  implicit none
  provide mo_integrals_erf_map mo_bielec_integral_erf_jj mo_bielec_integral_erf_jj_anti
  provide mo_bielec_integral_erf_jj_exchange mo_bielec_integrals_erf_in_map
  
end
