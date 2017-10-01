BEGIN_PROVIDER [double precision, core_energy_from_core]
 implicit none
 BEGIN_DOC
! energy of the block of doubly occupied electrons NOT TAKING INTO ACCOUNT THE ACTIVE ELECTRONS 
 END_DOC
 core_energy_from_core = nuclear_repulsion
 integer :: i,j

 do j=1,ao_num
   do i=1,ao_num
     core_energy_from_core += 0.5d0 * (                                          &
         (ao_mono_elec_integral(i,j) + Fock_matrix_alpha_ao(i,j) ) *  HF_density_matrix_ao_alpha_core_inact(i,j) +&
         (ao_mono_elec_integral(i,j) + Fock_matrix_beta_ao (i,j) ) *  HF_density_matrix_ao_beta_core_inact(i,j) )
   enddo
 enddo

END_PROVIDER 


 BEGIN_PROVIDER [double precision, coulomb_core_energy_from_act]
&BEGIN_PROVIDER [double precision, alpha_exch_core_energy_from_act]
&BEGIN_PROVIDER [double precision, beta_exch_core_energy_from_act]
 implicit none
 BEGIN_DOC
! coulomb_core_energy_from_act = \sum_{m,n,p,q} (mn|pq) * HF_density_matrix_ao_alpha_core_inact(m,n) *2 * HF_density_matrix_ao_act(p,q)
! alpha_exch_core_energy_from_act = \sum_{m,n,p,q} (mn|pq) * HF_density_matrix_ao_alpha_core_inact(m,p)  * HF_density_matrix_ao_act_alpha(n,q)
! beta_exch_core_energy_from_act = \sum_{m,n,p,q} (mn|pq) * HF_density_matrix_ao_beta_core_inact(m,p)  * HF_density_matrix_ao_act_beta(n,q)
 END_DOC

 integer                        :: i,j,k,l,k1
 double precision               :: integral

   PROVIDE ao_bielec_integrals_in_map 
           
   integer*8                      :: i8
   integer                        :: ii(8), jj(8), kk(8), ll(8), k2
   integer(cache_map_size_kind)   :: n_elements_max, n_elements
   integer(key_kind), allocatable :: keys(:)
   double precision, allocatable  :: values(:)

   call get_cache_map_n_elements_max(ao_integrals_map,n_elements_max)
   allocate(keys(n_elements_max), values(n_elements_max))

  !!$OMP PARALLEL DEFAULT(NONE)                                      &
  !!OMP PRIVATE (i8,n_elements,keys,values,kk,ii,ll,jj,i,j,k,l,k2,k1) & 
  !!OMP SHARED (ao_integrals_map,HF_density_matrix_ao_alpha_core_inact,HF_density_matrix_ao_alpha_core_inact,HF_density_matrix_ao_beta_core_inact, & 
  !!OMP SHARED  HF_density_matrix_ao_act,HF_density_matrix_ao_act_alpha,HF_density_matrix_ao_act_beta,n_elements_max) &
  !!$OMP REDUCTION (+:coulomb_core_energy_from_act)        &
  !!$OMP REDUCTION (+:alpha_exch_core_energy_from_act)      &  
  !!$OMP REDUCTION (+:beta_exch_core_energy_from_act)       
  !!$OMP DO SCHEDULE(dynamic) &
   do i8=0_8,ao_integrals_map%map_size
     n_elements = n_elements_max
     call get_cache_map(ao_integrals_map,i8,keys,values,n_elements)
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
         coulomb_core_energy_from_act += HF_density_matrix_ao_alpha_core_inact(k,l) * 2.d0 * HF_density_matrix_ao_act(i,j) * values(k1)
         alpha_exch_core_energy_from_act -= HF_density_matrix_ao_alpha_core_inact(i,k) * HF_density_matrix_ao_act_alpha(j,l) * values(k1)
         beta_exch_core_energy_from_act -= HF_density_matrix_ao_beta_core_inact(i,k) * HF_density_matrix_ao_act_beta(j,l) * values(k1)
       enddo
     enddo
   enddo
  !!$OMP END DO NOWAIT
   deallocate(keys,values)
  !!$OMP END PARALLEL
END_PROVIDER 


BEGIN_PROVIDER [double precision, ref_energy_core_core_and_core_act]
 implicit none
 ref_energy_core_core_and_core_act = core_energy_from_core + coulomb_core_energy_from_act +  alpha_exch_core_energy_from_act + beta_exch_core_energy_from_act


END_PROVIDER 

BEGIN_PROVIDER [double precision, ref_energy_act_act]
 implicit none
 integer :: i,j,k,l,h1,h2,p1,p2
 double precision  :: get_mo_bielec_integral
 ref_energy_act_act = 0.d0
 ! purely active part of the two body dm 
 do l = 1, n_act_orb  ! p2 
  p2 = list_act(l)
  do k = 1, n_act_orb  ! h2 
   h2 = list_act(k)
   do j = 1, n_act_orb  ! p1 
    p1 = list_act(j)
    do i = 1,n_act_orb   ! h1 
     h1 = list_act(i)
     if (dabs(transformed_act1_act1_act2_act2(i,j,k,l) - get_mo_bielec_integral(h1,h2,p1,p2,mo_integrals_map)).gt.1.d-10)then
      print*,i,j,k,l
      print*, get_mo_bielec_integral(h1,h2,p1,p2,mo_integrals_map),transformed_act1_act1_act2_act2(i,j,k,l)
     endif
     ref_energy_act_act += two_body_dm_ab_big_array_act(i,j,k,l) * get_mo_bielec_integral(h1,h2,p1,p2,mo_integrals_map)
    enddo
   enddo
  enddo
 enddo

END_PROVIDER 
