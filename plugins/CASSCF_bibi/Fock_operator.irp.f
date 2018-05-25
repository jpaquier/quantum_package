BEGIN_PROVIDER [double precision, MR_Fock_mo_core_virt_from_core_inact, (mo_tot_num,mo_tot_num)]
 implicit none
 integer :: i,j,k,l,iorb,jorb,korb
 double precision :: get_mo_bielec_integral
 MR_Fock_mo_core_virt_from_core_inact = 0.d0
 do i = 1, mo_tot_num
  do j = 1, mo_tot_num 
   MR_Fock_mo_core_virt_from_core_inact(i,j) = mo_mono_elec_integral(i,j)
  enddo
 enddo

 do k = 1, n_core_inact_orb
  korb = list_core_inact(k)
  do i = 1, n_core_inact_orb 
   iorb = list_core_inact(i)
   do j = 1, n_core_inact_orb
    jorb = list_core_inact(j)
    MR_Fock_mo_core_virt_from_core_inact(iorb,jorb) += 2.d0 * get_mo_bielec_integral(iorb,korb,jorb,korb) -  get_mo_bielec_integral(iorb,korb,korb,jorb) 
   enddo
   do j = 1, n_virt_orb
    jorb = list_virt(j)
    MR_Fock_mo_core_virt_from_core_inact(iorb,jorb) += 2.d0 * get_mo_bielec_integral(iorb,korb,jorb,korb) -  get_mo_bielec_integral(iorb,korb,korb,jorb) 
   enddo
  enddo
  do j = 1, n_virt_orb
   jorb = list_virt(j)
   do i = 1, n_virt_orb
    iorb = list_virt(i)
    MR_Fock_mo_core_virt_from_core_inact(iorb,jorb) += 2.d0 * get_mo_bielec_integral(iorb,korb,jorb,korb) -  get_mo_bielec_integral(iorb,korb,korb,jorb) 
   enddo
  enddo
 enddo

END_PROVIDER 

 BEGIN_PROVIDER [double precision, MR_Fock_matrix_alpha_ao, (ao_num, ao_num,N_states)]
&BEGIN_PROVIDER [double precision, MR_Fock_matrix_beta_ao, (ao_num, ao_num,N_states)]
 implicit none
 integer                        :: i,j,k
 do k = 1, N_states
  do j=1,ao_num
    !DIR$ VECTOR ALIGNED
    do i=1,ao_num
      MR_Fock_matrix_alpha_ao(i,j,k) = ao_mono_elec_integral(i,j) + ao_bi_elec_integral_alpha_core_inact(i,j) + ao_bi_elec_integral_alpha_act(i,j,k) 
      MR_Fock_matrix_beta_ao(i,j,k) = ao_mono_elec_integral(i,j) + ao_bi_elec_integral_beta_core_inact(i,j) + ao_bi_elec_integral_beta_act(i,j,k)
    enddo
  enddo
 enddo


 END_PROVIDER 

 BEGIN_PROVIDER [ double precision, ao_bi_elec_integral_alpha_core_inact, (ao_num, ao_num) ]
&BEGIN_PROVIDER [ double precision, ao_bi_elec_integral_beta_core_inact ,  (ao_num, ao_num) ]
&BEGIN_PROVIDER [ double precision, ao_bi_elec_integral_alpha_act,  (ao_num, ao_num,N_states) ]
&BEGIN_PROVIDER [ double precision, ao_bi_elec_integral_beta_act,  (ao_num, ao_num,N_states) ]
&BEGIN_PROVIDER [ double precision, repulsion_elec_core_core ]
 use map_module
 implicit none
 BEGIN_DOC
 ! Alpha MR_Fock matrix in AO basis set
 END_DOC
 
 integer                        :: i,j,k,l,k1,r,s
 integer                        :: i0,j0,k0,l0
 integer*8                      :: p,q
 double precision               :: integral, c0, c1, c2
 double precision               :: ao_bielec_integral, local_threshold
 double precision, allocatable  :: ao_bi_elec_integral_alpha_core_inact_tmp(:,:)
 double precision, allocatable  :: ao_bi_elec_integral_beta_core_inact_tmp(:,:)
 double precision, allocatable  :: ao_bi_elec_integral_alpha_act_tmp(:,:,:)
 double precision, allocatable  :: ao_bi_elec_integral_beta_act_tmp(:,:,:)
 double precision :: repulsion_elec_core_core_tmp
 print*, 'PROVIDING THE FOCK MATRIX ...'
 !DIR$ ATTRIBUTES ALIGN : $IRP_ALIGN :: ao_bi_elec_integral_beta_core_inact_tmp
 !DIR$ ATTRIBUTES ALIGN : $IRP_ALIGN :: ao_bi_elec_integral_alpha_core_inact_tmp
 !DIR$ ATTRIBUTES ALIGN : $IRP_ALIGN :: ao_bi_elec_integral_alpha_act_tmp
 !DIR$ ATTRIBUTES ALIGN : $IRP_ALIGN :: ao_bi_elec_integral_beta_act_tmp

 repulsion_elec_core_core = 0.d0
 ao_bi_elec_integral_alpha_core_inact = 0.d0
 ao_bi_elec_integral_beta_core_inact  = 0.d0
 ao_bi_elec_integral_alpha_act = 0.d0
 ao_bi_elec_integral_beta_act = 0.d0
   PROVIDE ao_bielec_integrals_in_map 
           
   integer(omp_lock_kind) :: lck(ao_num)
   integer*8                      :: i8
   integer                        :: ii(8), jj(8), kk(8), ll(8), k2,m
   integer(cache_map_size_kind)   :: n_elements_max, n_elements
   integer(key_kind), allocatable :: keys(:)
   double precision, allocatable  :: values(:)

   !$OMP PARALLEL DEFAULT(NONE)                                      &
   !$OMP PRIVATE(i,j,l,k1,k,integral,ii,jj,kk,ll,i8,keys,values,n_elements_max,m, &
   !$OMP  n_elements,ao_bi_elec_integral_alpha_core_inact_tmp,ao_bi_elec_integral_beta_core_inact_tmp,ao_bi_elec_integral_alpha_act_tmp,ao_bi_elec_integral_beta_act_tmp,repulsion_elec_core_core_tmp)&
   !$OMP SHARED(ao_num,HF_density_matrix_ao_alpha_core_inact,HF_density_matrix_ao_beta_core_inact,&
   !$OMP  ao_integrals_map, ao_bi_elec_integral_alpha_core_inact, ao_bi_elec_integral_beta_core_inact, ao_bi_elec_integral_alpha_act, ao_bi_elec_integral_beta_act,density_matrix_ao_act_alpha,density_matrix_ao_act_beta,repulsion_elec_core_core,density_matrix_ao_act,N_states) 

   call get_cache_map_n_elements_max(ao_integrals_map,n_elements_max)
   allocate(keys(n_elements_max), values(n_elements_max))
   allocate(ao_bi_elec_integral_alpha_core_inact_tmp(ao_num,ao_num), &
            ao_bi_elec_integral_beta_core_inact_tmp(ao_num,ao_num))
   ao_bi_elec_integral_alpha_core_inact_tmp = 0.d0
   ao_bi_elec_integral_beta_core_inact_tmp  = 0.d0
   allocate(ao_bi_elec_integral_alpha_act_tmp(ao_num,ao_num,N_states))
   allocate(ao_bi_elec_integral_beta_act_tmp(ao_num,ao_num,N_states))
   ao_bi_elec_integral_alpha_act_tmp  = 0.d0
   ao_bi_elec_integral_beta_act_tmp  = 0.d0
   repulsion_elec_core_core_tmp = 0.d0

   !$OMP DO SCHEDULE(dynamic)
   !DIR$ NOVECTOR
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
         integral = (HF_density_matrix_ao_alpha_core_inact(k,l)+HF_density_matrix_ao_beta_core_inact(k,l)) * values(k1)
         repulsion_elec_core_core_tmp += integral * HF_density_matrix_ao_alpha_core_inact(i,j)
         ao_bi_elec_integral_alpha_core_inact_tmp(i,j) += integral
         ao_bi_elec_integral_beta_core_inact_tmp (i,j) += integral
         integral = values(k1)
         repulsion_elec_core_core_tmp -= integral * HF_density_matrix_ao_alpha_core_inact(k,i) * HF_density_matrix_ao_alpha_core_inact(j,l)
         ao_bi_elec_integral_alpha_core_inact_tmp(l,j) -= HF_density_matrix_ao_alpha_core_inact(k,i) * integral
         ao_bi_elec_integral_beta_core_inact_tmp (l,j) -= HF_density_matrix_ao_beta_core_inact (k,i) * integral

         do m = 1, N_states
          integral = density_matrix_ao_act(k,l,m) * values(k1)
          ao_bi_elec_integral_alpha_act_tmp(i,j,m) += integral
          ao_bi_elec_integral_beta_act_tmp(i,j,m) += integral
          integral = values(k1)
          ao_bi_elec_integral_alpha_act_tmp(l,j,m) -= density_matrix_ao_act_alpha(k,i,m) * integral
          ao_bi_elec_integral_beta_act_tmp(l,j,m) -= density_matrix_ao_act_beta(k,i,m) * integral
         enddo
       enddo
     enddo
   enddo
   !$OMP END DO NOWAIT
   !$OMP CRITICAL
   ao_bi_elec_integral_alpha_core_inact += ao_bi_elec_integral_alpha_core_inact_tmp
   !$OMP END CRITICAL
   !$OMP CRITICAL
   ao_bi_elec_integral_beta_core_inact  += ao_bi_elec_integral_beta_core_inact_tmp
   !$OMP END CRITICAL

   !$OMP CRITICAL
   ao_bi_elec_integral_alpha_act  += ao_bi_elec_integral_alpha_act_tmp
   !$OMP END CRITICAL

   !$OMP CRITICAL
   ao_bi_elec_integral_beta_act  += ao_bi_elec_integral_beta_act_tmp
   !$OMP END CRITICAL

   !$OMP CRITICAL
   repulsion_elec_core_core += repulsion_elec_core_core_tmp
   !$OMP END CRITICAL
   deallocate(keys,values,ao_bi_elec_integral_alpha_core_inact_tmp,ao_bi_elec_integral_beta_core_inact_tmp)
   deallocate(ao_bi_elec_integral_beta_act_tmp,ao_bi_elec_integral_alpha_act_tmp)
   !$OMP END PARALLEL


END_PROVIDER

 BEGIN_PROVIDER [double precision, MR_Fock_matrix_alpha_from_act_mo_bis, (mo_tot_num, mo_tot_num,N_states)]
 implicit none
 integer :: i,j,k,l,korb,lorb,m
 double precision :: integral,get_mo_bielec_integral
 MR_Fock_matrix_alpha_from_act_mo_bis = 0.d0
 do m = 1, N_states
  do i = 1, mo_tot_num
   do j = 1, mo_tot_num
    do k = 1, n_act_orb
     korb = list_act(k)
     do l = 1, n_act_orb
      lorb = list_act(l)
      integral = get_mo_bielec_integral(i,korb,j,lorb,mo_integrals_map)
      MR_Fock_matrix_alpha_from_act_mo_bis(j,i,m) += integral * density_matrix_mo_act(korb,lorb,m)
      integral = get_mo_bielec_integral(i,j,korb,lorb,mo_integrals_map)
      MR_Fock_matrix_alpha_from_act_mo_bis(j,i,m) -= integral * density_matrix_mo_act_beta(korb,lorb,m)
     enddo
    enddo
   enddo
  enddo
 enddo

 END_PROVIDER 
