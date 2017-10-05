BEGIN_PROVIDER [double precision, mono_elec_core_energy]
 implicit none
 integer :: p,q
 mono_elec_core_energy = 0.d0
 do p = 1, ao_num
  do q = 1, ao_num
   mono_elec_core_energy += 2.d0 * HF_density_matrix_ao_alpha_core_inact(q,p) * ao_mono_elec_integral(q,p)
  enddo
 enddo
 

END_PROVIDER 


BEGIN_PROVIDER [double precision, core_energy_from_core]
 implicit none
 BEGIN_DOC
! energy of the block of doubly occupied electrons NOT TAKING INTO ACCOUNT THE ACTIVE ELECTRONS 
 END_DOC
 core_energy_from_core = 0.d0 ! nuclear_repulsion
 integer :: i,j

 do j=1,ao_num
   do i=1,ao_num
     core_energy_from_core += 0.5d0 * (                                          &
         (ao_mono_elec_integral(i,j) + Fock_matrix_alpha_ao(i,j) ) *  HF_density_matrix_ao_alpha_core_inact(i,j) +&
         (ao_mono_elec_integral(i,j) + Fock_matrix_beta_ao (i,j) ) *  HF_density_matrix_ao_beta_core_inact(i,j) )
   enddo
 enddo

END_PROVIDER 


 BEGIN_PROVIDER [double precision, coulomb_aa_core_energy_from_act]
&BEGIN_PROVIDER [double precision, coulomb_bb_core_energy_from_act]
&BEGIN_PROVIDER [double precision, coulomb_ab_core_energy_from_act]
&BEGIN_PROVIDER [double precision, coulomb_core_energy_from_act_total]
&BEGIN_PROVIDER [double precision, alpha_exch_core_energy_from_act]
&BEGIN_PROVIDER [double precision, beta_exch_core_energy_from_act]
 implicit none
 BEGIN_DOC
! coulomb_core_energy_from_act = \sum_{m,n,p,q} (mn|pq) * HF_density_matrix_ao_alpha_core_inact(m,n) *2 * density_matrix_ao_act(p,q)
! alpha_exch_core_energy_from_act = \sum_{m,n,p,q} (mn|pq) * HF_density_matrix_ao_alpha_core_inact(m,p)  * density_matrix_ao_act_alpha(n,q)
! beta_exch_core_energy_from_act = \sum_{m,n,p,q} (mn|pq) * HF_density_matrix_ao_beta_core_inact(m,p)  * density_matrix_ao_act_beta(n,q)
 END_DOC

 integer                        :: i,j,k,l,k1
 double precision               :: integral
 coulomb_aa_core_energy_from_act = 0.d0
 coulomb_bb_core_energy_from_act = 0.d0
 coulomb_ab_core_energy_from_act = 0.d0
 alpha_exch_core_energy_from_act = 0.d0
 beta_exch_core_energy_from_act = 0.d0

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
  !!OMP SHARED  density_matrix_ao_act,density_matrix_ao_act_alpha,density_matrix_ao_act_beta,n_elements_max) &
  !!$OMP REDUCTION (+:coulomb_aa_core_energy_from_act)        &
  !!$OMP REDUCTION (+:coulomb_bb_core_energy_from_act)        &
  !!$OMP REDUCTION (+:coulomb_ab_core_energy_from_act)        &
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
         ! a-a
         coulomb_aa_core_energy_from_act += HF_density_matrix_ao_alpha_core_inact(k,l) *  density_matrix_ao_act_alpha(i,j) * values(k1)
         coulomb_bb_core_energy_from_act += HF_density_matrix_ao_beta_core_inact(k,l) *  density_matrix_ao_act_beta(i,j) * values(k1)
         coulomb_ab_core_energy_from_act += HF_density_matrix_ao_beta_core_inact(k,l) *  density_matrix_ao_act_alpha(i,j) * values(k1)
         coulomb_ab_core_energy_from_act += HF_density_matrix_ao_alpha_core_inact(k,l) *  density_matrix_ao_act_beta(i,j) * values(k1)
         alpha_exch_core_energy_from_act -= HF_density_matrix_ao_alpha_core_inact(i,k) * density_matrix_ao_act_alpha(j,l) * values(k1)
         beta_exch_core_energy_from_act -= HF_density_matrix_ao_beta_core_inact(i,k) * density_matrix_ao_act_beta(j,l) * values(k1)
       enddo
     enddo
   enddo
  !!$OMP END DO NOWAIT
   deallocate(keys,values)
  !!$OMP END PARALLEL
  coulomb_core_energy_from_act_total = coulomb_aa_core_energy_from_act + coulomb_bb_core_energy_from_act + coulomb_ab_core_energy_from_act
END_PROVIDER 

 BEGIN_PROVIDER [double precision, coulomb_energy_core_act_mo_total]
&BEGIN_PROVIDER [double precision, coulomb_energy_core_act_mo_aa]
&BEGIN_PROVIDER [double precision, coulomb_energy_core_act_mo_bb]
&BEGIN_PROVIDER [double precision, coulomb_energy_core_act_mo_ab]
&BEGIN_PROVIDER [double precision, exch_energy_core_act_mo_alpha]
&BEGIN_PROVIDER [double precision, exch_energy_core_act_mo_beta]
 implicit none
 integer :: i,j,iorb,jorb,k,korb
 double precision :: integral,get_mo_bielec_integral
 coulomb_energy_core_act_mo_aa = 0.d0
 coulomb_energy_core_act_mo_bb = 0.d0
 coulomb_energy_core_act_mo_ab = 0.d0
 exch_energy_core_act_mo_alpha = 0.d0
 exch_energy_core_act_mo_beta = 0.d0
 double precision :: accu
 accu = 0.d0
 do i = 1, n_core_inact_orb
  iorb = list_core_inact(i) 
! accu += 2.d0 * mo_mono_elec_integral(iorb,iorb)
  do j = 1, n_core_inact_orb
   jorb = list_core_inact(j) 
   integral = mo_bielec_integral_jj(iorb,jorb)
   accu += 2.d0 * integral 
   integral = mo_bielec_integral_jj_exchange(iorb,jorb)
   accu -= integral 
  enddo
  do j = 1, n_act_orb
   jorb = list_act(j)
   do k = 1, n_act_orb
    korb = list_act(k)
    integral = get_mo_bielec_integral(iorb,jorb,iorb,korb,mo_integrals_map)
    if(density_matrix_mo_act(jorb,korb).lt.1.d-10)cycle
    coulomb_energy_core_act_mo_aa += 1.d0 * density_matrix_mo_act_alpha(jorb,korb) * integral 
    coulomb_energy_core_act_mo_bb += 1.d0 * density_matrix_mo_act_beta(jorb,korb) * integral 
    coulomb_energy_core_act_mo_ab += 1.d0 * density_matrix_mo_act_beta(jorb,korb) * integral + 1.d0 * density_matrix_mo_act_alpha(jorb,korb) * integral 
    integral = get_mo_bielec_integral(iorb,jorb,korb,iorb,mo_integrals_map)
    exch_energy_core_act_mo_alpha -= density_matrix_mo_act_alpha(jorb,korb) * integral
    exch_energy_core_act_mo_beta  -= density_matrix_mo_act_alpha(jorb,korb) * integral
   enddo
  enddo
 enddo
 coulomb_energy_core_act_mo_total = coulomb_energy_core_act_mo_aa + coulomb_energy_core_act_mo_bb + coulomb_energy_core_act_mo_ab
 
 
END_PROVIDER 


BEGIN_PROVIDER [double precision, ref_energy_core_core_and_core_act]
 implicit none
 ref_energy_core_core_and_core_act   = core_energy_from_core 
 ref_energy_core_core_and_core_act  += 0.5d0*(coulomb_aa_core_energy_from_act + coulomb_bb_core_energy_from_act + coulomb_ab_core_energy_from_act)
 ref_energy_core_core_and_core_act  += 0.5d0*(alpha_exch_core_energy_from_act + beta_exch_core_energy_from_act)
END_PROVIDER 

BEGIN_PROVIDER [double precision, ref_energy_act_act]
 BEGIN_DOC 
 ! purely active part of the two body dm 
 END_DOC
 implicit none
 integer :: i,j,k,l
 double precision  :: get_mo_bielec_integral
 ref_energy_act_act = 0.d0

 ! Diag part of the active two body dm
 do i = 1, n_act_orb
  do j = 1, n_act_orb
   double precision :: integral
   integral = transformed_act1_act1_act2_act2(i,i,j,j)
   ref_energy_act_act += (two_body_dm_ab_diag_act(j,i) + two_body_dm_aa_diag_act(j,i) + two_body_dm_bb_diag_act(j,i) ) *  integral
   integral = transformed_act1_act1_act2_act2(i,j,i,j)
   ref_energy_act_act += (two_body_dm_bb_diag_exchange_act(j,i) + two_body_dm_aa_diag_exchange_act(j,i) ) * integral
  enddo
 enddo

 ! Extra diagonal part of the active part of the two body dm 
 do l = 1, n_act_orb  ! p2 
  do k = 1, n_act_orb  ! h2 
   do j = 1, n_act_orb  ! p1 
    do i = 1,n_act_orb   ! h1 
     integral = transformed_act1_act1_act2_act2(i,j,k,l)
     ref_energy_act_act += (two_body_dm_ab_big_array_act(i,j,k,l) + two_body_dm_aa_big_array_act(i,j,k,l) + two_body_dm_bb_big_array_act(i,j,k,l)) * integral
    enddo
   enddo
  enddo
 enddo
 ref_energy_act_act = ref_energy_act_act * 0.5d0

END_PROVIDER 

BEGIN_PROVIDER [double precision, ref_energy_act_mono]
 implicit none
 integer :: p,q
 ref_energy_act_mono = 0.d0
 do p = 1, ao_num
  do q = 1, ao_num
   ref_energy_act_mono += density_matrix_ao_act(q,p) * ao_mono_elec_integral(q,p)
  enddo
 enddo
END_PROVIDER 

BEGIN_PROVIDER [double precision,reference_energy_superci ]
implicit none
 if(n_act_orb.gt.0)then
  reference_energy_superci = coulomb_core_energy_from_act_total +alpha_exch_core_energy_from_act + beta_exch_core_energy_from_act + & 
                             mono_elec_core_energy + repulsion_elec_core_core+ ref_energy_act_mono + ref_energy_act_act + nuclear_repulsion
 else
  reference_energy_superci = mono_elec_core_energy + repulsion_elec_core_core+  nuclear_repulsion
 endif

END_PROVIDER 

