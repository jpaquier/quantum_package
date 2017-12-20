BEGIN_PROVIDER [ double precision, MR_Fock_matrix_alpha_mo, (mo_tot_num,mo_tot_num,N_states) ]
   implicit none
   BEGIN_DOC
   ! MR_Fock matrix for a ALPHA excitation in the MO basis
   END_DOC
   integer :: m
   double precision, allocatable  :: T(:,:)
   allocate ( T(ao_num,mo_tot_num) )
   !DIR$ ATTRIBUTES ALIGN : $IRP_ALIGN :: T
   do m = 1, N_states
    call dgemm('N','N', ao_num, mo_tot_num, ao_num,                   &
        1.d0, MR_Fock_matrix_alpha_ao(1,1,m),size(MR_Fock_matrix_alpha_ao,1),      &
        mo_coef, size(mo_coef,1),                                     &
        0.d0, T, ao_num)
    call dgemm('T','N', mo_tot_num, mo_tot_num, ao_num,               &
        1.d0, mo_coef,size(mo_coef,1),                                &
        T, size(T,1),                                                 &
        0.d0, MR_Fock_matrix_alpha_mo(1,1,m), mo_tot_num)
   enddo
   deallocate(T)
END_PROVIDER


BEGIN_PROVIDER [ double precision, MR_Fock_matrix_beta_mo, (mo_tot_num,mo_tot_num,N_States) ]
   implicit none
   BEGIN_DOC
   ! MR_Fock matrix for a BETA excitation in the MO basis
   END_DOC
   double precision, allocatable  :: T(:,:)
   integer :: m
   allocate ( T(ao_num,mo_tot_num) )
   !DIR$ ATTRIBUTES ALIGN : $IRP_ALIGN :: T
   do m = 1, N_states
    call dgemm('N','N', ao_num, mo_tot_num, ao_num,                   &
        1.d0, MR_Fock_matrix_beta_ao(1,1,m),size(MR_Fock_matrix_beta_ao,1),      &
        mo_coef, size(mo_coef,1),                                     &
        0.d0, T, ao_num)
    call dgemm('T','N', mo_tot_num, mo_tot_num, ao_num,               &
        1.d0, mo_coef,size(mo_coef,1),                                &
        T, size(T,1),                                                 &
        0.d0, MR_Fock_matrix_beta_mo(1,1,m), mo_tot_num)
   enddo
   deallocate(T)
END_PROVIDER

BEGIN_PROVIDER [double precision, MR_Fock_canonical_MO, (mo_tot_num, mo_tot_num, N_states)]
   implicit none
   BEGIN_DOC
   ! Fock matrix on the MO basis.
   ! For open shells, the ROHF Fock Matrix is
   !
   !  |   F-K    |  F + K/2  |    F     |
   !  |---------------------------------|
   !  | F + K/2  |     F     |  F - K/2 |
   !  |---------------------------------|
   !  |    F     |  F - K/2  |  F + K   |
   !
   ! F = 1/2 (Fa + Fb)
   !
   ! K = Fb - Fa
   !
   END_DOC
  integer :: i,j,k,l,iorb,jorb,m
  MR_Fock_canonical_MO = 0.d0
 do m = 1, N_states
  do i = 1, n_core_inact_orb
   iorb = list_core_inact(i)
   do j = 1, n_core_inact_orb
    jorb = list_core_inact(j)
    MR_Fock_canonical_MO(jorb,iorb,m) = 0.5d0 * (MR_Fock_matrix_alpha_mo(iorb,jorb,m) + MR_Fock_matrix_beta_mo(iorb,jorb,m)) & 
                                      + (MR_Fock_matrix_alpha_mo(iorb,jorb,m) - MR_Fock_matrix_beta_mo(iorb,jorb,m))
    MR_Fock_canonical_MO(iorb,jorb,m) = 0.5d0 * (MR_Fock_matrix_alpha_mo(iorb,jorb,m) + MR_Fock_matrix_beta_mo(iorb,jorb,m)) & 
                                      + (MR_Fock_matrix_alpha_mo(iorb,jorb,m) - MR_Fock_matrix_beta_mo(iorb,jorb,m))
   enddo

   do j = 1, n_act_orb
    jorb = list_act(j)
    MR_Fock_canonical_MO(jorb,iorb,m) = 0.5d0 * (MR_Fock_matrix_alpha_mo(iorb,jorb,m) + MR_Fock_matrix_beta_mo(iorb,jorb,m)) & 
                                      - (MR_Fock_matrix_alpha_mo(iorb,jorb,m) - MR_Fock_matrix_beta_mo(iorb,jorb,m))
    MR_Fock_canonical_MO(iorb,jorb,m) = 0.5d0 * (MR_Fock_matrix_alpha_mo(iorb,jorb,m) + MR_Fock_matrix_beta_mo(iorb,jorb,m)) & 
                                      - (MR_Fock_matrix_alpha_mo(iorb,jorb,m) - MR_Fock_matrix_beta_mo(iorb,jorb,m))
   enddo
  

   do j = 1, n_virt_orb
    jorb = list_virt(j)
    MR_Fock_canonical_MO(jorb,iorb,m) = 0.5d0 * (MR_Fock_matrix_alpha_mo(iorb,jorb,m) + MR_Fock_matrix_beta_mo(iorb,jorb,m)) 
    MR_Fock_canonical_MO(iorb,jorb,m) = 0.5d0 * (MR_Fock_matrix_alpha_mo(iorb,jorb,m) + MR_Fock_matrix_beta_mo(iorb,jorb,m))
   enddo
  enddo


  do i = 1, n_act_orb
   iorb = list_act(i)
   do j = 1, n_act_orb
    jorb = list_act(j)
    MR_Fock_canonical_MO(jorb,iorb,m) = 0.5d0 * (MR_Fock_matrix_alpha_mo(iorb,jorb,m) + MR_Fock_matrix_beta_mo(iorb,jorb,m)) 
    MR_Fock_canonical_MO(iorb,jorb,m) = 0.5d0 * (MR_Fock_matrix_alpha_mo(iorb,jorb,m) + MR_Fock_matrix_beta_mo(iorb,jorb,m)) 
   enddo
   do j = 1, n_virt_orb
    jorb = list_virt(j)
    MR_Fock_canonical_MO(jorb,iorb,m) = 0.5d0 * (MR_Fock_matrix_alpha_mo(iorb,jorb,m) + MR_Fock_matrix_beta_mo(iorb,jorb,m)) & 
                                      + (MR_Fock_matrix_alpha_mo(iorb,jorb,m) - MR_Fock_matrix_beta_mo(iorb,jorb,m))
    MR_Fock_canonical_MO(iorb,jorb,m) = 0.5d0 * (MR_Fock_matrix_alpha_mo(iorb,jorb,m) + MR_Fock_matrix_beta_mo(iorb,jorb,m)) & 
                                      + (MR_Fock_matrix_alpha_mo(iorb,jorb,m) - MR_Fock_matrix_beta_mo(iorb,jorb,m))
   enddo
  enddo

  do i = 1, n_virt_orb
   iorb = list_virt(i)
   do j = 1, n_virt_orb
    jorb = list_virt(j)
    MR_Fock_canonical_MO(jorb,iorb,m) = 0.5d0 * (MR_Fock_matrix_alpha_mo(iorb,jorb,m) + MR_Fock_matrix_beta_mo(iorb,jorb,m)) & 
                                      - (MR_Fock_matrix_alpha_mo(iorb,jorb,m) - MR_Fock_matrix_beta_mo(iorb,jorb,m))
    MR_Fock_canonical_MO(iorb,jorb,m) = 0.5d0 * (MR_Fock_matrix_alpha_mo(iorb,jorb,m) + MR_Fock_matrix_beta_mo(iorb,jorb,m)) & 
                                      - (MR_Fock_matrix_alpha_mo(iorb,jorb,m) - MR_Fock_matrix_beta_mo(iorb,jorb,m))
   enddo
  enddo
 enddo



END_PROVIDER 

BEGIN_PROVIDER [ double precision, MR_Fock_matrix_alpha_beta_spin_average_mo, (mo_tot_num,mo_tot_num,N_states) ]
   implicit none
   BEGIN_DOC
   ! Average alpha/beta MR_Fock matrix in the MO basis
   END_DOC

  MR_Fock_matrix_alpha_beta_spin_average_mo = MR_Fock_matrix_beta_mo + MR_Fock_matrix_alpha_mo
  MR_Fock_matrix_alpha_beta_spin_average_mo = MR_Fock_matrix_alpha_beta_spin_average_mo * 0.5d0

END_PROVIDER 

BEGIN_PROVIDER [double precision, MR_Fock_matrix_spin_and_state_average_mo, (mo_tot_num,mo_tot_num)]
  implicit none
  integer :: i
  BEGIN_DOC
 ! State average and spin average MR_Fock matrix in the MO basis
  END_DOC
  MR_Fock_matrix_spin_and_state_average_mo = 0.d0
  do i = 1, N_states
   MR_Fock_matrix_spin_and_state_average_mo(:,:) += MR_Fock_matrix_alpha_beta_spin_average_mo(:,:,i) * state_average_weight(i)
  enddo
END_PROVIDER 


BEGIN_PROVIDER [ double precision, MR_Fock_matrix_alpha_from_act_mo, (mo_tot_num,mo_tot_num,N_states) ]
   implicit none
   BEGIN_DOC
   ! MR_Fock matrix for a ALPHA excitation in the MO basis
   END_DOC
   integer :: m
   double precision, allocatable  :: T(:,:)
   allocate ( T(ao_num,mo_tot_num) )
 do m = 1, N_states
   !DIR$ ATTRIBUTES ALIGN : $IRP_ALIGN :: T
   call dgemm('N','N', ao_num, mo_tot_num, ao_num,                   &
       1.d0, ao_bi_elec_integral_alpha_act(1,1,m),size(ao_bi_elec_integral_alpha_act,1),      &
       mo_coef, size(mo_coef,1),                                     &
       0.d0, T, ao_num)
   call dgemm('T','N', mo_tot_num, mo_tot_num, ao_num,               &
       1.d0, mo_coef,size(mo_coef,1),                                &
       T, size(T,1),                                                 &
       0.d0, MR_Fock_matrix_alpha_from_act_mo(1,1,m), mo_tot_num)
 enddo
   deallocate(T)
END_PROVIDER

BEGIN_PROVIDER [ double precision, MR_Fock_matrix_beta_from_act_mo, (mo_tot_num,mo_tot_num,N_states) ]
   implicit none
   BEGIN_DOC
   ! MR_Fock matrix for a ALPHA excitation in the MO basis
   END_DOC
   integer :: m
   double precision, allocatable  :: T(:,:)
   allocate ( T(ao_num,mo_tot_num) )
 do m = 1, N_states
   !DIR$ ATTRIBUTES ALIGN : $IRP_ALIGN :: T
   call dgemm('N','N', ao_num, mo_tot_num, ao_num,                   &
       1.d0, ao_bi_elec_integral_beta_act(1,1,m),size(ao_bi_elec_integral_beta_act,1),      &
       mo_coef, size(mo_coef,1),                                     &
       0.d0, T, ao_num)
   call dgemm('T','N', mo_tot_num, mo_tot_num, ao_num,               &
       1.d0, mo_coef,size(mo_coef,1),                                &
       T, size(T,1),                                                 &
       0.d0, MR_Fock_matrix_beta_from_act_mo(1,1,m), mo_tot_num)
 enddo
   deallocate(T)
END_PROVIDER
