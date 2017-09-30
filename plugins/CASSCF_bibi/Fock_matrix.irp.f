BEGIN_PROVIDER [ double precision, Fock_matrix_alpha_mo, (mo_tot_num_align,mo_tot_num) ]
   implicit none
   BEGIN_DOC
   ! Fock matrix for a ALPHA excitation in the MO basis
   END_DOC
   double precision, allocatable  :: T(:,:)
   allocate ( T(ao_num_align,mo_tot_num) )
   !DIR$ ATTRIBUTES ALIGN : $IRP_ALIGN :: T
   call dgemm('N','N', ao_num, mo_tot_num, ao_num,                   &
       1.d0, Fock_matrix_ao,size(Fock_matrix_alpha_ao,1),      &
       mo_coef, size(mo_coef,1),                                     &
       0.d0, T, ao_num_align)
   call dgemm('T','N', mo_tot_num, mo_tot_num, ao_num,               &
       1.d0, mo_coef,size(mo_coef,1),                                &
       T, size(T,1),                                                 &
       0.d0, Fock_matrix_alpha_mo, mo_tot_num_align)
   deallocate(T)
END_PROVIDER


BEGIN_PROVIDER [ double precision, Fock_matrix_beta_mo, (mo_tot_num_align,mo_tot_num) ]
   implicit none
   BEGIN_DOC
   ! Fock matrix for a BETA excitation in the MO basis
   END_DOC
   double precision, allocatable  :: T(:,:)
   allocate ( T(ao_num_align,mo_tot_num) )
   !DIR$ ATTRIBUTES ALIGN : $IRP_ALIGN :: T
   call dgemm('N','N', ao_num, mo_tot_num, ao_num,                   &
       1.d0, Fock_matrix_ao,size(Fock_matrix_beta_ao,1),      &
       mo_coef, size(mo_coef,1),                                     &
       0.d0, T, ao_num_align)
   call dgemm('T','N', mo_tot_num, mo_tot_num, ao_num,               &
       1.d0, mo_coef,size(mo_coef,1),                                &
       T, size(T,1),                                                 &
       0.d0, Fock_matrix_beta_mo, mo_tot_num_align)
   deallocate(T)
END_PROVIDER


BEGIN_PROVIDER [ double precision, Fock_matrix_alpha_beta_average_mo, (mo_tot_num_align,mo_tot_num) ]
   implicit none
   BEGIN_DOC
   ! Average alpha/beta Fock matrix in the MO basis
   END_DOC

  Fock_matrix_alpha_beta_average_mo = Fock_matrix_beta_mo + Fock_matrix_alpha_mo
  Fock_matrix_alpha_beta_average_mo = Fock_matrix_alpha_beta_average_mo * 0.5d0

END_PROVIDER 
