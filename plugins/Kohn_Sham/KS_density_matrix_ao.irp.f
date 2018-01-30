BEGIN_PROVIDER [ double precision, KS_density_matrix_ao_alpha, (ao_num,ao_num) ]
   implicit none
   BEGIN_DOC
   ! S^-1 x Alpha density matrix in the ao basis x S^-1
   END_DOC
   
   call dgemm('N','T',ao_num,ao_num,elec_alpha_num,1.d0, &
        mo_coef, size(mo_coef,1), &
        mo_coef, size(mo_coef,1), 0.d0, &
        KS_density_matrix_ao_alpha, size(KS_density_matrix_ao_alpha,1))

END_PROVIDER

BEGIN_PROVIDER [ double precision, KS_density_matrix_ao_beta,  (ao_num,ao_num) ]
   implicit none
   BEGIN_DOC
   ! S^-1 Beta density matrix in the ao basis x S^-1
   END_DOC
   
   call dgemm('N','T',ao_num,ao_num,elec_beta_num,1.d0, &
        mo_coef, size(mo_coef,1), &
        mo_coef, size(mo_coef,1), 0.d0, &
        KS_density_matrix_ao_beta, size(KS_density_matrix_ao_beta,1))

END_PROVIDER
 
BEGIN_PROVIDER [ double precision, KS_density_matrix_ao, (ao_num,ao_num) ]
   implicit none
   BEGIN_DOC
   ! S^-1 Density matrix in the ao basis S^-1
   END_DOC
   ASSERT (size(KS_density_matrix_ao,1) == size(KS_density_matrix_ao_alpha,1))
   if (elec_alpha_num== elec_beta_num) then
     KS_density_matrix_ao = KS_density_matrix_ao_alpha + KS_density_matrix_ao_alpha
   else
     ASSERT (size(KS_density_matrix_ao,1) == size(KS_density_matrix_ao_beta ,1))
     KS_density_matrix_ao = KS_density_matrix_ao_alpha + KS_density_matrix_ao_beta
   endif
   
END_PROVIDER
 
