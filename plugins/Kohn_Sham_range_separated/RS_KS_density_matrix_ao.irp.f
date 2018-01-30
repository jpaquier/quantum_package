BEGIN_PROVIDER [ double precision, RS_KS_density_matrix_ao_alpha, (ao_num,ao_num) ]
   implicit none
   BEGIN_DOC
   ! S^-1 x Alpha density matrix in the AO basis x S^-1
   END_DOC
   
   call dgemm('N','T',ao_num,ao_num,elec_alpha_num,1.d0, &
        mo_coef, size(mo_coef,1), &
        mo_coef, size(mo_coef,1), 0.d0, &
        RS_KS_density_matrix_ao_alpha, size(RS_KS_density_matrix_ao_alpha,1))

END_PROVIDER

BEGIN_PROVIDER [ double precision, RS_KS_density_matrix_ao_beta,  (ao_num,ao_num) ]
   implicit none
   BEGIN_DOC
   ! S^-1 Beta density matrix in the AO basis x S^-1
   END_DOC
   
   call dgemm('N','T',ao_num,ao_num,elec_beta_num,1.d0, &
        mo_coef, size(mo_coef,1), &
        mo_coef, size(mo_coef,1), 0.d0, &
        RS_KS_density_matrix_ao_beta, size(RS_KS_density_matrix_ao_beta,1))

END_PROVIDER
 
BEGIN_PROVIDER [ double precision, RS_KS_density_matrix_ao, (ao_num,ao_num) ]
   implicit none
   BEGIN_DOC
   ! S^-1 Density matrix in the AO basis S^-1
   END_DOC
   ASSERT (size(RS_KS_density_matrix_ao,1) == size(RS_KS_density_matrix_ao_alpha,1))
   if (elec_alpha_num== elec_beta_num) then
     RS_KS_density_matrix_ao = RS_KS_density_matrix_ao_alpha + RS_KS_density_matrix_ao_alpha
   else
     ASSERT (size(RS_KS_density_matrix_ao,1) == size(RS_KS_density_matrix_ao_beta ,1))
     RS_KS_density_matrix_ao = RS_KS_density_matrix_ao_alpha + RS_KS_density_matrix_ao_beta
   endif
   
END_PROVIDER
 
