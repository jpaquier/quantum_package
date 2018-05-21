 BEGIN_PROVIDER [double complex, dirac_SCF_density_matrix_ao, (2*dirac_ao_num,2*dirac_ao_num) ]
   implicit none
   BEGIN_DOC
   ! S^{-1}.P.S^{-1}  where P = C.C^t
   END_DOC
   call zgemm('N','C',2*dirac_ao_num,2*dirac_ao_num,elec_num,(1.d0,0.d0), &
        dirac_mo_coef, size(dirac_mo_coef,1), &
        dirac_mo_coef, size(dirac_mo_coef,1), (0.d0,0.d0), &
        dirac_SCF_density_matrix_ao, size(dirac_SCF_density_matrix_ao,1))
 END_PROVIDER

