BEGIN_PROVIDER [ double precision, HF_density_matrix_ao_alpha_core_inact, (ao_num_align,ao_num) ]
   implicit none
  BEGIN_DOC
 ! HF_density_matrix_ao_alpha_core_inact(i,j) = \sum_{k = 1, n_core_inact_orb} mo_coef(i,k) * mo_coef(j,k)
  END_DOC
   call dgemm('N','T',ao_num,ao_num,n_core_inact_orb,1.d0, &
        mo_coef, size(mo_coef,1), &
        mo_coef, size(mo_coef,1), 0.d0, &
        HF_density_matrix_ao_alpha_core_inact, size(HF_density_matrix_ao_alpha,1))
END_PROVIDER

BEGIN_PROVIDER [ double precision, HF_density_matrix_ao_beta_core_inact,  (ao_num_align,ao_num) ]
   implicit none
  BEGIN_DOC
 ! HF_density_matrix_ao_beta_core_inact(i,j) = \sum_{k = 1, n_core_inact_orb} mo_coef(i,k) * mo_coef(j,k)
  END_DOC
   call dgemm('N','T',ao_num,ao_num,n_core_inact_orb,1.d0, &
        mo_coef, size(mo_coef,1), &
        mo_coef, size(mo_coef,1), 0.d0, &
        HF_density_matrix_ao_beta_core_inact, size(HF_density_matrix_ao_beta,1))
END_PROVIDER

 BEGIN_PROVIDER [ double precision, density_matrix_ao_act_alpha, (ao_num_align,ao_num) ]
&BEGIN_PROVIDER [ double precision, density_matrix_ao_act_beta, (ao_num_align,ao_num) ]
&BEGIN_PROVIDER [ double precision, density_matrix_ao_act, (ao_num_align,ao_num) ]
   implicit none
  BEGIN_DOC
 ! density_matrix_ao_act_alpha(i,j) = \sum_{k_act,l_act = 1, n_act_orb} mo_coef(i,l_act) * mo_coef(j,k_act) * one_body_dm_mo_alpha_average(l_act,k_act)
 ! density_matrix_ao_act_beta(i,j) = \sum_{k_act,l_act = 1, n_act_orb} mo_coef(i,l_act) * mo_coef(j,k_act) * one_body_dm_mo_beta_average(l_act,k_act)
  END_DOC
   integer :: i,j,k,l,korb,lorb
   density_matrix_ao_act_alpha = 0.d0
   density_matrix_ao_act_beta = 0.d0
   do k = 1, n_act_orb
    korb = list_act(k)
    do l = 1, n_act_orb
     lorb = list_act(l)
     do i = 1, ao_num
      do j = 1, ao_num
        density_matrix_ao_act_alpha(i,j) += mo_coef(i,korb) * mo_coef(j,lorb) * one_body_dm_mo_alpha_average(korb,lorb)
        density_matrix_ao_act_beta(i,j) += mo_coef(i,korb) * mo_coef(j,lorb) * one_body_dm_mo_beta_average(korb,lorb)
      enddo
     enddo
    enddo
   enddo
  density_matrix_ao_act = density_matrix_ao_act_beta + density_matrix_ao_act_alpha

END_PROVIDER

 BEGIN_PROVIDER [double precision, density_matrix_mo_act, (mo_tot_num, mo_tot_num)]
&BEGIN_PROVIDER [double precision, density_matrix_mo_act_alpha, (mo_tot_num, mo_tot_num)]
&BEGIN_PROVIDER [double precision, density_matrix_mo_act_beta, (mo_tot_num, mo_tot_num)]
 implicit none
 integer :: i,j,k,l,korb,lorb
 density_matrix_mo_act = 0.d0
 density_matrix_mo_act_alpha = 0.d0
 density_matrix_mo_act_beta = 0.d0
   do k = 1, n_act_orb
    korb = list_act(k)
    do l = 1, n_act_orb
     lorb = list_act(l)
     density_matrix_mo_act(korb,lorb) = (one_body_dm_mo_alpha_average(korb,lorb) + one_body_dm_mo_beta_average(korb,lorb))
     density_matrix_mo_act_alpha(korb,lorb) = one_body_dm_mo_alpha_average(korb,lorb) 
     density_matrix_mo_act_beta(korb,lorb) = one_body_dm_mo_beta_average(korb,lorb) 
    enddo
   enddo

END_PROVIDER 

BEGIN_PROVIDER [ double precision, HF_density_matrix_ao, (ao_num_align,ao_num) ]
   implicit none
   BEGIN_DOC
   ! S^-1 Density matrix in the AO basis S^-1
   END_DOC
   ASSERT (size(HF_density_matrix_ao,1) == size(HF_density_matrix_ao_alpha,1))
   if (elec_alpha_num== elec_beta_num) then
     HF_density_matrix_ao = HF_density_matrix_ao_alpha + HF_density_matrix_ao_alpha
   else
     ASSERT (size(HF_density_matrix_ao,1) == size(HF_density_matrix_ao_beta ,1))
     HF_density_matrix_ao = HF_density_matrix_ao_alpha + HF_density_matrix_ao_beta
   endif
   integer :: i
   do i = 1, ao_num
    write(33,'(1000(F16.10,X))') HF_density_matrix_ao(i,:)
   enddo
   
END_PROVIDER
 
