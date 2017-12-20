BEGIN_PROVIDER [double precision, mo_general_density_alpha, (mo_tot_num,mo_tot_num)]
 implicit none
 integer  :: i,j,k,l
 mo_general_density_alpha = 0.d0
 do k = 1, N_states
  do j = 1, mo_tot_num
   do i = 1, mo_tot_num
    mo_general_density_alpha(i,j) +=  one_body_dm_mo_alpha_generators_restart(i,j,k) * state_average_weight(k) / norm_generators_restart(k)
   enddo
  enddo
 enddo

END_PROVIDER


BEGIN_PROVIDER [double precision, mo_general_density_beta, (mo_tot_num,mo_tot_num)]
 implicit none
 integer  :: i,j,k,l
 mo_general_density_beta = 0.d0
 do k = 1, N_states
  do j = 1, mo_tot_num
   do i = 1, mo_tot_num
    mo_general_density_beta(i,j) +=  one_body_dm_mo_beta_generators_restart(i,j,k) * state_average_weight(k) / norm_generators_restart(k)
   enddo
  enddo
 enddo

END_PROVIDER


