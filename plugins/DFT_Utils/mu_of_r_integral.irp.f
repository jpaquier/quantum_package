double precision function integral_of_mu_of_r_on_HF(r,mu_in)
 implicit none
 double precision, intent(in) :: r(3), mu_in
 double precision :: integrals_ao(ao_num,ao_num), NAI_pol_mult_erf_ao
 double precision :: mos_array(mo_tot_num)
 integer :: i,j,l,k,m
 double precision :: integrals_mo(elec_alpha_num)
 do k = 1, ao_num
  do m = 1, ao_num
   integrals_ao(m,k) = NAI_pol_mult_erf_ao(m,k,mu_in,r)
  enddo
 enddo
 integrals_mo = 0.d0
 do i = 1, elec_alpha_num
  do k = 1, ao_num
   do m = 1, ao_num
    integrals_mo(i) += mo_coef(m,i) * mo_coef(k,i) * integrals_ao(m,k)
   enddo
  enddo
 enddo
 call give_all_mos_at_r(r,mos_array)
 double precision :: density_beta
 density_beta = 0.d0
 do j = 1, elec_beta_num
  density_beta += mos_array(j)**2 
 enddo
 integral_of_mu_of_r_on_HF = 0.d0
 do i = 1, elec_alpha_num
  integral_of_mu_of_r_on_HF += integrals_mo(i)
 enddo
 integral_of_mu_of_r_on_HF = integral_of_mu_of_r_on_HF * density_beta  


end

 double precision function mu_integral_old(integral_f,r)
 implicit none 
 double precision, intent(in) :: integral_f, r(3)
 double precision :: mos_array(mo_tot_num),local_potential,mu_0,mu_in,tmp0,tmp_in
 double precision :: threshold,mu_tmp,integrals_mo(mo_tot_num,mo_tot_num)
 integer :: i,j
 threshold = 1.d-2
 call give_all_mos_at_r(r,mos_array)
 call local_r12_operator_on_hf(r,r,local_potential)
 mu_0 =  local_potential * dsqrt(dacos(-1.d0)) * 0.5d0
 tmp0 = 0.d0
 call give_all_erf_mu_of_r_kl_mo(integrals_mo,mu_0,r)
 do i = 1, elec_alpha_num
  do j = 1, elec_beta_num
   tmp0+= mos_array(j)**2 * integrals_mo(i,i)
  enddo
 enddo

 mu_in = 50.d0
 tmp_in = 0.d0
 call give_all_erf_mu_of_r_kl_mo(integrals_mo,mu_in,r)
 do i = 1, elec_alpha_num
  do j = 1, elec_beta_num
   tmp_in += mos_array(j)**2 * integrals_mo(i,i)
  enddo
 enddo
 if(integral_f.gt.tmp_in)then
  mu_integral_old = mu_0
! mu_integral = mu_in 
  return
 else
  do while (dabs(mu_in -mu_0).gt.threshold)
   mu_tmp = 0.5d0 * (mu_in + mu_0)
   call give_all_erf_mu_of_r_kl_mo(integrals_mo,mu_tmp,r)
   tmp_in = 0.d0
   do i = 1, elec_alpha_num
    do j = 1, elec_beta_num
     tmp_in += mos_array(j)**2 * integrals_mo(i,i)
    enddo
   enddo
   if((tmp0 - integral_f) * (tmp_in - integral_f).le.0.d0)then
    mu_in = mu_tmp
    else 
    mu_0 = mu_tmp
   endif
  enddo
 endif
 mu_integral_old = (mu_0 + mu_in)*0.5d0
 end

 double precision function mu_integral(integral_f,r)
 implicit none 
 double precision, intent(in) :: integral_f, r(3)
 double precision :: mos_array(mo_tot_num),local_potential,mu_0,mu_in,tmp0,tmp_in
 double precision :: threshold,mu_tmp,integrals_mo(mo_tot_num,mo_tot_num)
 integer :: i,j
 double precision :: integral_of_mu_of_r_on_HF
 threshold = 1.d-2
 call give_all_mos_at_r(r,mos_array)
 call local_r12_operator_on_hf(r,r,local_potential)
 mu_0 =  local_potential * dsqrt(dacos(-1.d0)) * 0.5d0
 tmp0 = 0.d0
 call give_all_erf_mu_of_r_kl_mo(integrals_mo,mu_0,r)
 do i = 1, elec_alpha_num
  do j = 1, elec_beta_num
   tmp0+= mos_array(j)**2 * integrals_mo(i,i)
  enddo
 enddo
 mu_in = 50.d0
 tmp_in = integral_of_mu_of_r_on_HF(r,mu_in)
 if(integral_f.gt.tmp_in)then
! mu_integral = mu_0
  mu_integral = 1000.d0
  return
 else
  do while (dabs(mu_in -mu_0).gt.threshold)
   mu_tmp = 0.5d0 * (mu_in + mu_0)
   tmp_in = integral_of_mu_of_r_on_HF(r,mu_tmp)
   if((tmp0 - integral_f) * (tmp_in - integral_f).le.0.d0)then
    mu_in = mu_tmp
    else 
    mu_0 = mu_tmp
   endif
  enddo
 endif
 mu_integral = (mu_0 + mu_in)*0.5d0
 end

BEGIN_PROVIDER [double precision, HF_mu_of_r_bielec_energy]
 implicit none
 integer :: i,j,k
 double precision :: r(3), weight,tmp,integral_of_mu_of_r_on_HF
 double precision, allocatable :: integrals_mo(:,:),mos_array(:)
 allocate(integrals_mo(mo_tot_num,mo_tot_num),mos_array(mo_tot_num))
 HF_mu_of_r_bielec_energy = 0.d0
 do i = 1, n_points_final_grid
  r(1:3) = final_grid_points(1:3,i)  
  weight = final_weight_functions_at_final_grid_points(i)

  tmp = integral_of_mu_of_r_on_HF(r,mu_of_r_vector(i))
  HF_mu_of_r_bielec_energy += weight * tmp
 enddo

END_PROVIDER 


BEGIN_PROVIDER [double precision, energy_c_LDA_mu_of_r]
 implicit none
 integer :: i,j,k
 double precision :: r(3), weight,tmp,rho_a,rho_b,e_lda
 energy_c_LDA_mu_of_r = 0.d0
 do i = 1, n_points_final_grid
  r(1:3) = final_grid_points(1:3,i)  
  weight = final_weight_functions_at_final_grid_points(i)
  call dm_dft_alpha_beta_at_r(r,rho_a,rho_b)
  call ec_only_lda_sr(mu_of_r_vector(i),rho_a,rho_b,e_lda)
  energy_c_LDA_mu_of_r += weight * e_lda
 enddo
END_PROVIDER 
