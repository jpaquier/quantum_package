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

 double precision function mu_integral(integral_f,r)
 implicit none 
 double precision, intent(in) :: integral_f, r(3)
 double precision :: mos_array(mo_tot_num),local_potential,mu_0,mu_in,tmp0,tmp_in
 double precision :: threshold,mu_tmp,integrals_mo(mo_tot_num,mo_tot_num)
 integer :: i,j
 double precision :: integral_of_mu_of_r_on_HF
 threshold = 1.d-4
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


 BEGIN_PROVIDER [double precision, mu_of_r_integral_hf_vector, (n_points_final_grid) ]
 implicit none 
 BEGIN_DOC
 ! mu_of_r and mu_average computation 
 END_DOC
 integer :: i_point
 double precision :: r(3)
 double precision :: local_potential,two_body_dm
 double precision :: cpu0,cpu1,integral_f,mu_integral
 print*,'providing the mu_of_r ...'
 call wall_time(cpu0)
 r = 0.d0
 call integral_of_f_12_on_hf(r,integral_f)
 !$OMP PARALLEL DO &
 !$OMP DEFAULT (NONE)  &
 !$OMP PRIVATE (i_point,r,integral_f) & 
 !$OMP shARED (n_points_final_grid,final_grid_points,mu_of_r_integral_hf_vector) 
 do i_point = 1, n_points_final_grid
  r(1) = final_grid_points(1,i_point)
  r(2) = final_grid_points(2,i_point)
  r(3) = final_grid_points(3,i_point)
  call integral_of_f_12_on_hf(r,integral_f)
  mu_of_r_integral_hf_vector(i_point) = mu_integral(integral_f,r)
 enddo
 !$OMP END PARALLEL DO
 call wall_time(cpu1)
 print*,'Time to provide mu_of_r = ',cpu1-cpu0
 END_PROVIDER 

