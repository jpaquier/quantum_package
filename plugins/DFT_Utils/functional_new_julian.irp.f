 BEGIN_PROVIDER [double precision, energy_Hxc, (N_states)]
 implicit none
 include 'Utils/constants.include.F'
 integer :: istate
 do istate = 1, N_states
  energy_Hxc(istate) =  pi/(2.d0 * mu_erf**2) * E_cor_tot(istate) 
 enddo
 END_PROVIDER

subroutine give_epsilon_h(r,eps_h)
  implicit none
  double precision, intent(in)  :: r(3)
  double precision, intent(out) :: eps_h(N_states)
  double precision :: mu,rho,accu
  double precision :: NAI_pol_mult_erf_ao
  integer :: istate,i,j
  double precision, allocatable :: aos_array(:)
  allocate(aos_array(ao_num)) 
  mu=1.d6
  eps_h = 0.d0
  call give_all_aos_at_r(r,aos_array)
  do istate = 1, N_states
   accu = 0d0
   rho = 0.d0
   do i = 1, ao_num
    do j = 1, ao_num
     rho += (one_body_dm_alpha_ao_for_dft(j,i,istate)+one_body_dm_beta_ao_for_dft(j,i,istate))*aos_array(i)*aos_array(j)
     accu += (one_body_dm_alpha_ao_for_dft(j,i,istate)+one_body_dm_beta_ao_for_dft(j,i,istate))*NAI_pol_mult_erf_ao(j,i,mu,r)
    enddo
   enddo
   eps_h(istate) = 0.5*accu*rho
  enddo
  deallocate(aos_array)
end

subroutine give_epsilon_h_SR(r,mu,eps_h_SR)
  implicit none
  double precision, intent(in)  :: r(3)
  double precision, intent(in)  :: mu
  double precision, intent(out) :: eps_h_SR(N_states)
  double precision :: rho,accu
  double precision :: NAI_pol_mult_erf_ao
  integer :: istate,i,j
  double precision, allocatable :: aos_array(:),eps_h(:)
  allocate(aos_array(ao_num),eps_h(N_states))
  eps_h_SR = 0.d0
  call give_all_aos_at_r(r,aos_array)
  do istate = 1, N_states
   accu = 0d0
   rho = 0.d0
   do i = 1, ao_num
    do j = 1, ao_num
     rho += (one_body_dm_alpha_ao_for_dft(j,i,istate)+one_body_dm_beta_ao_for_dft(j,i,istate))*aos_array(i)*aos_array(j)
     accu += (one_body_dm_alpha_ao_for_dft(j,i,istate)+one_body_dm_beta_ao_for_dft(j,i,istate))*NAI_pol_mult_erf_ao(j,i,mu,r)
    enddo
   enddo
   call give_epsilon_h(r,eps_h)
   eps_h_SR(istate) = eps_h(istate) - 0.5*accu*rho
  enddo
  deallocate(aos_array,eps_h)
end

subroutine give_epsilon_hxc_PBE_mu_corrected(mu,r,eps_hxc_on_top_PBE)
  implicit none
  double precision, intent(in)  :: mu , r(3)
  double precision, intent(out) :: eps_hxc_on_top_PBE(N_states)
  double precision :: two_dm_in_r, pi,e_PBE_hxc(N_states),ec_PBE(N_states),ex_PBE(N_states),eps_h(N_states),beta(N_states),on_top_two_dm_in_r_mu_corrected
  double precision :: aos_array(ao_num), grad_aos_array(3,ao_num)
  double precision :: rho_a(N_states),rho_b(N_states)
  double precision :: grad_rho_a(3,N_states),grad_rho_b(3,N_states)
  double precision :: grad_rho_a_2(N_states),grad_rho_b_2(N_states),grad_rho_a_b(N_states)
  double precision :: rhoc,rhoo,sigmacc,sigmaco,sigmaoo,vrhoc,vrhoo,vsigmacc,vsigmaco,vsigmaoo
  integer :: m, istate
  pi = 4d0 * datan(1d0)

  eps_hxc_on_top_PBE = 0d0
  call density_and_grad_alpha_beta_and_all_aos_and_grad_aos_at_r(r,rho_a,rho_b,grad_rho_a, grad_rho_b, aos_array, grad_aos_array)
  grad_rho_a_2 = 0.d0
  grad_rho_b_2 = 0.d0
  grad_rho_a_b = 0.d0
  do istate = 1, N_states
   do m = 1, 3
    grad_rho_a_2(istate) += grad_rho_a(m,istate)*grad_rho_a(m,istate)
    grad_rho_b_2(istate) += grad_rho_b(m,istate)*grad_rho_b(m,istate)
    grad_rho_a_b(istate) += grad_rho_a(m,istate)*grad_rho_b(m,istate)
   enddo
  enddo
  do istate = 1, N_states
   ! convertion from (alpha,beta) formalism to (closed, open) formalism
   call rho_ab_to_rho_oc(rho_a(istate),rho_b(istate),rhoo,rhoc)
   call grad_rho_ab_to_grad_rho_oc(grad_rho_a_2(istate),grad_rho_b_2(istate),grad_rho_a_b(istate),sigmaoo,sigmacc,sigmaco)
   call Ec_sr_PBE(0d0,rhoc,rhoo,sigmacc,sigmaco,sigmaoo,ec_PBE(istate))
   call Energy_x_pbe_sr(0d0,rho_a,rho_b,grad_rho_a_2(istate),grad_rho_b_2(istate),grad_rho_a_b(istate),ex_PBE(istate))
   call give_epsilon_h(r,eps_h)
   e_PBE_hxc(istate)=eps_h(istate)+ec_PBE(istate)+ex_PBE(istate) 
   beta(istate) = (2d0*e_PBE_hxc(istate))/(pi*on_top_two_dm_in_r_mu_corrected(mu,r,istate))
   eps_hxc_on_top_PBE(istate)=e_PBE_hxc(istate)/(1d0+beta(istate)*mu**2d0)
  enddo
 end


subroutine give_epsilon_hxc_PBE_SR_mu_corrected(mu,r,eps_hxc_on_top_PBE_SR)
  implicit none
  double precision, intent(in)  :: mu , r(3)
  double precision, intent(out) :: eps_hxc_on_top_PBE_SR(N_states)
  double precision :: two_dm_in_r,pi,e_PBE_hxc(N_states),ec_PBE(N_states),ex_PBE(N_states),eps_h(N_states),beta(N_states),on_top_two_dm_in_r_mu_corrected
  double precision :: aos_array(ao_num), grad_aos_array(3,ao_num)
  double precision :: rho_a(N_states),rho_b(N_states)
  double precision :: grad_rho_a(3,N_states),grad_rho_b(3,N_states)
  double precision :: grad_rho_a_2(N_states),grad_rho_b_2(N_states),grad_rho_a_b(N_states)
  double precision :: rhoc,rhoo,sigmacc,sigmaco,sigmaoo,vrhoc,vrhoo,vsigmacc,vsigmaco,vsigmaoo
  integer :: m, istate
  pi = 4d0 * datan(1d0)

  eps_hxc_on_top_PBE_SR = 0d0
  call density_and_grad_alpha_beta_and_all_aos_and_grad_aos_at_r(r,rho_a,rho_b,grad_rho_a,grad_rho_b, aos_array, grad_aos_array)
  grad_rho_a_2 = 0.d0
  grad_rho_b_2 = 0.d0
  grad_rho_a_b = 0.d0
  do istate = 1, N_states
   do m = 1, 3
    grad_rho_a_2(istate) += grad_rho_a(m,istate)*grad_rho_a(m,istate)
    grad_rho_b_2(istate) += grad_rho_b(m,istate)*grad_rho_b(m,istate)
    grad_rho_a_b(istate) += grad_rho_a(m,istate)*grad_rho_b(m,istate)
   enddo
  enddo
  do istate = 1, N_states
   ! convertion from (alpha,beta) formalism to (closed, open) formalism
   call rho_ab_to_rho_oc(rho_a(istate),rho_b(istate),rhoo,rhoc)
   call grad_rho_ab_to_grad_rho_oc(grad_rho_a_2(istate),grad_rho_b_2(istate),grad_rho_a_b(istate),sigmaoo,sigmacc,sigmaco)
   call Ec_sr_PBE(mu,rhoc,rhoo,sigmacc,sigmaco,sigmaoo,ec_PBE(istate))
   call Energy_x_pbe_sr(mu,rho_a,rho_b,grad_rho_a_2(istate),grad_rho_b_2(istate),grad_rho_a_b(istate),ex_PBE(istate))
   call give_epsilon_h_SR(r,mu,eps_h)
   e_PBE_hxc(istate)=eps_h(istate)+ec_PBE(istate)+ex_PBE(istate)
   beta(istate) = (2d0*e_PBE_hxc(istate))/(pi*on_top_two_dm_in_r_mu_corrected(mu,r,istate))
   eps_hxc_on_top_PBE_SR(istate)=e_PBE_hxc(istate)/(1d0+beta(istate)*mu**2d0)
  enddo
 end



 BEGIN_PROVIDER [double precision, energy_Hxc_bis, (N_states)]
 implicit none
 integer :: j,k,l,istate
 double precision, allocatable :: r(:), eps_hxc_on_top_PBE(:)
 energy_Hxc_bis = 0.d0
 allocate(r(3),eps_hxc_on_top_PBE(N_states))
  do j = 1, nucl_num
   do k = 1, n_points_radial_grid  -1
    do l = 1, n_points_integration_angular
     r(1) = grid_points_per_atom(1,l,k,j)
     r(2) = grid_points_per_atom(2,l,k,j)
     r(3) = grid_points_per_atom(3,l,k,j)
     do istate = 1, N_states
      call give_epsilon_hxc_PBE_mu_corrected(mu_erf,r,eps_hxc_on_top_PBE) 
      energy_Hxc_bis(istate) += final_weight_functions_at_grid_points(l,k,j) * eps_hxc_on_top_PBE(istate)
     enddo
    enddo
   enddo
  enddo
 deallocate(r,eps_hxc_on_top_PBE)
END_PROVIDER

 BEGIN_PROVIDER [double precision, energy_Hxc_ter, (N_states)]
 implicit none
 integer :: j,k,l,istate
 double precision, allocatable :: r(:), eps_hxc_on_top_PBE_SR(:)
 energy_Hxc_ter = 0.d0
 allocate(r(3),eps_hxc_on_top_PBE_SR(N_states))
  do j = 1, nucl_num
   do k = 1, n_points_radial_grid  -1
    do l = 1, n_points_integration_angular
     r(1) = grid_points_per_atom(1,l,k,j)
     r(2) = grid_points_per_atom(2,l,k,j)
     r(3) = grid_points_per_atom(3,l,k,j)
     do istate = 1, N_states
      call give_epsilon_hxc_PBE_SR_mu_corrected(mu_erf,r,eps_hxc_on_top_PBE_SR)
      energy_Hxc_ter(istate) += final_weight_functions_at_grid_points(l,k,j) * eps_hxc_on_top_PBE_SR(istate)
     enddo
    enddo
   enddo
  enddo
 deallocate(r,eps_hxc_on_top_PBE_SR)
END_PROVIDER

!!!!!!!!!!!new test Julian!!!!!!!!

subroutine give_epsilon_h_xc_PBE_SR_mu_corrected(mu,r,eps_h_xc_on_top_PBE_SR)
  implicit none
  double precision, intent(in)  :: mu , r(3)
  double precision, intent(out) :: eps_h_xc_on_top_PBE_SR(N_states)
  double precision :: two_dm_in_r,pi,e_PBE_xc(N_states),ec_PBE(N_states),ex_PBE(N_states),eps_h(N_states),beta(N_states),on_top_two_dm_in_r_mu_corrected
  double precision :: aos_array(ao_num), grad_aos_array(3,ao_num)
  double precision :: rho_a(N_states),rho_b(N_states)
  double precision :: grad_rho_a(3,N_states),grad_rho_b(3,N_states)
  double precision :: grad_rho_a_2(N_states),grad_rho_b_2(N_states),grad_rho_a_b(N_states)
  double precision :: rhoc,rhoo,sigmacc,sigmaco,sigmaoo,vrhoc,vrhoo,vsigmacc,vsigmaco,vsigmaoo
  integer :: m, istate
  pi = 4d0 * datan(1d0)

  eps_h_xc_on_top_PBE_SR = 0d0
  call density_and_grad_alpha_beta_and_all_aos_and_grad_aos_at_r(r,rho_a,rho_b,grad_rho_a,grad_rho_b,aos_array, grad_aos_array)
  grad_rho_a_2 = 0.d0
  grad_rho_b_2 = 0.d0
  grad_rho_a_b = 0.d0
  do istate = 1, N_states
   do m = 1, 3
    grad_rho_a_2(istate) += grad_rho_a(m,istate)*grad_rho_a(m,istate)
    grad_rho_b_2(istate) += grad_rho_b(m,istate)*grad_rho_b(m,istate)
    grad_rho_a_b(istate) += grad_rho_a(m,istate)*grad_rho_b(m,istate)
   enddo
  enddo
  do istate = 1, N_states
   ! convertion from (alpha,beta) formalism to (closed, open) formalism
   call rho_ab_to_rho_oc(rho_a(istate),rho_b(istate),rhoo,rhoc)
   call grad_rho_ab_to_grad_rho_oc(grad_rho_a_2(istate),grad_rho_b_2(istate),grad_rho_a_b(istate),sigmaoo,sigmacc,sigmaco)
   call Ec_sr_PBE(0d0,rhoc,rhoo,sigmacc,sigmaco,sigmaoo,ec_PBE(istate))
   call Energy_x_pbe_sr(0d0,rho_a,rho_b,grad_rho_a_2(istate),grad_rho_b_2(istate),grad_rho_a_b(istate),ex_PBE(istate))
   call give_epsilon_h_SR(r,mu,eps_h)
   e_PBE_xc(istate)=ec_PBE(istate)+ex_PBE(istate)
   beta(istate) =(2d0*e_PBE_xc(istate))/(pi*(on_top_two_dm_in_r_mu_corrected(mu,r,istate)-(rho_a(istate)+rho_b(istate))**2d0))
   eps_h_xc_on_top_PBE_SR(istate)=eps_h(istate)+(e_PBE_xc(istate)/(1d0+beta(istate)*mu**2d0))
  enddo
 end

 BEGIN_PROVIDER [double precision, energy_Hxc_4, (N_states)]
 implicit none
 integer :: j,k,l,istate
 double precision, allocatable :: r(:), eps_h_xc_on_top_PBE_SR(:)
 energy_Hxc_4 = 0.d0
 allocate(r(3),eps_h_xc_on_top_PBE_SR(N_states))
  do j = 1, nucl_num
   do k = 1, n_points_radial_grid  -1
    do l = 1, n_points_integration_angular
     r(1) = grid_points_per_atom(1,l,k,j)
     r(2) = grid_points_per_atom(2,l,k,j)
     r(3) = grid_points_per_atom(3,l,k,j)
     do istate = 1, N_states
      call give_epsilon_h_xc_PBE_SR_mu_corrected(mu_erf,r,eps_h_xc_on_top_PBE_SR)
      energy_Hxc_4(istate) += final_weight_functions_at_grid_points(l,k,j) * eps_h_xc_on_top_PBE_SR(istate)
     enddo
    enddo
   enddo
  enddo
 deallocate(r,eps_h_xc_on_top_PBE_SR)
END_PROVIDER

subroutine give_epsilon_hx_c_PBE_SR_mu_corrected(mu,r,eps_hx_c_on_top_PBE_SR)
  implicit none
  double precision, intent(in)  :: mu , r(3)
  double precision, intent(out) :: eps_hx_c_on_top_PBE_SR(N_states)
  double precision :: two_dm_in_r,pi,e_PBE_xc(N_states),ec_PBE(N_states),ex_PBE(N_states),eps_h(N_states),beta(N_states),on_top_two_dm_in_r_mu_corrected
  double precision :: aos_array(ao_num), grad_aos_array(3,ao_num)
  double precision :: rho_a(N_states),rho_b(N_states)
  double precision :: grad_rho_a(3,N_states),grad_rho_b(3,N_states)
  double precision :: grad_rho_a_2(N_states),grad_rho_b_2(N_states),grad_rho_a_b(N_states)
  double precision :: rhoc,rhoo,sigmacc,sigmaco,sigmaoo,vrhoc,vrhoo,vsigmacc,vsigmaco,vsigmaoo
  integer :: m, istate
  pi = 4d0 * datan(1d0)

  eps_hx_c_on_top_PBE_SR = 0d0
  call density_and_grad_alpha_beta_and_all_aos_and_grad_aos_at_r(r,rho_a,rho_b,grad_rho_a,grad_rho_b,aos_array,grad_aos_array)
  grad_rho_a_2 = 0.d0
  grad_rho_b_2 = 0.d0
  grad_rho_a_b = 0.d0
  do istate = 1, N_states
   do m = 1, 3
    grad_rho_a_2(istate) += grad_rho_a(m,istate)*grad_rho_a(m,istate)
    grad_rho_b_2(istate) += grad_rho_b(m,istate)*grad_rho_b(m,istate)
    grad_rho_a_b(istate) += grad_rho_a(m,istate)*grad_rho_b(m,istate)
   enddo
  enddo
  do istate = 1, N_states
   ! convertion from (alpha,beta) formalism to (closed, open) formalism
   call rho_ab_to_rho_oc(rho_a(istate),rho_b(istate),rhoo,rhoc)
   call grad_rho_ab_to_grad_rho_oc(grad_rho_a_2(istate),grad_rho_b_2(istate),grad_rho_a_b(istate),sigmaoo,sigmacc,sigmaco)
   call Ec_sr_PBE(0d0,rhoc,rhoo,sigmacc,sigmaco,sigmaoo,ec_PBE(istate))
   call Energy_x_pbe_sr(mu,rho_a,rho_b,grad_rho_a_2(istate),grad_rho_b_2(istate),grad_rho_a_b(istate),ex_PBE(istate))
   call give_epsilon_h_SR(r,mu,eps_h)
   beta(istate)=(2d0*ec_PBE(istate))/(pi*(on_top_two_dm_in_r_mu_corrected(mu,r,istate)-0.5d0*(rho_a(istate)+rho_b(istate))**2d0))
   eps_hx_c_on_top_PBE_SR(istate)=eps_h(istate)+ex_PBE(istate)+(ec_PBE(istate)/(1d0+beta(istate)*mu**2d0))
  enddo
 end





 BEGIN_PROVIDER [double precision, energy_Hxc_5, (N_states)]
 implicit none
 integer :: j,k,l,istate
 double precision, allocatable :: r(:), eps_hx_c_on_top_PBE_SR(:)
 energy_Hxc_5 = 0.d0
 allocate(r(3),eps_hx_c_on_top_PBE_SR(N_states))
  do j = 1, nucl_num
   do k = 1, n_points_radial_grid  -1
    do l = 1, n_points_integration_angular
     r(1) = grid_points_per_atom(1,l,k,j)
     r(2) = grid_points_per_atom(2,l,k,j)
     r(3) = grid_points_per_atom(3,l,k,j)
     do istate = 1, N_states
      call give_epsilon_hx_c_PBE_SR_mu_corrected(mu_erf,r,eps_hx_c_on_top_PBE_SR)
      energy_Hxc_5(istate) += final_weight_functions_at_grid_points(l,k,j)*eps_hx_c_on_top_PBE_SR(istate)
     enddo
    enddo
   enddo
  enddo
 deallocate(r,eps_hx_c_on_top_PBE_SR)
END_PROVIDER
