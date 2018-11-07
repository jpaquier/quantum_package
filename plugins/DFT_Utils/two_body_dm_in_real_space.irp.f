double precision function on_top_two_dm_in_r_mu_corrected(mu,r,istate)
 implicit none
 integer, intent(in) :: istate
 double precision, intent(in) :: r(3),mu
 double precision :: mos_array_r(mo_tot_num),pi
 integer :: i,j,k,l
 double precision :: accu1,accu2,accu3,accu4,threshold
 threshold = 1.d-10
 pi = 4d0 * datan(1d0)
 call give_all_mos_at_r(r,mos_array_r) 

 on_top_two_dm_in_r_mu_corrected = 0.d0
 do i = 1, mo_tot_num
  do j = 1, mo_tot_num
   do k = 1, mo_tot_num
    do l = 1, mo_tot_num
     on_top_two_dm_in_r_mu_corrected += two_bod_alpha_beta_mo_transposed(l,k,j,i,istate) * mos_array_r(i) * mos_array_r(l) * mos_array_r(k) * mos_array_r(j)
    enddo
   enddo
  enddo
 enddo
 on_top_two_dm_in_r_mu_corrected = on_top_two_dm_in_r_mu_corrected / ( 1d0 + 2d0/(dsqrt(pi)*mu) )
 on_top_two_dm_in_r_mu_corrected = max(on_top_two_dm_in_r_mu_corrected,1.d-15)
end



double precision function on_top_two_dm_in_r_mu_corrected_UEG(mu,r,istate)
 implicit none
 integer, intent(in) :: istate
 double precision, intent(in) :: r(3),mu
 double precision :: mos_array_r(mo_tot_num),correction_to_on_top_from_UEG
 integer :: i,j,k,l
 double precision :: accu1,accu2,accu3,accu4,threshold
 threshold = 1.d-10
 call give_all_mos_at_r(r,mos_array_r) 

 on_top_two_dm_in_r_mu_corrected_UEG = 0.d0
 do i = 1, mo_tot_num
  do j = 1, mo_tot_num
   do k = 1, mo_tot_num
    do l = 1, mo_tot_num
     on_top_two_dm_in_r_mu_corrected_UEG += two_bod_alpha_beta_mo_transposed(l,k,j,i,istate) * mos_array_r(i) * mos_array_r(l) * mos_array_r(k) * mos_array_r(j)
    enddo
   enddo
  enddo
 enddo
 on_top_two_dm_in_r_mu_corrected_UEG = on_top_two_dm_in_r_mu_corrected_UEG * correction_to_on_top_from_UEG(mu,r,istate)
 on_top_two_dm_in_r_mu_corrected_UEG = max(on_top_two_dm_in_r_mu_corrected_UEG,1.d-15)
end



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


double precision function correction_to_on_top_from_UEG(mu,r,istate)
 implicit none
 integer, intent(in) :: istate
 double precision, intent(in) :: mu,r(3)
 double precision :: rho_a(N_states),rho_b(N_states)
 double precision :: g0_UEG_mu_inf, g0_UEG_mu 
 call dm_dft_alpha_beta_at_r(r,rho_a,rho_b)

 correction_to_on_top_from_UEG = g0_UEG_mu_inf(rho_a(istate),rho_b(istate)) / g0_UEG_mu(mu,rho_a(istate),rho_b(istate)) 

end




double precision function g0_UEG_mu_inf(rho_a,rho_b)
 implicit none
 double precision, intent(in) :: rho_a,rho_b
 double precision :: rho,pi,x
 double precision :: B, C, D, E, d2, rs, ahd
 rho = rho_a+rho_b
 pi = 4d0 * datan(1d0)
 ahd = -0.36583d0
 d2 = 0.7524d0
 B = -2d0 * ahd - d2
 C = 0.08193d0
 D = -0.01277d0
 E = 0.001859d0
 rs = (3d0 / 4d0*pi*rho)**(1d0/3d0)
 x = -d2*rs

 g0_UEG_mu_inf= 0.5d0 * (1d0- B*rs + C*rs**2d0 + D*rs**3d0 + E*rs**4d0)*exp(x)
end



double precision function g0_UEG_mu(mu,rho_a,rho_b)
 implicit none
 double precision, intent(in) :: rho_a,rho_b,mu
 double precision :: zeta,pi,rho,x,alpha
 double precision :: B, C, D, E, d2, rs, ahd, h, kf
 pi = 4d0 * datan(1d0)
 rho = rho_a+rho_b
!zeta = (rho_a-rho_b)/(rho_a+rho_b)
 alpha = (4d0/(9d0*pi))**(1d0/3d0)
 ahd = -0.36583d0
 d2 = 0.7524d0
 B = -2d0 * ahd - d2
 C = 0.08193d0
 D = -0.01277d0
 E = 0.001859d0
 rs = (3d0 / 4d0*pi*rho)**(1d0/3d0)
 kf = (alpha*rs)**(-1d0)
 zeta = mu / kf
 x = -d2*rs*h(zeta)/ahd 
 g0_UEG_mu = (exp(x)/2d0) * (1d0- B*(h(zeta)/ahd)*rs + C*((h(zeta)**2d0)/(ahd**2d0))*(rs**2d0) + D*((h(zeta)**3d0)/(ahd**3d0))*(rs**3d0) + E*((h(zeta)**4d0)/(ahd**4d0))*(rs**4d0) )
 
end



double precision function h(zeta)
 implicit none
 double precision, intent(in) :: zeta
 double precision :: pi
 double precision :: a1, a2, b1, b2, b3, ahd, alpha
 pi = 4d0 * datan(1d0)
 ahd = -0.36583d0
 alpha = (4d0/(9d0*pi))**(1d0/3d0)
 a1 = -(6d0*alpha/pi)*(1d0-log(2d0))
 b1 = 1.4919d0
 b3 = 1.91528d0
 a2 = ahd * b3
 b2 = (a1 - (b3*alpha/sqrt(pi)))/ahd

 h = (a1*zeta**2d0 + a2*zeta**3d0) / (1d0 + b1*zeta + b2*zeta**2d0 + b3*zeta**3d0)
end










!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

double precision function on_top_dm_integral_with_mu_correction(mu,istate)
 implicit none
 integer, intent(in) :: istate
 double precision, intent(in) :: mu
 double precision :: two_dm_in_r, pi
 double precision :: weight
 integer :: j,k,l
 pi = 4d0 * datan(1d0)
 on_top_dm_integral_with_mu_correction = 0d0

 do j = 1, nucl_num
  do k = 1, n_points_radial_grid  -1
   do l = 1, n_points_integration_angular 
    weight = final_weight_functions_at_grid_points(l,k,j) 
    on_top_dm_integral_with_mu_correction += on_top_of_r(l,k,j,istate) * weight
   enddo
  enddo
 enddo
 on_top_dm_integral_with_mu_correction = 2d0 * on_top_dm_integral_with_mu_correction / ( 1d0 + 2d0/(dsqrt(pi)*mu) )

end


 BEGIN_PROVIDER [double precision, Energy_c_md_on_top, (N_states)]
 BEGIN_DOC
  ! Give the Ec_md energy with a good large mu behaviour in function of the on top pair density.
  ! Ec_md_on_top = (alpha/mu**3) * int n2(r,r) dr  where alpha = (sqrt(2pi)*(-2+sqrt(2)))/(3mu**3) 
 END_DOC
 implicit none 
 integer :: istate
 double precision :: pi,mu
 double precision :: on_top_dm_integral_with_mu_correction 
 mu = mu_erf
 pi = 4d0 * datan(1d0)
 do istate = 1, N_states
  Energy_c_md_on_top(istate) = ((-2d0+sqrt(2d0))*sqrt(2d0*pi)/(3d0*(mu**3)))*on_top_dm_integral_with_mu_correction(mu,istate)
 enddo
 END_PROVIDER


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


subroutine give_epsilon_c_md_on_top_PBE(mu,r,eps_c_md_on_top_PBE)
  implicit none
  double precision, intent(in)  :: mu , r(3)
  double precision, intent(out) :: eps_c_md_on_top_PBE(N_states)
  double precision :: two_dm_in_r, pi, e_pbe(N_states),beta(N_states)
  double precision :: aos_array(ao_num), grad_aos_array(3,ao_num)
  double precision :: rho_a(N_states),rho_b(N_states)
  double precision :: grad_rho_a(3,N_states),grad_rho_b(3,N_states)
  double precision :: grad_rho_a_2(N_states),grad_rho_b_2(N_states),grad_rho_a_b(N_states)
  double precision :: rhoc,rhoo,sigmacc,sigmaco,sigmaoo,vrhoc,vrhoo,vsigmacc,vsigmaco,vsigmaoo
  integer :: m, istate
  pi = 4d0 * datan(1d0)

  eps_c_md_on_top_PBE = 0d0
  call density_and_grad_alpha_beta_and_all_aos_and_grad_aos_at_r(r,rho_a,rho_b, grad_rho_a, grad_rho_b, aos_array, grad_aos_array)
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
   call Ec_sr_PBE(0d0,rhoc,rhoo,sigmacc,sigmaco,sigmaoo,e_PBE(istate))
   beta(istate) = (3d0*e_PBE(istate))/( (-2d0+sqrt(2d0))*sqrt(2d0*pi)*2d0*two_dm_in_r(r,r,istate) )
   eps_c_md_on_top_PBE(istate)=e_PBE(istate)/(1d0+beta(istate)*mu**3d0)
  enddo
 end



subroutine give_epsilon_c_md_on_top_PBE_and_corrected(mu,r,on_top,eps_c_md_on_top_PBE,eps_c_md_on_top_PBE_corrected)
  implicit none
  double precision, intent(in)  :: mu , r(3)
  double precision, intent(in)  :: on_top(N_states)
  double precision, intent(out) :: eps_c_md_on_top_PBE(N_states),eps_c_md_on_top_PBE_corrected(N_states)
  double precision              :: on_top_corrected(N_states)
  double precision :: two_dm_in_r, pi, e_pbe(N_states),beta(N_states)
  double precision :: aos_array(ao_num), grad_aos_array(3,ao_num)
  double precision :: rho_a(N_states),rho_b(N_states)
  double precision :: grad_rho_a(3,N_states),grad_rho_b(3,N_states)
  double precision :: grad_rho_a_2(N_states),grad_rho_b_2(N_states),grad_rho_a_b(N_states)
  double precision :: rhoc,rhoo,sigmacc,sigmaco,sigmaoo,vrhoc,vrhoo,vsigmacc,vsigmaco,vsigmaoo
  integer :: m, istate
  pi = 4d0 * datan(1d0)
  double precision :: on_top_tmp(N_states)
  on_top_tmp = max(on_top,1.d-15)
  eps_c_md_on_top_PBE = 0d0
  eps_c_md_on_top_PBE_corrected = 0d0
  call density_and_grad_alpha_beta_and_all_aos_and_grad_aos_at_r(r,rho_a,rho_b, grad_rho_a, grad_rho_b, aos_array, grad_aos_array)
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
   call Ec_sr_PBE(0d0,rhoc,rhoo,sigmacc,sigmaco,sigmaoo,e_PBE(istate))
   beta(istate) = (3d0*e_PBE(istate))/( (-2d0+sqrt(2d0))*sqrt(2d0*pi)*2d0*on_top_tmp(istate) )
   eps_c_md_on_top_PBE(istate)=e_PBE(istate)/(1d0+beta(istate)*mu**3d0)
   on_top_corrected(istate) = on_top_tmp(istate) / ( 1d0 + 2d0/(dsqrt(pi)*mu) )
   beta(istate) = (3d0*e_PBE(istate))/( (sqrt(2d0)-2d0)*sqrt(2d0*pi)*2d0*on_top_corrected(istate) )
   eps_c_md_on_top_PBE_corrected(istate)=e_PBE(istate)/(1d0+beta(istate)*mu**3d0)
  enddo
 end




 BEGIN_PROVIDER [double precision, Energy_c_md_on_top_PBE, (N_states)]
 BEGIN_DOC
  ! Give the Ec_md energy with a good large mu behaviour in function of the on top pair density coupled to the PBE correlation energy at mu=0
  ! Ec_md_on_top_PBE = Int epsilon_c_PBE_mu=0 / ( 1 + beta*mu**3 ) = Int eps_c_md_on_top_PBE  with beta chosen to recover the good large mu behaviour of the Energy_c_md_on_top functional
 END_DOC
 implicit none
 double precision :: eps_c_md_on_top_PBE(N_states)
 double precision :: two_dm_in_r, r(3)
 double precision :: weight,mu
 integer :: j,k,l,istate
 double precision :: wall1,wall0  
!call cpu_time(wall0)
 mu = mu_erf
 Energy_c_md_on_top_PBE = 0d0
  
 do j = 1, nucl_num
  do k = 1, n_points_radial_grid  -1
   do l = 1, n_points_integration_angular 
    r(:) = grid_points_per_atom(:,l,k,j)
    weight = final_weight_functions_at_grid_points(l,k,j) 
    call give_epsilon_c_md_on_top_PBE(mu,r,eps_c_md_on_top_PBE)
    do istate = 1, N_states
     Energy_c_md_on_top_PBE(istate) += eps_c_md_on_top_PBE(istate) * weight
    enddo
   enddo
  enddo
 enddo
!call cpu_time(wall1)
!print*,'cpu time for Energy_c_md_on_top_PBE       '
!print*,wall1 - wall0
 
 END_PROVIDER


 BEGIN_PROVIDER [double precision, Energy_c_md_on_top_PBE_cycle, (N_states)]
 BEGIN_DOC
  ! Give the Ec_md energy with a good large mu behaviour in function of the on top pair density coupled to the PBE correlation energy at mu=0
  ! Ec_md_on_top_PBE = Int epsilon_c_PBE_mu=0 / ( 1 + beta*mu**3 ) = Int eps_c_md_on_top_PBE  with beta chosen to recover the good large mu behaviour of the Energy_c_md_on_top functional
 END_DOC
 implicit none
 double precision :: eps_c_md_on_top_PBE(N_states)
 double precision :: two_dm_in_r, r(3)
 double precision :: weight,mu
 integer :: j,k,l,istate
 double precision :: dm_a,dm_b
 double precision :: wall1,wall0  
 call cpu_time(wall0)
 mu = mu_erf

 Energy_c_md_on_top_PBE_cycle = 0d0
 do j = 1, nucl_num
  do k = 1, n_points_radial_grid  -1
   do l = 1, n_points_integration_angular 
    r(:) = grid_points_per_atom(:,l,k,j)
    weight = final_weight_functions_at_grid_points(l,k,j) 
    call dm_dft_alpha_beta_at_r(r,dm_a,dm_b)
    if(weight * one_body_dm_mo_alpha_at_grid_points(l,k,j,1).lt.threshold_grid_dft)cycle
    call give_epsilon_c_md_on_top_PBE(mu,r,eps_c_md_on_top_PBE)
    do istate = 1, N_states
     Energy_c_md_on_top_PBE_cycle(istate) += eps_c_md_on_top_PBE(istate) * weight
    enddo
   enddo
  enddo
 enddo

 call cpu_time(wall1)
 print*,'cpu time for Energy_c_md_on_top_PBE_cycle '
 print*,wall1 - wall0
 
 END_PROVIDER


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


subroutine give_epsilon_c_md_on_top_PBE_mu_corrected(mu,r,eps_c_md_on_top_PBE)
  implicit none
  double precision, intent(in)  :: mu , r(3)
  double precision, intent(out) :: eps_c_md_on_top_PBE(N_states)
  double precision :: two_dm_in_r, pi, e_pbe(N_states),beta(N_states),on_top_two_dm_in_r_mu_corrected
  double precision :: aos_array(ao_num), grad_aos_array(3,ao_num)
  double precision :: rho_a(N_states),rho_b(N_states)
  double precision :: grad_rho_a(3,N_states),grad_rho_b(3,N_states)
  double precision :: grad_rho_a_2(N_states),grad_rho_b_2(N_states),grad_rho_a_b(N_states)
  double precision :: rhoc,rhoo,sigmacc,sigmaco,sigmaoo,vrhoc,vrhoo,vsigmacc,vsigmaco,vsigmaoo
  integer :: m, istate
  pi = 4d0 * datan(1d0)

  eps_c_md_on_top_PBE = 0d0
  call density_and_grad_alpha_beta_and_all_aos_and_grad_aos_at_r(r,rho_a,rho_b, grad_rho_a, grad_rho_b, aos_array, grad_aos_array)
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
   call Ec_sr_PBE(0d0,rhoc,rhoo,sigmacc,sigmaco,sigmaoo,e_PBE(istate))
   beta(istate) = (3d0*e_PBE(istate))/( (-2d0+sqrt(2d0))*sqrt(2d0*pi)*2d0*on_top_two_dm_in_r_mu_corrected(mu,r,istate) )
   eps_c_md_on_top_PBE(istate)=e_PBE(istate)/(1d0+beta(istate)*mu**3d0)
  enddo
 end
 

subroutine give_epsilon_c_md_on_top_PBE_mu_corrected_from_two_dm(mu,r,two_dm,eps_c_md_on_top_PBE)
  implicit none
  double precision, intent(in)  :: mu , r(3), two_dm(N_states)
  double precision, intent(out) :: eps_c_md_on_top_PBE(N_states)
  double precision :: two_dm_in_r, pi, e_pbe(N_states),beta(N_states),on_top_two_dm_in_r_mu_corrected_from_two_dm
  double precision :: aos_array(ao_num), grad_aos_array(3,ao_num)
  double precision :: rho_a(N_states),rho_b(N_states)
  double precision :: grad_rho_a(3,N_states),grad_rho_b(3,N_states)
  double precision :: grad_rho_a_2(N_states),grad_rho_b_2(N_states),grad_rho_a_b(N_states)
  double precision :: rhoc,rhoo,sigmacc,sigmaco,sigmaoo,vrhoc,vrhoo,vsigmacc,vsigmaco,vsigmaoo
  integer :: m, istate

  pi = 4.d0 * datan(1.d0)

  eps_c_md_on_top_PBE = 0.d0
  call density_and_grad_alpha_beta_and_all_aos_and_grad_aos_at_r(r,rho_a,rho_b, grad_rho_a, grad_rho_b, aos_array, grad_aos_array)
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
   call Ec_sr_PBE(0.d0,rhoc,rhoo,sigmacc,sigmaco,sigmaoo,e_PBE(istate))
   beta(istate) = (3.d0*e_PBE(istate))/( (-2.d0+sqrt(2d0))*sqrt(2.d0*pi)*2.d0*on_top_two_dm_in_r_mu_corrected_from_two_dm(mu,r,istate,two_dm) )
   eps_c_md_on_top_PBE(istate)=e_PBE(istate)/(1.d0+beta(istate)*mu**3.d0)
  enddo
 end
 
 double precision function on_top_two_dm_in_r_mu_corrected_from_two_dm(mu,r,istate,two_dm)
 implicit none
 double precision, intent(in) :: mu,r(3),two_dm(N_states)
 integer, intent(in) :: istate
 double precision :: pi
 pi = 4.d0 * datan(1.d0)
 on_top_two_dm_in_r_mu_corrected_from_two_dm = two_dm(istate)  / ( 1.d0 + 2.d0/(dsqrt(pi)*mu) )
 on_top_two_dm_in_r_mu_corrected_from_two_dm = max(on_top_two_dm_in_r_mu_corrected_from_two_dm ,1.d-15)
 end

 double precision function on_top_two_dm_in_r_mu_corrected_UEG_from_two_dm(mu,r,istate,two_dm)
 implicit none
 double precision, intent(in) :: mu,r(3),two_dm(N_states)
 integer, intent(in) :: istate
 
 double precision :: correction_to_on_top_from_UEG
 on_top_two_dm_in_r_mu_corrected_UEG_from_two_dm = two_dm(istate) * correction_to_on_top_from_UEG(mu,r,istate)
 on_top_two_dm_in_r_mu_corrected_UEG_from_two_dm = max(on_top_two_dm_in_r_mu_corrected_UEG_from_two_dm,1.d-15)
 end

 BEGIN_PROVIDER [double precision, Energy_c_md_on_top_PBE_mu_corrected, (N_states)]
 BEGIN_DOC
  ! Give the Ec_md energy with a good large mu behaviour in function of the on top pair density with mu correction coupled to the PBE correlation energy at mu=0
  ! Ec_md_on_top_PBE = Int epsilon_c_PBE_mu=0 / ( 1 + beta*mu**3 ) = Int eps_c_md_on_top_PBE  with beta chosen to recover the good large mu behaviour of the Energy_c_md_on_top functional
 END_DOC
 implicit none
 double precision :: eps_c_md_on_top_PBE(N_states)
 double precision :: two_dm_in_r, r(3)
 double precision :: weight,mu
 integer :: j,k,l,istate
 mu = mu_erf
 Energy_c_md_on_top_PBE_mu_corrected = 0d0
  
 do j = 1, nucl_num
  do k = 1, n_points_radial_grid  -1
   do l = 1, n_points_integration_angular 
    r(:) = grid_points_per_atom(:,l,k,j)
    weight = final_weight_functions_at_grid_points(l,k,j) 
    call give_epsilon_c_md_on_top_PBE_mu_corrected(mu,r,eps_c_md_on_top_PBE)
    do istate = 1, N_states
     Energy_c_md_on_top_PBE_mu_corrected(istate) += eps_c_md_on_top_PBE(istate) * weight
    enddo
   enddo
  enddo
 enddo
 END_PROVIDER




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


subroutine give_epsilon_c_md_on_top_PBE_mu_corrected_UEG(mu,r,eps_c_md_on_top_PBE)
  implicit none
  double precision, intent(in)  :: mu , r(3)
  double precision, intent(out) :: eps_c_md_on_top_PBE(N_states)
  double precision :: two_dm_in_r, pi, e_pbe(N_states),beta(N_states),on_top_two_dm_in_r_mu_corrected_UEG
  double precision :: aos_array(ao_num), grad_aos_array(3,ao_num)
  double precision :: rho_a(N_states),rho_b(N_states)
  double precision :: grad_rho_a(3,N_states),grad_rho_b(3,N_states)
  double precision :: grad_rho_a_2(N_states),grad_rho_b_2(N_states),grad_rho_a_b(N_states)
  double precision :: rhoc,rhoo,sigmacc,sigmaco,sigmaoo,vrhoc,vrhoo,vsigmacc,vsigmaco,vsigmaoo
  integer :: m, istate
  pi = 4d0 * datan(1d0)

  eps_c_md_on_top_PBE = 0d0
  call density_and_grad_alpha_beta_and_all_aos_and_grad_aos_at_r(r,rho_a,rho_b, grad_rho_a, grad_rho_b, aos_array, grad_aos_array)
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
   call Ec_sr_PBE(0d0,rhoc,rhoo,sigmacc,sigmaco,sigmaoo,e_PBE(istate))
   beta(istate) = (3d0*e_PBE(istate))/( (-2d0+sqrt(2d0))*sqrt(2d0*pi)*2d0*on_top_two_dm_in_r_mu_corrected_UEG(mu,r,istate) )
   eps_c_md_on_top_PBE(istate)=e_PBE(istate)/(1d0+beta(istate)*mu**3d0)
  enddo
 end
 

subroutine give_epsilon_c_md_on_top_PBE_mu_corrected_UEG_from_two_dm(mu,r,two_dm,eps_c_md_on_top_PBE)
  implicit none
  double precision, intent(in)  :: mu , r(3), two_dm(N_states)
  double precision, intent(out) :: eps_c_md_on_top_PBE(N_states)
  double precision :: two_dm_in_r, pi, e_pbe(N_states),beta(N_states),on_top_two_dm_in_r_mu_corrected_UEG_from_two_dm
  double precision :: aos_array(ao_num), grad_aos_array(3,ao_num)
  double precision :: rho_a(N_states),rho_b(N_states)
  double precision :: grad_rho_a(3,N_states),grad_rho_b(3,N_states)
  double precision :: grad_rho_a_2(N_states),grad_rho_b_2(N_states),grad_rho_a_b(N_states)
  double precision :: rhoc,rhoo,sigmacc,sigmaco,sigmaoo,vrhoc,vrhoo,vsigmacc,vsigmaco,vsigmaoo
  integer :: m, istate
  pi = 4d0 * datan(1d0)

  eps_c_md_on_top_PBE = 0d0
  call density_and_grad_alpha_beta_and_all_aos_and_grad_aos_at_r(r,rho_a,rho_b, grad_rho_a, grad_rho_b, aos_array, grad_aos_array)
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
   call Ec_sr_PBE(0d0,rhoc,rhoo,sigmacc,sigmaco,sigmaoo,e_PBE(istate))
   beta(istate) = (3d0*e_PBE(istate))/( (-2d0+sqrt(2d0))*sqrt(2d0*pi)*2d0*on_top_two_dm_in_r_mu_corrected_UEG_from_two_dm(mu,r,istate,two_dm) )
   eps_c_md_on_top_PBE(istate)=e_PBE(istate)/(1d0+beta(istate)*mu**3d0)
  enddo
 end
 


 BEGIN_PROVIDER [double precision, Energy_c_md_on_top_PBE_mu_corrected_UEG, (N_states)]
 BEGIN_DOC
  ! Give the Ec_md energy with a good large mu behaviour in function of the on top pair density with mu correction based on the on-top of the UEG coupled to the PBE correlation energy at mu=0
  ! Ec_md_on_top_PBE = Int epsilon_c_PBE_mu=0 / ( 1 + beta*mu**3 ) = Int eps_c_md_on_top_PBE  with beta chosen to recover the good large mu behaviour of the Energy_c_md_on_top functional
 END_DOC
 implicit none
 double precision :: two_dm_in_r, r(3)
 double precision :: weight,mu
 double precision,allocatable  :: eps_c_md_on_top_PBE(:)
 allocate(eps_c_md_on_top_PBE(N_states))
 integer :: j,k,l,istate
 mu = mu_erf
 Energy_c_md_on_top_PBE_mu_corrected_UEG = 0d0
  
 do j = 1, nucl_num
  do k = 1, n_points_radial_grid  -1
   do l = 1, n_points_integration_angular 
    r(:) = grid_points_per_atom(:,l,k,j)
    weight = final_weight_functions_at_grid_points(l,k,j) 
    call give_epsilon_c_md_on_top_PBE_mu_corrected_UEG(mu,r,eps_c_md_on_top_PBE)
    do istate = 1, N_states
     Energy_c_md_on_top_PBE_mu_corrected_UEG(istate) += eps_c_md_on_top_PBE(istate) * weight
    enddo
   enddo
  enddo
 enddo
 END_PROVIDER


 BEGIN_PROVIDER [double precision, Energy_c_md_on_top_PBE_mu_UEG_vector, (N_states)]
 BEGIN_DOC
  ! Give the Ec_md energy with a good large mu behaviour in function of the on top pair density with mu correction based on the on-top of the UEG coupled to the PBE correlation energy at mu=0
  ! Ec_md_on_top_PBE = Int epsilon_c_PBE_mu=0 / ( 1 + beta*mu**3 ) = Int eps_c_md_on_top_PBE  with beta chosen to recover the good large mu behaviour of the Energy_c_md_on_top functional
 END_DOC
 implicit none
 double precision ::  r(3)
 double precision :: weight,mu
 integer :: i,istate
 double precision,allocatable  :: eps_c_md_on_top_PBE(:),two_dm(:)
 allocate(eps_c_md_on_top_PBE(N_states),two_dm(N_states))
 mu = mu_erf
 Energy_c_md_on_top_PBE_mu_UEG_vector = 0.d0
  
 do i = 1, n_points_final_grid
  r(1) = final_grid_points(1,i)
  r(2) = final_grid_points(2,i)
  r(3) = final_grid_points(3,i)
  weight=final_weight_functions_at_final_grid_points(i)
  two_dm(:) = on_top_of_r_vector(i,:)
  call give_epsilon_c_md_on_top_PBE_mu_corrected_UEG_from_two_dm(mu,r,two_dm,eps_c_md_on_top_PBE)
  do istate = 1, N_states
   Energy_c_md_on_top_PBE_mu_UEG_vector(istate) += eps_c_md_on_top_PBE(istate) * weight
  enddo
 enddo
 END_PROVIDER

 BEGIN_PROVIDER [double precision, Energy_c_md_on_top_PBE_mu_vector, (N_states)]
 BEGIN_DOC
  ! Give the Ec_md energy with a good large mu behaviour in function of the on top pair density with mu correction based on the on-top of the UEG coupled to the PBE correlation energy at mu=0
  ! Ec_md_on_top_PBE = Int epsilon_c_PBE_mu=0 / ( 1 + beta*mu**3 ) = Int eps_c_md_on_top_PBE  with beta chosen to recover the good large mu behaviour of the Energy_c_md_on_top functional
 END_DOC
 implicit none
 double precision ::  r(3)
 double precision :: weight,mu
 integer :: i,istate
 double precision,allocatable  :: eps_c_md_on_top_PBE(:),two_dm(:)
 allocate(eps_c_md_on_top_PBE(N_states),two_dm(N_states))
 mu = mu_erf
 Energy_c_md_on_top_PBE_mu_vector = 0.d0
  
 do i = 1, n_points_final_grid
  r(1) = final_grid_points(1,i)
  r(2) = final_grid_points(2,i)
  r(3) = final_grid_points(3,i)
  weight=final_weight_functions_at_final_grid_points(i)
  two_dm(:) = on_top_of_r_vector(i,:)
  call give_epsilon_c_md_on_top_PBE_mu_corrected_from_two_dm(mu,r,two_dm,eps_c_md_on_top_PBE)
  do istate = 1, N_states
   Energy_c_md_on_top_PBE_mu_vector(istate) += eps_c_md_on_top_PBE(istate) * weight
  enddo
 enddo
 END_PROVIDER





!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!BONJOUR BONJOUR!!!!!!!!

subroutine give_epsilon_c_on_top_PBE_mu_corrected(mu,r,eps_c_on_top_PBE)
  implicit none
  double precision, intent(in)  :: mu , r(3)
  double precision, intent(out) :: eps_c_on_top_PBE(N_states)
  double precision :: two_dm_in_r, pi,e_pbe(N_states),beta2(N_states),on_top_two_dm_in_r_mu_corrected,on_top_two_dm_in_r_c_mu_corrected,on_top_hf
  double precision :: aos_array(ao_num), grad_aos_array(3,ao_num)
  double precision :: rho_a(N_states),rho_b(N_states)
  double precision :: grad_rho_a(3,N_states),grad_rho_b(3,N_states)
  double precision :: grad_rho_a_2(N_states),grad_rho_b_2(N_states),grad_rho_a_b(N_states)
  double precision ::rhoc,rhoo,sigmacc,sigmaco,sigmaoo,vrhoc,vrhoo,vsigmacc,vsigmaco,vsigmaoo
  double precision :: beta, alpha, on_top_corrected
  integer :: m, istate
  pi = 4d0 * datan(1d0)

  eps_c_on_top_PBE = 0d0
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
   call Ec_sr_PBE(0d0,rhoc,rhoo,sigmacc,sigmaco,sigmaoo,e_PBE(istate))
   on_top_corrected = on_top_two_dm_in_r_mu_corrected(mu,r,istate)
   on_top_two_dm_in_r_c_mu_corrected = on_top_corrected - on_top_hf(r,istate)

   beta  =  (3d0*e_PBE(istate))/( dsqrt(2.d0 * pi) * 4.d0 * on_top_corrected)
   alpha =  (2d0*e_PBE(istate))/( pi * 2.d0 * on_top_two_dm_in_r_c_mu_corrected ) 

!  eps_c_on_top_PBE(istate)=  e_PBE(istate)/(2d0+alpha*mu**2d0) + e_PBE(istate)/(2d0+beta*mu**3d0) 
!  if(mu.lt.1.d0)then
!   print*,mu,rho_a,rho_b 
!  endif
!  if(mu.lt.0.1d0)then
!   eps_c_on_top_PBE(istate) = 0.d0
!   print*,mu,rho_a,rho_b 
!  else
!   eps_c_on_top_PBE(istate)=  pi * on_top_two_dm_in_r_c_mu_corrected / mu**2 + 2.d0/3.d0 * dsqrt(2.d0 * pi) * 2.d0 * on_top_corrected / mu**3
!  endif
   if(mu.lt.0.001d0)then
    eps_c_on_top_PBE(istate) = 0.d0
   else
!   eps_c_on_top_PBE(istate) = erf(mu**1) * ( pi * on_top_two_dm_in_r_c_mu_corrected / mu**2 + 2.d0/3.d0 * dsqrt(2.d0 * pi) * 2.d0 * on_top_corrected / mu**3 )  &  
!                            + (1.d0 - erf(mu**1)) * e_PBE(istate)
!   eps_c_on_top_PBE(istate) = e_PBE(istate)
    eps_c_on_top_PBE(istate) = 4.d0 * dsqrt(pi) * (1.d0 - dsqrt(2.d0))/(3.d0 * mu**3) * on_top_corrected
   endif
   if(eps_c_on_top_PBE(istate).gt.0.d0)then
    print*,''
    print*,mu,rho_a,rho_b 
    print*,on_top_two_dm_in_r_c_mu_corrected,on_top_corrected
   
   endif
!  print*,'erf = ',erf(0.d0),erf(5.d0)
!  eps_c_on_top_PBE(istate)= e_PBE(istate)/(1d0+alpha*mu**2d0)
!  eps_c_on_top_PBE(istate)= e_PBE(istate)/(1d0+beta*mu**3d0)

!  eps_c_on_top_PBE(istate)= 0.5d0 * ( e_PBE(istate)/(1d0+beta*mu**3d0)  + e_PBE(istate)/(1d0+alpha*mu**2d0) )

!  eps_c_on_top_PBE(istate)=e_PBE(istate)*(1d0/(2d0+(3d0*mu**3*e_PBE(istate)/(2d0*sqrt(2d0*pi)*on_top_corrected*2))))
  enddo
 end


double precision function on_top_hf(r,istate)
 implicit none
 integer, intent(in) :: istate
 double precision, intent(in) :: r(3)
 double precision :: mos_array_r(mo_tot_num),rho_alpha,rho_beta
 integer :: i,j
 call give_all_mos_at_r(r,mos_array_r)
 rho_alpha = 0.d0
 rho_beta = 0.d0
 on_top_hf = 0.d0

 do i = 1, elec_alpha_num
  rho_alpha += mos_array_r(i) **2 
 enddo

 do i = 1, elec_beta_num
  rho_beta += mos_array_r(i) **2 
 enddo

 on_top_hf= rho_beta*rho_alpha
end

 BEGIN_PROVIDER [double precision, Energy_c_on_top_PBE_mu_corrected,(N_states)]
 BEGIN_DOC
  ! Ecrire doc 
 END_DOC
 implicit none
 double precision :: eps_c_on_top_PBE(N_states)
 double precision :: two_dm_in_r, r(3)
 double precision :: weight,mu
 integer :: j,k,l,istate
 mu = mu_erf
 Energy_c_on_top_PBE_mu_corrected = 0d0

 do j = 1, nucl_num
  do k = 1, n_points_radial_grid  -1
   do l = 1, n_points_integration_angular
    r(:) = grid_points_per_atom(:,l,k,j)
    weight = final_weight_functions_at_grid_points(l,k,j)
    call give_epsilon_c_on_top_PBE_mu_corrected(mu,r,eps_c_on_top_PBE)
    do istate = 1, N_states
     Energy_c_on_top_PBE_mu_corrected(istate) += eps_c_on_top_PBE(istate) * weight
    enddo
   enddo
  enddo
 enddo
 END_PROVIDER
