
 BEGIN_PROVIDER [double precision,Energy_c_md_LDA_mu_of_r_grille_sphe_of_r,(N_states,n_points_radial_grid_spherical)]
&BEGIN_PROVIDER [double precision,Energy_c_md_LDA_mu_of_r_grille_sphe_of_r_integrated,(N_states)]
&BEGIN_PROVIDER [double precision, rho_r,(N_states,n_points_radial_grid_spherical)]
&BEGIN_PROVIDER [double precision,mu_r,(N_states,n_points_radial_grid_spherical)]
&BEGIN_PROVIDER [double precision,rhomu_r,(N_states,n_points_radial_grid_spherical)]
&BEGIN_PROVIDER [double precision,rhomu_r_integrated, (N_states)]
&BEGIN_PROVIDER [double precision,Energy_c_md_LDA_mu_of_r_integrated_core,(N_states)]
&BEGIN_PROVIDER [double precision,Energy_c_md_LDA_mu_of_r_integrated_valence,(N_states)]
include 'Utils/constants.include.F'
 implicit none
 integer :: i,j,istate
 double precision, allocatable :: aos_array(:), r(:), rho_a(:), rho_b(:),ec(:)
 logical :: dospin
 double precision :: r2(3),dr2(3),local_potential,r12,dx2,mu,mu_coulomb,coulomb,two_body_dm
 double precision :: threshold,two_dm,two_dm_laplacian,total_dm,two_dm_HF,two_dm_laplacian_HF,total_dm_HF,dpi
 double precision :: cpu0,cpu1,integral_f,mu_integral
 dpi = 1.5d0 * dsqrt(dacos(-1.d0))
 dospin = .True. ! JT dospin have to be set to true for open shell
 threshold = 1.d-07
 !mu_average_grille_sphe = 0.d0
 allocate(aos_array(ao_num),r(3), rho_a(N_states),rho_b(N_states),ec(N_states))
 call cpu_time(cpu0)
!!!!!!spherical integrale!!!!!
 rhomu_r_integrated = 0.d0
 Energy_c_md_LDA_mu_of_r_integrated_core(:) = 0.d0
 Energy_c_md_LDA_mu_of_r_integrated_valence(:) = 0.d0
 Energy_c_md_LDA_mu_of_r_grille_sphe_of_r_integrated = 0.d0 
do i = 1, n_points_radial_grid_spherical
   Energy_c_md_LDA_mu_of_r_grille_sphe_of_r(:,i) = 0.d0 
   rho_r(:,i) = 0.d0
   mu_r(:,i) = 0.d0
   rhomu_r(:,i) = 0.d0
   do j = 1, n_points_integration_angular
    r(1) = spherical_grid_of_r(1,j,i)
    r(2) = spherical_grid_of_r(2,j,i)
    r(3) = spherical_grid_of_r(3,j,i)

    call dm_dft_alpha_beta_and_all_aos_at_r(r,rho_a,rho_b,aos_array)
    if(mu_of_r_potential.EQ."cusp_condition")then
     istate = 1
     call spherical_averaged_two_dm_at_second_order(r,0.d0,istate,two_dm,two_dm_laplacian,total_dm)
     call spherical_averaged_two_dm_HF_at_second_order(r,0.d0,istate,two_dm_HF,two_dm_laplacian_HF,total_dm_HF)
     two_dm = max(two_dm,1.d-15)
     two_dm_HF = max(two_dm_HF,1.d-15)
     mu =  dpi * (two_dm_laplacian / two_dm - two_dm_laplacian_HF / two_dm_HF)
     mu = max(mu,1.d-15)
    else if(mu_of_r_potential.EQ."integral_hf")then
     call integral_of_f_12_on_hf(r,integral_f)
     mu = mu_integral(integral_f,r)
    else if(mu_of_r_potential.EQ."hf_coallescence")then
     call local_r12_operator_on_hf(r,r,local_potential)
     mu =  local_potential * dsqrt(dacos(-1.d0)) * 0.5d0
    else if(mu_of_r_potential.EQ."psi_coallescence")then
     call expectation_value_in_real_space(r,r,local_potential,two_body_dm)
     mu =  local_potential * dsqrt(dacos(-1.d0)) * 0.5d0
    else 
      print*,'you requested the following mu_of_r_potential'
      print*,mu_of_r_potential
      print*,'which does not correspond to any of the options for such keyword'
      stop
    endif
    if(mu.lt.0.d0)then
     print*,r
     print*,mu
    endif
    do istate = 1, N_states
!!!!!!!!!! CORRELATION PART
     call ESRC_MD_LDAERF (mu,rho_a(istate),rho_b(istate),dospin,ec(istate))
     Energy_c_md_LDA_mu_of_r_grille_sphe_of_r(istate,i) += weights_rad_of_r(j,i) * ec(istate) 
     rho_r(istate,i) += (rho_a(istate)+rho_b(istate)) * weights_rad_of_r(j,i)
     mu_r(istate,i) +=  mu * weights_rad_of_r(j,i) 
     rhomu_r(istate,i) +=  (rho_a(istate)+rho_b(istate)) * mu * weights_rad_of_r(j,i)       
    enddo
  enddo
  do istate = 1, N_states
   rhomu_r(istate,i) = rhomu_r(istate,i) / dble(elec_num)
   mu_r(istate,i) = mu_r(istate,i) /  max(list_r(i)*list_r(i) * 4.d0 * pi,1.d-10)
   rhomu_r_integrated(istate) += rhomu_r(istate,i) * pas_grid_sphe
   Energy_c_md_LDA_mu_of_r_grille_sphe_of_r_integrated(istate) += Energy_c_md_LDA_mu_of_r_grille_sphe_of_r(istate,i) * pas_grid_sphe
   if(list_r(i) < r_core_sphe) then
     Energy_c_md_LDA_mu_of_r_integrated_core(istate) += Energy_c_md_LDA_mu_of_r_grille_sphe_of_r(istate,i) * pas_grid_sphe   
   Else
     Energy_c_md_LDA_mu_of_r_integrated_valence(istate) += Energy_c_md_LDA_mu_of_r_grille_sphe_of_r(istate,i) * pas_grid_sphe 
   endif
      
  enddo
enddo
 deallocate(aos_array,r,rho_a,rho_b, ec)
 call cpu_time(cpu1)
 print*,'Time for the ec_md integration :',cpu1-cpu0


END_PROVIDER

