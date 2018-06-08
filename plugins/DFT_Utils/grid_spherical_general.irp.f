
BEGIN_PROVIDER [integer, n_points_radial_grid_spherical]
&BEGIN_PROVIDER [integer, n_points_total_shperical]
 implicit none
 BEGIN_DOC
! number of radial points per atom for 3d numerical integration, needed for DFT
! for example
 END_DOC
 n_points_radial_grid_spherical= 10000
 n_points_total_shperical = n_points_radial_grid_spherical*n_points_integration_angular
END_PROVIDER 



 BEGIN_PROVIDER [double precision,spherical_grid_of_r,(3,n_points_integration_angular,n_points_radial_grid_spherical)]
&BEGIN_PROVIDER [double precision,weights_rad_of_r,(n_points_integration_angular,n_points_radial_grid_spherical)]
&BEGIN_PROVIDER [double precision, r0_sphere_grid,(3)]
&BEGIN_PROVIDER [double precision, pas_grid_sphe]
 implicit none
  double precision :: r_max, r_localoca
  integer :: i,j
 r_max=4
 pas_grid_sphe=r_max/n_points_radial_grid_spherical
 r_localoca=0
 r0_sphere_grid(1)=0.0
 r0_sphere_grid(2)=0.0
 r0_sphere_grid(3)=0.0

 do i = 1,n_points_radial_grid_spherical
  do j=1,n_points_integration_angular
   spherical_grid_of_r(1,j,i)= (r_localoca - r0_sphere_grid(1))* angular_quadrature_points(j,1)
   spherical_grid_of_r(2,j,i)= (r_localoca - r0_sphere_grid(2)) *angular_quadrature_points(j,2)
   spherical_grid_of_r(3,j,i)= (r_localoca - r0_sphere_grid(3)) * angular_quadrature_points(j,3)

   weights_rad_of_r(j,i)= weights_angular_points(j)  * r_localoca * r_localoca
  enddo
  
  r_localoca += pas_grid_sphe
 enddo

 END_PROVIDER


 BEGIN_PROVIDER [double precision, Energy_c_md_LDA_mu_of_r_grille_sphe, (N_states)]
&BEGIN_PROVIDER [double precision, mu_average_grille_sphe, (N_states)]
 implicit none
 BEGIN_DOC
 ! Corelation energy for the multi determinent short range LDA. PRB 73 155111
 ! 2006
 END_DOC
 integer :: j,i,istate
 double precision, allocatable :: aos_array(:), r(:), rho_a(:), rho_b(:), ec(:)
 logical :: dospin
 double precision :: r2(3),dr2(3),local_potential,r12,dx2,mu,mu_coulomb,coulomb,two_body_dm
 double precision :: threshold
 double precision :: cpu0,cpu1
 dospin = .True. ! JT dospin have to be set to true for open shell
 threshold = 1.d-07
 Energy_c_md_LDA_mu_of_r_grille_sphe = 0.d0
 mu_average_grille_sphe = 0.d0
 allocate(aos_array(ao_num),r(3), rho_a(N_states), rho_b(N_states),ec(N_states))
 call cpu_time(cpu0)
 do i = 1, n_points_radial_grid_spherical
   do j = 1, n_points_integration_angular 
    r(1) = spherical_grid_of_r(1,j,i)
    r(2) = spherical_grid_of_r(2,j,i)
    r(3) = spherical_grid_of_r(3,j,i)
    call dm_dft_alpha_beta_and_all_aos_at_r(r,rho_a,rho_b,aos_array)
    if(dabs(weights_rad_of_r(j,i) *pas_grid_sphe * (rho_a(1)+rho_b(1))).lt.threshold)cycle
    if(basis_set_hf_potential)then
     call local_r12_operator_on_hf(r,r,local_potential)
    else
     call expectation_value_in_real_space(r,r,local_potential,two_body_dm)
    endif
    mu = local_potential * dsqrt(dacos(-1.d0)) * 0.5d0
    do istate = 1, N_states
!!!!!!!!! CORRELATION PART
     call ESRC_MD_LDAERF (mu,rho_a(istate),rho_b(istate),dospin,ec(istate))
     Energy_c_md_LDA_mu_of_r_grille_sphe(istate) += weights_rad_of_r(j,i) * pas_grid_sphe * ec(istate)
     mu_average_grille_sphe(istate) += weights_rad_of_r(j,i) * pas_grid_sphe * mu
    enddo
   enddo
 enddo
 deallocate(aos_array,r,rho_a,rho_b, ec)
 call cpu_time(cpu1)
 print*,'Time for the ec_md integration :',cpu1-cpu0
END_PROVIDER

 BEGIN_PROVIDER [double precision,Energy_c_md_LDA_mu_of_r_grille_sphe_of_r,(N_states,n_points_radial_grid_spherical)]
&BEGIN_PROVIDER [double precision, rho_r,(N_states,n_points_radial_grid_spherical)]
&BEGIN_PROVIDER [double precision,mu_r,(N_states,n_points_radial_grid_spherical)]
&BEGIN_PROVIDER [double precision,rhomu_r,(N_states,n_points_radial_grid_spherical)]
include 'Utils/constants.include.F'
 implicit none
 integer :: i,j,istate
 double precision, allocatable :: aos_array(:), r(:), rho_a(:), rho_b(:),ec(:),list_r(:)
 logical :: dospin
 double precision :: r2(3),dr2(3),local_potential,r12,dx2,mu,mu_coulomb,coulomb,two_body_dm
 double precision :: threshold
 double precision :: cpu0,cpu1
 dospin = .True. ! JT dospin have to be set to true for open shell
 threshold = 1.d-07
 !mu_average_grille_sphe = 0.d0
 allocate(aos_array(ao_num),r(3), rho_a(N_states),rho_b(N_states),ec(N_states),list_r(n_points_radial_grid_spherical))
 call cpu_time(cpu0)
!!!!!!spherical integrale!!!!!

do i = 1, n_points_radial_grid_spherical
    
   do j = 1, n_points_integration_angular
    r(1) = spherical_grid_of_r(1,j,i)
    r(2) = spherical_grid_of_r(2,j,i)
    r(3) = spherical_grid_of_r(3,j,i)

    call dm_dft_alpha_beta_and_all_aos_at_r(r,rho_a,rho_b,aos_array)
    if(dabs(weights_rad_of_r(j,i) * (rho_a(1)+rho_b(1))).lt.threshold)cycle
    if(basis_set_hf_potential)then
     call local_r12_operator_on_hf(r,r,local_potential)
    else
     call expectation_value_in_real_space(r,r,local_potential,two_body_dm)
    endif
    mu = local_potential * dsqrt(dacos(-1.d0)) * 0.5d0
    do istate = 1, N_states
!!!!!!!!!! CORRELATION PART
     call ESRC_MD_LDAERF (mu,rho_a(istate),rho_b(istate),dospin,ec(istate))
     Energy_c_md_LDA_mu_of_r_grille_sphe_of_r(istate,i) += weights_rad_of_r(j,i) * ec(istate) 
     rho_r(istate,i) += (rho_a(istate)+rho_b(istate)) * weights_rad_of_r(j,i)
     mu_r(istate,i) +=  mu * weights_rad_of_r(j,i) 
     rhomu_r(istate,i) +=  (rho_a(istate)+rho_b(istate)) * mu * weights_rad_of_r(j,i)       
    enddo
  enddo
  list_r(i)= sqrt(r(1)*r(1)+r(2)*r(2)+r(3)*r(3))
  do istate = 1, N_states
   rhomu_r(istate,i) = rhomu_r(istate,i) / dble(elec_num)
   mu_r(istate,i) = mu_r(istate,i) /  (list_r(i)*list_r(i) * 4.d0 * pi)
  enddo
enddo
 deallocate(aos_array,r,rho_a,rho_b, ec)
 call cpu_time(cpu1)
 print*,'Time for the ec_md integration :',cpu1-cpu0


END_PROVIDER

