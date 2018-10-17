
 BEGIN_PROVIDER [integer, n_points_radial_grid_spherical]
&BEGIN_PROVIDER [integer, n_points_total_shperical]
 implicit none
 BEGIN_DOC
! number of radial points per atom for 3d numerical integration, needed for DFT
! for example
 END_DOC
 n_points_radial_grid_spherical= 100
 n_points_total_shperical = n_points_radial_grid_spherical*n_points_integration_angular
END_PROVIDER 

 BEGIN_PROVIDER [double precision, r_max_grid_spherical]
 r_max_grid_spherical=6.d0
 implicit none
 END_PROVIDER 

 BEGIN_PROVIDER [double precision, r_core_sphe]
 r_core_sphe=0.00001d0
 implicit none
 END_PROVIDER


 BEGIN_PROVIDER [double precision, pas_grid_sphe]
 implicit none
 pas_grid_sphe=r_max_grid_spherical/n_points_radial_grid_spherical
 END_PROVIDER 

 BEGIN_PROVIDER [double precision, r0_sphere_grid,(3)]
 implicit none
 r0_sphere_grid(1)=0.0
 r0_sphere_grid(2)=0.0
 r0_sphere_grid(3)=0.0
 END_PROVIDER 

 BEGIN_PROVIDER [double precision,spherical_grid_of_r,(3,n_points_integration_angular,n_points_radial_grid_spherical)]
&BEGIN_PROVIDER [double precision,weights_rad_of_r,(n_points_integration_angular,n_points_radial_grid_spherical)]
&BEGIN_PROVIDER [double precision, list_r ,(n_points_radial_grid_spherical)]
 implicit none
 double precision :: r_localoca
 integer :: i,j
 r_localoca=pas_grid_sphe

 do i = 1,n_points_radial_grid_spherical
  do j=1,n_points_integration_angular
   spherical_grid_of_r(1,j,i)= (r_localoca - r0_sphere_grid(1))* angular_quadrature_points(j,1)
   spherical_grid_of_r(2,j,i)= (r_localoca - r0_sphere_grid(2)) *angular_quadrature_points(j,2)
   spherical_grid_of_r(3,j,i)= (r_localoca - r0_sphere_grid(3)) * angular_quadrature_points(j,3)
   weights_rad_of_r(j,i)= weights_angular_points(j)  * r_localoca * r_localoca
  enddo
  list_r(i) = r_localoca
  r_localoca += pas_grid_sphe
 enddo

 END_PROVIDER

