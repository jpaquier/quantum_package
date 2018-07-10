program test
 read_wf = .True.
 touch read_wf
!call densitymap
!call correlation_hole
 call normalisation_on_top
!call print_weight
!call test_tbdm_with_symmetry
end

subroutine densitymap
 implicit none
 double precision ::r(3),accu(n_points_radial_grid-1),w,nr(n_points_radial_grid-1)
 double precision :: rho_a, rho_b, rho,  aos_array(ao_num)
 integer :: l,k
 accu = 0d0

  do k = 1, n_points_radial_grid  -1
   do l = 1, n_points_integration_angular 
    r(:) = grid_points_per_atom(:,l,k,1)
    call dm_dft_alpha_beta_and_all_aos_at_r(r,rho_a,rho_b,aos_array)
    w =  final_weight_functions_at_grid_points(l,k,1)
    rho = rho_a + rho_b
    accu(k) += w*rho
   enddo
   nr(k) = sqrt(r(1)**2+r(2)**2+r(3)**2)
  enddo

  
  do k = 1, n_points_radial_grid  -1
   print*, nr(k), accu(k)
  enddo
end
 
 


subroutine correlation_hole
 implicit none
 integer :: i,x1,y1,z1,x2,y2,z2, npas
 double precision :: two_dm_in_r
 double precision :: r1(3),r2(3), integral, theta, dtheta, pi,r,tbdm
 r=0.85d0
 pi = 3.14159265359d0
 npas = 100 
 dtheta = 360d0/ dfloat(npas)

  theta = pi
  r2(1) = r * dcos(theta)
  r2(2) = r * dsin(theta)
  r2(3) = 0d0

 do i = 0,npas
  theta = i*dtheta
  theta = theta * 2d0*pi/360d0
  r1(1) = r * dcos(theta)
  r1(2) = r * dsin(theta)
  r1(3) = 0.d0
  tbdm = two_dm_in_r(r1,r2,1)
  print*, theta/pi , tbdm
 enddo
 
print*, 'hahah'
end


