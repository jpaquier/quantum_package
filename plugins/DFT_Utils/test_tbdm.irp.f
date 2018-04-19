program test
 read_wf = .True.
 touch read_wf
!call correlation_hole
 call normalisation_on_top
end



subroutine correlation_hole
 implicit none
 integer :: i,x1,y1,z1,x2,y2,z2, npas
 double precision :: two_dm_in_r
 double precision :: r1(3),r2(3), integral, theta, dtheta, pi,r

 r=0.5d0
 pi = 3.14159265359d0
 npas = 1000
 dtheta = 360d0/ dfloat(npas)

  theta = pi
  r2(1) = r * cos(theta)
  r2(2) = r * sin(theta)
  r2(3) = 0d0

 do i = 0,npas
  theta = i*dtheta
  theta = theta * 2d0*pi/360d0
  r1(1) = r * cos(theta)
  r1(2) = r * sin(theta)
  r1(3) = 0d0
  print *, theta/pi , two_dm_in_r(r1,r2)
 enddo
 
print*, 'hahah'
end


subroutine normalisation_on_top
 implicit none
 integer :: i,j,k,l 
 double precision :: r(3), rho_a, rho_b, rho, aos_array(ao_num)
 double precision :: two_dm_in_r,dif,tdm

 do j = 1, nucl_num
  do k = 1, n_points_radial_grid  -1
   do l = 1, n_points_integration_angular 
    r(:) = grid_points_per_atom(:,l,k,j)
    call dm_dft_alpha_beta_and_all_aos_at_r(r,rho_a,rho_b,aos_array)
    rho = rho_a + rho_b
    tdm = two_dm_in_r(r,r)
    dif = (rho**2d0)/4d0  - tdm    
!    if( abs(dif) > 1d-8 )then
     print*,(rho**2d0)/4d0,tdm,dif
!    endif
   enddo
  enddo
 enddo
end
