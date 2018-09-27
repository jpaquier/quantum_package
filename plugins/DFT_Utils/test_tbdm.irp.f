program test
 read_wf = .True.
 touch read_wf
!call correlation_hole_bis
!call correlation_hole
!call test_mu
 call test_corr
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

! r2(1) = 0.d0
! r2(2) = 0.d0
! r2(3) = 0.d0

 do i = 0,npas
  theta = i*dtheta
  theta = theta * 2d0*pi/360d0
  r1(1) = r * dcos(theta)
  r1(2) = r * dsin(theta)
  r1(3) = 0.d0
  tbdm = two_dm_in_r(r1,r2,1)
  double precision :: r12,two_dm,two_dm_laplacian,total_dm
  r12 = dsqrt((r1(1) - r2(1))**2 + (r1(2) - r2(2))**2 + (r1(3) - r2(3))**2)
  call spherical_averaged_two_dm_at_second_order(r2,r12,1,two_dm,two_dm_laplacian,total_dm)
  print*,two_dm,two_dm_in_r(r2,r2,1)
  write(33,'(100(F16.10,X))')theta , tbdm, two_dm, two_dm_laplacian, total_dm
 enddo
 
print*, 'hahah'
end

subroutine test_mu
 implicit none
 provide on_top_of_r_vector_parallel


end


subroutine correlation_hole_bis
 implicit none
 integer :: i,x1,y1,z1,x2,y2,z2, npas
 double precision :: two_dm_in_r
 double precision :: r1(3),r2(3), integral, theta, dtheta, pi,r,tbdm,dr
 double precision :: r12,two_dm,two_dm_laplacian,total_dm
 npas = 100 
 dr = 3.d0/dfloat(npas)

 r2(1) = 0.d0
 r2(2) = 0.d0
 r2(3) = 0.d0
 r1 = 0.d0
 r = 0.d0
 do i = 0,npas
  r += dr
  r1(1) = r
  tbdm = two_dm_in_r(r1,r2,1)
  r12 = dsqrt((r1(1) - r2(1))**2 + (r1(2) - r2(2))**2 + (r1(3) - r2(3))**2)
  call spherical_averaged_two_dm_at_second_order(r2,r12,1,two_dm,two_dm_laplacian,total_dm)
  print*,two_dm,two_dm_in_r(r2,r2,1)
  write(33,'(100(F16.10,X))')r12 , tbdm, two_dm, two_dm_laplacian, total_dm
 enddo
 
print*, 'hahah'
end

subroutine test_corr
 implicit none
 integer :: j,k,l 
 double precision :: r(3),two_dm,two_dm_laplacian,total_dm
 j = 1
 print*,'n_points_radial_grid = ',n_points_radial_grid
  do k = 1, n_points_radial_grid  -1
   l = 1
   r(1) = grid_points_per_atom(1,l,k,j)
   r(2) = grid_points_per_atom(2,l,k,j)
   r(3) = grid_points_per_atom(3,l,k,j)
   call spherical_averaged_two_dm_at_second_order(r,0.d0,1,two_dm,two_dm_laplacian,total_dm)
  !if(two_dm_laplacian.lt.0.d0)then
    print*,k
    print*,'r = '
    print*,r
    print*,two_dm,two_dm_laplacian
  !endif
  enddo




end
