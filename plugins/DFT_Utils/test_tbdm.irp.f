program test
 read_wf = .True.
 touch read_wf
!call correlation_hole
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


