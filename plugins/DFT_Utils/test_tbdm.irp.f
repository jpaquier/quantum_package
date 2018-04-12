program test
 read_wf = .True.
 touch read_wf
!call correlation_hole
 call test_Ec_md_mu_inf_corected
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

subroutine test_Ec_md_mu_inf_corected  
 implicit none 
 double precision :: Ec_md_mu_inf_corected
 print*, "Ec_md_mu_inf_corected =   ", Ec_md_mu_inf_corected(mu_erf)
end
