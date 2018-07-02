program test
 read_wf = .True.
 touch read_wf
!call densitymap
!call correlation_hole
 call normalisation_on_top
!call print_weight
!call test_opt_on_top
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
 double precision :: two_dm_in_r_new_cycle
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
  tbdm = two_dm_in_r_new_cycle(r1,r2,1)
  print*, theta/pi , tbdm
 enddo
 
print*, 'hahah'
end


subroutine normalisation_on_top
 implicit none
 integer :: i,j,k,l 
 double precision :: r(3), rho_a, rho_b, rho, aos_array(ao_num)
 double precision :: two_dm_in_r,dif,tdm
 double precision :: two_dm_in_r_new_cycle

 do j = 1, nucl_num
  do k = 1, n_points_radial_grid  -1
   do l = 1, n_points_integration_angular 
    r(:) = grid_points_per_atom(:,l,k,j)
    call dm_dft_alpha_beta_and_all_aos_at_r(r,rho_a,rho_b,aos_array)
    rho = rho_a + rho_b
    tdm = two_dm_in_r_new_cycle(r,r,1)
    dif = (rho**2d0)/4d0  - tdm    
!    if( abs(dif) > 1d-8 )then
     print*,(rho**2d0)/4d0,tdm,dif
!    endif
   enddo
  enddo
 enddo
end

subroutine print_weight
 implicit none
 integer :: i,j,k
 double precision :: r(3),distance
 double precision :: dm_a(N_states),dm_b(N_states)
 double precision :: accu_dm(N_states),accu_weight
 do i = 1, n_points_radial_grid - 1
  accu_dm = 0.d0
  accu_weight = 0.d0
  do j = 1, n_points_integration_angular
   r(:) = grid_points_per_atom(:,1,i,1)
   distance = 0.d0
   do k = 1, 3
    distance += r(k)**2
   enddo
   distance = dsqrt(distance)
   call dm_dft_alpha_beta_at_r(r,dm_a,dm_b)
   accu_weight += final_weight_functions_at_grid_points(j,i,1) 
   accu_dm += dm_a * final_weight_functions_at_grid_points(j,i,1)
  enddo
  write(33,*) distance, accu_dm,accu_weight, final_weight_functions_at_grid_points(1,i,1), dm_a+dm_b
 enddo

end

subroutine test_opt_on_top
 implicit none
 integer :: i,j,k,l
 double precision :: two_dm_in_r, on_top_two_dm_in_r_with_symmetry , on_top_two_dm_in_r, on_top_two_dm_in_r_new, r(3), wall1, wall2, accu1, accu2, accu3, weight
 double precision :: on_top_two_dm_in_r_mu_corrected,on_top_two_dm_in_r_mu_corrected_UEG
 double precision :: accu4,two_dm_in_r_new
 double precision :: accu5,on_top_two_dm_in_r_sym,accu6


 accu3 = 0d0
 call cpu_time(wall1)
 do j = 1, nucl_num
  do k = 1, n_points_radial_grid  -1
   do l = 1, n_points_integration_angular 
    r(:) = grid_points_per_atom(:,l,k,j)
    weight = final_weight_functions_at_grid_points(l,k,j)
    accu3 +=  max(two_dm_in_r(r,r,1) * weight,1d-20)
   enddo
  enddo
 enddo
 call cpu_time(wall2)
 print*, "cpu time tbdm :", wall2-wall1 


!accu4 = 0d0
!call cpu_time(wall1)
!do j = 1, nucl_num
! do k = 1, n_points_radial_grid  -1
!  do l = 1, n_points_integration_angular 
!   r(:) = grid_points_per_atom(:,l,k,j)
!   weight = final_weight_functions_at_grid_points(l,k,j)
!   accu4 +=  max(two_dm_in_r_new(r,r,1) * weight,1d-20)
!  enddo
! enddo
!enddo
!call cpu_time(wall2)
!print*, "cpu time tbdm new :", wall2-wall1 




!accu6 = 0d0
!call cpu_time(wall1)
!do j = 1, nucl_num
! do k = 1, n_points_radial_grid  -1
!  do l = 1, n_points_integration_angular 
!   r(:) = grid_points_per_atom(:,l,k,j)
!   weight = final_weight_functions_at_grid_points(l,k,j)
!   accu6 += max(on_top_two_dm_in_r_new(r,1) * weight,1d-20)
!  enddo
! enddo
!enddo
!call cpu_time(wall2)
!print*, "cpu time new without symmetry with cycle :", wall2-wall1 


!accu5 = 0d0
!call cpu_time(wall1)
!do j = 1, nucl_num
! do k = 1, n_points_radial_grid  -1
!  do l = 1, n_points_integration_angular 
!   r(:) = grid_points_per_atom(:,l,k,j)
!   weight = final_weight_functions_at_grid_points(l,k,j)
!   accu5 +=  max(on_top_two_dm_in_r_sym(r,1) * weight,1d-20)
!  enddo
! enddo
!enddo
!call cpu_time(wall2)
!print*, "cpu time tbdm new with symmetry and cycle :", wall2-wall1 


!accu1 = 0d0
!call cpu_time(wall1)
!do j = 1, nucl_num
! do k = 1, n_points_radial_grid  -1
!  do l = 1, n_points_integration_angular 
!   r(:) = grid_points_per_atom(:,l,k,j)
!   weight = final_weight_functions_at_grid_points(l,k,j)
!   accu1 += max(on_top_two_dm_in_r(r,1) * weight,1d-20)
!  enddo
! enddo
!enddo
!call cpu_time(wall2)
!print*, "cpu time without symmetry with cycle :", wall2-wall1 


 accu2 = 0d0
 call cpu_time(wall1)
 do j = 1, nucl_num
  do k = 1, n_points_radial_grid  -1
   do l = 1, n_points_integration_angular 
    r(:) = grid_points_per_atom(:,l,k,j)
    weight = final_weight_functions_at_grid_points(l,k,j)
    accu2 += max(on_top_two_dm_in_r_mu_corrected(mu_erf,r,1) * weight,1d-20)
   enddo
  enddo
 enddo
 call cpu_time(wall2)
 print*, "cpu time with symmetry :", wall2-wall1 


 accu4 = 0d0
 call cpu_time(wall1)
 do j = 1, nucl_num
  do k = 1, n_points_radial_grid  -1
   do l = 1, n_points_integration_angular 
    r(:) = grid_points_per_atom(:,l,k,j)
    weight = final_weight_functions_at_grid_points(l,k,j)
    accu4 += max(on_top_two_dm_in_r_mu_corrected_UEG(mu_erf,r,1) * weight,1d-20)
   enddo
  enddo
 enddo
 call cpu_time(wall2)
 print*, "cpu time with symmetry :", wall2-wall1 



 print*, 'ref : ',accu3
!print*, accu1
!print*, accu6
 print*, 'norm: ',accu2
 print*, 'ueg : ',accu4
 print*, 'n+r : ',accu2-accu3 
 print*, 'u+r : ',accu4-accu3 
 print*, 'u+n : ',accu4-accu2
!print*, accu5
!print*, accu3-accu1 ,accu3-accu6 ,accu3-accu2 ,accu3-accu5
    
end

       

subroutine test_tbdm_with_symmetry
 implicit none
 integer :: i,j,k,l
 double precision :: two_dm_in_r, on_top_two_dm_in_r, on_top_two_dm_in_r_with_symmetry,r(3),on_top_two_dm_in_r_sym
 double precision :: pouet , toto, tata
 
 do j = 1, nucl_num
  do k = 1, n_points_radial_grid  -1
   do l = 1, n_points_integration_angular 
    r(:) = grid_points_per_atom(:,l,k,j)
  ! pouet = two_dm_in_r(r,r,1)
  ! toto = on_top_two_dm_in_r(r,1)
  ! tata = on_top_two_dm_in_r_with_symmetry(r,1)
    tata = on_top_two_dm_in_r_sym(r,1)
  ! if (abs(toto-tata)>1d-10) then
     print*,"pts :", j,k,l
     print*,'r :' , r 
     print*,'on top : ',tata
     print*," "
!   endif
   enddo
  enddo
 enddo

end






















