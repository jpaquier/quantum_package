program pouet
 read_wf = .True.
 touch read_wf
 integer :: nx
 double precision :: rmax
 nx = 100
 rmax = 3.d0
 call test_mu_erf
!call test_erf_coulomb_oprerator
!call test
!call test_inv_erf
!call test_grad
!call test_r12_psi(nx,rmax)
!call test_one_dm_mo
!call routine3
!call test_rho2
!call test_one_dm_ao
!call test_coulomb_oprerator
!call test_expectation_value
!call test_erf_coulomb_oprerator
!call test_nuclear_coulomb_oprerator
!call test_integratio_mo
!call test_naive_grid
end

subroutine test
implicit none
 print*,'e_c = ',energy_c
 print*,'e_x = ',energy_x
end

subroutine routine3
 implicit none
 integer :: i,j,k,l
 double precision :: accu
 accu = 0.d0
 print*, 'energy_x      =',energy_x
 print*, 'energy_c      =',energy_c
 print*, 'energy_c_md   =',energy_c_md
end

subroutine test_rho2
 implicit none
 integer :: j,k,l 
 double precision ::r(3)
 double precision :: test,integral_f

 test = 0.d0
 do j = 1, nucl_num
  do k = 1, n_points_radial_grid  -1
   do l = 1, n_points_integration_angular 
    
    r(1) = grid_points_per_atom(1,l,k,j)
    r(2) = grid_points_per_atom(2,l,k,j)
    r(3) = grid_points_per_atom(3,l,k,j)
    call integral_of_f_12_on_hf(r,integral_f)
    test += integral_f * final_weight_functions_at_grid_points(l,k,j) 
    enddo
   enddo
  enddo
 print*,'test = ',test
! print*,'test Barth = ',test_bart

end




subroutine test_one_dm_mo
 implicit none
 double precision :: r(3)
 double precision, allocatable :: mos_array(:), aos_array(:), aos_array_bis(:), mos_array_bis(:)
 allocate(mos_array(mo_tot_num),aos_array(ao_num),aos_array_bis(ao_num),mos_array_bis(mo_tot_num))
 integer :: m,n,p
 integer :: j,k,l 
 double precision :: ao,mo
 ao = 0.d0
 mo = 0.d0
  do j = 1, nucl_num
   do k = 1, n_points_radial_grid  -1
    do l = 1, n_points_integration_angular 
     r(1) = grid_points_per_atom(1,l,k,j)
     r(2) = grid_points_per_atom(2,l,k,j)
     r(3) = grid_points_per_atom(3,l,k,j)
     call give_all_mos_at_r(r,mos_array) 
     call give_all_aos_at_r(r,aos_array)
     do m = 1, mo_tot_num
      do n = 1, mo_tot_num
       mo += mos_array(n) * mos_array(m) * (one_body_dm_mo_alpha_average(n,m)+one_body_dm_mo_beta_average(n,m)) & 
                                         * final_weight_functions_at_grid_points(l,k,j)
      enddo
     enddo
     enddo
    enddo
   enddo
  print*,'mo = ',mo
end

subroutine test_integratio_mo
 implicit none
 double precision :: r(3)
 double precision, allocatable :: mos_array(:), aos_array(:), aos_array_bis(:), mos_array_bis(:)
 allocate(mos_array(mo_tot_num),aos_array(ao_num),aos_array_bis(ao_num),mos_array_bis(mo_tot_num))
 integer :: m,n,p
 integer :: j,k,l 
  mos_array_bis = 0.d0
  do j = 1, nucl_num
   do k = 1, n_points_radial_grid  -1
!   print*,k
    do l = 1, n_points_integration_angular 
     r(1) = grid_points_per_atom(1,l,k,j)
     r(2) = grid_points_per_atom(2,l,k,j)
     r(3) = grid_points_per_atom(3,l,k,j)
     call give_all_mos_at_r(r,mos_array) 
     do m = 1, mo_tot_num
      mos_array_bis(m) += mos_array(m) * final_weight_functions_at_grid_points(l,k,j)
     enddo
    enddo
   enddo
  enddo
  do m = 1, mo_tot_num
   if(dabs(mos_array_bis(m) - integrals_mo_over_r(m)).gt.1.d-10)then
    print*,m
    print*,dabs(mos_array_bis(m) - integrals_mo_over_r(m)),mos_array_bis(m) , integrals_mo_over_r(m)
   endif
  enddo
end

subroutine test_r12(nx,rmax)
 implicit none
 double precision, intent(in) :: rmax
 integer, intent(in) :: nx
 double precision :: r1(3),r2(3),dr(3),dx,r12
 character*(128) :: output
 character*(128) :: filename
 integer :: i_unit_output,getUnitAndOpen

 output=trim(ezfio_filename)//'.r12_delta'
 output=trim(output)
 print*,'output = ',trim(output)
 i_unit_output = getUnitAndOpen(output,'w')
 dx = rmax/dble(nx)
 dr = 0.d0
 dr(1) = dx
 dr(2) = dx

 r1 = 0.d0
 r2 = r1
 double precision :: integral
 integer :: i
 do i = 1, nx 
  r2 += dr
  r12 = dsqrt((r1(1)-r2(1))**2 + (r1(2)-r2(2))**2 +(r1(3)-r2(3))**2 )
  call numerical_delta_function_coulomb(r1,r2,integral)
  write(i_unit_output,'(100(F16.10,X))')r12,1.d0/r12,integral
 enddo
end

subroutine test_inv_erf
 implicit none
 integer :: ny
 double precision :: y,dy,ymax, inverse_erf,x,thr
 thr = 1.d-10
 ny = 1000
 ymax = 1.d0
 dy = dble(ymax/ny) 
 integer :: i
 y = 0.d0
 do i = 1, ny
  x = inverse_erf(y,thr)
  write(33,*)y,x,erf(x) 
  y+=dy
 enddo


end

subroutine test_r12_psi(nx,rmax)
  include 'constants.include.F'
 implicit none
 double precision, intent(in) :: rmax
 integer, intent(in) :: nx
 double precision :: r1(3),r2(3),dr1(3),dr2(3),dx,r12
 character*(128) :: output
 character*(128) :: filename
 integer :: i_unit_output,getUnitAndOpen

 output=trim(ezfio_filename)//'.r12_delta_psi'
 output=trim(output)
 print*,'output = ',trim(output)
 i_unit_output = getUnitAndOpen(output,'w')
 dx = rmax/dble(nx)
 dr1 = 0.d0
 dr2 = 0.d0
 dr1(1) = dx
 dr2(2) = dx
!dr(2) = dx

!r1(1) = 0.5d0
 double precision :: integral
 double precision :: integral_erf
 integer :: i,j
 double precision :: mos_array_r1(mo_tot_num)
 double precision :: mos_array_r2(mo_tot_num)
 double precision :: dtheta,r12max,dr12,theta,r12_test
 integer :: ntheta,nr12
 ntheta = 10
 r12max = 5.d0 
 nr12 = 500
 dr12 = r12max/dble(nr12) 
 dtheta = 2.d0 * pi /dble(ntheta)

 r1 = 0.d0
 r1(1) = 0.8d0
 character*(128) :: filename_theta
 character*(128) :: output_array(1000)
 do i = 1, ntheta
  if (i.lt.10)then
   write (filename_theta, "(I1)")i
  else
   write (filename_theta, "(I2)")i
  endif
  print*,filename_theta
  output_array(i)=trim(ezfio_filename)//'.'//trim(filename_theta)//'.r12_theta'
  output_array(i)=trim(output_array(i))
  print*,'output = ',trim(output_array(i))
  i_unit_output_array(i) = getUnitAndOpen(output_array(i),'w')
 enddo
 integer :: i_unit_output_array(1000)

 theta = 0.d0
 do i = 1, ntheta
  r12 = 0.d0
  do j = 1, nr12
   r12 += dr12
   r2 = r1
   r2(1) += r12 * dcos(theta)
   r2(2) += r12 * dsin(theta)
   r12_test = dsqrt((r1(1)-r2(1))**2 + (r1(2)-r2(2))**2 +(r1(3)-r2(3))**2 )
   call give_all_mos_at_r(r1,mos_array_r1)
   call give_all_mos_at_r(r2,mos_array_r2)
   if(dabs(r12-r12_test).gt.1.d-10)then
    print*,'error' 
    print*,r12,r12_test
   endif
   call local_r12_operator_on_hf(r1,r2,integral)
   if(integral.le.0.d0)then
    print*,integral,r12,mos_array_r1(1)*mos_array_r2(1) 
    pause
   endif
   write(i_unit_output_array(i),'(100(F16.10,X))')theta,r12,1.d0/r12,integral,mos_array_r1(1)*mos_array_r2(1)
  enddo
  theta += dtheta
 enddo
end




subroutine test_erf_coulomb_oprerator
 implicit none
 double precision :: r1(3),r2(3)
 integer :: i,j,nx
 double precision :: dx, xmax
 double precision :: rinit(3)
 double precision :: coulomb
 double precision :: r12
 character*(128) :: output
 character*(128) :: filename
 integer :: i_unit_output,getUnitAndOpen
 provide ezfio_filename 
 output=trim(ezfio_filename)//'.r12_erf'
 output=trim(output)
 print*,'output = ',trim(output)
 i_unit_output = getUnitAndOpen(output,'w')
 xmax = 5.d0
 nx = 100
 dx = xmax/dble(nx)
 rinit = 0.d0
 r1 = rinit

 r2 = r1
 do i = 1, nx 
  r2(1) += dx
  r2(2) += dx
  r12 = dsqrt((r1(1)-r2(1))**2 + (r1(2)-r2(2))**2 +(r1(3)-r2(3))**2 )
! call erf_coulomb_operator_in_real_space(r1,r2,coulomb) 
  call local_r12_operator_on_hf(r1,r2,coulomb)
  write(i_unit_output,*)r12,1.d0/r12 * (1.d0 - erfc(mu_erf * r12)) ,coulomb
 enddo
end

subroutine test_coulomb_oprerator
 implicit none
 double precision :: r1(3),r2(3)
 integer :: i,j,nx
 double precision :: dx, xmax
 double precision :: rinit(3)
 double precision :: coulomb
 double precision :: r12
 character*(128) :: output
 character*(128) :: filename
 integer :: i_unit_output,getUnitAndOpen
 provide ezfio_filename 
 output=trim(ezfio_filename)//'.r12'
 output=trim(output)
 print*,'output = ',trim(output)
 i_unit_output = getUnitAndOpen(output,'w')
 xmax = 2.d0
 nx = 100
 dx = xmax/dble(nx)
 rinit = 0.d0
 rinit(1) = nucl_coord(1,1)
 rinit(2) = nucl_coord(1,2)
 rinit(3) = nucl_coord(1,3)
 r1 = rinit

 r2 = r1
 do i = 1, nx 
  r2(1) += dx
  r2(2) += dx
  r12 = dsqrt((r1(1)-r2(1))**2 + (r1(2)-r2(2))**2 +(r1(3)-r2(3))**2 )
  call coulomb_operator_in_real_space(r1,r2,coulomb) 
  write(i_unit_output,*)r12,1.d0/r12,coulomb
 enddo
end

subroutine test_expectation_value_and_coulomb(nx,xmax)
 implicit none
 double precision, intent(in) :: xmax
 integer, intent(in) :: nx
 double precision :: r1(3),r2(3)
 integer :: i,j
 double precision :: dx,  dr(3)
 double precision :: rinit(3)
 double precision :: coulomb,coulomb_bis
 double precision :: r12
 character*(128) :: output
 character*(128) :: filename
 integer :: i_unit_output,getUnitAndOpen
 provide ezfio_filename 
 output=trim(ezfio_filename)//'.r12_psi'
 output=trim(output)
 print*,'output = ',trim(output)
 i_unit_output = getUnitAndOpen(output,'w')
 dx = xmax/dble(nx)
 rinit = 0.d0
 rinit(1) = nucl_coord(1,1)
 rinit(2) = nucl_coord(1,2)
 rinit(3) = nucl_coord(1,3)
 r1 = rinit
!r1(1) = 0.5d0

 dr = 0.d0
 dr(1) = dx
 dr(2) = dx

 r2 = r1
!do i = 1, nx 
! r2 += dr
  r12 = dsqrt((r1(1)-r2(1))**2 + (r1(2)-r2(2))**2 +(r1(3)-r2(3))**2 )
  call expectation_value_in_real_space_old(r1,r2,coulomb_bis,two_dm) 
  double precision :: two_dm,aos_array(ao_num),rho_a,rho_b,two_body_dm
! call dm_dft_alpha_beta_and_all_aos_at_r(r2,rho_a,rho_b,aos_array)
! call coulomb_operator_in_real_space(r1,r2,coulomb) 
  write(i_unit_output,'(100(F32.20,X))')r12,1.d0/r12,coulomb,two_body_dm,coulomb/two_body_dm
  write(i_unit_output,'(100(F32.20,X))')r12,1.d0/r12,coulomb_bis,two_dm,coulomb_bis/two_dm
!enddo
end


subroutine test_nuclear_coulomb_oprerator
 implicit none
 double precision :: r1(3),r2(3)
 integer :: i,j,nx
 double precision :: dx, xmax
 double precision :: rinit(3)
 double precision :: coulomb
 double precision :: r12
 character*(128) :: output
 character*(128) :: filename
 integer :: i_unit_output,getUnitAndOpen
 provide ezfio_filename 
 output=trim(ezfio_filename)//'.nuclear'
 output=trim(output)
 print*,'output = ',trim(output)
 i_unit_output = getUnitAndOpen(output,'w')
 xmax = 10.d0
 nx = 10000
 dx = xmax/dble(nx)
 rinit = 0.d0
 r1 = rinit

 r2 = r1
 do i = 1, nx 
  r2(1) += dx
  r12 = dsqrt((r1(1)-r2(1))**2 + (r1(2)-r2(2))**2 +(r1(3)-r2(3))**2 )
  call nuclear_coulomb_operator_in_real_space(r2,coulomb) 
  write(i_unit_output,*)r12,-2.* 1.d0/r12,coulomb
 enddo
end

subroutine test_grad
 implicit none
 integer :: i,j,k,l,m
 double precision :: r(3), rdx_plus(3),dr(3),rdx_minus(3)
 double precision :: dm_a(N_states),dm_b(N_states)
 double precision :: grad_dm_a(3,N_states),grad_dm_b(3,N_states)
 double precision :: grad_aos_array(3,ao_num)
 double precision :: grad_aos_array_bis(3,ao_num)

 double precision :: aos_array(ao_num),aos_array_plus(ao_num),aos_array_minus(ao_num)

 double precision :: xmax , dx_2
 double precision :: dx
 xmax = 5.d0
 integer :: nx
! variables to explore the r3 space
 nx = 1000
 dx_2 = dble(xmax/nx)
 dr(1) = dx_2
 dr(2) = dx_2
 dr(3) = dx_2

! dx for differentiation 
 do l = 1, 16
  dx = 1.d0/(10.d0**l)

  r = 0.d0
  double precision :: accu(3)
  accu = 0.d0
  do i = 1, nx
   call give_all_aos_and_grad_at_r(r,aos_array,grad_aos_array)
   do m = 1, 3
    rdx_plus = r
    rdx_plus(m) +=  dx
    rdx_minus = r
    rdx_minus(m) -=  dx
    call give_all_aos_at_r(rdx_plus,aos_array_plus)
    call give_all_aos_at_r(rdx_minus,aos_array_minus)
    do k = 1, ao_num
     grad_aos_array_bis(m,k) = (aos_array_plus(k) - aos_array_minus(k)) / (2.d0 * dx)
     accu(m) += dabs(grad_aos_array_bis(m,k) - grad_aos_array(m,k)) 
     if(isnan(dabs(grad_aos_array_bis(m,k) - grad_aos_array(m,k))))then
      print*,'AHAHAH !!' 
      print*,r
      print*,grad_aos_array_bis(m,k),grad_aos_array(m,k)
      pause 
     endif
    !if(dabs(grad_aos_array_bis(m,k) - grad_aos_array(m,k)).gt.1.d-8)then
    ! print*,'AHAHAH !!' 
    ! print*,'gradients '
    ! print*,dabs(grad_aos_array_bis(m,k) - grad_aos_array(m,k)),grad_aos_array_bis(m,k),grad_aos_array(m,k)
    ! print*,'positions'
    ! print*,r
    ! print*,rdx_plus
    ! print*,rdx_minus
    ! pause
    !endif
    enddo
   enddo
   r += dr
  enddo
  accu *= dr
 !print*,'accu = ',accu(1)
 !print*,'accu = ',accu(2)
 !print*,'accu = ',accu(3)
  print*,dx,(accu(1)+accu(2)+accu(3))/3.d0
 enddo
end

subroutine test_naive_grid
 implicit none
 integer :: i,j
 double precision :: r(3)
 double precision :: rmax,dx
 
 double precision, allocatable :: mos_array(:)
 rmax = 10.d0
 integer :: nx
 nx = 10000
 allocate(mos_array(mo_tot_num))
 dx = rmax/dble(nx)
 r = 0.d0
 r(3) = -rmax * 0.5d0
 do i = 1, nx
  call give_all_mos_at_r(r,mos_array)
  write(33,*)r(3),mos_array(1),mos_array(2)
  r(3) += dx
 enddo
!do j = 1, n_radial
! do i = 1, n_points_integration_angular
!  r = r_points(:,i,j) 
!  
! enddo
!enddo



end
