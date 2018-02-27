program pouet
 read_wf = .True.
 touch read_wf
 integer :: nx
 double precision :: rmax
 nx = 100
 rmax = 5.d0
 call test_expectation_value_and_coulomb(nx,rmax)
!call routine3
!call test_rho2
!call test_one_dm_ao
!call test_coulomb_oprerator
!call test_expectation_value
!call test_erf_coulomb_oprerator
!call test_nuclear_coulomb_oprerator
!call test_integratio_mo
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
 double precision ::r(3),rho2
 double precision :: test

 r(1) = 0.d0
 r(2) = 0.d0
 r(3) = 0.d0
 call  on_top_pair_density_in_real_space(r,rho2)
 print*,'rho2(0) = ',rho2
!stop
 test = 0.d0
  do j = 1, nucl_num
   do k = 1, n_points_radial_grid  -1
    print*,k
    do l = 1, n_points_integration_angular 
     
     r(1) = grid_points_per_atom(1,l,k,j)
     r(2) = grid_points_per_atom(2,l,k,j)
     r(3) = grid_points_per_atom(3,l,k,j)
     call  on_top_pair_density_in_real_space(r,rho2)
!    call  on_top_pair_density_in_real_space_from_ao(r,rho2)
     test += rho2 * final_weight_functions_at_grid_points(l,k,j) 
     enddo
    enddo
   enddo
 print*,'test = ',test

end


subroutine test_one_dm_ao
 implicit none
 double precision :: r(3)
 double precision, allocatable :: mos_array(:), aos_array(:), aos_array_bis(:), mos_array_bis(:)
 allocate(mos_array(mo_tot_num),aos_array(ao_num),aos_array_bis(ao_num),mos_array_bis(mo_tot_num))
 integer :: m,n,p
 integer :: j,k,l 
  do j = 1, nucl_num
   do k = 1, n_points_radial_grid  -1
    print*,k
    do l = 1, n_points_integration_angular 
     r(1) = grid_points_per_atom(1,l,k,j)
     r(2) = grid_points_per_atom(2,l,k,j)
     r(3) = grid_points_per_atom(3,l,k,j)
     call give_all_mos_at_r(r,mos_array) 
     call give_all_aos_at_r(r,aos_array)
     aos_array_bis = 0.d0
     do m = 1, ao_num
      print*,ao_overlap(m,m)
      do n = 1, mo_tot_num 
       do p = 1, ao_num
        aos_array_bis(m) += mos_array(n) * mo_coef_transp(n,p) * ao_overlap(p,m)
       !if(dabs(mos_array(n) * mo_coef_transp(n,p) * ao_overlap(p,m)).gt.1.d-10)then
       ! print*,n,p,m
       ! print*,mos_array(n) , mo_coef_transp(n,m) , ao_overlap(p,m)
       !endif
       enddo
      enddo
      if(dabs(aos_array_bis(m) - aos_array(m)).gt.1.d-10)then
       print*,'PB AO !!' 
       print*,r
       print*,m
       print*,ao_overlap(m,m)
       print*,dabs(aos_array_bis(m) - aos_array(m)),aos_array(m),aos_array_bis(m)
!      stop
      endif
!     print*,dabs(aos_array_bis(m) - aos_array(m))!,aos_array(m),aos_array_bis(m)
     enddo


    !do m = 1, ao_num
    ! do n = 1, mo_tot_num 
    !   aos_array_bis(m) += mos_array(n) * mo_to_ao_matrix(n,m)
    ! enddo
    ! if(dabs(aos_array_bis(m) - aos_array(m)).gt.1.d-10)then
    !  print*,'PB AO !!' 
    !  print*,r
    !  print*,dabs(aos_array_bis(m) - aos_array(m)),aos_array(m),aos_array_bis(m)
    !  stop
    ! endif
    ! print*,dabs(aos_array_bis(m) - aos_array(m))!,aos_array(m),aos_array_bis(m)
    !enddo
    stop
    enddo
   enddo
  enddo
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
 xmax = 10.d0
 nx = 100
 dx = xmax/dble(nx)
 rinit = 0.d0
 r1 = rinit

 r2 = r1
 do i = 1, nx 
  r2(1) += dx
  r2(2) += dx
  r12 = dsqrt((r1(1)-r2(1))**2 + (r1(2)-r2(2))**2 +(r1(3)-r2(3))**2 )
  call erf_coulomb_operator_in_real_space(r1,r2,coulomb) 
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

 dr = 0.d0
 dr(3) = dx

 r2 = r1
 do i = 1, nx 
  r2 += dr
  r12 = dsqrt((r1(1)-r2(1))**2 + (r1(2)-r2(2))**2 +(r1(3)-r2(3))**2 )
  call expectation_value_in_real_space(r1,r2,coulomb_bis) 
  call coulomb_operator_in_real_space(r1,r2,coulomb) 
  write(i_unit_output,'(F10.7,X,F16.8,X,F16.8,X,F16.8)')r12,1.d0/r12,coulomb_bis,coulomb
 enddo
end


subroutine test_expectation_value
 implicit none
 double precision :: r1(3),r2(3)
 integer :: i,j,nx
 double precision :: dx, xmax
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
  call expectation_value_in_real_space(r1,r2,coulomb) 
  write(i_unit_output,*)r12,1.d0/r12,coulomb
 enddo
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
