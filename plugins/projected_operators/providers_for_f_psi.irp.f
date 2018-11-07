BEGIN_PROVIDER [double precision, V_kl_contracted, (mo_tot_num,mo_tot_num,n_points_final_grid)]
 implicit none
 BEGIN_DOC
! V_kl_contracted(k,l,ipoint) = \sum_{ij} V_{ij}^{kl} phi_i(r_ipoint) phi_j(r_ipoint)
 END_DOC
 integer :: i,j,k,l
 double precision :: ipoint
 double precision, allocatable :: integrals_array(:,:), mos_array_r(:),r(:)
 allocate(integrals_array(mo_tot_num,mo_tot_num),mos_array_r(mo_tot_num),r(3))
 do ipoint = 1, n_points_final_grid
  r(1) = final_grid_points(1,ipoint)
  r(2) = final_grid_points(2,ipoint)
  r(3) = final_grid_points(3,ipoint)
  call give_all_mos_at_r(r,mos_array_r) 
  do l = 1, mo_tot_num ! 2 
   do k = 1, mo_tot_num ! 1 
    V_kl_contracted(k,l,ipoint) = 0.d0
    call get_mo_bielec_integrals_ij(k,l,mo_tot_num,integrals_array,mo_integrals_map) 
    do j = 1, mo_tot_num
     do i = 1, mo_tot_num
                     !1 2                            1 2
      V_kl_contracted(k,l,ipoint) += integrals_array(i,j) * mos_array_r(i) * mos_array_r(j)
     enddo
    enddo
   enddo
  enddo
 enddo

END_PROVIDER 

BEGIN_PROVIDER [double precision, rho2_mn_contracted, (mo_tot_num,mo_tot_num,n_points_final_grid)]
 implicit none
 BEGIN_DOC
! rho2_kl_contracted(k,l,ipoint) = \sum_{ij} rho2_{ij}^{kl} phi_i(r_ipoint) phi_j(r_ipoint)
 END_DOC
 integer :: i,j,k,l
 double precision :: ipoint
 double precision, allocatable :: integrals_array(:,:), mos_array_r(:),r(:)
  allocate(integrals_array(mo_tot_num,mo_tot_num),mos_array_r(mo_tot_num),r(3))
 do ipoint = 1, n_points_final_grid
  r(1) = final_grid_points(1,ipoint)
  r(2) = final_grid_points(2,ipoint)
  r(3) = final_grid_points(3,ipoint)
  call give_all_mos_at_r(r,mos_array_r) 
  do l = 1, mo_tot_num ! 2 
   do k = 1, mo_tot_num ! 1 
    rho2_mn_contracted(k,l,ipoint) = 0.d0
    do j = 1, mo_tot_num
     do i = 1, mo_tot_num
                        !1 2                                            1 2 1 2 
      rho2_mn_contracted(k,l,ipoint) += two_bod_alpha_beta_mo_physician(i,j,k,l,1) * mos_array_r(i) * mos_array_r(j)
     enddo
    enddo
   enddo
  enddo
 enddo

END_PROVIDER 


BEGIN_PROVIDER [double precision, f_psi_B, (n_points_final_grid)]
 implicit none
 double precision :: ipoint
 integer :: k,l 
 do ipoint = 1, n_points_final_grid
  f_psi_B(ipoint) = 0.d0
  do l = 1, mo_tot_num ! 2 
   do k = 1, mo_tot_num ! 1
    f_psi_B(ipoint) += V_kl_contracted(k,l,ipoint) * rho2_mn_contracted(k,l,ipoint)
   enddo
  enddo
 enddo
END_PROVIDER 


BEGIN_PROVIDER [double precision, f_psi_B_old, (n_points_final_grid)]
 implicit none
 double precision :: ipoint
 double precision :: r(3),coulomb,two_body_dm
 integer :: k,l 
 do ipoint = 1, n_points_final_grid
  f_psi_B_old(ipoint) = 0.d0
  r(1) = final_grid_points(1,ipoint)
  r(2) = final_grid_points(2,ipoint)
  r(3) = final_grid_points(3,ipoint)
  call expectation_value_in_real_space_no_divide(r,r,coulomb,two_body_dm)
  f_psi_B_old(ipoint) = coulomb
 enddo
END_PROVIDER 

