
BEGIN_PROVIDER [double precision, V_kl_contracted_sequential, (mo_tot_num,mo_tot_num,n_points_final_grid)]
 implicit none
 BEGIN_DOC
! V_kl_contracted(k,l,ipoint) = \sum_{ij} V_{ij}^{kl} phi_i(r_ipoint) phi_j(r_ipoint)
 END_DOC
 integer :: i,j,k,l
 integer :: ipoint
 double precision, allocatable :: integrals_array(:,:), mos_array_r(:),r(:)
 allocate(integrals_array(mo_tot_num,mo_tot_num),mos_array_r(mo_tot_num),r(3))
 do ipoint = 1, n_points_final_grid
  r(1) = final_grid_points(1,ipoint)
  r(2) = final_grid_points(2,ipoint)
  r(3) = final_grid_points(3,ipoint)
  call give_all_mos_at_r(r,mos_array_r) 
  do l = 1, mo_tot_num ! 2 
   do k = 1, mo_tot_num ! 1 
    V_kl_contracted_sequential(k,l,ipoint) = 0.d0
    call get_mo_bielec_integrals_ij(k,l,mo_tot_num,integrals_array,mo_integrals_map) 
    do j = 1, mo_tot_num
     do i = 1, mo_tot_num
                     !1 2                            1 2
      V_kl_contracted_sequential(k,l,ipoint) += integrals_array(i,j) * mos_array_r(i) * mos_array_r(j)
     enddo
    enddo
   enddo
  enddo
 enddo
 deallocate(integrals_array)

END_PROVIDER 


BEGIN_PROVIDER [double precision, rho2_kl_contracted_sequential, (mo_tot_num,mo_tot_num,n_points_final_grid)]
 implicit none
 BEGIN_DOC
! rho2_kl_contracted(k,l,ipoint) = \sum_{ij} rho2_{ij}^{kl} phi_i(r_ipoint) phi_j(r_ipoint)
 END_DOC
 integer :: i,j,k,l
 integer :: ipoint
 double precision, allocatable :: mos_array_r(:),r(:)
 allocate(mos_array_r(mo_tot_num),r(3))
 do ipoint = 1, n_points_final_grid
  r(1) = final_grid_points(1,ipoint)
  r(2) = final_grid_points(2,ipoint)
  r(3) = final_grid_points(3,ipoint)
  call give_all_mos_at_r(r,mos_array_r) 
  do l = 1, mo_tot_num ! 2 
   do k = 1, mo_tot_num ! 1 
    rho2_kl_contracted_sequential(k,l,ipoint) = 0.d0
    do j = 1, mo_tot_num
     do i = 1, mo_tot_num
                        !1 2                                            1 2 1 2 
      rho2_kl_contracted_sequential(k,l,ipoint) += two_bod_alpha_beta_mo_physician(i,j,k,l,1) * mos_array_r(i) * mos_array_r(j)
     enddo
    enddo
   enddo
  enddo
 enddo
 deallocate(mos_array_r,r)

END_PROVIDER 

