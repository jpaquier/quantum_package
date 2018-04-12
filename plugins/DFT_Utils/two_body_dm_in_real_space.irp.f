double precision function two_dm_in_r(r1,r2)
 implicit none
 double precision, intent(in) :: r1(3),r2(3)
 double precision :: mos_array_r1(mo_tot_num)
 double precision :: mos_array_r2(mo_tot_num)
 integer :: i,j,k,l
 call give_all_mos_at_r(r1,mos_array_r1) 
 call give_all_mos_at_r(r2,mos_array_r2) 
 two_dm_in_r = 0.d0
 do i = 1, mo_tot_num
  do j = 1, mo_tot_num
   do k = 1, mo_tot_num
    do l = 1, mo_tot_num
     two_dm_in_r += two_bod_alpha_beta_mo_transposed(l,k,j,i,1) * mos_array_r1(i) * mos_array_r1(l)  * mos_array_r2(k) * mos_array_r2(j)
    enddo
   enddo
  enddo
 enddo
end

double precision function on_top_dm_integral_with_mu_correction(mu)
 implicit none
 double precision, intent(in) :: mu
 double precision :: two_dm_in_r, pi, r(3)
 double precision :: weight
 integer :: j,k,l
 pi = 4d0 * datan(1d0)
 on_top_dm_integral_with_mu_correction = 0d0

 do j = 1, nucl_num
  do k = 1, n_points_radial_grid  -1
   do l = 1, n_points_integration_angular 
    r(:) = grid_points_per_atom(:,l,k,j)
    weight = final_weight_functions_at_grid_points(l,k,j) 
    on_top_dm_integral_with_mu_correction += two_dm_in_r(r,r)* weight**2d0
   enddo
  enddo
 enddo
 on_top_dm_integral_with_mu_correction = on_top_dm_integral_with_mu_correction / ( 1d0 + 2d0/(dsqrt(pi)*mu) )

end


double precision function Ec_md_mu_inf_corected(mu)
 implicit none
 double precision, intent(in) :: mu
 double precision :: on_top_dm_integral_with_mu_correction 
 double precision :: pi
 pi = 4d0 * datan(1d0)
 Ec_md_mu_inf_corected = ((-2d0+sqrt(2d0))*sqrt(2d0*pi)/(3d0*(mu**3)))*on_top_dm_integral_with_mu_correction(mu)
end
