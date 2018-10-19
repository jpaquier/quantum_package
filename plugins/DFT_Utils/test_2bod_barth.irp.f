program pouet
 implicit none
 integer :: nktotto
 read_wf = .true.
 touch read_wf
 !provide k_sorted_order_couple 
 call routine2
! call number_of_k(nktotto)
end


subroutine routine2
double precision :: test_ter2,test,test_ter,wall_1,wall_2,wall_3,wall_4

 call wall_time(wall_1)
 call test_rho2_bourrin(test)
 call wall_time(wall_2)

 !call test_rho2_selec_k_num_inte_sorted(test_ter)

 call wall_time(wall_3)
 call test_rho2_selec_k_num_inte_sortedi2(test_ter2)
 call wall_time(wall_4)
 print*,'*******************************'
 print*,'Bourin Numerical integrals                  = ',test 
 print*,'Selected k sorted new numerical integrals 2 = ',test_ter2
 print*,'Err sorted/bourrin  = ',dabs(test_ter2 - test)/dabs(test)
 print*,'  '
 print*,'  '
 print*,'wall time tot bourrin  = ',wall_2 - wall_1
 print*,'wall time tot Sorted   = ',wall_4 - wall_3
end

subroutine routine 
 implicit none
 double precision :: test_bbis,test_bart,test_ter,test_ter2,test,test_bart2,test_inte,rho2_ana(1)
 call on_top_pair_density_thresh_ec(rho2_ana)
 call test_rho2_selec_k_num_inte(test_bbis)
 call test_rho2_integral(test_inte)
 double precision :: wall_1,wall_2,wall_3,wall_4
 call wall_time(wall_1)
 call test_rho2_bourrin(test)
 call wall_time(wall_2) 

 call test_rho2_selec_k_num_inte_sorted(test_ter)
 
 call wall_time(wall_3)
 call test_rho2_selec_k_num_inte_sortedi2(test_ter2)
 call wall_time(wall_4)
 print*,' '
 print*,' '
 print*,'***********Error*******'
 print*,'E_cor_tot                                   = ',E_cor_tot
 print*,'Bourin Numerical integrals                  = ',test
 print*,'Bourin analatical integrals                 = ',test_inte
 print*,'Selected k analatical integrals             = ',rho2_ana
 print*,'Selected k numerical integrals              = ',test_bbis
 !print*,'Selected k sorted new numerical integrals   = ',test_ter 
 print*,'Selected k sorted new numerical integrals 2 = ',test_ter2
 print*,'*******************'
 print*,'*******************'
 print*,'*******************'
 print*,'*******************'
!print*,'Total          = ',E_cor_tot 
!print*,'couple         = ',E_cor_couple
!print*,'rho2_ana       = ',rho2_ana
 print*,'**'
 print*,'**'
 print*,'**'
 print*,'Err couple/tot = ',dabs(E_cor_tot(1) - E_cor_couple(1))/dabs(E_cor_tot(1))
 print*,'thr couple     = ',thr_couple_2dm
 print*,'**                             '
 print*,'Err couple/eig = ',dabs(E_cor_couple(1) - rho2_ana(1))/dabs(E_cor_couple(1))
 print*,'thr eigen      = ',thr_eig_2dm
 print*,'**'
 print*,'Err tot/final  = ',dabs(E_cor_tot(1) - rho2_ana(1))/dabs(E_cor_tot(1))
 print*,'**'
 print*,'**'
 print*,'********Sorted matrix elemets'
 print*,'**'
 print*,'Err sorted/bourrin  = ',dabs(test_ter2 - test)/dabs(test)
 print*,'  '
 print*,'  '
 print*,'wall time tot bourrin  = ',wall_2 - wall_1
 print*,'wall time tot Sorted   = ',wall_4 - wall_3




end




subroutine test_rho2_bourrin(test)
 implicit none
 double precision, intent(out) :: test
 integer :: j,k,l,istate
 double precision :: r(3),rho2
 double precision :: two_dm_in_r
 double precision :: wall_1, wall_2

do istate = 1, N_states
 test = 0.d0
 r(1) = 0.d0
 r(2) = 0.d0
 r(3) = 0.d0
 call wall_time(wall_1)
  do j = 1, nucl_num
   do k = 1, n_points_radial_grid  -1
    do l = 1, n_points_integration_angular
     r(1) = grid_points_per_atom(1,l,k,j)
     r(2) = grid_points_per_atom(2,l,k,j)
     r(3) = grid_points_per_atom(3,l,k,j)
     rho2 = two_dm_in_r(r,r,istate)
     test += rho2 * final_weight_functions_at_grid_points(l,k,j)
    enddo
   enddo
  enddo
 enddo
 call wall_time(wall_2)
 print*,'wall time bourrin num inte = ',wall_2 - wall_1
end

subroutine test_rho2_selec_k_num_inte(test_bbis)
 implicit none
 double precision, intent(out) :: test_bbis
 integer :: j,k,l,istate
 double precision :: r(3),rho2
 double precision :: two_dm_in_r_k_selected
 double precision :: wall_1, wall_2

 do istate = 1, N_states
 test_bbis = 0.d0
 r(1) = 0.d0
 r(2) = 0.d0
 r(3) = 0.d0
 call wall_time(wall_1)
  do j = 1, nucl_num
   do k = 1, n_points_radial_grid  -1
    do l = 1, n_points_integration_angular
     r(1) = grid_points_per_atom(1,l,k,j)
     r(2) = grid_points_per_atom(2,l,k,j)
     r(3) = grid_points_per_atom(3,l,k,j)
     rho2 = two_dm_in_r_k_selected(r,r,istate)
     test_bbis += rho2 * final_weight_functions_at_grid_points(l,k,j)
    enddo
   enddo
  enddo
 enddo
 call wall_time(wall_2)
 print*,'wall time selected k num integral = ',wall_2 - wall_1
end

subroutine test_rho2_selec_k_num_inte_sorted(test_ter)
 implicit none
 double precision, intent(out) :: test_ter
 integer :: j,k,l,istate
 double precision :: r(3),rho2 
 double precision :: two_dm_in_r_k_selected_sorted
 double precision :: wall_1, wall_2
 
 do istate = 1, N_states
 test_ter  = 0.d0
 r(1) = 0.d0
 r(2) = 0.d0
 r(3) = 0.d0
 call wall_time(wall_1)
  do j = 1, nucl_num
   do k = 1, n_points_radial_grid  -1
    do l = 1, n_points_integration_angular
     r(1) = grid_points_per_atom(1,l,k,j)
     r(2) = grid_points_per_atom(2,l,k,j)
     r(3) = grid_points_per_atom(3,l,k,j)
     rho2 = two_dm_in_r_k_selected_sorted(r,r,istate)
     test_ter += rho2 * final_weight_functions_at_grid_points(l,k,j)
    enddo
   enddo
  enddo
 enddo
 call wall_time(wall_2)
 print*,'wall time sorted k num integral = ',wall_2 - wall_1
end


subroutine test_rho2_selec_k_num_inte_sortedi2(test_ter2)
 implicit none
 double precision, intent(out) :: test_ter2
 integer :: j,k,l,istate
 double precision :: r(3),rho2
 double precision :: two_dm_in_r_k_sorted_couple 
 double precision :: wall_1, wall_2

 do istate = 1, N_states
 test_ter2  = 0.d0
 r(1) = 0.d0
 r(2) = 0.d0
 r(3) = 0.d0
 call wall_time(wall_1)
  do j = 1, nucl_num
   do k = 1, n_points_radial_grid  -1
    do l = 1, n_points_integration_angular
     r(1) = grid_points_per_atom(1,l,k,j)
     r(2) = grid_points_per_atom(2,l,k,j)
     r(3) = grid_points_per_atom(3,l,k,j)
     rho2 = two_dm_in_r_k_sorted_couple(r,r,istate)
     test_ter2 += rho2 * final_weight_functions_at_grid_points(l,k,j)
    enddo
   enddo
  enddo
 enddo
 call wall_time(wall_2)
 print*,'wall time sorted k num integral 2 = ',wall_2 - wall_1
end




subroutine test_rho2_integral(test_inte)
 implicit none
 double precision, intent(out) :: test_inte
 integer :: i,j,k,l,istate
 double precision :: wall_1, wall_2,get_mo_bielec_integral_ijkl_r3
 do istate = 1, N_states
 test_inte = 0.d0
 call wall_time(wall_1)
  do i = 1, mo_tot_num
   do j = 1, mo_tot_num
    do k = 1, mo_tot_num
     do l = 1, mo_tot_num
      test_inte += two_bod_alpha_beta_mo_transposed(l,k,j,i,istate) * get_mo_bielec_integral_ijkl_r3(l,k,j,i,mo_integrals_ijkl_r3_map)
     enddo
    enddo
   enddo
  enddo
 enddo
 call wall_time(wall_2)
 print*,'wall time bourrin analitycal integrals = ',wall_2 - wall_1
end

