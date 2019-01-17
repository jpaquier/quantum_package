 BEGIN_PROVIDER[complex*16, dirac_aos_vc_LDA_w, (n_points_final_grid,2*dirac_ao_num,N_states)]
 &BEGIN_PROVIDER[complex*16, dirac_aos_vx_LDA_w, (n_points_final_grid,2*dirac_ao_num,N_states)]
 implicit none
 BEGIN_DOC
 ! dirac_aos_vxc_LDA_w(j,i) = dirac_ao_i(r_j) * (v^x_alpha(r_j) + v^c_alpha(r_j)) * W(r_j)
 END_DOC
 integer :: istate,i,j
 double precision :: r(3)
 double precision :: mu,weight
 double precision :: e_c,v_c,e_x,v_x
 double precision, allocatable :: rho(:)
 allocate(rho(N_states))
 do istate = 1, N_states
  do i = 1, n_points_final_grid
   r(1) = final_grid_points(1,i)
   r(2) = final_grid_points(2,i)
   r(3) = final_grid_points(3,i)
   weight=final_weight_functions_at_final_grid_points(i)
   rho(istate) = dirac_one_body_dm_at_r(i,istate)
   call dirac_ec_LDA_sr(mu_erf,rho(istate),e_c,v_c)
   call dirac_ex_LDA_sr(mu_erf,rho(istate),e_x,v_x)
   do j =1, 2*dirac_ao_num
    dirac_aos_vc_LDA_w(i,j,istate) = (1.d0,0.d0)* v_c * dirac_aos_in_r_array(j,i)*weight
    dirac_aos_vx_LDA_w(i,j,istate) = (1.d0,0.d0)* v_x * dirac_aos_in_r_array(j,i)*weight
   enddo
  enddo
 enddo
 END_PROVIDER 

 BEGIN_PROVIDER[double precision, dirac_energy_x_LDA, (N_states) ]
 &BEGIN_PROVIDER[double precision, dirac_energy_c_LDA, (N_states) ]
 implicit none
 BEGIN_DOC
 ! exchange/correlation energy with the relativistic short range LDA functional
 END_DOC
 integer :: istate,i,j
 double precision :: r(3)
 double precision :: mu,weight
 double precision :: e_c,v_c,e_x,v_x
 double precision, allocatable :: rho(:)
 allocate(rho(N_states))
 dirac_energy_x_LDA = 0.d0
 dirac_energy_c_LDA = 0.d0
 do istate = 1, N_states
  do i = 1000, 10000
! do i = 1, n_points_final_grid
   r(1) = final_grid_points(1,i)
   r(2) = final_grid_points(2,i)
   r(3) = final_grid_points(3,i)
   weight=final_weight_functions_at_final_grid_points(i)
   rho(istate) = dirac_one_body_dm_at_r(i,istate)
   call dirac_ec_LDA_sr(mu_erf,rho(istate),e_c,v_c)
   call dirac_ex_LDA_sr(mu_erf,rho(istate),e_x,v_x)

   print*,r(1),r(2),r(3), rho(istate),e_x,v_x

   dirac_energy_x_LDA(istate) += weight * e_x
   dirac_energy_c_LDA(istate) += weight * e_c
  enddo
 enddo
 END_PROVIDER 


 BEGIN_PROVIDER [complex*16, dirac_potential_x_ao_LDA,(2*dirac_ao_num,2*dirac_ao_num,N_states)]
 &BEGIN_PROVIDER [complex*16, dirac_potential_c_ao_LDA,(2*dirac_ao_num,2*dirac_ao_num,N_states)]
 implicit none
 BEGIN_DOC 
 ! short range exchange/correlation alpha/beta potentials with LDA functional on the AO basis
 END_DOC
 integer :: istate
 double precision :: wall_1,wall_2
 call wall_time(wall_1)
 do istate = 1, N_states 
  call zgemm('N','C',2*dirac_ao_num,2*dirac_ao_num,n_points_final_grid,(1.d0,0.d0),dirac_aos_in_r_array,2*dirac_ao_num,dirac_aos_vc_LDA_w(1,1,istate),n_points_final_grid,(0.d0,0.d0),dirac_potential_c_ao_LDA(1,1,istate),2*dirac_ao_num)
  call zgemm('N','C',2*dirac_ao_num,2*dirac_ao_num,n_points_final_grid,(1.d0,0.d0),dirac_aos_in_r_array,2*dirac_ao_num,dirac_aos_vx_LDA_w(1,1,istate),n_points_final_grid,(0.d0,0.d0),dirac_potential_x_ao_LDA(1,1,istate),2*dirac_ao_num)
 enddo
 call wall_time(wall_2)
 print*,'time to provide dirac_potential_x/c_ao_LDA = ',wall_2 - wall_1
 END_PROVIDER 

