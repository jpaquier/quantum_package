
 BEGIN_PROVIDER [double precision, on_top_of_r,(n_points_integration_angular,n_points_radial_grid,nucl_num,N_states) ]
 implicit none
 BEGIN_DOC
 ! on_top computation 
 END_DOC
 integer :: j,k,l,istate
 double precision :: two_dm_in_r,on_top_in_r_sorted
 double precision, allocatable :: r(:)
 double precision :: cpu0,cpu1
 allocate(r(3))

 r = 0.d0
 istate = 1
 ! initialization for parallel loops
 If (ontop_approx .EQV. .TRUE.) Then
  on_top_of_r(1,1,1,1) = on_top_in_r_sorted(r,istate) 
 else      
  double precision :: two_dm,two_dm_laplacian,total_dm
  on_top_of_r(1,1,1,1) = two_dm_in_r(r,r,istate)
  call spherical_averaged_two_dm_at_second_order(r,0.d0,istate,two_dm,two_dm_laplacian,total_dm)
 endif

 print*,'providing the on_top_of_r ...'
 call wall_time(cpu0)
 If (ontop_approx .EQV. .TRUE.) Then
  do j = 1, nucl_num
   do k = 1, n_points_radial_grid  -1
    do l = 1, n_points_integration_angular
     r(1) = grid_points_per_atom(1,l,k,j)
     r(2) = grid_points_per_atom(2,l,k,j)
     r(3) = grid_points_per_atom(3,l,k,j)
     do istate = 1, N_states
!!!!!!!!!!!! CORRELATION PART
       on_top_of_r(l,k,j,istate) = on_top_in_r_sorted(r,istate) 
     enddo
    enddo
   enddo
  enddo
 else
  integer :: i_point,i
  do i_point = 1, n_points_final_grid
   k = index_final_points(1,i_point)
   i = index_final_points(2,i_point)
   j = index_final_points(3,i_point)
   on_top_of_r(k,i,j,istate) = on_top_of_r_vector(i_point,1)
  enddo  
 endif
 call wall_time(cpu1)
 print*,'Time to provide on_top_of_r = ',cpu1-cpu0
 deallocate(r)
 END_PROVIDER

 BEGIN_PROVIDER [double precision, on_top_of_r_vector,(n_points_final_grid,N_states) ]
&BEGIN_PROVIDER [double precision, mu_of_r_cusp_condition_vector,(n_points_final_grid,N_states) ]
 implicit none
 integer :: i_point,istate
 double precision :: two_dm_in_r_selected_points,dpi,r(3),two_dm,two_dm_laplacian,total_dm
 double precision :: two_dm_HF,two_dm_laplacian_HF,total_dm_HF
 double precision :: corr_hole_2,alpha,mu0
 istate = 1
 double precision :: wall_0,wall_1
 double precision :: alpha_bis,delta,beta,mu
 print*,'providing the on_top_of_r_vector'
 i_point = 1
 on_top_of_r_vector(i_point,istate) = two_dm_in_r_selected_points(i_point,istate)
 dpi = 1.5d0 * dsqrt(dacos(-1.d0))
 call wall_time(wall_0)
 !$OMP PARALLEL DO &
 !$OMP DEFAULT (NONE)  &
 !$OMP PRIVATE (i_point,r,two_dm,two_dm_laplacian,total_dm,two_dm_HF,two_dm_laplacian_HF,total_dm_HF,alpha,mu0, & 
 !$OMP          alpha_bis,delta,beta,mu) &
 !$OMP SHARED(on_top_of_r_vector,istate,n_points_final_grid,dpi,mu_of_r_cusp_condition_vector,final_grid_points)
 do i_point = 1, n_points_final_grid
  r(1) = final_grid_points(1,i_point)
  r(2) = final_grid_points(2,i_point)
  r(3) = final_grid_points(3,i_point)
  call spherical_averaged_two_dm_at_second_order(r,0.d0,istate,two_dm,two_dm_laplacian,total_dm)
  call spherical_averaged_two_dm_HF_at_second_order(r,0.d0,istate,two_dm_HF,two_dm_laplacian_HF,total_dm_HF)
 !corr_hole_2 = (two_dm_laplacian - two_dm_laplacian_HF * two_dm/total_dm_HF)/two_dm_HF
  two_dm = max(two_dm,1.d-15)
  two_dm_HF = max(two_dm_HF,1.d-15)
  on_top_of_r_vector(i_point,istate) = two_dm
!! approximated polynom 
  mu0 =  dpi * (two_dm_laplacian / two_dm - two_dm_laplacian_HF / two_dm_HF)
  mu0 = max(mu0,1.d-15)
  alpha = 1.d0 + 2.d0/(dsqrt(dacos(-1.d0)) * mu0)
  mu = mu0*alpha
!! new version with exact polynom 
  alpha_bis = (two_dm_laplacian / two_dm - two_dm_laplacian_HF / two_dm_HF)
  alpha_bis = max(alpha_bis,1.d-15)
  beta = 2.d0/(3.d0*dacos(-1.d0))
  delta = 2.d0/dsqrt(dacos(-1.d0))
  mu = 1.d0/(2.d0*beta)*(alpha_bis + dsqrt(alpha_bis*alpha_bis + 4.d0 * alpha_bis * beta * delta))

  mu_of_r_cusp_condition_vector(i_point,istate) = mu
 enddo
 !$OMP END PARALLEL DO
 call wall_time(wall_1)
 print*,'provided the on_top_of_r_vector'
 print*,'Time to provide :',wall_1 - wall_0
 END_PROVIDER 
 

 BEGIN_PROVIDER [double precision, integral_two_body]
 implicit none
 integer :: i_point,istate
 double precision :: cpu0,cpu1
 istate = 1
 integral_two_body = 0.d0
 call wall_time(cpu0)
 do i_point = 1, n_points_final_grid
  integral_two_body += on_top_of_r_vector(i_point,istate)
 enddo
 call wall_time(cpu1)
 print*,'Time to provide on_top_of_r = ',cpu1-cpu0
 END_PROVIDER 
 


double precision function on_top_in_r_sorted(r,istate)
 implicit none
 BEGIN_DOC
 !Compute the on top locally  
 END_DOC
 integer, intent(in) :: istate
 double precision, intent(in) :: r(3)
 double precision, allocatable :: mos_array_r(:)
 integer :: i,j,k,l,m,t,icouple,n,jcouple
 double precision :: tmp,psi_temp_l,psi_temp_r
 double precision :: u_dot_v
 allocate(mos_array_r(mo_tot_num))
 call give_all_mos_at_r(r,mos_array_r)
 on_top_in_r_sorted = 0.d0
 tmp = 0.d0
 do m = 1, n_coupleij_bis(istate)
  icouple = id_coupleij(m,istate)
  i = couple_to_array_reverse(icouple,1)
  j = couple_to_array_reverse(icouple,2) 
  tmp=0.d0
  do n = 1, n_kl(m,istate)
   jcouple = n_kl_id(n,m,istate)
   k = couple_to_array_reverse(jcouple,1)
   l = couple_to_array_reverse(jcouple,2)
   tmp += gamma_ijkl(n,m,istate) * mos_array_r(k)*mos_array_r(l)
  enddo
  on_top_in_r_sorted += tmp *mos_array_r(i) * mos_array_r(j)
 enddo
 deallocate(mos_array_r)
end

double precision function two_dm_in_r(r1,r2,istate)
 implicit none
 integer, intent(in) :: istate
 double precision, intent(in) :: r1(3),r2(3)
 double precision, allocatable :: mos_array_r1(:), mos_array_r2(:)
 integer :: i,j,k,l
 allocate(mos_array_r2(mo_tot_num), mos_array_r1(mo_tot_num))
 call give_all_mos_at_r(r1,mos_array_r1) 
 call give_all_mos_at_r(r2,mos_array_r2) 
 two_dm_in_r = 0.d0
 do i = 1, mo_tot_num
  do j = 1, mo_tot_num
   do k = 1, mo_tot_num
    do l = 1, mo_tot_num
     two_dm_in_r += two_bod_alpha_beta_mo_transposed(l,k,j,i,istate) * mos_array_r1(i) * mos_array_r1(l) * mos_array_r2(k) * mos_array_r2(j)
    enddo
   enddo
  enddo
 enddo
 two_dm_in_r = max(two_dm_in_r,1.d-15)
end


double precision function two_dm_in_r_selected_points(i_point,istate)
 implicit none
 integer, intent(in) :: istate,i_point
 integer :: i,j,k,l
 two_dm_in_r_selected_points = 0.d0
 do i = 1, mo_tot_num
  do j = 1, mo_tot_num
   do k = 1, mo_tot_num
    do l = 1, mo_tot_num
     two_dm_in_r_selected_points += two_bod_alpha_beta_mo_transposed(l,k,j,i,istate) * mos_in_r_array(i,i_point) * mos_in_r_array(l,i_point) * mos_in_r_array(k,i_point) * mos_in_r_array(j,i_point)
    enddo
   enddo
  enddo
 enddo
 two_dm_in_r_selected_points = max(two_dm_in_r_selected_points,1.d-15)
end

 BEGIN_PROVIDER [double precision, on_top_of_r_vector_simple,(n_points_final_grid,N_states) ]
 implicit none
 integer :: ipoint
 double precision :: on_top_of_r_from_provider
 double precision :: wall0,wall1
 on_top_of_r_vector_simple = 0.d0
 provide two_bod_alpha_beta_mo_physician mos_in_r_array 
 
 call wall_time(wall0)
 !$OMP PARALLEL        &
 !$OMP DEFAULT (NONE)  &
 !$OMP PRIVATE (ipoint) & 
 !$OMP SHARED  (n_points_final_grid, on_top_of_r_vector_simple )
 !$OMP DO              
 do ipoint = 1, n_points_final_grid
  on_top_of_r_vector_simple(ipoint,1) = on_top_of_r_from_provider(ipoint)
 enddo
 !$OMP END DO
 !$OMP END PARALLEL
 call wall_time(wall1)
 print*,'Time to provide on_top_of_r_vector_simple : ',wall1-wall0

 END_PROVIDER 

 double precision function on_top_of_r_from_provider(ipoint)
 implicit none
 integer, intent(in) :: ipoint
 integer :: i,j,k,l 
 on_top_of_r_from_provider = 0.d0
 do l = 1, mo_tot_num ! 2 
  do k = 1, mo_tot_num ! 1 
    do j = 1, mo_tot_num
     do i = 1, mo_tot_num
     on_top_of_r_from_provider += two_bod_alpha_beta_mo_physician(i,j,k,l,1) * mos_in_r_array(j,ipoint) * mos_in_r_array(i,ipoint) * mos_in_r_array(l,ipoint) * mos_in_r_array(k,ipoint)
    enddo
   enddo
  enddo
 enddo
 end
