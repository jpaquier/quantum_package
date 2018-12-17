 subroutine dirac_dm_dft_at_r(r,dm)
 implicit none
 BEGIN_DOC
 ! input: r(1) ==> r(1) = x, r(2) = y, r(3) = z
 ! output : dm = density evaluated at r(3)
 END_DOC
 double precision, intent(in) :: r(3)
 double precision, intent(out) :: dm(N_states)
 integer :: istate
 double precision  :: aos_array(ao_num),aos_array_bis(ao_num),u_dot_v
 call give_all_dirac_aos_at_r(r,dirac_aos_array)
 do istate = 1, N_states
 !aos_array_bis = aos_array
 !! alpha density
 !call dgemv('N',ao_num,ao_num,1.d0,one_body_dm_alpha_ao_for_dft(1,1,istate),ao_num,aos_array,1,0.d0,aos_array_bis,1)
 !dm_a(istate) = u_dot_v(aos_array,aos_array_bis,ao_num)
 !! beta density
 !aos_array_bis = aos_array
 !call dgemv('N',ao_num,ao_num,1.d0,one_body_dm_beta_ao_for_dft(1,1,istate),ao_num,aos_array,1,0.d0,aos_array_bis,1)
 !dm_b(istate) = u_dot_v(aos_array,aos_array_bis,ao_num)
 enddo
 end


 BEGIN_PROVIDER [double precision, dirac_one_body_dm_at_r, (n_points_final_grid,N_states) ]
 implicit none
 BEGIN_DOC
 ! dirac_one_body_dm_at_r(i,istate) = n(r_i,istate)
 ! where r_i is the ith point of the grid and istate is the state number
 END_DOC
 integer :: i,istate
 double precision :: r(3)
 double precision, allocatable :: dm(:)
 allocate(dm(N_states))
 do istate = 1, N_states
  do i = 1, n_points_final_grid
   r(1) = final_grid_points(1,i)
   r(2) = final_grid_points(2,i)
   r(3) = final_grid_points(3,i)
   call dirac_dm_dft_at_r(r,dm)
   dirac_one_body_dm_at_r(i,istate) = dm(istate)
  enddo
 enddo

END_PROVIDER 

