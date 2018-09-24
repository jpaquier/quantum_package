
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
 If (ontop_approx .EQV. .TRUE.) Then
  on_top_of_r(1,1,1,1) = on_top_in_r_sorted(r,istate) 
 else      
  on_top_of_r(1,1,1,1) = two_dm_in_r(r,r,istate)
 endif

 print*,'providing the on_top_of_r ...'
 call cpu_time(cpu0)
  do j = 1, nucl_num
   do k = 1, n_points_radial_grid  -1
    do l = 1, n_points_integration_angular
     r(1) = grid_points_per_atom(1,l,k,j)
     r(2) = grid_points_per_atom(2,l,k,j)
     r(3) = grid_points_per_atom(3,l,k,j)
     do istate = 1, N_states
!!!!!!!!!!!! CORRELATION PART
      If (ontop_approx .EQV. .TRUE.) Then
       on_top_of_r(l,k,j,istate) = on_top_in_r_sorted(r,istate) 
      else      
       on_top_of_r(l,k,j,istate) = two_dm_in_r(r,r,istate)
      endif
     enddo
    enddo
   enddo
 enddo
 call cpu_time(cpu1)
 print*,'Time to provide on_top_of_r = ',cpu1-cpu0
 deallocate(r)
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

