program pouet
 implicit none
 read_wf = .True.
 call test_rho2
end

subroutine routine
 implicit none
 integer :: i,j,k,l
 double precision, allocatable :: mat(:,:),eigvec(:,:),eigval(:)
 double precision :: trace_1,trace_2,thr
 allocate(mat(mo_tot_num,mo_tot_num),eigvec(mo_tot_num,mo_tot_num),eigval(mo_tot_num))
 thr = 1.d-15
 do i = 1, mo_tot_num ! loop over the first electron 
  do j = 1, mo_tot_num ! loop over the second electron 
   ! get matrix 
   trace_1 = 0.d0
   do k = 1, mo_tot_num
    do l = 1, mo_tot_num
     mat(l,k) = two_bod_alpha_beta_mo_transposed(l,k,j,i,1)
    enddo
    trace_1 += mat(k,k)
   enddo
   ! diagonalize the matrix 
   call lapack_diagd(eigval,eigvec,mat,mo_tot_num,mo_tot_num)
   trace_2 = 0.d0
   do k = 1, mo_tot_num
   !print*,'eigval(k) = ',eigval(k)
    trace_2 += dabs(eigval(k))
   enddo
   if(dabs(trace_2).gt.1.d-10)then
    print*,''
    print*,'i,j',i,j
    print*,'trace = ',trace_2
    write(*,'(100(F10.5,X))')eigval(1),eigval(2),eigval(mo_tot_num-1),eigval(mo_tot_num)
   endif
  enddo
 enddo


 deallocate(mat,eigvec,eigval)



end


subroutine test_rho2
 implicit none
 integer :: j,k,l,istate,n_k_loc,m
 double precision :: r(3),rho2
 double precision, allocatable :: rho2_ap(:)
 double precision :: test,test_bart,two_dm_in_r
 allocate(rho2_ap(N_states))

do istate = 1, N_states
 r(1) = 0.d0
 r(2) = 0.d0
 r(3) = 0.d0
 call  on_top_pair_density_in_real_space(r,rho2)
 print*,'rho2(0) = ',rho2
 call  on_top_pair_density_approx(r,rho2_ap)
 print*,'rho2(0) barth= ',rho2_ap(istate)
 test = 0.d0
 test_bart= 0.d0 
  do j = 1, nucl_num
   do k = 1, n_points_radial_grid  -1
    print*,k
    do l = 1, n_points_integration_angular
     r(1) = grid_points_per_atom(1,l,k,j)
     r(2) = grid_points_per_atom(2,l,k,j)
     r(3) = grid_points_per_atom(3,l,k,j)
     rho2 = two_dm_in_r(r,r,istate)
     call on_top_pair_density_approx(r,rho2_ap)
     if(dabs(rho2-rho2_ap(istate)).gt.1.d-10)then
      print*,'' 
      print*,'r = '
      print*,r
      print*,'Problem!!, rho2(r)-rho2_bart(r) ='
      print*,rho2-rho2_ap(istate),rho2,rho2_ap(istate)
      print*,'' 
     endif 
!    call  on_top_pair_density_in_real_space_from_ao(r,rho2)
     test += rho2 * final_weight_functions_at_grid_points(l,k,j)
     test_bart += rho2_ap(istate) * final_weight_functions_at_grid_points(l,k,j) 
    enddo
   enddo
  enddo
 print*,'test = ',test
 print*,'test Barth = ',test_bart
 enddo
end

