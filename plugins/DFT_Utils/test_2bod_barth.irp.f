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

subroutine on_top_pair_density_approx(r,rho2_ap)
 implicit none
 double precision, intent(in)  :: r(3)
 double precision, intent(out) :: rho2_ap(N_states)
 double precision :: psi_temp
 psi_temp = 0.d0
 double precision :: mos_array(mo_tot_num)
 double precision :: aos_array(mo_tot_num)
 call give_all_mos_at_r(r,mos_array)
 integer :: i,j,k,l,t,m,istate
 do istate = 1, N_states
 rho2_ap(istate) = 0.d0
 do m = 1, n_couple(istate) ! loop over the first electron 
   ! get matrix 
   i=list_couple(m,1,istate)
   j=list_couple(m,2,istate)
   print*,'i,j',i,j
    do k = 1, n_k_dens(m,istate) 
     do t=1,mo_tot_num
      psi_temp+=psi_k_couple(t,k,m,istate)*mos_array(t)
     enddo
     rho2_ap(istate) +=mos_array(i)*mos_array(j)*lambda_k(k,m,istate)*psi_temp*psi_temp
     print*,'rho2(0) barth= ',rho2_ap(istate)
   enddo
 enddo
 enddo
end


subroutine test_rho2
 implicit none
 integer :: j,k,l,istate,n_k_loc,m
 double precision :: r(3),rho2
 double precision :: rho2_ap(N_states)
 double precision :: test,test_bart

!do istate = 1, N_states
! do m = 1, n_couple(istate)
!  do n_k_loc=1, n_k_dens(m,istate)

!  print*,'lambda: ',lambda_k(n_k_loc,m,istate)
!  do j=1,mo_tot_num
!  print*,'j pote: ',psi_k_couple(j,n_k_loc,m,istate)
!  enddo
!  enddo
! enddo
!enddo

call  on_top_pair_density_approx(r,rho2_ap)
!do istate = 1, N_states
! r(1) = 0.d0
! r(2) = 0.d0
! r(3) = 0.d0
 print*,'rho2(0) barth= ',rho2_ap(1)
!enddo


!do istate = 1, N_states
! r(1) = 0.d0
! r(2) = 0.d0
! r(3) = 0.d0
! call  on_top_pair_density_in_real_space(r,rho2)
! print*,'rho2(0) = ',rho2
! call  on_top_pair_density_approx(r,rho2_ap)
! print*,'rho2(0) barth= ',rho2_ap(istate)
!stop
! test = 0.d0

!  do j = 1, nucl_num
!   do k = 1, n_points_radial_grid  -1
!    print*,k
!    do l = 1, n_points_integration_angular
!
!     r(1) = grid_points_per_atom(1,l,k,j)
!     r(2) = grid_points_per_atom(2,l,k,j)
!     r(3) = grid_points_per_atom(3,l,k,j)
!     call  on_top_pair_density_in_real_space(r,rho2)
!     call  on_top_pair_density_approx(r,rho2_ap(istate))
!     print*,'rho2(r) normal = ',rho2
!!     print*,'rho2(r) barth =  ',rho2_ap(istate)
!    call  on_top_pair_density_in_real_space_from_ao(r,rho2)
!     test += rho2 * final_weight_functions_at_grid_points(l,k,j)
!     test_bart += rho2_ap(istate) * final_weight_functions_at_grid_points(l,k,j) 
!     enddo
!    enddo
!   enddo
! print*,'test = ',test
! print*,'test Barth = ',test_bart
! enddo
end

