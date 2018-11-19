program projected_operators
  implicit none
  BEGIN_DOC
! TODO
  END_DOC
  read_wf = .True.
  touch read_wf
  call routine_v
  call routine_rho 
  call routine_final

end

subroutine routine_rho
 implicit none
 provide rho2_kl_contracted_transposed

end
subroutine routine_v
 implicit none
 integer :: ipoint,k,l
 double precision :: accu, weight
 accu = 0.d0
 do ipoint  = 1, n_points_final_grid
  weight=final_weight_functions_at_final_grid_points(ipoint)
  do l = 1, mo_tot_num
   do k = 1, mo_tot_num
    accu += dabs(V_kl_contracted(k,l,ipoint) ) * weight
   enddo
  enddo
 enddo
 print*,'accu = ',accu
end
 

subroutine routine_final
 implicit none
 integer :: ipoint
 double precision :: accu, weight
 accu = 0.d0
 do ipoint  = 1, n_points_final_grid
  weight=final_weight_functions_at_final_grid_points(ipoint)
  accu += (f_psi_B(ipoint)) *weight
 enddo
 print*,'*******************'
 print*,'*******************'
 print*,'*******************'
 print*,'*******************'
 print*,'accu = ',accu


end
