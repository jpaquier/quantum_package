program projected_operators
  implicit none
  BEGIN_DOC
! TODO
  END_DOC
  read_wf = .True.
  touch read_wf
  call routine

end

subroutine routine
 implicit none
 integer :: ipoint
 double precision :: accu, weight
 do ipoint  = 1, n_points_final_grid
  weight=final_weight_functions_at_final_grid_points(ipoint)
  accu += dabs(f_psi_B_old(ipoint) - f_psi_B(ipoint)) *weight
 enddo
 print*,'accu = ',accu


end
