program osoci_program
implicit none
 double precision :: norm_total(N_States)
!  do_it_perturbative = .True.
!  touch do_it_perturbative
!  threshold_perturbative = 1.d0
!  touch threshold_perturbative
   call FOBOCI_lmct_mlct_old_thr(1,norm_total)
   call provide_all_the_rest
end
subroutine provide_all_the_rest
implicit none
integer :: i
   call update_one_body_dm_mo
   call set_lmct_mlct_to_psi_det
   call diagonalize_CI
   call save_wavefunction


end

