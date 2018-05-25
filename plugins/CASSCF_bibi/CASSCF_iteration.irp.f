subroutine iteration_scf_state_average(delta_e_sci)
 implicit none
 double precision, intent(out) :: delta_e_sci
 call initialize_mo_coef_begin_iteration
 print*, 'Delta E SUPERCI  = ',eigenvalues_sci_state_average
 delta_e_sci = eigenvalues_sci_state_average
 
 call set_superci_natural_mos_state_average
 call reorder_active_orb
end


subroutine iteration_scf_state_specific(delta_e_sci)
 implicit none
 double precision, intent(out) :: delta_e_sci
 call initialize_mo_coef_begin_iteration
 print*, 'Delta E SUPERCI  = ',eigenvalues_sci(1)
 delta_e_sci = eigenvalues_sci(1)
 
 call set_superci_natural_mos_state_specific
 call reorder_active_orb
end


subroutine casscf_routine
 implicit none
 integer :: i,niter,m,j
 double precision :: energy(N_states),thresh_casscf,delta_e(N_states),delta_e_sci
 energy(1) = 0.d0
 thresh_casscf = 1.d-6
 do i = 1, 100
  print*, 'Iteration  = ',i
  do m = 1, N_states
   print*, 'State ',m
   print*, 'Reference Energy = ',m,reference_energy_superci(m)
  enddo

! call iteration_scf_state_specific(delta_e_sci)
  call iteration_scf_state_average(delta_e_sci)


  touch mo_coef
  do m = 1, N_states
   delta_e(m) = reference_energy_superci(m) - energy(m)
   print*, 'delta_e(m) ',delta_e(m),reference_energy_superci(m),energy(m)
  enddo
   if (dabs(delta_e(1)).lt.thresh_casscf)then
    exit
   endif
  do m = 1, N_states
   energy(m) = reference_energy_superci(m)
  enddo
 enddo
 niter = i
 
 print*, '*******************'
 print*, 'SUPER CI converged in ',niter
 do m = 1, N_states
  print*,  'Final Energy     = ',reference_energy_superci(m)
 enddo
 call diag_inactive_virt_and_update_mos_MR_Fock

!call save_mos

end
