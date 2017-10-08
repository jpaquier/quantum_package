subroutine iteration_scf(delta_e)
 implicit none
 integer :: m
 double precision, intent(out) :: delta_e(N_states)
 call initialize_mo_coef_begin_iteration
 do m = 1, N_states
  print*, 'State ',m
  print*, 'SUPERCI   Energy = ',eigenvalues_sci(m)+reference_energy_superci(m)
  print*, 'Delta E SUPERCI  = ',eigenvalues_sci(m)
 enddo
 delta_e = eigenvalues_sci
 
 call set_superci_natural_mos 
 touch mo_coef
 call reorder_active_orb
end


subroutine casscf_routine
 implicit none
 integer :: i,niter,m
 double precision :: energy(N_states),thresh_casscf,delta_e(N_states)
 energy(1) = 0.d0
 thresh_casscf = 1.d-10
 do i = 1, 100
  print*, 'Iteration  = ',i
  do m = 1, N_states
   print*, 'State ',m
   print*, 'Reference Energy = ',i,reference_energy_superci(m)
  enddo
  call iteration_scf(delta_e)
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
 call save_mos

end
