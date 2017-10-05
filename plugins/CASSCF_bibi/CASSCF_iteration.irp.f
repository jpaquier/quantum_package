subroutine iteration_scf(delta_e)
 implicit none
 double precision, intent(out) :: delta_e(N_states)
 print*, 'SUPERCI   Energy = ',eigenvalues_sci(1)+reference_energy_superci
 print*, 'Delta E SUPERCI  = ',eigenvalues_sci(1)
 delta_e = eigenvalues_sci
 
 call set_superci_natural_mos 
 touch mo_coef
end


subroutine casscf_routine
 implicit none
 integer :: i,niter
 double precision :: energy(N_states),thresh_casscf,delta_e(N_states)
 energy(1) = 0.d0
 thresh_casscf = 1.d-10
 do i = 1, 100
  print*, 'Iteration  = ',i
  print*, 'Reference Energy = ',i,reference_energy_superci
  call iteration_scf(delta_e)
  if (dabs(delta_e(1)).lt.thresh_casscf)then
   niter = i
   exit
  endif
  energy(1) = reference_energy_superci
 enddo
 
 print*, '*******************'
 print*, 'SUPER CI converged in ',niter
 print*,  'Final Energy     = ',reference_energy_superci
 call save_mos

end
