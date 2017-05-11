program print_energy
 implicit none
 read_wf = .true.
 touch read_wf
 call routine
end

subroutine routine
 print*, 'psi_energy_bielec   = ',psi_energy_bielec
 print*, 'psi_energy_monoelec = ',psi_energy_monoelec
 print*, 'total energy        = ',psi_energy_monoelec + psi_energy_bielec
 print*, 'psi_energy          = ',psi_energy
end
