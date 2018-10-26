program print_mono_energy
 implicit none
 read_wf = .True.
 touch read_wf
 call routine
end

subroutine routine
 implicit none
 integer :: i,j
 double precision :: average
 call get_average(mo_spread_x,one_body_dm_mo,average)
 print*, 'density mono energy = ',average
 print*, '**********************'
!print*, 'psi_energy_bielec   = ',psi_energy_bielec(1)
!print*, 'psi_energy_monoelec = ',psi_energy_monoelec(1)
!print*, 'total energy        = ',psi_energy_monoelec(1) + psi_energy_bielec(1)
!print*, 'total energy + nucl = ',psi_energy_monoelec(1) + psi_energy_bielec(1) + nuclear_repulsion
!print*, 'psi_energy          = ',psi_energy(1)
!print*, 'psi_energy + nucl   = ',psi_energy(1) + nuclear_repulsion
 

end
