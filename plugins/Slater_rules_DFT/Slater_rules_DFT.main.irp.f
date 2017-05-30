program Slater_rules_DFT
  implicit none
  print*, 'psi_energy_erf      ',psi_energy_erf
  print*, 'psi_energy_core     ',psi_energy_core
  print*, 'psi_energy_hartree  ',psi_energy_hartree
  print*, 'energy_x            ',energy_x
  print*, 'energy_c            ',energy_c
  print*, 'energy total        ',energy_c+energy_x+psi_energy_hartree+psi_energy_core+psi_energy_erf

end
