program Slater_rules_DFT
  implicit none
  print*, 'psi_energy_erf      ',psi_energy_erf
  print*, 'psi_energy_core     ',psi_energy_core
  print*, 'short_range_Hartree ',short_range_Hartree
  print*, 'energy_x            ',energy_x
  print*, 'energy_c            ',energy_c
  print*, 'energy total        ',energy_c+energy_x+short_range_Hartree+psi_energy_core+psi_energy_erf

end
