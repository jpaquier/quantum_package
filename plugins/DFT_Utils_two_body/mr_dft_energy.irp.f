BEGIN_PROVIDER [double precision, electronic_energy_mr_dft, (N_states)]
 implicit none
 BEGIN_DOC
 ! Energy for the multi determinantal DFT calculation
 END_DOC
 
 if(projected_wft_for_dft)then
  print*,'You are using a projected WFT which uses an eigenvalue stored in the EZFIO folder' 
  electronic_energy_mr_dft = data_energy_proj + short_range_Hartree + energy_x + energy_c - Trace_v_Hxc
 else
  print*,'You are using a variational method which uses the wave function stored in the EZFIO folder'
  electronic_energy_mr_dft = total_range_separated_electronic_energy
 endif


END_PROVIDER 

 subroutine print_variational_energy_dft
 implicit none
 print*,'/////////////////////////'
  print*,  '****************************************'
  print*,'///////////////////'
  print*,  ' Regular range separated DFT energy '
  write(*, '(A22,X,F32.10)') 'mu_erf              = ',mu_erf          
  write(*, '(A22,X,F16.10)') 'TOTAL ENERGY        = ',electronic_energy_mr_dft+nuclear_repulsion
  print*, ''
  print*, 'Component of the energy ....'
  print*, ''
  write(*, '(A22,X,F16.10)') 'nuclear_repulsion   = ',nuclear_repulsion
  write(*, '(A22,X,F16.10)') 'psi_energy_erf      = ',psi_energy_erf      
  write(*, '(A22,X,F16.10)') 'psi_energy_core     = ',psi_energy_core    
  write(*, '(A22,X,F16.10)') 'short_range_Hartree = ',short_range_Hartree
  write(*, '(A22,X,F16.10)') 'two_elec_energy     = ',two_elec_energy_dft
  write(*, '(A22,X,F16.10)') 'energy_x            = ',energy_x         
  write(*, '(A22,X,F16.10)') 'energy_c            = ',energy_c          
  print*, ''
  print*,  '****************************************'
  print*, ''
  write(*, '(A22,X,F16.10)') 'Approx eigenvalue   = ',Fock_matrix_expectation_value + psi_energy_erf
  write(*, '(A22,X,F16.10)') 'Trace_v_xc          = ',Trace_v_xc
 
  write(*, '(A28,X,F16.10)') 'Variational energy of Psi = ',psi_energy
  write(*, '(A28,X,F16.10)') 'psi_energy_bielec         = ',psi_energy_bielec
  write(*, '(A28,X,F16.10)') 'psi_energy_monoelec       = ',psi_kinetic_energy+psi_nuclear_elec_energy
 end

subroutine print_projected_energy_dft
 implicit none
 print*,  '****************************************'
 write(*, '(A22,X,F16.10)') 'TOTAL ENERGY        = ',electronic_energy_mr_dft 
 print*, ''
 print*, 'Component of the energy ....'
 print*, ''
 write(*, '(A22,X,F16.10)') 'nuclear_repulsion   = ',nuclear_repulsion
 write(*, '(A22,X,F16.10)') 'projected energy    = ',data_energy_proj + short_range_Hartree + energy_x + energy_c - Trace_v_Hxc
 write(*, '(A22,X,F16.10)') 'short_range_Hartree = ',short_range_Hartree 
 write(*, '(A22,X,F16.10)') 'energy_x            = ',energy_x 
 write(*, '(A22,X,F16.10)') 'energy_c            = ',energy_c 
 write(*, '(A22,X,F16.10)') '- Trace_v_Hxc       = ',- Trace_v_Hxc

end

