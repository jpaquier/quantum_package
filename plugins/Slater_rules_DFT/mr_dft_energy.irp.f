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

 print*,  '****************************************'
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

 print*,  '****************************************'
 print*,  '****************************************'
 print*,  ' MR DFT energy with pure correlation part for the DFT '
 write(*, '(A22,X,F16.10)') 'TOTAL ENERGY CORR   = ',elec_energy_dft_pure_corr_funct+nuclear_repulsion
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 write(*, '(A22,X,F16.10)') 'CORRECTED E_TOT CORR= ',Energy_c_md_on_top(1)+psi_energy+nuclear_repulsion
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 print*, ''
 print*, 'Component of the energy ....'
 print*, ''
 write(*, '(A28,X,F16.10)') 'nuclear_repulsion         = ',nuclear_repulsion
 write(*, '(A28,X,F16.10)') 'Variational energy of Psi = ',psi_energy
 write(*, '(A28,X,F16.10)') 'psi_energy_bielec         = ',psi_energy_bielec
 write(*, '(A28,X,F16.10)') 'psi_energy_monoelec       = ',psi_energy_monoelec
 write(*, '(A28,X,F16.10)') 'DFT Multi-det correlation = ',Energy_c_md
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 write(*, '(A28,X,F16.10)') 'corrected Multi-det correl= ',Energy_c_md_on_top(1)
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


end


subroutine print_variational_energy_dft_mu_of_r
 implicit none

 print*,  '****************************************'
 print*, 'Functional used = ',md_correlation_functional
 print*,  '****************************************'
 print*,  ' MR DFT energy with pure correlation part for the DFT '
 if(md_correlation_functional.EQ."basis_set_LDA")then
   write(*, '(A28,X,F16.10)') 'TOTAL ENERGY CORR         = ',psi_energy+Energy_c_md_mu_of_r_LDA+nuclear_repulsion
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   print*, ''
   write(*, '(A28,X,F16.10)') 'Variational energy of Psi = ',psi_energy
   print*, 'Component of the energy ....'
   print*, ''
   write(*, '(A28,X,F16.10)') 'nuclear_repulsion         = ',nuclear_repulsion
   write(*, '(A28,X,F16.10)') 'psi_energy_bielec         = ',psi_energy_bielec
   write(*, '(A28,X,F16.10)') 'psi_energy_monoelec       = ',psi_energy_monoelec
   write(*, '(A28,X,F16.10)') 'DFT mu(r)     correlation = ',Energy_c_md_mu_of_r_LDA
   print*, ''
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 else if(md_correlation_functional.EQ."basis_set_on_top_PBE")then
   write(*, '(A28,X,F16.10)') 'TOTAL ENERGY CORR         = ',psi_energy+Energy_c_md_mu_of_r_PBE_on_top+nuclear_repulsion
   write(*, '(A28,X,F16.10)') 'TOTAL ENERGY CORR CORRECTED=',psi_energy+Energy_c_md_mu_of_r_PBE_on_top_corrected+nuclear_repulsion
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   print*, ''
   write(*, '(A28,X,F16.10)') 'Variational energy of Psi = ',psi_energy
   print*, 'Component of the energy ....'
   print*, ''
   write(*, '(A28,X,F16.10)') 'nuclear_repulsion         = ',nuclear_repulsion
   write(*, '(A28,X,F16.10)') 'psi_energy_bielec         = ',psi_energy_bielec
   write(*, '(A28,X,F16.10)') 'psi_energy_monoelec       = ',psi_energy_monoelec
   write(*, '(A28,X,F16.10)') 'DFT mu(r) correlation     = ',Energy_c_md_mu_of_r_PBE_on_top
   write(*, '(A28,X,F16.10)') 'DFT mu(r) correlation corr= ',Energy_c_md_mu_of_r_PBE_on_top_corrected
   print*, ''
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


 endif 
   write(*, '(A28,X,F16.10)') 'mu_average for basis set  = ',mu_average


end



subroutine print_variational_energy_dft_no_ecmd
 implicit none

 print*,  '****************************************'
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
 print*, ''
 print*,  '****************************************'


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
