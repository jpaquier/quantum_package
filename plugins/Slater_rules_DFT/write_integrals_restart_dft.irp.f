program write_integrals
 implicit none
 read_wf = .true.
 touch read_wf
 disk_access_mo_one_integrals = "None"
 touch disk_access_only_mo_one_integrals
 disk_access_mo_integrals = "None"
 touch disk_access_mo_integrals
 disk_access_ao_integrals = "None"
 touch disk_access_ao_integrals
 call routine
 call routine2

end

subroutine routine
 implicit none
 call save_one_e_effective_potential
 call save_erf_bi_elec_integrals_ao
 call save_erf_bi_elec_integrals_mo

end

subroutine routine2
 implicit none

 print*,  '****************************************'
 write(*, '(A22,X,F16.10)') 'TOTAL ENERGY        = ',total_electronic_energy+nuclear_repulsion
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


end
