program pouet
 implicit none
 read_wf = .true.
 touch read_wf
 disk_access_mo_one_integrals = "None"
 touch disk_access_only_mo_one_integrals
 disk_access_mo_integrals = "None"
 touch disk_access_mo_integrals
 disk_access_ao_integrals = "None"
 touch disk_access_ao_integrals

 call routine2

end
subroutine routine 
 implicit none
 write(*, '(A31,X,F16.10)') 'psi_energy_core              = ',psi_energy_core   

end

subroutine routine2
 implicit none

 print*, '****************************************'
 write(*, '(A31,X,F16.10)') 'mu_erf                       = ',mu_erf
 print*, '****************************************'

 write(*, '(A31,X,F16.10)') 'TOTAL ENERGY                 = ',total_electronic_energy + nuclear_repulsion
 print*, ''
 print*, 'Component of the energy ....'
 print*, ''
 write(*, '(A31,X,F16.10)') 'nuclear_repulsion            = ',nuclear_repulsion
 write(*, '(A31,X,F16.10)') 'psi_energy_erf               = ',psi_energy_erf      
 write(*, '(A31,X,F16.10)') 'psi_energy_core              = ',psi_energy_core    
 write(*, '(A31,X,F16.10)') 'short_range_Hartree          = ',short_range_Hartree
 write(*, '(A31,X,F16.10)') 'two_elec_energy              = ',two_elec_energy_dft
 write(*, '(A31,X,F16.10)') 'psi_energy_core_and_hartree  = ',short_range_Hartree+psi_energy_core
 write(*, '(A31,X,F16.10)') 'energy_x                     = ',energy_x         
 write(*, '(A31,X,F16.10)') 'energy_c                     = ',energy_c          
 print*, ''
 print*,  '****************************************'
 print*, ''
 write(*, '(A31,X,F16.10)') 'Approx eigenvalue            = ',Fock_matrix_expectation_value + psi_energy_erf
 write(*, '(A31,X,F16.10)') 'Trace_v_xc                   = ',Trace_v_xc


end
