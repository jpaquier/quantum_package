program test_RS_KS
 implicit none
 read_wf = .False.
 touch read_wf

 call test_energy 
 call test_bis
end

subroutine test_energy
 implicit none
 print*, 'exchange_functional    =  ',exchange_functional
 print*, 'correlation_functional =  ',correlation_functional
 print*, 'SCF_energy          = ',SCF_energy
 print*, 'RS_KS_energy        = ',RS_KS_energy
 print*, 'one_electron_energy = ',one_electron_energy
 print*, 'two_electron_energy = ',two_electron_energy
 print*, 'e_exchange_dft      = ',e_exchange_dft
 print*, 'e_correlation_dft   = ',e_correlation_dft
 print*, 'Fock_matrix_energy  = ',Fock_matrix_energy
 print*, 'trace_potential_xc  = ',trace_potential_xc
 print*, 'extra_energy_contrib=  ',extra_energy_contrib_from_density
end

subroutine test_bis

implicit none
integer :: i,j
 print*,'energy_c    = ',energy_c
 print*,'energy_x    = ',energy_x
 do i = 1, ao_num
  do j = 1, ao_num
    print*,Fock_matrix_ao_alpha(i,j),Fock_matrix_ao_beta(i,j)
  enddo
 enddo


end
