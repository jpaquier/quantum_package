 BEGIN_PROVIDER [complex*16, dirac_ao_potential_xc, (2*dirac_ao_num, 2*dirac_ao_num)]
 implicit none
 integer :: i,j,k,l
 dirac_ao_potential_xc = (0.d0,0.d0)
  do i = 1, 2*dirac_ao_num
   do j = 1, 2*dirac_ao_num
 !  dirac_ao_potential_xc(i,j) =  dirac_potential_c_ao(i,j,1) + dirac_potential_x_ao(i,j,1)
   enddo
  enddo
 END_PROVIDER

 BEGIN_PROVIDER [double precision, dirac_e_exchange_dft]
 implicit none
  dirac_e_exchange_dft = dirac_energy_x(1)
 END_PROVIDER

 BEGIN_PROVIDER [double precision, dirac_e_correlation_dft]
 implicit none
  dirac_e_correlation_dft = dirac_energy_c(1)
 END_PROVIDER

