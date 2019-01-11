 BEGIN_PROVIDER [complex*16, dirac_potential_x_ao,(2*dirac_ao_num,2*dirac_ao_num,N_states)]
 &BEGIN_PROVIDER [complex*16, dirac_potential_c_ao,(2*dirac_ao_num,2*dirac_ao_num,N_states)]
 implicit none
 BEGIN_DOC
 ! relativistic exchange/correlation potentials on the AO basis
 END_DOC

  if(dirac_exchange_functional.EQ."None")then
   dirac_potential_x_ao = (0.d0,0.d0) 
  else
   dirac_potential_x_ao = dirac_potential_x_ao_LDA
  endif
  if(dirac_correlation_functional.EQ."None")then
   dirac_potential_c_ao = (0.d0,0.d0) 
  else
   dirac_potential_c_ao = dirac_potential_c_ao_LDA
  endif
 END_PROVIDER 

 BEGIN_PROVIDER [complex*16, dirac_potential_x_mo,(2*dirac_mo_tot_num,2*dirac_mo_tot_num,N_states)] 
 &BEGIN_PROVIDER [complex*16, dirac_potential_c_mo,(2*dirac_mo_tot_num,2*dirac_mo_tot_num,N_states)]
 implicit none
 BEGIN_DOC
 !  exchange/correlation potentials on the MO basis
 END_DOC
 integer :: istate 
  do istate = 1, N_states
    call dirac_ao_to_mo(                                                   &
        dirac_potential_x_ao(1,1,istate),                                 &
        size(dirac_potential_x_ao,1),                                &
        dirac_potential_x_mo(1,1,istate),                                 &
        size(dirac_potential_x_mo,1)                                 &
        )
    call dirac_ao_to_mo(                                                   &
        dirac_potential_c_ao(1,1,istate),                                 &
        size(dirac_potential_c_ao,1),                                &
        dirac_potential_c_mo(1,1,istate),                                 &
        size(dirac_potential_c_mo,1)                                 &
        )
  enddo
 END_PROVIDER 


 BEGIN_PROVIDER [double precision, dirac_energy_x, (N_states)]
 &BEGIN_PROVIDER [double precision, dirac_energy_c, (N_states)]
 implicit none
  if(dirac_exchange_functional.EQ."None")then
   dirac_energy_x = 0.d0
  else
   dirac_energy_x = dirac_energy_x_LDA
  endif
  if(dirac_correlation_functional.EQ."None")then
   dirac_energy_c = 0.d0
  else
   dirac_energy_c = dirac_energy_c_LDA
  endif  
 END_PROVIDER 
