
 BEGIN_PROVIDER [double precision, potential_x_alpha_ao,(ao_num,ao_num,N_states)]
&BEGIN_PROVIDER [double precision, potential_x_beta_ao,(ao_num,ao_num,N_states)]
&BEGIN_PROVIDER [double precision, potential_c_alpha_ao,(ao_num,ao_num,N_states)]
&BEGIN_PROVIDER [double precision, potential_c_beta_ao,(ao_num,ao_num,N_states)]
 implicit none
 BEGIN_DOC
! alpha/beta exchange/correlation potentials on the AO basis
 END_DOC

  if(trim(exchange_functional)=="short_range_LDA")then
   potential_x_alpha_ao = potential_x_alpha_ao_LDA
   potential_x_beta_ao = potential_x_beta_ao_LDA
  else if(exchange_functional.EQ."short_range_PBE")then
   potential_x_alpha_ao = potential_x_alpha_ao_PBE
   potential_x_beta_ao = potential_x_beta_ao_PBE
  else if(exchange_functional.EQ."None")then
   potential_x_alpha_ao = 0.d0 
   potential_x_beta_ao = 0.d0 
  else 
   print*, 'Exchange functional required does not exist ...'
   print*,'exchange_functional',exchange_functional
   stop
  endif

  if(trim(correlation_functional)=="short_range_LDA")then
   potential_c_alpha_ao = potential_c_alpha_ao_LDA
   potential_c_beta_ao = potential_c_beta_ao_LDA
  else if(correlation_functional.EQ."short_range_PBE")then
   potential_c_alpha_ao = potential_c_alpha_ao_PBE
   potential_c_beta_ao = potential_c_beta_ao_PBE
  else if(correlation_functional.EQ."None")then
   potential_c_alpha_ao = 0.d0 
   potential_c_beta_ao = 0.d0 
  else 
   print*, 'Correlation functional required does not ecist ...'
   print*,'correlation_functional',correlation_functional
   stop
  endif


END_PROVIDER 





 BEGIN_PROVIDER [double precision, potential_x_alpha_mo,(mo_tot_num,mo_tot_num,N_states)]
&BEGIN_PROVIDER [double precision, potential_x_beta_mo,(mo_tot_num,mo_tot_num,N_states)]
&BEGIN_PROVIDER [double precision, potential_c_alpha_mo,(mo_tot_num,mo_tot_num,N_states)]
&BEGIN_PROVIDER [double precision, potential_c_beta_mo,(mo_tot_num,mo_tot_num,N_states)]
 implicit none
 BEGIN_DOC
! alpha/beta exchange/correlation potentials on the MO basis
 END_DOC
 integer :: istate 
 do istate = 1, N_states
    call ao_to_mo(                                                   &
        potential_x_alpha_ao(1,1,istate),                                 &
        size(potential_x_alpha_ao,1),                                &
        potential_x_alpha_mo(1,1,istate),                                 &
        size(potential_x_alpha_mo,1)                                 &
        )

    call ao_to_mo(                                                   &
        potential_x_beta_ao(1,1,istate),                                  &
        size(potential_x_beta_ao,1),                                 &
        potential_x_beta_mo(1,1,istate),                                  &
        size(potential_x_beta_mo,1)                                  &
        )


    call ao_to_mo(                                                   &
        potential_c_alpha_ao(1,1,istate),                                 &
        size(potential_c_alpha_ao,1),                                &
        potential_c_alpha_mo(1,1,istate),                                 &
        size(potential_c_alpha_mo,1)                                 &
        )

    call ao_to_mo(                                                   &
        potential_c_beta_ao(1,1,istate),                                  &
        size(potential_c_beta_ao,1),                                 &
        potential_c_beta_mo(1,1,istate),                                  &
        size(potential_c_beta_mo,1)                                  &
        )

 enddo

END_PROVIDER 


  BEGIN_PROVIDER [double precision, energy_x, (N_states)]
 &BEGIN_PROVIDER [double precision, energy_c, (N_states)]
 implicit none
  if(trim(exchange_functional)=="short_range_LDA")then
   energy_x = energy_x_LDA
   energy_x = energy_x_LDA
  else if(exchange_functional.EQ."short_range_PBE")then
   energy_x = energy_x_PBE
   energy_x = energy_x_PBE
  else if(exchange_functional.EQ."None")then
   energy_x = 0.d0 
   energy_x = 0.d0 
  else 
   print*, 'Exchange functional required does not exist ...'
   print*,'exchange_functional',exchange_functional
   stop
  endif

  if(trim(correlation_functional)=="short_range_LDA")then
   energy_c = energy_c_LDA
   energy_c = energy_c_LDA
  else if(correlation_functional.EQ."short_range_PBE")then
   energy_c = energy_c_PBE
   energy_c = energy_c_PBE
  else if(correlation_functional.EQ."None")then
   energy_c = 0.d0 
   energy_c = 0.d0 
  else 
   print*, 'Correlation functional required does not ecist ...'
   print*,'correlation_functional',correlation_functional
   stop
  endif
 
END_PROVIDER 
