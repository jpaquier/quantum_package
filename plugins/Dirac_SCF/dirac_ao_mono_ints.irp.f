 BEGIN_PROVIDER [ complex*16, dirac_ao_mono_elec_integral,(2*(dirac_ao_num),2*(dirac_ao_num))]
  implicit none
  integer          :: i,j
  BEGIN_DOC
 ! array of the mono electronic hamiltonian on the AOs basis
 ! in the 4x4 component formalism with cartesian basis and 
 ! the unrestricted kinetic-balance scheme  
  END_DOC
  do i = 1, 2*(dirac_ao_num)
   if (i .le. ao_num) then 
    do j = 1, 2*(dirac_ao_num)
     if (j .le. ao_num) then
      dirac_ao_mono_elec_integral(i,j) = (1.d0,0.d0)*ao_nucl_elec_integral(i,j)  
     elseif (j .gt. ao_num .and. j .le. 2*ao_num) then
      dirac_ao_mono_elec_integral(i,j) = (1.d0,0.d0)*0.d0
     elseif (j .gt. 2*ao_num .and. j .le. (2*ao_num+small_ao_num)) then
      dirac_ao_mono_elec_integral(i,j) = Conjg(dirac_ao_kinetic_integral_z(d_L(j),d_L(i)))
     elseif (j .gt. (2*ao_num+small_ao_num)) then
      dirac_ao_mono_elec_integral(i,j) = Conjg(dirac_ao_kinetic_integral_plus(d_L(j),d_L(i)))  
     endif
    enddo
   elseif (i .gt. ao_num .and. i .le. 2*ao_num) then
    do j = 1, 2*(dirac_ao_num)
     if (j .le. ao_num) then 
      dirac_ao_mono_elec_integral(i,j) = (1.d0,0.d0)*0.d0  
     elseif (j .gt. ao_num .and. j .le. 2*ao_num) then
      dirac_ao_mono_elec_integral(i,j) = (1.d0,0.d0)*ao_nucl_elec_integral(d_L(i),d_L(j))
     elseif (j .gt. 2*ao_num .and. j .le. (2*ao_num+small_ao_num)) then
      dirac_ao_mono_elec_integral(i,j) = Conjg(dirac_ao_kinetic_integral_minus(d_L(j),d_L(i)))
     elseif (j > (2*ao_num+small_ao_num)) then
      dirac_ao_mono_elec_integral(i,j) = - Conjg(dirac_ao_kinetic_integral_z(d_L(j),d_L(i)))  
     endif
    enddo
   elseif (i .gt. 2*ao_num .and. i .le. (2*ao_num+small_ao_num)) then 
    do j = 1, 2*(dirac_ao_num)
     if (j .le. ao_num) then
      dirac_ao_mono_elec_integral(i,j) =  dirac_ao_kinetic_integral_z(d_L(i),d_L(j)) 
     elseif (j .gt. ao_num .and. i .le. 2*ao_num) then
      dirac_ao_mono_elec_integral(i,j) =  dirac_ao_kinetic_integral_minus(d_L(i),d_L(j))
     elseif (j .gt. 2*ao_num .and. j .le. (2*ao_num+small_ao_num)) then
      dirac_ao_mono_elec_integral(i,j) = (1.d0,0.d0)*(small_ao_nucl_elec_integral(d_L(i),d_L(j)) + small_ao_mass_energy(d_L(i),d_L(j)))
     elseif (j > (2*ao_num+small_ao_num)) then
      dirac_ao_mono_elec_integral(i,j) = (1.d0,0.d0)*0.d0
     endif
    enddo
   elseif (i .gt. (2*ao_num+small_ao_num)) then
    do j = 1, 2*(dirac_ao_num)
     if (j .le. ao_num) then
      dirac_ao_mono_elec_integral(i,j) =  dirac_ao_kinetic_integral_plus(d_L(i),d_L(j)) 
     elseif (j .gt. ao_num .and. j .le. 2*ao_num) then
      dirac_ao_mono_elec_integral(i,j) = - dirac_ao_kinetic_integral_z(d_L(i),d_L(j))
     elseif (j .gt. 2*ao_num .and. j .le. (2*ao_num+small_ao_num)) then
      dirac_ao_mono_elec_integral(i,j) = (1.d0,0.d0)*0.d0
     elseif  (j > (2*ao_num+small_ao_num)) then
      dirac_ao_mono_elec_integral(i,j) = (1.d0,0.d0)*(small_ao_nucl_elec_integral(d_L(i),d_L(j)) + small_ao_mass_energy(d_L(i),d_L(j)))
     endif
    enddo
   endif    
  enddo
  END_PROVIDER

