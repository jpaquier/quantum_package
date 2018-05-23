 BEGIN_PROVIDER [ complex*16, dirac_ao_mono_elec_integral,(2*(dirac_ao_num),2*(dirac_ao_num))]
  implicit none
  integer :: i,j,k,l
  BEGIN_DOC
 ! array of the mono electronic hamiltonian on the AOs basis
 ! in the 4x4 component formalism with cartesian basis and 
 ! the unrestricted kinetic-balance scheme  
  END_DOC
  do j = 1, 2*(dirac_ao_num)
   if (j .le. ao_num) then 
    l = j - 0
    do i = 1, 2*(dirac_ao_num)
     if (i .le. ao_num) then 
      k = i - 0
      dirac_ao_mono_elec_integral(i,j) = (1.d0,0.d0)*ao_nucl_elec_integral(k,l)  
     elseif (i .gt. ao_num .and. i .le. 2*ao_num) then
      k = i - ao_num
      dirac_ao_mono_elec_integral(i,j) = (1.d0,0.d0)*0.d0
     elseif (i .gt. 2*ao_num .and. i .le. (2*ao_num+small_ao_num)) then
      k = i - 2*ao_num
      dirac_ao_mono_elec_integral(i,j) = Conjg(dirac_ao_kinetic_integral_z(k,l))
     else
      k = i - (2*ao_num+small_ao_num)
      dirac_ao_mono_elec_integral(i,j) = Conjg(dirac_ao_kinetic_integral_plus(k,l))  
     endif
    enddo
   elseif (j .gt. ao_num .and. j .le. 2*ao_num) then
    l = j - ao_num
    do i = 1, 2*(ao_num + small_ao_num)
     if (i .le. ao_num) then 
      k = i - 0
      dirac_ao_mono_elec_integral(i,j) = (1.d0,0.d0)*0.d0  
     elseif (i .gt. ao_num .and. i .le. 2*ao_num) then
      k = i - ao_num
      dirac_ao_mono_elec_integral(i,j) = (1.d0,0.d0)*ao_nucl_elec_integral(k,l)
     elseif (i .gt. 2*ao_num .and. i .le. (2*ao_num+small_ao_num)) then
      k = i - 2*ao_num
      dirac_ao_mono_elec_integral(i,j) = Conjg(dirac_ao_kinetic_integral_minus(k,l))
     else
      k = i - (2*ao_num+small_ao_num)
      dirac_ao_mono_elec_integral(i,j) = - Conjg(dirac_ao_kinetic_integral_z(k,l))  
     endif
    enddo
   elseif (j .gt. 2*ao_num .and. j .le. (2*ao_num+small_ao_num)) then 
    l = j - 2*ao_num
    do i = 1, 2*(ao_num + small_ao_num)
     if (i .le. ao_num) then
      k = i - 0 
      dirac_ao_mono_elec_integral(i,j) =  dirac_ao_kinetic_integral_z(l,k) 
     elseif (i .gt. ao_num .and. i .le. 2*ao_num) then
      k = i - ao_num
      dirac_ao_mono_elec_integral(i,j) =  dirac_ao_kinetic_integral_minus(l,k)
     elseif (i .gt. 2*ao_num .and. i .le. (2*ao_num+small_ao_num)) then
      k = i - 2*ao_num
      dirac_ao_mono_elec_integral(i,j) = (1.d0,0.d0)*(small_ao_nucl_elec_integral(k,l) + small_ao_mass_energy(k,l))
     else
      k = i - (2*ao_num+small_ao_num)
      dirac_ao_mono_elec_integral(i,j) = (1.d0,0.d0)*0.d0
     endif
    enddo
   else 
    l = j - (2*ao_num+small_ao_num)
    do i = 1, 2*(ao_num + small_ao_num)
     if (i .le. ao_num) then
      k = i - 0 
      dirac_ao_mono_elec_integral(i,j) =  dirac_ao_kinetic_integral_plus(l,k) 
     elseif (i .gt. ao_num .and. i .le. 2*ao_num) then
      k = i - ao_num
      dirac_ao_mono_elec_integral(i,j) = - dirac_ao_kinetic_integral_z(l,k)
     elseif (i .gt. 2*ao_num .and. i .le. (2*ao_num+small_ao_num)) then
      k = i - 2*ao_num
      dirac_ao_mono_elec_integral(i,j) = (1.d0,0.d0)*0.d0
     else
      k = i - (2*ao_num+small_ao_num)
      dirac_ao_mono_elec_integral(i,j) = (1.d0,0.d0)*(small_ao_nucl_elec_integral(k,l) + small_ao_mass_energy(k,l))
     endif
    enddo
   endif    
  enddo
  END_PROVIDER

