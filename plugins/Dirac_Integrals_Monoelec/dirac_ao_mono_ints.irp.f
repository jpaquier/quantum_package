 BEGIN_PROVIDER [ complex*16, dirac_ao_mono_elec_nucl_integral,(2*dirac_ao_num,2*dirac_ao_num)]
  implicit none
  integer          :: i,i_minus,j,j_minus
  BEGIN_DOC
  !Array of the mono electronic nucleus-electron
  ! hamiltonian on the AOs basis in the 4x4 component
  ! formalism with cartesian basis and  the unrestricted kinetic-balance scheme  
  END_DOC
  dirac_ao_mono_elec_nucl_integral = 0.d0
  do j = 1, 2*(dirac_ao_num)
   if (j .le. large_ao_num) then 
    do i = 1, 2*(dirac_ao_num)
     if (i .le. large_ao_num) then
      dirac_ao_mono_elec_nucl_integral(i,j) = (1.d0,0.d0)*large_ao_nucl_elec_integral(i,j)  
     endif
    enddo
   elseif (j .gt. large_ao_num .and. j .le. 2*large_ao_num) then
    j_minus = j - large_ao_num
    do i = 1, 2*(dirac_ao_num)
     if (i .gt. large_ao_num .and. i .le. 2*large_ao_num) then
     i_minus = i - large_ao_num
      dirac_ao_mono_elec_nucl_integral(i,j) = (1.d0,0.d0)*large_ao_nucl_elec_integral(i_minus,j_minus)
     endif
    enddo
   elseif (j .gt. 2*large_ao_num .and. j .le. (2*large_ao_num+small_ao_num)) then 
    j_minus = j - 2*large_ao_num
    do i = 1, 2*(dirac_ao_num)
     if (i .gt. 2*large_ao_num .and. i .le. (2*large_ao_num+small_ao_num)) then
      i_minus = i - 2*large_ao_num 
      dirac_ao_mono_elec_nucl_integral(i,j) = (1.d0,0.d0)*small_ao_nucl_elec_integral(i_minus,j_minus)
     endif
    enddo
   elseif (j .gt. (2*large_ao_num+small_ao_num)) then
    j_minus = j - (2*large_ao_num + small_ao_num)
    do i = 1, 2*(dirac_ao_num)
     if  (i > (2*large_ao_num+small_ao_num)) then
      i_minus = i - (2*large_ao_num + small_ao_num)
      dirac_ao_mono_elec_nucl_integral(i,j) = (1.d0,0.d0)*small_ao_nucl_elec_integral(i_minus,j_minus)
     endif
    enddo
   endif    
  enddo
 END_PROVIDER


 BEGIN_PROVIDER [ complex*16, dirac_ao_mono_elec_mass_integral,(2*dirac_ao_num,2*dirac_ao_num)]
  implicit none
  integer          :: i,i_minus,j,j_minus
  BEGIN_DOC
  !Array of the mono electronic mass hamiltonian on the AOs basis
  ! in the 4x4 component formalism with cartesian basis and 
  ! the unrestricted kinetic-balance scheme  
  END_DOC
  dirac_ao_mono_elec_mass_integral = 0
  do j = 1, 2*(dirac_ao_num)
   if (j .gt. 2*large_ao_num .and. j .le. (2*large_ao_num+small_ao_num)) then 
    j_minus = j - 2*large_ao_num
    do i = 1, 2*(dirac_ao_num)
     if (i .gt. 2*large_ao_num .and. i .le. (2*large_ao_num+small_ao_num)) then
      i_minus = i - 2*large_ao_num 
      dirac_ao_mono_elec_mass_integral(i,j) = (1.d0,0.d0)* small_ao_mass_energy(i_minus,j_minus)
     endif
    enddo
   elseif (j .gt. (2*large_ao_num+small_ao_num)) then
    j_minus = j - (2*large_ao_num + small_ao_num)
    do i = 1, 2*(dirac_ao_num)
     if  (i > (2*large_ao_num+small_ao_num)) then
      i_minus = i - (2*large_ao_num + small_ao_num)
      dirac_ao_mono_elec_mass_integral(i,j) = (1.d0,0.d0)* small_ao_mass_energy(i_minus,j_minus)
     endif
    enddo
   endif    
  enddo
 END_PROVIDER


 BEGIN_PROVIDER [ complex*16, dirac_ao_mono_elec_kinetic_integral,(2*dirac_ao_num,2*dirac_ao_num)]
  implicit none
  integer          :: i,i_minus,j,j_minus
  BEGIN_DOC
  !Array of the mono electronic kinetic hamiltonian on the AOs basis
  ! in the 4x4 component formalism with cartesian basis and 
  ! the unrestricted kinetic-balance scheme  
  END_DOC
  dirac_ao_mono_elec_kinetic_integral = 0
  do j = 1, 2*(dirac_ao_num)
   if (j .le. large_ao_num) then 
    do i = 1, 2*(dirac_ao_num)
     if (i .gt. 2*large_ao_num .and. i .le. (2*large_ao_num+small_ao_num)) then
      i_minus = i - 2*large_ao_num
      dirac_ao_mono_elec_kinetic_integral(i,j) = dirac_ao_kinetic_integral_z(i_minus,j)
     elseif (i .gt. (2*large_ao_num+small_ao_num)) then
      i_minus = i - (2*large_ao_num + small_ao_num)
      dirac_ao_mono_elec_kinetic_integral(i,j) = dirac_ao_kinetic_integral_plus(i_minus,j)  
     endif
    enddo
   elseif (j .gt. large_ao_num .and. j .le. 2*large_ao_num) then
    j_minus = j - large_ao_num
    do i = 1, 2*(dirac_ao_num)
     if (i .gt. 2*large_ao_num .and. i .le. (2*large_ao_num+small_ao_num)) then
     i_minus = i - 2*large_ao_num
      dirac_ao_mono_elec_kinetic_integral(i,j) = dirac_ao_kinetic_integral_minus(i_minus,j_minus)
     elseif (i > (2*large_ao_num+small_ao_num)) then
      i_minus = i - (2*large_ao_num + small_ao_num)
      dirac_ao_mono_elec_kinetic_integral(i,j) = - dirac_ao_kinetic_integral_z(i_minus,j_minus)  
     endif
    enddo
   elseif (j .gt. 2*large_ao_num .and. j .le. (2*large_ao_num+small_ao_num)) then 
        j_minus = j - 2*large_ao_num
    do i = 1, 2*(dirac_ao_num)
     if (i .le. large_ao_num) then
      dirac_ao_mono_elec_kinetic_integral(i,j) =  Conjg(dirac_ao_kinetic_integral_z(j_minus,i))
     elseif (i .gt. large_ao_num .and. i .le. 2*large_ao_num) then
      i_minus = i - large_ao_num
      dirac_ao_mono_elec_kinetic_integral(i,j) =  Conjg(dirac_ao_kinetic_integral_minus(j_minus,i_minus))
     endif
    enddo
   elseif (j .gt. (2*large_ao_num+small_ao_num)) then
    j_minus = j - (2*large_ao_num + small_ao_num)
    do i = 1, 2*(dirac_ao_num)
     if (i .le. large_ao_num) then
      dirac_ao_mono_elec_kinetic_integral(i,j) =  Conjg(dirac_ao_kinetic_integral_plus(j_minus,i)) 
     elseif (i .gt. large_ao_num .and. i .le. 2*large_ao_num) then
      i_minus = i - large_ao_num
      dirac_ao_mono_elec_kinetic_integral(i,j) = - Conjg(dirac_ao_kinetic_integral_z(j_minus,i_minus))
     endif
    enddo
   endif    
  enddo
 END_PROVIDER


 BEGIN_PROVIDER [ complex*16, dirac_ao_mono_elec_integral,(2*dirac_ao_num,2*dirac_ao_num)]
 &BEGIN_PROVIDER [ complex*16, dirac_ao_mono_elec_integral_diag, (2*dirac_ao_num) ] 
 implicit none
  integer          :: i,j
  BEGIN_DOC
  !Array of the mono electronic hamiltonian on the dirac AO basis
  ! in the 4x4 component formalism with cartesian basis and 
  ! the unrestricted kinetic-balance scheme  
  END_DOC
  print*,'Computing the mono-electronic Fock matrix'
  dirac_ao_mono_elec_integral = (0.d0,0.d0)
  do j = 1, 2*(dirac_ao_num)
   do i = 1, 2*(dirac_ao_num)
    dirac_ao_mono_elec_integral(i,j) += (dirac_ao_mono_elec_nucl_integral(i,j) + dirac_ao_mono_elec_mass_integral(i,j) + dirac_ao_mono_elec_kinetic_integral(i,j) )
   enddo
  enddo
  do j = 1, 2*dirac_ao_num
   dirac_ao_mono_elec_integral_diag(j) = dirac_ao_mono_elec_integral(j,j)
  enddo
 END_PROVIDER
