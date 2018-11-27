 BEGIN_PROVIDER[complex*16, dirac_fock_matrix_ao, (2*dirac_ao_num,2*dirac_ao_num) ]
  implicit none
  BEGIN_DOC
  !Dirac Fock matrix in AO basis set 
  END_DOC
  integer                        ::i,j
  if (dirac_interaction == "Coulomb") then
   dirac_fock_matrix_ao = dirac_Fock_matrix_C_ao
  elseif (dirac_interaction == "Coulomb_Gaunt") then
   dirac_fock_matrix_ao = dirac_Fock_matrix_C_G_ao
  else
   print *,  'Unrecognized dirac_interaction : '//dirac_interaction
   stop 1
  endif
 END_PROVIDER

 BEGIN_PROVIDER[complex*16, eigenvectors_dirac_fock_matrix_ao, (2*dirac_mo_tot_num,2*dirac_mo_tot_num) ]
  implicit none
  BEGIN_DOC
  !Dirac Fock eigenvectors in AO basis set
  END_DOC
  integer                        ::i,j
  if (dirac_interaction == "Coulomb") then
   eigenvectors_dirac_fock_matrix_ao = eigenvectors_dirac_Fock_matrix_C_ao
  elseif (dirac_interaction == "Coulomb_Gaunt") then
   eigenvectors_dirac_fock_matrix_ao = eigenvectors_dirac_Fock_matrix_C_G_ao
  else
   print *,  'Unrecognized dirac_interaction : '//dirac_interaction
   stop 1
  endif
 END_PROVIDER

 BEGIN_PROVIDER[double precision, dirac_SCF_energy ]
  implicit none
  BEGIN_DOC
  !Dirac_SCF energy
  END_DOC
  integer                        ::i,j
  if (dirac_interaction == "Coulomb") then
   dirac_SCF_energy = dirac_SCF_C_energy
  elseif (dirac_interaction == "Coulomb_Gaunt") then
   dirac_SCF_energy = dirac_SCF_C_G_energy
  else
   print *,  'Unrecognized dirac_interaction : '//dirac_interaction
   stop 1
  endif
 END_PROVIDER

