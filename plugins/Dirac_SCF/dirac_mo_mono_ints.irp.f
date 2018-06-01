 BEGIN_PROVIDER [complex*16, dirac_mo_mono_elec_integral,(2*(dirac_mo_tot_num),2*(dirac_mo_tot_num))]
  implicit none
  BEGIN_DOC
  ! array of the mono electronic hamiltonian on the MOs basis 
  ! obtained from the canonical orthonormalisation of the AOs 
  ! basis, in the 4x4 component formalism with cartesian basis 
  ! and the unrestricted kinetic-balance scheme  
  END_DOC
    call dirac_ao_to_mo(                                                     &
        dirac_ao_mono_elec_integral,                                         &
        size(dirac_ao_mono_elec_integral,1),                                 &
        dirac_mo_mono_elec_integral,                                         &
        size(dirac_mo_mono_elec_integral,1)                                  &
        )
 END_PROVIDER


 BEGIN_PROVIDER [double precision,dirac_mo_mono_elec_eigenvalues,(2*(dirac_mo_tot_num))]
 &BEGIN_PROVIDER [complex*16, dirac_mo_mono_elec_eigenvectors,(2*(dirac_mo_tot_num),2*(dirac_mo_tot_num))]
 implicit none
 integer :: n,nmax
 double precision :: eigenvalues( 2*(dirac_mo_tot_num))
 complex*16       :: eigenvectors(2*(dirac_mo_tot_num),2*(dirac_mo_tot_num))
 n = 2*(dirac_mo_tot_num)
 nmax = n
 call lapack_diag_complex(eigenvalues,eigenvectors,dirac_mo_mono_elec_integral,nmax,n)
 dirac_mo_mono_elec_eigenvalues = eigenvalues
 dirac_mo_mono_elec_eigenvectors = eigenvectors
 END_PROVIDER

 
