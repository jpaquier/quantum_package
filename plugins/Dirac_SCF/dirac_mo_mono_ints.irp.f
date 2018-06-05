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



 BEGIN_PROVIDER [ complex*16, dirac_mo_overlap,(2*dirac_mo_tot_num,2*dirac_mo_tot_num)]
  implicit none
  integer :: i,j,n,l
  double precision :: f
  integer :: lmax
  lmax = (2*dirac_ao_num/4) * 4
  do j=1,2*dirac_mo_tot_num
   do i= 1,2*dirac_mo_tot_num
    dirac_mo_overlap(i,j) = 0.d0
    do n = 1, lmax,4
     do l = 1, 2*dirac_ao_num
      dirac_mo_overlap(i,j) += dirac_mo_coef(l,i) * &
                              ( dirac_mo_coef(n  ,j) * dirac_ao_overlap(l,n  )  &
                              + dirac_mo_coef(n+1,j) * dirac_ao_overlap(l,n+1)  &
                              + dirac_mo_coef(n+2,j) * dirac_ao_overlap(l,n+2)  &
                              + dirac_mo_coef(n+3,j) * dirac_ao_overlap(l,n+3)  )
     enddo
    enddo
    do n = lmax+1, 2*dirac_ao_num
     do l = 1, 2*dirac_ao_num
      dirac_mo_overlap(i,j) +=  dirac_mo_coef(n,j) * dirac_mo_coef(l,i) * dirac_ao_overlap(l,n)
     enddo
    enddo
   enddo
  enddo
 END_PROVIDER
 
 BEGIN_PROVIDER [complex*16, dirac_mo_overlap_bis,(2*(dirac_mo_tot_num),2*(dirac_mo_tot_num))]
  implicit none
  BEGIN_DOC
  ! array of the mono electronic hamiltonian on the MOs basis 
  ! obtained from the canonical orthonormalisation of the AOs 
  ! basis, in the 4x4 component formalism with cartesian basis 
  ! and the unrestricted kinetic-balance scheme  
  END_DOC
    call dirac_ao_to_mo(                                          &
        dirac_ao_overlap,                                         &
        size(dirac_ao_overlap,1),                                 &
        dirac_mo_overlap_bis,                                         &
        size(dirac_mo_overlap_bis,1)                                  &
        )
 END_PROVIDER
  



 BEGIN_PROVIDER [double precision, eigenvalues_dirac_mono_elec_mo, (2*(dirac_mo_tot_num))]
 &BEGIN_PROVIDER [complex*16, eigenvectors_dirac_mono_elec_mo, (2*(dirac_mo_tot_num),2*(dirac_mo_tot_num))]
 implicit none
 integer :: n,nmax
 double precision :: eigenvalues( 2*(dirac_mo_tot_num))
 complex*16       :: eigenvectors(2*(dirac_mo_tot_num),2*(dirac_mo_tot_num))
 n = 2*(dirac_mo_tot_num)
 nmax = n
 call lapack_diag_complex(eigenvalues,eigenvectors,dirac_mo_mono_elec_integral,nmax,n)
 eigenvalues_dirac_mono_elec_mo = eigenvalues
 eigenvectors_dirac_mono_elec_mo = eigenvectors
 END_PROVIDER

 
 BEGIN_PROVIDER [complex*16, eigenvectors_dirac_mono_elec_ao, (2*(dirac_mo_tot_num),2*(dirac_mo_tot_num))]
 implicit none
 integer :: n,nmax
  call zgemm('N','N', 2*(dirac_ao_num), 2*(dirac_mo_tot_num), 2*(dirac_ao_num),              &
      (1.d0,0.d0), dirac_mo_coef,size(dirac_mo_coef,1),                                      &
      eigenvectors_dirac_mono_elec_mo, size(eigenvectors_dirac_mono_elec_mo,1),              &
      (0.d0,0.d0), eigenvectors_dirac_mono_elec_ao, size(eigenvectors_dirac_mono_elec_ao,1)) 
 END_PROVIDER



