program print_ecmd_energy
 implicit none
 read_wf = .true.
 touch read_wf
 disk_access_mo_one_integrals = "None"
 touch disk_access_only_mo_one_integrals
 disk_access_mo_integrals = "None"
 touch disk_access_mo_integrals
 disk_access_ao_integrals = "None"
 touch disk_access_ao_integrals
 call ecmd_energy_printer
!call pouet 
end


subroutine pouet
 implicit none

!write(*, '(A28,X,F16.10)') 'on_top Multi-det correl   = ',Energy_c_md_on_top(1)
!write(*, '(A22,X,F16.10)') 'EC_MD_ON_TOP_PBE_cor= ',Energy_c_md_on_top_PBE_mu_UEG(1)+psi_energy+nuclear_repulsion
!write(*, '(A28,X,F16.10)') 'on_top_PBE_cor MD correl  = ',Energy_c_md_on_top_PBE_mu_UEG(1)
 write(*, '(A28,X,F16.10)') 'on_top_PBE_cor MD bis     = ',Energy_c_md_on_top_PBE_mu_UEG_vector(1)
!write(*, '(A28,X,F16.10)') 'on_top_PBE_cor MD correl  = ',Energy_c_md_on_top_PBE_mu(1)
 write(*, '(A28,X,F16.10)') 'on_top_PBE_cor MD bis     = ',Energy_c_md_on_top_PBE_mu_vector(1)

end

subroutine ecmd_energy_printer
 implicit none
 
 print*,  '****************************************'
 write(*, '(A22,X,F32.10)') 'mu_erf              = ',mu_erf          
 print*,  ' MR DFT energy with pure correlation part for the DFT '
 write(*, '(A22,X,F16.10)') 'EC_MD_LDA           = ',Energy_c_md+psi_energy+nuclear_repulsion
 write(*, '(A22,X,F16.10)') 'EC_MD_ON_TOP_UEG    = ',Energy_c_md_on_top_PBE_mu_UEG_vector(1)+psi_energy+nuclear_repulsion
 write(*, '(A22,X,F16.10)') 'EC_MD_ON_TOP_NO UEG = ',Energy_c_md_on_top_PBE_mu_vector(1)+psi_energy+nuclear_repulsion
 print*, ''
 print*, 'Component of the energy ....'
 print*, ''
 write(*, '(A28,X,F16.10)') 'nuclear_repulsion         = ',nuclear_repulsion
 write(*, '(A28,X,F16.10)') 'Variational energy of Psi = ',psi_energy
 write(*, '(A28,X,F16.10)') 'psi_energy_bielec         = ',psi_energy_bielec
 write(*, '(A28,X,F16.10)') 'LDA Multi-det correlation = ',Energy_c_md
 write(*, '(A28,X,F16.10)') 'on_top_PBE_cor UEG        = ',Energy_c_md_on_top_PBE_mu_UEG_vector(1)
 write(*, '(A28,X,F16.10)') 'on_top_PBE_cor NO UEG     = ',Energy_c_md_on_top_PBE_mu_vector(1)

end
