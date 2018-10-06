program test
 read_wf = .True.
 touch read_wf
 call test_subou
!call test_potepote
end


subroutine test_subou
 implicit none
 double precision :: someee,nbE,nbE_core,nbE_valence
 integer :: j,i
 character*(128) :: output
 character*(128) :: filename
 double precision :: Ecmd, mu_moy
 integer :: i_unit_output,getUnitAndOpen
 output=trim(ezfio_filename)//'_mu_of_r_info'
 output=trim(output)
 print*,'*****************************'
 print*,'*****************************'
 print*,'*****************************'
 print*,'*****************************'
 print*,'output = ',trim(output)
 i_unit_output = getUnitAndOpen(output,'w')
 print*,'i_unit_output',i_unit_output
 print*,'*****************************'
 print*,'*****************************'
 print*,'*****************************'


 nbE=0.0
 nbE_core=0.0
 nbE_valence=0.0
 do i = 1, n_points_radial_grid_spherical
  write(i_unit_output,'(100(F16.8,X))')list_r(i),rho_r(1,i),mu_r(1,i),Energy_c_md_LDA_mu_of_r_grille_sphe_of_r(1,i)
  nbE += rho_r(1,i) *pas_grid_sphe 
  if(list_r(i)< r_core_sphe) then
    nbE_core += rho_r(1,i) *pas_grid_sphe 
  else
    nbE_valence += rho_r(1,i) *pas_grid_sphe
  endif 
 enddo
print*,'nbre delectrons =  ',nbE
print*,'nbre delectrons core=  ',nbE_core
print*,'nbre delectrons valence=  ',nbE_valence

write(*, '(A28,X,F16.10)') 'DFT mu(r) correlation new grid =',Energy_c_md_LDA_mu_of_r_grille_sphe_of_r_integrated
!write(*, '(A28,X,F16.10)') 'DFT mu(r) correlation old grid =',Energy_c_md_mu_of_r_LDA

!write(*, '(A28,X,F16.10)') 'DFT mu correlation core =',Energy_c_md_LDA_mu_of_r_integrated_core
!write(*, '(A28,X,F16.10)') 'DFT mu correlation valence=',Energy_c_md_LDA_mu_of_r_integrated_valence

write(*, '(A28,X,F16.10)') 'Mu moyen new    = ',rhomu_r_integrated
write(*, '(A28,X,F16.10)') 'mu_average old  = ',mu_average

end
