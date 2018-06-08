program test
 read_wf = .True.
 touch read_wf
call test_subou
end



subroutine test_subou
 implicit none
double precision :: someee,nbE, nbE_a,nbE_b
integer :: j,i

 nbE=0.0
double precision :: Ecmd, mu_moy
 do i = 1, n_points_radial_grid_spherical
  write(33,'(100(F10.5,X))')list_r(i),rho_r(1,i),rhomu_r(1,i),Energy_c_md_LDA_mu_of_r_grille_sphe_of_r(1,i)
  nbE += rho_r(1,i) *pas_grid_sphe  
  Ecmd += Energy_c_md_LDA_mu_of_r_grille_sphe_of_r(1,i)*pas_grid_sphe
  mu_moy +=rhomu_r(1,i)*pas_grid_sphe
 enddo
print*,'nbre delectrons =  ',nbE

write(*, '(A28,X,F16.10)') 'DFT mu(r)correlation nouvelle grille=',Ecmd
write(*, '(A28,X,F16.10)') 'DFT mu(r) correlation ancienne grille=',Energy_c_md_mu_of_r_LDA

write(*, '(A28,X,F16.10)') 'Mu moyen nouveau= ',mu_moy 
write(*, '(A28,X,F16.10)') 'mu_average old=',mu_average

end
