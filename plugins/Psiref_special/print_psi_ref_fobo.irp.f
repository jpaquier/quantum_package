program pouet
 implicit none
!read_wf = .True.
!touch read_wf
 use bitmasks
 integer(bit_kind), allocatable :: key_tmp(:,:)
 allocate(key_tmp(N_int,2))
 integer :: i,j,k,l
 double precision :: h00,hii
 do k = 1, N_int
  key_tmp(k,1) = ref_bitmask(k,1)
  key_tmp(k,2) = ref_bitmask(k,2)
 enddo
 integer :: sigma,sigma_twisted,pi_star_1,pi_star_2
 integer :: i_ok
 sigma = 71 
 sigma_twisted = 68 
 pi_star_1 = 72 
 pi_star_2 = 73
 print*, 'J_sigma/pi*_1        ',mo_bielec_integral_jj(sigma,pi_star_1)
 print*, 'J_sigma/pi*_2        ',mo_bielec_integral_jj(sigma,pi_star_2)
 print*, 'J_sigma_twisted/pi*_1',mo_bielec_integral_jj(sigma_twisted,pi_star_1)
 print*, 'J_sigma_twisted/pi*_2',mo_bielec_integral_jj(sigma_twisted,pi_star_2)
 print*, 'ref bitmask'
 call debug_det(ref_bitmask,N_int)
 call do_mono_excitation(key_tmp,69,pi_star_1,1,i_ok)
 print*, ' LMCT      '
 call debug_det(key_tmp,N_int)
 call i_H_j(ref_bitmask,ref_bitmask,N_int,h00) 
 call i_H_j(key_tmp,key_tmp,N_int,hii) 
 print*, 'h00 - hii',h00-hii
! integer :: i
!print*, 'n_list_lmct = ',n_list_lmct
!do i = 1, n_list_lmct
! print*, list_lmct(i)
!enddo
!print*, ''
!print*, 'n_list_mlct = ',n_list_mlct
!do i = 1, n_list_mlct
! print*, list_mlct(i)
!enddo
!call routine

end

subroutine routine
 implicit none
 integer :: i
 print*, 'N_det_ref_fobo     = ',N_det_ref_fobo
 print*, 'N_det_non_ref_fobo = ',N_det_non_ref_fobo
 print*, 'N_det_lmct         = ',N_det_lmct
 print*, 'N_det_mlct         = ',N_det_mlct
 print*, '-----'
 print*, '-----'
 print*, '     '
 print*, 'psi_det fobo'
 do i = 1, N_det_ref_fobo
 print*, '-----'
  print*, 'i',i
  print*, 'coef = ',psi_ref_fobo_coef(i,:)
  call debug_det(psi_ref_fobo(1,1,i),N_int)
 print*, '-----'
 enddo

end
