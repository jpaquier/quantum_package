program save_fock_inactiv_virt_mos
 implicit none
 read_wf = .True.
 touch read_wf
 call routine
!call print_fock

end


subroutine print_fock
 implicit none
 integer :: i,j
!print*, ''
!print*, 'SR FOCK MATRIX'
!print*, ''
!do i = 1, n_core_inact_orb
! j = list_core_inact(i)
! write(*,'(500(F16.10,X))') Fock_matrix_mo(j,:)
! write(*,'(500(F16.10,X))') Fock_matrix_alpha_mo(j,:)
! write(*,'(500(F16.10,X))') Fock_matrix_beta_mo(j,:)
!enddo
 print*, ''
 print*, ''
 print*, 'MR FOCK MATRIX'
 print*, ''
 do i = 1, n_core_inact_orb
  j = list_core_inact(i)
  write(*,'(500(F16.10,X))') MR_Fock_canonical_MO(j,:,1)
! write(*,'(500(F16.10,X))') MR_Fock_matrix_alpha_mo(j,:,1)
! write(*,'(500(F16.10,X))') MR_Fock_matrix_beta_mo(j,:,1)
 enddo


end

subroutine routine
implicit none
!if(N_det.gt.1)then

  print*, 'Using MR FOCK matrix'
! call diag_inactive_virt_and_update_mos_MR_Fock
!else
! print*, 'Using SR FOCK matrix'
  call diag_inactive_virt_and_update_mos_SR_Fock
!endif
!call save_mos


end
