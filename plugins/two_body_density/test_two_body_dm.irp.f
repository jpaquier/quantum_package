program test_dm
 implicit none
 read_wf = .True.
 touch read_wf 
!call test_energy_dm
 call test_2dm
end

subroutine test_2dm
 implicit none
 provide two_bod_alpha_beta_mo


end

!subroutine test_energy_dm
!integer :: i,j,k,l
!double precision :: energy_2,get_mo_bielec_integral
!energy_2 = 0.d0
!do istate = 1, N_states
! do i = 1, mo_tot_num
!  do j = 1, mo_tot_num 
!   do k = 1, mo_tot_num
!    do l = 1, mo_tot_num
!     energy_2 += get_mo_bielec_integral(k,i,l,j,mo_integrals_map) * two_bod_alpha_beta_mo(k,l,j,i,istate)
!    enddo
!   enddo
!  enddo
! enddo
!enddo
!integer :: sze
!sze = mo_tot_num
!double precision, allocatable :: out_array(:,:)
!double precision :: integral
!allocate(out_array(sze,sze))
!do istate = 1, N_states
! do i = 1, mo_tot_num
!  do j = 1, mo_tot_num 
!   call get_mo_bielec_integrals_ij(i,j,sze,out_array,mo_map_integrals) 
!   do k = 1, mo_tot_num
!    do l = 1, mo_tot_num
!     integral = get_mo_bielec_integral(i,j,k,l,mo_integrals_map)
!     if(dabs(integral - out_array(k,l)).gt.1.d-10)then
!      print*,'i,j,k,l',i,j,k,l 
!      print*,dabs(integral - out_array(k,l)), integral, out_array(k,l)
!     endif
!    enddo
!   enddo
!  enddo
! enddo
!enddo
!
!print*,'psi_energy_bielec = ',psi_energy_bielec 
!print*,'test              = ',energy_2


!end
