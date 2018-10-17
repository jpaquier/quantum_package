program test_dm
 implicit none
 read_wf = .True.
 touch read_wf 
 call routine
end

subroutine routine
 implicit none
 integer :: i,j,k,l,i_state
 i_state = 1
 do i = 1, mo_tot_num
  do j = 1, mo_tot_num
   do k = 1, mo_tot_num
    do l = 1, mo_tot_num
     if(dabs(two_bod_alpha_beta_mo_contracted(l,k,j,i,i_state) - two_bod_alpha_beta_mo_contracted_parallel(l,k,j,i,i_state)).gt.1.d-10)then
      print*,i,j,k,l
      print*,dabs(two_bod_alpha_beta_mo_contracted(l,k,j,i,i_state) - two_bod_alpha_beta_mo_contracted_parallel(l,k,j,i,i_state)),two_bod_alpha_beta_mo_contracted(l,k,j,i,i_state) , two_bod_alpha_beta_mo_contracted_parallel(l,k,j,i,i_state)
     endif
    enddo
   enddo
  enddo
 enddo


end
