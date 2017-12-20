program print_bitmask
 implicit none
 read_wf = .True. 
 touch read_wf
 print*,'core'
 call debug_det(core_bitmask,N_int) 
 print*,'inact'
 call debug_det(inact_bitmask,N_int) 
 print*,'virt'
 call debug_det(virt_bitmask,N_int) 
 call routine

end

subroutine routine 
 implicit none
 integer :: i
 do i = 1, mo_tot_num
  write(*,'(1000(F16.10,X))')one_body_dm_mo_beta(i,:,1)
 enddo
end
