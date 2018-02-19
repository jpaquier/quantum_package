program test_dm
 implicit none
 read_wf = .True.
 touch read_wf 
 call test_energy_dm
end

subroutine test_energy_dm
 integer :: i,j,k,l
 double precision :: energy_2,get_mo_bielec_integral
 energy_2 = 0.d0
 do istate = 1, N_states
  do i = 1, mo_tot_num
   do j = 1, mo_tot_num 
    do k = 1, mo_tot_num
     do l = 1, mo_tot_num
      energy_2 += get_mo_bielec_integral(k,i,l,j,mo_integrals_map) * two_bod_alpha_beta(k,l,j,i,istate)
     enddo
    enddo
   enddo
  enddo
 enddo
 print*,'psi_energy_bielec = ',psi_energy_bielec 
 print*,'test              = ',energy_2


end
