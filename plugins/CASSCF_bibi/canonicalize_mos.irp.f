
subroutine diag_inactive_virt_and_update_mos_MR_Fock
 implicit none
 integer :: i,j,i_inact,j_inact,i_virt,j_virt
 double precision :: tmp(mo_tot_num,mo_tot_num)
 character*(64) :: label
 print*,'Diagonalizing the core, inact and virt Fock operator'
 tmp = 0.d0
 do i = 1, mo_tot_num
  do j = 1, mo_tot_num
   tmp(i,j) = MR_Fock_canonical_MO(i,j,1)
  enddo
 enddo
 
 do i = 1, n_act_orb
  integer :: i_act 
  i_act = list_act(i)
  do j = 1, n_virt_orb
   j_virt = list_virt(j)
   tmp(i_act,j_virt) = 0.d0
   tmp(j_virt,i_act) = 0.d0
  enddo
  do j = 1, n_core_inact_orb
   j_inact = list_core_inact(j)
   tmp(i_act,j_inact) = 0.d0
   tmp(j_inact,i_act) = 0.d0
  enddo
 enddo
 
 if(n_core_inact_orb.gt.0)then
  call diag_matrix_mo(tmp, mo_tot_num, list_core_inact, n_core_inact_orb, size(mo_coef,1),mo_coef)
 endif
 if(n_virt_orb.gt.0)then
  call diag_matrix_mo(tmp, mo_tot_num, list_virt, n_virt_orb, size(mo_coef,1),mo_coef)
 endif


 label = "Canonical"
!call mo_as_eigvectors_of_mo_matrix(tmp,size(tmp,1),size(tmp,2),label,1)
 touch mo_coef
end

