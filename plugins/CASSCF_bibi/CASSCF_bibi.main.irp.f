program CASSCF_bibi
  implicit none
  integer :: i,iorb,j,jorb,k,korb,l,lorb
! do i = 1, n_core_inact_orb
!  iorb = list_core_inact(i)
!  do j = 1, n_core_inact_orb
!   jorb = list_core_inact(j)
!    write(*,*)iorb,jorb, Fock_matrix_mo(iorb,jorb)
!  enddo
! enddo
!print*, '*****************'
!print*, '*****************'
!print*, '*****************'
!print*, '*****************'
!print*, '*****************'
!print*, '*****************'
!do i = 1, ao_num
! do j = 1, ao_num
!  write(*,*)Fock_matrix_ao(i,j)
! enddo
!enddo
!call iteration_scf
 double precision, allocatable :: matrix_integrals(:,:)
 double precision :: get_mo_bielec_integral,integral
 allocate(matrix_integrals(n_core_inact_orb, n_virt_orb))
 do i = 1, n_core_inact_orb
  iorb = list_core_inact(i)
  do j = 1, n_virt_orb
   jorb = list_virt(j)
   call get_all_core_inact_virt_integrals(i,j,matrix_integrals)
   do k = 1, n_core_inact_orb
    korb = list_core_inact(k)
    do l = 1, n_virt_orb
     lorb = list_virt(l)
     integral = get_mo_bielec_integral(iorb,korb,jorb,lorb,mo_integrals_map)
     if(dabs(matrix_integrals(k,l) - integral ).gt.1.d-10)then
      print*, iorb,jorb,korb,lorb
      print*, integral,matrix_integrals(k,l),dabs(matrix_integrals(k,l) - integral )
     endif
    enddo
   enddo
  enddo
 enddo
end
