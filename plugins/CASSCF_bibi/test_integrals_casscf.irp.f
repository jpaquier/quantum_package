program CASSCF_bibi
  implicit none
  call check_core1_virt2_virt2
  call check_virt1_core2_core2 
  call check_core1_core2_virt2
  call check_core1_virt1_core2_virt2
  call check_core1_core1_virt2_virt2
end 

subroutine check_core1_virt2_virt2
 implicit none

 double precision :: get_mo_bielec_integral,integral
 integer :: i,j,k,l
 integer :: iorb,jorb,korb,lorb
 do i = 1, n_virt_orb
  iorb = list_virt(i)
  do j = 1, n_virt_orb
   jorb = list_virt(j)
   do k = 1, n_core_inact_orb
    korb = list_core_inact(k)
     integral = get_mo_bielec_integral(iorb,korb,jorb,korb,mo_integrals_map)
     if(dabs(transformed_occ1_virt2_virt2(k,i,j) - integral ).gt.1.d-10)then
      print*, 'pb in fully transformed'
      print*, iorb,jorb,korb
      print*, integral,transformed_occ1_virt2_virt2(k,i,j),dabs(transformed_occ1_virt2_virt2(k,i,j) - integral )
     endif
    enddo
   enddo
  enddo
 print*, 'passed check_core1_virt2_virt2'

end

subroutine check_virt1_core2_core2
 implicit none

 double precision :: get_mo_bielec_integral,integral
 integer :: i,j,k,l
 integer :: iorb,jorb,korb,lorb
 do i = 1, n_core_inact_orb
  iorb = list_core_inact(i)
  do j = 1, n_core_inact_orb
   jorb = list_core_inact(j)
   do k = 1, n_virt_orb
    korb = list_virt(k)
     integral = get_mo_bielec_integral(iorb,korb,jorb,korb,mo_integrals_map)
     if(dabs(transformed_virt1_occ2_occ2(k,i,j) - integral ).gt.1.d-10)then
      print*, 'pb in fully transformed'
      print*, iorb,jorb,korb
      print*, integral,transformed_virt1_occ2_occ2(k,i,j),dabs(transformed_virt1_occ2_occ2(k,i,j) - integral )
     endif
    enddo
   enddo
  enddo
 print*, 'passed check_virt1_core2_core2'

end

subroutine check_core1_core2_virt2
 implicit none

 double precision :: get_mo_bielec_integral,integral
 integer :: i,j,k,l
 integer :: iorb,jorb,korb,lorb
 do i = 1, n_core_inact_orb
  iorb = list_core_inact(i)
  do j = 1, n_virt_orb
   jorb = list_virt(j)
   do k = 1, n_core_inact_orb
    korb = list_core_inact(k)
     integral = get_mo_bielec_integral(iorb,korb,jorb,korb,mo_integrals_map)
     if(dabs(transformed_occ1_occ2_virt2(k,i,j) - integral ).gt.1.d-10)then
      print*, 'pb in fully transformed'
      print*, iorb,jorb,korb
      print*, integral,transformed_occ1_occ2_virt2(k,i,j),dabs(transformed_occ1_occ2_virt2(k,i,j) - integral )
     endif
    enddo
   enddo
  enddo
 print*, 'passed check_core1_core2_virt2'

end

subroutine check_core1_virt1_core2_virt2
 implicit none

 double precision :: get_mo_bielec_integral,integral
 integer :: i,j,k,l
 integer :: iorb,jorb,korb,lorb
 do i = 1, n_core_inact_orb
  iorb = list_core_inact(i)
  do j = 1, n_virt_orb
   jorb = list_virt(j)
   do k = 1, n_core_inact_orb
    korb = list_core_inact(k)
    do l = 1, n_virt_orb
     lorb = list_virt(l)
     integral = get_mo_bielec_integral(iorb,korb,jorb,lorb,mo_integrals_map)
     if(dabs(transformed_occ1_virt1_occ2_virt2(k,l,i,j) - integral ).gt.1.d-10)then
      print*, 'pb in fully transformed'
      print*, iorb,jorb,korb
      print*, integral,transformed_occ1_virt1_occ2_virt2(k,l,i,j),dabs(transformed_occ1_virt1_occ2_virt2(k,l,i,j) - integral )
     endif
    enddo
   enddo
  enddo
 enddo
 
 print*, 'passed check_core1_virt1_core2_virt2'

end

subroutine check_core1_core1_virt2_virt2
 implicit none

 double precision :: get_mo_bielec_integral,integral
 integer :: i,j,k,l
 integer :: iorb,jorb,korb,lorb
 do i = 1, n_core_inact_orb
  iorb = list_core_inact(i)
  do j = 1, n_core_inact_orb
   jorb = list_core_inact(j)
   do k = 1, n_virt_orb
    korb = list_virt(k)
    do l = 1, n_virt_orb
     lorb = list_virt(l)
     integral = get_mo_bielec_integral(iorb,korb,jorb,lorb,mo_integrals_map)
     if(dabs(transformed_virt1_virt1_occ2_occ2(k,l,i,j) - integral ).gt.1.d-10)then
      print*, 'pb in fully transformed'
      print*, iorb,jorb,korb
      print*, integral,transformed_virt1_virt1_occ2_occ2(k,l,i,j),dabs(transformed_virt1_virt1_occ2_occ2(k,l,i,j) - integral )
     endif
    enddo
   enddo
  enddo
 enddo
 
 print*, 'passed check_core1_core1_virt2_virt2'

end
