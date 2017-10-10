program CASSCF_bibi
  implicit none
  call check_core1_virt2_virt2
  call check_virt1_core2_core2 
  call check_core1_core2_virt2
  call check_core1_virt1_core2_virt2
  call check_core1_core1_virt2_virt2

!!! VERY SLOW 
! call check_semi_trans_virt_virt
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
      print*, integral,transformed_occ1_virt2_virt2(k,i,j),(transformed_occ1_virt2_virt2(k,i,j) - integral )
      pause
     endif
    enddo
   enddo
  enddo
 print*, 'passed check_core1_virt2_virt2'

end

subroutine check_semi_trans_virt_virt
 implicit none

 double precision :: get_mo_bielec_integral,integral
 real(integral_kind), allocatable :: bielec_tmp_0(:,:),matrix_tmp_1(:,:),matrix_final(:,:)
 integer :: i,j,k,l,p,q,m,n

 allocate(bielec_tmp_0(ao_num,ao_num))
 do i = 1, n_virt_orb
  print*, 'i',i
  do j = 1, n_virt_orb
   print*, 'j',j
    do p = 1,ao_num
      do q = 1,ao_num
      integral = 0.d0
      print*, i,j,p,q
       do m = 1, ao_num
        call get_ao_bielec_integrals(p,m,q,ao_num,bielec_tmp_0(1,m)) ! k,l :: r1, m :: r2
       enddo
       do m = 1, ao_num
        do n = 1, ao_num
         integral += mo_coef_virt(n,j) * mo_coef_virt(m,i) * bielec_tmp_0(n,m)
        enddo
       enddo
      if(dabs(integral - semi_transformed_virt_virt(j,i,q,p)).gt.1.d-10)then
       print*, j,i,q,p
       print*, integral, semi_transformed_virt_virt(j,i,q,p), dabs(integral - semi_transformed_virt_virt(j,i,q,p))
       stop 
      endif
      enddo
    enddo
    
  enddo
 enddo
 deallocate(bielec_tmp_0)

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
      print*, iorb,jorb,korb,lorb
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
