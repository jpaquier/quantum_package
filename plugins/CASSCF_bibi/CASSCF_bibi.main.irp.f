program CASSCF_bibi
  implicit none
  call routine1


end
subroutine routine1
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

 integer :: n,m
 double precision :: accu,ao_bielec_integral,get_ao_bielec_integral
 double precision :: integral_map,integral_ao
 double precision :: matrix_tmp_1(n_core_inact_orb,ao_num),matrix_final(n_core_inact_orb,n_virt_orb)
 double precision :: vec_1(ao_num),vec_2(ao_num),vec_3(ao_num),accu_2,bielec_tmp_0(ao_num,ao_num)
 double precision :: vec_4(ao_num),accu_3
 do k = 1, ao_num
  do l = 1, ao_num
   matrix_tmp_1 = 0.d0
   matrix_final = 0.d0
   do m = 1, ao_num
    call get_ao_bielec_integrals(k,m,l,ao_num,bielec_tmp_0(1,m)) ! all integrals for a given l1, k1
   enddo
   call dgemm('N','N',n_core_inact_orb,ao_num,ao_num,1.d0,mo_coef_core_inact_transp,n_core_inact_orb,bielec_tmp_0,ao_num,0.d0,matrix_tmp_1,n_core_inact_orb)
   call dgemm('N','N',n_core_inact_orb,n_virt_orb,ao_num,1.d0,matrix_tmp_1,n_core_inact_orb,mo_coef_virt,ao_num,0.d0,matrix_final,n_core_inact_orb)
   do i = 1, n_core_inact_orb
    iorb = list_core_inact(i)
    do j = 1, n_virt_orb
     jorb = list_virt(j)
      accu_2 = 0.d0
      do m = 1, ao_num
       call get_ao_bielec_integrals(k,m,l,ao_num,bielec_tmp_0(1,m)) ! all integrals for a given l1, k1
      enddo
      do m = 1, ao_num
       do n = 1, ao_num
        accu_2 += bielec_tmp_0(m,n) * mo_coef(n,jorb) * mo_coef(m,iorb)
       enddo
      enddo
      if(dabs(accu_2-semi_transformed_occ_virt(k,l,i,j)).gt.1.d-10)then
       print*, 'pb !!'
       print*, k,l,i,j
       print*, accu_2,semi_transformed_occ_virt(k,l,i,j)
       stop
      endif
    enddo
   enddo
  enddo
 enddo


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


 do i = 1, n_core_inact_orb
  iorb = list_core_inact(i)
  do j = 1, n_virt_orb
   jorb = list_virt(j)
   do k = 1, n_core_inact_orb
    korb = list_core_inact(k)
    do l = 1, n_virt_orb
     lorb = list_virt(l)
     integral = get_mo_bielec_integral(iorb,korb,jorb,lorb,mo_integrals_map)
     if(dabs(transformed_occ_virt(k,l,i,j) - integral ).gt.1.d-10)then
      print*, 'pb in fully transformed'
      print*, iorb,jorb,korb,lorb
      print*, integral,transformed_occ_virt(k,l,i,j),dabs(transformed_occ_virt(k,l,i,j) - integral )
     endif
    enddo
   enddo
  enddo
 enddo

end

subroutine routine2
  integer :: i,iorb,j,jorb,k,korb,l,lorb
  double precision, allocatable :: matrix_integrals(:,:)
  double precision :: get_mo_bielec_integral,integral
  allocate(matrix_integrals(n_core_inact_orb, n_virt_orb))
  iorb = 1
  jorb = 6 
  i = list_core_inact_reverse(iorb)
  j = list_virt_reverse(jorb)
  print*, '------'
  print*, i,j
  print*, '---------'
  call get_all_core_inact_virt_integrals(i,j,matrix_integrals)
end
