BEGIN_PROVIDER [ double precision, mo_coef_begin_iteration, (ao_num,mo_tot_num) ]
   implicit none
   BEGIN_DOC
   ! Alpha and beta one-body density matrix that will be used for the 1h1p approach
   END_DOC
END_PROVIDER

subroutine initialize_mo_coef_begin_iteration
 implicit none
 mo_coef_begin_iteration = mo_coef
end

subroutine reorder_active_orb
 implicit none
 integer :: i,j,iorb
 integer :: k,l
 double precision, allocatable :: accu(:)
 integer, allocatable :: index_active_orb(:),iorder(:)
 double precision, allocatable :: mo_coef_tmp(:,:)
 allocate(accu(mo_tot_num),index_active_orb(n_act_orb),iorder(mo_tot_num))
 allocate(mo_coef_tmp(ao_num,mo_tot_num))
 
 do i = 1, n_act_orb
  iorb = list_act(i)
  do j = 1, mo_tot_num
   accu(j) = 0.d0
   iorder(j) = j
   do k = 1, ao_num
    do l = 1, ao_num
     accu(j) += mo_coef_begin_iteration(k,iorb) * mo_coef(l,j) * ao_overlap(k,l)
    enddo
   enddo
   accu(j) = -dabs(accu(j))
  enddo
  call dsort(accu,iorder,mo_tot_num)
  index_active_orb(i) = iorder(1) 
! print*, i,iorder(1),accu(1)
 enddo

 double precision :: x
 integer :: i1,i2
 print*, 'swapping the active MOs'
 do j = 1, n_act_orb
  i1 = list_act(j)
  i2 = index_active_orb(j)
  print*, i1,i2
  do i=1,ao_num
    x = mo_coef(i,i1)
    mo_coef(i,i1) = mo_coef(i,i2)
    mo_coef(i,i2) = x
  enddo
 enddo
 call loc_cele_routine

 deallocate(accu,index_active_orb, iorder)
end


subroutine reorder_all_orb
 implicit none
 integer :: i,j,iorb
 integer :: k,l
 integer :: index_orb(mo_tot_num)
 double precision :: mo_coef_tmp(ao_num,mo_tot_num)
 integer :: i1,i2
 
 mo_coef_tmp = mo_coef
 
 call reorder_all_mos(mo_coef_begin_iteration,index_orb)

 print*, 'swapping the MOs'
 do i1 = 1, mo_tot_num
  i2 = index_orb(i1)
  print*, i1,i2
  do i=1,ao_num
    mo_coef(i,i1) = mo_coef_tmp(i,i2)
    mo_coef(i,i2) = mo_coef_tmp(i,i1)
  enddo
 enddo

!call loc_cele_routine
!touch mo_coef

end

subroutine reorder_set_of_mos(mo_coef_before,list_orb,n_orb,index_orb)
 implicit none
 integer, intent(in) :: n_orb
 integer, intent(in) :: list_orb(n_orb)
 double precision, intent(in) :: mo_coef_before(ao_num, mo_tot_num)
 integer, intent(out) :: index_orb(mo_tot_num)

 double precision :: mo_coef_tmp
 double precision :: accu(mo_tot_num)
 logical          :: is_chosen(mo_tot_num)
 integer          :: iorder(mo_tot_num)
 integer          :: i,j,k,l
 
 integer :: iorb
 
 accu = 0.d0
 is_chosen = .False.
 do i = 1, n_orb
  iorb = list_orb(i)
  accu = 0.d0
  if(.not.is_chosen(iorb))then
  do j = 1, mo_tot_num
   accu(j) = 0.d0
   iorder(j) = j
   if(.not.is_chosen(j))then
    do k = 1, ao_num
     do l = 1, ao_num
      accu(j) += mo_coef_before(k,iorb) * mo_coef(l,j) * ao_overlap(k,l)
     enddo
    enddo
   endif
   accu(j) = -dabs(accu(j))
  enddo
   call dsort(accu,iorder,mo_tot_num)
   is_chosen(iorder(1)) = .True.
   is_chosen(iorb) = .True.
   index_orb(iorb) = iorder(1) 
   index_orb(iorder(1)) = iorb
  endif
 enddo


end

subroutine reorder_all_mos(mo_coef_before,index_orb)
 implicit none
 double precision, intent(in) :: mo_coef_before(ao_num, mo_tot_num)
 integer, intent(out) :: index_orb(mo_tot_num)

 double precision :: mo_coef_tmp
 double precision :: accu(mo_tot_num)
 logical          :: is_chosen(mo_tot_num)
 integer          :: iorder(mo_tot_num)
 integer          :: i,j,k,l
 
 integer :: iorb
 
 accu = 0.d0
 is_chosen = .False.


 do i = 1, n_core_orb
  iorb = list_core(i)
  accu = 0.d0
  if(.not.is_chosen(iorb))then
  do j = 1, mo_tot_num
   accu(j) = 0.d0
   iorder(j) = j
   if(.not.is_chosen(j))then
    do k = 1, ao_num
     do l = 1, ao_num
      accu(j) += mo_coef_before(k,iorb) * mo_coef(l,j) * ao_overlap(k,l)
     enddo
    enddo
   endif
   accu(j) = -dabs(accu(j))
  enddo
   call dsort(accu,iorder,mo_tot_num)
   is_chosen(iorder(1)) = .True.
   is_chosen(iorb) = .True.
   index_orb(iorb) = iorder(1) 
   index_orb(iorder(1)) = iorb
  endif
 enddo


 do i = 1, n_inact_orb
  iorb = list_inact(i)
  accu = 0.d0
  if(.not.is_chosen(iorb))then
  do j = 1, mo_tot_num
   accu(j) = 0.d0
   iorder(j) = j
   if(.not.is_chosen(j))then
    do k = 1, ao_num
     do l = 1, ao_num
      accu(j) += mo_coef_before(k,iorb) * mo_coef(l,j) * ao_overlap(k,l)
     enddo
    enddo
   endif
   accu(j) = -dabs(accu(j))
  enddo
   call dsort(accu,iorder,mo_tot_num)
   is_chosen(iorder(1)) = .True.
   is_chosen(iorb) = .True.
   index_orb(iorb) = iorder(1) 
   index_orb(iorder(1)) = iorb
  endif
 enddo

 do i = 1, n_act_orb
  iorb = list_act(i)
  accu = 0.d0
  if(.not.is_chosen(iorb))then
  do j = 1, mo_tot_num
   accu(j) = 0.d0
   iorder(j) = j
   if(.not.is_chosen(j))then
    do k = 1, ao_num
     do l = 1, ao_num
      accu(j) += mo_coef_before(k,iorb) * mo_coef(l,j) * ao_overlap(k,l)
     enddo
    enddo
   endif
   accu(j) = -dabs(accu(j))
  enddo
   call dsort(accu,iorder,mo_tot_num)
   is_chosen(iorder(1)) = .True.
   is_chosen(iorb) = .True.
   index_orb(iorb) = iorder(1) 
   index_orb(iorder(1)) = iorb
  endif
 enddo

 do i = 1, n_virt_orb
  iorb = list_virt(i)
  accu = 0.d0
  if(.not.is_chosen(iorb))then
  do j = 1, mo_tot_num
   accu(j) = 0.d0
   iorder(j) = j
   if(.not.is_chosen(j))then
    do k = 1, ao_num
     do l = 1, ao_num
      accu(j) += mo_coef_before(k,iorb) * mo_coef(l,j) * ao_overlap(k,l)
     enddo
    enddo
   endif
   accu(j) = -dabs(accu(j))
  enddo
   call dsort(accu,iorder,mo_tot_num)
   is_chosen(iorder(1)) = .True.
   is_chosen(iorb) = .True.
   index_orb(iorb) = iorder(1) 
   index_orb(iorder(1)) = iorb
  endif
 enddo



end
