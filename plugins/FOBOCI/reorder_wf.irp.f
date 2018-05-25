subroutine reorder_wf
  use bitmasks
 implicit none
 integer(bit_kind), allocatable :: psi_det_tmp(:,:,:)
 double precision, allocatable  :: psi_coef_tmp(:,:)
 integer :: i,j,k,l
 integer :: index_ref_generators_restart(N_states), degree
 logical, allocatable :: is_a_ref_det(:) 

 allocate(psi_det_tmp(N_int,2,N_det), is_a_ref_det(N_det), psi_coef_tmp(N_det, N_states))
 do i = 1, N_det
  is_a_ref_det(i) = .False.
  do k = 1, N_int
   psi_det_tmp(k,1,i) = psi_det(k,1,i)
   psi_det_tmp(k,2,i) = psi_det(k,2,i)
  enddo
  do k = 1, N_states
   psi_coef_tmp(i,k) = psi_coef(i,k)
  enddo
   call get_excitation_degree(ref_generators_restart(1,1,1),psi_det(1,1,i),degree,N_int)   
   if(degree == 0)then
    index_ref_generators_restart = i
   endif
  do j = 1, N_det_generators_restart
   call get_excitation_degree(psi_det(1,1,i),psi_det_generators_restart(1,1,j),degree,N_int)  
   if(degree == 0)then
    is_a_ref_det(i) = .True.
    exit
   endif
  enddo
 enddo
 
 j = 0
 do i = 1, N_det
  if(is_a_ref_det(i))then
   j+=1 
  endif
 enddo
 if(j.ne.N_det_generators_restart)then
  print*, 'pb in reorder_wf !!'
  print*, 'pv in is_a_ref_det'
  print*, 'j .ne. N_det_generators_restart'
  print*, j,N_det_generators_restart
 
 
 endif

 ! ref determinant in first position
 do k = 1, N_int
  psi_det(k,1,1) = psi_det_tmp(k,1,index_ref_generators_restart(1)) 
  psi_det(k,2,1) = psi_det_tmp(k,2,index_ref_generators_restart(1)) 
 enddo
 do k = 1, N_states
  psi_coef(1,k) = psi_coef_tmp(index_ref_generators_restart(1),k)
 enddo
 

 ! other CAS determinants after
 j = 1
 do i = 1, N_det
  if(.not.is_a_ref_det(i))cycle
  if(i==index_ref_generators_restart(1))cycle
  j+=1 
  do k = 1, N_int
   psi_det(k,1,j) = psi_det_tmp(k,1,i) 
   psi_det(k,2,j) = psi_det_tmp(k,2,i) 
  enddo
  do k = 1, N_states
   psi_coef(j,k) = psi_coef_tmp(i,k)
  enddo
 enddo
 if(j.ne.N_det_generators_restart)then
  print*, 'pb in reorder_wf !!'
  print*, 'j .ne. N_det_generators_restart'
  print*, j,N_det_generators_restart
  do i = 1, j
   print*, ''
   print*, i
   write(*, '(100(I16,X))') psi_det(:,:,i)
   print*, ''
  enddo
  stop
 endif

 j = N_det_generators_restart
 do i = 1, N_det
  if(is_a_ref_det(i))cycle
  j+=1 
  do k = 1, N_int
   psi_det(k,1,j) = psi_det_tmp(k,1,i) 
   psi_det(k,2,j) = psi_det_tmp(k,2,i) 
  enddo
  do k = 1, N_states
   psi_coef(j,k) = psi_coef_tmp(i,k)
  enddo
 enddo
 if(j.ne.N_det)then
  print*, 'pb in reorder_wf !!'
  print*, 'j .ne. N_det'
  print*, j,N_det
  stop
 endif
 
 deallocate(psi_det_tmp,is_a_ref_det,psi_coef_tmp)
 SOFT_TOUCH N_det psi_det psi_coef

end
