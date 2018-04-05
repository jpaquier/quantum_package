BEGIN_PROVIDER [ integer, number_of_uniq_expo_per_l_shell_per_atom, (0:7,nucl_num)]
 implicit none
 BEGIN_DOC
! number_of_uniq_expo_per_l_shell_per_atom(i,j) = number of unique primitive exponents in the shell "i" of the nucleus "j"
 END_DOC
 integer :: i,j,l_type,index_ao,i_count,k
 integer :: iorder(ao_prim_num_max*ao_num)
 double precision :: expo_primitives(ao_prim_num_max*ao_num)
 double precision :: alpha
 do i = 1, nucl_num
   print*,''
   print*,'*******************'
  print*,'nucleus = ',i
  do l_type = 0,7  
   print*,"l_type = ",l_type
   i_count = 0
   do j = 1, Nucl_num_l_type_Aos(l_type,i)
    ! "physical" index of all the AO corresponding to the shell type l_type attached to atom i 
    index_ao = Nucl_list_l_type_Aos(j,l_type,i) 
    ! you browse all the primitives of the AO index_ao
    do k = 1, ao_prim_num(index_ao)
     ! exponent of the kth primitive of the index_aoth AO
     i_count += 1
     alpha = ao_expo_ordered_transp(k,index_ao)
     expo_primitives(i_count) = -alpha
     iorder(i_count) = i_count
    enddo
   enddo 
   print*,'number of expo = ',i_count
   call dsort(expo_primitives,iorder,i_count) 
   print*,'*************'
   do k = 1, i_count
    print*,k,expo_primitives(k),iorder(k)
   enddo
   print*,'*************'
   print*,''
  enddo
 enddo

END_PROVIDER 
