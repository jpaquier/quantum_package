 BEGIN_PROVIDER [ integer, number_of_uniq_expo_per_l_shell_per_atom, (0:7,nucl_num)]
&BEGIN_PROVIDER [ integer, nmax_of_uniq_expo_per_l_shell_per_atom]
 implicit none
 BEGIN_DOC
! number_of_uniq_expo_per_l_shell_per_atom(i,j) = number of unique primitive exponents in the shell "i" of the nucleus "j"
 END_DOC
 integer :: i,j,l_type,index_ao,i_count,k,l,m,index_uniq
 integer :: iorder(ao_prim_num_max*ao_num),imax
 double precision :: expo_primitives(ao_prim_num_max*ao_num)
 double precision :: alpha
 imax = 0
 do i = 1, nucl_num
 !print*,''
 !print*,'**************************************************'
 !print*,'nucleus = ',i
  do l_type = 0,7  
 ! print*,"l_type = ",l_type
   i_count = 0
   number_of_uniq_expo_per_l_shell_per_atom(l_type,i) = 0
   if(Nucl_num_l_type_Aos(l_type,i) == 0)cycle
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
 ! print*,'number of expo = ',i_count
 ! print*,'*************'
   call dsort(expo_primitives,iorder,i_count) 
 ! print*,'*************'
   index_uniq = 1
   do l = 2, i_count
    if (expo_primitives(l) /= expo_primitives(l-1)) then
 !   print*,l,expo_primitives(l),iorder(l) 
     index_uniq += 1
    endif
   enddo
   if(imax .lt.index_uniq)then
    imax = index_uniq
   endif
 ! print*,'number of non-redondant expo = ',index_uniq
 ! print*,'*************'
 ! print*,''
   number_of_uniq_expo_per_l_shell_per_atom(l_type,i) = index_uniq
  enddo
 enddo
 nmax_of_uniq_expo_per_l_shell_per_atom = imax
 print*,'nmax_of_uniq_expo_per_l_shell_per_atom=',nmax_of_uniq_expo_per_l_shell_per_atom

END_PROVIDER


BEGIN_PROVIDER [ double precision, uniq_expo_per_l_shell_per_atom,(nmax_of_uniq_expo_per_l_shell_per_atom,0:7,nucl_num)]
 implicit none
 BEGIN_DOC
!Arrays containing the values of the exponents for the uncontracted gaussians
 END_DOC
integer :: i,j,l_type,index_ao,i_count,k,l,m,index_uniq
 integer :: iorder(ao_prim_num_max*ao_num)
 double precision :: expo_primitives(ao_prim_num_max*ao_num)
 double precision :: alpha
 do i = 1, nucl_num
  do l_type = 0,7
   i_count = 0
   if(Nucl_num_l_type_Aos(l_type,i) == 0)cycle
   do j = 1, Nucl_num_l_type_Aos(l_type,i)
    ! "physical" index of all the AO corresponding to the shell type l_type
    ! attached to atom i 
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
   call dsort(expo_primitives,iorder,i_count)
   index_uniq = 1
   uniq_expo_per_l_shell_per_atom(index_uniq,l_type,i) = expo_primitives(index_uniq)
   do l = 2, i_count
    if (expo_primitives(l) /= expo_primitives(l-1)) then
     index_uniq += 1
     uniq_expo_per_l_shell_per_atom(index_uniq,l_type,i) = expo_primitives(l)
    endif
   enddo
  enddo
 enddo



 
END_PROVIDER


