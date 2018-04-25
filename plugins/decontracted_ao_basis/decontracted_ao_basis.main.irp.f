program decontracted_ao_basis
  implicit none
  BEGIN_DOC
! TODO
  END_DOC
  integer :: i,j,l,l_type,index_uniq
 
!Print each aos on each nuclei, and labem them.
 
 !do i = 1, nucl_num
 ! print*,'nucleus number ',i
 ! do l_type = 0,7
 !  print*,'l_type = ', l_type
 !  print*,'Nucl_num_l_type_Aos    ',Nucl_num_l_type_Aos(l_type,i)
 !  do j = 1, Nucl_num_l_type_Aos(l_type,i)
 !   print*,'Nucl_list_l_type_Aos = ',Nucl_list_l_type_Aos(j,l_type,i) 
 !  enddo
 ! enddo 
 !enddo

 do i = 1, nucl_num
  print*,''
  print*,'**************************************************'
  print*,'nucleus = ',i
  do l_type = 0,7
   print*,"l_type = ",l_type
   print*,'number_of_uniq_expo_per_l_shell_per_atom = ', number_of_uniq_expo_per_l_shell_per_atom(l_type,i)
   do index_uniq = 1, number_of_uniq_expo_per_l_shell_per_atom(l_type,i)
    print*, uniq_expo_per_l_shell_per_atom(index_uniq,l_type,i)
   enddo
  enddo
 enddo

end
