 BEGIN_PROVIDER [integer, n_shell_in_basis, (0:ao_l_max)]
&BEGIN_PROVIDER [integer, n_shell_max_in_basis]
 implicit none
 integer :: i
 n_shell_in_basis = 0
 do i = 1, ao_num
  n_shell_in_basis(ao_l(i)) += 1
 enddo
 n_shell_max_in_basis = maxval(n_shell_in_basis) 
END_PROVIDER 


 BEGIN_PROVIDER [integer, n_spherical_shell_in_basis, (0:ao_l_max)]
 implicit none
 integer :: i
 !n_spherical_shell_in_basis(i) = number of spherical AOS of order "i" in the basis
 ! 6  "D-like" GTOs turns into 5  spherical harmonics, which are accounted as 1 shell
 ! 10 "F-like" GTOs turns into 7  spherical harmonics, which are accounted as 1 shell
 ! 15 "G-like" GTOs turns into 9  spherical harmonics, which are accounted as 1 shell
 ! 21 "H-like" GTOs turns into 11 spherical harmonics, which are accounted as 1 shell
 ! 28 "I-like" GTOs turns into 13 spherical harmonics, which are accounted as 1 shell
 ! 36 "J-like" GTOs turns into 15 spherical harmonics, which are accounted as 1 shell
 END_DOC
 n_spherical_shell_in_basis(2) = n_shell_in_basis(2)/6
 n_spherical_shell_in_basis(3) = n_shell_in_basis(3)/10
 n_spherical_shell_in_basis(4) = n_shell_in_basis(4)/15
 n_spherical_shell_in_basis(5) = n_shell_in_basis(5)/21
 n_spherical_shell_in_basis(6) = n_shell_in_basis(6)/28
 n_spherical_shell_in_basis(7) = n_shell_in_basis(7)/36
END_PROVIDER 

BEGIN_PROVIDER [integer, n_spherical_AOs_in_basis]
 implicit none
 BEGIN_DOC 
 ! total number of spherical harmonics present in the basis
 END_DOC
 integer :: i,num_shell
 n_spherical_AOs_in_basis = 0
 do i = 0, ao_l_max
  if(i==0)then
   num_shell = n_shell_in_basis(i)
  else if(i==1)then
   num_shell = n_shell_in_basis(i)
  else if(i==2)then
   num_shell = n_spherical_shell_in_basis(i) * 5
  else if(i==3)then
   num_shell = n_spherical_shell_in_basis(i) * 7
  else if(i==4)then
   num_shell = n_spherical_shell_in_basis(i) * 9
  else if(i==5)then
   num_shell = n_spherical_shell_in_basis(i) * 11
  else if(i==6)then
   num_shell = n_spherical_shell_in_basis(i) * 13
  else if(i==7)then
   num_shell = n_spherical_shell_in_basis(i) * 15
  endif
  n_spherical_AOs_in_basis += num_shell
 enddo
END_PROVIDER 


!BEGIN_PROVIDER [double precision, cartesian_to_spherical_matrix, (ao_num,n_spherical_AOs_in_basis)]
!BEGIN_PROVIDER [integer, list_shell_in_basis, (n_shell_max_in_basis,0:ao_l_max)]
!BEGIN_PROVIDER [integer, ao_spherical_l, (n_spherical_AOs_in_basis)]
!BEGIN_PROVIDER [integer, ao_spherical_nucl, (n_spherical_AOs_in_basis)]
!implicit none
!BEGIN_DOC 
!! carthesian to spherical transformation matrix containing overlap !!!
!! cartesian_to_spherical_matrix(j,i) = <AO_j_cartesian|AO_i_spherical>
!! |AO_i_spherical> = \sum_{j = 1, ao_num} cartesian_to_spherical_matrix(j,i) |AO_i_spherical>
!! list_shell_in_basis = list of all spherical AOs by shell
!! list_shell_in_basis(j,i) = "jth" spherical AO of angular momentum "i" 
!END_DOC
!integer :: i,k,l,i_start,i_start_0,n_shell_passed,index_ao,n_shell_passed_0
!cartesian_to_spherical_matrix = 0.d0
!i_start = 1
!n_shell_passed = 1 
!do while (i_start .le. ao_num)
! if      (ao_l(i_start)==0)then
!  cartesian_to_spherical_matrix(i_start, n_shell_passed) = 1.d0
!  i_start += 1
!  ao_spherical_l(0) += 1
!  index_ao = ao_spherical_l(0) 
!  list_shell_in_basis(index_ao,0) = n_shell_passed
!  ao_spherical_nucl(n_shell_passed) = ao_nucl(i_start)
!  n_shell_passed += 1
! else if (ao_l(i_start)==1)then
!  i_start_0 = i_start
!  n_shell_passed_0 = n_shell_passed
!  do k = 1, 3 
!   cartesian_to_spherical_matrix(i_start, n_shell_passed) = 1.d0
!   i_start += 1
!   do i = 1, 3
!    do l = i,i 
!     cartesian_to_spherical_matrix(i_start_0+i-1,n_shell_passed-1+k) += cart_to_sphe_1(l,k) !* ao_overlap(i+i_start_0-1,l+i_start_0-1)
!    enddo
!   enddo
!   ao_spherical_l(1) += 1
!   index_ao = ao_spherical_l(1) 
!   list_shell_in_basis(index_ao,1) = n_shell_passed_0
!   n_shell_passed_0 += 1
!  enddo
!  n_shell_passed += 3
! else if (ao_l(i_start)==2)then
!  i_start_0 = i_start
!  n_shell_passed_0 = n_shell_passed
!  do k = 1, 5
!   do i = 1, 6
!    do l = 1, 6
!     cartesian_to_spherical_matrix(i_start_0+i-1,n_shell_passed-1+k) += cart_to_sphe_2(l,k) * ao_overlap(i+i_start_0-1,l+i_start_0-1)
!    enddo
!   enddo
!   ao_spherical_l(2) += 1
!   index_ao = ao_spherical_l(1) 
!   list_shell_in_basis(index_ao,2) = n_shell_passed_0
!   n_shell_passed_0 += 1
!  enddo
!  n_shell_passed += 5
!  i_start += 6
! else if (ao_l(i_start)==3)then
!  i_start_0 = i_start
!  n_shell_passed_0 = n_shell_passed
!  do k = 1, 7
!   do i = 1, 10
!    do l = 1, 10
!     cartesian_to_spherical_matrix(i_start_0+i-1,n_shell_passed-1+k) += cart_to_sphe_3(l,k) * ao_overlap(i+i_start_0-1,l+i_start_0-1)
!    enddo
!   enddo
!   ao_spherical_l(2) += 1
!   index_ao = ao_spherical_l(3) 
!   list_shell_in_basis(index_ao,3) = n_shell_passed_0
!   n_shell_passed_0 += 1
!  enddo
!  n_shell_passed += 7
!  i_start += 10
! else if (ao_l(i_start)==4)then
!  i_start_0 = i_start
!  n_shell_passed_0 = n_shell_passed
!  do k = 1, 9
!   do i = 1, 15
!    do l = 1, 15
!     cartesian_to_spherical_matrix(i_start_0+i-1,n_shell_passed-1+k) += cart_to_sphe_4(l,k) * ao_overlap(i+i_start_0-1,l+i_start_0-1)
!    enddo
!   enddo
!   ao_spherical_l(4) += 1
!   index_ao = ao_spherical_l(1) 
!   list_shell_in_basis(index_ao,4) = n_shell_passed_0
!   n_shell_passed_0 += 1
!  enddo
!  n_shell_passed += 9
!  i_start += 15
! else if (ao_l(i_start)==5)then
!  i_start_0 = i_start
!  n_shell_passed_0 = n_shell_passed
!  do k = 1, 11
!   do i = 1, 21
!    do l = 1, 21
!     cartesian_to_spherical_matrix(i_start_0+i-1,n_shell_passed-1+k) += cart_to_sphe_5(l,k) * ao_overlap(i+i_start_0-1,l+i_start_0-1)
!    enddo
!   enddo
!   ao_spherical_l(5) += 1
!   index_ao = ao_spherical_l(1) 
!   list_shell_in_basis(index_ao,5) = n_shell_passed_0
!   n_shell_passed_0 += 1
!  enddo
!  n_shell_passed += 11
!  i_start += 21
! else if (ao_l(i_start)==6)then
!  i_start_0 = i_start
!  n_shell_passed_0 = n_shell_passed
!  do k = 1, 13
!   do i = 1, 28
!    do l = 1, 28
!     cartesian_to_spherical_matrix(i_start_0+i-1,n_shell_passed-1+k) += cart_to_sphe_6(l,k) * ao_overlap(i+i_start_0-1,l+i_start_0-1)
!    enddo
!   enddo
!   ao_spherical_l(6) += 1
!   index_ao = ao_spherical_l(1) 
!   list_shell_in_basis(index_ao,6) = n_shell_passed
!   n_shell_passed_0 += 1
!  enddo
!  n_shell_passed += 13
!  i_start += 28
! else if (ao_l(i_start)==7)then
!  i_start_0 = i_start
!  n_shell_passed_0 = n_shell_passed
!  do k = 1, 15
!   do i = 1, 36
!    do l = 1, 36
!     cartesian_to_spherical_matrix(i_start_0+i-1,n_shell_passed-1+k) += cart_to_sphe_7(l,k) * ao_overlap(i+i_start_0-1,l+i_start_0-1)
!    enddo
!   enddo
!   ao_spherical_l(7) += 1
!   index_ao = ao_spherical_l(1) 
!   list_shell_in_basis(index_ao,7) = n_shell_passed_0
!   n_shell_passed_0 += 1
!  enddo
!  n_shell_passed += 15
!  i_start += 36
! endif
!enddo

!ND_PROVIDER


 BEGIN_PROVIDER [double precision, cartesian_to_spherical_matrix, (ao_num,n_spherical_AOs_in_basis)]
&BEGIN_PROVIDER [integer, list_shell_in_basis, (n_shell_max_in_basis,0:ao_l_max)]
&BEGIN_PROVIDER [integer, list_shell_in_basis_per_atom, (n_shell_max_in_basis,0:ao_l_max,nucl_num)]
&BEGIN_PROVIDER [integer, n_shell_in_basis_per_atom, (0:ao_l_max,nucl_num)]
&BEGIN_PROVIDER [integer, ao_spherical_l, (n_spherical_AOs_in_basis)]
&BEGIN_PROVIDER [integer, ao_spherical_nucl, (n_spherical_AOs_in_basis)]
 implicit none
 BEGIN_DOC 
 ! carthesian to spherical transformation matrix containing overlap !!!
 ! cartesian_to_spherical_matrix(j,i) = <AO_j_cartesian|AO_i_spherical>
 ! |AO_i_spherical> = \sum_{j = 1, ao_num} cartesian_to_spherical_matrix(j,i) |AO_i_spherical>
 ! list_shell_in_basis = list of all spherical AOs by shell
 ! list_shell_in_basis(j,i) = "jth" spherical AO of angular momentum "i" 
 END_DOC
 integer :: i,k,l,i_start,i_start_0,n_shell_passed,index_ao,n_shell_passed_0
 integer :: i_nucl,n_shell
 cartesian_to_spherical_matrix = 0.d0
 i_start = 1
 n_shell_passed = 1 
 do while (i_start .le. ao_num)
  if      (ao_l(i_start)==0)then
   cartesian_to_spherical_matrix(i_start, n_shell_passed) = 1.d0
   i_start += 1
   ao_spherical_l(0) += 1
   index_ao = ao_spherical_l(0) 
   list_shell_in_basis(index_ao,0) = n_shell_passed
   i_nucl = ao_nucl(i_start)
   ao_spherical_nucl(n_shell_passed) = i_nucl
   n_shell_in_basis_per_atom(0,i_nucl) += 1
   n_shell = n_shell_in_basis_per_atom(0,i_nucl)
   list_shell_in_basis_per_atom(n_shell,0,i_nucl) = n_shell_passed
   n_shell_passed += 1
  else if (ao_l(i_start)==1)then
   i_start_0 = i_start
   do i = 1, 3
    do k = 1, 3
     do l = 1,3 
      cartesian_to_spherical_matrix(i_start_0+i-1,n_shell_passed-1+k) += cart_to_sphe_1(l,k) * ao_overlap(i+i_start_0-1,l+i_start_0-1)
     enddo
    enddo
   enddo
   do i = 1, 3
    ao_spherical_l(1) += 1
    index_ao = ao_spherical_l(1) 
    list_shell_in_basis(index_ao,1) = n_shell_passed
    i_nucl = ao_nucl(i_start)
    n_shell_in_basis_per_atom(1,i_nucl) += 1
    n_shell = n_shell_in_basis_per_atom(1,i_nucl)
    list_shell_in_basis_per_atom(n_shell,1,i_nucl) = n_shell_passed
    n_shell_passed += 1
    i_start +=1 
   enddo
  else if (ao_l(i_start)==2)then
   i_start_0 = i_start
   do i = 1, 6
    do k = 1, 5
     do l = 1, 6
      cartesian_to_spherical_matrix(i_start_0+i-1,n_shell_passed-1+k) += cart_to_sphe_2(l,k) * ao_overlap(i+i_start_0-1,l+i_start_0-1)
     enddo
    enddo
   enddo
   do i = 1, 5
    ao_spherical_l(2) += 1
    index_ao = ao_spherical_l(2) 
    list_shell_in_basis(index_ao,2) = n_shell_passed
    i_nucl = ao_nucl(i_start)
    n_shell_in_basis_per_atom(2,i_nucl) += 1
    n_shell = n_shell_in_basis_per_atom(2,i_nucl)
    list_shell_in_basis_per_atom(n_shell,2,i_nucl) = n_shell_passed
    n_shell_passed += 1
   enddo
   i_start += 6
  else if (ao_l(i_start)==3)then
   i_start_0 = i_start
   do i = 1, 10
    do k = 1, 7
     do l = 1, 10
      cartesian_to_spherical_matrix(i_start_0+i-1,n_shell_passed-1+k) += cart_to_sphe_3(l,k) * ao_overlap(i+i_start_0-1,l+i_start_0-1)
     enddo
    enddo
   enddo
   do i = 1, 7
    ao_spherical_l(3) += 1
    index_ao = ao_spherical_l(3) 
    list_shell_in_basis(index_ao,3) = n_shell_passed
    i_nucl = ao_nucl(i_start)
    n_shell_in_basis_per_atom(3,i_nucl) += 1
    n_shell = n_shell_in_basis_per_atom(3,i_nucl)
    list_shell_in_basis_per_atom(n_shell,3,i_nucl) = n_shell_passed
    n_shell_passed += 1
   enddo
   i_start += 10
  else if (ao_l(i_start)==4)then
   i_start_0 = i_start
   do i = 1, 15
    do k = 1, 9
     do l = 1, 15
      cartesian_to_spherical_matrix(i_start_0+i-1,n_shell_passed-1+k) += cart_to_sphe_4(l,k) * ao_overlap(i+i_start_0-1,l+i_start_0-1)
     enddo
    enddo
   enddo
   do i = 1, 9
    ao_spherical_l(4) += 1
    index_ao = ao_spherical_l(4) 
    list_shell_in_basis(index_ao,4) = n_shell_passed
    i_nucl = ao_nucl(i_start)
    n_shell_in_basis_per_atom(4,i_nucl) += 1
    n_shell = n_shell_in_basis_per_atom(4,i_nucl)
    list_shell_in_basis_per_atom(n_shell,4,i_nucl) = n_shell_passed
    n_shell_passed += 1
   enddo
   i_start += 15
  else if (ao_l(i_start)==5)then
   i_start_0 = i_start
   do i = 1, 21
    do k = 1, 11
     do l = 1, 21
      cartesian_to_spherical_matrix(i_start_0+i-1,n_shell_passed-1+k) += cart_to_sphe_5(l,k) * ao_overlap(i+i_start_0-1,l+i_start_0-1)
     enddo
    enddo
   enddo
   do i = 1, 11
    ao_spherical_l(5) += 1
    index_ao = ao_spherical_l(5) 
    list_shell_in_basis(index_ao,5) = n_shell_passed
    i_nucl = ao_nucl(i_start)
    n_shell_in_basis_per_atom(5,i_nucl) += 1
    n_shell = n_shell_in_basis_per_atom(5,i_nucl)
    list_shell_in_basis_per_atom(n_shell,5,i_nucl) = n_shell_passed
    n_shell_passed += 1
   enddo
   i_start += 21
  else if (ao_l(i_start)==6)then
   i_start_0 = i_start
   do i = 1, 28
    do k = 1, 13
     do l = 1, 28
      cartesian_to_spherical_matrix(i_start_0+i-1,n_shell_passed-1+k) += cart_to_sphe_6(l,k) * ao_overlap(i+i_start_0-1,l+i_start_0-1)
     enddo
    enddo
   enddo
   do i = 1, 13
    ao_spherical_l(6) += 1
    index_ao = ao_spherical_l(6) 
    list_shell_in_basis(index_ao,6) = n_shell_passed
    i_nucl = ao_nucl(i_start)
    n_shell_in_basis_per_atom(6,i_nucl) += 1
    n_shell = n_shell_in_basis_per_atom(6,i_nucl)
    list_shell_in_basis_per_atom(n_shell,6,i_nucl) = n_shell_passed
    n_shell_passed += 1
   enddo
   i_start += 28
  else if (ao_l(i_start)==7)then
   i_start_0 = i_start
   do i = 1, 36
    do k = 1, 15
     do l = 1, 36
      cartesian_to_spherical_matrix(i_start_0+i-1,n_shell_passed-1+k) += cart_to_sphe_7(l,k) * ao_overlap(i+i_start_0-1,l+i_start_0-1)
     enddo
    enddo
   enddo
   do i = 1, 15
    ao_spherical_l(7) += 1
    index_ao = ao_spherical_l(7) 
    list_shell_in_basis(index_ao,7) = n_shell_passed
    i_nucl = ao_nucl(i_start)
    n_shell_in_basis_per_atom(7,i_nucl) += 1
    n_shell = n_shell_in_basis_per_atom(7,i_nucl)
    list_shell_in_basis_per_atom(n_shell,7,i_nucl) = n_shell_passed
    n_shell_passed += 1
   enddo
   i_start += 36
  endif
 enddo

END_PROVIDER


