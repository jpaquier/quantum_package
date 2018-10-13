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
 ! number of spherical shell in the basis
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


BEGIN_PROVIDER [double precision, cartesian_to_spherical_matrix, (ao_num,n_spherical_AOs_in_basis)]
 implicit none
 BEGIN_DOC 
 ! carthesian to spherical transformation matrix containing overlap !!!
 END_DOC
 integer :: i,k,l,i_start,i_start_0,n_shell_passed
 cartesian_to_spherical_matrix = 0.d0
 i_start = 1
 n_shell_passed = 1 
 do while (i_start .le. ao_num)
  if      (ao_l(i_start)==0)then
   cartesian_to_spherical_matrix(i_start, n_shell_passed) = 1.d0
   i_start += 1
   n_shell_passed += 1
  else if (ao_l(i_start)==1)then
   do i = 1, 3 
    cartesian_to_spherical_matrix(i_start, n_shell_passed) = 1.d0
    i_start += 1
    n_shell_passed += 1
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
   n_shell_passed += 5
   i_start += 6
  else if (ao_l(i_start)==3)then
   i_start_0 = i_start
   do i = 1, 10
    do k = 1, 7
     do l = 1, 10
      cartesian_to_spherical_matrix(i_start_0+i-1,n_shell_passed-1+k) += cart_to_sphe_2(l,k) * ao_overlap(i+i_start_0-1,l+i_start_0-1)
     enddo
    enddo
   enddo
   n_shell_passed += 7
   i_start += 10
  else if (ao_l(i_start)==4)then
   i_start_0 = i_start
   do i = 1, 15
    do k = 1, 9
     do l = 1, 15
      cartesian_to_spherical_matrix(i_start_0+i-1,n_shell_passed-1+k) += cart_to_sphe_2(l,k) * ao_overlap(i+i_start_0-1,l+i_start_0-1)
     enddo
    enddo
   enddo
   n_shell_passed += 9
   i_start += 15
  else if (ao_l(i_start)==5)then
   i_start_0 = i_start
   do i = 1, 21
    do k = 1, 11
     do l = 1, 21
      cartesian_to_spherical_matrix(i_start_0+i-1,n_shell_passed-1+k) += cart_to_sphe_2(l,k) * ao_overlap(i+i_start_0-1,l+i_start_0-1)
     enddo
    enddo
   enddo
   n_shell_passed += 9
   i_start += 21
  else if (ao_l(i_start)==6)then
   i_start_0 = i_start
   do i = 1, 28
    do k = 1, 13
     do l = 1, 28
      cartesian_to_spherical_matrix(i_start_0+i-1,n_shell_passed-1+k) += cart_to_sphe_2(l,k) * ao_overlap(i+i_start_0-1,l+i_start_0-1)
     enddo
    enddo
   enddo
   n_shell_passed += 13
   i_start += 28
  else if (ao_l(i_start)==7)then
   i_start_0 = i_start
   do i = 1, 36
    do k = 1, 15
     do l = 1, 36
      cartesian_to_spherical_matrix(i_start_0+i-1,n_shell_passed-1+k) += cart_to_sphe_2(l,k) * ao_overlap(i+i_start_0-1,l+i_start_0-1)
     enddo
    enddo
   enddo
   n_shell_passed += 15
   i_start += 36
  endif
 enddo

END_PROVIDER
