 BEGIN_PROVIDER [ integer, dirac_ao_num ]
 &BEGIN_PROVIDER [ integer, dirac_ao_prim_num, (ao_num + small_ao_num) ]
 implicit none
  BEGIN_DOC
  ! Concatenation of the large and small components orbital properties
  ! in general arrays, for use in the bi-electronic integrals
  END_DOC
  integer                        :: i,j,k
  dirac_ao_num = (ao_num + small_ao_num)
  do i = 1, dirac_ao_num
   if (i .le. ao_num) then
   dirac_ao_prim_num (i) = ao_prim_num(i)
   else
   j = i - ao_num
   dirac_ao_prim_num (i) = small_ao_prim_num(j)
   endif 
  enddo
  END_PROVIDER

 BEGIN_PROVIDER [ integer, dirac_ao_prim_num_max ]
 implicit none
 BEGIN_DOC
 ! max number of primitives
 END_DOC
 dirac_ao_prim_num_max = maxval(dirac_ao_prim_num)
 END_PROVIDER

 BEGIN_PROVIDER [ integer, dirac_mo_tot_num ]
  implicit none
  BEGIN_DOC
  ! Number of small component MOs
  END_DOC
  dirac_mo_tot_num = mo_tot_num + small_mo_tot_num
  ASSERT (small_mo_tot_num > 0)
 END_PROVIDER

 BEGIN_PROVIDER [ integer, dirac_ao_nucl, (dirac_ao_num) ]
 &BEGIN_PROVIDER [ double precision, dirac_ao_coef_normalized_ordered_transp, (dirac_ao_prim_num_max,dirac_ao_num) ]
 &BEGIN_PROVIDER [ double precision, dirac_ao_expo_ordered_transp, (dirac_ao_prim_num_max,dirac_ao_num) ]
 &BEGIN_PROVIDER [ integer, dirac_ao_power, (dirac_ao_num,3) ]
 &BEGIN_PROVIDER [ integer, dirac_ao_l, (dirac_ao_num) ]
  implicit none
  BEGIN_DOC
  ! Concatenation of the large and small components orbital properties
  ! in general arrays, for use in the bi-electronic integrals
  END_DOC
  integer                        :: i,j,k,l
  do i = 1, dirac_ao_num
   if (i <= ao_num) then
    dirac_ao_nucl(i) = ao_nucl(i)
    dirac_ao_l(i) = ao_l(i)
    do k = 1, ao_prim_num(i)
     dirac_ao_coef_normalized_ordered_transp(k,i) = ao_coef_normalized_ordered_transp(k,i)
     dirac_ao_expo_ordered_transp(k,i) = ao_expo_ordered_transp(k,i)
    enddo
    do l = 1, 3
     dirac_ao_power(i,l) = ao_power(i,l)
    enddo
   else
    j = i - ao_num
    dirac_ao_nucl(i) = small_ao_nucl(j)
    dirac_ao_l(i) = small_ao_l(j)
    do k = 1, small_ao_prim_num(j)
     dirac_ao_coef_normalized_ordered_transp(k,i) = small_ao_coef_normalized_ordered_transp(k,j)
     dirac_ao_expo_ordered_transp(k,i) = small_ao_expo_ordered_transp(k,j)
    enddo
    do l = 1, 3
     dirac_ao_power(i,l) = small_ao_power(j,l)
    enddo
   endif
  enddo
 END_PROVIDER

 BEGIN_PROVIDER [ double precision, dirac_ao_integrals_threshold  ]
  implicit none
  BEGIN_DOC
  ! If |<pq|rs>| < ao_integrals_threshold then <pq|rs> is zero
  END_DOC
 dirac_ao_integrals_threshold =  1.0E-015
 END_PROVIDER

 BEGIN_PROVIDER [double precision, dirac_ao_overlap_abs, (dirac_ao_num,dirac_ao_num) ]
 implicit none
  BEGIN_DOC
  ! Concatenation of the large and small component
  ! overlap_abs 
  END_DOC
  integer                        :: i,j,k,l
  do i = 1, dirac_ao_num
   if (i <= ao_num) then
    do j = 1, dirac_ao_num
     if (j <= ao_num) then
      dirac_ao_overlap_abs(i,j) = ao_overlap_abs(i,j)
     else
      dirac_ao_overlap_abs(i,j) = 0.d0
     endif
    enddo
   else
    k = i - ao_num
    do j = 1, dirac_ao_num
     if (j <= ao_num) then
     dirac_ao_overlap_abs(i,j) = 0.d0
     else
     l = j - ao_num
     dirac_ao_overlap_abs(i,j) = small_ao_overlap_abs(k,l)
     endif
    enddo
   endif
  enddo
 END_PROVIDER


