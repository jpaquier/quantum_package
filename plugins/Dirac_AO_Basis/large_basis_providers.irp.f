 BEGIN_PROVIDER [ integer, large_ao_num ]
  implicit none
  BEGIN_DOC
  !Number of large component AO
  END_DOC  
 large_ao_num = ao_num
 END_PROVIDER

 BEGIN_PROVIDER [ integer, large_ao_prim_num, (large_ao_num) ]
 implicit none
  BEGIN_DOC
  !Number of large component primitives
  END_DOC  
  large_ao_prim_num = 1
 END_PROVIDER

 BEGIN_PROVIDER [ integer, large_ao_prim_num_max ]
 implicit none
  BEGIN_DOC
  !max number of primitives of the large component
  END_DOC  
  large_ao_prim_num_max = maxval(large_ao_prim_num)
 END_PROVIDER

 BEGIN_PROVIDER [ integer, large_ao_power, (large_ao_num,3) ]
 &BEGIN_PROVIDER [ integer, large_ao_l, (large_ao_num) ]
 &BEGIN_PROVIDER [ integer, large_ao_nucl, (large_ao_num) ]
 &BEGIN_PROVIDER [ double precision, large_ao_expo_ordered_transp, (large_ao_prim_num_max,large_ao_num) ]
 implicit none
 BEGIN_DOC
 !Nucleus on which the AOs are centered
 !large_ao_power is the AO power for x,y,z for the AOs of 
 ! the large component ao basis
 !large_ao_l = l value of : (a+b+c) in x^a y^b z^c for the AOs of 
 ! the large component basis
 !large_ao_nucl is the nucleus on which the large component AO is located
 !large_ao_expo_ordered_transp is the transposed ordered large_ao_expo (which
 ! are not defined given that they are the non-relativistic  ao_expo ) 
 END_DOC  
  large_ao_nucl = ao_nucl
  large_ao_power = ao_power
  large_ao_l = ao_l
  large_ao_expo_ordered_transp = ao_expo_ordered_transp
 END_PROVIDER


 BEGIN_PROVIDER [ double precision, large_ao_coef, (large_ao_num,large_ao_prim_num_max)]
 &BEGIN_PROVIDER [ double precision, large_ao_coef_normalized_ordered_transp, (large_ao_prim_num_max,large_ao_num)]
 BEGIN_DOC
  !Normalized and ordered coefficient of the large component AOs 
  END_DOC
  double precision               :: large_norm,large_overlap_x,large_overlap_y,large_overlap_z,C_A(3), c
  integer                        :: l, powA(3), nz
  integer                        :: i,j,k
  do i = 1, large_ao_num
   large_ao_coef(i,1) = 1
  enddo
  nz=100
  C_A(1) = 0.d0
  C_A(2) = 0.d0
  C_A(3) = 0.d0
  large_ao_coef_normalized_ordered_transp = 0.d0
  do i=1, large_ao_num
   powA(1) = large_ao_power(i,1)
   powA(2) = large_ao_power(i,2)
   powA(3) = large_ao_power(i,3)
   do j = 1, large_ao_prim_num(i)
    call overlap_gaussian_xyz(C_A,C_A,large_ao_expo_ordered_transp(j,i),large_ao_expo_ordered_transp(j,i),powA,powA,large_overlap_x,large_overlap_y,large_overlap_z,large_norm,nz)
    large_ao_coef_normalized_ordered_transp(j,i) = large_ao_coef(i,j)/sqrt(large_norm)
   enddo
  enddo
 END_PROVIDER



