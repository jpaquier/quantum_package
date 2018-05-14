 BEGIN_PROVIDER [ integer, number_of_small_component_expo_per_shell_per_atom,(0:7,nucl_num)] 
&BEGIN_PROVIDER [ integer, nmax_of_small_component_expo] 
&BEGIN_PROVIDER [ integer, number_of_small_component_ao_per_shell_per_atom,(0:7,nucl_num)]
&BEGIN_PROVIDER [ integer, number_of_small_component_ao_per_atom,(nucl_num)]
&BEGIN_PROVIDER [ integer, small_ao_num]
 implicit none
 BEGIN_DOC
! number_of_small_component_expo_per_shell_per_atom = number of ao exponents in
! the shell "i" of the nucleus "j" for the small component basis 
 END_DOC                                                                                                 
 integer :: i,i_count,l_type,j,index_ao,k,l,l_count                                                      
 integer :: imax
 imax = 0
 number_of_small_component_expo_per_shell_per_atom(:,:) = 0
 number_of_small_component_ao_per_shell_per_atom(:,:) = 0
 number_of_small_component_ao_per_atom(:) = 0
 small_ao_num = 0
 do i = 1, nucl_num
  do l_type = 0,7
   if(Nucl_num_l_type_Aos(l_type,i) == 0)cycle                                                           
   do j = 1, Nucl_num_l_type_Aos(l_type,i)
    ! "physical" index of all the AO corresponding to the shell type l_type                              
    ! attached to atom i                                                                                 
    index_ao = Nucl_list_l_type_Aos(j,l_type,i)                                                          
    if (l_type == 0) then
     number_of_small_component_expo_per_shell_per_atom(l_type+1,i) += 1
    else
     if (j == 1) then
      number_of_small_component_expo_per_shell_per_atom(l_type-1,i) += 1
      number_of_small_component_expo_per_shell_per_atom(l_type+1,i) += 1
     elseif ( j /= 1 .and. ao_expo_ordered_transp(1,index_ao) /= ao_expo_ordered_transp(1,Nucl_list_l_type_Aos(j-1,l_type,i))) then
      number_of_small_component_expo_per_shell_per_atom(l_type-1,i) += 1
      number_of_small_component_expo_per_shell_per_atom(l_type+1,i) += 1
     endif                                  
    endif
   enddo
  enddo 
  do l_type = 0,7  
   if(imax .lt. number_of_small_component_expo_per_shell_per_atom(l_type,i)) then
    imax =  number_of_small_component_expo_per_shell_per_atom(l_type,i)
   endif 
   number_of_small_component_ao_per_shell_per_atom(l_type,i) += number_of_small_component_expo_per_shell_per_atom(l_type,i)*(l_type+1)*(l_type+2)/2 
   number_of_small_component_ao_per_atom(i) += number_of_small_component_ao_per_shell_per_atom(l_type,i) 
  enddo
 small_ao_num += number_of_small_component_ao_per_atom(i)
 enddo
 nmax_of_small_component_expo = imax
 END_PROVIDER

 BEGIN_PROVIDER [ double precision, small_component_expo_per_shell_per_atom,(nmax_of_small_component_expo,0:7,nucl_num)]
 implicit none
 BEGIN_DOC
! Arrays containing the values of the exponents for the small component in 
! Unrestricted Kinetic Balance, taken from the uncontracted gaussian forming
! the large component basis
 END_DOC                                                                       
 integer :: i,l_type,j,index_ao,k,l,l_count
 integer :: expo_count(0:7) 
 integer :: iorder(nmax_of_small_component_expo)
 double precision :: small_component_expo_per_shell(nmax_of_small_component_expo,0:7)
 double precision :: small_component_expo(nmax_of_small_component_expo)
 small_component_expo_per_shell_per_atom(:,:,:) = 0
 small_component_expo_per_shell(:,:) = 0
 do i = 1, nucl_num
  iorder(:) = 0
  expo_count(:) = 0
  do l_type = 0,7
   if(Nucl_num_l_type_Aos(l_type,i) == 0)cycle
   do j = 1, Nucl_num_l_type_Aos(l_type,i)
    ! "physical" index of all the AO corresponding to the shell type l_type
    ! attached to atom i 
    index_ao = Nucl_list_l_type_Aos(j,l_type,i)
    if (l_type == 0) then
     expo_count(l_type+1) += 1
     small_component_expo_per_shell(expo_count(l_type+1),l_type+1) = ao_expo_ordered_transp(1,index_ao)
    else
     if (j == 1) then
      expo_count(l_type-1) += 1                              
      small_component_expo_per_shell(expo_count(l_type-1),l_type-1) = ao_expo_ordered_transp(1,index_ao) 
      expo_count(l_type+1) += 1
      small_component_expo_per_shell(expo_count(l_type+1),l_type+1) = ao_expo_ordered_transp(1,index_ao) 
     elseif ( j /= 1 .and. ao_expo_ordered_transp(1,index_ao) /=ao_expo_ordered_transp(1,Nucl_list_l_type_Aos(j-1,l_type,i))) then
      expo_count(l_type-1) += 1
      small_component_expo_per_shell(expo_count(l_type-1),l_type-1) = ao_expo_ordered_transp(1,index_ao)
      expo_count(l_type+1) += 1
      small_component_expo_per_shell(expo_count(l_type+1),l_type+1) = ao_expo_ordered_transp(1,index_ao) 
     endif
    endif
   enddo
  enddo 
 ! Second loop to sort the ao exponent in the right order
  do l_type = 0,7
   small_component_expo(:)=0
   do k = 1, expo_count(l_type)
    iorder(k) = k
    small_component_expo(k) = small_component_expo_per_shell(k,l_type)
   enddo
    call dsort(small_component_expo,iorder,expo_count(l_type))
   do l =1, expo_count(l_type)
    small_component_expo_per_shell_per_atom(l,l_type,i) = small_component_expo(expo_count(l_type)+1-l)
   enddo 
  enddo 
 enddo
 END_PROVIDER


 BEGIN_PROVIDER [ integer, small_ao_l, (small_ao_num) ]
&BEGIN_PROVIDER [ integer, small_ao_l_max  ]
&BEGIN_PROVIDER [ integer, small_ao_nucl, (small_ao_num) ]  
&BEGIN_PROVIDER [ double precision, small_ao_expo, (small_ao_num) ]
 implicit none
 BEGIN_DOC
! small_ao_l = l value of the small component AO: a+b+c in x^a y^b z^c
 END_DOC
 integer :: i,l,l_type,j,j_count,k,k_count
 j_count = 0
 k_count = 0
 do i = 1,nucl_num
  do l_type = 0,7
   do j = 1, number_of_small_component_ao_per_shell_per_atom(l_type,i)  
    j_count += 1
    small_ao_l(j_count) = l_type
    small_ao_nucl(j_count) = i 
   enddo
   do k = 1, number_of_small_component_expo_per_shell_per_atom(l_type,i)
    do l =1, (l_type+1)*(l_type+2)/2
     k_count += 1
     small_ao_expo(k_count) = small_component_expo_per_shell_per_atom(k,l_type,i)
    enddo
   enddo
  enddo
 enddo
 small_ao_l_max = maxval(small_ao_l)
 END_PROVIDER


 BEGIN_PROVIDER [ integer, small_ao_power, (small_ao_num,3)]
 implicit none
 BEGIN_DOC
  The AO power for x,y,z for the small component ao basis
 END_DOC
 integer :: i,j,j_count,k,k_count,k_j,l,l_type,l_count
 small_ao_power(:,:) = 0
 k_count = 0
 k_j = 0
 do i = 1,nucl_num
  do l_type = 0,7
   do k = 1, number_of_small_component_expo_per_shell_per_atom(l_type,i)
    do j =1, l_type + 1
     if (j == 1) then
      k_count += 1
      k_j = k_count
      small_ao_power(k_count,1) = l_type
     else
      do l = 1, j
       k_count += 1
       l_count += 1  
       small_ao_power(k_count,1) = small_ao_power(k_j,1) - j +1
       small_ao_power(k_count,2) = small_ao_power(k_j,2) + j -l 
       small_ao_power(k_count,3) = small_ao_power(k_j,3) +l -1 
      enddo
     endif
    enddo
   enddo
  enddo 
 enddo  
 END_PROVIDER

 
 BEGIN_PROVIDER [ double precision, small_ao_coef_normalized, (small_ao_num)]
&BEGIN_PROVIDER [ double precision, small_ao_coef, (small_ao_num)]
 BEGIN_DOC
  ! Normalized and ordered coefficient of the small component AOs in the
  ! unrestricted kinetic balance
  END_DOC
  double precision               :: small_norm,small_overlap_x,small_overlap_y,small_overlap_z,C_A(3), c
  integer                        :: l, powA(3), nz
  integer                        :: i,j,k
  do i = 1, small_ao_num
   small_ao_coef(i) = 1.d0
  enddo
  nz=100
  C_A(1) = 0.d0
  C_A(2) = 0.d0
  C_A(3) = 0.d0
  small_ao_coef_normalized = 0.d0
  do i=1, small_ao_num
   powA(1) = small_ao_power(i,1)
   powA(2) = small_ao_power(i,2)
   powA(3) = small_ao_power(i,3)
   call overlap_gaussian_xyz(C_A,C_A,small_ao_expo(i),small_ao_expo(i),powA,powA,small_overlap_x,small_overlap_y,small_overlap_z,small_norm,nz)
   small_ao_coef_normalized(i) = small_ao_coef(i)/sqrt(small_norm)
  enddo
 END_PROVIDER

 BEGIN_PROVIDER [ double precision, small_ao_ortho_canonical_coef, (small_ao_num,small_ao_num)]
&BEGIN_PROVIDER [ integer, small_ao_ortho_canonical_num ]
  implicit none
  BEGIN_DOC
! matrix of the coefficients of the small component mos generated by the 
! orthonormalization by the S^{-1/2} canonical transformation of the small
! component aos. small_ao_ortho_canonical_coef(i,j) = coefficient of the ith ao on the jth
! small_ao_ortho_canonical orbital
  END_DOC
  integer :: i
  small_ao_ortho_canonical_coef = 0.d0
  do i=1, small_ao_num
    small_ao_ortho_canonical_coef(i,i) = 1.d0
  enddo
  small_ao_ortho_canonical_num = small_ao_num
  call ortho_canonical(small_ao_overlap,size(small_ao_overlap,1), small_ao_num,small_ao_ortho_canonical_coef, size(small_ao_ortho_canonical_coef,1), small_ao_ortho_canonical_num)
 END_PROVIDER
  
 BEGIN_PROVIDER [double precision, small_ao_ortho_canonical_overlap, (small_ao_ortho_canonical_num, small_ao_ortho_canonical_num)]
  implicit none
  BEGIN_DOC
 ! overlap matrix of the small_ao_ortho_canonical.
 ! Expected to be the Identity
  END_DOC
  integer                        :: i,j,k,l
  double precision               :: c
  do j=1, small_ao_ortho_canonical_num
    do i=1, small_ao_ortho_canonical_num
      small_ao_ortho_canonical_overlap(i,j) = 0.d0
    enddo
  enddo
  do j=1, small_ao_ortho_canonical_num
    do k=1, small_ao_num
      c = 0.d0
      do l=1, small_ao_num
        c +=  small_ao_ortho_canonical_coef(l,j) * small_ao_overlap(l,k)
      enddo
      do i=1, small_ao_ortho_canonical_num
        small_ao_ortho_canonical_overlap(i,j) += small_ao_ortho_canonical_coef(k,i) * c
      enddo
    enddo
  enddo
 END_PROVIDER

 BEGIN_PROVIDER [ integer, small_mo_tot_num ]
  implicit none
  BEGIN_DOC
  ! Number of small component MOs
  END_DOC
  small_mo_tot_num = small_ao_ortho_canonical_num
  ASSERT (small_mo_tot_num > 0)
 END_PROVIDER


 BEGIN_PROVIDER [ double precision, small_mo_coef, (small_ao_num,small_mo_tot_num) ]
  implicit none
  BEGIN_DOC
  ! Molecular orbital coefficients on small component AO basis set
  ! small_mo_coef(i,j) = coefficient of the ith small component ao on the jth
  ! small component mo
  ! mo_label : Label characterizing the MOS (local, canonical, natural, etc)
  END_DOC
  integer                        :: i, j
  double precision, allocatable  :: buffer(:,:)
  logical                        :: exists
  PROVIDE ezfio_filename
  ! Orthonormalized AO basis
  do i=1,small_mo_tot_num
    do j=1,small_ao_num
      small_mo_coef(j,i) = small_ao_ortho_canonical_coef(j,i)
    enddo
  enddo
 END_PROVIDER

 BEGIN_PROVIDER [ double precision, dirac_mo_coef,(2*(ao_num+small_ao_num),2*(mo_tot_num+small_mo_tot_num))
 implicit none
  BEGIN_DOC
  ! Molecular orbital coefficients on AO basis set
  ! dirac_mo_coef(i,j) = coefficient of the ith ao on the jth mo
  ! mo_label : Label characterizing the MOS (local, canonical, natural, etc)
  END_DOC
  integer                        :: i, j, k, l
  double precision, allocatable  :: buffer(:,:)
  logical                        :: exists
  PROVIDE ezfio_filename
  do i=1, 2*(mo_tot_num+small_mo_tot_num)
   if (i .le. mo_tot_num) then
    l = i - 0
    do j=1, 2*(ao_num+small_ao_num)
     if (j .le. ao_num) then
      k = j - 0
      dirac_mo_coef(j,i) = ao_ortho_canonical_coef(k,l)
     elseif (j .gt. mo_tot_num) then
      dirac_mo_coef(j,i) = 0
     endif
    enddo
   elseif (i.gt. mo_tot_num .and. i .le. 2*mo_tot_num) then
    l = i - mo_tot_num
    do j=1, 2*(ao_num+small_ao_num)
     if (j .le. ao_num) then
      dirac_mo_coef(j,i) = 0
     elseif (j .gt. ao_num .and. j .le. 2*ao_num) then
      k = j - ao_num
      dirac_mo_coef(j,i) = ao_ortho_canonical_coef(k,l)
     elseif (j .gt. 2*ao_num) then
      dirac_mo_coef(j,i) = 0
     endif
    enddo
   elseif (i.gt. 2*mo_tot_num .and. i .le. (2*mo_tot_num+small_mo_tot_num)) then
    l = i - 2*mo_tot_num
    do j=1, 2*(ao_num+small_ao_num)
     if (j .le. 2*ao_num) then
      dirac_mo_coef(j,i) = 0
     elseif (j .gt. 2*ao_num .and. j .le. (2*ao_num+small_ao_num)) then
      k = j - 2*ao_num
      dirac_mo_coef(j,i) = small_ao_ortho_canonical_coef(k,l)
     elseif (j .gt. (2*ao_num+small_ao_num)) then
      dirac_mo_coef(j,i) = 0
     endif
    enddo
   else 
    l = i - (2*mo_tot_num+small_mo_tot_num)
    do j=1, 2*(ao_num+small_ao_num)
     if (j .le. (2*ao_num+small_ao_num)) then
      dirac_mo_coef(j,i) = 0
     else
      k = j - (2*ao_num+small_ao_num)
      dirac_mo_coef(j,i) = small_ao_ortho_canonical_coef(k,l)
     endif
    enddo
   endif
  enddo
 END_PROVIDER

 subroutine dirac_ao_to_mo(A_ao,LDA_ao,A_mo,LDA_mo)
  implicit none
  BEGIN_DOC
  ! Transform A from the AO basis to the MO basis
  !
  ! Ct.A_ao.C
  END_DOC
  integer, intent(in)            :: LDA_ao,LDA_mo
  double precision, intent(in)   :: A_ao(LDA_ao,2*(ao_num+small_ao_num))
  double precision, intent(out)  :: A_mo(LDA_mo,2*(mo_tot_num+small_mo_tot_num))
  double precision, allocatable  :: T(:,:)

  allocate ( T(2*(ao_num+small_ao_num),2*(mo_tot_num+small_mo_tot_num)) )
  !DIR$ ATTRIBUTES ALIGN : $IRP_ALIGN :: T

  call dgemm('N','N', 2*(ao_num+small_ao_num), 2*(mo_tot_num+small_mo_tot_num), 2*(ao_num+small_ao_num), &
      1.d0, A_ao,LDA_ao,                                             &
      dirac_mo_coef, size(dirac_mo_coef,1),                          &
      0.d0, T, size(T,1))                                            

  call dgemm('T','N', 2*(mo_tot_num+small_mo_tot_num), 2*(mo_tot_num+small_mo_tot_num), 2*(ao_num+small_ao_num), &
      1.d0, dirac_mo_coef,size(dirac_mo_coef,1),                     &
      T, 2*(ao_num+small_ao_num),                                                     &
      0.d0, A_mo, size(A_mo,1))

  deallocate(T)
 end


 BEGIN_PROVIDER [double complex, dirac_fock_matrix_eigenvalues,(2*(mo_tot_num+small_mo_tot_num))]
 &BEGIN_PROVIDER [double complex, dirac_fock_matrix_eigenvectors, (2*(mo_tot_num+small_mo_tot_num),2*(mo_tot_num+small_mo_tot_num))]
 implicit none 
 integer :: n,nmax
 double complex :: eigenvectors(2*(mo_tot_num+small_mo_tot_num),2*(mo_tot_num+small_mo_tot_num)), eigenvalues( 2*(mo_tot_num+small_mo_tot_num))
 n = 2*(mo_tot_num+small_mo_tot_num)
 nmax = n
 dirac_fock_matrix_eigenvalues = 0.d0
 dirac_fock_matrix_eigenvectors = 0.d0
 eigenvalues = 0.d0 
 eigenvectors = 0.d0
 call lapack_diag_complex(eigenvalues,eigenvectors,dirac_mo_mono_elec_integral,nmax,n)
 dirac_fock_matrix_eigenvalues = eigenvalues
 dirac_fock_matrix_eigenvectors = eigenvectors
 END_PROVIDER 



