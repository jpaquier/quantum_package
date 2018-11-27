BEGIN_PROVIDER [logical, active_guess]
 implicit none
 active_guess = .True.
END_PROVIDER 

BEGIN_PROVIDER [integer, n_ligand_loc]
 implicit none
 n_ligand_loc = 1
END_PROVIDER 

BEGIN_PROVIDER [integer, n_orb_metal_atoms_loc]
 implicit none
 if(active_guess)then
  n_orb_metal_atoms_loc = n_act_orb
 else
  n_orb_metal_atoms_loc = n_metal_atoms
 endif
END_PROVIDER 

BEGIN_PROVIDER [integer, index_metal_atom_orb_loc, (n_orb_metal_atoms_loc)]
 implicit none
 integer :: i,j
 if(active_guess)then
  do i = 1, n_act_orb
   index_metal_atom_orb_loc(i) = list_act(i)
  enddo
 else
  j = 0
  do i = elec_beta_num+1,elec_alpha_num
   j += 1
   index_metal_atom_orb_loc(j) = i
  enddo
 endif
END_PROVIDER 


 BEGIN_PROVIDER [integer, n_orb_ligand_loc]
&BEGIN_PROVIDER [integer, index_ligand_orb_loc, (mo_tot_num)]
&BEGIN_PROVIDER [integer, index_ligand_orb_loc_sorted, (mo_tot_num)]
 implicit none
 BEGIN_DOC
 ! n_orb_ligand_loc = number of ligand-like orbitals that have been selected by the maxoverlap procedure
 ! index_ligand_orb_loc(i) = index of the ith most overlapping ligand-like orbital 
 ! index_ligand_orb_loc_sorted(i) = same that index_ligand_orb_loc but sorted per increasing number of orbitals 
 END_DOC
 double precision :: thr_loc
 thr_loc = 0.1d0
 print*,''
 print*,'Selecting the interesting doubly occupied orbitals ...'
 call find_good_orb(index_ligand_orb_loc, n_orb_ligand_loc,threshold_overlap,list_inact,n_inact_orb)
 integer :: i
 integer, allocatable :: iorder(:)
 allocate(iorder(n_orb_ligand_loc))
 
 do i = 1, n_orb_ligand_loc
  iorder(i) = i
  index_ligand_orb_loc_sorted(i) = index_ligand_orb_loc(i)
 enddo
 call isort(index_ligand_orb_loc_sorted,iorder,n_orb_ligand_loc)
END_PROVIDER 

 BEGIN_PROVIDER [integer, n_orb_ligand_virt_loc]
&BEGIN_PROVIDER [integer, index_ligand_virt_orb_loc, (mo_tot_num)]
&BEGIN_PROVIDER [integer, index_ligand_virt_orb_loc_sorted, (mo_tot_num)]
 implicit none
 BEGIN_DOC
 ! n_orb_ligand_virt_loc = number of ligand_virt-like orbitals that have been selected by the maxoverlap procedure
 ! index_ligand_virt_orb_loc(i) = index of the ith most overlapping ligand_virt-like orbital 
 ! index_ligand_virt_orb_loc_sorted(i) = same that index_ligand_virt_orb_loc but sorted per increasing number of orbitals 
 END_DOC
 double precision :: thr_loc
 
 print*,''
 print*,'Selecting the interesting virtuals orbitals ...'
 call find_good_orb(index_ligand_virt_orb_loc, n_orb_ligand_virt_loc,threshold_overlap,list_virt,n_virt_orb)
 integer :: i
 integer, allocatable :: iorder(:)
 allocate(iorder(n_orb_ligand_virt_loc))
 
 do i = 1, n_orb_ligand_virt_loc
  iorder(i) = i
  index_ligand_virt_orb_loc_sorted(i) = index_ligand_virt_orb_loc(i)
 enddo
 call isort(index_ligand_virt_orb_loc_sorted,iorder,n_orb_ligand_virt_loc)
 print*,'interesting virtual orbitals selected !'
END_PROVIDER 


 BEGIN_PROVIDER [integer, n_metal_atoms]
&BEGIN_PROVIDER [integer, n_non_metal_atoms]
 implicit none
 integer :: i
 n_metal_atoms = 0
 n_non_metal_atoms = 0
 do i = 1, nucl_num
  if(nucl_charge(i).gt.20)then
   n_metal_atoms +=1 
  else 
   n_non_metal_atoms += 1
  endif
 enddo
END_PROVIDER 

 BEGIN_PROVIDER [integer, index_metal_atoms, (n_metal_atoms)]
&BEGIN_PROVIDER [integer, index_non_metal_atoms, (n_non_metal_atoms)]
 implicit none
 integer :: i,j_metal,j_non_metal
 j_metal = 0
 j_non_metal = 0
 do i = 1, nucl_num
  if(nucl_charge(i).gt.20)then
   j_metal+=1 
   index_metal_atoms(j_metal) = i
  else 
   j_non_metal+=1 
   index_non_metal_atoms(j_non_metal) = i
  endif
 enddo
END_PROVIDER 

BEGIN_PROVIDER [double precision, guess_loc_metal_mos,(ao_num,n_orb_metal_atoms_loc)]
 implicit none
 integer :: i,j,k,m,l
 guess_loc_metal_mos = 0.d0
 do i = 1, n_orb_metal_atoms_loc
  j = index_metal_atom_orb_loc(i)
  k = index_metal_atoms(1)
   do l = 1, Nucl_N_Aos(k)
    m = Nucl_Aos_transposed(l,k)
    guess_loc_metal_mos(m,i) =  mo_coef(m,j)
   enddo
 enddo
END_PROVIDER 

BEGIN_PROVIDER [double precision, guess_loc_ligand_mos,(ao_num,n_orb_ligand_loc)]
 implicit none
 integer :: i,j,k,m,l
 guess_loc_ligand_mos = 0.d0
 do i = 1, n_orb_ligand_loc
  j = index_ligand_orb_loc(i)
  do k = 1, nucl_num
   if(k == index_metal_atoms(1))cycle
   do l = 1, Nucl_N_Aos(k)
    m = Nucl_Aos_transposed(l,k)
    guess_loc_ligand_mos(m,i) =  mo_coef(m,j)
   enddo
  enddo
 enddo
END_PROVIDER 

subroutine find_good_orb(list_orb, n_orb,thr_loc,list_orb_in,n_orb_in)
 implicit none
 use bitmasks
 double precision, intent(in) :: thr_loc
 integer, intent(in)  :: list_orb_in(n_orb_in),n_orb_in
 integer, intent(out)  :: n_orb
 integer, intent(out) :: list_orb(mo_tot_num)
 integer  :: list_orb_tmp((N_int * bit_kind_size))
 integer :: i,j,k,l,m,n,p,jj,k1,l1,m1,kk,kk1
 double precision :: ovrlp(mo_tot_num)
 double precision :: accu
 integer :: iorder(mo_tot_num)
 logical :: is_selected(mo_tot_num)
 integer(bit_kind) :: key(N_int)
 double precision :: mo_coef_metal(ao_num,mo_tot_num)
 double precision :: mo_coef_list(ao_num,mo_tot_num)

 
 ! normalization of the METAL-CENTERED ACTIVE ORBITALS :
 ! THESE ORBITALS ARE NOTHING BUT THE ACTIVE ORBITALS WITH AO COEFFICIENTS ONLY ON THE METALS
!print*,'METALLIC ORBITALS NORMALIZATION '
 mo_coef_metal = 0.d0
 do i = 1, n_act_orb
  j = list_act(i)
  accu = 0.d0
  do kk = 1, n_metal_atoms ! you run on the metallic atoms 
   k = index_metal_atoms(kk)
   do l = 1, Nucl_N_Aos(k) ! you run on the AO attached to each metallic atoms 
    m = Nucl_Aos_transposed(l,k) ! m = AO attached to the Lth AO of the  Kth metalic atom 
    mo_coef_metal(m,j) = mo_coef(m,j)
    do kk1 = 1, n_metal_atoms ! you run on the metallic atoms 
     k1 = index_metal_atoms(kk1)
     do l1 = 1, Nucl_N_Aos(k1) ! you run on the AO attached to each metallic atoms 
      m1 = Nucl_Aos_transposed(l1,k1) ! m = AO attached to the Lth AO of the  Kth metalic atom 
      ! you compute the overlap
      accu += mo_coef(m,j) * ao_overlap(m1,m) * mo_coef(m1,j)
     enddo
    enddo
   enddo
  enddo
  accu = 1.d0/dsqrt(accu)
! print*,'accu = ',accu
  do m = 1, ao_num
  ! you normalize
   mo_coef_metal(m,j) = mo_coef_metal(m,j)/accu
  enddo
 enddo
 
!print*,''
!print*,''
!print*,'NON METALLIC ORBITALS NORMALIZATION '
!! YOU DO THE OPPOSIT FOR THE INPUT ORBITALS : ANYTHING BUT THE NON-METALLIC ATOMS
!mo_coef_list = 0.d0
!do i = 1, n_orb_in
! j = list_orb_in(i)
! accu = 0.d0
! do kk = 1, n_non_metal_atoms ! you run on the non_metallic atoms 
!  k = index_non_metal_atoms(kk)
!  do l = 1, Nucl_N_Aos(k) ! you run on the AO attached to each non_metallic atoms 
!   m = Nucl_Aos_transposed(l,k) ! m = AO attached to the Lth AO of the  Kth non_metalic atom 
!   mo_coef_list(m,j) =  mo_coef(m,j)
!   do kk1 = 1, n_non_metal_atoms ! you run on the non_metallic atoms 
!    k1 = index_non_metal_atoms(kk1)
!    do l1 = 1, Nucl_N_Aos(k1) ! you run on the AO attached to each non_non_metallic atoms 
!     m1 = Nucl_Aos_transposed(l1,k1) ! m = AO attached to the Lth AO of the  Kth non_metalic atom 
!     ! you compute the overlap
!     accu += mo_coef(m,j) * ao_overlap(m1,m) * mo_coef(m1,j)
!    enddo
!   enddo
!  enddo
! enddo
! print*,'accu = ',i,j,accu
! accu = 1.d0/dsqrt(accu)
! do m = 1, ao_num
! ! you normalize
!  mo_coef_list(m,j) = mo_coef_list(m,j)/accu
! enddo
!enddo

 mo_coef_list = mo_coef


 is_selected = .False.
! FOR EACH STRONGLY METALLIC ACTIVE ORBITALS 
 do i = 1, n_act_orb
  print*, 'ACTIVE ORBITAL ',i,index_metal_atom_orb_loc(i)
  print*, 'OVERLAP ...'
  ovrlp = 0.d0
  do j = 1, mo_tot_num
   iorder(j) = j
  enddo
! YOU COMPUTE THE OVERLAP WITH ALL THE INPUT ORBITALS 
  do jj = 1, n_orb_in
   j = list_orb_in(jj)
   do n = 1, ao_num
    do m = 1, ao_num
     ovrlp(j) += mo_coef_list(n,j) * ao_overlap(n,m) * mo_coef_metal(m,list_act(i))
    enddo
   enddo
   ovrlp(j) = -dabs(ovrlp(j))
   print*,j,dabs(ovrlp(j)),thr_loc
  enddo
!! YOU SORT THE OVERLAPS 
  call dsort(ovrlp,iorder,mo_tot_num)
  print*, 'MAXIMUM OVERLAPS AND CORRESPONDING ORBITALS '
 !print*,iorder(1),ovrlp(1)
 !print*,iorder(2),ovrlp(2)
!! YOU SELECT 
  do jj = 1, n_orb_in
    print*,iorder(jj),dabs(ovrlp(jj)),thr_loc
   if(dabs(ovrlp(jj)).gt.thr_loc)then
    j = iorder(jj)
    is_selected(j) = .True.
   endif
  enddo
 enddo
 key = 0_bit_kind
 do j = 1, mo_tot_num
  if(is_selected(j))then
   call set_bit_to_integer(j,key,N_int)
  endif
 enddo
 call bitstring_to_list( key, list_orb_tmp, n_orb, N_int)
 print*,'n_orb = ',n_orb
 do i = 1, n_orb
  list_orb(i) = list_orb_tmp(i)
  print*,list_orb(i)
 enddo
end
