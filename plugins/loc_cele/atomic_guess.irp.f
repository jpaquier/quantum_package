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
 implicit none
 double precision :: thr_loc
 thr_loc = 0.05d0
 call find_good_orb(index_ligand_orb_loc, n_orb_ligand_loc,thr_loc)
!n_orb_ligand_loc = 2 
!index_ligand_orb_loc(1) = 70
!index_ligand_orb_loc(2) = 65
!index_ligand_orb_loc(1) = 75
!index_ligand_orb_loc(2) = 71
!index_ligand_orb_loc(3) = 79
!index_ligand_orb_loc(4) = 80
 print*, 'n_orb_ligand_loc ',n_orb_ligand_loc
END_PROVIDER 

BEGIN_PROVIDER [integer, n_metal_atoms]
 implicit none
 integer :: i
 n_metal_Atoms = 0
 do i = 1, nucl_num
  if(nucl_charge(i).gt.20)then
   n_metal_atoms +=1 
  endif
 enddo
END_PROVIDER 

BEGIN_PROVIDER [integer, index_metal_atoms, (n_metal_atoms)]
 implicit none
 integer :: i,j
 j = 0
 do i = 1, nucl_num
  if(nucl_charge(i).gt.20)then
   j+=1 
   index_metal_atoms(j) = i
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

subroutine find_good_orb(list_orb, n_orb,thr_loc)
 implicit none
 double precision, intent(in) :: thr_loc
 integer, intent(out)  :: n_orb
 integer, intent(out) :: list_orb(mo_tot_num)
 integer :: i,j,k,l,m,n,p,jj
 double precision :: ovrlp(elec_beta_num)
 double precision :: accu
 integer :: iorder(elec_beta_num)
 logical :: is_selected(mo_tot_num)
 double precision :: ovrlp_selected(mo_tot_num)
 is_selected = .False.
 n_orb = 0
 if(active_guess)then
  do i = 1, n_act_orb
   print*, 'ACTIVE ORBITAL ',i,index_metal_atom_orb_loc(i)
   print*, 'OVERLAP ...'
   ovrlp = 0.d0
   do jj = 1, n_core_inact_orb
    j = list_core_inact(jj)
    iorder(j) = j
    do k = 1, n_metal_atoms
      p = index_metal_atoms(k)
      do l = 1, Nucl_N_Aos(k)
       m = Nucl_Aos_transposed(l,k)
       do n = 1, ao_num
        ovrlp(j) += mo_coef(n,j) * ao_overlap(m,n) * mo_coef(m,index_metal_atom_orb_loc(i))
       enddo
      enddo
     ovrlp(j) = -dabs(ovrlp(j))
    enddo
    print*, 'ovrlp(j)',ovrlp(j),j
   enddo
   call dsort(ovrlp,iorder,elec_beta_num)
   print*, 'MAXIMUM OVERLAPS AND CORRESPONDING ORBITALS '
   print*,iorder(1),ovrlp(1),is_selected(iorder(1)) 
   print*,iorder(2),ovrlp(2),is_selected(iorder(2)) 
   if (is_selected(iorder(1)))cycle
   is_selected(iorder(1)) = .True.
   n_orb += 1
   list_orb(n_orb) = iorder(1)
   ovrlp_selected(n_orb) = dabs(ovrlp(1))
   do j = 2, n_core_inact_orb
    if (is_selected(iorder(j)))cycle
    if(ovrlp(j)/ovrlp(1).gt.thr_loc)then
     n_orb +=1
     list_orb(n_orb) = iorder(j)
     is_selected(iorder(j)) = .True.
     ovrlp_selected(n_orb) = dabs(ovrlp(j))
    endif
   enddo
  enddo
  print*, 'n_orb (loc) = ',n_orb
  
 else
  do i = elec_beta_num+1, elec_alpha_num
   print*, 'SINGLY OCCUPIED ORBITAL ',i,index_metal_atom_orb_loc(i)
   print*, 'OVERLAP ...'
   do j = 1, elec_beta_num
    ovrlp(j) = 0.d0
    iorder(j) = j
    do k = 1, n_metal_atoms
      p = index_metal_atoms(k)
      do l = 1, Nucl_N_Aos(k)
       m = Nucl_Aos_transposed(l,k)
       do n = 1, ao_num
        ovrlp(j) += mo_coef(n,j) * ao_overlap(m,n) * mo_coef(m,index_metal_atom_orb_loc(i))
       enddo
      enddo
     ovrlp(j) = -dabs(ovrlp(j))
    enddo
    print*, 'ovrlp(j)',ovrlp(j),j
   enddo
   call dsort(ovrlp,iorder,elec_beta_num)
   print*, 'MAXIMUM OVERLAPS AND CORRESPONDING ORBITALS '
   print*,iorder(1),ovrlp(1),is_selected(iorder(1)) 
   print*,iorder(2),ovrlp(2),is_selected(iorder(2)) 
   if (is_selected(iorder(1)))cycle
   n_orb += 1
   list_orb(n_orb) = iorder(1)
   is_selected(iorder(1)) = .True.
   ovrlp_selected(n_orb) = dabs(ovrlp(1))
   do j = 2, elec_beta_num
    if (is_selected(iorder(j)))cycle
    if(ovrlp(j)/ovrlp(1).gt.thr_loc)then
     n_orb +=1
     list_orb(n_orb) = iorder(j)
     is_selected(iorder(j)) = .True.
     ovrlp_selected(n_orb) = dabs(ovrlp(j))
    endif
   enddo
  enddo
  print*, 'n_orb (loc) = ',n_orb
 endif
 
 do i = 1, n_orb
  print*, 'list_orb(i) == ',list_orb(i),ovrlp_selected(i)
 enddo
end
