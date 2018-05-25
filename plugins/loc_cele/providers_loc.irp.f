BEGIN_PROVIDER [integer, n_orb_loc ]
 implicit none
 n_orb_loc = n_orb_metal_atoms_loc + n_orb_ligand_loc
END_PROVIDER 

BEGIN_PROVIDER [integer, index_loc_orb, (n_orb_loc)]
 implicit none
 integer :: i,j
 j = 0
 do i = 1, n_orb_ligand_loc
  j += 1
  index_loc_orb(j) = index_ligand_orb_loc(i)
 enddo
 do i = 1, n_orb_metal_atoms_loc
  j+=1
  index_loc_orb(j) = index_metal_atom_orb_loc(i)
 enddo
END_PROVIDER 

BEGIN_PROVIDER [double precision, guess_loc_mos, (ao_num, n_orb_loc)]
 implicit none
 integer :: i,j,k
 j = 0

 do i  = 1, n_orb_ligand_loc
  j+=1 
  do k = 1, ao_num
   guess_loc_mos(k,j) = guess_loc_ligand_mos(k,i)
  enddo
 enddo

 do i  = 1, n_orb_metal_atoms_loc
  j+=1 
  do k = 1, ao_num
   guess_loc_mos(k,j) = guess_loc_metal_mos(k,i)
  enddo
 enddo



END_PROVIDER 
