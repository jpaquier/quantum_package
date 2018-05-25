logical function is_the_hole_in_det(key_in,ispin,i_hole)
  use bitmasks
  ! returns true if the electron ispin is absent from i_hole
 implicit none
 integer, intent(in) :: i_hole,ispin
 integer(bit_kind), intent(in) :: key_in(N_int,2)
 integer(bit_kind) :: key_tmp(N_int)
 integer(bit_kind) :: itest(N_int)
 integer :: i,j,k
 do i = 1, N_int
  itest(i) = 0_bit_kind
 enddo
 k = ishft(i_hole-1,-bit_kind_shift)+1
 j = i_hole-ishft(k-1,bit_kind_shift)-1
 itest(k) = ibset(itest(k),j)
 j = 0
 do i = 1, N_int
  key_tmp(i) = iand(itest(i),key_in(i,ispin))
  j += popcnt(key_tmp(i))
 enddo
 if(j==0)then
  is_the_hole_in_det = .True.
 else
  is_the_hole_in_det = .False.
 endif

end

logical function is_the_particl_in_det(key_in,ispin,i_particl)
  use bitmasks
  ! returns true if the electron ispin is absent from i_particl
 implicit none
 integer, intent(in) :: i_particl,ispin
 integer(bit_kind), intent(in) :: key_in(N_int,2)
 integer(bit_kind) :: key_tmp(N_int)
 integer(bit_kind) :: itest(N_int)
 integer :: i,j,k
 do i = 1, N_int
  itest(i) = 0_bit_kind
 enddo
 k = ishft(i_particl-1,-bit_kind_shift)+1
 j = i_particl-ishft(k-1,bit_kind_shift)-1
 itest(k) = ibset(itest(k),j)
 j = 0
 do i = 1, N_int
  key_tmp(i) = iand(itest(i),key_in(i,ispin))
  j += popcnt(key_tmp(i))
 enddo
 if(j==0)then
  is_the_particl_in_det = .False.
 else
  is_the_particl_in_det = .True.
 endif

end

subroutine find_hole_in_det(key_in,i_hole)
 use bitmasks
 implicit none
 integer(bit_kind), intent(in) :: key_in(N_int,2)
 integer, intent(out) :: i_hole(0:mo_tot_num,2)
 integer :: i
 integer(bit_kind) :: key_tmp(N_int,2)
 integer :: occ(N_int*bit_kind_size), itest
 do i = 1, N_int
  key_tmp(i,1) = xor(reunion_of_core_inact_bitmask(i,1),iand(reunion_of_core_inact_bitmask(i,1),key_in(i,1)))
  key_tmp(i,2) = xor(reunion_of_core_inact_bitmask(i,1),iand(reunion_of_core_inact_bitmask(i,2),key_in(i,2)))
 enddo
 
 call bitstring_to_list(key_tmp(1,1), occ, itest, N_int)
 i_hole(0,1) = itest
 do i = 1, itest
  i_hole(i,1) = occ(i)
 enddo

 call bitstring_to_list(key_tmp(1,2), occ, itest, N_int)
 i_hole(0,2) = itest
 do i = 1, itest
  i_hole(i,2) = occ(i)
 enddo

end

subroutine find_particle_in_det(key_in,i_part)
 use bitmasks
 implicit none
 integer(bit_kind), intent(in) :: key_in(N_int,2)
 integer, intent(out) :: i_part(0:mo_tot_num,2)
 integer :: i
 integer(bit_kind) :: key_tmp(N_int,2)
 integer :: occ(N_int*bit_kind_size), itest
 do i = 1, N_int
  key_tmp(i,1) = iand(virt_bitmask(i,1),key_in(i,1))
  key_tmp(i,2) = iand(virt_bitmask(i,2),key_in(i,2))
 enddo
 
 call bitstring_to_list(key_tmp(1,1), occ, itest, N_int)
 i_part(0,1) = itest
 do i = 1, itest
  i_part(i,1) = occ(i)
 enddo

 call bitstring_to_list(key_tmp(1,2), occ, itest, N_int)
 i_part(0,2) = itest
 do i = 1, itest
  i_part(i,2) = occ(i)
 enddo

end
