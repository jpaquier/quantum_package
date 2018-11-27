program print_good_orbs
 implicit none
 character(len=100):: str_final_core
 character(len=100):: str_final_act
 character(len=100):: str_final_del
!provide index_ligand_orb_loc
!call routine_core(str_final_core)
!call routine_act(str_final_act)
!print*,'core = ',trim(str_final_core)
!print*,'act  = ',trim(str_final_act)
 print*,'n_orb_ligand_virt_loc = ',n_orb_ligand_virt_loc
end

subroutine routine_core(final_str)
 character(len=100), intent(out) :: final_str
 use bitmasks
 implicit none
 integer :: i 
 integer, allocatable :: list_core_overlap(:)
 allocate(list_core_overlap(N_int * bit_kind_size))
 integer :: n_lmct,n_core_orb_tmp
 integer(bit_kind), allocatable :: key(:),key2(:)
 allocate(key(N_int),key2(N_int))
 key = 0_bit_kind
 do i = 1, n_orb_ligand_loc
  call set_bit_to_integer(index_ligand_orb_loc(i),key,N_int)
 enddo
 do i = 1, N_int
  key2(i) = xor(key(i),inact_bitmask(i,1)) 
 enddo
 call bitstring_to_list( key2, list_core_overlap, n_core_orb_tmp, N_int)
  print*,'i,list_core_overlap(i) '
 do i = 1, n_core_orb_tmp
  print*,i,list_core_overlap(i) 
 enddo
 integer :: n_start(mo_tot_num), n_end(mo_tot_num), n_couple
 n_couple = 1
 n_start(n_couple) = list_core_overlap(n_couple)
 n_end(n_couple) = index_ligand_orb_loc_sorted(n_couple) - 1
 do i = 1, n_orb_ligand_loc
  if(index_ligand_orb_loc_sorted(i)+1.lt.index_metal_atom_orb_loc(1))then
   n_couple +=1
   n_start(n_couple) = index_ligand_orb_loc_sorted(i)+1
   if(n_couple.le.n_orb_ligand_loc)then
    n_end(n_couple) = index_ligand_orb_loc_sorted(i+1)-1
   else 
    n_end(n_couple) = index_metal_atom_orb_loc(1)-1
   endif
  endif
 enddo
 print*,'n_couple = ',n_couple
 do i = 1, n_couple
  print*,i,n_start(i),n_end(i) 
 enddo

 character(len=20) :: str,i1,i2,test,print_core_begin,print_core_end
 character(len=20), allocatable :: print_core(:)
 allocate(print_core(n_couple))
 do i = 1, n_couple
  i1 =  trim(str(n_start(i)))
  i2 =  trim(str(n_end(i)))
  print_core(i)=trim(trim(i1)//'-'//trim(i2)//',')
 enddo
 i = 1
 i1 =  trim(str(n_start(i)))
 i2 =  trim(str(n_end(i)))
 print_core_begin=trim('"['//trim(i1)//'-'//trim(i2)//',')
 i = n_couple
 i1 =  trim(str(n_start(i)))
 i2 =  trim(str(n_end(i)))
 print_core_end=trim(trim(i1)//'-'//trim(i2)//']"')
 final_str = trim(print_core_begin)
 do i = 2, n_couple -1 
  final_str = trim(trim(final_str)//trim(print_core(i)))
 enddo
 final_str = trim(trim(final_str)//trim(print_core_end))
 integer:: nchar
 nchar = len(trim(final_str))
 print*,'nchar = ',nchar
 write(*,*)trim(trim(final_str))
end

subroutine routine_act(final_str)
 character(len=100), intent(out) :: final_str
 
 use bitmasks
 implicit none
 integer :: i 
 integer, allocatable :: list_act_overlap(:)
 allocate(list_act_overlap(N_int * bit_kind_size))
 integer :: n_lmct,n_act_orb_tmp
 integer(bit_kind), allocatable :: key(:),key2(:)
 allocate(key(N_int),key2(N_int))
 key = 0_bit_kind
 do i = 1, n_orb_ligand_loc
  call set_bit_to_integer(index_ligand_orb_loc(i),key,N_int)
 enddo
 do i = 1, N_int
  key2(i) = ior(key(i),act_bitmask(i,1)) 
 enddo
 call bitstring_to_list( key2, list_act_overlap, n_act_orb_tmp, N_int)

 character(len=20) :: str,i1,i2,test,print_act_begin,print_act_end
 character(len=20), allocatable :: print_act(:)
 allocate(print_act(n_act_orb_tmp))
 do i = 1, n_act_orb_tmp
  i1 =  trim(str(list_act_overlap(i)))
  print_act(i)=trim(trim(i1)//',')
 enddo
 i = 1
 i1 =  trim(str(list_act_overlap(i)))
 print_act_begin=trim('"['//trim(i1)//',')
 i = n_act_orb_tmp
 i1 =  trim(str(list_act_overlap(i)))
 print_act_end=trim(trim(i1)//']"')

 final_str = trim(print_act_begin)
 do i = 2, n_act_orb_tmp -1 
  final_str = trim(trim(final_str)//trim(print_act(i)))
 enddo
 final_str = trim(trim(final_str)//trim(print_act_end))
 integer:: nchar
 nchar = len(trim(final_str))
 print*,'nchar = ',nchar
 write(*,*)trim(trim(final_str))
end

subroutine routine_del(final_str)
 character(len=100), intent(out) :: final_str
 
 use bitmasks
 implicit none
 integer :: i 
 integer, allocatable :: list_act_overlap(:)
 allocate(list_act_overlap(N_int * bit_kind_size))
 integer :: n_lmct,n_act_orb_tmp
 integer(bit_kind), allocatable :: key(:),key2(:)
 allocate(key(N_int),key2(N_int))
 key = 0_bit_kind
 do i = 1, n_orb_ligand_loc
  call set_bit_to_integer(index_ligand_orb_loc(i),key,N_int)
 enddo
 do i = 1, N_int
  key2(i) = ior(key(i),act_bitmask(i,1)) 
 enddo
 call bitstring_to_list( key2, list_act_overlap, n_act_orb_tmp, N_int)

 character(len=20) :: str,i1,i2,test,print_act_begin,print_act_end
 character(len=20), allocatable :: print_act(:)
 allocate(print_act(n_act_orb_tmp))
 do i = 1, n_act_orb_tmp
  i1 =  trim(str(list_act_overlap(i)))
  print_act(i)=trim(trim(i1)//',')
 enddo
 i = 1
 i1 =  trim(str(list_act_overlap(i)))
 print_act_begin=trim('"['//trim(i1)//',')
 i = n_act_orb_tmp
 i1 =  trim(str(list_act_overlap(i)))
 print_act_end=trim(trim(i1)//']"')

 final_str = trim(print_act_begin)
 do i = 2, n_act_orb_tmp -1 
  final_str = trim(trim(final_str)//trim(print_act(i)))
 enddo
 final_str = trim(trim(final_str)//trim(print_act_end))
 integer:: nchar
 nchar = len(trim(final_str))
 print*,'nchar = ',nchar
 write(*,*)trim(trim(final_str))
end

