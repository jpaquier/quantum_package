
  use map_module
BEGIN_PROVIDER [integer, n_couple_ao_ao ]
 implicit none
 n_couple_ao_ao = ao_num * (ao_num  + 1) /2 
END_PROVIDER 

 BEGIN_PROVIDER [integer, index_ao_ao, (ao_num, ao_num)]
&BEGIN_PROVIDER [integer, index_ao_ao_reverse, (n_couple_ao_ao,2)]
 implicit none
 integer :: i,j,n_couple
 n_couple = 0
 do i = 1, ao_num
  do j = i, ao_num
   n_couple +=1 
   index_ao_ao(i,j) = n_couple
   index_ao_ao(j,i) = n_couple
   index_ao_ao_reverse(n_couple,1) = i
   index_ao_ao_reverse(n_couple,2) = j
  enddo
 enddo
 if(n_couple .ne.n_couple_ao_ao)then
  print*, 'PB !!! n_couple .ne.n_couple_ao_ao' 
 endif
 
END_PROVIDER 

BEGIN_PROVIDER [integer, n_couple_core_inact_virt]
 implicit none
 n_couple_core_inact_virt = n_core_inact_orb * n_virt_orb
END_PROVIDER 

 BEGIN_PROVIDER [integer, index_core_inact_virt, (n_core_orb, n_virt_orb)]
&BEGIN_PROVIDER [integer, index_core_inact_virt_reverse, (n_couple_core_inact_virt,2)]
 implicit none
 integer :: i,j,n_couple
 n_couple = 0
 do i = 1, n_core_orb
  do j = 1, n_virt_orb 
   n_couple +=1 
   index_core_inact_virt(i,j) = n_couple
   index_core_inact_virt(j,i) = n_couple
   index_ao_ao_reverse(n_couple,1) = i
   index_ao_ao_reverse(n_couple,2) = j
  enddo
 enddo
 if(n_couple .ne.n_couple_core_inact_virt)then
  print*, 'PB !!! n_couple .ne.n_couple_core_inact_virt' 
 endif
 
END_PROVIDER 


 BEGIN_PROVIDER [real(integral_kind), semi_transformed_occ_virt, (ao_num,ao_num,n_core_inact_orb,n_virt_orb)]
&BEGIN_PROVIDER [real(integral_kind), semi_transformed_occ_virt_bis, (ao_num,ao_num,n_core_inact_orb,n_virt_orb)]
  use map_module
 implicit none
 semi_transformed_occ_virt = 0.d0

 semi_transformed_occ_virt_bis = 0.d0

 integer :: i,j,k,l,iorb,jorb,korb,lorb,m,n
 double precision :: c
 real(integral_kind), allocatable :: bielec_tmp_0(:,:),matrix_tmp_1(:,:),matrix_final(:,:)
 real(integral_kind) :: integral,ao_bielec_integral

 allocate(bielec_tmp_0(ao_num,ao_num),matrix_tmp_1(n_core_inact_orb,ao_num),matrix_final(n_core_inact_orb,n_virt_orb))

 do k = 1,ao_num
   do l = 1,ao_num
    matrix_tmp_1 = 0.d0
    matrix_final = 0.d0
    do m = 1, ao_num
     call get_ao_bielec_integrals(k,m,l,ao_num,bielec_tmp_0(1,m)) ! k,l :: r1, m :: r2
    enddo
    call dgemm('N','N',n_core_inact_orb,ao_num,ao_num,1.d0,mo_coef_core_inact_transp,n_core_inact_orb,bielec_tmp_0,ao_num,0.d0,matrix_tmp_1,n_core_inact_orb)
    call dgemm('N','N',n_core_inact_orb,n_virt_orb,ao_num,1.d0,matrix_tmp_1,n_core_inact_orb,mo_coef_virt,ao_num,0.d0,matrix_final,n_core_inact_orb)
    do j = 1, n_virt_orb
     do i = 1, n_core_inact_orb
      semi_transformed_occ_virt(k,l,i,j) = matrix_final(i,j)
     enddo
    enddo
   enddo
  enddo
 
END_PROVIDER 

 BEGIN_PROVIDER [double precision, mo_coef_core_inact, (ao_num, n_core_inact_orb)]
&BEGIN_PROVIDER [double precision, mo_coef_core_inact_transp, (n_core_inact_orb,ao_num)]
 implicit none
 integer :: i,j,k,iorb,jorb
 do i = 1, n_core_inact_orb
  iorb = list_core_inact(i) 
  do k = 1, ao_num
   mo_coef_core_inact(k,i) = mo_coef(k,iorb)
   mo_coef_core_inact_transp(i,k) = mo_coef(k,iorb)
  enddo
 enddo

END_PROVIDER 

 BEGIN_PROVIDER [double precision, mo_coef_virt, (ao_num, n_virt_orb)]
&BEGIN_PROVIDER [double precision, mo_coef_virt_transp, (n_virt_orb,ao_num)]
 implicit none
 integer :: i,j,k,iorb,jorb
 do i = 1, n_virt_orb
  iorb = list_virt(i) 
  do k = 1, ao_num
   mo_coef_virt(k,i) = mo_coef(k,iorb)
   mo_coef_virt_transp(i,k) = mo_coef(k,iorb)
  enddo
 enddo

END_PROVIDER 

subroutine get_all_core_inact_virt_integrals(iorb,vorb,matrix_integrals)
  use map_module
 implicit none
 integer, intent(in) :: iorb,vorb
 real(integral_kind), intent(out) :: matrix_integrals(n_core_inact_orb,n_virt_orb)
 real(integral_kind) :: integrals_matrix(ao_num,ao_num)
 real(integral_kind) :: matrix_tmp_1(n_core_inact_orb,ao_num)
 integer :: i,j
 integer :: m,n

 do m = 1, ao_num
  do n = 1, ao_num
   integrals_matrix(m,n) = semi_transformed_occ_virt(m,n,iorb,vorb)
  enddo
 enddo
 matrix_integrals = 0.d0
 
 call dgemm('N','N',n_core_inact_orb,ao_num,ao_num,1.d0,mo_coef_core_inact_transp,n_core_inact_orb,integrals_matrix,ao_num,0.d0,matrix_tmp_1,n_core_inact_orb)
 call dgemm('N','N',n_core_inact_orb,n_virt_orb,ao_num,1.d0,matrix_tmp_1,n_core_inact_orb,mo_coef_virt,ao_num,0.d0,matrix_integrals,n_core_inact_orb)

end
