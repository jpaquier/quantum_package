
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
  use map_module
 implicit none
 semi_transformed_occ_virt = 0.d0
 integer :: i,j,k,l,iorb,jorb,korb,lorb,m
 integer :: i1,j1,k1,l1,kmax
 double precision :: c
 real(integral_kind), allocatable :: bielec_tmp_0(:,:)
 integer, allocatable           :: bielec_tmp_0_idx(:)
 double precision, allocatable :: mo_coef_one_column(:),mo_coef_one_column_bis(:),contraction_tmp(:)
 double precision :: u_dot_v
 double precision :: tmp(ao_num)
 allocate(bielec_tmp_0(ao_num,ao_num),bielec_tmp_0_idx(ao_num),mo_coef_one_column(ao_num_align),contraction_tmp(ao_num_align),mo_coef_one_column_bis(ao_num_align))

 do l1 = 1,ao_num
   !DEC$ VECTOR ALIGNED
   do k1 = 1,ao_num
     !DEC$ VECTOR ALIGNED
     do j1 = 1,ao_num
       call get_ao_bielec_integrals(k1,j1,l1,ao_num,bielec_tmp_0(1,j1)) ! all integrals for a given l1, k1
     enddo
     do i = 1, n_core_inact_orb
      iorb = list_core_inact(i)
      do k = 1, ao_num
       mo_coef_one_column(k) = mo_coef(k,iorb)
       contraction_tmp(k) = 0.d0
      enddo
!     tmp = 0.d0
!     do m = 1, ao_num
!      do k = 1, ao_num
!       tmp(m) += bielec_tmp_0(m,k) * mo_coef_one_column(k)
!      enddo
!     enddo
      call matrix_vector_product(mo_coef_one_column,contraction_tmp,bielec_tmp_0,ao_num,ao_num)
!     do m = 1, ao_num
!      print*, m,'m'
!      if(dabs(tmp(m) - contraction_tmp(m)).gt.1.d-10)then
!       print*, tmp(m),contraction_tmp(m),dabs(tmp(m) - contraction_tmp(m))
!      endif
!     enddo
      
      do j = 1, n_virt_orb 
       jorb = list_virt(i)
       double precision :: accu 
       accu = 0.d0
       do k = 1, ao_num
        mo_coef_one_column_bis(k) = mo_coef(k,jorb)
        accu += mo_coef(k,jorb) * contraction_tmp(k)
       enddo
       
       semi_transformed_occ_virt(l1,k1,i,j) = u_dot_v(mo_coef_one_column_bis,contraction_tmp,ao_num)
       if(dabs(accu - semi_transformed_occ_virt(l1,k1,i,j)).gt.1.d-10)then
        print*, accu,semi_transformed_occ_virt(l1,k1,i,j),dabs(accu - semi_transformed_occ_virt(l1,k1,i,j))
       endif
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
 integer :: iorb_i,vorb_v
 real(integral_kind), intent(out) :: matrix_integrals(n_core_inact_orb,n_virt_orb)
 real(integral_kind) :: integrals_matrix(ao_num,ao_num)
 real(integral_kind) :: integrals_tmp(n_core_inact_orb,ao_num)
 integer :: i,j,k,l
!iorb_i = list_core_inact_reverse(iorb_i)
!vorb_v = list_virt_reverse(vorb_v)
 do i = 1, ao_num
  do j = 1, ao_num
   integrals_matrix(i,j) = semi_transformed_occ_virt(j,i,iorb,vorb)
  enddo
 enddo
 
 ! integrals_tmp = C integrals_matrix  
 call dgemm( 'N','N',n_core_inact_orb,ao_num,ao_num, 1.d0,mo_coef_core_inact_transp,n_core_inact_orb,integrals_matrix,ao_num,1.d0,integrals_tmp,n_core_inact_orb)
 ! matrix_integrals = integrals_tmp C^+
 call dgemm( 'N','N',n_core_inact_orb,n_virt_orb,ao_num, 1.d0,integrals_tmp,n_core_inact_orb,mo_coef_virt,ao_num,1.d0,matrix_integrals,n_virt_orb)

end
