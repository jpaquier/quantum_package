

double precision function compute_extra_diag_two_body_dm_ab(r1,r2) 
 implicit none
 BEGIN_DOC
! compute the extra diagonal contribution to the alpha/bet two body density at r1, r2
 END_DOC 
 double precision :: r1(3), r2(3)
 double precision :: compute_extra_diag_two_body_dm_ab_act,compute_extra_diag_two_body_dm_ab_core_act
 compute_extra_diag_two_body_dm_ab = compute_extra_diag_two_body_dm_ab_act(r1,r2)+compute_extra_diag_two_body_dm_ab_core_act(r1,r2)
end

double precision function compute_extra_diag_two_body_dm_ab_act(r1,r2)
 implicit none
 BEGIN_DOC
! compute the extra diagonal contribution to the two body density at r1, r2
! involving ONLY THE ACTIVE PART, which means that the four index of the excitations 
! involved in the two body density matrix are ACTIVE 
 END_DOC 
 PROVIDE n_act_orb
 double precision, intent(in) :: r1(3),r2(3)
 integer :: i,j,k,l
 double precision :: mos_array_r1(n_act_orb),mos_array_r2(n_act_orb)
 double precision :: contrib
 double precision :: contrib_tmp
!print*,'n_act_orb = ',n_act_orb
 compute_extra_diag_two_body_dm_ab_act = 0.d0
 call give_all_act_mos_at_r(r1,mos_array_r1)
 call give_all_act_mos_at_r(r2,mos_array_r2)
 do l = 1, n_act_orb  ! p2 
  do k = 1, n_act_orb  ! h2 
   do j = 1, n_act_orb  ! p1 
    do i = 1,n_act_orb   ! h1 
     contrib_tmp = mos_array_r1(i) * mos_array_r1(j) * mos_array_r2(k) * mos_array_r2(l)
     compute_extra_diag_two_body_dm_ab_act += two_body_dm_ab_big_array_act(i,j,k,l,1) * contrib_tmp
    enddo
   enddo
  enddo
 enddo

end

double precision function compute_extra_diag_two_body_dm_ab_core_act(r1,r2)
 implicit none
 BEGIN_DOC
! compute the extra diagonal contribution to the two body density at r1, r2
! involving ONLY THE ACTIVE PART, which means that the four index of the excitations 
! involved in the two body density matrix are ACTIVE 
 END_DOC 
 double precision, intent(in) :: r1(3),r2(3)
 integer :: i,j,k,l
 double precision :: mos_array_act_r1(n_act_orb),mos_array_act_r2(n_act_orb)
 double precision :: mos_array_core_r1(n_core_orb),mos_array_core_r2(n_core_orb)
 double precision :: contrib_core_1,contrib_core_2
 double precision :: contrib_act_1,contrib_act_2
 double precision :: contrib_tmp
 compute_extra_diag_two_body_dm_ab_core_act = 0.d0
 call give_all_act_mos_at_r(r1,mos_array_act_r1)
 call give_all_act_mos_at_r(r2,mos_array_act_r2)
 call give_all_core_mos_at_r(r1,mos_array_core_r1)
 call give_all_core_mos_at_r(r2,mos_array_core_r2)
  do i = 1, n_act_orb  ! h1 
   do j = 1, n_act_orb  ! p1 
    contrib_act_1 = mos_array_act_r1(i) * mos_array_act_r1(j) 
    contrib_act_2 = mos_array_act_r2(i) * mos_array_act_r2(j) 
    do k = 1,n_core_orb  ! h2 
     contrib_core_1 = mos_array_core_r1(k)  * mos_array_core_r1(k)
     contrib_core_2 = mos_array_core_r2(k)  * mos_array_core_r2(k)
     contrib_tmp = 0.5d0 * (contrib_act_1 * contrib_core_2 + contrib_act_2 * contrib_core_1)
     compute_extra_diag_two_body_dm_ab_core_act += two_body_dm_ab_big_array_core_act(k,i,j,1) * contrib_tmp
    enddo
   enddo
  enddo

end

double precision function compute_diag_two_body_dm_ab_core(r1,r2)
 implicit none
 double precision :: r1(3),r2(3)
 integer :: i,j,k,l
 double precision :: mos_array_r1(n_core_orb_allocate),mos_array_r2(n_core_orb_allocate)
 double precision :: contrib,contrib_tmp
 compute_diag_two_body_dm_ab_core = 0.d0
 call give_all_core_mos_at_r(r1,mos_array_r1)
 call give_all_core_mos_at_r(r2,mos_array_r2)
 do l = 1, n_core_orb  ! 
  contrib = mos_array_r2(l)*mos_array_r2(l)
! if(dabs(contrib).lt.threshld_two_bod_dm)cycle
  do k = 1, n_core_orb  ! 
    contrib_tmp = contrib * mos_array_r1(k)*mos_array_r1(k)
!  if(dabs(contrib).lt.threshld_two_bod_dm)cycle
     compute_diag_two_body_dm_ab_core += two_body_dm_ab_diag_core(k,l) * contrib_tmp
  enddo
 enddo

end


double precision function compute_diag_two_body_dm_ab_act(r1,r2)
 implicit none
 double precision :: r1(3),r2(3)
 integer :: i,j,k,l
 double precision :: mos_array_r1(n_act_orb),mos_array_r2(n_act_orb)
 double precision :: contrib,contrib_tmp
 compute_diag_two_body_dm_ab_act = 0.d0
 call give_all_act_mos_at_r(r1,mos_array_r1)
 call give_all_act_mos_at_r(r2,mos_array_r2)
 do l = 1, n_act_orb  ! 
  contrib = mos_array_r2(l)*mos_array_r2(l)
! if(dabs(contrib).lt.threshld_two_bod_dm)cycle
  do k = 1, n_act_orb  ! 
    contrib_tmp = contrib * mos_array_r1(k)*mos_array_r1(k)
!  if(dabs(contrib).lt.threshld_two_bod_dm)cycle
     compute_diag_two_body_dm_ab_act += two_body_dm_ab_diag_act(k,l,1) * contrib_tmp
  enddo
 enddo
end

double precision function compute_diag_two_body_dm_ab_core_act(r1,r2)
 implicit none
 double precision :: r1(3),r2(3)
 integer :: i,j,k,l
 double precision :: mos_array_core_r1(n_core_orb_allocate),mos_array_core_r2(n_core_orb_allocate)
 double precision :: mos_array_act_r1(n_act_orb),mos_array_act_r2(n_act_orb)
 double precision :: contrib_core_1,contrib_core_2
 double precision :: contrib_act_1,contrib_act_2
 double precision :: contrib_tmp
 compute_diag_two_body_dm_ab_core_act = 0.d0
 call give_all_act_mos_at_r(r1,mos_array_act_r1)
 call give_all_act_mos_at_r(r2,mos_array_act_r2)
 call give_all_core_mos_at_r(r1,mos_array_core_r1)
 call give_all_core_mos_at_r(r2,mos_array_core_r2)
! if(dabs(contrib).lt.threshld_two_bod_dm)cycle
 do k = 1, n_act_orb  ! 
    contrib_act_1 = mos_array_act_r1(k) * mos_array_act_r1(k)
    contrib_act_2 = mos_array_act_r2(k) * mos_array_act_r2(k)
    contrib_tmp = 0.5d0 * (contrib_act_1 * contrib_act_2 + contrib_act_2 * contrib_act_1)
!  if(dabs(contrib).lt.threshld_two_bod_dm)cycle
    do l = 1, n_core_orb  ! 
     contrib_core_1 = mos_array_core_r1(l) * mos_array_core_r1(l)
     contrib_core_2 = mos_array_core_r2(l) * mos_array_core_r2(l)
     compute_diag_two_body_dm_ab_core_act += two_body_dm_diag_core_act(l,k,1) * contrib_tmp
    enddo
 enddo
end

double precision function compute_diag_two_body_dm_ab(r1,r2)
 implicit none
 double precision,intent(in) :: r1(3),r2(3)
 double precision :: compute_diag_two_body_dm_ab_act,compute_diag_two_body_dm_ab_core
 double precision :: compute_diag_two_body_dm_ab_core_act
 compute_diag_two_body_dm_ab  =   compute_diag_two_body_dm_ab_act(r1,r2)+compute_diag_two_body_dm_ab_core(r1,r2) & 
                                + compute_diag_two_body_dm_ab_core_act(r1,r2)
end
