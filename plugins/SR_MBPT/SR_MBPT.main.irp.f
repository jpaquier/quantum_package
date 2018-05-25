program SR_MBPT
 implicit none
 read_wf = .True.
 touch read_wf
 call routine

end

subroutine routine
  use bitmasks
  implicit none
  BEGIN_DOC
! TODO
  END_DOC

  integer :: i,j,k,l
  integer          :: exc(0:2,2,2)
  integer          :: degree
  double precision :: phase,hij,hjj
  integer          :: h1,h2,p1,p2,s1,s2
  integer          :: hole,particl
  provide Fock_matrix_diag_mo
  double precision, allocatable :: psi_coef_pt(:,:), e_corr(:), delta_e0(:), delta_v(:),h0(:)
  double precision :: v1,href,vref,h00,e_corr_total
  allocate(psi_coef_pt(N_det, n_order_sr_pt),e_corr(n_order_sr_pt),delta_e0(N_det), delta_v(N_Det),h0(N_det))
  psi_coef_pt = 0.d0
  delta_e0 = 0.d0 
  delta_v = 0.d0
  e_corr = 0.d0
  e_corr_total = 0.d0
  call i_H_j(psi_det(1,1,1),psi_det(1,1,1),N_int,h00)
  call get_excitation_degree(psi_det(1,1,1),psi_det(1,1,2),degree,N_int)
  call get_excitation(psi_det(1,1,1),psi_det(1,1,2),exc,degree,phase,N_int)
  call decode_exc(exc,degree,hole,particl,h2,p2,s1,s2)
  psi_coef_pt(1,1) = 1.d0  ! RHF coef
  vref = 0.d0
  call bitstring_to_list( psi_det(1,1,1), list, n_elements, N_int)
  do k = 1, n_elements
   vref += Fock_matrix_diag_mo(list(k))
  enddo
  call bitstring_to_list( psi_det(1,2,1), list, n_elements, N_int)
  do k = 1, n_elements
   vref += Fock_matrix_diag_mo(list(k))
  enddo
  call i_H_j(psi_det(1,1,1),psi_det(1,1,1),N_int,href)
  v1 = href - vref  
  do j = 2, N_det
   if(sr_mbpt_h0==0)then ! Molle Plesset H0
     ! Diagonal matrix informations 
     call get_excitation_degree(psi_det(1,1,1),psi_det(1,1,j),degree,N_int)
     call get_excitation(psi_det(1,1,1),psi_det(1,1,j),exc,degree,phase,N_int)
     call decode_exc(exc,degree,h1,p1,h2,p2,s1,s2)
     if(degree==1)then
      delta_e0(j) = Fock_matrix_diag_mo(h1) - Fock_matrix_diag_mo(p1)
     else if (degree==2)then
      delta_e0(j) = Fock_matrix_diag_mo(h1) - Fock_matrix_diag_mo(p1) + Fock_matrix_diag_mo(h2) - Fock_matrix_diag_mo(p2)
     endif
     integer  :: list(N_int*bit_kind_size),n_elements
     call bitstring_to_list( psi_det(1,1,j), list, n_elements, N_int)
     h0(j) = 0.d0
     do k = 1, n_elements
      h0(j) += Fock_matrix_diag_mo(list(k))
     enddo
     call bitstring_to_list( psi_det(1,2,j), list, n_elements, N_int)
     do k = 1, n_elements
      h0(j) += Fock_matrix_diag_mo(list(k))
     enddo
     call i_H_j(psi_det(1,1,j),psi_det(1,1,j),N_int,hjj)
     delta_v(j) = hjj - (h0(j) + v1)
     
   else ! Epstein Nesbet H0
     call i_H_j(psi_det(1,1,j),psi_det(1,1,j),N_int,hjj)
     h0(j) = hjj
     delta_e0(j) = h00 - h0(j)
     delta_v(j) = 0.d0
   endif

     call i_H_j(psi_det(1,1,1),psi_det(1,1,j),N_int,hij)
     psi_coef_pt(j,1) = hij / delta_e0(j)
  enddo
  do i = 2, n_order_sr_pt
   do j = 2, N_det
    call i_H_j(psi_det(1,1,1),psi_det(1,1,j),N_int,hij)
    if(dabs(hij).gt.1.d0)then
     print*, hij,j,i
    endif
    e_corr(i) += psi_coef_pt(j,i-1) * hij
   enddo
   e_corr_total += e_corr(i)
   write(*,'(A16,I3,A5,F16.10)') ' partial e_corr(',i,')  = ',e_corr(i)
   write(*,'(A16,I3,A5,F16.10)') ' total   e_corr(',i,')  = ',e_corr_total
   do j = 2, N_det
    ! RENORMALIZATION TERM
    do l = 1, i-1
     psi_coef_pt(j,i) += - psi_coef_pt(j,i-l) * e_corr(l) / delta_e0(j)
    enddo
    ! DIAGONAL POTENTIAL TERM 
    psi_coef_pt(j,i) += psi_coef_pt(j,i-1) * delta_v(j) / delta_e0(j)
    
    ! INTERACTION TERM 
    do k = j+1, N_det
     call i_H_j(psi_det(1,1,k),psi_det(1,1,j),N_int,hij)
     if(verbose_large_interactions)then
      if(dabs(hij).gt.0.1d0)then
       print*, 'Large interaction found '
       print*, hij,delta_e0(j),psi_coef_pt(k,i-1)
       call get_excitation_degree(psi_det(1,1,j),psi_det(1,1,k),degree,N_int)
       call get_excitation(psi_det(1,1,j),psi_det(1,1,k),exc,degree,phase,N_int)
       call decode_exc(exc,degree,h1,p1,h2,p2,s1,s2)
       print*, 'degree = ',degree
       if(degree==2)then
        print*, h1,p1,h2,p2
        print*, s1,s2
       else 
        print*, h1,p1
        print*, s1
       endif
       double precision :: get_mo_bielec_integral,hole_integral,particl_integral
       hole_integral = get_mo_bielec_integral(h1,hole,p1,hole,mo_integrals_map)
       particl_integral = get_mo_bielec_integral(h1,particl,p1,particl,mo_integrals_map)
       print*, hole_integral,particl_integral,particl_integral - hole_integral
      endif
     endif
     psi_coef_pt(j,i) += hij / delta_e0(j) * psi_coef_pt(k,i-1)
     psi_coef_pt(k,i) += hij / delta_e0(k) * psi_coef_pt(j,i-1)
!    if(dabs(hij).gt.1.d0)then
!     print*, hij,j,k 
!    endif
    enddo
   enddo
   if(verbose_large_interactions)then
    print*, 'psi_coef ',psi_coef_pt(2,i) 
   endif
  enddo
  print*, 'TOTAL ECORR        =',e_corr_total 
  print*, 'VARIATIONAL ECORR  =',psi_energy(1) -  h00
  print*, 'TOTAL ENERGY       =',e_corr_total + href  + nuclear_repulsion
  print*, 'VARIATIONAL ENERGY =',psi_energy(1) + nuclear_repulsion
  print*, 'ERROR TO VAR CALC  =',psi_energy(1) -  h00 - e_corr_total 
end
