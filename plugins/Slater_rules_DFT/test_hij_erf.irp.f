program test_hij
 implicit none
 read_wf = .true.
 touch read_wf 
 call routine_2
end

subroutine routine
!!!!!!!!! see slater_rules_general.irp.f for the hij with H_mu
 implicit none 
 integer :: i
 integer :: degree
 double precision :: hij
 integer          :: exc(0:2,2,2)
 double precision :: phase,hii_ref,hii
 integer :: h1,p1,h2,p2,s1,s2
 double precision :: get_mo_bielec_integral
 do i = 1, N_det
  print*,''
  print*,'i = ',i
  call debug_det(psi_det(1,1,i),N_int)
  call get_excitation_degree(psi_det(1,1,i),psi_det(1,1,1),degree,N_int)
  print*,'degree = ',degree
  if(degree == 0)then
   print*,'Reference determinant '
   call i_H_j(psi_det(1,1,i),psi_det(1,1,i),N_int,hii_ref)
  else 
   call i_H_j(psi_det(1,1,i),psi_det(1,1,i),N_int,hii)
   call i_H_j(psi_det(1,1,1),psi_det(1,1,i),N_int,hij)
   call get_excitation(psi_det(1,1,1),psi_det(1,1,i),exc,degree,phase,N_int)
   call decode_exc(exc,degree,h1,p1,h2,p2,s1,s2)
   print*,'phase = ',phase
   print*, 'hij',hij
   if(degree == 1)then
    print*,'s1',s1
    print*,'h1,p1 = ',h1,p1
    double precision :: hmono,hdouble
    call  i_H_j_verbose(psi_det(1,1,1),psi_det(1,1,i),N_int,hij,hmono,hdouble)
   write(*,'(A17,X,F16.12)')'hmono          = ',hmono
   write(*,'(A17,X,F16.12)')'hdouble        = ',hdouble
   write(*,'(A17,X,F16.12)')'hmono+hdouble  = ',(hmono+hdouble) * phase
   else
    print*,'s1',s1
    print*,'h1,p1 = ',h1,p1
    print*,'s2',s2
    print*,'h2,p2 = ',h2,p2
    write(*,'(A17,X,F16.12)'), 'hij            = ',get_mo_bielec_integral(h1,h2,p1,p2,mo_integrals_map)
    if(s1==s2)then
     write(*,'(A17,X,F16.12)'), 'exchange term  = ',get_mo_bielec_integral(h1,h2,p2,p1,mo_integrals_map)
     write(*,'(A17,X,F16.12)'), 'total          = ',phase * (get_mo_bielec_integral(h1,h2,p1,p2,mo_integrals_map) - get_mo_bielec_integral(h1,h2,p2,p1,mo_integrals_map))
    endif
   endif
   
   write(*,'(A17,X,F16.12)') ' <Ref| H |D_I> = ',hij
   write(*,'(A17,X,F16.12)') ' Delta E (Ref) = ',-(hii - hii_ref)
   write(*,'(A17,X,F16.12)')' amplitude (1) = ',hij/(hii_ref - hii)
  endif
   write(*,'(A17,X,F16.12)')' amplitude     = ',psi_coef(i,1)/psi_coef(1,1)

 enddo


 print*,''
 print*,''
 print*,''

end

subroutine routine_2
 implicit none
 double precision :: accu
 integer :: i
 accu = 0.d0
 do i = 1, N_det
  accu += psi_coef(i,1)**2
 enddo
 print*,'accu = ',accu

end
