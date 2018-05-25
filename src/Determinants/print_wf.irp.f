program printwf
 implicit none 
 read_wf = .True.
 touch read_wf
 print*,'ref_bitmask_energy = ',ref_bitmask_energy
 call routine
!call print_cas

end

subroutine routine
 implicit none 
 integer :: i
 integer :: degree
 double precision :: hij
 integer          :: exc(0:2,2,2)
 double precision :: phase,hii_ref,hii
 integer :: h1,p1,h2,p2,s1,s2
 double precision :: get_mo_bielec_integral
 double precision :: norm_mono_a,norm_mono_b
 norm_mono_a = 0.d0
 norm_mono_b = 0.d0
 do i = 1, N_det
!do i = 1, min(500,N_det)
  print*,''
  print*,'i = ',i
  call debug_det(psi_det(1,1,i),N_int)
  call get_excitation_degree(psi_det(1,1,i),psi_det(1,1,1),degree,N_int)
  print*,'degree = ',degree
  if(degree == 0)then
   print*,'Reference determinant '
   call i_H_j(psi_det(1,1,i),psi_det(1,1,i),N_int,hii_ref)
  else 
   integer :: nh,np,number_of_holes,number_of_particles
   nh = number_of_holes(psi_det(1,1,i))
   np = number_of_particles(psi_det(1,1,i))
   print*, 'nh,np',nh,np
   call i_H_j(psi_det(1,1,i),psi_det(1,1,i),N_int,hii)
   call i_H_j(psi_det(1,1,1),psi_det(1,1,i),N_int,hij)
   call get_excitation(psi_det(1,1,1),psi_det(1,1,i),exc,degree,phase,N_int)
   call decode_exc(exc,degree,h1,p1,h2,p2,s1,s2)
   print*,'phase = ',phase
   print*, 'hij',hij
   if(degree == 1)then
    print*,'s1',s1
    print*,'h1,p1 = ',h1,p1
    if(s1 == 1)then
     norm_mono_a += dabs(psi_coef(i,1)/psi_coef(1,1))
    else
     norm_mono_b += dabs(psi_coef(i,1)/psi_coef(1,1))
    endif
   double precision :: hmono,hdouble
   call  i_H_j_verbose(psi_det(1,1,1),psi_det(1,1,i),N_int,hij,hmono,hdouble)
   write(*,'(A17,X,F16.12)')'hmono          = ',hmono
   write(*,'(A17,X,F16.12)')'hdouble        = ',hdouble
   write(*,'(A17,X,F16.12)')'hmono+hdouble  = ',hmono+hdouble
   write(*,'(A17,X,F16.12)')'hij            = ',hij
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
   write(*,'(A17,X,F16.12)') ' Delta E (Ref) = ',hii - hii_ref
   write(*,'(A17,X,F16.12)')' amplitude (1) = ',hij/(hii_ref - hii)
  endif
   write(*,'(A17,X,F16.12)')' amplitude     = ',psi_coef(i,1)/psi_coef(1,1)

 enddo


 print*,''
 print*,''
 print*,''
 print*,'mono alpha = ',norm_mono_a
 print*,'mono beta  = ',norm_mono_b

end

subroutine print_cas
 implicit none
 integer :: i,j
 do i = 1, N_det_cas
  call debug_det(psi_cas(1,1,i),N_int)
 enddo

end
