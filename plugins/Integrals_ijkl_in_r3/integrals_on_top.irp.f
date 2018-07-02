program test
double precision :: inte_ana,inte_num,inte_mo_bour
!print*,'coucou'
!print*,ao_integrals_threshold
!print*,'coucou'
 call test_inte_num_bart
 call test_inte_mo_bart
end

subroutine test_inte_mo_bart
 implicit none
 double precision:: inte_mo_bour
 double precision :: inte_mo_bg, inte_mo_bg_loc, inte_mo_bour_loc
 integer :: m,n,p,q,i,j,k,l
 double precision :: integral_bourrin_mo, get_mo_bielec_integral_ijkl_r3 
 inte_mo_bour_loc = 0.d0
 inte_mo_bg_loc = 0.d0
 inte_mo_bour = 0.d0
 inte_mo_bg = 0.d0

 do l = 1, mo_tot_num
  do k = 1, mo_tot_num
   do j = 1, mo_tot_num
    do i = 1, mo_tot_num
     inte_mo_bour_loc = integral_bourrin_mo(i,j,k,l)
     inte_mo_bg_loc = get_mo_bielec_integral_ijkl_r3(i,j,k,l,mo_integrals_ijkl_r3_map)
     inte_mo_bour += inte_mo_bour_loc
     inte_mo_bg += inte_mo_bg_loc 
     if(dabs(inte_mo_bour_loc - inte_mo_bg_loc)/dabs(inte_mo_bour_loc + inte_mo_bg_loc).gt.1.d-6.and.dabs(inte_mo_bg_loc).gt.1.d-10)then
      print*,'Bourrin = ',inte_mo_bour_loc
      print*,'On map  = ',inte_mo_bg_loc
      print*,dabs(inte_mo_bour_loc - inte_mo_bg_loc)/dabs(inte_mo_bour_loc + inte_mo_bg_loc)
     endif  
    enddo
   enddo
  enddo
 enddo
 print*,' '
 print*,' '
 print*,'********************'
 print*,'Bourrin = ',inte_mo_bour
 print*,'On map  = ',inte_mo_bg

end


subroutine test_inte_ana_bart(inte_ana)
 implicit none
 double precision, intent(out) :: inte_ana 
 integer :: m,n,p,q
 double precision :: intete,ao_bielec_integral_ijkl_r3 
 inte_ana = 0 
 do m = 1, ao_num
  do n = 1, ao_num 
   do p = 1, ao_num 
    do q = 1, ao_num 
     inte_ana += ao_bielec_integral_ijkl_r3(q,p,n,m)
     print*,'Analytical integral  = ',inte_ana        
    enddo
   enddo
  enddo
 enddo
end


subroutine test_inte_num_bart
 implicit none
 double precision:: inte_num
 integer :: i,j,k,l,m,n,p,q
 double precision :: intete
 double precision :: r(3), ao_bielec_integral_ijkl_r3
 double precision,allocatable :: aos_array(:)
 double precision :: ao_bielec_integral_ijkl_r3_test, integrals_in_map
 double precision :: get_ao_bielec_integral_ijkl_r3
 allocate(aos_array(ao_num))
 
  r(1) = 0.d0 
  r(2) = 0.d0
  r(3) = 0.d0
!print*,ao_integrals_threshold
 do m= 1, ao_num
  do n= 1, ao_num
   do p= 1, ao_num
    do q= 1, ao_num
!     intete = 0.d0
!     do j = 1, nucl_num
!      do k = 1, n_points_radial_grid  -1
!       do l = 1, n_points_integration_angular
!        r(1) = grid_points_per_atom(1,l,k,j)
!        r(2) = grid_points_per_atom(2,l,k,j)
!        r(3) = grid_points_per_atom(3,l,k,j)
!        call give_all_aos_at_r(r,aos_array)
!        intete += aos_array(q)*aos_array(p)*aos_array(n)*aos_array(m) * final_weight_functions_at_grid_points(l,k,j) 
!       enddo
!      enddo
!     enddo
       inte_num = ao_bielec_integral_ijkl_r3(q,p,n,m)
       integrals_in_map = get_ao_bielec_integral_ijkl_r3(q,p,n,m,ao_integrals_ijkl_r3_map)
       if(dabs(inte_num - integrals_in_map)/dabs(inte_num + integrals_in_map).gt.1.d-6.and.dabs(integrals_in_map).gt.1.d-10)then
        print*,'analytical = ',inte_num
        print*,'numerical  = ',integrals_in_map
        print*,dabs(inte_num - integrals_in_map)/dabs(inte_num + integrals_in_map)
       endif
     enddo 
    enddo
   enddo
  enddo

end



