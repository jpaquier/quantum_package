program foboscf
 implicit none
!if(disk_access_ao_integrals == "None" .or. disk_access_ao_integrals == "Read" )then
! disk_access_ao_integrals = "Write"
! touch disk_access_ao_integrals
!endif
!print*, 'disk_access_ao_integrals',disk_access_ao_integrals
 no_oa_or_av_opt = .True.
 touch no_oa_or_av_opt
 call run_prepare
 call routine_fobo_scf
!call save_mos

end

subroutine run_prepare
 implicit none
 call opt_orb
end

subroutine routine_fobo_scf
 implicit none
 integer :: i,j
 double precision :: norm_total(N_States)
 print*,''
 print*,''
 character*(64) :: label
 label = "Natural"
 do i = 1, 5
  print*,'*******************************************************************************'
  print*,'*******************************************************************************'
  print*,'FOBO-SCF Iteration ',i
  print*, 'ao_bielec_integrals_in_map = ',ao_bielec_integrals_in_map
  print*,'*******************************************************************************'
  print*,'*******************************************************************************'
  if(speed_up_convergence_foboscf)then
   if(i==3)then
    threshold_lmct = max(threshold_lmct,0.001)
    threshold_mlct = max(threshold_mlct,0.05)
    soft_touch threshold_lmct threshold_mlct
   endif
   if(i==4)then
    threshold_lmct = max(threshold_lmct,0.005)
    threshold_mlct = max(threshold_mlct,0.07)
    soft_touch threshold_lmct threshold_mlct
   endif
   if(i==5)then
    threshold_lmct = max(threshold_lmct,0.01)
    threshold_mlct = max(threshold_mlct,0.1)
    soft_touch threshold_lmct threshold_mlct
   endif
  endif
  call initialize_mo_coef_begin_iteration
  call print_mos(i)
  call FOBOCI_lmct_mlct_old_thr(i,norm_total)
  call save_osoci_natural_mos(norm_total)
  call reorder_active_orb
! touch mo_coef
  call opt_orb

  call clear_mo_map
  call provide_properties
  call save_mos 
 enddo



end
