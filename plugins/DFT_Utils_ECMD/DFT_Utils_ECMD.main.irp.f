program DFT_Utils_ECMD
  implicit none
  BEGIN_DOC
! TODO
  END_DOC
  read_wf = .True.
  touch read_wf
  call routine
end

subroutine routine
  implicit none
! print*,'integral_on_top    = ',integral_on_top
  print*,'Energy_c_md_on_top = ',Energy_c_md_on_top
! print*,'Energy_c_md_LDA    = ',Energy_c_md_LDA
! print*,'Energy_c_md_on_top_PBE_mu_UEG_vector = ',Energy_c_md_on_top_PBE_mu_UEG_vector
! print*,'Energy_c_md_on_top_PBE_mu_vector     = ',Energy_c_md_on_top_PBE_mu_vector
end
