 ! 1 
 BEGIN_PROVIDER [ complex*16, dirac_ao_bi_elec_integral_L_alpha_L_alpha, (ao_num, ao_num) ]
  implicit none
  BEGIN_DOC
  ! L_alpha L_alpha bloc of the bi-electronic Fock matrix in dirac AO basis set
   END_DOC
  integer                        :: i,j,k,l 
  complex*16                     :: D
  double precision               :: dirac_ao_bielec_integral   
  dirac_ao_bi_elec_integral_L_alpha_L_alpha = (0.d0,0.d0)
  do i = 1, ao_num
   do j = 1, ao_num
    do k = 1, 2*dirac_ao_num
     do l = 1, 2*dirac_ao_num
      D = dirac_SCF_density_matrix_ao(k,l)     
      if (k .le. ao_num .and. l .le. ao_num) then
       dirac_ao_bi_elec_integral_L_alpha_L_alpha(i,j) += D*(dirac_ao_bielec_integral(i,j,d_L(k),d_L(l))- dirac_ao_bielec_integral(i,d_L(l),d_L(k),j))
      elseif ((k .gt. ao_num .and. k .le. (2*ao_num) .and. l.gt. ao_num .and. l.le. (2*ao_num)) .or. &
              (k .gt. (2*ao_num) .and. k .le. (2*ao_num+small_ao_num) .and. l .gt. (2*ao_num) .and. l .le. (2*ao_num+small_ao_num)) .or. &
              (k .gt. (2*ao_num+small_ao_num) .and. l.gt.(2*ao_num+small_ao_num))) then
       dirac_ao_bi_elec_integral_L_alpha_L_alpha(i,j) += D*dirac_ao_bielec_integral(i,j,d_L(k),d_L(l))
      endif
     enddo
    enddo
   enddo
  enddo
 END_PROVIDER 

 !2
 BEGIN_PROVIDER [ complex*16, dirac_ao_bi_elec_integral_L_beta_L_beta, (ao_num, ao_num) ]
  implicit none
  BEGIN_DOC
  ! L_beta L_beta bloc of the bi-electronic Fock matrix in dirac AO basis set
   END_DOC
  integer                        :: i,j,k,l 
  complex*16                     :: D
  double precision               :: dirac_ao_bielec_integral   
  dirac_ao_bi_elec_integral_L_beta_L_beta = (0.d0,0.d0)
  do i = 1, ao_num
   do j = 1, ao_num
    do k = 1, 2*dirac_ao_num
     do l = 1, 2*dirac_ao_num
      D = dirac_SCF_density_matrix_ao(k,l)
      if (k .gt. ao_num .and. k .le. (2*ao_num) .and. l.gt. ao_num .and. l.le. (2*ao_num)) then
       dirac_ao_bi_elec_integral_L_beta_L_beta(i,j) += D*(dirac_ao_bielec_integral(i,j,d_L(k),d_L(l)) - dirac_ao_bielec_integral(i,d_L(l),d_L(k),j))
      elseif ((k .le. ao_num .and. l .le. ao_num) .or. &
              (k > (2*ao_num) .and. k <= (2*ao_num+small_ao_num) .and. l> (2*ao_num) .and. l <= (2*ao_num+small_ao_num)) .or. &
              (k .gt. (2*ao_num+small_ao_num) .and. l.gt.(2*ao_num+small_ao_num))) then
       dirac_ao_bi_elec_integral_L_beta_L_beta(i,j) += D*dirac_ao_bielec_integral(i,j,d_L(k),d_L(l))
      endif
     enddo
    enddo
   enddo
  enddo
 END_PROVIDER 

 !3
 BEGIN_PROVIDER [ complex*16, dirac_ao_bi_elec_integral_L_beta_L_alpha, (ao_num, ao_num) ]
  implicit none
  BEGIN_DOC
  ! L_beta L_alpha bloc of the bi-electronic Fock matrix in dirac AO basis set
   END_DOC
  integer                        :: i,j,k,l 
  complex*16                     :: D
  double precision               :: dirac_ao_bielec_integral   
  dirac_ao_bi_elec_integral_L_beta_L_alpha = (0.d0,0.d0)
  do i = 1, ao_num
   do j = 1, ao_num
    do k = 1, 2*dirac_ao_num
     do l = 1, 2*dirac_ao_num
      if (k .le. ao_num .and. l .gt. ao_num .and. l .le. 2*ao_num) then
       D = dirac_SCF_density_matrix_ao(k,l)
       dirac_ao_bi_elec_integral_L_beta_L_alpha(i,j) += D*(- dirac_ao_bielec_integral(i,d_L(l),d_L(k),j))
      endif
     enddo
    enddo
   enddo
  enddo
 END_PROVIDER 

 !4
 BEGIN_PROVIDER [ complex*16, dirac_ao_bi_elec_integral_L_alpha_L_beta, (ao_num, ao_num) ]
  implicit none
  BEGIN_DOC
  ! L_beta L_alpha bloc of the bi-electronic Fock matrix in dirac AO basis set
   END_DOC
  integer                        :: i,j,k,l 
  complex*16                     :: D
  double precision               :: dirac_ao_bielec_integral   
  dirac_ao_bi_elec_integral_L_alpha_L_beta = (0.d0,0.d0)
  do i = 1, ao_num
   do j = 1, ao_num
    dirac_ao_bi_elec_integral_L_alpha_L_beta(i,j) = Conjg(dirac_ao_bi_elec_integral_L_beta_L_alpha(j,i))
   !do k = 1, 2*dirac_ao_num
   ! do l = 1, 2*dirac_ao_num
   !  D = dirac_SCF_density_matrix_ao(k,l)
   !   if (k .gt. ao_num .and. k .le. 2*ao_num .and. l .le. ao_num) then
   !    dirac_ao_bi_elec_integral_L_beta_L_alpha(i,j) += D*(- dirac_ao_bielec_integral(i,d_L(l),d_L(k),j))
   !  endif
   ! enddo
   !enddo
   enddo
  enddo
 END_PROVIDER 

 !5
 BEGIN_PROVIDER [ complex*16, dirac_ao_bi_elec_integral_S_alpha_S_alpha, (small_ao_num, small_ao_num) ]
  implicit none
  BEGIN_DOC
  ! S_alpha S_alpha bloc of the bi-electronic Fock matrix in dirac AO basis set
   END_DOC
  integer                        :: i,i_plus,j,j_plus,k,l 
  complex*16                     :: D
  double precision               :: dirac_ao_bielec_integral   
  dirac_ao_bi_elec_integral_S_alpha_S_alpha = (0.d0,0.d0)
  do i = 1, small_ao_num
   i_plus = i + large_ao_num
   do j = 1, small_ao_num
    j_plus = j + large_ao_num
    do k = 1, 2*dirac_ao_num
     do l = 1, 2*dirac_ao_num
      D = dirac_SCF_density_matrix_ao(k,l)
      if (k .gt. (2*large_ao_num) .and. k .le. (2*large_ao_num+small_ao_num) .and. l .gt. (2*large_ao_num) .and. l .le. (2*large_ao_num+small_ao_num)) then
       dirac_ao_bi_elec_integral_S_alpha_S_alpha(i,j) += D*(dirac_ao_bielec_integral(i_plus,j_plus,d_L(k),d_L(l)) - dirac_ao_bielec_integral(i_plus,d_L(l),d_L(k),j_plus))
      elseif ((k .le. ao_num .and. l .le. ao_num) .or. &
              (k .gt. ao_num .and. k .le. (2*ao_num) .and. l.gt. ao_num .and. l.le. (2*ao_num)) .or. &
              (k .gt. (2*ao_num+small_ao_num) .and. l.gt.(2*ao_num+small_ao_num))) then
       dirac_ao_bi_elec_integral_S_alpha_S_alpha(i,j) += D*dirac_ao_bielec_integral(i_plus,j_plus,d_L(k),d_L(l))
      endif
     enddo
    enddo
   enddo
  enddo
 END_PROVIDER 

 !6
 BEGIN_PROVIDER [ complex*16, dirac_ao_bi_elec_integral_S_beta_S_beta, (small_ao_num, small_ao_num) ]
  implicit none
  BEGIN_DOC
  ! S_beta S_beta bloc of the bi-electronic Fock matrix in dirac AO basis set
   END_DOC
  integer                        :: i,i_plus,j,j_plus,k,l 
  complex*16                     :: D
  double precision               :: dirac_ao_bielec_integral   
  dirac_ao_bi_elec_integral_S_beta_S_beta = (0.d0,0.d0)
  do i = 1, small_ao_num
   i_plus = i + large_ao_num
   do j = 1, small_ao_num 
    j_plus = j + large_ao_num
    do k = 1, 2*dirac_ao_num
     do l = 1, 2*dirac_ao_num
      D = dirac_SCF_density_matrix_ao(k,l)
      if (k .gt. (2*large_ao_num+small_ao_num) .and. l .gt. (2*large_ao_num+small_ao_num)) then
       dirac_ao_bi_elec_integral_S_beta_S_beta(i,j) += D*(dirac_ao_bielec_integral(i_plus,j_plus,d_L(k),d_L(l)) - dirac_ao_bielec_integral(i_plus,d_L(l),d_L(k),j_plus))
      elseif ((k .le. ao_num .and. l .le. ao_num) .or. &
              (k .gt. ao_num .and. k .le. (2*ao_num) .and. l.gt. ao_num .and. l.le. (2*ao_num)) .or. &
              (k .gt. (2*ao_num) .and. k .le. (2*ao_num+small_ao_num) .and. l .gt. (2*ao_num) .and. l .le. (2*ao_num+small_ao_num))) then
       dirac_ao_bi_elec_integral_S_beta_S_beta(i,j) += D*dirac_ao_bielec_integral(i_plus,j_plus,d_L(k),d_L(l))
      endif
     enddo
    enddo
   enddo
  enddo
 END_PROVIDER 

 !7
 BEGIN_PROVIDER [ complex*16, dirac_ao_bi_elec_integral_S_beta_S_alpha, (small_ao_num, small_ao_num) ]
  implicit none
  BEGIN_DOC
  ! S_beta S_alpha bloc of the bi-electronic Fock matrix in dirac AO basis set
   END_DOC
  integer                        :: i,i_plus,j,j_plus,k,l 
  complex*16                     :: D
  double precision               :: dirac_ao_bielec_integral   
  dirac_ao_bi_elec_integral_S_beta_S_alpha = (0.d0,0.d0)
  do i = 1, small_ao_num
   i_plus = i + large_ao_num
   do j = 1, small_ao_num
   j_plus = j + large_ao_num
    do k = 1, 2*dirac_ao_num
     do l = 1, 2*dirac_ao_num
      if (k .gt. (2*ao_num) .and. k .le. (2*ao_num+small_ao_num) .and. l .gt. (2*ao_num+small_ao_num)) then
       D = dirac_SCF_density_matrix_ao(k,l)
       dirac_ao_bi_elec_integral_S_beta_S_alpha(i,j) += D*(- dirac_ao_bielec_integral(i_plus,d_L(l),d_L(k),j_plus))
      endif
     enddo
    enddo
   enddo
  enddo
 END_PROVIDER 

 !8
 BEGIN_PROVIDER [ complex*16, dirac_ao_bi_elec_integral_S_alpha_S_beta, (small_ao_num, small_ao_num) ]
  implicit none
  BEGIN_DOC
  ! S_alpha S_beta bloc of the bi-electronic Fock matrix in dirac AO basis set
   END_DOC
  integer                        :: i,i_plus,j,j_plus,k,l 
  complex*16                     :: D
  double precision               :: dirac_ao_bielec_integral   
  dirac_ao_bi_elec_integral_S_alpha_S_beta = (0.d0,0.d0)
  do i = 1, small_ao_num 
   i_plus = i + large_ao_num
   do j = 1, small_ao_num
   j_plus = j + large_ao_num
    dirac_ao_bi_elec_integral_S_alpha_S_beta(i,j) = Conjg(dirac_ao_bi_elec_integral_S_beta_S_alpha(j,i))
   !do k = 1, 2*dirac_ao_num
   ! do l = 1, 2*dirac_ao_num
   !  D = dirac_SCF_density_matrix_ao(k,l)
   !   if (k .gt. (2*ao_num+small_ao_num) .and. l .gt. (2*ao_num) .and. l .le. (2*ao_num+small_ao_num)) then
   !    dirac_ao_bi_elec_integral_S_beta_S_alpha(i,j) += D*(- dirac_ao_bielec_integral(i_plus,d_L(l),d_L(k),j_plus))
   !  endif
   ! enddo
   !enddo
   enddo
  enddo
 END_PROVIDER 

 !9
 BEGIN_PROVIDER [ complex*16, dirac_ao_bi_elec_integral_S_alpha_L_alpha, (small_ao_num, large_ao_num) ]
  implicit none
  BEGIN_DOC
  ! S_alpha L_alpha bloc of the bi-electronic Fock matrix in dirac AO basis set
   END_DOC
  integer                        :: i,i_plus,j,k,l 
  complex*16                     :: D
  double precision               :: dirac_ao_bielec_integral   
  dirac_ao_bi_elec_integral_S_alpha_L_alpha = (0.d0,0.d0)
  do i = 1, small_ao_num
   i_plus = i + large_ao_num
   do j = 1, ao_num
    do k = 1, 2*dirac_ao_num
     do l = 1, 2*dirac_ao_num
      if (k.le. ao_num .and. l .gt. 2*ao_num .and. l .le. (2*ao_num+small_ao_num)) then
       D = dirac_SCF_density_matrix_ao(k,l)
       dirac_ao_bi_elec_integral_S_alpha_L_alpha(i,j) += D*(- dirac_ao_bielec_integral(i_plus,d_L(l),d_L(k),j))
      endif
     enddo
    enddo
   enddo
  enddo
 END_PROVIDER

 !10
 BEGIN_PROVIDER [ complex*16, dirac_ao_bi_elec_integral_L_alpha_S_alpha, (large_ao_num, small_ao_num) ]
  implicit none
  BEGIN_DOC
  ! L_alpha S_alpha bloc of the bi-electronic Fock matrix in dirac AO basis set
   END_DOC
  integer                        :: i,j,j_plus,k,l 
  complex*16                     :: D
  double precision               :: dirac_ao_bielec_integral   
  dirac_ao_bi_elec_integral_L_alpha_S_alpha = (0.d0,0.d0)
  do i = 1, ao_num
   do j = 1, small_ao_num
   j_plus = j + large_ao_num
    dirac_ao_bi_elec_integral_L_alpha_S_alpha(i,j) = Conjg(dirac_ao_bi_elec_integral_S_alpha_L_alpha(j,i))
   !do k = 1, 2*dirac_ao_num
   ! do l = 1, 2*dirac_ao_num
   !  if (k .gt. 2*ao_num .and. k .le. (2*ao_num+small_ao_num) .and. l .le. ao_num) then
   !   D = dirac_SCF_density_matrix_ao(k,l)
   !   dirac_ao_bi_elec_integral_S_alpha_L_alpha(i,j) += D*(- dirac_ao_bielec_integral(i,d_L(l),d_L(k),j_plus))
   !  endif
   ! enddo
   !enddo
   enddo
  enddo
 END_PROVIDER
 
 !11
 BEGIN_PROVIDER [ complex*16, dirac_ao_bi_elec_integral_S_beta_L_beta, (small_ao_num, large_ao_num) ]
  implicit none
  BEGIN_DOC
  ! S_beta L_beta bloc of the bi-electronic Fock matrix in dirac AO basis set
   END_DOC
  integer                        :: i,i_plus,j,k,l 
  complex*16                     :: D
  double precision               :: dirac_ao_bielec_integral   
  dirac_ao_bi_elec_integral_S_beta_L_beta = (0.d0,0.d0)
  do i = 1, small_ao_num
   i_plus = i + large_ao_num
   do j = 1, ao_num
    do k = 1, 2*dirac_ao_num
     do l = 1, 2*dirac_ao_num
      if (k .gt. ao_num .and. k .le. 2*ao_num .and. l .gt. (2*ao_num+small_ao_num)) then
       D = dirac_SCF_density_matrix_ao(k,l)
       dirac_ao_bi_elec_integral_S_alpha_L_alpha(i,j) += D*(- dirac_ao_bielec_integral(i_plus,d_L(l),d_L(k),j))
      endif
     enddo
    enddo
   enddo
  enddo
 END_PROVIDER

 !12
 BEGIN_PROVIDER [ complex*16, dirac_ao_bi_elec_integral_L_beta_S_beta, (large_ao_num, small_ao_num) ]
  implicit none
  BEGIN_DOC
  ! L_beta S_beta bloc of the bi-electronic Fock matrix in dirac AO basis set
   END_DOC
  integer                        :: i,j,j_plus,k,l 
  complex*16                     :: D
  double precision               :: dirac_ao_bielec_integral   
  dirac_ao_bi_elec_integral_L_beta_S_beta = (0.d0,0.d0)
  do i = 1, ao_num
   do j = 1, small_ao_num
   j_plus = j + large_ao_num
    dirac_ao_bi_elec_integral_L_beta_S_beta(i,j) = Conjg(dirac_ao_bi_elec_integral_S_beta_L_beta(j,i))
   !do k = 1, 2*dirac_ao_num
   ! do l = 1, 2*dirac_ao_num
   !  if (k .gt. (2*ao_num+small_ao_num) .and. l .gt. ao_num .and. l .le. 2*ao_num) then
   !   D = dirac_SCF_density_matrix_ao(k,l)
   !   dirac_ao_bi_elec_integral_S_alpha_L_alpha(i,j) += D*(- dirac_ao_bielec_integral(i,d_L(l),d_L(k),j_plus))
   !  endif
   ! enddo
   !enddo
   enddo
  enddo
 END_PROVIDER

 !13
 BEGIN_PROVIDER [ complex*16, dirac_ao_bi_elec_integral_S_alpha_L_beta, (small_ao_num, large_ao_num) ]
  implicit none
  BEGIN_DOC
  ! S_alpha L_beta bloc of the bi-electronic Fock matrix in dirac AO basis set
   END_DOC
  integer                        :: i,i_plus,j,k,l 
  complex*16                     :: D
  double precision               :: dirac_ao_bielec_integral   
  dirac_ao_bi_elec_integral_S_alpha_L_beta = (0.d0,0.d0)
  do i = 1, small_ao_num
   i_plus = i + large_ao_num
   do j = 1, ao_num
    do k = 1, 2*dirac_ao_num
     do l = 1, 2*dirac_ao_num
      if (k .gt. ao_num .and. k .le. 2*ao_num .and. l .gt. (2*ao_num) .and. l .le. (2*ao_num+small_ao_num)) then
       D = dirac_SCF_density_matrix_ao(k,l)
       dirac_ao_bi_elec_integral_S_alpha_L_alpha(i,j) += D*(- dirac_ao_bielec_integral(i_plus,d_L(l),d_L(k),j))
      endif
     enddo
    enddo
   enddo
  enddo
 END_PROVIDER

 !14
 BEGIN_PROVIDER [ complex*16, dirac_ao_bi_elec_integral_L_beta_S_alpha, (large_ao_num, small_ao_num) ]
  implicit none
  BEGIN_DOC
  ! L_beta S_alpha bloc of the bi-electronic Fock matrix in dirac AO basis set
   END_DOC
  integer                        :: i,j,j_plus,k,l 
  complex*16                     :: D
  double precision               :: dirac_ao_bielec_integral   
  dirac_ao_bi_elec_integral_L_beta_S_alpha = (0.d0,0.d0)
  do i = 1, ao_num
   do j = 1, small_ao_num
   j_plus = j + large_ao_num
    dirac_ao_bi_elec_integral_L_beta_S_alpha(i,j) = Conjg(dirac_ao_bi_elec_integral_S_alpha_L_beta(j,i))
   !do k = 1, 2*dirac_ao_num
   ! do l = 1, 2*dirac_ao_num
   !  if (k .gt. (2*ao_num) .and. k .le. (2*ao_num+small_ao_num) .and. l .gt. ao_num .and. l .le. 2*ao_num) then
   !   D = dirac_SCF_density_matrix_ao(k,l)
   !   dirac_ao_bi_elec_integral_S_alpha_L_alpha(i,j) += D*(- dirac_ao_bielec_integral(i,d_L(l),d_L(k),j_plus))
   !  endif
   ! enddo
   !enddo
   enddo
  enddo
 END_PROVIDER

 !15
 BEGIN_PROVIDER [ complex*16, dirac_ao_bi_elec_integral_S_beta_L_alpha, (small_ao_num, large_ao_num) ]
  implicit none
  BEGIN_DOC
  ! S_beta L_alpha bloc of the bi-electronic Fock matrix in dirac AO basis set
   END_DOC
  integer                        :: i,i_plus,j,k,l 
  complex*16                     :: D
  double precision               :: dirac_ao_bielec_integral   
  dirac_ao_bi_elec_integral_S_beta_L_alpha = (0.d0,0.d0)
  do i = 1, small_ao_num
   i_plus = i + large_ao_num
   do j = 1, ao_num
    do k = 1, 2*dirac_ao_num
     do l = 1, 2*dirac_ao_num
      if (k .le. ao_num .and. l .gt. (2*ao_num+small_ao_num)) then
       D = dirac_SCF_density_matrix_ao(k,l)
       dirac_ao_bi_elec_integral_S_alpha_L_alpha(i,j) += D*(- dirac_ao_bielec_integral(i_plus,d_L(l),d_L(k),j))
      endif
     enddo
    enddo
   enddo
  enddo
 END_PROVIDER

 !16
 BEGIN_PROVIDER [ complex*16, dirac_ao_bi_elec_integral_L_alpha_S_beta, (large_ao_num, small_ao_num) ]
  implicit none
  BEGIN_DOC
  ! L_alpha S_beta bloc of the bi-electronic Fock matrix in dirac AO basis set
   END_DOC
  integer                        :: i,j,j_plus,k,l 
  complex*16                     :: D
  double precision               :: dirac_ao_bielec_integral   
  dirac_ao_bi_elec_integral_L_alpha_S_beta = (0.d0,0.d0)
  do i = 1, ao_num
   do j = 1, small_ao_num
   j_plus = j + large_ao_num
    dirac_ao_bi_elec_integral_L_alpha_S_beta(i,j) = Conjg(dirac_ao_bi_elec_integral_S_beta_L_alpha(j,i))
   !do k = 1, 2*dirac_ao_num
   ! do l = 1, 2*dirac_ao_num
   !  if (k .gt. (2*ao_num+small_ao_num) .and. l .le. ao_num) then
   !   D = dirac_SCF_density_matrix_ao(k,l)
   !   dirac_ao_bi_elec_integral_S_alpha_L_alpha(i,j) += D*(- dirac_ao_bielec_integral(i,d_L(l),d_L(k),j_plus))
   !  endif
   ! enddo
   !enddo
   enddo
  enddo
 END_PROVIDER



