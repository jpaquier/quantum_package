program pouet
 read_wf = .True.
 touch read_wf
 call routine3
end

subroutine routine3
 implicit none
 integer :: i,j,k,l
 double precision :: accu
 accu = 0.d0
 print*, 'energy_x      =',energy_x
 print*, 'energy_c      =',energy_c
 print*, 'energy_c_md   =',energy_c_md
end
