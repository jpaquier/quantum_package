
BEGIN_PROVIDER [integer, type_of_superci]
implicit none
 BEGIN_DOC 
 ! if 0 then use the simplest FOCK matrix for the SUPERCI , if 1 use the integrals with three index to correct it, if 2 use the (ia|jb) integrals also 
 END_DOC  
type_of_superci = 2
END_PROVIDER 

BEGIN_PROVIDER [logical, state_following]
implicit none
 BEGIN_DOC
! if .True. then performs a state following during the superci CASSCF optimization
 END_DOC 


END_PROVIDER 
