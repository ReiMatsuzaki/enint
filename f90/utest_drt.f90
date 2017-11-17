program main
  use Mod_DRT
  implicit none
  type(Obj_DRT) :: drt
  call DRT_new(drt, 4,0,4)
  call DRT_dump(drt)
 
end program main
