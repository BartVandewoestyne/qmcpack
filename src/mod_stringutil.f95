module mod_stringutil

!  implicit none

  public :: trimlr

  contains

  function trimlr(mystring) result(outstring)
    character(len=*), intent(in) :: mystring
    character(len=len(mystring)) :: outstring

    outstring = trim(adjustl(mystring))

  end function trimlr

end module mod_stringutil
