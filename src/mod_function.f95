module mod_function

  use numeric_kinds

  private

  public :: allocate_function_parameters

  type, public :: functionparams
    real(kind=qp), dimension(:), pointer   :: a
    real(kind=qp), dimension(:), pointer   :: b
    real(kind=qp), dimension(:), pointer   :: c
    real(kind=qp), dimension(:), pointer   :: u
    real(kind=qp), dimension(:,:), pointer :: x
    integer(kind=i4b), dimension(:), pointer :: k
  end type functionparams

  contains

    subroutine allocate_function_parameters(params, s)
      type(functionparams), intent(inout) :: params
      integer(kind=i4b), intent(in)       :: s
      
      allocate(params%a(s))
      allocate(params%c(s))
      allocate(params%k(s))
      allocate(params%u(s))

    end subroutine allocate_function_parameters

end module mod_function
