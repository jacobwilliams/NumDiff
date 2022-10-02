!*******************************************************************************
!> author: Jacob Williams
!  date: October 29, 2016
!
!  Numerical differentiation module for computing the Jacobian matrix
!  (the derivative matrix of `m` functions w.r.t. `n` variables) using
!  finite differences.

    module numerical_differentiation_module

    use numdiff_kinds_module
    use numdiff_utilities_module
    use iso_fortran_env,      only: error_unit
    use dsm_module,           only: dsm
    use diff_module,          only: diff_func
    use numdiff_cache_module, only: function_cache

    implicit none

    private

    real(wp),parameter :: zero = 0.0_wp

    type,public :: finite_diff_method

        !! defines the finite difference method
        !! used to compute the Jacobian.
        !!
        !! See: [[get_finite_difference_method]]
        !! for the different methods.

        private

        integer :: id = 0 !! unique ID for the method
        character(len=:),allocatable :: name  !! the name of the method
        integer :: class = 0 !! 2=backward diffs, 3=central diffs, etc...
        real(wp),dimension(:),allocatable :: dx_factors  !! multiplicative factors for `dx` perturbation
        real(wp),dimension(:),allocatable :: df_factors  !! multiplicative factors for accumulating function evaluations
        real(wp)                          :: df_den_factor = zero  !! denominator factor for finite difference equation (times `dx`)
    contains
        private
        procedure,public :: get_formula
        procedure,public :: print => print_finite_difference_method
    end type finite_diff_method
    interface finite_diff_method
        !! constructor
        module procedure initialize_finite_difference_method
    end interface

    type,public :: meth_array
        !! to store an array of [[finite_diff_method]] types
        !! this is used when the `mode=2` option is used
        !! in [[numdiff_type]]
        private
        type(finite_diff_method),dimension(:),allocatable :: meth
    end type meth_array

    type,public :: sparsity_pattern

        !! A sparsity pattern

        private

        logical :: sparsity_computed = .false. !! has the sparsity pattern already been computed?
        integer :: num_nonzero_elements = 0 !! number of nonzero elements in the jacobian
                                            !! (will be the dimension of `irow` and `icol`)
        integer,dimension(:),allocatable :: irow  !! sparsity pattern - rows of non-zero elements
        integer,dimension(:),allocatable :: icol  !! sparsity pattern - columns of non-zero elements

        logical :: linear_sparsity_computed = .false. !! if the linear pattern has been populated
        integer :: num_nonzero_linear_elements = 0 !! number of constant elements in the jacobian
        integer,dimension(:),allocatable  :: linear_irow  !! linear sparsity pattern - rows of non-zero elements
        integer,dimension(:),allocatable  :: linear_icol  !! linear sparsity pattern - columns of non-zero elements
        real(wp),dimension(:),allocatable :: linear_vals  !! linear elements of the jacobian

        integer,dimension(:),allocatable :: indices  !! index vector
                                                     !! `[1,2,...,num_nonzero_elements]`
                                                     !! for putting `df` into `jac`

        integer :: maxgrp = 0 !! the number of groups in the partition
                              !! of the columns of `a`.
        integer,dimension(:),allocatable :: ngrp !! specifies the partition of the columns of `a`.
                                                 !! column `jcol` belongs to group `ngrp(jcol)`.
                                                 !! `size(n)`

        contains
        private
        procedure :: dsm_wrapper
        procedure :: compute_indices
        procedure,public :: destroy => destroy_sparsity
        procedure,public :: print => print_sparsity
        procedure,public :: columns_in_partition_group
    end type sparsity_pattern

    type,public :: numdiff_type

        !! base type for sparsity and Jacobian computations.

        private

        logical :: exception_raised = .false. !! if true, an exception has been raised
        integer :: istat = 0   !! a non-zero value will cause all routine to exit.
                               !! this can be set to `-1` by calling [[terminate]].
        character(len=:),allocatable :: error_msg !! error message (if `istat/=0`)

        integer :: n = 0 !! number of `x` variables
        integer :: m = 0 !! number of `f` functions

        real(wp),dimension(:),allocatable :: xlow  !! lower bounds on `x`
        real(wp),dimension(:),allocatable :: xhigh !! upper bounds on `x`

        logical :: print_messages = .true. !! if true, warning messages are printed
                                           !! to the `error_unit` for any errors.

        integer :: chunk_size = 100  !! chuck size for allocating the arrays (>0)

        integer :: perturb_mode = 1  !! perturbation mode:
                                     !!
                                     !! **1** - perturbation is `dx=dpert`,
                                     !! **2** - perturbation is `dx=dpert*x`,
                                     !! **3** - perturbation is `dx=dpert*(1+x)`
        real(wp),dimension(:),allocatable :: dpert !! perturbation vector for `x`

        logical :: partition_sparsity_pattern = .false.  !! to partition the sparsity pattern using [[dsm]]
        type(sparsity_pattern) :: sparsity  !! the sparsity pattern
        real(wp),dimension(:),allocatable :: xlow_for_sparsity  !! lower bounds on `x` for computing sparsity (optional)
        real(wp),dimension(:),allocatable :: xhigh_for_sparsity !! upper bounds on `x` for computing sparsity (optional)
        real(wp),dimension(:),allocatable :: dpert_for_sparsity !! perturbation vector for `x` when computing
                                                                !! sparsity for `sparsity_mode=4`
        integer :: num_sparsity_points = 3    !! for `sparsity_mode=4`, the number of jacobian
                                              !! evaluations used to estimate the sparsity pattern.
                                              !! See [[compute_sparsity_random_2]]
        integer :: sparsity_perturb_mode = 1  !! perturbation mode (if `sparsity_mode=4`):
                                              !!
                                              !! **1** - perturbation is `dx=dpert`,
                                              !! **2** - perturbation is `dx=dpert*x`,
                                              !! **3** - perturbation is `dx=dpert*(1+x)`

        ! if computing the sparsity pattern, we also have an option to
        ! compute the linear pattern, which indicates constant elements
        ! of the jacobian. these elements don't need to be computed again.
        logical :: compute_linear_sparsity_pattern = .false.  !! to also compute the linear sparsity pattern
        real(wp) :: linear_sparsity_tol = epsilon(1.0_wp)    !! the equality tolerance for derivatives to
                                                             !! indicate a constant jacobian element (linear sparsity)
        real(wp) :: function_precision_tol = epsilon(1.0_wp) !! the function precision. two functions values
                                                             !! that are the within this tolerance are
                                                             !! considered the same value. This is used
                                                             !! when estimating the sparsity pattern when
                                                             !! `sparsity_mode=2` in [[compute_sparsity_random]]

        integer :: mode = 1 !! **1** = use `meth` (specified methods),
                            !! **2** = use `class` (specified class, method is selected on-the-fly).
        type(finite_diff_method),dimension(:),allocatable :: meth   !! the finite difference method to use
                                                                    !! compute the `n`th column of the Jacobian
                                                                    !! `size(n)`.  Either this or `class` is used
        integer,dimension(:),allocatable :: class  !! the class of method to use to
                                                   !! compute the `n`th column of the Jacobian
                                                   !! `size(n)`. Either this or `meth` is used
        type(meth_array),dimension(:),allocatable :: class_meths !! array of methods for the specified classes.
                                                                 !! used with `class` when `mode=2`

        ! parameters when using diff:
        real(wp) :: eps = 1.0e-9_wp !! tolerance parameter for [[diff]]
        real(wp) :: acc = 0.0_wp    !! tolerance parameter for [[diff]]

        type(function_cache) :: cache    !! if using the function cache

        logical :: use_diff = .false. !! if we are using the Neville's process method,
                                      !! rather than finite differences

        procedure(func),pointer :: problem_func => null()
            !! the user-defined function

        ! these are required to be defined by the user:
        procedure(func),pointer :: compute_function => null()
            !! compute the user-defined function
            !! this can point to the `problem_func` or, if using
            !! the cache, it points to [[compute_function_with_cache]].

        procedure(spars_f),pointer :: compute_sparsity => null()
            !! for computing the sparsity pattern

        procedure(info_f),pointer :: info_function => null()
            !! an optional function the user can define
            !! which is called when each column of the jacobian is computed.
            !! It can be used to perform any setup operations.

        procedure(jacobian_f),pointer :: jacobian_function => null()
            !! the low-level function called by [[compute_jacobian]]
            !! that actually computes the jacobian.

    contains

        private

        procedure,public :: initialize => initialize_numdiff !! initialize the class
        procedure,public :: diff_initialize => initialize_numdiff_for_diff !! initialize the class

        procedure,public :: compute_jacobian              !! main routine to compute the Jacobian
                                                          !! using the selected options. It
                                                          !! returns the sparse (vector) form.
        procedure,public :: compute_jacobian_dense        !! return the dense `size(m,n)`
                                                          !! matrix form of the Jacobian.
        procedure,public :: compute_jacobian_times_vector !! returns the product of the Jacobian
                                                          !! matrix and an input vector
        procedure,public :: destroy => destroy_numdiff_type  !! destroy the class
        procedure,public :: print_sparsity_pattern  !! print the sparsity pattern in vector form to a file
        procedure,public :: print_sparsity_matrix   !! print the sparsity pattern in matrix form to a file
        procedure,public :: set_sparsity_pattern    !! manually set the sparsity pattern
        procedure,public :: select_finite_diff_method  !! select a method in a specified class so
                                                       !! that the variable bounds are not violated
                                                       !! when by the perturbations.
        procedure,public :: select_finite_diff_method_for_partition_group  !! version of [[select_finite_diff_method]]
                                                                           !! for partitioned sparsity pattern.
        procedure,public :: set_numdiff_bounds  !! can be called to change the variable bounds.

        procedure,public :: compute_sparsity_pattern !! if the user needs to compute the sparsity pattern manually.
                                                     !! (otherwise it will be done the first time the Jacobian is computed)
        procedure,public :: get_sparsity_pattern     !! returns the sparsity pattern (if it is allocated)

        procedure,public :: terminate !! can be called by user to stop the computation

        procedure,public :: failed !! to check if an exception was raised.
        procedure,public :: get_error_status !! the status of error condition

        ! internal routines:
        procedure :: destroy_sparsity_pattern      !! destroy the sparsity pattern
        procedure :: compute_perturb_vector
        procedure :: compute_perturbation_vector   !! computes the variable perturbation factor
        procedure :: compute_sparsity_perturbation_vector
        procedure :: perturb_x_and_compute_f
        procedure :: perturb_x_and_compute_f_partitioned
        procedure :: set_numdiff_sparsity_bounds
        procedure :: set_sparsity_mode
        procedure :: generate_dense_sparsity_partition
        procedure :: compute_jacobian_for_sparsity
        procedure :: resize_sparsity_vectors
        procedure :: raise_exception
        procedure :: clear_exceptions

    end type numdiff_type

    abstract interface
        subroutine func(me,x,f,funcs_to_compute)
            !! The function (vector array of output functions `f`, computed
            !! from a vector of input variables `x`).
            !! This must be defined for all computations.
            import :: numdiff_type,wp
            implicit none
            class(numdiff_type),intent(inout) :: me
            real(wp),dimension(:),intent(in) :: x !! array of variables (size `n`)
            real(wp),dimension(:),intent(out) :: f !! array of functions (size `m`)
            integer,dimension(:),intent(in) :: funcs_to_compute !! the elements of the
                                                                !! function vector that need
                                                                !! to be computed (the other
                                                                !! are ignored)
        end subroutine func
        subroutine spars_f(me,x)
            !! The function to compute the sparsity pattern.
            !! It populates the `irow` and `icol` variables in the class.
            import :: numdiff_type,wp
            implicit none
            class(numdiff_type),intent(inout) :: me
            real(wp),dimension(:),intent(in) :: x !! vector of variables (size `n`)
        end subroutine spars_f
        subroutine info_f(me,column,i,x)
            !! User-defined info function (optional).
            !! Informs user what is being done during Jacobian computation.
            !! It can be used to perform any setup operations that need to
            !! done on the user's end.
            import :: numdiff_type,wp
            implicit none
            class(numdiff_type),intent(inout)  :: me
            integer,dimension(:),intent(in) :: column !! the columns being computed.
            integer,intent(in) :: i                   !! perturbing these columns for the `i`th time (1,2,...)
            real(wp),dimension(:),intent(in)  :: x    !! the nominal variable vector
        end subroutine info_f
        subroutine jacobian_f(me,x,dx,jac)
            !! Actual function for computing the Jacobian
            !! called by [[compute_jacobian]].
            import :: numdiff_type,wp
            implicit none
            class(numdiff_type),intent(inout)   :: me
            real(wp),dimension(:),intent(in)    :: x    !! vector of variables (size `n`)
            real(wp),dimension(me%n),intent(in) :: dx   !! absolute perturbation (>0) for each variable
            real(wp),dimension(:),intent(out)   :: jac  !! sparse jacobian vector (size `num_nonzero_elements`)
        end subroutine jacobian_f
    end interface

    ! other:
    public :: get_finite_diff_formula
    public :: get_all_methods_in_class

    contains
!*******************************************************************************

!*******************************************************************************
!>
!  Constructor for a [[finite_diff_method]].
!
!@note factors are input as integers for convenience, but are converted
!      to reals for the actual computations. (note: this means we can't
!      currently define methods that have non-integer factors).

    function initialize_finite_difference_method(id,name,class,dx_factors,&
                                                 df_factors,df_den_factor) result(me)

    implicit none

    type(finite_diff_method)        :: me
    integer,intent(in)              :: id            !! unique ID for the method
    character(len=*),intent(in)     :: name          !! the name of the method
    integer,intent(in)              :: class         !! 2=backward diffs, 3=central diffs, etc...
    integer,dimension(:),intent(in) :: dx_factors    !! multiplicative factors for dx perturbation
    integer,dimension(:),intent(in) :: df_factors    !! multiplicative factors for accumulating function evaluations
    integer,intent(in)              :: df_den_factor !! denominator factor for finite difference equation (times dx)

    if (size(dx_factors)/=size(df_factors)) then
        error stop 'Error in initialize_finite_difference_method: '//&
                   'dx_factors and df_factors arrays must be the same size.'
    else

        me%id     = id
        me%name   = trim(name)
        me%class  = class

        ! the following is not strictly necessary if the
        ! compiler fully supports auto-LHS allocation:
        allocate(me%dx_factors(size(dx_factors)))
        allocate(me%df_factors(size(df_factors)))

        me%dx_factors = real(dx_factors,wp)
        me%df_factors = real(df_factors,wp)

        me%df_den_factor = real(df_den_factor,wp)

    end if

    end function initialize_finite_difference_method
!*******************************************************************************

!*******************************************************************************
!>
!  Print the contents of a [[finite_diff_method]]. Used for debugging.

    subroutine print_finite_difference_method(me,iunit)

    implicit none

    class(finite_diff_method),intent(in) :: me
    integer,intent(in)                   :: iunit  !! file unit for printing
                                                   !! (assumed to be opened)

    write(iunit,'(A,1X,I5)')        'id            :', me%id
    write(iunit,'(A,1X,A)')         'name          :', me%name
    write(iunit,'(A,1X,I5)')        'class         :', me%class
    write(iunit,'(A,1X,*(I5,","))') 'dx_factors    :', int(me%dx_factors)
    write(iunit,'(A,1X,*(I5,","))') 'df_factors    :', int(me%df_factors)
    write(iunit,'(A,1X,I5)')        'df_den_factor :', int(me%df_den_factor)

    end subroutine print_finite_difference_method
!*******************************************************************************

!*******************************************************************************
!>
!  Wrapper for computing the function, using the cache.

    subroutine compute_function_with_cache(me,x,f,funcs_to_compute)

    implicit none

    class(numdiff_type),intent(inout) :: me
    real(wp),dimension(:),intent(in) :: x !! array of variables (size `n`)
    real(wp),dimension(:),intent(out) :: f !! array of functions (size `m`)
    integer,dimension(:),intent(in) :: funcs_to_compute !! the elements of the
                                                        !! function vector that need
                                                        !! to be computed (the other
                                                        !! are ignored)

    integer :: i !! index in the cache
    logical,dimension(size(funcs_to_compute)) :: ffound  !! functions found in the cache
    logical :: xfound  !! if `x` was found in the cache

    if (me%exception_raised) return ! check for exceptions

    call me%cache%get(x,funcs_to_compute,i,f,xfound,ffound)

    if (xfound .and. any(ffound)) then

        ! at least one of the functions was found
        if (all(ffound)) return ! all were found

        ! compute the ones that weren't found,
        ! and add them to the cache:
        call me%problem_func(x,f,pack(funcs_to_compute,mask=(.not. ffound)))
        call me%cache%put(i,x,f,pack(funcs_to_compute,mask=(.not. ffound)))

    else

        ! compute the function and add it to the cache:
        call me%problem_func(x,f,funcs_to_compute)
        call me%cache%put(i,x,f,funcs_to_compute)

    end if

    end subroutine compute_function_with_cache
!*******************************************************************************

!*******************************************************************************
!>
!  Return a string with the finite difference formula.
!
!### Example
!  * For 3-point backward: `dfdx = (f(x-2h)-4f(x-h)+3f(x)) / (2h)`

    subroutine get_formula(me,formula)

    class(finite_diff_method),intent(in) :: me
    character(len=:),allocatable,intent(out) :: formula

    integer :: i !! counter
    character(len=:),allocatable :: x !! temp variable for integer to string conversion
    character(len=:),allocatable :: f !! temp variable for integer to string conversion

    if (allocated(me%dx_factors) .and. allocated(me%df_factors)) then

        formula = 'dfdx = ('

        do i = 1, size(me%dx_factors)

            if (int(me%df_factors(i))==1) then
                if (i==1) then
                    formula = formula//'f('
                else
                    formula = formula//'+f('
                end if
            elseif (int(me%df_factors(i))==-1) then
                formula = formula//'-f('
            else
                if (i==1) then
                    f = integer_to_string(int(me%df_factors(i)))
                else
                    f = integer_to_string(int(me%df_factors(i)), with_sign = .true.)
                end if
                formula = formula//trim(adjustl(f))//'f('
            end if

            if (int(me%dx_factors(i))==0) then
                formula = formula//'x'
            elseif (int(me%dx_factors(i))==1) then
                formula = formula//'x+h'
            elseif (int(me%dx_factors(i))==-1) then
                formula = formula//'x-h'
            else
                x = integer_to_string(int(me%dx_factors(i)), with_sign = .true.)
                formula = formula//'x'//trim(adjustl(x))//'h'
            end if

            formula = formula//')'

        end do

        f = integer_to_string(int(me%df_den_factor))
        if (int(me%df_den_factor)==1) then
            formula = formula//') / h'
        else
            formula = formula//') / ('//trim(adjustl(f))//'h)'
        end if

    else
        formula = ''
    end if

    end subroutine get_formula
!*******************************************************************************

!*******************************************************************************
!>
!  Return a string with the finite difference formula.
!  Input is the method `id` code.
!
!###See also:
!  * [[get_formula]]

    subroutine get_finite_diff_formula(id,formula,name)

    implicit none

    integer,intent(in) :: id  !! the id code for the method
    character(len=:),allocatable,intent(out) :: formula !! the formula string
    character(len=:),allocatable,intent(out),optional :: name

    type(finite_diff_method) :: fd
    logical :: found

    call get_finite_difference_method(id,fd,found)
    if (found .and. fd%id/=0) then
        call get_formula(fd,formula)
        if (present(name)) name = fd%name
    else
        formula = ''
        if (present(name)) name = ''
    end if

    end subroutine get_finite_diff_formula
!*******************************************************************************

!*******************************************************************************
!>
!  Return a [[finite_diff_method]] given the `id` code.
!  (the `id` codes begin at 1, are sequential, and uniquely define the method).
!
!### Available methods
!
!   * \( (f(x+h)-f(x)) / h \)
!   * \( (f(x)-f(x-h)) / h \)
!   * \( (f(x+h)-f(x-h)) / (2h) \)
!   * \( (-3f(x)+4f(x+h)-f(x+2h)) / (2h) \)
!   * \( (f(x-2h)-4f(x-h)+3f(x)) / (2h) \)
!   * \( (-2f(x-h)-3f(x)+6f(x+h)-f(x+2h)) / (6h) \)
!   * \( (f(x-2h)-6f(x-h)+3f(x)+2f(x+h)) / (6h) \)
!   * \( (-11f(x)+18f(x+h)-9f(x+2h)+2f(x+3h)) / (6h) \)
!   * \( (-2f(x-3h)+9f(x-2h)-18f(x-h)+11f(x)) / (6h) \)
!   * \( (f(x-2h)-8f(x-h)+8f(x+h)-f(x+2h)) / (12h) \)
!   * \( (-3f(x-h)-10f(x)+18f(x+h)-6f(x+2h)+f(x+3h)) / (12h) \)
!   * \( (-f(x-3h)+6f(x-2h)-18f(x-h)+10f(x)+3f(x+h)) / (12h) \)
!   * \( (-25f(x)+48f(x+h)-36f(x+2h)+16f(x+3h)-3f(x+4h)) / (12h) \)
!   * \( (3f(x-4h)-16f(x-3h)+36f(x-2h)-48f(x-h)+25f(x)) / (12h) \)
!   * \( (3f(x-2h)-30f(x-h)-20f(x)+60f(x+h)-15f(x+2h)+2f(x+3h)) / (60h) \)
!   * \( (-2f(x-3h)+15f(x-2h)-60f(x-h)+20f(x)+30f(x+h)-3f(x+2h)) / (60h) \)
!   * \( (-12f(x-h)-65f(x)+120f(x+h)-60f(x+2h)+20f(x+3h)-3f(x+4h)) / (60h) \)
!   * \( (3f(x-4h)-20f(x-3h)+60f(x-2h)-120f(x-h)+65f(x)+12f(x+h)) / (60h) \)
!   * \( (-137f(x)+300f(x+h)-300f(x+2h)+200f(x+3h)-75f(x+4h)+12f(x+5h)) / (60h) \)
!   * \( (-12f(x-5h)+75f(x-4h)-200f(x-3h)+300f(x-2h)-300f(x-h)+137f(x)) / (60h) \)
!   * \( (-f(x-3h)+9f(x-2h)-45f(x-h)+45f(x+h)-9f(x+2h)+f(x+3h)) / (60h) \)
!   * \( (2f(x-2h)-24f(x-h)-35f(x)+80f(x+h)-30f(x+2h)+8f(x+3h)-f(x+4h)) / (60h) \)
!   * \( (f(x-4h)-8f(x-3h)+30f(x-2h)-80f(x-h)+35f(x)+24f(x+h)-2f(x+2h)) / (60h) \)
!   * \( (-10f(x-h)-77f(x)+150f(x+h)-100f(x+2h)+50f(x+3h)-15f(x+4h)+2f(x+5h)) / (60h) \)
!   * \( (-2f(x-5h)+15f(x-4h)-50f(x-3h)+100f(x-2h)-150f(x-h)+77f(x)+10f(x+h)) / (60h) \)
!   * \( (-147f(x)+360f(x+h)-450f(x+2h)+400f(x+3h)-225f(x+4h)+72f(x+5h)-10f(x+6h)) / (60h) \)
!   * \( (10f(x-6h)-72f(x-5h)+225f(x-4h)-400f(x-3h)+450f(x-2h)-360f(x-h)+147f(x)) / (60h) \)
!   * \( (-4f(x-3h)+42f(x-2h)-252f(x-h)-105f(x)+420f(x+h)-126f(x+2h)+28f(x+3h)-3f(x+4h)) / (420h) \)
!   * \( (3f(x-4h)-28f(x-3h)+126f(x-2h)-420f(x-h)+105f(x)+252f(x+h)-42f(x+2h)+4f(x+3h)) / (420h) \)
!   * \( (10f(x-2h)-140f(x-h)-329f(x)+700f(x+h)-350f(x+2h)+140f(x+3h)-35f(x+4h)+4f(x+5h)) / (420h) \)
!   * \( (-4f(x-5h)+35f(x-4h)-140f(x-3h)+350f(x-2h)-700f(x-h)+329f(x)+140f(x+h)-10f(x+2h)) / (420h) \)
!   * \( (-60f(x-h)-609f(x)+1260f(x+h)-1050f(x+2h)+700f(x+3h)-315f(x+4h)+84f(x+5h)-10f(x+6h)) / (420h) \)
!   * \( (10f(x-6h)-84f(x-5h)+315f(x-4h)-700f(x-3h)+1050f(x-2h)-1260f(x-h)+609f(x)+60f(x+h)) / (420h) \)
!   * \( (-1089f(x)+2940f(x+h)-4410f(x+2h)+4900f(x+3h)-3675f(x+4h)+1764f(x+5h)-490f(x+6h)+60f(x+7h)) / (420h) \)
!   * \( (-60f(x-7h)+490f(x-6h)-1764f(x-5h)+3675f(x-4h)-4900f(x-3h)+4410f(x-2h)-2940f(x-h)+1089f(x)) / (420h) \)
!   * \( (3f(x-4h)-32f(x-3h)+168f(x-2h)-672f(x-h)+672f(x+h)-168f(x+2h)+32f(x+3h)-3f(x+4h)) / (840h) \)
!   * \( (-5f(x-3h)+60f(x-2h)-420f(x-h)-378f(x)+1050f(x+h)-420f(x+2h)+140f(x+3h)-30f(x+4h)+3f(x+5h)) / (840h) \)
!   * \( (-3f(x-5h)+30f(x-4h)-140f(x-3h)+420f(x-2h)-1050f(x-h)+378f(x)+420f(x+h)-60f(x+2h)+5f(x+3h)) / (840h) \)
!   * \( (15f(x-2h)-240f(x-h)-798f(x)+1680f(x+h)-1050f(x+2h)+560f(x+3h)-210f(x+4h)+48f(x+5h)-5f(x+6h)) / (840h) \)
!   * \( (5f(x-6h)-48f(x-5h)+210f(x-4h)-560f(x-3h)+1050f(x-2h)-1680f(x-h)+798f(x)+240f(x+h)-15f(x+2h)) / (840h) \)
!   * \( (-105f(x-h)-1338f(x)+2940f(x+h)-2940f(x+2h)+2450f(x+3h)-1470f(x+4h)+588f(x+5h)-140f(x+6h)+15f(x+7h)) / (840h) \)
!   * \( (-15f(x-7h)+140f(x-6h)-588f(x-5h)+1470f(x-4h)-2450f(x-3h)+2940f(x-2h)-2940f(x-h)+1338f(x)+105f(x+h)) / (840h) \)
!   * \( (-2283f(x)+6720f(x+h)-11760f(x+2h)+15680f(x+3h)-14700f(x+4h)+9408f(x+5h)-3920f(x+6h)+960f(x+7h)-105f(x+8h)) / (840h) \)
!   * \( (105f(x-8h)-960f(x-7h)+3920f(x-6h)-9408f(x-5h)+14700f(x-4h)-15680f(x-3h)+11760f(x-2h)-6720f(x-h)+2283f(x)) / (840h) \)
!   * \( (-2f(x-5h)+25f(x-4h)-150f(x-3h)+600f(x-2h)-2100f(x-h)+2100f(x+h)-600f(x+2h)+150f(x+3h)-25f(x+4h)+2f(x+5h)) / (2520h) \)
!   * \( (5f(x-6h)-72f(x-5h)+495f(x-4h)-2200f(x-3h)+7425f(x-2h)-23760f(x-h)+23760f(x+h)-7425f(x+2h)+2200f(x+3h)-495f(x+4h)+
!        72f(x+5h)-5f(x+6h)) / (27720h) \)
!   * \( (-15f(x-7h)+245f(x-6h)-1911f(x-5h)+9555f(x-4h)-35035f(x-3h)+105105f(x-2h)-315315f(x-h)+315315f(x+h)-105105f(x+2h)+
!        35035f(x+3h)-9555f(x+4h)+1911f(x+5h)-245f(x+6h)+15f(x+7h)) / (360360h) \)
!   * \( (7f(x-8h)-128f(x-7h)+1120f(x-6h)-6272f(x-5h)+25480f(x-4h)-81536f(x-3h)+224224f(x-2h)-640640f(x-h)+640640f(x+h)-
!        224224f(x+2h)+81536f(x+3h)-25480f(x+4h)+6272f(x+5h)-1120f(x+6h)+128f(x+7h)-7f(x+8h)) / (720720h) \)
!
!  Where \(f(x)\) is the user-defined function of \(x\)
!  and \(h\) is a "small" perturbation.
!
!### References
!  * G. E. Mullges, F. Uhlig, "Numerical Algorithms with Fortran", Springer, 1996.
!  * G. W. Griffiths, "Numerical Analysis Using R", Cambridge University Press, 2016.
!
!@note This is the only routine that has to be changed if a new
!      finite difference method is added.
!
!@note The order within a class is assumed to be the order that we would prefer
!      to use them (e.g., central diffs are first, etc.) This is used in
!      the [[select_finite_diff_method]] routine.

    subroutine get_finite_difference_method(id,fd,found)

    implicit none

    integer,intent(in)                   :: id     !! the id code for the method
    type(finite_diff_method),intent(out) :: fd     !! this method (can be used
                                                   !! in [[compute_jacobian]])
    logical,intent(out)                  :: found  !! true if it was found

    found = .true.

    select case (id)

    ! 2 point methods
    case(1)
        ! (f(x+h)-f(x)) / h
        fd = finite_diff_method(id,'2-point forward 1',    2,[1,0],[1,-1],1)
    case(2)
        ! (f(x)-f(x-h)) / h
        fd = finite_diff_method(id,'2-point backward 1',   2,[0,-1],[1,-1],1)

    ! 3 point methods
    case(3)
        ! (f(x+h)-f(x-h)) / (2h)
        fd = finite_diff_method(id,'3-point central',    3,[1,-1],[1,-1],2)
    case(4)
        ! (-3f(x)+4f(x+h)-f(x+2h)) / (2h)
        fd = finite_diff_method(id,'3-point forward 2',    3,[0,1,2],[-3,4,-1],2)
    case(5)
        ! (f(x-2h)-4f(x-h)+3f(x)) / (2h)
        fd = finite_diff_method(id,'3-point backward 2',   3,[-2,-1,0],[1,-4,3],2)

    ! 4 point methods
    case(6)
        ! (-2f(x-h)-3f(x)+6f(x+h)-f(x+2h)) / (6h)
        fd = finite_diff_method(id,'4-point forward 2',  4,[-1,0,1,2],[-2,-3,6,-1],6)
    case(7)
        ! (f(x-2h)-6f(x-h)+3f(x)+2f(x+h)) / (6h)
        fd = finite_diff_method(id,'4-point backward 2', 4,[-2,-1,0,1],[1,-6,3,2],6)
    case(8)
        ! (-11f(x)+18f(x+h)-9f(x+2h)+2f(x+3h)) / (6h)
        fd = finite_diff_method(id,'4-point forward 3',  4,[0,1,2,3],[-11,18,-9,2],6)
    case(9)
        ! (-2f(x-3h)+9f(x-2h)-18f(x-h)+11f(x)) / (6h)
        fd = finite_diff_method(id,'4-point backward 3', 4,[-3,-2,-1,0],[-2,9,-18,11],6)

    ! 5 point methods
    case(10)
        ! (f(x-2h)-8f(x-h)+8f(x+h)-f(x+2h)) / (12h)
        fd = finite_diff_method(id,'5-point central',    5,[-2,-1,1,2],[1,-8,8,-1],12)
    case(11)
        ! (-3f(x-h)-10f(x)+18f(x+h)-6f(x+2h)+f(x+3h)) / (12h)
        fd = finite_diff_method(id,'5-point forward 3',  5,[-1,0,1,2,3],[-3,-10,18,-6,1],12)
    case(12)
        ! (-f(x-3h)+6f(x-2h)-18f(x-h)+10f(x)+3f(x+h)) / (12h)
        fd = finite_diff_method(id,'5-point backward 3', 5,[-3,-2,-1,0,1],[-1,6,-18,10,3],12)
    case(13)
        ! (-25f(x)+48f(x+h)-36f(x+2h)+16f(x+3h)-3f(x+4h)) / (12h)
        fd = finite_diff_method(id,'5-point forward 4',  5,[0,1,2,3,4],[-25,48,-36,16,-3],12)
    case(14)
        ! (3f(x-4h)-16f(x-3h)+36f(x-2h)-48f(x-h)+25f(x)) / (12h)
        fd = finite_diff_method(id,'5-point backward 4', 5,[-4,-3,-2,-1,0],[3,-16,36,-48,25],12)

    ! 6 point methods
    case(15)
        ! (3f(x-2h)-30f(x-h)-20f(x)+60f(x+h)-15f(x+2h)+2f(x+3h)) / (60h)
        fd = finite_diff_method(id,'6-point forward 3',  6,[-2,-1,0,1,2,3],[3,-30,-20,60,-15,2],60)
    case(16)
        ! (-2f(x-3h)+15f(x-2h)-60f(x-h)+20f(x)+30f(x+h)-3f(x+2h)) / (60h)
        fd = finite_diff_method(id,'6-point backward 3', 6,[-3,-2,-1,0,1,2],[-2,15,-60,20,30,-3],60)
    case(17)
        ! (-12f(x-h)-65f(x)+120f(x+h)-60f(x+2h)+20f(x+3h)-3f(x+4h)) / (60h)
        fd = finite_diff_method(id,'6-point forward 4',  6,[-1,0,1,2,3,4],[-12,-65,120,-60,20,-3],60)
    case(18)
        ! (3f(x-4h)-20f(x-3h)+60f(x-2h)-120f(x-h)+65f(x)+12f(x+h)) / (60h)
        fd = finite_diff_method(id,'6-point backward 4', 6,[-4,-3,-2,-1,0,1],[3,-20,60,-120,65,12],60)
    case(19)
        ! (-137f(x)+300f(x+h)-300f(x+2h)+200f(x+3h)-75f(x+4h)+12f(x+5h)) / (60h)
        fd = finite_diff_method(id,'6-point forward 5',  6,[0,1,2,3,4,5],[-137,300,-300,200,-75,12],60)
    case(20)
        ! (-12f(x-5h)+75f(x-4h)-200f(x-3h)+300f(x-2h)-300f(x-h)+137f(x)) / (60h)
        fd = finite_diff_method(id,'6-point backward 5', 6,[-5,-4,-3,-2,-1,0],[-12,75,-200,300,-300,137],60)

    ! 7 point methods
    case(21)
        ! (-f(x-3h)+9f(x-2h)-45f(x-h)+45f(x+h)-9f(x+2h)+f(x+3h)) / (60h)
        fd = finite_diff_method(id,'7-point central', 7,[-3,-2,-1,1,2,3],[-1,9,-45,45,-9,1],60)
    case(22)
        ! (2f(x-2h)-24f(x-h)-35f(x)+80f(x+h)-30f(x+2h)+8f(x+3h)-f(x+4h)) / (60h)
        fd = finite_diff_method(id,'7-point forward 4', 7,[-2,-1,0,1,2,3,4],[2,-24,-35,80,-30,8,-1],60)
    case(23)
        ! (f(x-4h)-8f(x-3h)+30f(x-2h)-80f(x-h)+35f(x)+24f(x+h)-2f(x+2h)) / (60h)
        fd = finite_diff_method(id,'7-point backward 4', 7,[-4,-3,-2,-1,0,1,2],[1,-8,30,-80,35,24,-2],60)
    case(24)
        ! (-10f(x-h)-77f(x)+150f(x+h)-100f(x+2h)+50f(x+3h)-15f(x+4h)+2f(x+5h)) / (60h)
        fd = finite_diff_method(id,'7-point forward 5', 7,[-1,0,1,2,3,4,5],[-10,-77,150,-100,50,-15,2],60)
    case(25)
        ! (-2f(x-5h)+15f(x-4h)-50f(x-3h)+100f(x-2h)-150f(x-h)+77f(x)+10f(x+h)) / (60h)
        fd = finite_diff_method(id,'7-point backward 5', 7,[-5,-4,-3,-2,-1,0,1],[-2,15,-50,100,-150,77,10],60)
    case(26)
        ! (-147f(x)+360f(x+h)-450f(x+2h)+400f(x+3h)-225f(x+4h)+72f(x+5h)-10f(x+6h)) / (60h)
        fd = finite_diff_method(id,'7-point forward 6', 7,[0,1,2,3,4,5,6],[-147,360,-450,400,-225,72,-10],60)
    case(27)
        ! (10f(x-6h)-72f(x-5h)+225f(x-4h)-400f(x-3h)+450f(x-2h)-360f(x-h)+147f(x)) / (60h)
        fd = finite_diff_method(id,'7-point backward 6', 7,[-6,-5,-4,-3,-2,-1,0],[10,-72,225,-400,450,-360,147],60)

    ! 8 point methods
    case(28)
        ! (-4f(x-3h)+42f(x-2h)-252f(x-h)-105f(x)+420f(x+h)-126f(x+2h)+28f(x+3h)-3f(x+4h)) / (420h)
        fd = finite_diff_method(id,'8-point forward 4',  8,[-3,-2,-1,0,1,2,3,4],[-4,42,-252,-105,420,-126,28,-3],420)
    case(29)
        ! (3f(x-4h)-28f(x-3h)+126f(x-2h)-420f(x-h)+105f(x)+252f(x+h)-42f(x+2h)+4f(x+3h)) / (420h)
        fd = finite_diff_method(id,'8-point backward 4', 8,[-4,-3,-2,-1,0,1,2,3],[3,-28,126,-420,105,252,-42,4],420)
    case(30)
        ! (10f(x-2h)-140f(x-h)-329f(x)+700f(x+h)-350f(x+2h)+140f(x+3h)-35f(x+4h)+4f(x+5h)) / (420h)
        fd = finite_diff_method(id,'8-point forward 5',  8,[-2,-1,0,1,2,3,4,5],[10,-140,-329,700,-350,140,-35,4],420)
    case(31)
        ! (-4f(x-5h)+35f(x-4h)-140f(x-3h)+350f(x-2h)-700f(x-h)+329f(x)+140f(x+h)-10f(x+2h)) / (420h)
        fd = finite_diff_method(id,'8-point backward 5', 8,[-5,-4,-3,-2,-1,0,1,2],[-4,35,-140,350,-700,329,140,-10],420)
    case(32)
        ! (-60f(x-h)-609f(x)+1260f(x+h)-1050f(x+2h)+700f(x+3h)-315f(x+4h)+84f(x+5h)-10f(x+6h)) / (420h)
        fd = finite_diff_method(id,'8-point forward 6',  8,[-1,0,1,2,3,4,5,6],[-60,-609,1260,-1050,700,-315,84,-10],420)
    case(33)
        ! (10f(x-6h)-84f(x-5h)+315f(x-4h)-700f(x-3h)+1050f(x-2h)-1260f(x-h)+609f(x)+60f(x+h)) / (420h)
        fd = finite_diff_method(id,'8-point backward 6', 8,[-6,-5,-4,-3,-2,-1,0,1],[10,-84,315,-700,1050,-1260,609,60],420)
    case(34)
        ! (-1089f(x)+2940f(x+h)-4410f(x+2h)+4900f(x+3h)-3675f(x+4h)+1764f(x+5h)-490f(x+6h)+60f(x+7h)) / (420h)
        fd = finite_diff_method(id,'8-point forward 7',  8,[0,1,2,3,4,5,6,7],[-1089,2940,-4410,4900,-3675,1764,-490,60],420)
    case(35)
        ! (-60f(x-7h)+490f(x-6h)-1764f(x-5h)+3675f(x-4h)-4900f(x-3h)+4410f(x-2h)-2940f(x-h)+1089f(x)) / (420h)
        fd = finite_diff_method(id,'8-point backward 7', 8,[-7,-6,-5,-4,-3,-2,-1,0],[-60,490,-1764,3675,-4900,4410,-2940,1089],420)

    ! 9 point methods
    case(36)
        ! (3f(x-4h)-32f(x-3h)+168f(x-2h)-672f(x-h)+672f(x+h)-168f(x+2h)+32f(x+3h)-3f(x+4h)) / (840h)
        fd = finite_diff_method(id,'9-point central',    &
                9,[-4,-3,-2,-1,1,2,3,4],[3,-32,168,-672,672,-168,32,-3],840)
    case(37)
        ! (-5f(x-3h)+60f(x-2h)-420f(x-h)-378f(x)+1050f(x+h)-420f(x+2h)+140f(x+3h)-30f(x+4h)+3f(x+5h)) / (840h)
        fd = finite_diff_method(id,'9-point forward 5',  &
                9,[-3,-2,-1,0,1,2,3,4,5],[-5,60,-420,-378,1050,-420,140,-30,3],840)
    case(38)
        ! (-3f(x-5h)+30f(x-4h)-140f(x-3h)+420f(x-2h)-1050f(x-h)+378f(x)+420f(x+h)-60f(x+2h)+5f(x+3h)) / (840h)
        fd = finite_diff_method(id,'9-point backward 5', &
                9,[-5,-4,-3,-2,-1,0,1,2,3],[-3,30,-140,420,-1050,378,420,-60,5],840)
    case(39)
        ! (15f(x-2h)-240f(x-h)-798f(x)+1680f(x+h)-1050f(x+2h)+560f(x+3h)-210f(x+4h)+48f(x+5h)-5f(x+6h)) / (840h)
        fd = finite_diff_method(id,'9-point forward 6',  &
                9,[-2,-1,0,1,2,3,4,5,6],[15,-240,-798,1680,-1050,560,-210,48,-5],840)
    case(40)
        ! (5f(x-6h)-48f(x-5h)+210f(x-4h)-560f(x-3h)+1050f(x-2h)-1680f(x-h)+798f(x)+240f(x+h)-15f(x+2h)) / (840h)
        fd = finite_diff_method(id,'9-point backward 6', &
                9,[-6,-5,-4,-3,-2,-1,0,1,2],[5,-48,210,-560,1050,-1680,798,240,-15],840)
    case(41)
        ! (-105f(x-h)-1338f(x)+2940f(x+h)-2940f(x+2h)+2450f(x+3h)-1470f(x+4h)+588f(x+5h)-140f(x+6h)+15f(x+7h)) / (840h)
        fd = finite_diff_method(id,'9-point forward 7',  &
                9,[-1,0,1,2,3,4,5,6,7],[-105,-1338,2940,-2940,2450,-1470,588,-140,15],840)
    case(42)
        ! (-15f(x-7h)+140f(x-6h)-588f(x-5h)+1470f(x-4h)-2450f(x-3h)+2940f(x-2h)-2940f(x-h)+1338f(x)+105f(x+h)) / (840h)
        fd = finite_diff_method(id,'9-point backward 7', &
                9,[-7,-6,-5,-4,-3,-2,-1,0,1],[-15,140,-588,1470,-2450,2940,-2940,1338,105],840)
    case(43)
        ! (-2283f(x)+6720f(x+h)-11760f(x+2h)+15680f(x+3h)-14700f(x+4h)+9408f(x+5h)-3920f(x+6h)+960f(x+7h)-105f(x+8h)) / (840h)
        fd = finite_diff_method(id,'9-point forward 8',  &
                9,[0,1,2,3,4,5,6,7,8],[-2283,6720,-11760,15680,-14700,9408,-3920,960,-105],840)
    case(44)
        ! (105f(x-8h)-960f(x-7h)+3920f(x-6h)-9408f(x-5h)+14700f(x-4h)-15680f(x-3h)+11760f(x-2h)-6720f(x-h)+2283f(x)) / (840h)
        fd = finite_diff_method(id,'9-point backward 8', &
                9,[-8,-7,-6,-5,-4,-3,-2,-1,0],[105,-960,3920,-9408,14700,-15680,11760,-6720,2283],840)

    ! 11 point methods
    case(500)
        ! (-2f(x-5h)+25f(x-4h)-150f(x-3h)+600f(x-2h)-2100f(x-h)+2100f(x+h)-600f(x+2h)+150f(x+3h)-25f(x+4h)+2f(x+5h)) / (2520h)
        fd = finite_diff_method(id,'11-point central', &
                11,[-5,-4,-3,-2,-1,1,2,3,4,5],[-2,25,-150,600,-2100,2100,-600,150,-25,2],2520)

    ! 13 point methods
    case(600)
        ! (5f(x-6h)-72f(x-5h)+495f(x-4h)-2200f(x-3h)+7425f(x-2h)-23760f(x-h)+23760f(x+h)-7425f(x+2h)+2200f(x+3h)-495f(x+4h)+72f(x+5h)-5f(x+6h)) / (27720h)
        fd = finite_diff_method(id,'13-point central', &
                13,[-6,-5,-4,-3,-2,-1,1,2,3,4,5,6],&
                   [5,-72,495,-2200,7425,-23760,23760,-7425,2200,-495,72,-5],27720)

    ! 15 point methods
    case(700)
        ! (-15f(x-7h)+245f(x-6h)-1911f(x-5h)+9555f(x-4h)-35035f(x-3h)+105105f(x-2h)-315315f(x-h)+315315f(x+h)-105105f(x+2h)+35035f(x+3h)-9555f(x+4h)+1911f(x+5h)-245f(x+6h)+15f(x+7h)) / (360360h)
        fd = finite_diff_method(id,'15-point central', &
                15,[-7,-6,-5,-4,-3,-2,-1,1,2,3,4,5,6,7],&
                   [-15,245,-1911,9555,-35035,105105,-315315,315315,-105105,35035,-9555,1911,-245,15],360360)

    ! 17 point methods
    case(800)
        ! (7f(x-8h)-128f(x-7h)+1120f(x-6h)-6272f(x-5h)+25480f(x-4h)-81536f(x-3h)+224224f(x-2h)-640640f(x-h)+640640f(x+h)-224224f(x+2h)+81536f(x+3h)-25480f(x+4h)+6272f(x+5h)-1120f(x+6h)+128f(x+7h)-7f(x+8h)) / (720720h)
        fd = finite_diff_method(id,'17-point central', &
                17,[-8,-7,-6,-5,-4,-3,-2,-1,1,2,3,4,5,6,7,8],&
                   [7,-128,1120,-6272,25480,-81536,224224,-640640,640640,-224224,81536,-25480,6272,-1120,128,-7],720720)

    case default
        found = .false.
    end select

    end subroutine get_finite_difference_method
!*******************************************************************************

!*******************************************************************************
!>
!  Returns all the methods with the given `class`
!  (i.e., number of points in the formula).

    function get_all_methods_in_class(class) result(list_of_methods)

    implicit none

    integer,intent(in) :: class  !! the `class` is the number of points in the
                                 !! finite different computation
    type(meth_array) :: list_of_methods  !! all the methods with the specified `class`

    type(finite_diff_method) :: fd  !! temp variable for getting a
                                    !! method from [[get_finite_difference_method]]
    integer :: id     !! method id counter
    logical :: found  !! status flag
    integer :: n      !! temp size variable
    type(finite_diff_method),dimension(:),allocatable :: tmp  !! for resizing `meth` array

    ! currently, the only way to do this is to call the
    ! get_finite_difference_method routine and see if there
    ! is one available.
    id = 0
    do
        id = id + 1
        call get_finite_difference_method(id,fd,found)
        if (found) then
            if (fd%class==class) then
                if (allocated(list_of_methods%meth)) then
                    ! add to the list
                    n = size(list_of_methods%meth)
                    allocate(tmp(n+1))
                    tmp(1:n) = list_of_methods%meth
                    tmp(n+1) = fd
                    call move_alloc(tmp,list_of_methods%meth)
                    ! ... this doesn't appear to work on Intel compiler ...
                    !list_of_methods%meth = [list_of_methods%meth,fd]  ! add to the list
                else
                    ! first element:
                    allocate(list_of_methods%meth(1))
                    list_of_methods%meth = fd
                end if
            elseif (fd%class>class) then ! we assume they are in increasing order
                exit
            end if
        else
            exit ! done
        end if
    end do

    end function get_all_methods_in_class
!*******************************************************************************

!*******************************************************************************
!>
!  Select a finite diff method of a given `class` so that the perturbations
!  of `x` will not violate the variable bounds.

    subroutine select_finite_diff_method(me,x,xlow,xhigh,dx,list_of_methods,fd,status_ok)

    implicit none

    class(numdiff_type),intent(inout)    :: me
    real(wp),intent(in)                  :: x               !! the variable value
    real(wp),intent(in)                  :: xlow            !! the variable lower bound
    real(wp),intent(in)                  :: xhigh           !! the variable upper bound
    real(wp),intent(in)                  :: dx              !! the perturbation value (>0)
    type(meth_array),intent(in)          :: list_of_methods !! list of available methods to choose from
    type(finite_diff_method),intent(out) :: fd              !! this method can be used
    logical,intent(out)                  :: status_ok       !! true if it really doesn't violate the bounds
                                                            !! (say, the bounds are very close to each other)
                                                            !! if `status_ok=False`, then the first method in
                                                            !! the given class is returned in `fd`.

    integer  :: i   !! counter
    integer  :: j   !! counter
    real(wp) :: xs  !! if `x` is outside the bounds, this is the value
                    !! on the nearest bound. otherwise equal to `x`.
    real(wp) :: xp  !! perturbed `x` value

    ! initialize:
    status_ok = .false.

    if (me%exception_raised) return ! check for exceptions

    ! make sure it is within the bounds
    xs = min(xhigh,max(xlow,x))

    ! try all the methods in the class:
    do i = 1, size(list_of_methods%meth)
        status_ok = .true. ! will be set to false if any
                           ! perturbation violates the bounds
        ! check each of the perturbations:
        do j = 1, size(list_of_methods%meth(i)%dx_factors)
            xp = xs + list_of_methods%meth(i)%dx_factors(j)*dx
            if (xp < xlow .or. xp > xhigh) then
                status_ok = .false.
                exit
            end if
        end do
        if (status_ok) then   ! this one is OK to use
            fd = list_of_methods%meth(i)
            exit
        end if
    end do

    if (.not. status_ok) then
        ! no method was found that doesn't violate the bounds,
        ! so just return the first one in the list.
        fd = list_of_methods%meth(1)
    end if

    end subroutine select_finite_diff_method
!*******************************************************************************

!*******************************************************************************
!>
!  Select a finite diff method of a given `class` so that the perturbations
!  of `x` will not violate the variable bounds for any variable in the group.
!
!  The `x` vector are only the variables in a group (not the full variable vector)

    subroutine select_finite_diff_method_for_partition_group(me,x,xlow,xhigh,dx,&
                                                             list_of_methods,fd,status_ok)

    implicit none

    class(numdiff_type),intent(inout)     :: me
    real(wp),dimension(:),intent(in)      :: x               !! the variable values
    real(wp),dimension(:),intent(in)      :: xlow            !! the variable lower bounds
    real(wp),dimension(:),intent(in)      :: xhigh           !! the variable upper bounds
    real(wp),dimension(:),intent(in)      :: dx              !! the perturbation values (>0)
    type(meth_array),intent(in)           :: list_of_methods !! list of available methods to choose from
    type(finite_diff_method),intent(out)  :: fd              !! this method can be used
    logical,intent(out)                   :: status_ok       !! true if it really doesn't violate the bounds
                                                             !! (say, the bounds are very close to each other)
                                                             !! if `status_ok=False`, then the first method in
                                                             !! the given class is returned in `fd`.

    integer  :: i   !! counter
    integer  :: j   !! counter
    real(wp),dimension(size(x)) :: xp  !! perturbed `x` values
    real(wp),dimension(size(x)) :: xs  !! if `x` is outside the bounds, this is the value
                                       !! on the nearest bound. otherwise equal to `x`.

    ! initialize:
    status_ok = .false.

    if (me%exception_raised) return ! check for exceptions

    ! make sure they are within the bounds
    xs = min(xhigh,max(xlow,x))

    ! try all the methods in the class:
    do i = 1, size(list_of_methods%meth)
        status_ok = .true. ! will be set to false if any
                           ! perturbation violates the bounds
        ! check each of the perturbations:
        do j = 1, size(list_of_methods%meth(i)%dx_factors)
            xp = xs + list_of_methods%meth(i)%dx_factors(j)*dx
            if (any(xp < xlow) .or. any(xp > xhigh)) then
                status_ok = .false.
                exit
            end if
        end do
        if (status_ok) then   ! this one is OK to use
            fd = list_of_methods%meth(i)
            exit
        end if
    end do

    if (.not. status_ok) then
        ! no method was found that doesn't violate the bounds,
        ! so just return the first one in the list.
        fd = list_of_methods%meth(1)
    end if

    end subroutine select_finite_diff_method_for_partition_group
!*******************************************************************************

!*******************************************************************************
!>
!  Alternate version of [[initialize_numdiff]] routine when
!  using [[diff]] to compute the Jacobian.

    subroutine initialize_numdiff_for_diff(me,n,m,xlow,xhigh,&
                                    problem_func,sparsity_mode,info,&
                                    chunk_size,eps,acc,cache_size,&
                                    xlow_for_sparsity,xhigh_for_sparsity,&
                                    dpert_for_sparsity,sparsity_perturb_mode,&
                                    linear_sparsity_tol,function_precision_tol,&
                                    num_sparsity_points,print_messages)

    implicit none

    class(numdiff_type),intent(inout):: me
    integer,intent(in)               :: n             !! number of `x` variables
    integer,intent(in)               :: m             !! number of `f` functions
    real(wp),dimension(n),intent(in) :: xlow          !! lower bounds on `x`
    real(wp),dimension(n),intent(in) :: xhigh         !! upper bounds on `x`
    procedure(func)                  :: problem_func  !! the user function that defines the problem
                                                      !! (returns `m` functions)
    integer,intent(in)               :: sparsity_mode !! the sparsity computation method:
                                                      !!
                                                      !! **1** - assume dense,
                                                      !! **2** - three-point simple method,
                                                      !! **3** - will be specified by the user in
                                                      !! a subsequent call to [[set_sparsity_pattern]].
                                                      !! **4** - computes a two-point jacobian
                                                      !! at `num_sparsity_points` points.
    procedure(info_f),optional       :: info          !! a function the user can define
                                                      !! which is called when each column
                                                      !! of the jacobian is computed.
                                                      !! It can be used to perform any
                                                      !! setup operations.
    integer,intent(in),optional      :: chunk_size    !! chunk size for allocating the arrays
                                                      !! (must be >0) [default is 100]
    real(wp),intent(in),optional     :: eps           !! tolerance parameter for [[diff]]
                                                      !! if not present, default is `1.0e-9_wp`
    real(wp),intent(in),optional     :: acc           !! tolerance parameter for [[diff]]
                                                      !! if not present, default is `0.0_wp`
    integer,intent(in),optional      :: cache_size    !! if present, this is the cache size
                                                      !! for the function cache
                                                      !! (default is not to enable cache)
    real(wp),dimension(n),intent(in),optional :: xlow_for_sparsity   !! lower bounds on `x` used for
                                                                     !! sparsity computation (when
                                                                     !! `sparsity_mode` is 2). If not
                                                                     !! present, then `xlow` is used.
    real(wp),dimension(n),intent(in),optional :: xhigh_for_sparsity  !! upper bounds on `x` used for
                                                                     !! sparsity computation (when
                                                                     !! `sparsity_mode` is 2). If not
                                                                     !! present, then `xhigh` is used.
    real(wp),dimension(n),intent(in),optional :: dpert_for_sparsity  !! required if `sparsity_mode=4`
    integer,intent(in),optional :: sparsity_perturb_mode   !! perturbation mode (required if `sparsity_mode=4`):
                                                           !!
                                                           !! **1** - perturbation is `dx=dpert`,
                                                           !! **2** - perturbation is `dx=dpert*x`,
                                                           !! **3** - perturbation is `dx=dpert*(1+x)`
    integer,intent(in),optional :: num_sparsity_points  !! for `sparsity_mode=4`, the number of jacobian
                                                        !! evaluations used to estimate the sparsity pattern.
    real(wp),intent(in),optional :: linear_sparsity_tol !! the equality tolerance for derivatives to
                                                        !! indicate a constant jacobian element (linear sparsity)
    real(wp),intent(in),optional :: function_precision_tol  !! the function precision. two functions values
                                                            !! that are the within this tolerance are
                                                            !! considered the same value. This is used
                                                            !! when estimating the sparsity pattern when
                                                            !! `sparsity_mode=2` in [[compute_sparsity_random]]
    logical,intent(in),optional :: print_messages !! if true, print error messages to `error_unit`.
                                                  !! default is True.

    logical :: cache  !! if the cache is to be used

    me%use_diff = .true.

    ! cache:
    cache = present(cache_size)
    if (cache) cache = cache_size>0
    if (cache) then
        me%compute_function  => compute_function_with_cache
        call me%cache%initialize(cache_size,n,m)
    else
        me%compute_function  => problem_func
        call me%cache%destroy()
    end if

    ! functions:
    me%problem_func      => problem_func
    me%jacobian_function => compute_jacobian_with_diff

    ! size of the problem:
    me%n = n
    me%m = m

    ! input variable bounds:
    call me%set_numdiff_bounds(xlow,xhigh)

    ! set the sparsity function
    call me%set_sparsity_mode(sparsity_mode,xlow_for_sparsity,xhigh_for_sparsity)

    if (sparsity_mode==4) then
        ! these must be present, since we don't have a dpert and perturb mode for the gradients:
        if (present(dpert_for_sparsity) .and. present(sparsity_perturb_mode)) then
            me%dpert_for_sparsity = abs(dpert_for_sparsity)
            ! perturbation options:
            select case (sparsity_perturb_mode)
            case(1:3)
                me%sparsity_perturb_mode = sparsity_perturb_mode
            case default
                call me%raise_exception(1,'initialize_numdiff_for_diff',&
                                          'sparsity_perturb_mode must be 1, 2, or 3.')
                return
            end select
        else
            call me%raise_exception(2,'initialize_numdiff_for_diff',&
                                      'missing required inputs for sparsity mode 4')
            return
        end if
    end if
    ! if these aren't present, they will just keep the defaults:
    if (present(linear_sparsity_tol))    me%linear_sparsity_tol    = linear_sparsity_tol
    if (present(function_precision_tol)) me%function_precision_tol = function_precision_tol
    if (present(num_sparsity_points))    me%num_sparsity_points    = num_sparsity_points

    ! optional:
    if (present(chunk_size))     me%chunk_size = abs(chunk_size)
    if (present(eps))            me%eps = eps
    if (present(acc))            me%acc = acc
    if (present(info))           me%info_function => info
    if (present(print_messages)) me%print_messages = print_messages

    end subroutine initialize_numdiff_for_diff
!*******************************************************************************

!*******************************************************************************
!>
!  Change the variable bounds in a [[numdiff_type]].
!
!### See also
!  * [[set_numdiff_sparsity_bounds]]
!
!@note The bounds must be set when the class is initialized,
!      but this routine can be used to change them later if required.

    subroutine set_numdiff_bounds(me,xlow,xhigh)

    implicit none

    class(numdiff_type),intent(inout) :: me
    real(wp),dimension(:),intent(in)  :: xlow    !! lower bounds on `x`
    real(wp),dimension(:),intent(in)  :: xhigh   !! upper bounds on `x`

    integer :: i  !! counter for error print
    character(len=:),allocatable :: error_info !! error message info
    character(len=:),allocatable :: istr   !! for integer to string
    character(len=30) :: xlow_str, xhigh_str !! for real to string

    if (me%exception_raised) return ! check for exceptions

    if (allocated(me%xlow))  deallocate(me%xlow)
    if (allocated(me%xhigh)) deallocate(me%xhigh)

    if (size(xlow)/=me%n .or. size(xhigh)/=me%n) then
        call me%raise_exception(3,'set_numdiff_bounds',&
                                  'invalid size of xlow or xhigh')
        return
    else if (any(xlow>=xhigh)) then
        error_info = 'all xlow must be < xhigh'
        do i = 1, size(xlow)
            if (xlow(i)>=xhigh(i)) then
                istr = integer_to_string(i)
                write(xlow_str,'(F30.16)') xlow(i)
                write(xhigh_str,'(F30.16)') xhigh(i)
                error_info = error_info//new_line('')//'  Error for optimization variable '//trim(adjustl(istr))//&
                                ': xlow='//trim(adjustl(xlow_str))//&
                                ' >= xhigh='//trim(adjustl(xhigh_str))
            end if
        end do
        call me%raise_exception(4,'set_numdiff_bounds',error_info)
        return
    else
        allocate(me%xlow(me%n))
        allocate(me%xhigh(me%n))
        me%xlow  = xlow
        me%xhigh = xhigh
    end if

    end subroutine set_numdiff_bounds
!*******************************************************************************

!*******************************************************************************
!>
!  Set sparsity mode.

    subroutine set_sparsity_mode(me,sparsity_mode,xlow_for_sparsity,xhigh_for_sparsity)

    implicit none

    class(numdiff_type),intent(inout) :: me
    integer,intent(in) :: sparsity_mode !! the sparsity computation method:
                                        !! **1** - assume dense,
                                        !! **2** - three-point simple method,
                                        !! **3** - will be specified by the user in
                                        !! a subsequent call to [[set_sparsity_pattern]].
                                        !! **4** - computes a two-point jacobian
                                        !! at `num_sparsity_points` points.
    real(wp),dimension(:),intent(in),optional :: xlow_for_sparsity   !! lower bounds on `x` used for
                                                                     !! sparsity computation (when
                                                                     !! `sparsity_mode` is 2). If not
                                                                     !! present, then `xlow` is used.
    real(wp),dimension(:),intent(in),optional :: xhigh_for_sparsity  !! upper bounds on `x` used for
                                                                     !! sparsity computation (when
                                                                     !! `sparsity_mode` is 2). If not
                                                                     !! present, then `xhigh` is used.

    if (me%exception_raised) return ! check for exceptions

    ! set the sparsity function
    select case (sparsity_mode)
    case(1)  ! dense
        me%compute_sparsity => compute_sparsity_dense
    case(3)  ! user defined
        me%compute_sparsity => null()
    case(2)  ! three-point simple method
        me%compute_sparsity => compute_sparsity_random
        ! in this case, we have the option of specifying
        ! separate bounds for computing the sparsity:
        call me%set_numdiff_sparsity_bounds(xlow_for_sparsity,xhigh_for_sparsity)
    case(4)  ! compute 2-point jacobian in specified number of points
        me%compute_sparsity => compute_sparsity_random_2
        ! in this case, we have the option of specifying
        ! separate bounds for computing the sparsity:
        call me%set_numdiff_sparsity_bounds(xlow_for_sparsity,xhigh_for_sparsity)
    case default
        call me%raise_exception(5,'set_sparsity_mode',&
                                  'sparsity_mode must be 1, 2, 3, or 4.')
        return
    end select

    end subroutine set_sparsity_mode
!*******************************************************************************

!*******************************************************************************
!>
!  Sets the variable bounds for sparsity in a [[numdiff_type]].
!  These are only used for `sparsity_mode=2`.
!
!### See also
!  * [[set_numdiff_bounds]]
!
!@note This routine assumes that `xlow` and `xhigh` have already
!      been set in the class.

    subroutine set_numdiff_sparsity_bounds(me,xlow,xhigh)

    implicit none

    class(numdiff_type),intent(inout) :: me
    real(wp),dimension(:),intent(in),optional  :: xlow  !! lower bounds on `x` to be used for
                                                        !! sparsity computation. If not present,
                                                        !! then then `xlow` values in the class are used.
    real(wp),dimension(:),intent(in),optional  :: xhigh !! upper bounds on `x` to be used for
                                                        !! sparsity computation. If not present,
                                                        !! then then `xhigh` values in the class are used.

    if (me%exception_raised) return ! check for exceptions

    if (allocated(me%xlow_for_sparsity)) deallocate(me%xlow_for_sparsity)
    if (allocated(me%xhigh_for_sparsity)) deallocate(me%xhigh_for_sparsity)

    if (.not. present(xlow)) then
        allocate(me%xlow_for_sparsity(size(me%xlow)))
        me%xlow_for_sparsity = me%xlow
    else
        allocate(me%xlow_for_sparsity(size(xlow)))
        me%xlow_for_sparsity = xlow
    end if

    if (.not. present(xhigh)) then
        allocate(me%xhigh_for_sparsity(size(me%xhigh)))
        me%xhigh_for_sparsity = me%xhigh
    else
        allocate(me%xhigh_for_sparsity(size(xhigh)))
        me%xhigh_for_sparsity = xhigh
    end if

    ! error checks:
    if (size(me%xlow_for_sparsity)/=me%n .or. size(me%xhigh_for_sparsity)/=me%n) then
        call me%raise_exception(6,'set_numdiff_sparsity_bounds',&
                                  'invalid size of xlow or xhigh')
    else if (any(me%xlow_for_sparsity>=me%xhigh_for_sparsity)) then
        call me%raise_exception(7,'set_numdiff_sparsity_bounds',&
                                  'all xlow must be < xhigh')
    end if

    end subroutine set_numdiff_sparsity_bounds
!*******************************************************************************

!*******************************************************************************
!>
!  Initialize a [[numdiff_type]] class. This must be called first.
!
!@note Only one of the following inputs can be used: `jacobian_method`,
!      `jacobian_methods`, `class`, or `classes`.

    subroutine initialize_numdiff(me,n,m,xlow,xhigh,perturb_mode,dpert,&
                        problem_func,sparsity_mode,jacobian_method,jacobian_methods,&
                        class,classes,info,chunk_size,partition_sparsity_pattern,&
                        cache_size,xlow_for_sparsity,xhigh_for_sparsity,&
                        dpert_for_sparsity,sparsity_perturb_mode,&
                        linear_sparsity_tol,function_precision_tol,&
                        num_sparsity_points)

    implicit none

    class(numdiff_type),intent(inout)   :: me
    integer,intent(in)                  :: n               !! number of `x` variables
    integer,intent(in)                  :: m               !! number of `f` functions
    real(wp),dimension(n),intent(in)    :: xlow            !! lower bounds on `x`
    real(wp),dimension(n),intent(in)    :: xhigh           !! upper bounds on `x`
    integer,intent(in)                  :: perturb_mode    !! perturbation mode:
                                                           !! **1** - perturbation is `dx=dpert`,
                                                           !! **2** - perturbation is `dx=dpert*x`,
                                                           !! **3** - perturbation is `dx=dpert*(1+x)`
    real(wp),dimension(n),intent(in)    :: dpert           !! perturbation vector for `x`
    procedure(func)                     :: problem_func    !! the user function that defines the problem
                                                           !! (returns `m` functions)
    integer,intent(in)                  :: sparsity_mode   !! the sparsity computation method:
                                                           !! **1** - assume dense,
                                                           !! **2** - three-point simple method,
                                                           !! **3** - will be specified by the user in
                                                           !! a subsequent call to [[set_sparsity_pattern]].
                                                           !! **4** - computes a two-point jacobian
                                                           !! at `num_sparsity_points` points.
    integer,intent(in),optional         :: jacobian_method !! `id` code for the finite difference method
                                                           !! to use for all `n` variables.
                                                           !! see [[get_finite_difference_method]]
    integer,dimension(n),intent(in),optional :: jacobian_methods !! `id` codes for the finite difference method
                                                                 !! to use for each variable.
                                                                 !! see [[get_finite_difference_method]]
    integer,intent(in),optional :: class  !! method class for the finite difference method
                                          !! to use for all `n` variables.
                                          !! see [[get_finite_difference_method]]
    integer,dimension(n),intent(in),optional :: classes   !! method class for the finite difference methods
                                                          !! to use for each variable.
                                                          !! see [[get_finite_difference_method]]
    procedure(info_f),optional :: info  !! a function the user can define
                                        !! which is called when each column
                                        !! of the jacobian is computed.
                                        !! It can be used to perform any
                                        !! setup operations.
    integer,intent(in),optional :: chunk_size  !! chunk size for allocating the arrays
                                               !! (must be >0) [default is 100]
    logical,intent(in),optional :: partition_sparsity_pattern  !! if the sparisty pattern is to
                                                               !! be partitioned using [[DSM]]
                                                               !! [default is False]
    integer,intent(in),optional :: cache_size    !! if present, this is the cache size
                                                 !! for the function cache
                                                 !! (default is not to enable cache)
    real(wp),dimension(n),intent(in),optional :: xlow_for_sparsity   !! lower bounds on `x` used for
                                                                     !! sparsity computation (when
                                                                     !! `sparsity_mode` is 2). If not
                                                                     !! present, then `xlow` is used.
    real(wp),dimension(n),intent(in),optional :: xhigh_for_sparsity  !! upper bounds on `x` used for
                                                                     !! sparsity computation (when
                                                                     !! `sparsity_mode` is 2). If not
                                                                     !! present, then `xhigh` is used.
    real(wp),dimension(n),intent(in),optional :: dpert_for_sparsity  !! for `sparsity_mode=4`, the perturbation
    integer,intent(in),optional :: sparsity_perturb_mode   !! perturbation mode (required if `sparsity_mode=4`):
                                                           !! **1** - perturbation is `dx=dpert`,
                                                           !! **2** - perturbation is `dx=dpert*x`,
                                                           !! **3** - perturbation is `dx=dpert*(1+x)`
    real(wp),intent(in),optional :: linear_sparsity_tol !! the equality tolerance for derivatives to
                                                        !! indicate a constant jacobian element (linear sparsity)
    real(wp),intent(in),optional :: function_precision_tol  !! the function precision. two functions values
                                                            !! that are the within this tolerance are
                                                            !! considered the same value. This is used
                                                            !! when estimating the sparsity pattern when
                                                            !! `sparsity_mode=2` in [[compute_sparsity_random]]
    integer,intent(in),optional :: num_sparsity_points  !! for `sparsity_mode=4`, the number of jacobian
                                                        !! evaluations used to estimate the sparsity pattern.

    integer :: i      !! counter
    logical :: found  !! flag for [[get_finite_difference_method]]
    logical :: cache  !! if the cache is to be used

    me%use_diff = .false.

    ! cache:
    cache = present(cache_size)
    if (cache) cache = cache_size>0
    if (cache) then
        me%compute_function  => compute_function_with_cache
        call me%cache%initialize(cache_size,n,m)
    else
        me%compute_function  => problem_func
        call me%cache%destroy()
    end if

    ! functions:
    me%problem_func => problem_func

    if (present(partition_sparsity_pattern)) then
        me%partition_sparsity_pattern = partition_sparsity_pattern
    else
        me%partition_sparsity_pattern = .false.
    end if

    ! method:
    if (allocated(me%meth)) deallocate(me%meth)
    if (allocated(me%class)) deallocate(me%class)
    if (allocated(me%class_meths)) deallocate(me%class_meths)

    if (      present(jacobian_method) .and. .not. present(jacobian_methods) .and. &
        .not. present(class) .and. .not. present(classes)) then
        ! use the same for all variable
        me%mode = 1
        allocate(me%meth(n))
        do i=1,n
            call get_finite_difference_method(jacobian_method,me%meth(i),found)
            if (.not. found) then
                call me%raise_exception(8,'initialize_numdiff',&
                                          'invalid jacobian_method')
                return
            end if
        end do
    elseif (.not. present(jacobian_method) .and. present(jacobian_methods) .and. &
            .not. present(class) .and. .not. present(classes)) then
        ! specify a separate method for each variable
        me%mode = 1
        allocate(me%meth(n))
        do i=1,n
            call get_finite_difference_method(jacobian_methods(i),me%meth(i),found)
            if (.not. found) then
                call me%raise_exception(9,'initialize_numdiff',&
                                          'invalid jacobian_methods')
                return
            end if
        end do
        if (me%partition_sparsity_pattern) then
            call me%raise_exception(10,'initialize_numdiff',&
                                        'when using partitioned sparsity pattern, '//&
                                        'all columns must use the same '//&
                                        'finite diff method.')
            return
        end if
    elseif (.not. present(jacobian_method) .and. .not. present(jacobian_methods) .and. &
                  present(class) .and. .not. present(classes)) then
        ! use the class for all variables
        me%mode = 2
        allocate(me%class(n))
        me%class = class
        allocate(me%class_meths(n))
        me%class_meths(1) = get_all_methods_in_class(class)
        if (n>1) me%class_meths(2:n) = me%class_meths(1)  ! just copy them over
    elseif (.not. present(jacobian_method) .and. .not. present(jacobian_methods) .and. &
            .not. present(class) .and. present(classes)) then
        ! specify a separate class for each variable
        me%mode = 2
        me%class = classes
        allocate(me%class_meths(n))
        do i=1,n
            me%class_meths(i) = get_all_methods_in_class(me%class(i))
        end do
        if (me%partition_sparsity_pattern) then
            call me%raise_exception(11,'initialize_numdiff',&
                                        'when using partitioned sparsity pattern, '//&
                                        'all columns must use the same '//&
                                        'finite diff method.')
        end if
    else
        call me%raise_exception(12,'initialize_numdiff',&
                                    'must specify one of either '//&
                                    'jacobian_method, jacobian_methods, class, or classes.')
    end if

    ! size of the problem:
    me%n = n
    me%m = m

    ! input variable bounds:
    call me%set_numdiff_bounds(xlow,xhigh)

    ! perturbation options:
    select case (perturb_mode)
    case(1:3)
        me%perturb_mode = perturb_mode
    case default
        call me%raise_exception(13,'initialize_numdiff',&
                                   'perturb_mode must be 1, 2, or 3.')
        return
    end select
    if (allocated(me%dpert)) deallocate(me%dpert)
    allocate(me%dpert(n))
    me%dpert = abs(dpert)

    ! optional:
    if (present(info))       me%info_function => info
    if (present(chunk_size)) me%chunk_size = abs(chunk_size)

    ! set the jacobian function, depending on the options:
    if (me%partition_sparsity_pattern) then
        me%jacobian_function => compute_jacobian_partitioned
    else
        me%jacobian_function => compute_jacobian_standard
    end if

    ! set the sparsity function
    call me%set_sparsity_mode(sparsity_mode,xlow_for_sparsity,xhigh_for_sparsity)

    if (present(linear_sparsity_tol))    me%linear_sparsity_tol    = linear_sparsity_tol
    if (present(function_precision_tol)) me%function_precision_tol = function_precision_tol
    if (present(num_sparsity_points))    me%num_sparsity_points    = num_sparsity_points

    if (present(dpert_for_sparsity)) then
        me%dpert_for_sparsity = abs(dpert_for_sparsity)
    else
        me%dpert_for_sparsity = dpert ! use the same dpert as the jacobian
    end if

    if (present(sparsity_perturb_mode)) then
        select case (sparsity_perturb_mode)
        case(1:3)
            me%sparsity_perturb_mode = sparsity_perturb_mode
        case default
            call me%raise_exception(14,'initialize_numdiff',&
                                       'sparsity_perturb_mode must be 1, 2, or 3.')
            return
        end select
    else
        me%sparsity_perturb_mode = perturb_mode ! use the same perturb mode as the jacobian
    end if

    end subroutine initialize_numdiff
!*******************************************************************************

!*******************************************************************************
!>
!  destroy the [[numdiff_type]] class.

    subroutine destroy_numdiff_type(me)

    implicit none

    class(numdiff_type),intent(out) :: me

    end subroutine destroy_numdiff_type
!*******************************************************************************

!*******************************************************************************
!>
!  destroy a [[sparsity_pattern]] type.

    subroutine destroy_sparsity(me)

    implicit none

    class(sparsity_pattern),intent(out) :: me

    end subroutine destroy_sparsity
!*******************************************************************************

!*******************************************************************************
!>
!  Wrapper for [[dsm]] to compute the sparsity pattern partition.

    subroutine dsm_wrapper(me,n,m,info)

    implicit none

    class(sparsity_pattern),intent(inout) :: me
    integer,intent(in)  :: n     !! number of columns of jacobian matrix
    integer,intent(in)  :: m     !! number of rows of jacobian matrix
    integer,intent(out) :: info  !! status output from [[dsm]]

    integer :: mingrp !! for call to [[dsm]]
    integer,dimension(:),allocatable :: ipntr  !! for call to [[dsm]]
    integer,dimension(:),allocatable :: jpntr  !! for call to [[dsm]]
    integer,dimension(:),allocatable :: irow   !! for call to [[dsm]]
                                               !! (temp copy since [[dsm]]
                                               !! will modify it)
    integer,dimension(:),allocatable :: icol   !! for call to [[dsm]]
                                               !! (temp copy since [[dsm]]
                                               !! will modify it)

    allocate(ipntr(m+1))
    allocate(jpntr(n+1))
    allocate(me%ngrp(n))
    irow = me%irow
    icol = me%icol

    call dsm(m,n,me%num_nonzero_elements,&
             irow,icol,&
             me%ngrp,me%maxgrp,&
             mingrp,info,ipntr,jpntr)

    end subroutine dsm_wrapper
!*******************************************************************************

!*******************************************************************************
!>
!  Returns the columns in a sparsity partition group.
!
!@note This is just a wrapper to get data from `ngrp`.

    subroutine columns_in_partition_group(me,igroup,n_cols,cols,nonzero_rows,indices,status_ok)

    implicit none

    class(sparsity_pattern),intent(in)           :: me
    integer,intent(in)                           :: igroup       !! group number. Should be `>0` and `<=me%mxgrp`
    integer,intent(out)                          :: n_cols       !! number of columns in the `igroup` group.
    integer,dimension(:),allocatable,intent(out) :: cols         !! the column numbers in the `igroup` group.
                                                                 !! (if none, then it is not allocated)
    integer,dimension(:),allocatable,intent(out) :: nonzero_rows !! the row numbers of all the nonzero
                                                                 !! Jacobian elements in this group
    integer,dimension(:),allocatable,intent(out) :: indices      !! nonzero indices in `jac` for a group
    logical,intent(out)                          :: status_ok    !! true if the partition is valid

    integer :: i  !! counter
    integer :: num_nonzero_elements_in_col          !! number of nonzero elements in a column
    integer :: num_nonzero_elements_in_group        !! number of nonzero elements in a group
    integer,dimension(:),allocatable :: col_indices !! nonzero indices in `jac` for a column

    if (me%maxgrp>0 .and. allocated(me%ngrp)) then

        status_ok = .true.

        n_cols = count(me%ngrp==igroup)
        if (n_cols>0) then
            allocate(cols(n_cols))
            cols = pack([(i,i=1,size(me%ngrp))],mask=me%ngrp==igroup)
        end if

        ! get all the non-zero elements in each column:
        num_nonzero_elements_in_group = 0  ! initialize
        do i = 1, n_cols
            num_nonzero_elements_in_col = count(me%icol==cols(i))
            if (num_nonzero_elements_in_col/=0) then ! there are functions to
                                                     ! compute in this column
                num_nonzero_elements_in_group = num_nonzero_elements_in_group + &
                                                num_nonzero_elements_in_col
                if (allocated(col_indices)) deallocate(col_indices)
                allocate(col_indices(num_nonzero_elements_in_col))
                col_indices = pack(me%indices,mask=me%icol==cols(i))
                if (allocated(nonzero_rows)) then
                    nonzero_rows = [nonzero_rows,me%irow(col_indices)]
                    indices = [indices,col_indices]
                else
                    nonzero_rows = me%irow(col_indices)
                    indices = col_indices
                end if
            end if
        end do

    else
        status_ok = .false.
    end if

    end subroutine columns_in_partition_group
!*******************************************************************************

!*******************************************************************************
!>
!  Destroy the sparsity pattern in the class.

    subroutine destroy_sparsity_pattern(me)

    implicit none

    class(numdiff_type),intent(inout) :: me

    call me%sparsity%destroy()

    end subroutine destroy_sparsity_pattern
!*******************************************************************************

!*******************************************************************************
!>
!  Computes the `indices` vector in the class.

    subroutine compute_indices(me)

    implicit none

    class(sparsity_pattern),intent(inout) :: me

    integer :: i !! counter

    allocate(me%indices(me%num_nonzero_elements))
    me%indices = [(i,i=1,me%num_nonzero_elements)]

    end subroutine compute_indices
!*******************************************************************************

!*******************************************************************************
!>
!  To specify the sparsity pattern directly if it is already known.
!
!@note If specifying the linear pattern, all three optional arguments
!      must be present.

    subroutine set_sparsity_pattern(me,irow,icol,linear_irow,linear_icol,linear_vals,maxgrp,ngrp)

    implicit none

    class(numdiff_type),intent(inout)         :: me
    integer,dimension(:),intent(in)           :: irow        !! sparsity pattern nonzero elements row indices
    integer,dimension(:),intent(in)           :: icol        !! sparsity pattern nonzero elements column indices
    integer,dimension(:),intent(in),optional  :: linear_irow !! linear sparsity pattern nonzero elements row indices
    integer,dimension(:),intent(in),optional  :: linear_icol !! linear sparsity pattern nonzero elements column indices
    real(wp),dimension(:),intent(in),optional :: linear_vals !! linear sparsity values (constant elements of the Jacobian)
    integer,intent(in),optional               :: maxgrp      !! DSM sparsity partition
                                                             !! [only used if `me%partition_sparsity_pattern=True`]
    integer,dimension(:),intent(in),optional  :: ngrp        !! DSM sparsity partition (size `n`)
                                                             !! [only used if `me%partition_sparsity_pattern=True`]

    integer :: info !! status output form [[dsm]]

    call me%destroy_sparsity_pattern()

    if (me%exception_raised) return ! check for exceptions

    if (size(irow)/=size(icol) .or. any(irow>me%m) .or. any(icol>me%n)) then
        call me%raise_exception(15,'set_sparsity_pattern',&
                                   'invalid inputs')
        return
    else

        me%sparsity%sparsity_computed = .true.
        me%sparsity%num_nonzero_elements = size(irow)
        me%sparsity%irow = irow
        me%sparsity%icol = icol

        call me%sparsity%compute_indices()
        if (me%partition_sparsity_pattern) then
            if (present(maxgrp) .and. present(ngrp)) then
                ! use the user-input partition:
                if (maxgrp>0 .and. all(ngrp>=1 .and. ngrp<=maxgrp) .and. size(ngrp)==me%n) then
                    me%sparsity%maxgrp = maxgrp
                    me%sparsity%ngrp   = ngrp
                else
                    call me%raise_exception(28,'set_sparsity_pattern',&
                                            'invalid sparsity partition inputs.')
                    return
                end if
            else
                call me%sparsity%dsm_wrapper(me%n,me%m,info)
                if (info/=1) then
                    call me%raise_exception(16,'set_sparsity_pattern',&
                                            'error partitioning sparsity pattern.')
                    return
                end if
            end if
        end if

    end if

    ! linear pattern:
    if (present(linear_irow) .and. present(linear_icol) .and. present(linear_vals)) then
        if (size(linear_irow)/=size(linear_icol) .or. &
            size(linear_vals)/=size(linear_icol) .or. &
            any(linear_irow>me%m) .or. &
            any(linear_icol>me%n)) then
            call me%raise_exception(17,'set_sparsity_pattern',&
                                       'invalid linear sparsity pattern')
            return
        else
            me%sparsity%linear_irow = linear_irow
            me%sparsity%linear_icol = linear_icol
            me%sparsity%linear_vals = linear_vals
            me%sparsity%linear_sparsity_computed = .true.
            me%sparsity%num_nonzero_linear_elements = size(linear_irow)
        end if
    end if

    end subroutine set_sparsity_pattern
!*******************************************************************************

!*******************************************************************************
!>
!  assume all elements of Jacobian are non-zero.

    subroutine compute_sparsity_dense(me,x)

    implicit none

    class(numdiff_type),intent(inout) :: me
    real(wp),dimension(:),intent(in)  :: x   !! vector of variables (size `n`)

    integer :: i !! counter
    integer :: r !! row counter
    integer :: c !! column counter

    call me%destroy_sparsity_pattern()

    if (me%exception_raised) return ! check for exceptions

    me%sparsity%num_nonzero_elements = me%m * me%n
    allocate(me%sparsity%irow(me%sparsity%num_nonzero_elements))
    allocate(me%sparsity%icol(me%sparsity%num_nonzero_elements))

    call me%sparsity%compute_indices()

    ! create the dense matrix:
    i = 0
    do c = 1, me%n
        do r = 1, me%m
            i = i + 1
            me%sparsity%irow(i) = r
            me%sparsity%icol(i) = c
        end do
    end do

    ! No real need for this, since it can't be partitioned (all elements are true)
    if (me%partition_sparsity_pattern) call me%generate_dense_sparsity_partition()

    me%sparsity%sparsity_computed = .true.

    end subroutine compute_sparsity_dense
!*******************************************************************************

!*******************************************************************************
!>
!  Generate a "dense" sparsity partition.

    subroutine generate_dense_sparsity_partition(me)

    implicit none

    class(numdiff_type),intent(inout) :: me

    integer :: i !! counter

    if (me%exception_raised) return ! check for exceptions

    me%sparsity%maxgrp = me%n
    allocate(me%sparsity%ngrp(me%n))
    me%sparsity%ngrp = [(i, i=1,me%n)]

    end subroutine generate_dense_sparsity_partition
!*******************************************************************************

!*******************************************************************************
!>
!  Compute the sparsity pattern by computing the function at three
!  "random" points in the [`xlow_for_sparsity`, `xhigh_for_sparsity`] interval
!  and checking if the function values are the same.

    subroutine compute_sparsity_random(me,x)

    implicit none

    class(numdiff_type),intent(inout) :: me
    real(wp),dimension(:),intent(in) :: x !! vector of variables (size `n`)
                                          !! (not used here)

    integer :: i !! column counter
    integer :: j !! row counter
    integer :: n_icol  !! `icol` size counter
    integer :: n_irow  !! `irow` size counter
    integer :: n_linear_icol  !! `linear_icol` size counter
    integer :: n_linear_irow  !! `linear_irow` size counter
    integer :: n_linear_vals  !! `linear_vals` size counter
    integer,dimension(me%m) :: idx !! indices to compute `[1,2,...,m]`
    real(wp),dimension(me%n) :: x1 !! perturbed variable vector
    real(wp),dimension(me%n) :: x2 !! perturbed variable vector
    real(wp),dimension(me%n) :: x3 !! perturbed variable vector
    real(wp),dimension(me%n) :: x4 !! perturbed variable vector
    real(wp),dimension(me%m) :: f1 !! function evaluation
    real(wp),dimension(me%m) :: f2 !! function evaluation
    real(wp),dimension(me%m) :: f3 !! function evaluation
    real(wp),dimension(me%m) :: f4 !! function evaluation
    real(wp) :: dfdx1 !! for linear sparsity estimation
    real(wp) :: dfdx2 !! for linear sparsity estimation
    real(wp) :: dfdx3 !! for linear sparsity estimation
    real(wp) :: dfdx  !! for linear sparsity estimation
    integer :: info !! status output form [[dsm]]

    real(wp),dimension(4),parameter :: coeffs = [0.20123456787654321_wp,&
                                                 0.40123456787654321_wp,&
                                                 0.60123456787654321_wp,&
                                                 0.80123456787654321_wp]
        !! Pick three pseudo-random points roughly equally spaced.
        !! (add some noise in attempt to avoid freak zeros)
        !!````
        !! xlow---|----|--x--|---xhigh
        !!        1    2     3
        !!````
        !!
        !! Also using an extra point to estimate the
        !! constant elements of the Jacobian.

    ! initialize:
    call me%destroy_sparsity_pattern()

    if (me%exception_raised) return ! check for exceptions

    ! error check:
    if (.not. allocated(me%xlow_for_sparsity) .or. .not. allocated(me%xhigh_for_sparsity)) then
        call me%raise_exception(18,'compute_sparsity_random',&
                                   'the x bounds have not been set.')
        return
    end if

    associate (xlow => me%xlow_for_sparsity, xhigh => me%xhigh_for_sparsity)

        ! we will compute all the functions:
        idx = [(i,i=1,me%m)]

        n_icol = 0  ! initialize vector size counters
        n_irow = 0
        n_linear_icol = 0
        n_linear_irow = 0
        n_linear_vals = 0

        ! define a nominal point roughly in the middle:
        x2 = xlow + (xhigh-xlow)*coeffs(2)
        call me%compute_function(x2,f2,idx)
        if (me%exception_raised) return ! check for exceptions

        if (me%compute_linear_sparsity_pattern) then
            ! we need another point where we perturb all the variables
            ! to check to make sure it is linear in only one variable.
            x4 = xlow + (xhigh-xlow)*coeffs(4)
            call me%compute_function(x4,f4,idx)
            if (me%exception_raised) return ! check for exceptions
        end if

        do i = 1, me%n  ! columns of Jacobian

            ! restore nominal:
            x1 = x2
            x3 = x2

            x1(i) = xlow(i) + (xhigh(i)-xlow(i))*coeffs(1)
            x3(i) = xlow(i) + (xhigh(i)-xlow(i))*coeffs(3)

            call me%compute_function(x1,f1,idx)
            if (me%exception_raised) return ! check for exceptions
            call me%compute_function(x3,f3,idx)
            if (me%exception_raised) return ! check for exceptions

            do j = 1, me%m ! each function (rows of Jacobian)
                if (equal_within_tol([f1(j),f2(j),f3(j)], me%function_precision_tol)) then
                    ! no change in the function, so no sparsity element here.
                    cycle
                else
                    if (me%compute_linear_sparsity_pattern) then
                        ! if computing the linear pattern separately.
                        ! do the three points lie on a line? (x1,f1) -> (x2, f2), -> (x3,f3)
                        ! [check that the slopes are equal within a tolerance]
                        dfdx1 = (f1(j)-f2(j)) / (x1(i)-x2(i)) ! slope of line from 1->2
                        dfdx2 = (f1(j)-f3(j)) / (x1(i)-x3(i)) ! slope of line from 1->3
                        dfdx3 = (f1(j)-f4(j)) / (x1(i)-x4(i)) ! slope of line from 1->4
                        if (equal_within_tol([dfdx1,dfdx2,dfdx3],me%linear_sparsity_tol)) then
                            ! this is a linear element (constant value)
                            dfdx = (dfdx1 + dfdx2 + dfdx3) / 3.0_wp ! just take the average and use that
                            call expand_vector(me%sparsity%linear_icol,n_linear_icol,me%chunk_size,val=i)
                            call expand_vector(me%sparsity%linear_irow,n_linear_irow,me%chunk_size,val=j)
                            call expand_vector(me%sparsity%linear_vals,n_linear_vals,me%chunk_size,val=dfdx)
                            cycle ! this element will not be added to the nonlinear pattern
                        end if
                    end if
                    ! otherwise, add it to the nonlinear pattern:
                    call expand_vector(me%sparsity%icol,n_icol,me%chunk_size,val=i)
                    call expand_vector(me%sparsity%irow,n_irow,me%chunk_size,val=j)
                end if
            end do

        end do

        ! resize to correct size:
        call me%resize_sparsity_vectors(n_icol,n_irow,n_linear_icol,&
                                            n_linear_irow,n_linear_vals)

    end associate

    call me%sparsity%compute_indices()
    if (me%partition_sparsity_pattern) then
        call me%sparsity%dsm_wrapper(me%n,me%m,info)
        if (info/=1) then
            call me%raise_exception(19,'compute_sparsity_random',&
                                       'error partitioning sparsity pattern.')
            return
        end if
    end if

    ! finished:
    me%sparsity%sparsity_computed = .true.

    end subroutine compute_sparsity_random
!*******************************************************************************

!*******************************************************************************
!>
!  Resize the sparsity arrays after accumulating them.

    subroutine resize_sparsity_vectors(me,n_icol,n_irow,n_linear_icol,&
                                            n_linear_irow,n_linear_vals)

    implicit none

    class(numdiff_type),intent(inout) :: me
    integer,intent(inout) :: n_icol
    integer,intent(inout) :: n_irow
    integer,intent(inout) :: n_linear_icol
    integer,intent(inout) :: n_linear_irow
    integer,intent(inout) :: n_linear_vals

    if (me%exception_raised) return ! check for exceptions

    ! resize to correct size:
    if (allocated(me%sparsity%icol)) then
        call expand_vector(me%sparsity%icol,n_icol,me%chunk_size,finished=.true.)
        call expand_vector(me%sparsity%irow,n_irow,me%chunk_size,finished=.true.)
        me%sparsity%num_nonzero_elements = size(me%sparsity%irow)
    else
        me%sparsity%num_nonzero_elements = 0
    end if

    ! linear pattern (note: there may be no linear elements):
    me%sparsity%linear_sparsity_computed = me%compute_linear_sparsity_pattern .and. &
                                            allocated(me%sparsity%linear_vals)
    if (me%sparsity%linear_sparsity_computed) then
        call expand_vector(me%sparsity%linear_icol,n_linear_icol,me%chunk_size,finished=.true.)
        call expand_vector(me%sparsity%linear_irow,n_linear_irow,me%chunk_size,finished=.true.)
        call expand_vector(me%sparsity%linear_vals,n_linear_vals,me%chunk_size,finished=.true.)
        me%sparsity%num_nonzero_linear_elements = n_linear_vals
    else
        me%sparsity%num_nonzero_linear_elements = 0
    end if

    end subroutine resize_sparsity_vectors
!*******************************************************************************

!*******************************************************************************
!>
!  Compute the sparsity pattern by computing a 2-point jacobian at a specified
!  number of "random" points (`num_sparsity_points`) in the
!  [`xlow_for_sparsity`, `xhigh_for_sparsity`] interval and checking if
!  they are the same.

    subroutine compute_sparsity_random_2(me,x)

    implicit none

    class(numdiff_type),intent(inout) :: me
    real(wp),dimension(:),intent(in)  :: x !! vector of variables (size `n`)
                                           !! (not used here)

    type :: jac_type
        !! so we can define an array of jacobian columns
        real(wp),dimension(:),allocatable :: jac !! a column of the jacobian
    end type jac_type

    real(wp),dimension(:),allocatable :: coeffs           !! coefficients for computing `xp`
    real(wp),dimension(:,:),allocatable :: xp             !! perturbed `x` array (size `n,num_sparsity_points`)
    type(jac_type),dimension(:),allocatable :: jac_array  !! array of jacobian columns
    type(sparsity_pattern) :: tmp_sparsity_pattern        !! will accumulate the pattern here
                                                          !! and copy it into the class when
                                                          !! finished.
    integer :: i              !! counter for the number of points
    integer :: j              !! counter
    integer :: icol           !! column counter
    integer :: irow           !! row counter
    real(wp) :: dfdx          !! value of a jacobian element
    integer :: n_icol         !! `icol` size counter
    integer :: n_irow         !! `irow` size counter
    integer :: n_linear_icol  !! `linear_icol` size counter
    integer :: n_linear_irow  !! `linear_irow` size counter
    integer :: n_linear_vals  !! `linear_vals` size counter
    type(meth_array) :: class_meths  !! set of finite diff methods to use
    real(wp),dimension(:),allocatable :: jac !! array of jacobian element values
    integer :: info !! status output form [[dsm]]

    ! initialize:
    call me%destroy_sparsity_pattern()

    if (me%exception_raised) return ! check for exceptions

    ! error check:
    if (.not. allocated(me%xlow_for_sparsity) .or. .not. allocated(me%xhigh_for_sparsity)) then
        call me%raise_exception(20,'compute_sparsity_random_2',&
                                   'the x bounds have not been set.')
        return
    end if

    ! compute the coefficients:
    coeffs = divide_interval(me%num_sparsity_points)

    ! save the initial settings:
    tmp_sparsity_pattern = me%sparsity

    ! we create a provisional sparsity pattern here,
    ! to compute one row at a time:
    me%sparsity%irow = [(irow,irow=1,me%m)]
    me%sparsity%sparsity_computed = .true.
    me%sparsity%num_nonzero_elements = me%m
    call me%sparsity%compute_indices()
    if (me%partition_sparsity_pattern) call me%generate_dense_sparsity_partition()
    n_icol = 0
    n_irow = 0
    n_linear_icol = 0
    n_linear_irow = 0
    n_linear_vals = 0
    allocate(jac(me%num_sparsity_points))

    ! the idea here is to compute the (dense) jacobian
    ! one column at at time, and keep a running track of
    ! the ones that have changed.

    ! first precompute the perturbation arrays for each point:
    allocate(xp(me%n,me%num_sparsity_points))
    do i = 1, me%num_sparsity_points
        xp(:,i) = me%xlow_for_sparsity + (me%xhigh_for_sparsity-me%xlow_for_sparsity)*coeffs(i)
    end do

    ! we will use 2-point methods (simple differences):
    class_meths = get_all_methods_in_class(2)

    do icol = 1, me%n  ! column loop

        ! size the array:
        if (allocated(jac_array)) deallocate(jac_array)
        allocate(jac_array(me%num_sparsity_points))

        ! compute the ith column of the jacobian:
        me%sparsity%icol = [(icol, j=1,me%m)]
        do i = 1, me%num_sparsity_points
            call me%compute_jacobian_for_sparsity( icol, class_meths, xp(:,i), jac_array(i)%jac )
            if (me%exception_raised) return ! check for exceptions
        end do

        ! check each row:
        do irow = 1, me%m

            ! get the jacobian values for this row,col for all the points:
            do j = 1, me%num_sparsity_points
                jac(j) = jac_array(j)%jac(irow)
            end do

            ! put the results into the tmp_sparsity_pattern
            if (equal_within_tol([0.0_wp,jac],me%linear_sparsity_tol)) then
                ! they are all zero
                cycle
            else
                if (me%compute_linear_sparsity_pattern) then
                    if (equal_within_tol(jac,me%linear_sparsity_tol)) then
                        ! this is a linear element (constant value)
                        dfdx = sum(jac) / me%num_sparsity_points ! just take the average and use that
                        call expand_vector(tmp_sparsity_pattern%linear_icol,n_linear_icol,me%chunk_size,val=icol)
                        call expand_vector(tmp_sparsity_pattern%linear_irow,n_linear_irow,me%chunk_size,val=irow)
                        call expand_vector(tmp_sparsity_pattern%linear_vals,n_linear_vals,me%chunk_size,val=dfdx)
                        cycle ! this element will not be added to the nonlinear pattern
                    end if
                end if
                ! otherwise, add it to the nonlinear pattern:
                call expand_vector(tmp_sparsity_pattern%icol,n_icol,me%chunk_size,val=icol)
                call expand_vector(tmp_sparsity_pattern%irow,n_irow,me%chunk_size,val=irow)
            end if
        end do

    end do

    ! copy over the computed pattern:
    me%sparsity = tmp_sparsity_pattern

    ! resize to correct size:
    call me%resize_sparsity_vectors(n_icol,n_irow,n_linear_icol,&
                                        n_linear_irow,n_linear_vals)

    call me%sparsity%compute_indices()
    if (me%partition_sparsity_pattern) then
        call me%sparsity%dsm_wrapper(me%n,me%m,info)
        if (info/=1) then
            call me%raise_exception(21,'compute_sparsity_random_2',&
                                       'error partitioning sparsity pattern.')
            return
        end if
    end if
    me%sparsity%sparsity_computed = .true.

    end subroutine compute_sparsity_random_2
!*******************************************************************************

!*******************************************************************************
!>
!  Computes the sparsity pattern and return it.
!  Uses the settings currently in the class.

    subroutine compute_sparsity_pattern(me,x,irow,icol,linear_irow,linear_icol,linear_vals)

    implicit none

    class(numdiff_type),intent(inout) :: me
    real(wp),dimension(:),intent(in)  :: x     !! vector of variables (size `n`)
    integer,dimension(:),allocatable,intent(out) :: irow  !! sparsity pattern nonzero elements row indices
    integer,dimension(:),allocatable,intent(out) :: icol  !! sparsity pattern nonzero elements column indices
    integer,dimension(:),allocatable,intent(out),optional  :: linear_irow !! linear sparsity pattern
                                                                          !! nonzero elements row indices
    integer,dimension(:),allocatable,intent(out),optional  :: linear_icol !! linear sparsity pattern nonzero
                                                                          !! elements column indices
    real(wp),dimension(:),allocatable,intent(out),optional :: linear_vals !! linear sparsity values (constant
                                                                          !! elements of the Jacobian)

    if (me%exception_raised) return ! check for exceptions

    if (associated(me%compute_sparsity)) then
        call me%compute_sparsity(x)
        if (me%exception_raised) return ! check for exceptions
        call me%get_sparsity_pattern(irow,icol,linear_irow,linear_icol,linear_vals)
    end if

    end subroutine compute_sparsity_pattern
!*******************************************************************************

!*******************************************************************************
!>
!  Returns the sparsity pattern from the class.
!  If it hasn't been computed, the output arrays will not be allocated.

    subroutine get_sparsity_pattern(me,irow,icol,&
                                    linear_irow,linear_icol,linear_vals,&
                                    maxgrp,ngrp)

    implicit none

    class(numdiff_type),intent(inout) :: me
    integer,dimension(:),allocatable,intent(out)  :: irow  !! sparsity pattern nonzero elements row indices
    integer,dimension(:),allocatable,intent(out)  :: icol  !! sparsity pattern nonzero elements column indices
    integer,dimension(:),allocatable,intent(out),optional  :: linear_irow !! linear sparsity pattern
                                                                          !! nonzero elements row indices
    integer,dimension(:),allocatable,intent(out),optional  :: linear_icol !! linear sparsity pattern nonzero
                                                                          !! elements column indices
    real(wp),dimension(:),allocatable,intent(out),optional :: linear_vals !! linear sparsity values (constant
                                                                          !! elements of the Jacobian)
    integer,intent(out),optional                           :: maxgrp      !! DSM sparsity partition
    integer,dimension(:),allocatable,intent(out),optional  :: ngrp        !! DSM sparsity partition

    if (me%exception_raised) return ! check for exceptions

    if (allocated(me%sparsity%irow) .and. allocated(me%sparsity%icol)) then
        irow = me%sparsity%irow
        icol = me%sparsity%icol
    end if

    ! optional linear pattern output:
    if (present(linear_irow)) then
        if (allocated(me%sparsity%linear_irow)) linear_irow = me%sparsity%linear_irow
    end if
    if (present(linear_icol)) then
        if (allocated(me%sparsity%linear_icol)) linear_icol = me%sparsity%linear_icol
    end if
    if (present(linear_vals)) then
        if (allocated(me%sparsity%linear_vals)) linear_vals = me%sparsity%linear_vals
    end if

    ! optional DSM partition:
    if (present(ngrp) .and. allocated(me%sparsity%ngrp)) ngrp = me%sparsity%ngrp
    if (present(maxgrp)) maxgrp = me%sparsity%maxgrp

    end subroutine get_sparsity_pattern
!*******************************************************************************

!*******************************************************************************
!>
!  just a wrapper for [[compute_jacobian]], that returns a dense (`m x n`) matrix.
!
!@note This one will include the constant elements if the linear pattern is available.

    subroutine compute_jacobian_dense(me,x,jac)

    implicit none

    class(numdiff_type),intent(inout) :: me
    real(wp),dimension(:),intent(in) :: x !! vector of variables (size `n`)
    real(wp),dimension(:,:),allocatable,intent(out) :: jac !! the jacobian matrix

    real(wp),dimension(:),allocatable :: jac_vec  !! sparse jacobian representation
    integer :: i !! counter

    if (me%exception_raised) return ! check for exceptions

    ! size output matrix:
    allocate(jac(me%m,me%n))

    ! convert to dense form:
    jac = zero

    ! compute sparse form of jacobian:
    call me%compute_jacobian(x,jac_vec)
    if (me%exception_raised) return ! check for exceptions

    if (allocated(jac_vec)) then
        ! add the nonlinear elements:
        do i = 1, me%sparsity%num_nonzero_elements
            jac(me%sparsity%irow(i),me%sparsity%icol(i)) = jac_vec(i)
        end do
        deallocate(jac_vec)
    end if

    ! add the constant elements if necessary:
    do i = 1, me%sparsity%num_nonzero_linear_elements
        jac(me%sparsity%linear_irow(i),me%sparsity%linear_icol(i)) = me%sparsity%linear_vals(i)
    end do

    end subroutine compute_jacobian_dense
!*******************************************************************************

!*******************************************************************************
!>
!  Perturb the specified optimization variable, and compute the function.
!  This routine is designed so that `df` is accumulated as each function
!  evaluation is done, to avoid having to allocate more temporary storage.

    subroutine perturb_x_and_compute_f(me,x,dx_factor,dx,&
                                       df_factor,column,idx,df)

    implicit none

    class(numdiff_type),intent(inout) :: me
    real(wp),dimension(:),intent(in)  :: x         !! nominal variable vector
    real(wp),intent(in)               :: dx_factor !! factor to multiply `dx`
    real(wp),dimension(:),intent(in)  :: dx        !! the perturbation value for this column
    real(wp),intent(in)               :: df_factor !! factor to multiply function value
    integer,intent(in)                :: column    !! the variable to perturb
    integer,dimension(:),intent(in)   :: idx       !! the elements in this
                                                   !! column of the Jacobian
                                                   !! to compute (passed to function)
    real(wp),dimension(me%m),intent(inout) :: df   !! the accumulated function value
                                                   !! note: for the first call, this
                                                   !! should be set to zero

    real(wp),dimension(me%n) :: xp  !! the perturbed variable vector
    real(wp),dimension(me%m) :: f   !! function evaluation

    if (me%exception_raised) return ! check for exceptions

    xp = x
    if (dx_factor/=zero) xp(column) = xp(column) + dx_factor * dx(column)
    call me%compute_function(xp,f,idx)
    if (me%exception_raised) return ! check for exceptions
    df(idx) = df(idx) + df_factor * f(idx)

    end subroutine perturb_x_and_compute_f
!*******************************************************************************

!*******************************************************************************
!>
!  Returns the product `J*v`, where `J` is the `m x n` Jacobian matrix
!  and `v` is an `n x 1` vector.
!
!@note This one will include the constant elements if the linear pattern is available.

    subroutine compute_jacobian_times_vector(me,x,v,z)

    implicit none

    class(numdiff_type),intent(inout) :: me
    real(wp),dimension(:),intent(in)  :: x    !! vector of variables (size `n`)
    real(wp),dimension(:),intent(in)  :: v    !! vector (size `n`)
    real(wp),dimension(:),intent(out) :: z    !! The product `J*v` (size `m`)

    real(wp),dimension(:),allocatable :: jac !! sparse jacobian vector
    integer :: i !! counter
    integer :: r !! row number in full jacobian
    integer :: c !! column number in full jacobian

    if (me%exception_raised) return ! check for exceptions

    ! first compute the jacobian in sparse vector form:
    call me%compute_jacobian(x,jac)
    if (me%exception_raised) return ! check for exceptions

    ! initialize output vector:
    z = 0.0_wp

    if (allocated(jac)) then

        !   J    v   z
        !  ---   -   -
        !  X0X   X   X
        !  0X0 * X = X
        !  00X   X   X
        !  00X       X

        ! multiplication by input v:
        do i = 1, me%sparsity%num_nonzero_elements
            r = me%sparsity%irow(i)
            c = me%sparsity%icol(i)
            z(r) = z(r) + jac(i)*v(c)
        end do

        deallocate(jac)

    end if

    ! linear elements if available:
    do i = 1, me%sparsity%num_nonzero_linear_elements
        r = me%sparsity%linear_irow(i)
        c = me%sparsity%linear_icol(i)
        z(r) = z(r) + me%sparsity%linear_vals(i)*v(c)
    end do

    end subroutine compute_jacobian_times_vector
!*******************************************************************************

!*******************************************************************************
!>
!  Compute the Jacobian.
!
!@note The output `jac` only includes the elements of the nonlinear Jacobian.
!      If the constant elements are being handled separately (if the linear
!      pattern is available), then those elements can be obtained by
!      calling `get_sparsity_pattern` if required.

    subroutine compute_jacobian(me,x,jac)

    implicit none

    class(numdiff_type),intent(inout)             :: me
    real(wp),dimension(:),intent(in)              :: x    !! vector of variables (size `n`)
    real(wp),dimension(:),allocatable,intent(out) :: jac  !! sparse jacobian vector

    real(wp),dimension(me%n) :: dx  !! absolute perturbation (>0) for each variable

    if (me%exception_raised) return ! check for exceptions

    ! if we don't have a sparsity pattern yet then compute it:
    ! [also computes the indices vector]
    if (.not. me%sparsity%sparsity_computed) call me%compute_sparsity(x)
    if (me%sparsity%num_nonzero_elements==0) return

    ! size the jacobian vector:
    allocate(jac(me%sparsity%num_nonzero_elements))

    ! compute dx vector:
    if (.not. me%use_diff) then
        ! only need this for the finite difference methods (not diff)
        call me%compute_perturbation_vector(x,dx)
        if (me%exception_raised) return ! check for exceptions
    end if

    ! compute the jacobian:
    if (associated(me%jacobian_function)) then
        call me%jacobian_function(x,dx,jac)
    else
        call me%raise_exception(22,'compute_jacobian',&
                                   'jacobian_function has not been associated.')
        return
    end if

    end subroutine compute_jacobian
!*******************************************************************************

!*******************************************************************************
!>
!  A separate version of [[compute_jacobian]] to be used only when
!  computing the sparsity pattern in [[compute_sparsity_random_2]].
!  It uses `class_meths` and the sparsity dperts and bounds.
!
!@note Based on [[compute_jacobian]]. The index manipulation here could be
!      greatly simplified, since we realdy know we are computed all the
!      elements in one column.

    subroutine compute_jacobian_for_sparsity(me,i,class_meths,x,jac)

    implicit none

    class(numdiff_type),intent(inout)             :: me
    integer,intent(in)                            :: i           !! the column being computed
    type(meth_array),intent(in)                   :: class_meths !! set of finite diff methods to use
    real(wp),dimension(:),intent(in)              :: x           !! vector of variables (size `n`)
    real(wp),dimension(:),allocatable,intent(out) :: jac         !! sparse jacobian vector

    real(wp),dimension(me%n) :: dx  !! absolute perturbation (>0) for each variable
    integer,dimension(:),allocatable :: nonzero_elements_in_col  !! the indices of the
                                                                 !! nonzero Jacobian
                                                                 !! elements in a column
    integer                  :: j   !! function evaluation counter
    real(wp),dimension(me%m) :: df  !! accumulated function
    type(finite_diff_method) :: fd  !! a finite different method (when
                                    !! specifying class rather than the method)
    logical :: status_ok                    !! error flag
    integer :: num_nonzero_elements_in_col  !! number of nonzero elements in a column

    if (me%exception_raised) return ! check for exceptions

    ! Note that a sparsity pattern has already been set

    ! size the jacobian vector:
    allocate(jac(me%sparsity%num_nonzero_elements))

    ! compute the perturbation vector (really we only need dx(i)):
    call me%compute_sparsity_perturbation_vector(x,dx)
    if (me%exception_raised) return ! check for exceptions

    ! initialize:
    jac = zero

    ! determine functions to compute for this column:
    num_nonzero_elements_in_col = count(me%sparsity%icol==i)
    if (num_nonzero_elements_in_col/=0) then ! there are functions to compute

        nonzero_elements_in_col = pack(me%sparsity%irow,mask=me%sparsity%icol==i)

        call me%select_finite_diff_method(x(i),me%xlow_for_sparsity(i),me%xhigh_for_sparsity(i),&
                                            dx(i),class_meths,fd,status_ok)
        if (.not. status_ok) then
            if (me%print_messages) then
                write(error_unit,'(A,1X,I5)') &
                'Error in compute_jacobian_for_sparsity: variable bounds violated for column: ',i
            end if
        end if

        ! compute this column of the Jacobian:
        df = zero
        do j = 1, size(fd%dx_factors)
            if (associated(me%info_function)) call me%info_function([i],j,x)
            call me%perturb_x_and_compute_f(x,fd%dx_factors(j),&
                                            dx,fd%df_factors(j),&
                                            i,nonzero_elements_in_col,df)
            if (me%exception_raised) return ! check for exceptions
        end do
        df(nonzero_elements_in_col) = df(nonzero_elements_in_col) / &
                                        (fd%df_den_factor*dx(i))

        ! put result into the output vector:
        jac(pack(me%sparsity%indices,mask=me%sparsity%icol==i)) = &
            df(nonzero_elements_in_col)

    end if

    end subroutine compute_jacobian_for_sparsity
!*******************************************************************************

!*******************************************************************************
!>
!  Compute the Jacobian using finite differences.
!  (one column at a time)

    subroutine compute_jacobian_standard(me,x,dx,jac)

    implicit none

    class(numdiff_type),intent(inout)   :: me
    real(wp),dimension(:),intent(in)    :: x    !! vector of variables (size `n`)
    real(wp),dimension(me%n),intent(in) :: dx   !! absolute perturbation (>0)
                                                !! for each variable
    real(wp),dimension(:),intent(out)   :: jac  !! sparse jacobian vector (size
                                                !! `num_nonzero_elements`)

    integer,dimension(:),allocatable :: nonzero_elements_in_col  !! the indices of the
                                                                 !! nonzero Jacobian
                                                                 !! elements in a column
    integer                  :: i   !! column counter
    integer                  :: j   !! function evaluation counter
    real(wp),dimension(me%m) :: df  !! accumulated function
    type(finite_diff_method) :: fd  !! a finite different method (when
                                    !! specifying class rather than the method)
    logical :: status_ok   !! error flag
    integer :: num_nonzero_elements_in_col  !! number of nonzero elements in a column

    if (me%exception_raised) return ! check for exceptions

    ! initialize:
    jac = zero

    ! compute Jacobian matrix column-by-column:
    do i=1,me%n

        ! determine functions to compute for this column:
        num_nonzero_elements_in_col = count(me%sparsity%icol==i)
        if (num_nonzero_elements_in_col/=0) then ! there are functions to compute

            nonzero_elements_in_col = pack(me%sparsity%irow,mask=me%sparsity%icol==i)

            select case (me%mode)
            case(1) ! use the specified methods

                ! compute this column of the Jacobian:
                df = zero
                do j = 1, size(me%meth(i)%dx_factors)
                    if (associated(me%info_function)) call me%info_function([i],j,x)
                    call me%perturb_x_and_compute_f(x,me%meth(i)%dx_factors(j),&
                                                    dx,me%meth(i)%df_factors(j),&
                                                    i,nonzero_elements_in_col,df)
                    if (me%exception_raised) return ! check for exceptions
                end do

                df(nonzero_elements_in_col) = df(nonzero_elements_in_col) / &
                                              (me%meth(i)%df_den_factor*dx(i))

            case(2) ! select the method from the class so as not to violate the bounds

                call me%select_finite_diff_method(x(i),me%xlow(i),me%xhigh(i),&
                                                  dx(i),me%class_meths(i),fd,status_ok)
                if (.not. status_ok) then
                    if (me%print_messages) then
                        write(error_unit,'(A,1X,I5)') &
                        'Error in compute_jacobian_standard: variable bounds violated for column: ',i
                    end if
                end if

                ! compute this column of the Jacobian:
                df = zero
                do j = 1, size(fd%dx_factors)
                    if (associated(me%info_function)) call me%info_function([i],j,x)
                    call me%perturb_x_and_compute_f(x,fd%dx_factors(j),&
                                                    dx,fd%df_factors(j),&
                                                    i,nonzero_elements_in_col,df)
                    if (me%exception_raised) return ! check for exceptions
                end do
                df(nonzero_elements_in_col) = df(nonzero_elements_in_col) / &
                                              (fd%df_den_factor*dx(i))

            case default
                call me%raise_exception(23,'compute_jacobian_standard',&
                                           'invalid mode')
                return
            end select

            ! put result into the output vector:
            jac(pack(me%sparsity%indices,mask=me%sparsity%icol==i)) = &
                df(nonzero_elements_in_col)

        end if

    end do

    end subroutine compute_jacobian_standard
!*******************************************************************************

!*******************************************************************************
!>
!  Compute the Jacobian one element at a time using the Neville's process
!  algorithm [[diff]]. This takes a very large number of function evaluations,
!  but should give a very accurate answer.

    subroutine compute_jacobian_with_diff(me,x,dx,jac)

    implicit none

    class(numdiff_type),intent(inout)   :: me
    real(wp),dimension(:),intent(in)    :: x    !! vector of variables (size `n`)
    real(wp),dimension(me%n),intent(in) :: dx   !! absolute perturbation (>0)
                                                !! for each variable
    real(wp),dimension(:),intent(out)   :: jac  !! sparse jacobian vector (size
                                                !! `num_nonzero_elements`)

    integer,parameter  :: iord  = 1  !! tells [[diff]] to compute first derivative

    integer                  :: i      !! counter for nonzero elements of jacobian
    type(diff_func)          :: d      !! the [[diff]] class to use
    real(wp)                 :: x0     !! value of ith variable for [[diff]]
    real(wp)                 :: xmin   !! ith variable lower bound
    real(wp)                 :: xmax   !! ith variable upper bound
    real(wp)                 :: deriv  !! derivative value df/dx from [[diff]]
    real(wp)                 :: error  !! estimated error from [[diff]]
    integer                  :: ifail  !! [[diff]] error flag
    real(wp),dimension(me%n) :: xp     !! used by [[dfunc]].
    real(wp),dimension(me%m) :: fvec   !! used by [[dfunc]].
    integer                  :: ir     !! row index of next non-zero element
    integer                  :: ic     !! column index of next non-zero element

    integer :: ic_prev  !! previous column perturbed
    integer :: icount   !! count of number of times a column has been perturbed
    logical :: use_info !! if we are reporting to the user

    if (me%exception_raised) return ! check for exceptions

    ! initialize:
    jac = zero
    use_info = associated(me%info_function)
    if (use_info) ic_prev = -1

    ! set the function for diff:
    call d%set_function(dfunc)

    ! each element is computed one by one
    do i = 1, me%sparsity%num_nonzero_elements

        ir   = me%sparsity%irow(i)
        ic   = me%sparsity%icol(i)
        x0   = x(ic)
        xmin = me%xlow(ic)
        xmax = me%xhigh(ic)

        ! if reporting to the user, have to keep track
        ! of which column is being perturbed, and how
        ! many times:
        if (use_info) then
            if (ic/=ic_prev) icount = 0
        end if

        call d%compute_derivative(iord,x0,xmin,xmax,me%eps,me%acc,deriv,error,ifail)

        if (ifail==0 .or. ifail==1) then
            jac(i) = deriv
        else if (ifail == -1) then
            ! this indicates an exception was
            ! already raised, so just return.
            return
        else
            call me%raise_exception(24, 'compute_jacobian_with_diff',&
                                        'error computing derivative with DIFF.')
            return
        end if

        if (use_info) ic_prev = ic

    end do

    contains

        function dfunc(this,xval) result(fx)

        !! interface to user function for [[diff]] routine.

        implicit none

        class(diff_func),intent(inout)  :: this
        real(wp),intent(in)             :: xval  !! input variable (`ic` variable)
        real(wp)                        :: fx    !! derivative of `ir` function
                                                 !! w.r.t. `xval` variable

        if (use_info) then
            icount = icount + 1
            call me%info_function([ic],icount,x)
        end if

        xp = x
        xp(ic) = xval
        call me%compute_function(xp,fvec,funcs_to_compute=[ir])

        if (me%exception_raised) then ! check for exceptions
            fx = 0.0_wp
            call me%terminate()
        else
            fx = fvec(ir)
        end if

        end function dfunc

    end subroutine compute_jacobian_with_diff
!*******************************************************************************

!*******************************************************************************
!>
!  Compute the Jacobian using finite differences,
!  (using the partitioned sparsity pattern to compute multiple columns
!  at a time).

    subroutine compute_jacobian_partitioned(me,x,dx,jac)

    implicit none

    class(numdiff_type),intent(inout)    :: me
    real(wp),dimension(:),intent(in)     :: x    !! vector of variables (size `n`)
    real(wp),dimension(me%n),intent(in)  :: dx   !! absolute perturbation (>0) for each variable
    real(wp),dimension(:),intent(out)    :: jac  !! sparse jacobian vector

    integer                          :: i             !! column counter
    integer                          :: j             !! function evaluation counter
    integer                          :: igroup        !! group number counter
    integer                          :: n_cols        !! number of columns in a group
    integer,dimension(:),allocatable :: cols          !! array of column indices in a group
    integer,dimension(:),allocatable :: nonzero_rows  !! the indices of the nonzero Jacobian
                                                      !! elementes (row numbers) in a group
    integer,dimension(:),allocatable :: indices       !! nonzero indices in `jac` for a group
    integer,dimension(:),allocatable :: col_indices   !! nonzero indices in `jac` for a column
    real(wp),dimension(me%m)         :: df            !! accumulated function
    type(finite_diff_method)         :: fd            !! a finite different method (when
                                                      !! specifying class rather than the method)
    logical                          :: status_ok     !! error flag
    integer                          :: num_nonzero_elements_in_col

    if (me%exception_raised) return ! check for exceptions

    ! initialize:
    jac = zero

    ! compute by group:
    do igroup = 1, me%sparsity%maxgrp

        ! get the columns in this group:
        call me%sparsity%columns_in_partition_group(igroup,n_cols,cols,&
                                                    nonzero_rows,indices,status_ok)
        if (.not. status_ok) then
            call me%raise_exception(25,'compute_jacobian_partitioned',&
                                       'the partition has not been computed.')
            return
        end if

        if (n_cols>0) then

            if (allocated(nonzero_rows)) then

                select case (me%mode)
                case(1) ! use the specified methods

                    ! note: all the methods must be the same within a group

                    ! compute the columns of the Jacobian in this group:
                    df = zero
                    do j = 1, size(me%meth(1)%dx_factors)
                         if (associated(me%info_function)) call me%info_function(cols,j,x)
                         call me%perturb_x_and_compute_f_partitioned(x,me%meth(1)%dx_factors(j),&
                                                         dx,me%meth(1)%df_factors(j),&
                                                         cols,nonzero_rows,df)
                         if (me%exception_raised) return ! check for exceptions
                    end do
                    ! divide by the denominator, which can be different for each column:
                    do i = 1, n_cols
                        num_nonzero_elements_in_col = count(me%sparsity%icol==cols(i))
                        if (allocated(col_indices)) deallocate(col_indices)
                        allocate(col_indices(num_nonzero_elements_in_col))
                        col_indices = pack(me%sparsity%indices,mask=me%sparsity%icol==cols(i))
                        df(me%sparsity%irow(col_indices)) = df(me%sparsity%irow(col_indices)) / &
                                                            (me%meth(1)%df_den_factor*dx(cols(i)))
                    end do

                case(2) ! select the method from the class so as not to violate
                        ! the bounds on *any* of the variables in the group

                    ! note: all the classes must be the same within a group

                    call me%select_finite_diff_method_for_partition_group( &
                                    x(cols),me%xlow(cols),me%xhigh(cols),&
                                    dx(cols),me%class_meths(1),fd,status_ok)

                    if (.not. status_ok) then
                        if (me%print_messages) then
                            ! will not consider this a fatal error for now:
                            write(error_unit,'(A,1X,I5,1X,A,1X,*(I5,1X))') &
                                'Error in compute_jacobian_partitioned: '//&
                                'variable bounds violated for group: ',&
                                igroup,'. columns: ',cols
                        end if
                    end if

                    ! compute the columns of the Jacobian in this group:
                    df = zero
                    do j = 1, size(fd%dx_factors)
                        if (associated(me%info_function)) call me%info_function(cols,j,x)
                        call me%perturb_x_and_compute_f_partitioned(x,fd%dx_factors(j),&
                                                        dx,fd%df_factors(j),&
                                                        cols,nonzero_rows,df)
                        if (me%exception_raised) return ! check for exceptions
                    end do
                    ! divide by the denominator, which can be different for each column:
                    do i = 1, n_cols
                        num_nonzero_elements_in_col = count(me%sparsity%icol==cols(i))
                        if (allocated(col_indices)) deallocate(col_indices)
                        allocate(col_indices(num_nonzero_elements_in_col))
                        col_indices = pack(me%sparsity%indices,mask=me%sparsity%icol==cols(i))
                        df(me%sparsity%irow(col_indices)) = df(me%sparsity%irow(col_indices)) / &
                                                            (fd%df_den_factor*dx(cols(i)))
                    end do

                case default
                    call me%raise_exception(26,'compute_jacobian_partitioned','invalid mode')
                    return
                end select

                ! put result into the output vector:
                jac(indices) = df(nonzero_rows)

            end if

        end if

    end do

    end subroutine compute_jacobian_partitioned
!*******************************************************************************

!*******************************************************************************
!>
!  Perturb the specified optimization variable, and compute the function.
!  This routine is designed so that `df` is accumulated as each function
!  evaluation is done, to avoid having to allocate more temporary storage.

    subroutine perturb_x_and_compute_f_partitioned(me,x,dx_factor,dx,&
                                       df_factor,columns,idx,df)

    implicit none

    class(numdiff_type),intent(inout) :: me
    real(wp),dimension(:),intent(in)  :: x         !! nominal variable vector
    real(wp),intent(in)               :: dx_factor !! factor to multiply `dx`
    real(wp),dimension(:),intent(in)  :: dx        !! the perturbation value for this column
    real(wp),intent(in)               :: df_factor !! factor to multiply function value
    integer,dimension(:),intent(in)   :: columns   !! the variables to perturb
    integer,dimension(:),intent(in)   :: idx       !! the elements in this
                                                   !! column of the Jacobian
                                                   !! to compute (passed to function)
    real(wp),dimension(me%m),intent(inout) :: df   !! the accumulated function value
                                                   !! note: for the first call, this
                                                   !! should be set to zero

    real(wp),dimension(me%n) :: xp  !! the perturbed variable vector
    real(wp),dimension(me%m) :: f   !! function evaluation

    if (me%exception_raised) return ! check for exceptions

    xp = x
    if (dx_factor/=zero) xp(columns) = xp(columns) + dx_factor * dx(columns)
    call me%compute_function(xp,f,idx)
    if (me%exception_raised) return ! check for exceptions
    df(idx) = df(idx) + df_factor * f(idx)

    end subroutine perturb_x_and_compute_f_partitioned
!*******************************************************************************

!*******************************************************************************
!>
!  Compute `dx`, the perturbation vector for `x`

    subroutine compute_perturb_vector(me,x,dpert,perturb_mode,dx)

    implicit none

    class(numdiff_type),intent(inout)    :: me
    real(wp),dimension(me%n),intent(in)  :: x     !! vector of variables (size `n`)
    real(wp),dimension(me%n),intent(in)  :: dpert
    integer,intent(in)                   :: perturb_mode
    real(wp),dimension(me%n),intent(out) :: dx  !! absolute perturbation (>0)
                                                !! for each variable

    real(wp),parameter :: eps = epsilon(1.0_wp) !! the smallest allowed absolute step

    if (me%exception_raised) return ! check for exceptions

    select case (perturb_mode)
    case(1)
        dx = abs(dpert)
    case(2)
        dx = abs(dpert * x)
    case(3)
        dx = abs(dpert) * (1.0_wp + abs(x))
    case default
        call me%raise_exception(27,'compute_perturb_vector',&
                                   'invalid value for perturb_mode (must be 1, 2, or 3)')
        return
    end select

    ! make sure none are too small:
    where (dx<eps)
        dx = eps
    end where

    end subroutine compute_perturb_vector
!*******************************************************************************

!*******************************************************************************
!>
!  Compute `dx`, the perturbation vector for `x` used
!  when computing the gradients.

    subroutine compute_perturbation_vector(me,x,dx)

    implicit none

    class(numdiff_type),intent(inout)    :: me
    real(wp),dimension(me%n),intent(in)  :: x  !! vector of variables (size `n`)
    real(wp),dimension(me%n),intent(out) :: dx !! absolute perturbation (>0)
                                               !! for each variable

    if (me%exception_raised) return ! check for exceptions

    call me%compute_perturb_vector(x,me%dpert,me%perturb_mode,dx)

    end subroutine compute_perturbation_vector
!*******************************************************************************

!*******************************************************************************
!>
!  Compute `dx`, the perturbation vector for `x` used
!  when computing the sparsity pattern.

    subroutine compute_sparsity_perturbation_vector(me,x,dx)

    implicit none

    class(numdiff_type),intent(inout)    :: me
    real(wp),dimension(me%n),intent(in)  :: x  !! vector of variables (size `n`)
    real(wp),dimension(me%n),intent(out) :: dx !! absolute perturbation (>0)
                                               !! for each variable

    if (me%exception_raised) return ! check for exceptions

    call me%compute_perturb_vector(x,me%dpert_for_sparsity,me%sparsity_perturb_mode,dx)

    end subroutine compute_sparsity_perturbation_vector
!*******************************************************************************

!*******************************************************************************
!>
!  Print the sparsity pattern.

    subroutine print_sparsity(me,n,m,iunit,dense)

    implicit none

    class(sparsity_pattern),intent(in) :: me
    integer,intent(in) :: n  !! number of variables (columns of jacobian)
    integer,intent(in) :: m  !! number of functions (rows of jacobian)
    integer,intent(in) :: iunit !! file unit to write to.
                                !! (assumed to be already opened)
    logical,intent(in),optional :: dense  !! if present and true, the matrix form
                                          !! of the sparsity pattern is printed
                                          !! (default is vector form)

    logical :: print_matrix !! if the matrix form is to be printed
    integer :: r   !! row counter
    character(len=1),dimension(n) :: row  !! a row of the sparsity matrix

    if (present(dense)) then
        print_matrix = dense
    else
        print_matrix = .false.  ! default
    end if

    write(iunit,'(A)') '---Sparsity pattern---'
    if (allocated(me%irow) .and. allocated(me%icol)) then

        if (print_matrix) then
            do r = 1,m    ! print by row
                row = '0'
                row(pack(me%icol,mask=me%irow==r)) = 'X'
                write(iunit,'(*(A1))') row
            end do
        else
            write(iunit,'(A,1X,*(I3,","))') 'irow:',me%irow
            write(iunit,'(A,1X,*(I3,","))') 'icol:',me%icol
        end if

        if (allocated(me%ngrp)) then
            write(iunit,'(A)') ''
            write(iunit,'(A)') '---Sparsity partition---'
            write(iunit,'(A,1x,I5)')       'Number of groups:',me%maxgrp
            write(iunit,'(A,1x,*(I5,1X))') 'Group array:     ',me%ngrp
        end if

    end if
    write(iunit,'(A)') ''

    ! print linear pattern if available
    if (me%linear_sparsity_computed) then
        write(iunit,'(A)') '---Linear sparsity pattern---'
        if (allocated(me%linear_icol) .and. allocated(me%linear_irow)) then
            if (print_matrix) then
                do r = 1,m    ! print by row
                    row = '0'
                    row(pack(me%linear_icol,mask=me%linear_irow==r)) = 'X'
                    write(iunit,'(*(A1))') row
                end do
            else
                write(iunit,'(A,1X,*(I3,","))') 'irow:',me%linear_irow
                write(iunit,'(A,1X,*(I3,","))') 'icol:',me%linear_icol
                write(iunit,'(A,1X,*(E30.16,","))') 'vals:',me%linear_vals
            end if
        end if
        write(iunit,'(A)') ''
    end if

    end subroutine print_sparsity
!*******************************************************************************

!*******************************************************************************
!>
!  Print the sparsity pattern in vector form (`irow`, `icol`).

    subroutine print_sparsity_pattern(me,iunit)

    implicit none

    class(numdiff_type),intent(inout) :: me
    integer,intent(in) :: iunit !! file unit to write to.
                                !! (assumed to be already opened)

    call me%sparsity%print(me%n,me%m,iunit,dense=.false.)

    end subroutine print_sparsity_pattern
!*******************************************************************************

!*******************************************************************************
!>
!  Print the sparsity pattern in matrix form.

    subroutine print_sparsity_matrix(me,iunit)

    implicit none

    class(numdiff_type),intent(inout) :: me
    integer,intent(in) :: iunit !! file unit to write to.
                                !! (assumed to be already opened)

    call me%sparsity%print(me%n,me%m,iunit,dense=.true.)

    end subroutine print_sparsity_matrix
!*******************************************************************************

!*******************************************************************************
!>
!  A user-callable routine. When called, it will terminate
!  all computations and return. The `istat` return code will be
!  set to `-1`. This can be called in the function or the info function.

    subroutine terminate(me)

    implicit none

    class(numdiff_type),intent(inout) :: me

    if (me%exception_raised) return ! check for existing exceptions

    me%istat = -1
    me%error_msg = 'Terminated by the user'
    me%exception_raised = .true.

    end subroutine terminate
!*******************************************************************************

!*******************************************************************************
!>
!  Raise an exception.

    subroutine raise_exception(me,istat,routine,error_msg)

    implicit none

    class(numdiff_type),intent(inout) :: me
    integer,intent(in)          :: istat      !! error code.
    character(len=*),intent(in) :: routine    !! the routine where the error was raised.
    character(len=*),intent(in) :: error_msg  !! error message string.

    me%istat = istat
    me%error_msg = 'Error in '//trim(routine)//' : '//trim(error_msg)
    me%exception_raised = .true.

    end subroutine raise_exception
!*******************************************************************************

!*******************************************************************************
!>
!  Clear all exceptions.

    subroutine clear_exceptions(me)

    implicit none

    class(numdiff_type),intent(inout) :: me

    me%istat = 0
    if (allocated(me%error_msg)) deallocate(me%error_msg)
    me%exception_raised = .false.

    end subroutine clear_exceptions
!*******************************************************************************

!*******************************************************************************
!>
!  Returns True if an exception has been raised.

    pure logical function failed(me)

    implicit none

    class(numdiff_type),intent(in) :: me

    failed = me%exception_raised

    end function failed
!*******************************************************************************

!*******************************************************************************
!>
!  Returns the current error code and message.

    subroutine get_error_status(me,istat,error_msg)

    implicit none

    class(numdiff_type),intent(in) :: me
    integer,intent(out),optional :: istat  !! error code (`istat=0` means no errors).
    character(len=:),allocatable,intent(out),optional :: error_msg !! error message string.

    if (present(istat)) istat = me%istat

    if (present(error_msg)) then
        if (allocated(me%error_msg)) then
            error_msg = me%error_msg
        else
            error_msg = ''
        end if
    end if

    end subroutine get_error_status
!*******************************************************************************

!*******************************************************************************
!>
!  Convert an integer to a string.

    function integer_to_string(i, with_sign) result(str)

    implicit none

    integer,intent(in)           :: i           !! the integer
    logical,intent(in),optional  :: with_sign   !! also include the sign (default is False)
    character(len=:),allocatable :: str         !! integer converted to a string

    integer :: istat !! `iostat` code for write statement
    character(len=100) :: tmp !!
    logical :: sgn !! local copy of `with_sign`

    if (present(with_sign)) then
        sgn = with_sign
    else
        sgn = .false.
    end if

    if (sgn) then
        write(tmp,'(SP,I100)',iostat=istat) i
    else
        write(tmp,'(I100)',iostat=istat) i
    end if

    if (istat == 0) then
        str = trim(adjustl(tmp))
    else
        str = '****'
    end if

    end function integer_to_string
!*******************************************************************************

!*******************************************************************************
    end module numerical_differentiation_module
!*******************************************************************************
