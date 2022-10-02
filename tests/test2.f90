!*******************************************************************************
!> author: Jacob Williams
!  date: October 1, 2022
!
!  Test2 for the numerical differentiation module.
!  Makes a plot of the errors using different methods and step sizes.

    program test2

    use iso_fortran_env
    use numerical_differentiation_module
    use numdiff_kinds_module, only: wp
    use pyplot_module

    implicit none

    integer,parameter :: n = 1 !! number of variables
    integer,parameter :: m = 1 !! number of functions
    real(wp),dimension(n),parameter :: x     = 1.0_wp !! point at which to compute the derivative
    real(wp),dimension(n),parameter :: xlow  = -100000.0_wp !! bounds not really needed
    real(wp),dimension(n),parameter :: xhigh =  100000.0_wp !!
    integer,parameter :: perturb_mode  = 1 !! absolute step
    integer,parameter :: cache_size    = 0 !! `0` indicates not to use cache
    integer,parameter :: sparsity_mode = 1 !! assume dense
    integer,dimension(*),parameter  :: methods = [1,3,10,21,36,500,600,700,800]
                                                !! array of method IDs:
                                                !! [1,3,10,21,36,500,600,700,800]
                                                !! forward + all the central diff ones
    integer,parameter :: exp_star = -ceiling(log10(epsilon(1.0_wp))) !! exponent to start with
    integer,parameter :: exp_stop = -2 !! exponent to end with
    integer,parameter :: exp_step = -1 !! ecponent step
    integer,parameter :: exp_scale = 5 !! number of substeps from one to the next

    type(numdiff_type)                  :: my_prob          !! main class to compute the derivatives
    integer                             :: i                !! method counter
    integer                             :: j                !! counter
    real(wp),dimension(:),allocatable   :: jac              !! jacobian
    integer                             :: func_evals       !! function evaluation counter
    character(len=:),allocatable        :: error_msg        !! error message string
    real(wp),dimension(n)               :: dpert            !! perturbation step size
    integer                             :: ipert            !! perturbation step size counter
    real(wp)                            :: error            !! diff from true derivative
    integer                             :: num_dperts       !! number of dperts to test
    integer                             :: num_methods      !! number of methods to test
    real(wp),dimension(:),allocatable   :: results_dpert    !! results array - dpert
    real(wp),dimension(:),allocatable   :: results_errors   !! results array - errors
    type(pyplot)                        :: plt              !! for plotting the results
    character(len=:),allocatable        :: formula          !! finite diff formula for the plot legend
    character(len=:),allocatable        :: name             !! finite diff name for the plot legend
    integer                             :: idx              !! index in results arrays
    character(len=:),allocatable        :: real_kind_str    !! real kind for the plot title
    real(wp),dimension(3)               :: color            !! line color array
    real(wp),dimension(2)               :: ylim             !! plot y limit array

    ! size the arrays:
    num_methods = size(methods)
    num_dperts = size([(i, i = exp_star*exp_scale, exp_stop*exp_scale, exp_step)])
    allocate(results_errors(num_dperts)); results_errors = -huge(1.0_wp)

    ! for plot title:
    select case(wp)
        case(REAL32);  real_kind_str = '[Single Precision]'
        case(REAL64);  real_kind_str = '[Double Precision]'
        case(REAL128); real_kind_str = '[Quad Precision]'
        case default; error stop 'Invalid real kind'
    end select

    ! initialize the plot:
    call plt%initialize(grid            = .true.,&
                        figsize         = [20,10],&
                        axes_labelsize  = 30, &
                        xtick_labelsize = 30, &
                        ytick_labelsize = 30, &
                        font_size       = 30, &
                        xlabel          = 'Finite Difference Perturbation Step Size $h$',&
                        ylabel          = 'Finite Difference Derivative Error',&
                        title           = 'Derivative of $x + \sin(x)$ at $x=1$ '//real_kind_str,&
                        legend          = .true., &
                        legend_fontsize = 10,&
                        usetex          = .true.)

    ! try different finite diff methods
    do j=1,size(methods)

        idx = 0
        i = methods(j)  ! method id
        call get_finite_diff_formula(i,formula,name)

        ! cycle through perturbation step sizes:
        do ipert = exp_star*exp_scale, exp_stop*exp_scale, exp_step

            idx = idx + 1 ! index for arrays
            dpert = 10.0_wp ** (-ipert/real(exp_scale,wp)) ! compute perturbation step size
            if (j==1) then
                if (.not. allocated(results_dpert)) allocate(results_dpert(0))
                results_dpert = [results_dpert, dpert(1)] ! save dpert
            end if

            func_evals = 0
            call my_prob%destroy()
            call my_prob%initialize(n,m,xlow,xhigh,perturb_mode,dpert,&
                                    problem_func=my_func,&
                                    sparsity_mode=sparsity_mode,&
                                    jacobian_method=i,&
                                    partition_sparsity_pattern=.false.,&
                                    cache_size=cache_size)
            if (my_prob%failed()) then
                call my_prob%get_error_status(error_msg=error_msg)
                error stop error_msg
            end if

            call my_prob%compute_jacobian(x,jac)
            if (my_prob%failed()) then
                call my_prob%get_error_status(error_msg=error_msg)
                error stop error_msg
            end if
            error = abs(deriv(x(1)) - jac(1))
            !write(output_unit,'(I5,1X,*(F27.16))') i , dpert, error
            if (error<=epsilon(1.0_wp)) error = 0.0_wp
            results_errors(idx) = error ! save result

        end do

        ! line color for the plot:
        select case (j)
        case(1); color = [1.0_wp, 0.0_wp, 0.0_wp] ! red
        case(2); color = [0.0_wp, 1.0_wp, 0.0_wp] ! green
        case default
            ! blue gradient for the others
            color = [real(j-2,wp)/num_methods, &
                     real(j-2,wp)/num_methods, &
                     0.9_wp]
        end select
        ylim = [10.0_wp**(ceiling(log10(epsilon(1.0_wp)))), 1.0_wp]

        ! plot for this method:
        call plt%add_plot(results_dpert,results_errors,&
                          xscale='log', yscale='log',&
                          label = formula,linestyle='.-',markersize=5,linewidth=2, &
                          color = color,&
                          xlim = ylim,&
                          ylim = ylim)

    end do

    ! save plot:
    call plt%savefig('results '//real_kind_str//'.pdf', pyfile='results.py')

contains

    function func(x)
        !! Problem function
        real(wp), intent(in) :: x
        real(wp) :: func
        func = x + sin(x)
    end function func
    function deriv(x)
        !! Problem function true derivative
        real(wp), intent(in) :: x
        real(wp) :: deriv
        deriv = 1.0_wp + cos(x)
    end function deriv

    subroutine my_func(me,x,f,funcs_to_compute)
        !! Problem function interface for numdiff
        class(numdiff_type),intent(inout) :: me
        real(wp),dimension(:),intent(in)  :: x
        real(wp),dimension(:),intent(out) :: f
        integer,dimension(:),intent(in)   :: funcs_to_compute
        if (any(funcs_to_compute==1)) f(1) = func(x(1))
        func_evals = func_evals + 1
    end subroutine my_func

!*******************************************************************************
    end program test2
!*******************************************************************************
