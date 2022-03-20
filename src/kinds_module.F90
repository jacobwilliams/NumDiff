!*******************************************************************************
!>
!  Numeric kinds for NumDiff.
!
!@note The default real kind (`wp`) can be
!      changed using optional preprocessor flags.
!      This library was built with real kind:
#ifdef REAL32
!      `real(kind=real32)` [4 bytes]
#elif REAL64
!      `real(kind=real64)` [8 bytes]
#elif REAL128
!      `real(kind=real128)` [16 bytes]
#else
!      `real(kind=real64)` [8 bytes]
#endif

    module numdiff_kinds_module

    use iso_fortran_env

    private

#ifdef REAL32
    integer,parameter,public :: wp = real32   !! default real kind [4 bytes]
#elif REAL64
    integer,parameter,public :: wp = real64   !! default real kind [8 bytes]
#elif REAL128
    integer,parameter,public :: wp = real128  !! default real kind [16 bytes]
#else
    integer,parameter,public :: wp = real64   !! default real kind [8 bytes]
#endif

    end module numdiff_kinds_module
!*******************************************************************************