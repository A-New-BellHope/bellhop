MODULE norms

CONTAINS
  FUNCTION NORM2b( x )

    ! Compute the 2-norm of a vector
    ! This is an intrinsic function in the newest Fortran (NORM2)
    ! Provided here for users of an older Fortran

    IMPLICIT NONE
    REAL (KIND=8), INTENT( IN ) :: x( : )
    REAL (KIND=8)               :: NORM2b

    NORM2b = SQRT( DOT_PRODUCT( x, x ) )

  END FUNCTION NORM2b
END MODULE norms
