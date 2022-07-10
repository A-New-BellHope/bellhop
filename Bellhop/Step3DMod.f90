MODULE Step3DMod

  USE bellhopMod
  USE Bdry3DMod
  USE sspMod
  IMPLICIT NONE
  
  REAL (KIND=8), PARAMETER, PRIVATE :: INFINITESIMAL_STEP_SIZE = 1.0d-6
  
CONTAINS

  SUBROUTINE Step3D( ray0, ray2, topRefl, botRefl )

    ! Does a single step along the ray
    ! x denotes the ray coordinate, ( x, y, z )
    ! t denotes the scaled tangent to the ray (previously (xi, eta, zeta) )
    ! c * t would be the unit tangent

    USE RayNormals

    ! rays
    TYPE( ray3DPt ) :: ray0, ray1, ray2
    LOGICAL, INTENT( OUT ) :: topRefl, botRefl
    INTEGER         :: iSegx0, iSegy0, iSegz0
    REAL  (KIND=8 ) :: gradc0( 3 ), gradc1( 3 ), gradc2( 3 ), &
         c0, cimag0, csq0, cxx0, cyy0, czz0, cxy0, cxz0, cyz0, cnn0, cmn0, cmm0, &
         c1, cimag1, csq1, cxx1, cyy1, czz1, cxy1, cxz1, cyz1, cnn1, cmn1, cmm1, &
         c2, cimag2,       cxx2, cyy2, czz2, cxy2, cxz2, cyz2, c_mat1( 2, 2 ), &
         urayt0( 3 ), urayt1( 3 ), urayt2( 3 ), h, halfh, hw0, hw1, w0, w1, rho
    COMPLEX( KIND=8 ) :: d_phi0, d_phi1
    REAL  (KIND=8 ) :: d_p_tilde0( 2 ), d_p_hat0( 2 ), d_q_tilde0( 2 ), d_q_hat0( 2 ), &
                       d_p_tilde1( 2 ), d_p_hat1( 2 ), d_q_tilde1( 2 ), d_q_hat1( 2 )
    REAL  (KIND=8 ) :: gradcjump( 3 ), csjump, cn1jump, cn2jump, tBdry( 3 ), nBdry( 3 ), RM, R1, R2, Tg, Th, &
         rayt( 3 ), rayn1( 3 ), rayn2( 3 )
    REAL  (KIND=8 ) :: e1( 3 ), e2( 3 )                              ! ray normals for ray-centered coordinates
    REAL  (KIND=8 ) :: p_tilde_in(  2 ), p_hat_in(  2 ), q_tilde_in(  2 ), q_hat_in(  2 ), &
         p_tilde_out( 2 ), p_hat_out( 2 ), RotMat( 2, 2 )
         
    WRITE( PRTFile, * )
    WRITE( PRTFile, * ) 'ray0 x t', ray0%x, ray0%t
    WRITE( PRTFile, * ) 'iSegx iSegy iSegz', iSegx, iSegy, iSegz
           
    ! The numerical integrator used here is a version of the polygon (a.k.a. midpoint, leapfrog, or Box method), and similar
    ! to the Heun (second order Runge-Kutta method).
    ! However, it's modified to allow for a dynamic step change, while preserving the second-order accuracy).

    ! *** Phase 1 (an Euler step)

    CALL EvaluateSSP3D( ray0%x, ray0%t, c0, cimag0, gradc0, cxx0, cyy0, czz0, cxy0, cxz0, cyz0, rho, freq, 'TAB' )
    CALL RayNormal( ray0%t, ray0%phi, c0, e1, e2 ) ! Compute ray normals e1 and e2
    CALL Get_c_partials( cxx0, cxy0, cxz0, cyy0, cyz0, czz0, e1, e2, cnn0, cmn0, cmm0 ) ! Compute second partials of c along ray normals
    CALL ComputeDeltaPQ( ray0, c0, gradc0, cnn0, cmn0, cmm0, d_phi0, d_p_tilde0, d_p_hat0, d_q_tilde0, d_q_hat0)


    iSegx0 = iSegx     ! make note of current layer
    iSegy0 = iSegy
    iSegz0 = iSegz

    csq0   = c0 * c0
    urayt0 = c0 * ray0%t  ! unit tangent
    h = Beam%deltas       ! initially set the step h, to the basic one, deltas

    CALL ReduceStep3D( ray0%x, urayt0, iSegx0, iSegy0, iSegz0, h ) ! reduce h to land on boundary
    halfh = 0.5 * h       ! first step of the modified polygon method is a half step

    ray1%x = ray0%x    + halfh * urayt0
    ray1%t = ray0%t    - halfh * gradc0 / csq0
    CALL UpdateRayPQ ( ray1, ray0, halfh, d_phi0, d_p_tilde0, d_p_hat0, d_q_tilde0, d_q_hat0 )

    ! *** Phase 2

    CALL EvaluateSSP3D( ray1%x, ray1%t, c1, cimag1, gradc1, cxx1, cyy1, czz1, cxy1, cxz1, cyz1, rho, freq, 'TAB' )
    ! LP: Fixed; should be ray1%phi; ray2%phi would be uninitialized memory or
    ! left over from the previous ray
    CALL RayNormal( ray1%t, ray1%phi, c1, e1, e2 )
    CALL Get_c_partials( cxx1, cxy1, cxz1, cyy1, cyz1, czz1, e1, e2, cnn1, cmn1, cmm1 ) ! Compute second partials of c along ray normals
    CALL ComputeDeltaPQ( ray1, c1, gradc1, cnn1, cmn1, cmm1, d_phi1, d_p_tilde1, d_p_hat1, d_q_tilde1, d_q_hat1)

    csq1   = c1 * c1
    urayt1 = c1 * ray1%t   ! unit tangent

    CALL ReduceStep3D( ray0%x, urayt1, iSegx0, iSegy0, iSegz0, h ) ! reduce h to land on boundary

    ! use blend of f' based on proportion of a full step used.
    w1  = h / ( 2.0d0 * halfh )
    w0  = 1.0d0 - w1
    urayt2 = w0 * urayt0 + w1 * urayt1
    ! Take the blended ray tangent ( urayt2 ) and find the minimum step size ( h )
    ! to put this on a boundary, and ensure that the resulting position
    ! ( ray2%x ) gets put precisely on the boundary.
    CALL StepToBdry3D( ray0%x, ray2%x, urayt2, iSegx0, iSegy0, iSegz0, h, topRefl, botRefl )
    
    ! Update other variables with this new h
    ! LP: Fixed: ray2%phi now depends on hw0 & hw1 like the other parameters,
    ! originally only depended on h and ray1 vars
    hw0 = h * w0
    hw1 = h * w1
    ray2%t   = ray0%t   - hw0 * gradc0 / csq0 - hw1 * gradc1 / csq1
    ray2%tau = ray0%tau + hw0 / CMPLX( c0, cimag0, KIND=8 ) + hw1 / CMPLX( c1, cimag1, KIND=8 )
    CALL UpdateRayPQ ( ray2, ray0, hw0, d_phi0, d_p_tilde0, d_p_hat0, d_q_tilde0, d_q_hat0 )
    CALL UpdateRayPQ ( ray2, ray2, hw1, d_phi1, d_p_tilde1, d_p_hat1, d_q_tilde1, d_q_hat1 ) ! Not a typo, accumulating into 2

    ray2%Amp       = ray0%Amp
    ray2%Phase     = ray0%Phase
    ray2%NumTopBnc = ray0%NumTopBnc
    ray2%NumBotBnc = ray0%NumBotBnc

    ! *** If we crossed an interface, apply jump condition ***

    CALL EvaluateSSP3D( ray2%x, ray2%t, c2, cimag2, gradc2, cxx2, cyy2, czz2, cxy2, cxz2, cyz2, rho, freq, 'TAB' )

    ray2%c = c2

    IF ( iSegx /= iSegx0 .OR. &
         iSegy /= iSegy0 .OR. &
         iSegz /= iSegz0 ) THEN

       gradcjump =  gradc2 - gradc0  ! this is precise only for c-linear layers

!!! what if we cross isegx, isegy, or isegz at the same time?
       IF ( iSegz /= iSegz0 ) THEN
          nBdry = [ 0D0, 0D0, -SIGN( 1.0D0, ray2%t( 3 ) ) ]   ! inward normal to layer
       ELSE IF ( iSegx /= iSegx0 ) THEN
          nBdry = [ -SIGN( 1.0D0, ray2%t( 1 ) ), 0D0, 0D0 ]   ! inward normal to x-segment
       ELSE
          nBdry = [ 0D0, -SIGN( 1.0D0, ray2%t( 2 ) ), 0D0 ]   ! inward normal to y-segment
       END IF
       CALL CurvatureCorrection

    END IF

  CONTAINS
    SUBROUTINE CurvatureCorrection

      ! correct p-q due to jumps in the gradient of the sound speed

      USE cross_products

      Th    = DOT_PRODUCT( ray2%t, nBdry )   ! component of ray tangent, normal to boundary
      tBdry = ray2%t - Th * nBdry            ! tangent, along the boundary, in the reflection plane
      tBdry = tBdry / NORM2( tBdry )         ! unit boundary tangent
      Tg    = DOT_PRODUCT( ray2%t, tBdry )   ! component of ray tangent, along the boundary

      rayt = c2 * ray2%t                     ! unit tangent to ray

      rayn2 = cross_product( rayt, nBdry )   ! ray tangent x boundary normal gives refl. plane normal
      rayn2 = rayn2 / NORM2( rayn2 )         ! unit normal
      rayn1 = -cross_product( rayt, rayn2 )  ! ray tangent x refl. plane normal is first ray normal

      ! normal and tangential derivatives of the sound speed
      cn1jump = DOT_PRODUCT( gradcjump, rayn1 )
      cn2jump = DOT_PRODUCT( gradcjump, rayn2 )
      csjump  = DOT_PRODUCT( gradcjump, rayt  )

      RM = Tg / Th   ! this is tan( alpha ) where alpha is the angle of incidence
      R1 = RM * ( 2 * cn1jump - RM * csjump ) / c2 ** 2
      R2 = RM * cn2jump / c2 ** 2

      ! *** curvature correction ***

      CALL RayNormal_unit( rayt, ray2%phi, e1, e2 )   ! Compute ray normals e1 and e2

      RotMat( 1, 1 ) = DOT_PRODUCT( rayn1, e1 )
      RotMat( 1, 2 ) = DOT_PRODUCT( rayn1, e2 )
      RotMat( 2, 1 ) = -RotMat( 1, 2 )             ! DOT_PRODUCT( rayn2, e1 )
      RotMat( 2, 2 ) = DOT_PRODUCT( rayn2, e2 )

      ! rotate p-q values in e1, e2 system, onto rayn1, rayn2 system

      p_tilde_in = RotMat( 1, 1 ) * ray2%p_tilde + RotMat( 1, 2 ) * ray2%p_hat
      p_hat_in   = RotMat( 2, 1 ) * ray2%p_tilde + RotMat( 2, 2 ) * ray2%p_hat

      q_tilde_in = RotMat( 1, 1 ) * ray2%q_tilde + RotMat( 1, 2 ) * ray2%q_hat
      q_hat_in   = RotMat( 2, 1 ) * ray2%q_tilde + RotMat( 2, 2 ) * ray2%q_hat

      ! here's the actual curvature change

      p_tilde_out = p_tilde_in - q_tilde_in * R1 - q_hat_in * R2  
      p_hat_out   = p_hat_in   - q_tilde_in * R2

      ! rotate p back to e1, e2 system, q does not change
      ! Note RotMat^(-1) = RotMat^T

      ray2%p_tilde = RotMat( 1, 1 ) * p_tilde_out + RotMat( 2, 1 ) * p_hat_out
      ray2%p_hat   = RotMat( 1, 2 ) * p_tilde_out + RotMat( 2, 2 ) * p_hat_out

    END SUBROUTINE CurvatureCorrection

    !**********************************************************************!

    SUBROUTINE Get_c_partials( cxx, cxy, cxz, cyy, cyz, czz, e1, e2, cnn, cmn, cmm )

      ! Computes the second partials of c along ray normals

      REAL (KIND=8), INTENT( IN  ) :: cxx, cxy, cxz, cyy, cyz, czz  ! curvature of sound speed (cartesian)
      REAL (KIND=8), INTENT( IN  ) :: e1( 3 ), e2( 3 )              ! principal normals
      REAL (KIND=8), INTENT( OUT ) :: cnn, cmn, cmm                 ! curvature of sound speed (ray-centered)

      cnn = cxx * e1( 1 )**2 + cyy * e1( 2 )**2 + czz * e1( 3 )**2 + 2.0 * cxy * e1( 1 ) * e1( 2 ) + &
           2.0 * cxz * e1( 1 ) * e1( 3 ) + 2.0 * cyz * e1( 2 ) * e1( 3 )

      cmn = cxx * e1( 1 ) * e2( 1 ) + cyy * e1( 2 ) * e2( 2 ) + czz * e1( 3 ) * e2( 3 ) + &
           cxy * ( e1( 1 ) * e2( 2 ) + e2( 1 ) * e1( 2 ) ) + cxz * ( e1( 1 ) * e2( 3 ) + e2( 1 ) * e1( 3 ) ) +  &
           cyz * ( e1( 2 ) * e2( 3 ) + e2( 2 ) * e1( 3 ) )

      cmm = cxx * e2( 1 )**2 + cyy * e2( 2 )**2 + czz * e2( 3 )**2 + 2.0 * cxy * e2( 1 ) * e2( 2 ) + &
           2.0 * cxz * e2( 1 ) * e2( 3 ) + 2.0 * cyz * e2( 2 ) * e2( 3 )

      RETURN
    END SUBROUTINE Get_c_partials
  END SUBROUTINE Step3D
  
  SUBROUTINE ComputeDeltaPQ( ray, c, gradc, cnn, cmn, cmm, d_phi, d_p_tilde, d_p_hat, d_q_tilde, d_q_hat )
    TYPE(    ray3DPt ), INTENT( IN )  :: ray
    REAL(    KIND=8  ), INTENT( IN )  :: c, gradc( 3 ), cnn, cmn, cmm
    COMPLEX( KIND=8  ), INTENT( OUT ) :: d_phi
    REAL(    KIND=8  ), INTENT( OUT ) :: d_p_tilde( 2 ), d_p_hat( 2 ), d_q_tilde( 2 ), d_q_hat( 2 )
    REAL(    KIND=8  )                :: c_mat( 2, 2 ), csq
    
    d_phi = ( 1.0 / c ) * ray%t( 3 ) * &
         ( ray%t( 2 ) * gradc( 1 ) - ray%t( 1 ) * gradc( 2 ) ) / &
         ( ray%t( 1 ) ** 2 + ray%t( 2 ) ** 2 )
    
    csq = c ** 2
    c_mat( 1, : ) = -[ cnn, cmn ] / csq
    c_mat( 2, : ) = -[ cmn, cmm ] / csq
    
    d_p_tilde = MATMUL( c_mat, ray%q_tilde )
    d_p_hat   = MATMUL( c_mat, ray%q_hat )

    d_q_tilde =         c    * ray%p_tilde 
    d_q_hat   =         c    * ray%p_hat 
    
    ! d_f    = ( c0 * ray0%DetP - cnn0 / csq0 * ray0%DetQ )
    ! d_g    = ( c0 * ray0%DetP - cmm0 / csq0 * ray0%DetQ )
    ! d_h    = (                - cmn0 / csq0 * ray0%DetQ )
    ! d_DetP = 1.0D0 / csq0 * ( -cmm0 * ray0%f - cnn0 * ray0%g + 2.0 * cmn0 * ray0%h )
    ! d_DetQ = c0 * ( ray0%f + ray0%g )
    
  END SUBROUTINE ComputeDeltaPQ
  
  SUBROUTINE UpdateRayPQ ( ray1, ray0, h, d_phi, d_p_tilde, d_p_hat, d_q_tilde, d_q_hat )
    TYPE(    ray3DPt ), INTENT( INOUT ) :: ray1
    TYPE(    ray3DPt ), INTENT( IN ) :: ray0
    REAL(    KIND=8  ), INTENT( IN ) :: h
    COMPLEX( KIND=8  ), INTENT( IN ) :: d_phi
    REAL(    KIND=8  ), INTENT( IN ) :: d_p_tilde( 2 ), d_p_hat( 2 ), d_q_tilde( 2 ), d_q_hat( 2 )
    
    ray1%phi = ray0%phi + h * d_phi
    ray1%p_tilde = ray0%p_tilde + h * d_p_tilde
    ray1%q_tilde = ray0%q_tilde + h * d_q_tilde 
    ray1%p_hat   = ray0%p_hat   + h * d_p_hat
    ray1%q_hat   = ray0%q_hat   + h * d_q_hat 
    
    ! LP: no longer missing the hw0 / hw1 blend
    ! ray1%f    = ray0%f    + h * d_f
    ! ray1%g    = ray0%g    + h * d_g
    ! ray1%h    = ray0%h    + h * d_h
    ! ray1%DetP = ray0%DetP + h * d_DetP
    ! ray1%DetQ = ray0%DetQ + h * d_DetQ

  END SUBROUTINE UpdateRayPQ

  ! **********************************************************************!

  SUBROUTINE ReduceStep3D( x0, urayt, iSegx0, iSegy0, iSegz0, h )

    INTEGER,       INTENT( IN    ) :: iSegx0, iSegy0, iSegz0
    REAL (KIND=8), INTENT( IN    ) :: x0( 3 ), urayt( 3 )  ! ray coordinate and tangent
    REAL (KIND=8), INTENT( INOUT ) :: h                    ! reduced step size
    REAL (KIND=8) :: d( 3 ), d0( 3 ), tri_n( 3 )
    REAL (KIND=8) :: x( 3 )                              ! ray coordinate after full trial step
    REAL (KIND=8) :: h1, h2, h3, h4, h5, h6, h7          ! step sizes
    REAL (KIND=8) :: xSeg( 2 ), ySeg( 2 )

    ! Detect interface or boundary crossing and reduce step, if necessary, to land on that crossing.
    ! Keep in mind possibility that user put source right on an interface
    ! and that multiple events can occur (crossing interface, top, and bottom in a single step).

    x = x0 + h * urayt ! make a trial step

    ! interface crossing in depth
    ! Step reduction is not done for the top or bottom layer
    ! Instead the SSP is extrapolated
    ! This prevents problems when the boundaries are outside the domain of the SSP
    h1 = huge( h1 )
    IF ( ABS( urayt( 3 ) ) > EPSILON( h1 ) ) THEN
       IF      ( SSP%z( iSegz0     ) > x(  3 ) .AND. iSegz0     > 1  ) THEN
          h1 = ( SSP%z( iSegz0     ) - x0( 3 ) ) / urayt( 3 )
          WRITE( PRTFile, * ) 'Shallower bound SSP Z > z; h', SSP%z( iSegz0     ), x( 3 ), h1
       ELSE IF ( SSP%z( iSegz0 + 1 ) < x(  3 ) .AND. iSegz0 + 1 < SSP%Nz ) THEN
          h1 = ( SSP%z( iSegz0 + 1 ) - x0( 3 ) ) / urayt( 3 )
          WRITE( PRTFile, * ) 'Deeper bound SSP Z < z; h', SSP%z( iSegz0 + 1 ), x( 3 ), h1
       END IF
    END IF

    ! top crossing
    h2 = huge( h2 )
    d  = x - Topx              ! vector from top to ray
    ! LP: Changed from EPSILON( h1 ) (should have been h2!) to 0 in conjunction 
    ! with StepToBdry3D change here.
    IF ( DOT_PRODUCT( Topn, d )  >= 0.0D0 ) THEN
       d0 = x0 - Topx   ! vector from top    node to ray origin
       h2 = -DOT_PRODUCT( d0, Topn ) / DOT_PRODUCT( urayt, Topn )
       WRITE( PRTFile, * ) 'Top crossing h', h2
    END IF

    ! bottom crossing
    h3 = huge( h3 )
    d  = x - Botx              ! vector from bottom to ray
    ! LP: Changed from EPSILON( h1 ) (should have been h3!) to 0 in conjunction 
    ! with StepToBdry3D change here.
    IF ( DOT_PRODUCT( Botn, d ) >= 0.0D0 ) THEN
       d0 = x0 - Botx   ! vector from bottom node to ray origin
       h3 = -DOT_PRODUCT( d0, Botn ) / DOT_PRODUCT( urayt, Botn )
       WRITE( PRTFile, * ) 'Bottom crossing h', h3
    END IF

    ! top/bottom segment crossing in x
    h4 = huge( h4 )
    xSeg( 1 ) = MAX( xTopSeg( 1 ), xBotSeg( 1 ) )
    xSeg( 2 ) = MIN( xTopSeg( 2 ), xBotSeg( 2 ) )
    
    IF ( SSP%Type == 'H' ) THEN   ! ocean segment
       xSeg( 1 ) = MAX( xSeg( 1 ), SSP%Seg%x( iSegx0     ) )
       xSeg( 2 ) = MIN( xSeg( 2 ), SSP%Seg%x( iSegx0 + 1 ) )
    END IF

    IF ( ABS( urayt( 1 ) ) > EPSILON( h1 ) ) THEN
       IF       ( x(  1 ) < xSeg( 1 ) ) THEN
          h4 = -( x0( 1 ) - xSeg( 1 ) ) / urayt( 1 )
          WRITE( PRTFile, * ) 'Min bound SSP X > x; h', xSeg( 1 ), x( 1 ), h4
       ELSE IF  ( x(  1 ) > xSeg( 2 ) ) THEN
          h4 = -( x0( 1 ) - xSeg( 2 ) ) / urayt( 1 )
          WRITE( PRTFile, * ) 'Max bound SSP X < x; h', xSeg( 2 ), x( 1 ), h4
       END IF
    END IF

    ! top/bottom segment crossing in y
    h5 = huge( h5 )
    ySeg( 1 ) = MAX( yTopSeg( 1 ), yBotSeg( 1 ) )
    ySeg( 2 ) = MIN( yTopSeg( 2 ), yBotSeg( 2 ) )

    IF ( SSP%Type == 'H' ) THEN   ! ocean segment
       ySeg( 1 ) = MAX( ySeg( 1 ), SSP%Seg%y( iSegy0     ) )
       ySeg( 2 ) = MIN( ySeg( 2 ), SSP%Seg%y( iSegy0 + 1 ) )
    END IF

    ! LP: removed 1000 * epsilon which mbp had comment questioning also
    IF ( ABS( urayt( 2 ) ) > EPSILON( h1 ) ) THEN
       IF       ( x(  2 ) < ySeg( 1 ) ) THEN
          h5 = -( x0( 2 ) - ySeg( 1 ) ) / urayt( 2 )
          WRITE( PRTFile, * ) 'Min bound SSP Y > y; h', ySeg( 1 ), x( 2 ), h5
       ELSE IF  ( x(  2 ) > ySeg( 2 ) ) THEN
          h5 = -( x0( 2 ) - ySeg( 2 ) ) / urayt( 2 )
          WRITE( PRTFile, * ) 'Max bound SSP Y < y; h', ySeg( 2 ), x( 2 ), h5
       END IF
    END IF

    ! triangle crossing within a top segment
    h6    = huge( h6 )
    d     = x  - Topx   ! vector from bottom node to ray end
    d0    = x0 - Topx   ! vector from bottom node to ray origin
    tri_n = [ -( yTopSeg( 2 ) - yTopSeg( 1 ) ), xTopSeg( 2 ) - xTopSeg( 1 ), 0.0d0 ]

    IF ( CheckDiagCrossing( tri_n, d0, d ) ) THEN
       h6 = -DOT_PRODUCT( d0, tri_n ) / DOT_PRODUCT( urayt, tri_n )
       WRITE( PRTFile, * ) 'Tri diag crossing h6', h6
    END IF

    ! triangle crossing within a bottom segment
    h7    = huge( h7 )
    d     = x  - Botx   ! vector from bottom node to ray end
    d0    = x0 - Botx   ! vector from bottom node to ray origin
    tri_n = [ -( yBotSeg( 2 ) - yBotSeg( 1 ) ), xBotSeg( 2 ) - xBotSeg( 1 ), 0.0d0 ]

    IF ( CheckDiagCrossing( tri_n, d0, d ) ) THEN
       h7 = -DOT_PRODUCT( d0, tri_n ) / DOT_PRODUCT( urayt, tri_n )
       WRITE( PRTFile, * ) 'Tri diag crossing h7', h7
    END IF

    h = MIN( h, h1, h2, h3, h4, h5, h6, h7 )  ! take limit set by shortest distance to a crossing

    IF ( h < INFINITESIMAL_STEP_SIZE * Beam%deltas ) THEN        ! is it taking an infinitesimal step?
       h = INFINITESIMAL_STEP_SIZE * Beam%deltas                 ! make sure we make some motion
       iSmallStepCtr = iSmallStepCtr + 1   ! keep a count of the number of sequential small steps
    ELSE
       iSmallStepCtr = 0                   ! didn't do a small step so reset the counter
    END IF
  END SUBROUTINE ReduceStep3D
  
  ! **********************************************************************!

  SUBROUTINE StepToBdry3D( x0, x2, urayt, iSegx0, iSegy0, iSegz0, h, topRefl, botRefl )

    INTEGER,       INTENT( IN    ) :: iSegx0, iSegy0, iSegz0
    REAL (KIND=8), INTENT( IN    ) :: x0( 3 ), urayt( 3 )  ! ray coordinate and tangent
    REAL (KIND=8), INTENT( INOUT ) :: x2( 3 ), h           ! output coord, reduced step size
    LOGICAL,       INTENT( OUT   ) :: topRefl, botRefl
    REAL (KIND=8) :: d( 3 ), d0( 3 ), tri_n( 3 )
    REAL (KIND=8) :: xSeg( 2 ), ySeg( 2 )                  ! boundary limits
    INTEGER                        :: k

    ! Original step due to maximum step size
    h = Beam%deltas
    x2 = x0 + h * urayt

    ! interface crossing in depth
    IF ( ABS( urayt( 3 ) ) > EPSILON( h ) ) THEN
       IF      ( SSP%z( iSegz0     ) > x2( 3 ) .AND. iSegz0     > 1  ) THEN
          h  = ( SSP%z( iSegz0     ) - x0( 3 ) ) / urayt( 3 )
          x2 = x0 + h * urayt
          x2( 3 ) = SSP%z( iSegz0 )
          WRITE( PRTFile, * ) 'StepToBdry3D shallower h to', h, x2
       ELSE IF ( SSP%z( iSegz0 + 1 ) < x2( 3 ) .AND. iSegz0 + 1 < SSP%Nz ) THEN
          h  = ( SSP%z( iSegz0 + 1 ) - x0( 3 ) ) / urayt( 3 )
          x2 = x0 + h * urayt
          x2( 3 ) = SSP%z( iSegz0 + 1 )
          WRITE( PRTFile, * ) 'StepToBdry3D deeper h to', h, x2
       END IF
    END IF

    ! top crossing
    d  = x2 - Topx ! vector from top to ray
    d0 = x0 - Topx ! vector from top node to ray origin
    ! LP: Originally, this value had to be > a small positive number, meaning the
    ! new step really had to be outside the boundary, not just to the boundary.
    ! Also, this is not missing a normalization factor, Topn is normalized so
    ! this is actually the distance above the top in meters.
    IF ( DOT_PRODUCT( Topn, d )  > -INFINITESIMAL_STEP_SIZE ) THEN
       d0 = x0 - Topx   ! vector from top    node to ray origin
       h  = -DOT_PRODUCT( d0, Topn ) / DOT_PRODUCT( urayt, Topn )
       x2 = x0 + h * urayt
       ! Snap to exact top depth value if it's flat
       IF ( ABS( Topn( 1 ) ) < EPSILON( Topn( 1 ) ) .AND. ABS( Topn( 2 ) ) < EPSILON( Topn( 2 ) ) ) THEN
          x2( 3 ) = Topx( 3 )
       END IF
       WRITE( PRTFile, * ) 'StepToBdry3D top crossing h to', h, x2
       topRefl = .TRUE.
    ELSE
       topRefl = .FALSE.
    END IF

    ! bottom crossing
    d  = x2 - Botx ! vector from bottom to ray
    d0 = x0 - Botx ! vector from bottom node to ray origin
    ! LP: See comment above for top case.
    IF ( DOT_PRODUCT( Botn, d ) > -INFINITESIMAL_STEP_SIZE ) THEN
       d0 = x0 - Botx   ! vector from bottom node to ray origin
       h  = -DOT_PRODUCT( d0, Botn ) / DOT_PRODUCT( urayt, Botn )
       x2 = x0 + h * urayt
       ! Snap to exact bottom depth value if it's flat
       IF ( ABS( Botn( 1 ) ) < EPSILON( Botn( 1 ) ) .AND. ABS( Botn( 2 ) ) < EPSILON( Botn( 2 ) ) ) THEN
          x2( 3 ) = Botx( 3 )
       END IF
       WRITE( PRTFile, * ) 'StepToBdry3D bottom crossing h to', h, x2
       botRefl = .TRUE.
       ! Should not ever be able to cross both, but in case it does, make sure
       ! only the crossing we exactly landed on is active
       topRefl = .FALSE.
    ELSE
       botRefl = .FALSE.
    END IF

    ! top/bottom segment crossing in x
    xSeg( 1 ) = MAX( xTopSeg( 1 ), xBotSeg( 1 ) ) ! LP: lower range bound (not an x value)
    xSeg( 2 ) = MIN( xTopSeg( 2 ), xBotSeg( 2 ) ) ! LP: upper range bound (not a y value)
    
    IF ( SSP%Type == 'H' ) THEN   ! ocean segment
       xSeg( 1 ) = MAX( xSeg( 1 ), SSP%Seg%x( iSegx0     ) )
       xSeg( 2 ) = MIN( xSeg( 2 ), SSP%Seg%x( iSegx0 + 1 ) )
    END IF

    IF ( ABS( urayt( 1 ) ) > EPSILON( h ) ) THEN
       IF       ( x2( 1 ) < xSeg( 1 ) ) THEN
          h  = -( x0( 1 ) - xSeg( 1 ) ) / urayt( 1 )
          x2 = x0 + h * urayt
          x2( 1 ) = xSeg( 1 )
          topRefl = .FALSE.
          botRefl = .FALSE.
          WRITE( PRTFile, * ) 'StepToBdry3D X min bound h to', h, x2
       ELSE IF  ( x2( 1 ) > xSeg( 2 ) ) THEN
          h  = -( x0( 1 ) - xSeg( 2 ) ) / urayt( 1 )
          x2 = x0 + h * urayt
          x2( 1 ) = xSeg( 2 )
          topRefl = .FALSE.
          botRefl = .FALSE.
          WRITE( PRTFile, * ) 'StepToBdry3D X max bound h to', h, x2
       END IF
    END IF

    ! top/bottom segment crossing in y
    ySeg( 1 ) = MAX( yTopSeg( 1 ), yBotSeg( 1 ) )
    ySeg( 2 ) = MIN( yTopSeg( 2 ), yBotSeg( 2 ) )

    IF ( SSP%Type == 'H' ) THEN   ! ocean segment
       ySeg( 1 ) = MAX( ySeg( 1 ), SSP%Seg%y( iSegy0     ) )
       ySeg( 2 ) = MIN( ySeg( 2 ), SSP%Seg%y( iSegy0 + 1 ) )
    END IF

    ! LP: removed 1000 * epsilon which mbp had comment questioning also
    IF ( ABS( urayt( 2 ) ) > EPSILON( h ) ) THEN
       IF       ( x2( 2 ) < ySeg( 1 ) ) THEN
          h  = -( x0( 2 ) - ySeg( 1 ) ) / urayt( 2 )
          x2 = x0 + h * urayt
          x2( 1 ) = ySeg( 1 )
          topRefl = .FALSE.
          botRefl = .FALSE.
          WRITE( PRTFile, * ) 'StepToBdry3D Y min bound h to', h, x2
       ELSE IF  ( x2( 2 ) > ySeg( 2 ) ) THEN
          h  = -( x0( 2 ) - ySeg( 2 ) ) / urayt( 2 )
          x2 = x0 + h * urayt
          x2( 1 ) = ySeg( 2 )
          topRefl = .FALSE.
          botRefl = .FALSE.
          WRITE( PRTFile, * ) 'StepToBdry3D Y max bound h to', h, x2
       END IF
    END IF

    ! triangle crossing within a top segment
    d     = x2 - Topx   ! vector from bottom node to ray end
    d0    = x0 - Topx   ! vector from bottom node to ray origin
    tri_n = [ -( yTopSeg( 2 ) - yTopSeg( 1 ) ), xTopSeg( 2 ) - xTopSeg( 1 ), 0.0d0 ]

    IF ( CheckDiagCrossing( tri_n, d0, d ) ) THEN
       h  = -DOT_PRODUCT( d0, tri_n ) / DOT_PRODUCT( urayt, tri_n )
       EnsureStepOverTriDiagTopBdry: DO k = 1, 100
          x2 = x0 + h * urayt
          ! LP: Since this is not an exact floating-point value to step
          ! to, make sure we have stepped over the boundary.
          IF ( CheckDiagCrossing( tri_n, d0, x2 - Topx ) ) EXIT EnsureStepOverTriDiagTopBdry
          ! LP: Slightly increase h if not.
          h = h * 1.000001
       END DO EnsureStepOverTriDiagTopBdry
       IF ( k >= 100 ) WRITE( PRTFile, * ) 'EnsureStepOverTriDiagBdry did not converge'
       topRefl = .FALSE.
       botRefl = .FALSE.
       WRITE( PRTFile, * ) 'StepToBdry3D diagonal crossing h to', h, x2
    END IF

    ! triangle crossing within a bottom segment
    d     = x2 - Botx   ! vector from bottom node to ray end
    d0    = x0 - Botx   ! vector from bottom node to ray origin
    tri_n = [ -( yBotSeg( 2 ) - yBotSeg( 1 ) ), xBotSeg( 2 ) - xBotSeg( 1 ), 0.0d0 ]

    IF ( CheckDiagCrossing( tri_n, d0, d ) ) THEN
       h  = -DOT_PRODUCT( d0, tri_n ) / DOT_PRODUCT( urayt, tri_n )
       EnsureStepOverTriDiagBotBdry: DO k = 1, 100
          x2 = x0 + h * urayt
          ! LP: Since this is not an exact floating-point value to step
          ! to, make sure we have stepped over the boundary.
          IF ( CheckDiagCrossing( tri_n, d0, x2 - Botx ) ) EXIT EnsureStepOverTriDiagBotBdry
          ! LP: Slightly increase h if not.
          h = h * 1.000001
       END DO EnsureStepOverTriDiagBotBdry
       IF ( k >= 100 ) WRITE( PRTFile, * ) 'EnsureStepOverTriDiagBdry did not converge'
       topRefl = .FALSE.
       botRefl = .FALSE.
       WRITE( PRTFile, * ) 'StepToBdry3D diagonal crossing h to', h, x2
    END IF

    IF ( h < INFINITESIMAL_STEP_SIZE * Beam%deltas ) THEN        ! is it taking an infinitesimal step?
       h = INFINITESIMAL_STEP_SIZE * Beam%deltas                 ! make sure we make some motion
       x2 = x0 + h * urayt
       WRITE( PRTFile, * ) 'StepToBdry3D small step forced h to ', h, x2
       ! Recheck reflection conditions
       d = x2 - Topx ! vector from top to ray
       IF ( DOT_PRODUCT( Topn, d ) > EPSILON( d( 1 ) ) ) THEN
          topRefl = .TRUE.
       ELSE
          topRefl = .FALSE.
       END IF
       d = x2 - Botx ! vector from bottom to ray
       IF ( DOT_PRODUCT( Botn, d ) > EPSILON( d( 1 ) ) ) THEN
          botRefl = .TRUE.
          topRefl = .FALSE.
       ELSE
          botRefl = .FALSE.
       END IF
    END IF
  END SUBROUTINE StepToBdry3D

  LOGICAL FUNCTION CheckDiagCrossing( tri_n, d0, d )
    REAL (KIND=8), INTENT( IN ) :: d( 3 ), d0( 3 ), tri_n( 3 )
    
    CheckDiagCrossing = &
       ( DOT_PRODUCT( tri_n, d0 ) > 0.0d0 .AND. DOT_PRODUCT( tri_n, d ) <= 0.0d0 ) .OR. &
       ( DOT_PRODUCT( tri_n, d0 ) < 0.0d0 .AND. DOT_PRODUCT( tri_n, d ) >= 0.0d0 )
  END FUNCTION CheckDiagCrossing

END MODULE Step3DMod
