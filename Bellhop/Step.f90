MODULE Step

  USE bellhopMod
  USE sspMod
  IMPLICIT NONE
  
  REAL (KIND=8), PARAMETER, PRIVATE :: REL_SNAP = 1.0d-6 ! Relative error/distance from a boundary to snap to it
CONTAINS

  SUBROUTINE Step2D( ray0, ray2, Topx, Topn, Botx, Botn )

    ! Does a single step along the ray
    ! x denotes the ray coordinate, (r,z)
    ! t denotes the scaled tangent to the ray (previously (rho, zeta))
    ! c * t would be the unit tangent

    USE BdryMod
    TYPE( ray2DPt )    :: ray0, ray1, ray2
    REAL (KIND=8 ), INTENT( IN ) :: Topx( 2 ), Topn( 2 ), Botx( 2 ), Botn( 2 )
    INTEGER            :: iSegz0, iSegr0
    REAL     (KIND=8 ) :: gradc0( 2 ), gradc1( 2 ), gradc2( 2 ), &
         c0, cimag0, crr0, crz0, czz0, csq0, cnn0_csq0, &
         c1, cimag1, crr1, crz1, czz1, csq1, cnn1_csq1, &
         c2, cimag2, crr2, crz2, czz2, urayt0( 2 ), urayt1( 2 ), &
         h, halfh, hw0, hw1, ray2n( 2 ), RM, RN, gradcjump( 2 ), cnjump, csjump, w0, w1, rho 
         
    WRITE( PRTFile, * )
    WRITE( PRTFile, * ) 'ray0 x t p q', ray0%x, ray0%t, ray0%p, ray0%q
    ! WRITE( PRTFile, * ) 'iSegz iSegr', iSegz, iSegr
    
    IF ( ray0%x( 1 ) > 420.0 ) THEN
       STOP 'Enough'
    END IF

    ! The numerical integrator used here is a version of the polygon (a.k.a. midpoint, leapfrog, or Box method), and similar
    ! to the Heun (second order Runge-Kutta method).
    ! However, it's modified to allow for a dynamic step change, while preserving the second-order accuracy).

    ! *** Phase 1 (an Euler step)

    CALL EvaluateSSP( ray0%x, c0, cimag0, gradc0, crr0, crz0, czz0, rho, freq, 'TAB' )

    csq0      = c0 * c0
    cnn0_csq0 = crr0 * ray0%t( 2 )**2 - 2.0 * crz0 * ray0%t( 1 ) * ray0%t( 2 ) + czz0 * ray0%t( 1 )**2
    iSegz0    = iSegz     ! make note of current layer
    iSegr0    = iSegr

    h = Beam%deltas       ! initially set the step h, to the basic one, deltas
    urayt0 = c0 * ray0%t  ! unit tangent

    CALL ReduceStep2D( ray0%x, urayt0, iSegz0, iSegr0, Topx, Topn, Botx, Botn, h ) ! reduce h to land on boundary
    halfh = 0.5 * h   ! first step of the modified polygon method is a half step

    ray1%x = ray0%x + halfh * urayt0
    ray1%t = ray0%t - halfh * gradc0 / csq0
    ray1%p = ray0%p - halfh * cnn0_csq0 * ray0%q
    ray1%q = ray0%q + halfh * c0        * ray0%p

    ! WRITE( PRTFile, * ) 'ray1 x t p q', ray1%x, ray1%t, ray1%p, ray1%q

    ! *** Phase 2

    CALL EvaluateSSP( ray1%x, c1, cimag1, gradc1, crr1, crz1, czz1, rho, freq, 'TAB' )
    csq1      = c1 * c1
    cnn1_csq1 = crr1 * ray1%t( 2 )**2 - 2.0 * crz1 * ray1%t( 1 ) * ray1%t( 2 ) + czz1 * ray1%t( 1 )**2

    ! The Munk test case with a horizontally launched ray caused problems.
    ! The ray vertexes on an interface and can ping-pong around that interface.
    ! Have to be careful in that case about big changes to the stepsize (that invalidate the leap-frog scheme) in phase II.
    ! A modified Heun or Box method could also work.

    urayt1 = c1 * ray1%t   ! unit tangent

    CALL ReduceStep2D( ray0%x, urayt1, iSegz0, iSegr0, Topx, Topn, Botx, Botn, h ) ! reduce h to land on boundary

    ! WRITE( PRTFile, * ) 'urayt1', urayt1, 'out h 2', h

    ! use blend of f' based on proportion of a full step used.
    w1  = h / ( 2.0d0 * halfh )
    w0  = 1.0d0 - w1
    hw0 = h * w0
    hw1 = h * w1
    
    ! WRITE( PRTFile, * ) 'w1 w0 hw0 hw1', w1, w0, hw0, hw1

    ray2%x   = ray0%x   + hw0 * urayt0              + hw1 * urayt1
    ray2%t   = ray0%t   - hw0 * gradc0 / csq0       - hw1 * gradc1 / csq1
    ray2%p   = ray0%p   - hw0 * cnn0_csq0 * ray0%q  - hw1 * cnn1_csq1 * ray1%q
    ray2%q   = ray0%q   + hw0 * c0        * ray0%p  + hw1 * c1        * ray1%p
    ray2%tau = ray0%tau + hw0 / CMPLX( c0, cimag0, KIND=8 ) + hw1 / CMPLX( c1, cimag1, KIND=8 )

    WRITE( PRTFile, * ) 'ray2 x t p q tau', ray2%x, ray2%t, ray2%p, ray2%q, ray2%tau

    ray2%Amp       = ray0%Amp
    ray2%Phase     = ray0%Phase
    ray2%NumTopBnc = ray0%NumTopBnc
    ray2%NumBotBnc = ray0%NumBotBnc
    
    ! ReduceStep2D aims to put the position exactly on an interface. However,
    ! since it is not modifying the position directly, but computing a step
    ! size which is then used to modify the position, rounding errors may cause
    ! the final position to be slightly before or after the interface. If
    ! before, the next (very small) step takes the ray over the interface; if
    ! after, no extra small step is needed. Since rounding errors are
    ! effectively random, the resulting trajectory receives or doesn't receive
    ! an extra step randomly. While the extra step is very small, it slightly
    ! changes all the values after it, potentially leading to divergence later.
    ! So, if there is an interface very close, snap to it.
    IF ( ABS( urayt0( 2 ) ) > EPSILON( urayt0( 2 ) ) .OR. ABS( urayt1( 2 ) ) > EPSILON( urayt1( 2 ) ) ) THEN
       ! Ray has some depth movement component
       CALL StepSnapTo( ray2%x( 2 ), SSP%z( iSegz0     ) )
       CALL StepSnapTo( ray2%x( 2 ), SSP%z( iSegz0 + 1 ) )
    END IF
    IF ( ABS( DOT_PRODUCT( Topn, ray2%x - Topx )) < REL_SNAP ) THEN
       ! On the top surface. If the top surface is slanted, doing a correction
       ! here won't necessarily be any more accurate than the original
       ! computation. But if it's flat, like usually, we want the exact top
       ! depth value.
       IF ( ABS( Topn( 1 ) ) < EPSILON( Topn( 1 ) ) ) THEN
          ray2%x( 2 ) = Topx( 2 )
          WRITE( PRTFile, * ) 'Snapped at top', ray2%x
       END IF
    END IF
    IF ( ABS( DOT_PRODUCT( Botn, ray2%x - Botx )) < REL_SNAP ) THEN
       ! On the bottom surface. If the bottom surface is slanted, doing a
       ! correction here won't necessarily be any more accurate than the
       ! original computation. But if it's flat, like usually, we want the
       ! exact bottom depth value.
       IF ( ABS( Botn( 1 ) ) < EPSILON( Botn( 1 ) ) ) THEN
          ray2%x( 2 ) = Topx( 2 )
          WRITE( PRTFile, * ) 'Snapped at bottom', ray2%x
       END IF
    END IF
    IF ( ABS( urayt0( 1 ) ) > EPSILON( urayt0( 1 ) ) .OR. ABS( urayt1( 1 ) ) > EPSILON( urayt1( 1 ) ) ) THEN
       ! Ray has some range movement component
       CALL StepSnapTo( ray2%x( 1 ), rTopSeg( 1 ) ) ! top seg minimum
       CALL StepSnapTo( ray2%x( 1 ), rTopSeg( 2 ) ) ! top seg maximum
       CALL StepSnapTo( ray2%x( 1 ), rBotSeg( 1 ) ) ! bot seg minimum
       CALL StepSnapTo( ray2%x( 1 ), rBotSeg( 2 ) ) ! bot seg maximum
       IF ( SSP%Type == 'Q' ) THEN
          CALL StepSnapTo( ray2%x( 1 ), SSP%Seg%r( iSegz0     ) )
          CALL StepSnapTo( ray2%x( 1 ), SSP%Seg%r( iSegz0 + 1 ) )
       END IF
    END IF

    ! If we crossed an interface, apply jump condition

    CALL EvaluateSSP( ray2%x, c2, cimag2, gradc2, crr2, crz2, czz2, rho, freq, 'TAB' )
    ray2%c = c2

    IF ( iSegz /= iSegz0 .OR. iSegr /= iSegr0 ) THEN
       gradcjump =  gradc2 - gradc0
       ray2n     = [ -ray2%t( 2 ), ray2%t( 1 ) ]   ! ray normal

       cnjump    = DOT_PRODUCT( gradcjump, ray2n  )
       csjump    = DOT_PRODUCT( gradcjump, ray2%t )

       IF ( iSegz /= iSegz0 ) THEN         ! crossing in depth
          RM = +ray2%t( 1 ) / ray2%t( 2 )  ! this is tan( alpha ) where alpha is the angle of incidence
       ELSE                                ! crossing in range
          RM = -ray2%t( 2 ) / ray2%t( 1 )  ! this is tan( alpha ) where alpha is the angle of incidence
       END IF

       RN     = RM * ( 2 * cnjump - RM * csjump ) / c2
       ray2%p = ray2%p - ray2%q * RN

    END IF

  END SUBROUTINE Step2D

  ! **********************************************************************!

  SUBROUTINE ReduceStep2D( x0, urayt, iSegz0, iSegr0, Topx, Topn, Botx, Botn, h )

    ! calculate a reduced step size, h, that lands on any points where the environment changes

    USE BdryMod
    INTEGER,       INTENT( IN    ) :: iSegz0, iSegr0                             ! SSP layer the ray is in
    REAL (KIND=8), INTENT( IN    ) :: x0( 2 ), urayt( 2 )                        ! ray coordinate and tangent
    REAL (KIND=8), INTENT( IN    ) :: Topx( 2 ), Topn( 2 ), Botx( 2 ), Botn( 2 ) ! Top, bottom coordinate and normal
    REAL (KIND=8), INTENT( INOUT ) :: h                                          ! reduced step size 
    REAL (KIND=8)                  :: x( 2 ), d( 2 ), d0( 2 ), h1, h2, h3, h4, rSeg( 2 )

    ! Detect interface or boundary crossing and reduce step, if necessary, to land on that crossing.
    ! Keep in mind possibility that user put source right on an interface
    ! and that multiple events can occur (crossing interface, top, and bottom in a single step).

    x = x0 + h * urayt ! make a trial step

    ! interface crossing in depth
    h1 = huge( h1 )
    IF ( ABS( urayt( 2 ) ) > EPSILON( h1 ) ) THEN
       IF      ( SSP%z( iSegz0     ) > x(  2 ) ) THEN
          h1 = ( SSP%z( iSegz0     ) - x0( 2 ) ) / urayt( 2 )
       ELSE IF ( SSP%z( iSegz0 + 1 ) < x(  2 ) ) THEN
          h1 = ( SSP%z( iSegz0 + 1 ) - x0( 2 ) ) / urayt( 2 )
       END IF
    END IF

    ! top crossing
    h2 = huge( h2 )
    d  = x - Topx              ! vector from top to ray
    IF ( DOT_PRODUCT( Topn, d ) > EPSILON( h2 ) ) THEN
       d0  = x0 - Topx         ! vector from top    node to ray origin
       h2 = -DOT_PRODUCT( d0, Topn ) / DOT_PRODUCT( urayt, Topn )
    END IF

    ! bottom crossing
    h3 = huge( h3 )
    d  = x - Botx              ! vector from bottom to ray
    IF ( DOT_PRODUCT( Botn, d ) > EPSILON( h2 ) ) THEN
       d0  = x0 - Botx         ! vector from bottom node to ray origin
       h3 = -DOT_PRODUCT( d0, Botn ) / DOT_PRODUCT( urayt, Botn )
    END IF

    ! top or bottom segment crossing in range
    rSeg( 1 ) = MAX( rTopSeg( 1 ), rBotSeg( 1 ) )
    rSeg( 2 ) = MIN( rTopSeg( 2 ), rBotSeg( 2 ) )

    IF ( SSP%Type == 'Q' ) THEN
       rSeg( 1 ) = MAX( rSeg( 1 ), SSP%Seg%r( iSegr0     ) )
       rSeg( 2 ) = MIN( rSeg( 2 ), SSP%Seg%r( iSegr0 + 1 ) )
    END IF

    h4 = huge( h4 )
    IF ( ABS( urayt( 1 ) )  > EPSILON( h4 ) ) THEN
       IF       ( x(  1 ) < rSeg( 1 ) ) THEN
          h4 = -( x0( 1 ) - rSeg( 1 ) ) / urayt( 1 )
       ELSE IF  ( x(  1 ) > rSeg( 2 ) ) THEN
          h4 = -( x0( 1 ) - rSeg( 2 ) ) / urayt( 1 )
       END IF
    END IF

    h = MIN( h, h1, h2, h3, h4 )           ! take limit set by shortest distance to a crossing
    IF ( h < 1.0d-4 * Beam%deltas ) THEN   ! is it taking an infinitesimal step?
       h = 1.0d-5 * Beam%deltas            ! make sure we make some motion
       iSmallStepCtr = iSmallStepCtr + 1   ! keep a count of the number of sequential small steps
    ELSE
       iSmallStepCtr = 0   ! didn't do a small step so reset the counter
    END IF

  END SUBROUTINE ReduceStep2D
  
  ! **********************************************************************!
  
  SUBROUTINE StepSnapTo( x, tgt )
    REAL (KIND=8), INTENT( INOUT ) :: x
    REAL (KIND=8), INTENT( IN    ) :: tgt
    
    IF ( ABS( ( x - tgt ) / x ) < REL_SNAP ) THEN
       WRITE( PRTFile, * ) 'Snapping to', tgt
       x = tgt
    END IF
    
  END SUBROUTINE StepSnapTo

END MODULE Step
