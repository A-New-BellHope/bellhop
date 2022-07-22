MODULE bdry3Dmod

  ! Loads
  ! altimetry (top bdry) and bathymetry (bottom bdry) data
  ! This version is for BELLHOP3D
  !
  ! x = coordinate of boundary
  ! t = tangent for a facet
  ! n = normal  for a facet (outward pointing)
  ! n1, n2 are normals for each of the triangles in a pair, n is selected from those
  ! Len = length of tangent (temporary variable to normalize tangent)

  USE SubTabulate
  USE monotonicMod
  USE FatalError

  IMPLICIT NONE
  SAVE
  INTEGER, PARAMETER :: ATIFile = 40, BTYFile = 41, Number_to_Echo = 21
  INTEGER            :: IsegTopx, IsegTopy, IsegBotx, IsegBoty, &
       NATIPts( 2 ), NBTYPts( 2 )
  INTEGER            :: ix, iy, IOStat, IAllocStat, iSmallStepCtr = 0
  REAL (KIND=8) :: xTopseg( 2 ), yTopseg( 2 ), xBotseg( 2 ), yBotseg( 2 ), &
       Topx( 3 ), Botx( 3 ), &   ! coordinates of corner of active rectangle
       Topn( 3 ), Botn( 3 ), &   ! tangent and normal    of active triangle
       Topxmid( 3 ), Botxmid( 3 ) ! coordinates of center of active rectangle
       ! because corners may be at big number and mess up floating point precision
  LOGICAL            :: Top_tridiag_pos, Bot_tridiag_pos ! whether in positive / n2 triangle
  REAL (KIND=8), PARAMETER :: TRIDIAG_THRESH = 3D-6
  REAL (KIND=8), PARAMETER :: big = 1E25                  ! large number used for domain termination when no altimetry/bathymetry given
  !big = sqrt( huge( Top( 1, 1 )%x ) ) / 1.0d5

  CHARACTER  (LEN=1) :: atiType, btyType

  TYPE BdryPt
     REAL (KIND=8) :: x( 3 ), t( 3 ), n( 3 ), n1( 3 ), n2( 3 ), Len, Noden( 3 ), Noden_unscaled( 3 ), z_xx, z_xy, z_yy, &
          phi_xx, phi_xy, phi_yy, kappa_xx, kappa_xy, kappa_yy
  END TYPE BdryPt

  TYPE(BdryPt), ALLOCATABLE :: Bot( :, : ), Top( :, : )

CONTAINS

  SUBROUTINE ReadATI3D( FileRoot, TopATI, DepthT, PRTFile )

    CHARACTER (LEN= 1), INTENT( IN ) :: TopATI        ! Set to '~' if altimetry is not flat
    INTEGER,            INTENT( IN ) :: PRTFile       ! unit number for print file
    REAL      (KIND=8), INTENT( IN ) :: DepthT        ! Nominal top depth
    CHARACTER (LEN=80), INTENT( IN ) :: FileRoot
    REAL (KIND=8), ALLOCATABLE :: Temp( : )
    REAL (KIND=8), ALLOCATABLE :: TopGlobalx( : ), TopGlobaly( : )

    SELECT CASE ( TopATI )
    CASE ( '~', '*' )
       WRITE( PRTFile, * ) '*********************************'
       WRITE( PRTFile, * ) 'Using top-altimetry file'

       OPEN( UNIT = ATIFile, FILE = TRIM( FileRoot ) // '.ati', STATUS = 'OLD', IOSTAT = IOStat, ACTION = 'READ' )
       IF ( IOsTAT /= 0 ) THEN
          WRITE( PRTFile, * ) 'ATIFile = ', TRIM( FileRoot ) // '.ati'
          CALL ERROUT( 'ReadATI', 'Unable to open altimetry file' )
       END IF

       READ(  ATIFile, * ) atiType
       SELECT CASE ( atiType )
       CASE ( 'R' )
          WRITE( PRTFile, * ) 'Regular grid for a 3D run'
       CASE ( 'C' )
          WRITE( PRTFile, * ) 'Regular grid for a 3D run (curvilinear)'
       CASE DEFAULT
          CALL ERROUT( 'ReadATI3D', 'Unknown option for selecting altimetry interpolation' )
       END SELECT

       ! x values
       READ(  ATIFile, * ) NatiPts( 1 )
       WRITE( PRTFile, * )
       WRITE( PRTFile, * ) 'Number of altimetry points in x-direction', NatiPts( 1 )

       ALLOCATE( TopGlobalx( MAX( NatiPts( 1 ), 3 ) ), Stat = IAllocStat )
       IF ( IAllocStat /= 0 ) &
          CALL ERROUT( 'BELLHOP3D:ReadATI3D', 'Insufficient memory for altimetry data: reduce # ati points' )

       TopGlobalx( 3 ) = -999.9
       READ(  ATIFile, * ) TopGlobalx( 1 : NatiPts( 1 ) )
       CALL SubTab( TopGlobalx, NatiPts( 1 ) )
       WRITE( PRTFile, "( 5G14.6 )" ) ( TopGlobalx( ix ), ix = 1, MIN( NatiPts( 1 ), Number_to_Echo ) )
       IF ( NatiPts( 1 ) > Number_to_Echo ) WRITE( PRTFile, "( G14.6 )" ) ' ... ', TopGlobalx( NatiPts( 1 ) )
       IF ( .NOT. monotonic( TopGlobalx, NatiPts( 1 ) ) ) THEN
          CALL ERROUT( 'BELLHOP3D:ReadATI3D', 'Altimetry X values are not monotonically increasing' )
       END IF

       ! y values
       READ(  ATIFile, * ) NatiPts( 2 )
       WRITE( PRTFile, * )
       WRITE( PRTFile, * ) 'Number of altimetry points in y-direction', NatiPts( 2 )

       ALLOCATE( TopGlobaly( MAX( NatiPts( 2 ), 3 ) ), Stat = IAllocStat )
       IF ( IAllocStat /= 0 ) &
          CALL ERROUT( 'BELLHOP3D:ReadATI3D', 'Insufficient memory for altimetry data: reduce # ati points' )

       TopGlobaly( 3 ) = -999.9
       READ(  ATIFile, * ) TopGlobaly( 1 : NatiPts( 2 ) )
       CALL SubTab( TopGlobaly, NatiPts( 2 ) )
       WRITE( PRTFile, "( 5G14.6 )" ) ( TopGlobaly( iy ), iy = 1, MIN( NatiPts( 2 ), Number_to_Echo ) )
       IF ( NatiPts( 2 ) > Number_to_Echo ) WRITE( PRTFile, "( G14.6 )"  ) ' ... ', TopGlobaly( NatiPts( 2 ) )
       IF ( .NOT. monotonic( TopGlobaly, NatiPts( 2 ) ) ) THEN
          CALL ERROUT( 'BELLHOP3D:ReadATI3D', 'Altimetry Y values are not monotonically increasing' )
       END IF

       TopGlobalx = 1000. * TopGlobalx   ! convert km to m
       TopGlobaly = 1000. * TopGlobaly

       ! z values
       ALLOCATE( Top( NatiPts( 1 ), NatiPts( 2 ) ), Temp( NatiPts( 1 ) ), Stat = IAllocStat )
       IF ( IAllocStat /= 0 ) &
            CALL ERROUT( 'BELLHOP3D:ReadATI3D', 'Insufficient memory for altimetry data: reduce # ati points' )

       WRITE( PRTFile, * )

       DO iy = 1, NatiPts( 2 )
          READ( ATIFile, * ) Top( :, iy )%x( 3 )   ! read a row of depths

          ! IF ( iy < Number_to_Echo .OR. iy == NatiPts( 2 ) ) THEN   ! echo some values
          !    WRITE( PRTFile, FMT = "(G11.3)" ) Top( :, iy )%x( 3 )
          ! END IF
          ! IF ( ANY( Top( :, iy )%x( 3 ) > DepthB ) ) THEN
          !    CALL ERROUT( 'BELLHOP3D:ReadATI3D', 'Altimetry drops below lowest point in the sound speed profile' )
          ! END IF
       END DO

       CLOSE( ATIFile )

       IF ( ANY( ISNAN( Top( :, : )%x( 3 ) ) ) ) THEN
          WRITE( PRTFile, * ) 'Warning in BELLHOP3D - ReadATI3D : The altimetry file contains a NaN'
       END IF
 
       DO ix = 1, NatiPts( 1 ) 
          DO iy = 1, NatiPts( 2 )
             Top( ix, iy )%x( 1 ) = TopGlobalx( ix )
             Top( ix, iy )%x( 2 ) = TopGlobaly( iy )
          END DO
       END DO
       
       DEALLOCATE( TopGlobalx )
       DEALLOCATE( TopGlobaly )

       CALL ComputeBdryTangentNormal( Top, 'Top' )

    CASE DEFAULT   ! no altimetry given, use SSP depth for flat top
       atiType = 'R'
       NatiPts = [ 2, 2 ]
       ALLOCATE( Top( 2, 2 ), Stat = IAllocStat )
       IF ( IAllocStat /= 0 ) CALL ERROUT( 'BELLHOP3D:ReadATI3D', 'Insufficient memory'  )

       Top( 1, 1 )%x = [ -big, -big, DepthT ]
       Top( 1, 2 )%x = [ -big,  big, DepthT ]
       Top( 2, 1 )%x = [  big, -big, DepthT ]
       Top( 2, 2 )%x = [  big,  big, DepthT ]

       Top( 1, 1 )%t  = [ 1.0, 0.0,  0.0 ]   ! tangent to top
       Top( 1, 1 )%n1 = [ 0.0, 0.0, -1.0 ]   ! outward-pointing normal
       Top( 1, 1 )%n2 = [ 0.0, 0.0, -1.0 ]   ! outward-pointing normal
       Top( 1, 2 )%t  = [ 1.0, 0.0,  0.0 ]   ! tangent to top
       Top( 1, 1 )%n1 = [ 0.0, 0.0, -1.0 ]   ! outward-pointing normal
       Top( 1, 1 )%n2 = [ 0.0, 0.0, -1.0 ]   ! outward-pointing normal
       Top( 2, 1 )%t  = [ 1.0, 0.0,  0.0 ]   ! tangent to top
       Top( 2, 1 )%n1 = [ 0.0, 0.0, -1.0 ]   ! outward-pointing normal
       Top( 2, 1 )%n2 = [ 0.0, 0.0, -1.0 ]   ! outward-pointing normal
       Top( 2, 2 )%t  = [ 1.0, 0.0,  0.0 ]   ! tangent to top
       Top( 2, 2 )%n1 = [ 0.0, 0.0, -1.0 ]   ! outward-pointing normal
       Top( 2, 2 )%n2 = [ 0.0, 0.0, -1.0 ]   ! outward-pointing normal
    END SELECT

    ! dummy TopSeg info to force GetTopSeg to search for the active segment on first call
    xTopSeg = [ +big, -big ]
    yTopSeg = [ +big, -big ]

  END SUBROUTINE ReadATI3D

  ! **********************************************************************!

  SUBROUTINE ReadBTY3D( FileRoot, BotBTY, DepthB, PRTFile )

    ! Reads in the bottom bathymetry

    CHARACTER (LEN= 1), INTENT( IN ) :: BotBTY        ! Set to '~' if bathymetry is not flat
    INTEGER,            INTENT( IN ) :: PRTFile       ! unit number for print file
    REAL      (KIND=8), INTENT( IN ) :: DepthB        ! Nominal bottom depth
    CHARACTER (LEN=80), INTENT( IN ) :: FileRoot
    REAL (KIND=8), ALLOCATABLE :: Temp( : )
    REAL (KIND=8), ALLOCATABLE :: BotGlobalx( : ), BotGlobaly( : )
 
    SELECT CASE ( BotBTY )
    CASE ( '~', '*' )
       WRITE( PRTFile, * ) '*********************************'
       WRITE( PRTFile, * ) 'Using bottom-bathymetry file'

       OPEN( UNIT = BTYFile, FILE = TRIM( FileRoot ) // '.bty', STATUS = 'OLD', IOSTAT = IOStat, ACTION = 'READ' )
       IF ( IOStat /= 0 ) THEN
         WRITE( PRTFile, * ) 'BTYFile = ', TRIM( FileRoot ) // '.bty'
         CALL ERROUT( 'ReadBTY3D', 'Unable to open bathymetry file' )
       END IF
 
       READ( BTYFile, * ) btyType

       SELECT CASE ( btyType )
       CASE ( 'R' )
          WRITE( PRTFile, * ) 'Regular grid for a 3D run'
       CASE ( 'C' )
          WRITE( PRTFile, * ) 'Regular grid for a 3D run (curvilinear)'
       CASE DEFAULT
          CALL ERROUT( 'ReadBTY3D', 'Unknown option for selecting bathymetry interpolation' )
       END SELECT

       ! x values
       READ(  BTYFile, * ) NbtyPts( 1 )
       WRITE( PRTFile, * )
       WRITE( PRTFile, * ) 'Number of bathymetry points in x-direction', NbtyPts( 1 )

       ALLOCATE( BotGlobalx( MAX( NbtyPts( 1 ), 3 ) ), Stat = IAllocStat )
       IF ( IAllocStat /= 0 ) &
          CALL ERROUT( 'BELLHOP3D:ReadBTY3D', 'Insufficient memory for bathymetry data: reduce # bty points' )

       BotGlobalx( 3 ) = -999.9
       READ(  BTYFile, * ) BotGlobalx( 1 : NbtyPts( 1 ) )
       CALL SubTab( BotGlobalx, NbtyPts( 1 ) )
       WRITE( PRTFile, "( 5G14.6 )" ) ( BotGlobalx( ix ), ix = 1, MIN( NbtyPts( 1 ), Number_to_Echo ) )
       IF ( NbtyPts( 1 ) > Number_to_Echo ) WRITE( PRTFile, "( G14.6 )" ) ' ... ', BotGlobalx( NbtyPts( 1 ) )
       IF ( .NOT. monotonic( BotGlobalx, NbtyPts( 1 ) ) ) THEN
          CALL ERROUT( 'BELLHOP3D:ReadATI3D', 'Bathymetry X values are not monotonically increasing' )
       END IF

       ! y values
       READ(  BTYFile, * ) NbtyPts( 2 )
       WRITE( PRTFile, * )
       WRITE( PRTFile, * ) 'Number of bathymetry points in y-direction', NbtyPts( 2 )

       ALLOCATE( BotGlobaly( MAX( NbtyPts( 2 ), 3 ) ), Stat = IAllocStat )
       IF ( IAllocStat /= 0 ) &
          CALL ERROUT( 'BELLHOP3D:ReadBTY3D', 'Insufficient memory for bathymetry data: reduce # bty points' )

       BotGlobaly( 3 ) = -999.9
       READ(  BTYFile, * ) BotGlobaly( 1 : NbtyPts( 2 ) )
       CALL SubTab( BotGlobaly, NbtyPts( 2 ) )
       WRITE( PRTFile, "( 5G14.6 )" ) ( BotGlobaly( iy ), iy = 1, MIN( NbtyPts( 2 ), Number_to_Echo ) )
       IF ( NbtyPts( 2 ) > Number_to_Echo ) WRITE( PRTFile, "( G14.6 )" ) ' ... ', BotGlobaly( NbtyPts( 2 ) )
       IF ( .NOT. monotonic( BotGlobaly, NbtyPts( 2 ) ) ) THEN
          CALL ERROUT( 'BELLHOP3D:ReadATI3D', 'Bathymetry Y values are not monotonically increasing' )
       END IF

       BotGlobalx = 1000. * BotGlobalx   ! convert km to m
       BotGlobaly = 1000. * BotGlobaly

       ! z values
       ALLOCATE( Bot( NbtyPts( 1 ), NbtyPts( 2 ) ), Temp( NbtyPts( 1 ) ), Stat = IAllocStat )
       IF ( IAllocStat /= 0 ) &
          CALL ERROUT( 'BELLHOP3D:ReadBTY3D', 'Insufficient memory for bathymetry data: reduce # bty points' )

       WRITE( PRTFile, * )

       DO iy = 1, NbtyPts( 2 )
          READ( BTYFile, * ) Bot( :, iy )%x( 3 )    ! read a row of depths
          ! IF ( iy < Number_to_Echo .OR. iy == NbtyPts( 2 ) ) THEN   ! echo some values
          !    WRITE( PRTFile, FMT = "(G11.3)" ) Bot( :, iy )%x( 3 )
          ! END IF
          ! IF ( ANY( Bot( :, iy )%x( 3 ) > DepthB ) ) THEN
          !    CALL ERROUT( 'BELLHOP3D:ReadBTY3D', 'Bathymetry drops below lowest point in the sound speed profile' )
          ! END IF
       END DO

       CLOSE( BTYFile )

       IF ( ANY( ISNAN( Bot( :, : )%x( 3 ) ) ) ) THEN
          WRITE( PRTFile, * ) 'Warning in BELLHOP3D - ReadBTY3D : The bathymetry file contains a NaN'
       END IF
 
       DO ix = 1, NbtyPts( 1 ) 
          DO iy = 1, NbtyPts( 2 )
             Bot( ix, iy )%x( 1 ) = BotGlobalx( ix )
             Bot( ix, iy )%x( 2 ) = BotGlobaly( iy )
          END DO
       END DO
       
       DEALLOCATE( BotGlobalx )
       DEALLOCATE( BotGlobaly )

       CALL ComputeBdryTangentNormal( Bot, 'Bot' )
    CASE DEFAULT   ! no bathymetry given, use SSP depth for flat bottom
       btyType = 'R'
       NbtyPts = [ 2, 2 ]
       ALLOCATE( Bot( 2, 2 ), Stat = IAllocStat )
       IF ( IAllocStat /= 0 ) CALL ERROUT( 'BELLHOP', 'Insufficient memory'  )

       Bot( 1, 1 )%x = [ -big, -big, DepthB ]
       Bot( 1, 2 )%x = [ -big,  big, DepthB ]
       Bot( 2, 1 )%x = [  big, -big, DepthB ]
       Bot( 2, 2 )%x = [  big,  big, DepthB ]

       Bot( 1, 1 )%t  = [ 1.0, 0.0, 0.0 ]   ! tangent to bottom
       Bot( 1, 1 )%n1 = [ 0.0, 0.0, 1.0 ]   ! outward-pointing normal
       Bot( 1, 1 )%n2 = [ 0.0, 0.0, 1.0 ]   ! outward-pointing normal
       Bot( 1, 2 )%t  = [ 1.0, 0.0, 0.0 ]   ! tangent to bottom
       Bot( 1, 2 )%n1 = [ 0.0, 0.0, 1.0 ]   ! outward-pointing normal
       Bot( 1, 2 )%n2 = [ 0.0, 0.0, 1.0 ]   ! outward-pointing normal
       Bot( 2, 1 )%t  = [ 1.0, 0.0, 0.0 ]   ! tangent to bottom
       Bot( 2, 1 )%n1 = [ 0.0, 0.0, 1.0 ]   ! outward-pointing normal
       Bot( 2, 1 )%n2 = [ 0.0, 0.0, 1.0 ]   ! outward-pointing normal
       Bot( 2, 2 )%t  = [ 1.0, 0.0, 0.0 ]   ! tangent to bottom
       Bot( 2, 2 )%n1 = [ 0.0, 0.0, 1.0 ]   ! outward-pointing normal
       Bot( 2, 2 )%n2 = [ 0.0, 0.0, 1.0 ]   ! outward-pointing normal

    END SELECT

    ! dummy BotSeg info to force GetBotSeg to search for the active segment on first call
    xBotSeg = [ +big, -big ]
    yBotSeg = [ +big, -big ]

  END SUBROUTINE ReadBTY3D

  ! **********************************************************************!

  SUBROUTINE GetTopSeg3D( x, t, isInit )

    ! Get the Top segment info (index and range interval) for XY position, x
    ! sets Topx and Topn
    ! LP: According to the original logic, if the ray escapes the box here,
    ! a warning is printed, and the other segment info is not updated. If
    ! the ray escaped the box to the negative side in either dimension, or
    ! the results are NaN below (should never happen), the ray is then
    ! terminated in TraceRay without the current step. However, if the ray
    ! escaped to the positive side, the ray is terminated as part of the
    ! normal stopping conditions, including the current step.
    ! On top of that weird asymmetry, of course the original logic allowed
    ! for steps to the edge of the same segment they are currently in,
    ! whereas we always step into the next segment according to the tangent.
    ! Finally, in 2D, the altimetry / bathymetry are extended to very large
    ! values so this never gets hit, whereas in Nx2D and 3D it is not so
    ! almost all rays hit this case.
    ! This has been changed to:
    ! - An override so being at the outer edge of a segment is OK for the
    !   last segment (just for this step)
    ! - If the position is actually outside the segment, print a warning
    !   and put it in the nearest valid segment.
    ! - Changed the stopping condition to stop if on the boundary and
    !   pointing outwards.

    INTEGER, PARAMETER :: PRTFile = 6
    REAL (KIND=8), INTENT( IN ) :: x( 3 ), t( 3 )
    LOGICAL,       INTENT( IN ) :: isInit
    INTEGER :: nx, ny
    REAL (KIND=8) :: Top_tri_n( 2 )   ! triangle normals
    REAL (KIND=8) :: over_diag_amount
    
    nx = NatiPts( 1 )
    ny = NatiPts( 2 )
    
    IsegTopx = MIN( MAX( IsegTopx, 1 ), nx - 1 )
    IsegTopy = MIN( MAX( IsegTopy, 1 ), ny - 1 )
    
    IF ( t( 1 ) >= 0.0 ) THEN
       DO WHILE ( IsegTopx >= 1 .AND. Top( IsegTopx, 1 )%x( 1 ) > x( 1 ) )
          IsegTopx = IsegTopx - 1
       END DO
       DO WHILE ( IsegTopx >= 1 .AND. IsegTopx < nx .AND. Top( IsegTopx + 1, 1 )%x( 1 ) <= x( 1 ) )
          IsegTopx = IsegTopx + 1
       END DO
    ELSE
      DO WHILE ( IsegTopx < nx .AND. Top( IsegTopx + 1, 1 )%x( 1 ) < x( 1 ) )
         IsegTopx = IsegTopx + 1
      END DO
      DO WHILE ( IsegTopx >= 1 .AND. IsegTopx < nx .AND. Top( IsegTopx, 1 )%x( 1 ) >= x( 1 ) )
         IsegTopx = IsegTopx - 1
      END DO
    END IF
    IF ( t( 2 ) >= 0.0 ) THEN
      DO WHILE ( IsegTopy >= 1 .AND. Top( 1, IsegTopy )%x( 2 ) > x( 2 ) )
         IsegTopy = IsegTopy - 1
      END DO
      DO WHILE ( IsegTopy >= 1 .AND. IsegTopy < ny .AND. Top( 1, IsegTopy + 1 )%x( 2 ) <= x( 2 ) )
         IsegTopy = IsegTopy + 1
      END DO
    ELSE
      DO WHILE ( IsegTopy < ny .AND. Top( 1, IsegTopy + 1 )%x( 2 ) < x( 2 ) )
         IsegTopy = IsegTopy + 1
      END DO
      DO WHILE ( IsegTopy >= 1 .AND. IsegTopy < ny .AND. Top( 1, IsegTopy )%x( 2 ) >= x( 2 ) )
         IsegTopy = IsegTopy - 1
      END DO
    END IF
    
    IF ( IsegTopx ==  0 .AND. Top(  1,  1 )%x( 1 ) == x( 1 ) ) IsegTopx = 1
    IF ( IsegTopx == nx .AND. Top( nx,  1 )%x( 1 ) == x( 1 ) ) IsegTopx = nx - 1
    IF ( IsegTopy ==  0 .AND. Top(  1,  1 )%x( 2 ) == x( 2 ) ) IsegTopy = 1
    IF ( IsegTopy == ny .AND. Top(  1, ny )%x( 2 ) == x( 2 ) ) IsegTopy = ny - 1
    
    IF ( IsegTopx <= 0 .OR. IsegTopx >= nx .OR. IsegTopy <= 0 .OR. IsegTopy >= ny ) THEN
       WRITE( PRTFile, * ) 'Warning: GetTopSeg3D: Top altimetry undefined above the ray, x', x
       IsegTopx = MIN( MAX( IsegTopx, 1 ), nx - 1 )
       IsegTopy = MIN( MAX( IsegTopy, 1 ), ny - 1 )
    END IF
    
    
    xTopSeg  = [ Top( IsegTopx, 1 )%x( 1 ), Top( IsegTopx + 1, 1 )%x( 1 ) ]   ! segment limits in range
    yTopSeg  = [ Top( 1, IsegTopy )%x( 2 ), Top( 1, IsegTopy + 1 )%x( 2 ) ]   ! segment limits in range
    
    Topx = Top( IsegTopx, IsegTopy )%x
    Topxmid = ( Topx + Top( IsegTopx + 1, IsegTopy + 1 )%x ) * 0.5D0
    
    ! WRITE( PRTFile, * ) 'IsegTop', IsegTopx, IsegTopy
    ! WRITE( PRTFile, * ) 'Topx x', Topx, x

    ! identify the normal based on the active triangle of a pair
    ! normal of triangle side pointing up and to the left
    Top_tri_n = [ -( yTopSeg( 2 ) - yTopSeg( 1 ) ), xTopSeg( 2 ) - xTopSeg( 1 ) ]
    Top_tri_n = Top_tri_n / NORM2( Top_tri_n )
    over_diag_amount = DOT_PRODUCT( x( 1 : 2 ) - Topxmid( 1 : 2 ), Top_tri_n )
    ! WRITE( PRTFile, * ) 'Top_tri_n over_diag_amount', Top_tri_n, over_diag_amount
    IF ( ABS( over_diag_amount ) > TRIDIAG_THRESH ) THEN
       Top_tridiag_pos = ( over_diag_amount >= 0.0D0 )
    ELSE IF ( isInit ) THEN
       Top_tridiag_pos = DOT_PRODUCT( t( 1 : 2 ), Top_tri_n ) >= 0.0D0
    END IF
    IF ( .NOT. Top_tridiag_pos ) THEN
       Topn = Top( IsegTopx, IsegTopy )%n1
    ELSE
       Topn = Top( IsegTopx, IsegTopy )%n2
    END IF

    ! if the Bot depth is bad (a NaN) then error out
    ! LP: Originally set segment to invalid and relied on later code to catch this.
    IF ( ISNAN( Topx( 3 ) ) .OR. ANY( ISNAN( Topn ) ) ) THEN
       WRITE( PRTFile, * ) 'Error: Boundary segment contains NaN!'
       CALL ERROUT( 'BELLHOP: GetTopSeg3D', 'Boundary segment contains NaN' )
    END IF
  END SUBROUTINE GetTopSeg3D

  ! **********************************************************************!

  SUBROUTINE GetBotSeg3D( x, t, isInit )
    
    ! Get the Bottom segment info (index and range interval) for XY position, x
    ! sets Botx and Botn
    ! LP: See comment in GetTopSeg3D.
    
    INTEGER, PARAMETER :: PRTFile = 6
    REAL (KIND=8), INTENT( IN ) :: x( 3 ), t( 3 )
    LOGICAL,       INTENT( IN ) :: isInit
    INTEGER :: nx, ny
    REAL (KIND=8) :: Bot_tri_n( 2 )   ! triangle normals
    REAL (KIND=8) :: over_diag_amount
    
    nx = NbtyPts( 1 )
    ny = NbtyPts( 2 )
    
    IsegBotx = MIN( MAX( IsegBotx, 1 ), nx - 1 )
    IsegBoty = MIN( MAX( IsegBoty, 1 ), ny - 1 )
    
    IF ( t( 1 ) >= 0.0 ) THEN
       DO WHILE ( IsegBotx >= 1 .AND. Bot( IsegBotx, 1 )%x( 1 ) > x( 1 ) )
          IsegBotx = IsegBotx - 1
       END DO
       DO WHILE ( IsegBotx >= 1 .AND. IsegBotx < nx .AND. Bot( IsegBotx + 1, 1 )%x( 1 ) <= x( 1 ) )
          IsegBotx = IsegBotx + 1
       END DO
    ELSE
      DO WHILE ( IsegBotx < nx .AND. Bot( IsegBotx + 1, 1 )%x( 1 ) < x( 1 ) )
         IsegBotx = IsegBotx + 1
      END DO
      DO WHILE ( IsegBotx >= 1 .AND. IsegBotx < nx .AND. Bot( IsegBotx, 1 )%x( 1 ) >= x( 1 ) )
         IsegBotx = IsegBotx - 1
      END DO
    END IF
    IF ( t( 2 ) >= 0.0 ) THEN
      DO WHILE ( IsegBoty >= 1 .AND. Bot( 1, IsegBoty )%x( 2 ) > x( 2 ) )
         IsegBoty = IsegBoty - 1
      END DO
      DO WHILE ( IsegBoty >= 1 .AND. IsegBoty < ny .AND. Bot( 1, IsegBoty + 1 )%x( 2 ) <= x( 2 ) )
         IsegBoty = IsegBoty + 1
      END DO
    ELSE
      DO WHILE ( IsegBoty < ny .AND. Bot( 1, IsegBoty + 1 )%x( 2 ) < x( 2 ) )
         IsegBoty = IsegBoty + 1
      END DO
      DO WHILE ( IsegBoty >= 1 .AND. IsegBoty < ny .AND. Bot( 1, IsegBoty )%x( 2 ) >= x( 2 ) )
         IsegBoty = IsegBoty - 1
      END DO
    END IF
    
    IF ( IsegBotx ==  0 .AND. Bot(  1,  1 )%x( 1 ) == x( 1 ) ) IsegBotx = 1
    IF ( IsegBotx == nx .AND. Bot( nx,  1 )%x( 1 ) == x( 1 ) ) IsegBotx = nx - 1
    IF ( IsegBoty ==  0 .AND. Bot(  1,  1 )%x( 2 ) == x( 2 ) ) IsegBoty = 1
    IF ( IsegBoty == ny .AND. Bot(  1, ny )%x( 2 ) == x( 2 ) ) IsegBoty = ny - 1
    
    IF ( IsegBotx <= 0 .OR. IsegBotx >= nx .OR. IsegBoty <= 0 .OR. IsegBoty >= ny ) THEN
       WRITE( PRTFile, * ) 'Warning: GetBotSeg3D: Bottom bathymetry undefined below the ray, x', x
       IsegBotx = MIN( MAX( IsegBotx, 1 ), nx - 1 )
       IsegBoty = MIN( MAX( IsegBoty, 1 ), ny - 1 )
    END IF
    
    ! WRITE( PRTFile, * ) 'IsegBot', IsegBotx, IsegBoty
    
    xBotSeg  = [ Bot( IsegBotx, 1 )%x( 1 ), Bot( IsegBotx + 1, 1 )%x( 1 ) ]   ! segment limits in range
    yBotSeg  = [ Bot( 1, IsegBoty )%x( 2 ), Bot( 1, IsegBoty + 1 )%x( 2 ) ]   ! segment limits in range
    
    Botx = Bot( IsegBotx, IsegBoty )%x
    Botxmid = ( Botx + Bot( IsegBotx + 1, IsegBoty + 1 )%x ) * 0.5D0
    
    ! identify the normal based on the active triangle of a pair
    ! normal of triangle side pointing up and to the left
    Bot_tri_n = [ -( yBotSeg( 2 ) - yBotSeg( 1 ) ), xBotSeg( 2 ) - xBotSeg( 1 ) ]
    Bot_tri_n = Bot_tri_n / NORM2( Bot_tri_n )
    over_diag_amount = DOT_PRODUCT( x( 1 : 2 ) - Botxmid( 1 : 2 ), Bot_tri_n )
    IF ( ABS( over_diag_amount ) > TRIDIAG_THRESH ) THEN
       Bot_tridiag_pos = ( over_diag_amount >= 0.0D0 )
    ELSE IF ( isInit ) THEN
       Bot_tridiag_pos = DOT_PRODUCT( t( 1 : 2 ), Bot_tri_n ) >= 0.0D0
    END IF
    IF ( .NOT. Bot_tridiag_pos ) THEN
       Botn = Bot( IsegBotx, IsegBoty )%n1
    ELSE
       Botn = Bot( IsegBotx, IsegBoty )%n2
    END IF

    ! if the Bot depth is bad (a NaN) then error out
    ! LP: Originally set segment to invalid and relied on later code to catch this.
    IF ( ISNAN( Botx( 3 ) ) .OR. ANY( ISNAN( Botn ) ) ) THEN
       WRITE( PRTFile, * ) 'Error: Boundary segment contains NaN!'
       CALL ERROUT( 'BELLHOP: GetBotSeg3D', 'Boundary segment contains NaN' )
    END IF
  END SUBROUTINE GetBotSeg3D

 ! **********************************************************************!

  SUBROUTINE ComputeBdryTangentNormal( Bdry, BotTop )

    ! Does some pre-processing on the boundary points to pre-compute segment
    ! lengths  (%Len),
    ! tangents (%t, %nodet),
    ! normals  (%n, %noden), and
    ! curvatures (%kappa)
    !
    ! The boundary is also extended with a constant depth to infinity to cover cases where the ray
    ! exits the domain defined by the user

    INTEGER                          :: NPts( 2 ) = [ 0, 0 ]
    REAL      (KIND=8)               :: p1( 3 ), p2( 3 ), p3( 3 ), p4( 3 ), U( 3 ), V( 3 )
    REAL      (KIND=8)               :: n1( 3 ), n2( 3 )      ! normal vectors to the pair of triangles
    REAL      (KIND=8)               :: tvec( 3 ), Len
    TYPE(BdryPt)                     :: Bdry( :, : )
    CHARACTER (LEN=3),  INTENT( IN ) :: BotTop           ! Flag indicating bottom or top reflection
    CHARACTER (LEN=2)                :: CurvilinearFlag = '-'
    REAL      (KIND=8)               :: mx, my, n( 3 )

    SELECT CASE ( BotTop )
    CASE ( 'Bot' )
       NPts = NbtyPts
       CurvilinearFlag = btyType
    CASE ( 'Top' )
       NPts = NatiPts
       CurvilinearFlag = atiType
    END SELECT

    ! normals on triangle faces
    DO ix = 1, NPts( 1 ) - 1
       DO iy = 1, NPts( 2 ) - 1
          ! coordinates of corner nodes, moving counter-clockwise around the rectangle
          p1 = Bdry( ix,     iy     )%x
          p2 = Bdry( ix + 1, iy     )%x
          p3 = Bdry( ix + 1, iy + 1 )%x
          p4 = Bdry( ix,     iy + 1 )%x

          ! edges for triangle 1
          U = p2 - p1   ! tangent along one edge
          V = p3 - p1   ! tangent along another edge

          ! normal vector is the cross-product of the edge tangents
          n1( 1 ) = U( 2 ) * V( 3 ) - U( 3 ) * V( 2 )
          n1( 2 ) = U( 3 ) * V( 1 ) - U( 1 ) * V( 3 )
          n1( 3 ) = U( 1 ) * V( 2 ) - U( 2 ) * V( 1 )
          IF ( BotTop == 'Top' ) n1 = -n1

          Bdry( ix, iy )%n1 = n1 / NORM2( n1 )   ! scale to make it a unit normal

          ! edges for triangle 2
          U = p3 - p1   ! tangent along one edge
          V = p4 - p1   ! tangent along another edge

          ! normal vector is the cross-product of the edge tangents
          n2( 1 ) = U( 2 ) * V( 3 ) - U( 3 ) * V( 2 )
          n2( 2 ) = U( 3 ) * V( 1 ) - U( 1 ) * V( 3 )
          n2( 3 ) = U( 1 ) * V( 2 ) - U( 2 ) * V( 1 )
          IF ( BotTop == 'Top' ) n2 = -n2
          
          Bdry( ix, iy )%n2 = n2 / NORM2( n2 )   ! scale to make it a unit normal
       
       END DO
    END DO

    ! normals at nodes
    ! use forward, centered, or backward difference formulas
    DO ix = 1, NPts( 1 )
       DO iy = 1, NPts( 2 )
          IF ( ix == 1 ) THEN
             mx = ( Bdry( ix + 1, iy     )%x( 3 ) - Bdry( ix    , iy     )%x( 3 ) ) / &
                  ( Bdry( ix + 1, iy     )%x( 1 ) - Bdry( ix    , iy     )%x( 1 ) )
          ELSE IF ( ix == Npts( 1 ) ) THEN
             mx = ( Bdry( ix    , iy     )%x( 3 ) - Bdry( ix - 1, iy     )%x( 3 ) ) / &
                  ( Bdry( ix    , iy     )%x( 1 ) - Bdry( ix - 1, iy     )%x( 1 ) )
          ELSE
             mx = ( Bdry( ix + 1, iy     )%x( 3 ) - Bdry( ix - 1, iy     )%x( 3 ) ) / &
                  ( Bdry( ix + 1, iy     )%x( 1 ) - Bdry( ix - 1, iy     )%x( 1 ) )
          END IF

          IF ( iy == 1 ) THEN
             my = ( Bdry( ix    , iy + 1 )%x( 3 ) - Bdry( ix    , iy     )%x( 3 ) ) / &
                  ( Bdry( ix    , iy + 1 )%x( 2 ) - Bdry( ix    , iy     )%x( 2 ) )
          ELSE IF ( iy == Npts( 2 ) ) THEN
             my = ( Bdry( ix    , iy     )%x( 3 ) - Bdry( ix    , iy - 1 )%x( 3 ) ) / &
                  ( Bdry( ix    , iy     )%x( 2 ) - Bdry( ix    , iy - 1 )%x( 2 ) )
          ELSE
             my = ( Bdry( ix    , iy + 1 )%x( 3 ) - Bdry( ix    , iy - 1 )%x( 3 ) ) / &
                  ( Bdry( ix    , iy + 1 )%x( 2 ) - Bdry( ix    , iy - 1 )%x( 2 ) )
          END IF

          n = [ -mx, -my, 1.0D0 ]   ! this a normal to the surface

          IF ( ix < NPts( 1 ) .AND. iy < NPts( 2 ) ) THEN
             ! xx term
             Bdry( ix, iy )%phi_xx = atan2( n( 3 ), n( 1 ) )   ! this is the angle at each node

             ! xy term
             tvec = Bdry( ix + 1, iy + 1 )%x - Bdry( ix, iy )%x
             Len  = SQRT( tvec( 1 ) ** 2 + tvec( 2 ) ** 2 )
             tvec = tvec / Len
             Bdry( ix, iy )%phi_xy = atan2( n( 3 ), n( 1 ) * tvec( 1 ) + n( 2 ) * tvec( 2 ) )   ! this is the angle at each node

             ! yy term
             Bdry( ix, iy )%phi_yy = atan2( n( 3 ), n( 2 ) )   ! this is the angle at each node
          END IF

          Bdry( ix, iy )%Noden_unscaled = n
          Bdry( ix, iy )%Noden = n / NORM2( n )          
       END DO
    END DO

    IF ( CurvilinearFlag( 1 : 1 ) == 'C' ) THEN ! curvilinear option: compute derivative as centered difference between two nodes

       ! compute curvatures in each segment
       ! ALLOCATE( phi( NPts ), Stat = IAllocStat )

       ! - sign below because the node normal = ( -mx, -my, 1 )
       DO ix = 1, NPts( 1 ) - 1
          DO iy = 1, NPts( 2 ) - 1
             ! z_xx (difference in x of z_x)

             Bdry( ix, iy )%z_xx = -( Bdry( ix + 1, iy     )%Noden_unscaled( 1 ) - Bdry( ix, iy )%Noden_unscaled( 1 ) ) / &
                                    ( Bdry( ix + 1, iy     )%x(              1 ) - Bdry( ix, iy )%x(              1 ) )

             tvec = Bdry( ix + 1, iy )%x - Bdry( ix, iy )%x
             Len  = SQRT( tvec( 1 ) ** 2 + tvec( 3 ) ** 2 )
             Bdry( ix, iy )%kappa_xx = ( Bdry( ix + 1, iy )%phi_xx - Bdry( ix, iy )%phi_xx ) / Len ! this is curvature = dphi/ds

             ! z_xy (difference in y of z_x)

             Bdry( ix, iy )%z_xy = -( Bdry( ix    , iy + 1 )%Noden_unscaled( 1 ) - Bdry( ix, iy )%Noden_unscaled( 1 ) ) / &
                                    ( Bdry( ix    , iy + 1 )%x(              2 ) - Bdry( ix, iy )%x(              2 ) )

             tvec = Bdry( ix + 1, iy + 1 )%x - Bdry( ix, iy )%x
             Len  = SQRT( tvec( 1 ) ** 2 + tvec( 2 ) ** 2 + tvec( 3 ) ** 2 )
             Bdry( ix, iy )%kappa_xy = ( Bdry( ix + 1, iy + 1 )%phi_xy - Bdry( ix, iy )%phi_xy ) / Len ! this is curvature = dphi/ds

             ! new
             tvec = Bdry( ix, iy + 1 )%x - Bdry( ix, iy )%x
             Len  = SQRT( tvec( 2 ) ** 2 + tvec( 3 ) ** 2 )
             Bdry( ix, iy )%kappa_xy = ( Bdry( ix, iy + 1 )%phi_xx - Bdry( ix, iy )%phi_xx ) / Len ! this is curvature = dphi/ds

             ! z_yy (difference in y of z_y)

             Bdry( ix, iy )%z_yy = -( Bdry( ix    , iy + 1 )%Noden_unscaled( 2 ) - Bdry( ix, iy )%Noden_unscaled( 2 ) ) / &
                                    ( Bdry( ix    , iy + 1 )%x(              2 ) - Bdry( ix, iy )%x(              2 ) )

             tvec = Bdry( ix, iy + 1 )%x - Bdry( ix, iy )%x
             Len  = SQRT( tvec( 2 ) ** 2 + tvec( 3 ) ** 2 )
             Bdry( ix, iy )%kappa_yy = ( Bdry( ix, iy + 1 )%phi_yy - Bdry( ix, iy )%phi_yy ) / Len ! this is curvature = dphi/ds

             ! introduce Len factor per Eq. 4.4.18 in Cerveny's book
             Len = NORM2( Bdry( ix, iy )%Noden_unscaled )
             Bdry( ix, iy )%z_xx = Bdry( ix, iy )%z_xx / Len
             Bdry( ix, iy )%z_xy = Bdry( ix, iy )%z_xy / Len
             Bdry( ix, iy )%z_yy = Bdry( ix, iy )%z_yy / Len
          END DO
       END DO
    ELSE
       Bdry%z_xx = 0
       Bdry%z_xy = 0
       Bdry%z_yy = 0

       Bdry%kappa_xx = 0
       Bdry%kappa_xy = 0
       Bdry%kappa_yy = 0
    END IF
!!$    write( *, * ) 'ix=1', Bdry( 1, : )%kappa_xx
!!$    write( *, * ) 'ix=1', Bdry( 1, : )%kappa_xy
!!$    write( *, * ) 'ix=1', Bdry( 1, : )%kappa_yy
!!$    write( *, * ) 'iy=1', Bdry( :, 1 )%kappa_xx
!!$    write( *, * ) 'iy=1', Bdry( :, 1 )%kappa_xy
!!$    write( *, * ) 'iy=1', Bdry( :, 1 )%kappa_yy
!!$    write( *, * ) 'D'
!!$    write( *, * ) Bdry( :, : )%z_xx
!!$    write( *, * ) Bdry( :, : )%z_xy
!!$    write( *, * ) Bdry( :, : )%z_yy
!!$
!!$    write( *, * ) 'kappa'
!!$    write( *, * ) Bdry( :, : )%kappa_xx
!!$    write( *, * ) Bdry( :, : )%kappa_xy
!!$    write( *, * ) Bdry( :, : )%kappa_yy

  END SUBROUTINE ComputeBdryTangentNormal

END MODULE bdry3Dmod

 
