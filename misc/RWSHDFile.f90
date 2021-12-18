MODULE RWSHDFile

  ! routines to read or write the SHDFile

  USE SourceReceiverPositions
  USE FatalError

  IMPLICIT NONE
  INTEGER, PARAMETER, PRIVATE :: SHDFile = 25
  INTEGER LRecl

  ! variables taken from SourceReceiverPositions:
  ! freqVec vector of frequencies
  ! theta   vector of bearing lines,   theta( 1 : Ntheta )
  ! Sz      vector of source   depths, Sz(    1 : NSz    )
  ! Rz      vector of receiver depths, Rz(    1 : NRz    )
  ! Rr      vector of receiver ranges, Rr(    1 : NRr    )

CONTAINS
  SUBROUTINE ReadHeader( FileName, Title, Atten, PlotType )

    ! Read header from disk file
    ! This is not used anywhere in the Fortran code for the Acoustics Toolbox, since Matlab is used to process SHDFiles

    ! FileName is a SHDFIL for complex pressure or a GRNFIL for a Green's function
    ! Title   arbitrary title

    INTEGER, PARAMETER                :: PRTFile = 6
    REAL,               INTENT( OUT ) :: Atten           ! stabilizing attenuation for SCOOTER FFP runs
    CHARACTER (LEN=80), INTENT( OUT ) :: Title, FileName
    CHARACTER (LEN=10), INTENT( OUT ) :: PlotType
    INTEGER                           :: IAllocStat, IOStat

    ! Open file, read header
    IF ( FileName( 1 : 1 ) == ' ' ) FileName = 'SHDFIL'

    ! INQUIRE( FILE = FileName, RECL = IRECL )
    OPEN( UNIT = SHDFile,   FILE = FileName, STATUS = 'OLD', ACCESS = 'DIRECT', FORM = 'UNFORMATTED', RECL = 4, &
         IOSTAT = IOStaT, ACTION = 'READ' )
    IF ( IOStat /= 0 ) CALL ERROUT( 'ReadHeader', 'Unable to open shade file' )

    READ( SHDFile, REC = 1 ) LRecl
    CLOSE( UNIT = SHDFile )
    OPEN(  UNIT = SHDFile,   FILE = FileName, STATUS = 'OLD', ACCESS = 'DIRECT', FORM = 'UNFORMATTED', RECL = 4 * LRecl )

    READ( SHDFile, REC = 1  ) LRecl, Title
    READ( SHDFile, REC = 2  ) PlotType
    READ( SHDFile, REC = 3  ) Nfreq, Pos%Ntheta, Pos%NSx, Pos%NSy, Pos%NSz, Pos%NRz, Pos%NRr, atten

    ALLOCATE( freqVec( Nfreq ), Pos%Sz( Pos%NSz ), Pos%Rz( Pos%NRz ), Pos%Rr( Pos%NRr ), Pos%theta( Pos%Ntheta ), Stat = IAllocStat)
    IF ( IAllocStat /= 0 ) CALL ERROUT( 'ReadHeader', 'Too many source/receiver combinations' )

    READ( SHDFile, REC = 4  ) freqVec
    READ( SHDFile, REC = 5  ) Pos%theta
    READ( SHDFile, REC = 6  ) Pos%Sx
    READ( SHDFile, REC = 7  ) Pos%Sy
    READ( SHDFile, REC = 8  ) Pos%Sz
    READ( SHDFile, REC = 9  ) Pos%Rz
    READ( SHDFile, REC = 10 ) Pos%Rr

    ! Pos%deltaR = Pos%r( Pos%NRr ) - Pos%r( Pos%NRr - 1 )

  END SUBROUTINE ReadHeader

  !**********************************************************************!

  SUBROUTINE WriteHeader( FileName, Title, freq0, Atten, PlotType )

    ! Write header to disk file

    REAL,      INTENT( IN ) :: freq0, Atten      ! Nominal frequency, stabilizing attenuation (for wavenumber integration only)
    CHARACTER, INTENT( IN ) :: FileName*( * )    ! Name of the file (could be a shade file or a Green's function file)
    CHARACTER, INTENT( IN ) :: Title*( * )       ! Arbitrary title
    CHARACTER, INTENT( IN ) :: PlotType*( 10 )   ! 

    ! receiver bearing angles
    IF ( .NOT. ALLOCATED( Pos%theta ) ) THEN
       ALLOCATE( Pos%theta( 1 ) )
       Pos%theta( 1 ) = 0   ! dummy bearing angle
       Pos%Ntheta     = 1
    END IF

    ! source x-coordinates
    IF ( .NOT. ALLOCATED( Pos%Sx ) ) THEN
       ALLOCATE( Pos%Sx( 1 ) )
       Pos%sx( 1 ) = 0      ! dummy x-coordinate
       Pos%NSx     = 1
    END IF

    ! source y-coordinates
    IF ( .NOT. ALLOCATED( Pos%Sy ) ) THEN
       ALLOCATE( Pos%Sy( 1 ) )
       Pos%sy( 1 ) = 0      ! dummy y-coordinate
       Pos%NSy     = 1
    END IF

    IF ( PlotType( 1 : 2 ) /= 'TL' ) THEN
       ! MAX( 41, ... ) below because Title is already 40 words (or 80 bytes)
       LRecl = MAX( 41, 2 * Nfreq, Pos%Ntheta, Pos%NSx, Pos%NSy, Pos%NSz, Pos%NRz, 2 * Pos%NRr )   ! words/record (NRr doubled for complex pressure storage)

       OPEN ( FILE = FileName, UNIT = SHDFile, STATUS = 'REPLACE', ACCESS = 'DIRECT', RECL = 4 * LRecl, FORM = 'UNFORMATTED')
       WRITE( SHDFile, REC = 1  ) LRecl, Title( 1 : 80 )
       WRITE( SHDFile, REC = 2  ) PlotType
       WRITE( SHDFile, REC = 3  ) Nfreq, Pos%Ntheta, Pos%NSx, Pos%NSy, Pos%NSz, Pos%NRz, Pos%NRr, freq0, atten
       WRITE( SHDFile, REC = 4  ) freqVec(   1 : Nfreq )
       WRITE( SHDFile, REC = 5  ) Pos%theta( 1 : Pos%Ntheta )

       WRITE( SHDFile, REC = 6  ) Pos%Sx( 1 : Pos%NSx )
       WRITE( SHDFile, REC = 7  ) Pos%Sy( 1 : Pos%NSy )
       WRITE( SHDFile, REC = 8  ) Pos%Sz( 1 : Pos%NSz )

       WRITE( SHDFile, REC = 9  ) Pos%Rz( 1 : Pos%NRz )
       WRITE( SHDFile, REC = 10 ) Pos%Rr( 1 : Pos%NRr )

    ELSE   ! compressed format for TL from FIELD3D
       LRecl = MAX( 41, 2 * Nfreq, Pos%Ntheta, Pos%NSz, Pos%NRz, 2 * Pos%NRr )   ! words/record (NR doubled for complex pressure storage)

       OPEN ( FILE = FileName, UNIT = SHDFile, STATUS = 'REPLACE', ACCESS = 'DIRECT', RECL = 4 * LRecl, FORM = 'UNFORMATTED')
       WRITE( SHDFile, REC = 1  ) LRecl, Title( 1 : 80 )
       WRITE( SHDFile, REC = 2  ) PlotType
       WRITE( SHDFile, REC = 3  ) Nfreq, Pos%Ntheta, Pos%NSx, Pos%NSy, Pos%NSz, Pos%NRz, Pos%NRr, freq0, atten
       WRITE( SHDFile, REC = 4  ) freqVec(   1 : Nfreq )
       WRITE( SHDFile, REC = 5  ) Pos%theta( 1 : Pos%Ntheta )

       WRITE( SHDFile, REC = 6  ) Pos%Sx( 1 ), Pos%Sx( Pos%NSx )
       WRITE( SHDFile, REC = 7  ) Pos%Sy( 1 ), Pos%Sy( Pos%NSy )
       WRITE( SHDFile, REC = 8  ) Pos%Sz( 1 : Pos%NSz )

       WRITE( SHDFile, REC = 9  ) Pos%Rz( 1 : Pos%NRz )
       WRITE( SHDFile, REC = 10 ) Pos%Rr( 1 : Pos%NRr )
    END IF

  END SUBROUTINE WriteHeader

  !**********************************************************************!

  SUBROUTINE WriteField( P, NRz, NRr, IRec )

    ! Write the field to disk

    INTEGER, INTENT( IN )    :: NRz, NRr        ! Number of receiver depths, ranges
    COMPLEX, INTENT( IN )    :: P( NRz, NRr )   ! Pressure field
    INTEGER, INTENT( INOUT ) :: iRec            ! last record read
    INTEGER                  :: iRz

    DO iRz = 1, NRz
       iRec = iRec + 1
       WRITE( SHDFile, REC = iRec ) P( iRz, : )
    END DO

  END SUBROUTINE WriteField

END MODULE RWSHDFile
