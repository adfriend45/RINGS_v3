!======================================================================!
MODULE STATE
!----------------------------------------------------------------------!
USE DOUBLE
USE CONTROL
!----------------------------------------------------------------------!
IMPLICIT NONE
!----------------------------------------------------------------------!
REAL(DP), DIMENSION (nfi, ncells_max) :: L     ! Radial cell length (�m)
REAL(DP), DIMENSION (nfi, ncells_max) :: M     ! Cell DM        (mg[DM])
REAL(DP), DIMENSION (nfi, ncells_max) :: rasym ! Asymmetry on enlarg.
REAL(DP), DIMENSION (nfi, ncells_max) :: I     ! Div. inhib. vol  (�m^3)
REAL(DP), DIMENSION (nfi, ncells_max) :: D ! Cell distance from phloem (�m)
!----------------------------------------------------------------------!
! Carbohydrate concentration in cell lumen (mg ml-1).
!----------------------------------------------------------------------!
REAL(DP), DIMENSION (nfi, ncells_max) :: SUC
!----------------------------------------------------------------------!
REAL(DP) :: rlat              ! Latitude of site               (radians)
REAL(DP) :: pz_min            ! Minimum width of proliferation zone (�m)
REAL(DP) :: pz                ! Width of proliferation zone         (�m)
REAL(DP) :: ez                ! Width of enlargement zone           (�m)
REAL(DP) :: tz                ! Width of wall thickening zone       (�m)
REAL(DP) :: T                 ! Mean daily temperature               (K)
REAL(DP) :: TC                ! Mean daily temperature              (oC)
REAL(DP) :: dd                ! Degree-days                         (oC)
!----------------------------------------------------------------------!
INTEGER, DIMENSION (nfi) :: ncells       ! No. cells/file            (n)
INTEGER :: ncells_alive                  ! No. cells in cambial zone (n)
!----------------------------------------------------------------------!
INTEGER :: kyr                                  ! Year              (CE)
INTEGER :: kday                                 ! Day of year      (day)
INTEGER :: cd                                   ! Chilling days   (days)
!----------------------------------------------------------------------!
LOGICAL :: dorm                                 ! Domancy state      (-)
!----------------------------------------------------------------------!
END MODULE STATE
!======================================================================!
