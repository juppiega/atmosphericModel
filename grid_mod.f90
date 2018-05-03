!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
! Grid module
!
! - Grid variables
!
! - Grid functions
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

MODULE grid_mod

!-----------------------------------------------------------------------------------------
! Load modules
!-----------------------------------------------------------------------------------------
USE parameters_mod

!-----------------------------------------------------------------------------------------
! Variable declaration
!-----------------------------------------------------------------------------------------
IMPLICIT NONE

PRIVATE
PUBLIC :: nz, z, zmid, updInd ! TODO: nz .> numLevels, z -> zLevels
integer :: i

!
! Number of height levels

!
integer, parameter :: topHeight = 3010
integer, parameter :: resolution = 10
INTEGER, PARAMETER :: nz = topHeight / resolution + 1
!integer, parameter :: nz = 50

integer, parameter :: updInd(nz-2) = (/ (i, i = 2, nz-1) /) ! UpdateIndices (excluding boundary points)

!
! Model height levels
!
!REAL(dp), PARAMETER, DIMENSION(nz) :: z = (/   0,   10,   20,   30,   40, &
!                                               50,   60,   70,   80,   90, &
!                                              100,  120,  140,  160,  180, &
!                                              200,  230,  260,  300,  350, &
!                                              400,  450,  500,  550,  600, &
!                                              650,  700,  800,  900, 1000, &
!                                             1100, 1200, 1300, 1400, 1500, &
!                                             1600, 1700, 1800, 1900, 2000, &
!                                             2100, 2200, 2300, 2400, 2500, &
!                                             2600, 2700, 2800, 2900, 3000 /)
real(dp), parameter, dimension(nz) :: z = (/ (i, i = 0, topHeight, resolution) /)
real(dp), parameter :: zmid(nz-1) = 0.5*(z(1:nz-1) + z(2:nz))

END MODULE grid_mod
