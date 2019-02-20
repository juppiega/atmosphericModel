!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
! Time module
!
! - Time variables
!
! - Time functions
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

MODULE Time_Mod

!-----------------------------------------------------------------------------------------
! Load modules
!-----------------------------------------------------------------------------------------
USE Parameters_Mod

!-----------------------------------------------------------------------------------------
! Variable declaration
!-----------------------------------------------------------------------------------------
IMPLICIT NONE

PUBLIC

INTEGER, PARAMETER :: one_hour = 60*60  ! [s], one hour in seconds

REAL(dp) :: time  ! [s], current time
REAL(dp) :: time_start, time_end  ! [s], start and end time

REAL(dp) :: dt  ! [s], time step for main loop, usually is equal to meteorology time step
REAL(dp) :: dt_chem  ! [s], time step for chemistry calculation
REAL(dp) :: dt_micro  ! [s], time step for aerosol processes
REAL(dp) :: dt_output  ! [s], time step for output

REAL(dp) :: time_start_chemistry  ! [s], time to start calculating chemistry
REAL(dp) :: time_start_aerosol ! [s], time to start calculating aerosol
REAL(dp) :: time_start_output ! [s], time to start calculating aerosol

INTEGER :: daynumber_start  ! [day], start Julian day
INTEGER :: daynumber  ! [day], current Julian day

INTEGER :: counter  ! [-], counter of time steps


CONTAINS


!-----------------------------------------------------------------------------------------
! Interface
!
SUBROUTINE Time_Init()
!
! Declaration
!

!
! Description
!
! Set initial values for time_start, time_end, time steps, day number, counter
!-----------------------------------------------------------------------------------------
  !
  ! Simulation time period
  !
  time_start = 0
  time_end = 3*24*3600 + 16*60!4.5*86400 + 10 + time_start

  !
  ! Time steps
  !
  dt = 0.1
  dt_chem = 60
  dt_micro = 10
  dt_output = 900

  !
  ! Get the Julian date of Aug. 10, 2011
  !
  daynumber_start = 31+28+31+30+31+30+31+10

  !
  ! Start to run chemistry module after 1 day to save computation time
  !
  if (box) then
    time_start_chemistry = 0*24*one_hour + time_start
  else
    time_start_chemistry = time_end
  end if
  time_start_aerosol = time_end

  time_start_output = dt_output + time_start

  !
  ! Current time and date
  !
  time = time_start
  daynumber = daynumber_start

  !
  ! Count how many time steps the code runs
  !
  counter = 0
END SUBROUTINE Time_Init

END MODULE Time_Mod
