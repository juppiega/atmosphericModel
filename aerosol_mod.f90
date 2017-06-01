
module Aerosol_mod
    use parameters_mod
    use time_mod
    IMPLICIT NONE
    private
    public test_aerosol

    !! ====================== Definition of variables =====================================================================
      ! so that numbers will be in 64bit floating point
    ! http://en.wikipedia.org/wiki/Double_precision_floating-point_format

    REAL(dp), PARAMETER :: ka = 0.4D0              ! von Karman constant, dimensionless
    INTEGER, PARAMETER ::  nr_bins = 100           ! Number of particle size bins
    INTEGER, PARAMETER ::  nr_cond = 2             ! Number of condensable vapours
    REAL(DP), PARAMETER :: Rg=8.3145D0             ! Universal gas constant J mol^-1 K^-1

    REAL(DP), DIMENSION(nr_bins) :: diameter  ,&    ! Diameter of each size bin
        particle_mass                        ,&    ! mass of one particle in each size bin
        size_distribution                        ,&    ! number concentration in each size bin
        particle_volume                      ,&    ! volume concentration in each size bin
        coag_loss                            ,&    ! coagulation loss rate of particles in each size bin
        v_dep                                      ! Dry deposition velocity of particles

    REAL(DP), DIMENSION(nr_cond) :: molecular_mass   ,&   ! molecular mass of the condensing vapours [kg/#]
        molecular_volume                            ,&   ! Molecule volume of condensable vapours [m^3]
        molecular_dia                               ,&   ! Molecule diameter of condensable vapours [m]
        molar_mass                                  ,&   ! Molar mass of condensable vapours [kg/m^3]
        cond_vapour                                      ! Concentration of condensable vapours [molec/m^3]

    REAL(DP), DIMENSION(nr_cond) :: Cond_sink             ! Condensation sink of vapours [s^-1]

    REAL(DP) ::         PN, PM, PV                  ! Total particle number and mass concentration [cm^-3]

    REAL(DP) :: vd_SO2,vd_O3, vd_HNO3           ! Dry deposition velocity of SO2, O3 & HNO3

    REAL(DP)    particle_density              ,&    ! [Kg]
        nucleation_coef               ,&    ! Nucleation coefficient
        mass_accomm                         ! mass accomodation coefficient

    REAL(dp) :: nucleation_rate                     ! #/(m3*s)

    real(dp) :: &
        temperature                   ,&    ! Temperature [K]
        Richards_nr10m                ,&    ! Richards number at the reference altitude 10 m
        wind_speed10m                 ,&    ! Wind speed at the reference altitude 10 m
        pressure                      ,&    ! Pressure [Pa]
        DSWF                          ,&    ! Downward Shortwave Radiation Flux (W/m^2)
        Mixing_height                       ! Boundary layer mixing height [m]

CONTAINS

    subroutine compute_PN()
        implicit none
        PN = SUM(size_distribution)*1D-6                 ! Total particle number concentration (cm^-3)
    end subroutine

    subroutine compute_PM()
        implicit none
        PM=SUM(size_distribution*particle_mass)*1D9    ! Total particle mass concentration (ug/m^3)
    end subroutine

    subroutine compute_PV()
        implicit none
        PV=SUM(size_distribution*particle_volume)*1D12    ! Total particle mass concentration (ug/m^3)
    end subroutine


    subroutine output_aerosols()
        implicit none
        logical :: first_call = .true.
        integer, save :: PM_unit, PN_unit, PV_unit, time_unit

        if (first_call) then
            OPEN(FILE = TRIM(ADJUSTL(outdir))//adjustl('/')//adjustl(trim('PM'))//adjustl('.dat')  , &
                STATUS = 'REPLACE', ACTION = 'WRITE', newunit=PM_unit)
            OPEN(FILE = TRIM(ADJUSTL(outdir))//adjustl('/')//adjustl(trim('PN'))//adjustl('.dat')  , &
                STATUS = 'REPLACE', ACTION = 'WRITE', newunit=PN_unit)
            OPEN(FILE = TRIM(ADJUSTL(outdir))//adjustl('/')//adjustl(trim('PV'))//adjustl('.dat')  , &
                STATUS = 'REPLACE', ACTION = 'WRITE', newunit=PV_unit)
            OPEN(FILE = TRIM(ADJUSTL(outdir))//adjustl('/')//adjustl(trim('time'))//adjustl('.dat')  , &
                STATUS = 'REPLACE', ACTION = 'WRITE', newunit=time_unit)
            first_call = .false.
!            print *, 'First call to aer. output'
        end if

        WRITE(PM_unit, *) PM
        WRITE(PN_unit, *) PN
        write(PV_unit, *) PV
        write(time_unit,*) time/86400

    end subroutine

    subroutine test_aerosol()
        implicit none
        real(kind=8) :: last_output, simu_hours
            !! ======================= Programe starts ===========================================================================

        ! Assign values to parameters and initialize the simulation

        CALL Aerosol_init(diameter, particle_mass, particle_volume, size_distribution, &
            particle_density, nucleation_coef, molecular_mass, molar_mass, &
            molecular_volume, molecular_dia, mass_accomm)

        simu_hours = 24D0 ! hours
        dt_aero = 1D1 ! s
        time_start = 0D0
        time_end = time_start + simu_hours*3600.
        time = time_start
        last_output = -dt_output

        DO WHILE (time .lt. time_end) ! Main program time step loop

            ! Meteorological parameters:
            temperature = 300D0  ! K
            pressure = 1D5       ! Pa

            ! In the 1D-blm you will instead use:
            ! dh = h(2:nz) - h(1:nz-1) ! The height of each model layer
            ! Cp = 1D3  ! Air specific heat at constant pressure, J/(kg��K)
            ! temperature = theta - (g/Cp)*h ! Temperature in K
            ! pressure(1)=1D5 ! Surface pressure (Pa)
            ! DO i=2,nz
            !    pressure(i) = pressure(i-1)*exp(-Mair*g/&
            !   (Rg*(temperature(i-1)+temperature(i))/2D0)*dh(i-1)) ! Barometric law
            ! END DO

            Richards_nr10m = 0D0 ! Richards number at 10 m above ground (0 for neutral atmosphere)
            wind_speed10m = 2D0  ! Wind speed at the reference altitude 10 m [m/s]
            Mixing_height = 1D3  ! Assumed mixing height of the box model (m)
            DSWF = 6D2 * sin(pi/86400 * time) ! Downward Shortwave radiation flux (diurnal cycle) (W/m^2)
            ! In the 1D-blm model we use:
            ! DSWF = 6D2 * exp_coszen ! Approximate downward shortwave rations flux (W/m^2)

            cond_vapour(1) = 1D13*sin(pi/86400 * time) ! H2SO4 molec / m^3
            cond_vapour(2) = 1D13                      ! ELVOC molec / m^3

            !!! Calculate particle dry deposition velocity and particle losses due to dry deposition here !!!

            !!! Calculate new 2nm particle formation (nucleation) here !!!
            call Nucleation(cond_vapour(1), size_distribution(1))

            !!! Calculate coagulation losses here !!!
            call Coagulation(dt_aero, size_distribution, diameter, &
            temperature,pressure,particle_mass)

            !!! Calculate condensation particle growth here !!!
            call Condensation(dt_aero, temperature, pressure, molecular_mass, &
            molecular_volume, molar_mass, molecular_dia, particle_mass, particle_volume, &
            size_distribution, diameter, cond_vapour)

            ! Compute statistics.

            call compute_PM()
            call compute_PN()
            call compute_PV()

            if (time - last_output >= dt_output) then
                call output_aerosols()
                last_output = time
                write(*,*) 'time', time/3600.
            end if

            time = time + dt_aero


        ENDDO

    end subroutine

    SUBROUTINE Aerosol_init(diameter, particle_mass, particle_volume, particle_conc, &
        particle_density, nucleation_coef, molecular_mass, molar_mass, &
        molecular_volume, molecular_dia, mass_accomm)

        !! ====================== Definition of variables =====================================================================

        REAL(DP), DIMENSION(nr_bins), INTENT(OUT) :: diameter       , &    ! diamter of each size bin
            particle_mass  , &    ! mass of one particle
            particle_volume  , &  ! volume of one particle
            particle_conc         ! number concentration

        REAL(DP), DIMENSION(nr_cond), INTENT(OUT) :: molecular_mass ,& ! molecular mass of the condensing vapours [kg/#]
            molecular_volume ,&     ! [m3]
            molecular_dia, &        ! [m]
            molar_mass              ! molar mass of the condensing vapours [kg/mol]

        REAL(DP), INTENT(OUT) :: nucleation_coef, mass_accomm

        REAL(DP), DIMENSION(nr_cond) :: density             ! Bulk density of condensing vapours [kg/m^3]
        real(dp), intent(inout) :: particle_density

        INTEGER :: i
        mass_accomm = 1D0   ! Mass accommodation coefficient

        nucleation_coef = 1D-20

        ! Particle diameters between 2D-9 and 2.5D-6 m:
        diameter(1)=2D-9
        DO i=2,nr_bins
            diameter(i)=diameter(i-1)*(2.5D-6/diameter(1))**(1D0/(nr_bins-1))
        END DO

        particle_conc = 1D0 ! Assume an initial particle number concentration of 1 m^-3
        where((abs(diameter-2D-7)-MINVAL(abs(diameter-2D-7)))<1D-20)  particle_conc=2D8 ! add 200 cm^-3 200 nm sized accumulation mode particles

        particle_density = 1.4D3                                        ! Assumed fixed particle density [kg/m^3]
        particle_volume = 1D0/6D0 * pi * diameter**3                      ! Single particle volume (m^3)
        particle_mass=  1D0/6D0 * pi * diameter**3 * particle_density     ! [kg]

        density = (/1.84D3, 1.4D3/)                                     ! density of sulphuric acid and SOA
        molar_mass = (/0.098D0, 0.3D0/)                                 ! H2SO4 and ELVOC
        molecular_mass = molar_mass / Na                                ! molecular mass [kg]
        molecular_volume = molecular_mass / density                     ! molecular volume [m^3]
        molecular_dia = (6D0 * molecular_volume / pi )**(1D0/3D0)       ! molecular diameter [m]

    END SUBROUTINE Aerosol_init

    SUBROUTINE Nucleation(H2SO4_conc, smallest_conc)
    real(kind = 8), intent(in) :: H2SO4_conc ! [molec/m^3]
    real(kind = 8), intent(inout) :: smallest_conc
    real(kind = 8) :: nucleation_rate
    real(kind = 8), parameter :: K = 1D-20 ! Nucleation coeff [m^3/molec/s]

    nucleation_rate = K * H2SO4_conc**2
    smallest_conc = smallest_conc + nucleation_rate*dt_aero

    END SUBROUTINE Nucleation

    SUBROUTINE Condensation(dt_aero, temperature, pressure, molecular_mass, &
        molecular_volume, molar_mass, molecular_dia, particle_mass, particle_volume, &
        particle_conc, diameter, cond_vapour) ! Add more variables if you need it

        REAL(DP), DIMENSION(nr_bins), INTENT(IN) :: diameter, particle_mass
        REAL(DP), DIMENSION(nr_cond), INTENT(IN) :: molecular_mass, molecular_dia, &
            molecular_volume, molar_mass
        REAL(DP), INTENT(IN) :: dt_aero, temperature, pressure

        REAL(DP), DIMENSION(nr_bins), INTENT(INOUT) :: particle_conc

        REAL(DP), DIMENSION(2), INTENT(IN) :: cond_vapour  ! condensing vapour concentrations, which is H2SO4 and organics (ELVOC) [#/m^3]

        REAL(DP), DIMENSION(nr_bins), INTENT(IN)   :: particle_volume

        REAL(DP), DIMENSION(nr_bins)   ::  slip_correction, diffusivity_p, speed_p, &
            particle_conc_new, particle_volume_new

        REAL(DP), DIMENSION(nr_cond)   ::  diffusivity_gas, speed_gas

        REAL(DP) :: dyn_visc, l_gas, dens_air, x1, x2

        real(dp),dimension(nr_bins) :: fuchs_sutugin, Kn, lambda, CR

        INTEGER :: j

        ! Add more variabels as you need it...

        dyn_visc = 1.8D-5*(temperature/298D0)**0.85D0  ! dynamic viscosity of air
        dens_air=Mair*pressure/(Rg*temperature)        ! Air density
        l_gas=2D0*dyn_visc/(pressure*SQRT(8D0*Mair/(pi*Rg*temperature))) ! Gas mean free path in air (m)

        slip_correction = 1D0+(2D0*l_gas/(diameter))*&
            (1.257D0+0.4D0*exp(-1.1D0/(2D0*l_gas/diameter))) ! Cunninghams slip correction factor (Seinfeld and Pandis eq 9.34)

        diffusivity_p = slip_correction*kb*temperature/(3D0*pi*dyn_visc*diameter)   ! Diffusivity for the different particle sizes m^2/s
        speed_p = SQRT(8D0*kb*temperature/(pi*particle_mass))                     ! speed of particles (m/s)

        diffusivity_gas=5D0/(16D0*Na*molecular_dia**2D0*dens_air)*&
            SQRT(Rg*temperature*Mair/(2D0*pi)*(molar_mass+Mair)/molar_mass)            ! Diffusivity of condensable vapours (m^2 s^-1)

        ! Thermal velocity of vapour molecule
        speed_gas=SQRT(8D0*kb*temperature/(pi*molecular_mass)) ! speed of H2SO4 molecule



        particle_volume_new = particle_volume
        do j = 1, size(cond_vapour) ! Loop over molecules
        ! Calculate the Fuchs-Sutugin correction factor:
        lambda = 3*(diffusivity_p + diffusivity_gas(j)) / sqrt(speed_p**2 + speed_gas(j)**2)
        Kn = 2*lambda / (diameter + molecular_dia(j))
        fuchs_sutugin = 0.75*(1+Kn) / (Kn**2 + Kn + 0.283*Kn + 0.75)

        ! Calculate the Collision rate (CR [m^3/2]) between gas molecules (H2SO4 and ELVOC) and the particles:
        CR = 2*pi * (diameter + molecular_dia(j)) * (diffusivity_gas(j) + diffusivity_p) * fuchs_sutugin

        ! Calculate the new single particle volume after condensation (particle_volume_new):
        particle_volume_new = particle_volume_new + CR*cond_vapour(j)*molecular_volume(j) *dt_aero
        end do

        particle_conc_new = 0.0
        particle_conc_new(nr_bins) = particle_conc(nr_bins)


        DO j = 1,nr_bins-1
            x1 = (particle_volume(j+1) - particle_volume_new(j)) / (particle_volume(j+1) - particle_volume(j))
            x2 = 1 - x1
            particle_conc_new(j) = particle_conc_new(j) + x1*particle_conc(j)
            particle_conc_new(j+1) = particle_conc_new(j+1) + x2*particle_conc(j)
        END DO
        ! Update the particle concentration in the particle_conc vector:
        particle_conc=particle_conc_new

    END SUBROUTINE Condensation

    SUBROUTINE Coagulation(dt_aero, particle_conc, diameter, &
        temperature,pressure,particle_mass) ! Add more variables if you need it

        REAL(DP), DIMENSION(nr_bins), INTENT(IN) :: diameter
        REAL(DP), DIMENSION(nr_bins), INTENT(INOUT) :: particle_conc
        REAL(DP), INTENT(IN) :: dt_aero
        REAL(DP), DIMENSION(nr_bins), INTENT(IN) :: particle_mass       ! mass of one particle
        REAL(DP), INTENT(IN) :: temperature, pressure

        REAL(DP), DIMENSION(nr_bins,nr_bins) :: coagulation_coef        ! coagulation coefficients [m^3/s]

        REAL(DP), DIMENSION(nr_bins) :: slip_correction, diffusivity, dist, speed_p, &
            Beta_Fuchs, free_path_p
        real(dp), dimension(nr_bins) :: loss_self, loss_larger

        REAL(DP) ::       dyn_visc, &                                   ! dynamic viscosity, kg/(m*s)
            l_gas                                         ! Gas mean free path in air
        INTEGER  :: i

        ! The Coagulation coefficient is calculated according to formula 13.56 in Seinfield and Pandis (2006), Page 603

        dyn_visc = 1.8D-5*(temperature/298.)**0.85                                              ! Dynamic viscosity of air

        l_gas=2D0*dyn_visc/(pressure*SQRT(8D0*Mair/(pi*Rg*temperature)))                        ! Gas mean free path in air (m)

        slip_correction = 1D0+(2D0*l_gas/(diameter))*&
            (1.257D0+0.4D0*exp(-1.1D0/(2D0*l_gas/diameter)))                                        ! Cunninghams slip correction factor (Seinfeld and Pandis eq 9.34)

        diffusivity = slip_correction*kb*temperature/(3D0*pi*dyn_visc*diameter)                 ! Diffusivity for the different particle sizes m^2/s

        speed_p = SQRT(8D0*kb*temperature/(pi*particle_mass))                                   ! Speed of particles (m/s)

        free_path_p = 8D0*diffusivity/(pi*speed_p)                                              ! Particle mean free path (m)

        dist = (1D0/(3D0*diameter*free_path_p))*((diameter+free_path_p)**3D0 &
            -(diameter**2D0+free_path_p**2D0)**(3D0/2D0))-diameter                    ! mean distance from the center of a sphere reached by particles leaving the sphere's surface (m)

        DO i = 1,nr_bins
            Beta_Fuchs = 1D0/((diameter+diameter(i))/(diameter+diameter(i)+&
                2D0*(dist**2D0+dist(i)**2D0)**0.5D0)+8D0*(diffusivity+diffusivity(i))/&
                (((speed_p**2D0+speed_p(i)**2D0)**0.5D0)*(diameter+diameter(i))))                    ! Fuchs correction factor from Seinfeld and Pandis, 2006, p. 600

            coagulation_coef(i,:) = 2D0*pi*Beta_Fuchs*(diameter*diffusivity(i)+&
                diameter*diffusivity+diameter(i)*diffusivity+diameter(i)*diffusivity(i))             ! coagulation rates between two particles of all size combinations  (m^3/s)
        END DO

    ! Write equations that considers how the particle number concentration in each size bin
    !(particle_conc) is influenced by the coagulation sink (loss of smaller particles when
    ! they collide with larger ones)


    ! You can first calculate the loss (loss1) do to self-coagulation between particles in the same size bin
    ! and then calculate the loss (loss2) due to coagulation with larger particles
    ! Then add the two loss terms together loss = loss1 + loss2

    do i = 1, nr_bins-1
        loss_self(i) = -0.5*coagulation_coef(i,i)*particle_conc(i)**2
        loss_larger(i) = -particle_conc(i)*sum(coagulation_coef(i, i+1 : nr_bins) * particle_conc(i+1 : nr_bins))
    end do
    loss_larger(nr_bins) = 0
    loss_self(nr_bins) = -0.5*coagulation_coef(nr_bins,nr_bins)*particle_conc(nr_bins)**2

    particle_conc = particle_conc + (loss_self + loss_larger) * dt_aero


    END SUBROUTINE Coagulation

    SUBROUTINE dry_dep_velocity(diameter,particle_density,temperature,pressure,DSWF, &
        Richards_nr10m,wind_speed10m) ! Add more variables if you need it

        REAL(DP), DIMENSION(nr_bins), INTENT(IN) :: diameter

        REAL(DP), INTENT(IN) :: temperature, pressure, Richards_nr10m, DSWF, &
            wind_speed10m, particle_density

        REAL(DP) :: z0, r_coll, a_landuse, j_landuse, v_kinematic,dyn_visc,Pr,beta,&
            gam,zr,u_friction,dens_air,L_Ob


        dens_air = Mair*pressure/(Rg*temperature)    ! Air density (kg/m^3)
        dyn_visc = 1.8D-5*(temperature/298.)**0.85   ! dynamic viscosity of air (kg/(m*s))
        v_kinematic = dyn_visc/dens_air              ! kinematic viscosity of air (m^2/s)

        zr=10D0                                      ! Reference height [m]
        L_Ob=zr/Richards_nr10m                       ! Monin-Obukhov length scale
        u_friction=ka*wind_speed10m/(log(zr/z0))     ! Friction velocity (Eq. 16.67 from Seinfeld and Pandis, 2006)

        ! Land use category paramaters from Seinfeld and Pandis, 2006 Table 19.2:
        z0 = 0.9D0 ! Roughness length evergreen, needleleaf trees (m)
        r_coll = 2D-3 ! radius of collector evergreen, needleleaf trees

        ! coefficients based on land use categories (evergreen, needleleaf trees)
        a_landuse = 1D0
        j_landuse = 0.56D0

        Pr = 0.95D0        ! Turbulent Prandtl number (when ka = 0.4 (Hogstrom, 1988))
        beta = 7.8D0       ! When ka = 0.4 (Hogstrom, 1988)
        gam = 11.6D0       ! When ka = 0.4 (Hogstrom, 1988)

       ! Calculate the particle sedimentation velocity:

       ! Calculation of aerodynamic resistance for particles for:
       ! stable boundary layer (Ri>1D-6)
       ! neutral boundary layer (abs(Ri)<1D-6
       ! unstable boundary layer Ri<-1D-6


       ! Calculate the quasi-laminar resistance (rb) for particles:

       ! Calculate the dry deposition velocity for particles:


    END SUBROUTINE dry_dep_velocity

END module Aerosol_mod
