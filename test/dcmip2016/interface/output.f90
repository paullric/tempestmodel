!=======================================================================
!
!   Interface output visualization
!
!
!
!=======================================================================

    SUBROUTINE output_baroclinic_viz(lon, lat, p, z, u, v, w, t, phis, ps, rho, q,zcoords )
    IMPLICIT NONE

    REAL(8), INTENT(IN)  :: &
    lon,        & ! Longitude (radians)
    lat,        & ! Latitude (radians)
    z,          & ! Altitude (m)
    p             ! Pressure (Pa)

    REAL(8), INTENT(IN) :: &
    u,          & ! Zonal wind (m s^-1)
    v,          & ! Meridional wind (m s^-1)
    w,          & ! Vertical Velocity (m s^-1)
    t,          & ! Temperature (K)
    phis,       & ! Surface Geopotential (m^2 s^-2)
    ps,         & ! Surface Pressure (Pa)
    rho,        & ! density (kg m^-3)
    q             ! water vapor mixing ratio (kg/kg)

    INTEGER, INTENT(IN) :: zcoords     ! 1 if z coordinates are specified
                                        ! 0 if p coordinates are specified





    !print *,lon,lat,z
    print *,"U Velocity",u
    print *,"V Velocity",v
    print *,"W Velocity",w
    print *,"Temperature",t
    print *,"Surf. Phi",phis
    print *,"Surf. Pres",ps
    print *,"Rho",rho
    print *,"Q",q
    
    if (zcoords .eq. 1) then
    print *,"Pressure",p
    else

    print *,"z, height in m",z
    end if

    END SUBROUTINE output_baroclinic_viz

    SUBROUTINE output_term_ini_value_viz(lat,lon,cl,cl2 )
    IMPLICIT NONE

    REAL(8), INTENT(IN)  :: &
    lon,        & ! Longitude (radians)
    lat,        & ! Latitude (radians)
    cl, cl2   ! molar mixing ratio of cl and cl2

    print *, " "
    print *, " Initial Values"
    print *, " Latitude: ",lat," Longitude: ",lon
    print *, " Cl: ",cl," Cl2: ",cl2

    END SUBROUTINE output_term_ini_value_viz

    SUBROUTINE output_term__tendency_viz(cl_f,cl_2f)
    IMPLICIT NONE

    REAL(8), INTENT(IN)   :: cl_f, cl_2f  ! time rate of change of cl and cl2

    print *, " "
    print *, " Tendency should be zero for the initial values"
    print *, " Cl rate ",cl_f," Cl2 rate: ",cl_2f
    print *, " "

    END SUBROUTINE output_term__tendency_viz




