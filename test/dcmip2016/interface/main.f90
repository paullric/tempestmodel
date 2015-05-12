


  SUBROUTINE tc_tropical(lon, lat, p, z, u,v,t,phis,ps,rho,q)
  USE tropical_cyclone_test

    REAL(8), INTENT(INOUT)  :: &
                z,          & ! Altitude (m)
                p             ! Pressure (Pa)

    REAL(8), INTENT(IN)  :: &
                lon,        & ! Longitude (radians)
                lat         ! Latitude (radians)
                


    REAL(8), INTENT(OUT)  :: &

                
                u,          & ! Zonal wind (m s^-1)
                v,          & ! Meridional wind (m s^-1)
                t,          & ! Temperature (K)
                phis,       & ! Surface Geopotential (m^2 s^-2)
                ps,         & ! Surface Pressure (Pa)
                rho,        & ! density (kg m^-3)
                q   		  ! Specific Humidity (kg/kg)


    REAL(8) ::X 
    Integer(4) :: deep, moist, pertt, zcoords

    deep=0
    moist=0
    pertt=0
    X=1.0d0
    zcoords=1

    CALL tc_initial_vortex(lon,lat,p,z,zcoords,u,v,t,phis,ps,rho,q)




  END SUBROUTINE tc_tropical

  SUBROUTINE tc_baroclinic(lon, lat, p, z,u,v,w,t,phis,ps,rho,q)
  USE test_baroclinic


    REAL(8), INTENT(IN)  :: &
                lon,        & ! Longitude (radians)
                lat           ! Latitude (radians)
    REAL(8), INTENT(INOUT)  :: &
                z,          & ! Altitude (m)
                p             ! Pressure (Pa)
    REAL(8), INTENT(OUT)  :: &


                u,          & ! Zonal wind (m s^-1)
                v,          & ! Meridional wind (m s^-1)
                w,          & ! Vertical Velocity (m s^-1)
                t,          & ! Temperature (K)
                phis,       & ! Surface Geopotential (m^2 s^-2)
                ps,         & ! Surface Pressure (Pa)
                rho,        & ! density (kg m^-3)
                q   		  ! Specific Humidity (kg/kg)

    REAL(8) ::  X
    Integer(4) :: deep, moist, pertt, zcoords

    deep=0
    moist=0
    pertt=0
    X=1.0d0
    zcoords=1

    CALL baroclinic_instability_alt(deep, moist, pertt, X, lon, lat, p, z, zcoords, u, v, w, t, phis, ps, rho, q)



  END SUBROUTINE tc_baroclinic



 !SUBROUTINE tc_terminator_test(lat,lon,dt,cl, cl2, cl_f, cl_2f)
 !USE Terminator, only: initial_value_Terminator, tendency_Terminator

  !  REAL(8), INTENT(IN)  :: &
   !             lon,        & ! Longitude (radians)
    !            lat,        & ! Latitude (radians)
        

   ! REAL(8), INTENT(OUT)  :: &
    !            cl,        &
     !           cl2,       &
      !          cl_f,      &
       !         cl_2f

    
   ! REAL(8) ::  dt=1800 ! seconds
    
   ! CALL initial_value_Terminator(lat, lon, cl, cl2)
   ! CALL tendency_Terminator(lat, lon, cl, cl2, dt, cl_f, cl_2f)

 ! END SUBROUTINE tc_terminator_test



