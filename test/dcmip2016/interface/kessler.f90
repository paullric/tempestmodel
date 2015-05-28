


!-----------------------------------------------------------------------
subroutine KESSLER( t, qv, qc, qr, rho, pk, dt, z, nz, rainnc )
!-----------------------------------------------------------------------
!
!     The KESSLER subroutine implements the Kessler (1969) microphysics
!     parameterization as described by Soong and Ogura (1973) and Klemp
!     and Wilhelmson (1978, KW). KESSLER is called at the end of each
!     time step and makes the final adjustments to the potential
!     temperature and moisture variables due to microphysical processes
!     occurring during that time step. KESSLER is called once for each
!     vertical column of grid cells. Increments are computed and added
!     into the respective variables. The Kessler scheme contains three
!     moisture categories: water vapor, cloud water (liquid water that
!     moves with the flow), and rain water (liquid water that falls
!     relative to the surrounding air). There  are no ice categories.
!     Variables in the column are ordered from the surface to the top.
! Input variables:
!     t      - potential temperature (K)
!     qv     - water vapor mixing ratio (gm/gm)
!     qc     - cloud water mixing ratio (gm/gm)
!     qr     - rain  water mixing ratio (gm/gm)
!     rho    - dry air density (not mean state as in KW) (kg/m^3)
!     pk     - Exner function  (not mean state as in KW) (p/p0)**(R/cp)
!     dt     - time step (s)
!     z      - heights of thermodynamic levels in the grid column (m)
!     nz     - number of thermodynamic levels in the column
!     rainnc - accumulated precip beneath the grid column (mm)
! Output variables:
!     Increments are added into t, qv, qc, qr, and rainnc which are
!     returned to the routine from which KESSLER was called. To obtain
!     the total precip qt, after calling the KESSLER routine, compute:
!       qt = sum over surface grid cells of (rainnc * cell area)  (kg)
!       [here, the conversion to kg uses (10^3 kg/m^3)*(10^-3 m/mm) = 1]
!
implicit none



!-----------------------------------------------------------------------
!     input/output params parameters at given location
!-----------------------------------------------------------------------

REAL(8), INTENT(INOUT) :: &
            t (nz),     & ! potential temperature (K)
            qv (nz),    & ! water vapor mixing ratio (gm/gm)
            qc (nz),    & ! cloud water mixing ratio (gm/gm)
            qr (nz),    & ! rain  water mixing ratio (gm/gm)
            rho (nz),   & ! dry air density (not mean state as in KW) (kg/m^3)
            rainnc        ! accumulated precip beneath the grid column (mm)


REAL(8), INTENT(IN) :: &
            z (nz),     & ! heights of thermodynamic levels in the grid column (m)
            pk(nz),         & ! Exner function  (not mean state as in KW) (p/p0)**(R/cp)
            dt         ! time step (s)

INTEGER, INTENT(IN) :: nz  	  ! number of thermodynamic levels in the column



!-----------------------------------------------------------------------


REAL, DIMENSION(nz) :: r, rhalf, velqr, sed, pc
real(8) :: f5, f2x, xk, ern, qrprod, prod, qvs, psl, rhoqr
integer :: k
!-----------------------------------------------------------------------
    f2x =17.27d0
    f5 = 237.3d0*f2x*2500000d0/1003.d0
    xk = .2875d0      !  r/cp
psl    = 1000.d0      !  mb
rhoqr  = 1000.d0    !  density of liquid water kg/m^3
do k=1,nz
r(k)     = 0.001d0*rho(k)
rhalf(k) = sqrt(rho(1)/rho(k))
pc(k)    = 3.8d0/(pk(k)**(1./xk)*psl)
!        liquid water terminal velocity (m/s) following KW eq. 2.15
velqr(k)  = 36.34d0*(qr(k)*r(k))**0.1364*rhalf(k)
end do
!     precipitation accumulated beneath the column
rainnc = rainnc + 1000.d0*rho(1)*qr(1)*velqr(1)*dt/rhoqr   ! mm rain
!    sedimentatiion term using upstream differencing
do k=1,nz-1
sed(k) = dt*(r(k+1)*qr(k+1)*velqr(k+1)-r(k)*qr(k)*velqr(k))/(r(k)*(z(k+1)-z(k)))
end do
sed(nz)  = -dt*qr(nz)*velqr(nz)/(.5*(z(nz)-z(nz-1)))
do k=1,nz
!        autoconversion and accretion rates following KW eq. 2.13a,b
qrprod = qc(k) - (qc(k)-dt*amax1(.001*(qc(k)-.001d0),0.))/(1.d0+dt*2.2d0*qr(k)**.875)
qc(k) = amax1(qc(k)-qrprod,0.)
qr(k) = amax1(qr(k)+qrprod+sed(k),0.)
!        saturation vapor mixing ratio (gm/gm) following KW eq. 2.11
qvs    = pc(k)*exp(f2x*(pk(k)*t(k)-273.d0)   &
/(pk(k)*t(k)- 36.d0))
prod    = (qv(k)-qvs)/(1.d0+qvs*f5/(pk(k)*t(k)-36.d0)**2)
!        evaporation rate following KW eq. 2.14a,b
ern     = amin1(dt*(((1.6d0+124.9d0*(r(k)*qr(k))**.2046)  &
*(r(k)*qr(k))**.525)/(2550000d0*pc(k)                &
/(3.8d0 *qvs)+540000d0))*(dim(qvs,qv(k))               &
/(r(k)*qvs)),amax1(-prod-qc(k),0.),qr(k))
!        saturation adjustment following KW eq. 3.10
t(k)= t(k) + 2500000d0/(1003.d0*pk(k))*(amax1( prod,-qc(k))-ern)
qv(k) = amax1(qv(k)-max(prod,-qc(k))+ern,0.)
qc(k) = qc(k)+max(prod,-qc(k))
qr(k) = qr(k)-ern
end do
end subroutine KESSLER
!-----------------------------------------------------------------------

