module glow_mod

use phys_consts, only : wp
use cglow,only: cglow_init   ! subroutine to allocate use-associated variables
use cglow,only: jmax,nbins,lmax,nmaj,nei,nex,nw,nc,nst
use cglow,only: idate,ut,glat,glong,f107a,f107,f107p,ap,ef,ec
use cglow,only: iscale,jlocal,kchem,xuvfac
use cglow,only: sza,dip,efrac,ierr
use cglow,only: zz,zo,zn2,zo2,zns,znd,zno,ztn,ze,zti,zte
use cglow,only: ener,del,phitop,wave1,wave2,sflux,pespec,sespec,uflx,dflx,sion
use cglow,only: photoi,photod,phono,aglw,ecalc,zxden,zeta,zceta,eheat,vcb
use cglow,only: data_dir
implicit none

logical :: first_call = .true.

contains

  subroutine glow_run(W0,PhiWmWm2,date_doy,UTsec,xf107,xf107a,xlat,xlon,alt,nn,Tn,ns,Ts,ionrate,eheating,iver)
  
  ! This software is part of the GLOW model.  Use is governed by the Open Source
  ! Academic Research License Agreement contained in the file glowlicense.txt.
  ! For more information see the file glow.txt.
  
  ! Stan Solomon and Ben Foster, 1/15
  ! Stan Solomon, 12/15, 1/16
  ! Stan Solomon, 3/16, MPI parallel version
  ! Guy Grubbs II, 4/17, modified for GEMINI integration
  
  ! Main multi-processor driver for the GLOW model.
  ! Uses TIE-GCM history files or MSIS/IRI for input.
  ! Requires MPI and netCDF libraries
  
  ! For definitions of use-associated variables, see subroutine GLOW and module CGLOW.
  ! For definitions of TGCM input variables see module READTGCM
  ! For definitions of output arrays see module OUTPUT
  
  ! Other definitions:
  ! f107p   Solar 10.7 cm flux for previous day
  ! ap      Ap index of geomagnetic activity
  ! z       altitude array, km
  
  ! Array dimensions:
  ! jmax    number of altitude levels
  ! nbins   number of energetic electron energy bins
  ! lmax    number of wavelength intervals for solar flux
  ! nmaj    number of major` species
  ! nst     number of states produced by photoionization/dissociation
  ! nei     number of states produced by electron impact
  ! nex     number of ionized/excited species
  ! nw      number of airglow emission wavelengths
  ! nc      number of component production terms for each emission
    
  
    real(wp), dimension(:), intent(in) :: W0,PhiWmWm2,alt,Tn
    real(wp), dimension(:,:), intent(in) :: nn,ns,Ts
    real(wp), dimension(:,:), intent(out) :: ionrate
    real(wp), dimension(:), intent(out) :: eheating, iver
    real(wp), intent(in) :: UTsec, xlat, xlon, xf107, xf107a
    integer, intent(in) :: date_doy

    real :: stl   
    real,allocatable :: z(:),zun(:),zvn(:)             ! glow height coordinate in km (jmax)

    real, dimension(nbins) :: phitoptmp = 0.0
    integer :: j
    character(len=1024) :: iri90_dir

  !
  ! Variables that need to be set for GLOW not sent from GEMINI
  !
    kchem = 4.
    jlocal = 0.
    iscale=1
    xuvfac=3.
  !
  ! Set data directories:
  !
    data_dir    ='ionization/glow/data/'
    iri90_dir   ='ionization/glow/data/iri90/'
  !
  ! Allocate arrays in other modules (formerly in common blocks):
  !
    if(first_call) then
      first_call = .false.
      jmax=102
      call cglow_init
    end if
    allocate(z(jmax))
    allocate(zun(jmax))
    allocate(zvn(jmax))
  !
  ! Set electron energy grid:
  !
    call egrid (ener, del, nbins)
  !
  ! Set variables given from GEMINI
  !
    glat = real(xlat,4)
    glong = real(xlon,4)
    idate = date_doy
    ut = real(UTsec,4)
    f107 = real(xf107,4)
    f107p = real(xf107,4)
    f107a = real(xf107a,4)
    ap = 31.5
  !
  ! Calculate local solar time:
  !
    stl = ut/3600. + glong/15.
    if (stl < 0.) stl = stl + 24.
    if (stl >= 24.) stl = stl - 24.
  !
  ! Call MZGRID to use MSIS/NOEM/IRI inputs on default altitude grid:
  !
    call mzgrid (jmax,nex,idate,ut,glat,glong,stl,f107a,f107,f107p,ap,iri90_dir, &
                 z,zo,zo2,zn2,zns,znd,zno,ztn,zun,zvn,ze,zti,zte,zxden)
    write(*,*) jmax,nex,idate,ut,glat,glong,stl,f107a,f107,f107p,ap,iri90_dir
  !
  ! Set Maxwellian distribution into phitop array
  !
    phitop=0.0; phitoptmp=0.0
    do j = 1, size(PhiWmWm2,1)
      write(*,*) real(PhiWmWm2(j),4), real(W0(j),4)
      call maxt(real(PhiWmWm2(j),4),real(W0(j),4),ener,del,nbins,0,0,0,phitoptmp)
      phitop=phitop+phitoptmp
    end do
  !
  ! Fill altitude array, converting from km to cm 
  !
    zz(:)  = z(:)*1.e5
  !
  ! Call GLOW to calculate ionized and excited species, and airglow emission rates:
  !
    call glow
    
    ionrate(:,1) = real(0.0,wp) !(real(SION(1,:),wp)+real(SION(2,:),wp)*0.3d0)*1.0d6 !O+
    ionrate(:,4) = real(0.0,wp) !(real(SION(2,:),wp)*0.7d0)*1.0d6 !O2+
    ionrate(:,3) = real(0.0,wp) !(real(SION(3,:),wp)*0.84d0)*1.0d6 !N2+
    ionrate(:,5) = real(0.0,wp) !(real(SION(3,:),wp)*0.16d0)*1.0d6 !N+
    ionrate(:,2) = real(0.0,wp) !real(0d0,wp)*1.0d6 !NO+
    ionrate(:,6) = real(0.0,wp) !real(0d0,wp)*1.0d6 !H+
    eheating = real(0.0,wp) !real(eheat,wp)*1.0d6
    iver = real(vcb,wp)

  ! ZXDEN   array of excited and and/or ionized state densities at each altitude:
  !           O+(2P), O+(2D), O+(4S), N+, N2+, O2+, NO+, N2(A), N(2P),
  !           N(2D), O(1S), O(1D); cm-3
  !ions (ns): 1=O+, 2=NO+, 3=N2+, 4=O2+, 5=N+, 6=H+
  !neutrals (nn): O,N2,O2,H,N,NO
  
  end subroutine glow_run

end module glow_mod
