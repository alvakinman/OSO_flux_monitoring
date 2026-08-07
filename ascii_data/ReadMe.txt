J/A+A/XXX/XXX  Flux density data of 361 active galactic nuclei  (Kinman+, 2026)
================================================================================
A multi-band radio flux density catalog of ICRF3 sources using the Onsala
Twin Telescopes
       Kinman A., Haas R., Le Bail K., Yang J., Varenius E.
================================================================================
ADC_Keywords: Active gal. nuclei ; Radio continuum ;
Keywords: catalogs - galaxies: active - radio continuum: galaxies -
          instrumentation: interferometers - reference systems -

Abstract:
    Context: The VLBI Global Observing System (VGOS) is the next generation
system for geodetic and astrometric Very Long Baseline Interferometry (VLBI).
VGOS is broadband and uses fast-moving antennas that are capable of observing
more scans per minute than its predecessors. To optimize the observing time for
each source in geodetic schedules, a flux density catalog is needed for the 
sources that are observed at the VGOS frequencies.
    Aim: The aim of this work is to monitor the flux densities of geodetic 
sources in the VGOS bands. The obtained flux density time series can be used for
more effective scheduling of geodetic and astrometric VLBI experiments, as well 
as probing active galactic nuclei (AGN) physics. 
   Method: The Onsala Twin Telescopes have been used as a single baseline 
interferometer to measure flux densities of AGN that are part of the
International Celestial Reference System (ICRF3). The telescopes observed four 
frequency bands simultaneously, centered at 3.2, 5.5, 6.6 and 10.4~GHz. Both 
locally planned flux monitoring sessions and international geodetic experiments 
have been analysed. The data has been correlated using the Distributed FX (DiFX)
software correlator and calibrated using the Common Astronomy Software 
Applications (CASA). The possibility of predicting geodetic signal-to-noise 
ratios (S/N) using the measured flux densities was also tested.
   Results: Simultaneous light curves in up to four frequencies have been
obtained for 361 sources, of which 269 sources have at least five observation 
epochs. The majority of the sources vary significantly in flux density during
the measurement period. Most sources have a flat or inverted spectrum, with only
6 % having a steep spectrum. There is also evidence of spectral index variation
over time for a number of sources. Furthermore, the flux densities from this 
work were shown to more precisely predict geodetic signal-to-noise ratios
compared to the standard VGOS flux density catalog, in particular when 
considering the most variable sources.
   Conclusions: The majority of the radio sources observed by VGOS show 
significant variations in flux density. Such variation needs to be taken into
account to obtain the most optimal VGOS schedules. The flux density catalog
presented here is expected to be of use for both astronomy and geodesy. We plan
to continue the monitoring program, potentially adding more sources and more
participating stations.


Description:
    fluxdensity.dat: Flux densities in four frequency bands for 361 ICRF3 
    sources. The frequency bands are: 3.00-3.48 GHz (Band 1), 5.24-5.72 GHz
    (Band 2), 6.36-6.84 GHz (Band 3) and 10.20-10.68 GHz (Band 4). These
    correspond to band A-D of the operational frequency setup used for geodetic
    VLBI with the VLBI Global Observing System (VGOS). Measurements were made by
    the Onsala Twin Telescopes, operating as a single-baseline interferometer.
    As the baseline length is 75 meters, the flux density is probed on
    arc minute scales.
    ivsnames.dat: Translation table matching IERS and IVS source designation.
    Only includes sources where the two designations are different.
	
File Summary:
--------------------------------------------------------------------------------
 FileName      Lrecl    Records    Explanations
--------------------------------------------------------------------------------
ReadMe            80      104     This file
fluxdensity.dat   79      8318    Flux densities in four frequency bands for
                                  361 ICRF3 sources
ivsnames.dat      16      32      Translation table between IERS and IVS source
                                  names
--------------------------------------------------------------------------------

Byte-by-byte Description of file: fluxdensity.dat
--------------------------------------------------------------------------------
   Bytes Format Units Label   Explanations
--------------------------------------------------------------------------------
   1-  8 A8     ---     source  IERS name of source (1)
   9- 15 F7.1   d       MJD     epoch (modified Julian days)
  16- 25 A10    ---     date    epoch (date)
  26- 31 A6     ---     exp     experiment name
  32- 37 F6.3   Jy      flux1   ? Flux density in 3.2 GHz band (VGOS-OP band A)
  38- 43 F6.3   Jy      flux2   ? Flux density in 5.5 GHz band (VGOS-OP band B)
  44- 49 F6.3   Jy      flux3   ? Flux density in 6.6 GHz band (VGOS-OP band C)
  50- 55 F6.3   Jy      flux4   ? Flux density in 10.4 GHz band (VGOS-OP band D)
  56- 61 F6.3   Jy      err1    ? Uncertainty in 3.2 GHz band (VGOS-OP band A)
  62- 67 F6.3   Jy      err2    ? Uncertainty in 5.5 GHz band (VGOS-OP band B)
  68- 73 F6.3   Jy      err3    ? Uncertainty in 6.6 GHz band (VGOS-OP band C)
  74- 79 F6.3   Jy      err4    ? Uncertainty in 10.4 GHz band (VGOS-OP band D)
--------------------------------------------------------------------------------
Note (1):   International Earth Rotation and Reference Systems Service. 
            Format is HHHH+DDD, describing B1950.0 coordinates.
            SIMBAD-searchable as "IERS B[HHHH+DDD]". The exception is 1409+524,
            Which is searchable as "3C295"
--------------------------------------------------------------------------------

Byte-by-byte Description of file: ivsnames.dat
--------------------------------------------------------------------------------
   Bytes Format  Units         Label     Explanations
--------------------------------------------------------------------------------
   1-  8  A8     ---           IERS_name Source name used by the IERS (1)
   9- 16  A8     ---           IVS_name  Source name used by the IVS (2)
-----
Note (1): International Earth Rotation and Reference Systems Service
Note (2): International VLBI Service for Geodesy and Astrometry
================================================================================
(End)                            Alva Kinman [OSO, Sweden]     17-Jul-2026
