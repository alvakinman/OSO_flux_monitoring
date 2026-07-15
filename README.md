# OSO_flux_monitoring
Data from AGN flux density monitoring at Onsala Space Observatory.

If you make use of this data for a publication, please cite Kinman et al. 2026 {... reference details to be added ...}.

The directory *lightcurves* contains flux density measurements dating back to 2023-01. The directory will be updated regularly.

Last included monitoring epochs: vo6119 (2026-04-29) and fm6120 (2026-04-30).

There are two experiment types.
VO - VGOS-operational sessions, coordinated by the IVS. Only the ONSA13NE-ONSA13SW baseline (75 m) is analysed.
FM - Flux monitoring sessions. Conducted locally at Onsala Space Observatory every month, targeting the ICRF3 defining sources. 

## Latest lightcurve data
Latest batch created at ...

Number of sources: ...

List of included experiments:

*end of experiment list

## Table structure
Data files are in csv format with the following columns:

|COLUMN    |UNIT    |DESCRIPTION|
|----------|--------|-----------|
|time      |MJD     |Observation time (average of all scans) in Modified Julian Days.|
|date      |---     |Date in YYYY-MM-DD format corresponding to the time column.|
|exp       |---     |Name of experiment|
|flux1     |Jy      |Flux density in VGOS-OP band A (3.2 GHz)|
|flux2     |Jy      |Flux density in VGOS-OP band B (5.5 GHz)|
|flux3     |Jy      |Flux density in VGOS-OP band C (6.6 GHz)|
|flux4     |Jy      |Flux density in VGOS-OP band D (10.4 GHz)|
|err1      |Jy      |Uncertainty of flux density in band A. Corresponds to the|
|         |       | error bars in Kinman et al. (2026).|
|err2      |Jy     | Uncertainty of flux density in band B.|
|err3      |Jy     | Uncertainty of flux density in band C.|
|err4      |Jy     | Uncertainty of flux density in band D.|


## Source name translation table
Data files are named {sourcename}.csv, where *sourcename* is the name in use by the IVS. For most sources, this corresponds to the IERS name/B1950 coordinates. For the sources where this is not the case, the translation between IVS and IERS name can be found in the table below.

*start of translation table

*end of translation table
