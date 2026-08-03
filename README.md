# OSO_flux_monitoring
Data from AGN flux density monitoring at Onsala Space Observatory.

If you make use of this data, please cite Kinman et al. 2026, "A multi-band radio flux density catalog of ICRF3 sources using the Onsala Twin Telescopes", arXiv:2607.20671. DOI: https://doi.org/10.48550/arXiv.2607.20671

The directory *lightcurves* contains flux density measurements dating back to 2023-01. The directory will be updated regularly.

### There are two experiment types:
VO - VGOS-operational sessions, coordinated by the IVS. Only the ONSA13NE-ONSA13SW baseline (75 m) is analysed.

FM - Flux monitoring sessions. Conducted locally at Onsala Space Observatory every month, targeting the ICRF3 defining sources. 

## Latest light curve data
Latest batch created 2026-07-15

Latest observation date included: 2026-04-30 (fm6120)

Number of sources: 361

List of included experiments:

fm3031, fm3134, fm3148, fm3182, fm3195, fm3293, vo3299, fm3306, fm4050, fm4074, vo4108, vo4115, fm4128, vo4150, fm4155, vo4248, vo4283, vo4311, fm4312, vo4339, fm4340, fm4341, vo4346, vo5022, vo5029, vo5043, fm5052, fm5053, fm5087, vo5099, fm5114, vo5134, vo5176, fm5183, fm5197, vo5232, fm5264, vo5281, vo5302, fm5306, fm5331, vo5351, fm5352, vo6014, fm6029, vo6049, fm6054, vo6056, vo6091, fm6092, vo6119, fm6120

*end of experiment list

## Table structure
Data files in the directory *lightcurves* are in csv format with the following columns:

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

Note: In the combined file *fluxdensity_{date}.csv*, the first column is a "source" column with IERS names (see translation table below).


## Source name translation table
Data files are named {sourcename}.csv, where *sourcename* is the name in use by the IVS. For most sources, this corresponds to the IERS name/B1950 coordinates. For the sources where this is not the case, the translation between IVS and IERS name can be found in the table below.

*start of translation table
|IERS name |IVS name|
|--------|--------|
|0007+106|IIIZW2|
|0046+316|NGC0262|
|0305+039|NGC1218|
|0316+413|3C84|
|0336-019|CTA26|
|0355+508|NRAO150|
|0538+498|3C147|
|0851+202|OJ287|
|0923+392|4C39.25|
|1222+131|M84|
|1226+023|3C273B|
|1228+126|3C274|
|1253-055|3C279|
|1328+307|3C286|
|1409+524|3C295|
|1633+382|1633+38|
|1637+826|NGC6251|
|1638+398|NRAO512|
|1652+398|DA426|
|1730-130|NRAO530|
|1807+698|3C371|
|1901+319|3C395|
|2017+745|2017+743|
|2037+511|3C418|
|2200+420|VR422201|
|2223-052|3C446|
|2230+114|CTA102|
|2250+190|2250+194|
|2251+158|3C454.3|

*end of translation table
