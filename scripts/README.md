# OSO_flux_monitoring
Scripts for AGN flux density monitoring at Onsala Space Observatory.

If you make use of this data for a publication, please reference Kinman et al. 2026 {... reference details to be added ...}.

cal_fm_kinman2026.py: CASA pipeline for obtaining flux densitites from a "FM" experiment involving only Oe and Ow.

cal_vo_kinman2026.py: CASA pipeline for obtaining flux densities from a VGOS-OP (VO) experiment. Flux densities will only be for the Oe-Ow baseline, but the fitsidi file may contain other antennas. A separate flux calibrator fitsidi-file is needed.
