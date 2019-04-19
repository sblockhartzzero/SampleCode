#ifndef L1A_HAWKEYE_H
#define L1A_HAWKEYE_H

#ifdef __cplusplus
extern "C" {
#endif
    
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <netcdf.h>
#include "mfhdf.h"
#include "passthebuck.h"
#include "l1_struc.h"
#include "filehandle.h"
    
// Call sequence
// openl1a_hawkeye      > qc_hawkeye_CCD_T      > nan_wmean                     // once
// readl1a_hawkeye      > read_cal_hawkeye                                      // once
// readl1a_hawkeye      > interp_hawkeye_CCD_T  > prep_for_interp_double        // per scan
// readl1a_hawkeye      > calibrate_hawkeye     > prep_for_interp_double        // per scan

// Open, read, close
int openl1a_hawkeye(filehandle *file);
int readl1a_hawkeye(filehandle *file, int32 recnum, l1str *l1rec);
int closel1a_hawkeye(filehandle *file);
// Data cleansing (one-time)
void qc_hawkeye_CCD_T();
// Calibration-related (one-time)
void read_cal_hawkeye(char *cal_path);
// Calibration-related (per scan)
void interp_hawkeye_CCD_T(double delta_time_ms_i);
int calibrate_hawkeye(double current_julian_date);
// Misc utilities
double nan_wmean(double *weights, double *data, size_t n, double fill_value);
void prep_for_interp_double(double *xi, double *x, int N);
double ***alloc3d_double(size_t m, size_t n, size_t o);
void free3d_all_double(double ***p, size_t m, size_t n, size_t o);

    
#ifdef __cplusplus
}
#endif
    
#endif