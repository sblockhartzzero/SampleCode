/* ============================================================================ */
/* module l1a_hawkeye.c - functions to read HAWKEYE L1A for MSL12               */
/* Written By:  Steve Lockhart SAIC May 2018                                    */
/*      -Started with l1_viirs_nc.c and made changes to support hawkeye.        */
/*                                                                              */
/* ============================================================================ */

/* Issues:
    open:
       -Why is file->nscan = nline = nscan * MBAND_NUM_DETECTORS? Setting MBAND_NUM_DETECTORS = 1 for hawkeye.    
       -What is spatial resolution? 120m?
    missing fields:
        -Commenting out HAM stuff, extract_pixel_start/stop stuff, as it is currently not in L1A. 
        -Setting orbit number to 0, as it is currently not in L1A.
        -scanQualityId skipped for now
        -Also skipping att angle, so ripples to "Compute polarization rotation angles"
    general: 
        -Clean up commented out sections (i.e. the ones that are VIIRS-specific)
        -What are finder_pixels, finder_lines in L1A file?
        -What are dark_pixels in L1A file? 
        -Does last line get run twice?                
*/

// Includes from l1_viirs_nc.c
#include <nc4utils.h>
#include "l12_proto.h"
#include "libnav.h"
#include <productInfo.h>

// Includes from demo_cal.c
#include <stdio.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <stdlib.h>
#include <fcntl.h>
#include <string.h>
#include <errno.h>
#include <unistd.h>
#include <libgen.h>
#include <math.h>
#include <timeutils.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_sort.h>
#include <gsl/gsl_statistics.h>
#include <stdbool.h>

// New includes
#include "l1a_hawkeye.h"


// MBAND_NUM_DETECTORS changed to 1 for hawkeye
#define MBAND_NUM_DETECTORS 1


// Declarations from l1_viirs_nc.c
static int32 prevScan = -1;

static int geoFileId;
static int geoNavigationGrp;
static int geoGeolocationGrp;
static int geoScanLineGrp;
static int l1aScanLineGrp;
static int l1aTelemetryGrp;
static int l1aEarthViewGrp;
static int lonId, latId, senzId, senaId, solzId, solaId, angId, posId, velId, HAMSideId, scanQualityId, pixelQualityId;

static double ***f_cal_corr = NULL; // f table correction [band][det][ms] 

static float *Fobar; // reflectance to radiance conversion factors
static int extract_pixel_start = 0;
static int extract_pixel_stop = 0;

static short *tmpShort;
static unsigned char *tmpByte;
static int nline;
static size_t num_scans, num_pixels, num_bands;

static int firstCall = 1;
static double starttime;
static double lastvalidtime;
static int lastvalidscan = 0;
static double time_interval;

static float latGeoFillValue = -999.9;
static float lonGeoFillValue = -999.9;
static short senzGeoFillValue = -32768;
static short senaGeoFillValue = -32768;
static short solzGeoFillValue = -32768;
static short solaGeoFillValue = -32768;

static float latL2FillValue = -999.0;
static float lonL2FillValue = -999.0;
static float senzL2FillValue = -32767;
static float senaL2FillValue = -32767;
static float solzL2FillValue = -32767;
static float solaL2FillValue = -32767;

// Declarations added for hawkeye L1A.
static size_t num_ccd_temps, num_tlm_blocks, num_tlm_blocks_cleansed;
static double **CCD_temperatures;                                   //[num_tlm_blocks][num_ccd_temps]
static double **CCD_temperatures_cleansed;                          //up to [num_tlm_blocks][num_ccd_temps]
static double CCD_temperature_default = 35.0;                       // same as K3T[3], used iff all CCD_temperatures are FILL
static double *tlm_delta_time_ms;                                   //[num_tlm_blocks], ms since start of image
static double *tlm_delta_time_ms_cleansed;                          // up to [num_tlm_blocks]
static float fill_value_in_CCD_T;
static int *dn_varid;                                               //[num_bands]
static double *CCD_temperatures_this_scan;                          //[num_ccd_temps] i.e. interpolated in time
static double **dn;                                                 //[num_bands][num_pixels], This scan, before cal
static double **Lt;                                                 //[num_bands][num_pixels], This scan, after cal

// Declarations added for hawkeye cal coeffs.
static size_t cal_num_dates;
static size_t cal_num_temperatures;
static double *time_ref;                                            //[cal_num_dates]
static double *temperature_ref;                                     //[cal_num_temperatures]
static double *K1;                                                  //Gain at nadir pixel, [num_bands]   
static double ***K2;                                                //Temporal correction, BEFORE interpolation, [num_bands][num_pixels][cal_num_dates] 
static double **K3;                                                 //Temperature correction BEFORE interpolation, [num_bands][cal_num_temperatures] 
static float **K4;                                                  //RVS, [num_bands][num_pixels]
static double **K5;                                                 //Non-linearity (quadratic term), [num_bands][num_pixels]

/**
 * Open the hawkeye L1A file and perform some one-time tasks (as opposed to tasks that
 * are per scan), including:
 *   -Get L1A dimensions num_scans, num_bands, num_pixels, num_ccd_temps, num_tlm_blocks (static).
 *    Allocate memory for some static arrays based upon these dimensions.
 *   -Get L1A group ids e.g. l1aScanLineGrp, l1aEarthViewGrp and it's 8 "band_%d" var ids (static)
 *   -Get L1A CCD_temperatures and tlm_delta_time_ms, and call qc_hawkeye_CCD_T to create cleansed 
 *    versions of these arrays (static).
 *   -Get L1A "time_coverage_start" and "time_coverage_end" to derive time_interval (static) etc.
 *

 * Get   
 * @param file
 * @return 
 */
int openl1a_hawkeye(filehandle * file) {
    char *fltime;

    size_t tmpSize;
    int ncid_L1A, varid, dimid, status;
    size_t att_len;
    int orbit_number;
    int band_num;
    char nc_search_string[6] = "";              //band_X
    
    // Initialize vars for reading arrays from L1A file
    size_t start[] = { 0, 0, 0 };
    size_t count[] = { 1, 1, 1 };
    
    // Open the netcdf4 input file
    printf("Opening hawkeye l1a file\n");
    status = nc_open(file->name, NC_NOWRITE, &ncid_L1A);
    if (status == FAIL) {
        fprintf(stderr, "-E- %s line %d: nc_open(%s) failed.\n",
                __FILE__, __LINE__, file->name);
        return (1);
    }
    
    // Get dims from L1A file: num_scans, num_bands, num_pixels, num_ccd_temps, num_tlm_blocks
    // num_scans
    status = nc_inq_dimid(ncid_L1A, "number_of_scans", &dimid);
    if (status != NC_NOERR) {
        fprintf(stderr, "-E- Error reading num_scans.\n");
        exit(EXIT_FAILURE);
    };
    nc_inq_dimlen(ncid_L1A, dimid, &num_scans);
    // num_bands
    status = nc_inq_dimid(ncid_L1A, "number_of_bands", &dimid);
    if (status != NC_NOERR) {
        fprintf(stderr, "-E- Error reading num_bands.\n");
        exit(EXIT_FAILURE);
    };
    nc_inq_dimlen(ncid_L1A, dimid, &num_bands);
    // num_pixels
    status = nc_inq_dimid(ncid_L1A, "number_of_pixels", &dimid);
    if (status != NC_NOERR) {
        fprintf(stderr, "-E- Error reading num_pixels.\n");
        exit(EXIT_FAILURE);
    };
    nc_inq_dimlen(ncid_L1A, dimid, &num_pixels);
    // num_pixels
    status = nc_inq_dimid(ncid_L1A, "ccd_temps", &dimid);
    if (status != NC_NOERR) {
        fprintf(stderr, "-E- Error reading num_ccd_temps.\n");
        exit(EXIT_FAILURE);
    };
    nc_inq_dimlen(ncid_L1A, dimid, &num_ccd_temps);
    // num_tlm_blocks
    status = nc_inq_dimid(ncid_L1A, "number_of_tlm_blocks", &dimid);
    if (status != NC_NOERR) {
        fprintf(stderr, "-E- Error reading num_tlm_blocks.\n");
        exit(EXIT_FAILURE);
    };
    nc_inq_dimlen(ncid_L1A, dimid, &num_tlm_blocks);
    
    // Derive nline (See issue above.) This line is from l1_viirs_nc.c.
    nline = num_scans * MBAND_NUM_DETECTORS;
    
    // Now that we know dims, prep for additional hawkeye one-time reads by allocating memory
    // for arrays (declared static above). (Memory is freed in closel1a_hawkeye.)
    CCD_temperatures = alloc2d_double(num_ccd_temps,num_tlm_blocks);
    CCD_temperatures_cleansed = alloc2d_double(num_ccd_temps,num_tlm_blocks);
    tlm_delta_time_ms = (double *) calloc(num_tlm_blocks, sizeof (double));
    tlm_delta_time_ms_cleansed = (double *) calloc(num_tlm_blocks, sizeof (double));
    dn_varid = (int *) calloc(num_bands, sizeof (int));
    
    // Get group id from L1A file for GROUP scan_line_attributes.
    if ((nc_inq_grp_ncid(ncid_L1A, "scan_line_attributes", &l1aScanLineGrp)) == NC_NOERR) {
    } else {
        fprintf(stderr, "-E- Error finding scan_line_attributes.\n");
        exit(EXIT_FAILURE);
    }   
    
    // get mirror side (from l1_viirs_nc.c) from this GROUP
    //status = nc_inq_varid(l1aScanLineGrp, "HAM_side", &HAMSideId);
    //check_err(status, __LINE__, __FILE__);    
    
    // Get CCD_temperatures, tlm_delta_time_ms from L1A file from GROUP parameters_telemetry_data
    if ((nc_inq_grp_ncid(ncid_L1A, "parameters_telemetry_data", &l1aTelemetryGrp)) == NC_NOERR) {
    } else {
        fprintf(stderr, "-E- Error finding parameters_telemetry_data.\n");
        exit(EXIT_FAILURE);
    }
    // Get CCD_temperatures from this GROUP 
    count[0] = num_tlm_blocks;                                                                  // 1 scan at a time
    count[1] = num_ccd_temps;
    count[2] = 0;
    status = nc_inq_varid(l1aTelemetryGrp, "CCD_temperatures", &varid);
    if (status != NC_NOERR) {
        fprintf(stderr, "-E- Error finding CCD_temperatures.\n");
        exit(EXIT_FAILURE);
    } 
    status = nc_get_vara_double(l1aTelemetryGrp, varid, start, count, &CCD_temperatures[0][0]);
    if (status != NC_NOERR) {
        fprintf(stderr, "-E- Error reading CCD_temperatures.\n");
        exit(EXIT_FAILURE);
    };
    // Determine what fill value was used for CCD_temperatures
    status = nc_inq_varid(l1aTelemetryGrp, "CCD_temperatures", &varid);
    if (status != NC_NOERR) {
        printf("Lost CCD_temperatures; status = %d\n", status);
        exit(EXIT_FAILURE);
    }
    status = nc_inq_var_fill(l1aTelemetryGrp, varid, NULL, &fill_value_in_CCD_T);
    if (status != NC_NOERR) {
        fprintf(stderr, "-W- Could not get FILL value for CCD_temperatures, so assuming it is -999.0.\n");
        fill_value_in_CCD_T = -999.0; 
    }
    // Get tlm_time_stamp from this GROUP
    count[0] = num_tlm_blocks;                                                                  // 1 scan at a time
    count[1] = 0;
    count[2] = 0;
    status = nc_inq_varid(l1aTelemetryGrp, "tlm_time_stamp", &varid);
    if (status != NC_NOERR) {
        fprintf(stderr, "-E- Error finding tlm_delta_time_ms.\n");
        exit(EXIT_FAILURE);
    } 
    status = nc_get_vara_double(l1aTelemetryGrp, varid, start, count, &tlm_delta_time_ms[0]);
    if (status != NC_NOERR) {
        fprintf(stderr, "-E- Error reading tlm_delta_time_ms.\n");
        exit(EXIT_FAILURE);
    };
    
    // Cleanse CCD_temperatures (and corresponding tlm_delta_time_ms), handling fill values
    // This sets static vars CCD_temperatures_cleansed, tlm_delta_time_ms_cleansed, num_tlm_blocks_cleansed
    qc_hawkeye_CCD_T();   

    // Get attribute values (from l1_viirs_nc.c) e.g. "time_coverage_start" and "time_coverage_end" to derive
    // time_interval etc. Note that time_coverage_end is the START of the last scan.
    nc_type vr_type;            // attribute type 
    size_t vr_len;              // attribute length 
    // Comment out extract_pixel_start/stop stuff
    /*
    if ((nc_inq_att(fileID, NC_GLOBAL, "extract_pixel_start", &vr_type, &vr_len) == 0)) {
        status = nc_get_att_int(fileID, NC_GLOBAL, "extract_pixel_start", &extract_pixel_start);
        check_err(status, __LINE__, __FILE__);
        extract_pixel_start--; // Attribute is one-based
        status = nc_get_att_int(fileID, NC_GLOBAL, "extract_pixel_stop", &extract_pixel_stop);
        check_err(status, __LINE__, __FILE__);
        extract_pixel_stop--; // Attribute is one-based
        if (npix != (extract_pixel_stop - extract_pixel_start + 1)) {
            fprintf(stderr, "-E- Problem with the extracted L1A file pixel dimension.\n");
            printf("    npix(%d), extract_pixel_stop(%d), extract_pixel_start(%d) do not work together.\n",
                    npix, extract_pixel_stop, extract_pixel_start);
            exit(EXIT_FAILURE);
        }
    }
    */

    if (want_verbose) {
        printf("Hawkeye L1A Npix  :%d Nlines:%d\n", num_pixels, nline);
    } // want_verbose

    // get start and end time
    status = nc_inq_attlen(ncid_L1A, NC_GLOBAL, "time_coverage_start", &att_len);
    check_err(status, __LINE__, __FILE__);

    // allocate required space before retrieving values 
    fltime = (char *) malloc(att_len + 1); // + 1 for trailing null 

    // get attribute values 
    status = nc_get_att_text(ncid_L1A, NC_GLOBAL, "time_coverage_start", fltime);
    check_err(status, __LINE__, __FILE__);
    fltime[att_len] = '\0';
    //    isodate2ydmsec(fltime, (int32_t *) &year,(int32_t *) &day, (int32_t *) &msec);

    // Convert "time_coverage_start" ISO string to unix (seconds since 1/1/1970)
    starttime = lastvalidtime = isodate2unix(fltime);
    free(fltime);

    status = nc_inq_attlen(ncid_L1A, NC_GLOBAL, "time_coverage_end", &att_len);
    check_err(status, __LINE__, __FILE__);

    // allocate required space before retrieving values 
    fltime = (char *) malloc(att_len + 1); // + 1 for trailing null

    // get attribute values 
    status = nc_get_att_text(ncid_L1A, NC_GLOBAL, "time_coverage_end", fltime);
    check_err(status, __LINE__, __FILE__);
    fltime[att_len] = '\0';

    // Convert "time_coverage_stop" ISO string to unix (seconds since 1/1/1970)
    double stoptime = isodate2unix(fltime);
    free(fltime);

    // time_interval may be used in readl1a_hawkeye if there is not a good scan time (per scan)
    time_interval = (stoptime - starttime)/(num_scans-1); // secs per scan

    if ((nc_inq_att(ncid_L1A, NC_GLOBAL, "orbit_number", &vr_type, &vr_len) == 0)) {
        status = nc_get_att_int(ncid_L1A, NC_GLOBAL, "orbit_number", &orbit_number);
        check_err(status, __LINE__, __FILE__);
    } else {
        orbit_number = 0;
    }

    
    // Identify the "earth_view_data" GROUP and its "band_%d" vars, to be used later by readl1a_hawkeye
    // Store the ids in static variables so we don't have to do nc_inq_grp_ncid per scan.
    if ((nc_inq_grp_ncid(ncid_L1A, "earth_view_data", &l1aEarthViewGrp)) == NC_NOERR) {
        for (band_num=0; band_num<num_bands; band_num++) {
            sprintf(nc_search_string, "band_%d", band_num+1);
            status = nc_inq_varid(l1aEarthViewGrp, nc_search_string, &dn_varid[band_num]);
            if (status != NC_NOERR) {
                fprintf(stderr, "-E- Error finding %s.\n", nc_search_string);
                exit(EXIT_FAILURE);
            }
        }
    } else {
        fprintf(stderr, "-E- Error finding earth_view_data.\n");
        exit(EXIT_FAILURE);
    }

    file->sd_id = ncid_L1A;
    file->nbands = num_bands;
    file->npix = num_pixels;
    file->nscan = nline;
    file->ndets = MBAND_NUM_DETECTORS;
    file->terrain_corrected = 1; // presumed.
    file->orbit_number = orbit_number;
    strcpy(file->spatialResolution, "120 m");

    rdsensorinfo(file->sensorID, input->evalmask,
            "Fobar", (void **) &Fobar);

    if (want_verbose)
        printf("file->nbands = %d\n", (int) file->nbands);

    // Setup geofile pointers
    if (file->geofile[0] != '\0') {

        status = nc_open(file->geofile, NC_NOWRITE, &geoFileId);
        if (status != NC_NOERR) {
            fprintf(stderr, "-E- Could not open GEO file \"%s\"\n", file->geofile);
            exit(EXIT_FAILURE);
        }

        status = nc_inq_grp_ncid(geoFileId, "geolocation_data", &geoGeolocationGrp);
        check_err(status, __LINE__, __FILE__);
        status = nc_inq_varid(geoGeolocationGrp, "longitude", &lonId);
        check_err(status, __LINE__, __FILE__);
        status = nc_inq_var_fill(geoGeolocationGrp, lonId, NULL, &lonGeoFillValue);
        check_err(status, __LINE__, __FILE__);
        status = nc_inq_varid(geoGeolocationGrp, "latitude", &latId);
        check_err(status, __LINE__, __FILE__);
        status = nc_inq_var_fill(geoGeolocationGrp, latId, NULL, &latGeoFillValue);
        check_err(status, __LINE__, __FILE__);
        status = nc_inq_varid(geoGeolocationGrp, "sensor_zenith", &senzId);
        check_err(status, __LINE__, __FILE__);
        status = nc_inq_var_fill(geoGeolocationGrp, senzId, NULL, &senzGeoFillValue);
        check_err(status, __LINE__, __FILE__);
        status = nc_inq_varid(geoGeolocationGrp, "sensor_azimuth", &senaId);
        check_err(status, __LINE__, __FILE__);
        status = nc_inq_var_fill(geoGeolocationGrp, senaId, NULL, &senaGeoFillValue);
        check_err(status, __LINE__, __FILE__);
        status = nc_inq_varid(geoGeolocationGrp, "solar_zenith", &solzId);
        check_err(status, __LINE__, __FILE__);
        status = nc_inq_var_fill(geoGeolocationGrp, solzId, NULL, &solzGeoFillValue);
        check_err(status, __LINE__, __FILE__);
        status = nc_inq_varid(geoGeolocationGrp, "solar_azimuth", &solaId);
        check_err(status, __LINE__, __FILE__);
        status = nc_inq_var_fill(geoGeolocationGrp, solaId, NULL, &solaGeoFillValue);
        check_err(status, __LINE__, __FILE__);
        status = nc_inq_varid(geoGeolocationGrp, "quality_flag", &pixelQualityId);
        check_err(status, __LINE__, __FILE__);

        status = nc_inq_grp_ncid(geoFileId, "navigation_data", &geoNavigationGrp);
        check_err(status, __LINE__, __FILE__);
        //angId not yet in GEO file
        //status = nc_inq_varid(geoNavigationGrp, "att_ang_mid", &angId);
        //check_err(status, __LINE__, __FILE__);
        status = nc_inq_varid(geoNavigationGrp, "orb_pos", &posId);
        check_err(status, __LINE__, __FILE__);
        status = nc_inq_varid(geoNavigationGrp, "orb_vel", &velId);
        check_err(status, __LINE__, __FILE__);

        status = nc_inq_grp_ncid(geoFileId, "scan_line_attributes", &geoScanLineGrp);
        check_err(status, __LINE__, __FILE__);
        // The following field is not yet in the GEO file
        //status = nc_inq_varid(geoScanLineGrp, "scan_quality", &scanQualityId);
        //check_err(status, __LINE__, __FILE__);
    } // geofile

    // Setup the fill values for the geo products
    // Replaced VIIRSN with HAWKEYE. 
    productInfo_t* info = allocateProductInfo();
    status = findProductInfo("lat", HAWKEYE, info);
    if (status)
        latL2FillValue = info->fillValue;
    status = findProductInfo("lon", HAWKEYE, info);
    if (status)
        lonL2FillValue = info->fillValue;
    status = findProductInfo("sena", HAWKEYE, info);
    if (status)
        senaL2FillValue = info->fillValue;
    status = findProductInfo("senz", HAWKEYE, info);
    if (status)
        senzL2FillValue = info->fillValue;
    status = findProductInfo("sola", HAWKEYE, info);
    if (status)
        solaL2FillValue = info->fillValue;
    status = findProductInfo("solz", HAWKEYE, info);
    if (status)
        solzL2FillValue = info->fillValue;
    freeProductInfo(info);

    return (LIFE_IS_GOOD);
}


/**
 * Read the specified scan line from the specified L1A file. If this is the first call, 
 * call read_cal_hawkeye to read the calibration file, which stores calibration coefficients
 * (for instrument calibration) in static arrays.
 * For each scan, get scan_delta_time_ms and use it to call interp_hawkeye_CCD_T, which sets 
 * CCD_temperatures_this_scan[num_ccd_temps]. Then, apply the instrument cal by calling calibrate_hawkeye.
 * Store Lt in l1rec->Lt.
 * Read GEO file.
 * 
 * 
 * @param file
 * @param line
 * @param l1rec
 * @return 
 */
int readl1a_hawkeye(filehandle *file, int32 line, l1str *l1rec) {

    // Declarations from l1_viirs_nc.c
    int32 ip, ib, ipb;
    int i;
    double scan_sec;
    int16 scan_year, scan_day;
    double f_corr;

    int status;
    size_t start[] = { 0, 0, 0 };
    size_t count[] = { 1, 1, 1 };

    int32 scan = line / MBAND_NUM_DETECTORS;
    
    // Additional declarations for hawkeye
    int band_num, pixel_num, varid;
    int16 scan_month, scan_dom;
    double current_julian_date;
    double scan_delta_time_ms;                                      // ms since start of image per scan

    //printf("reading hawkeye l1a file\n");
    for (ip = 0; ip < num_pixels; ip++) {
        l1rec->pixnum[ip] = ip + extract_pixel_start;
    }

    // If first call, 
    if (firstCall) {
        firstCall = 0;  
            
        // Get hawkeye calibration coefficients from cal file
        // This sets static cal coeffs: K1-K5 as well as time_ref and temperature_ref 
        // It should be done only once (i.e. for firstCall)
        read_cal_hawkeye(input->calfile);

        // One-time memory allocations
        tmpShort = (short *) malloc(num_pixels * sizeof(short));
        tmpByte = (unsigned char *) malloc(num_pixels);
        dn = alloc2d_double(num_pixels,num_bands);                        
        Lt = alloc2d_double(num_pixels,num_bands);
        CCD_temperatures_this_scan = (double *) calloc(num_ccd_temps, sizeof (double));

    }
    
    //    l1rec->sensorID = file->sensorID;
    l1rec->npix = file->npix;

    // Time
    // Get delta_time 
    start[0] = line;
    count[0] = 1;                                                                       // 1 scan at a time
    count[1] = 0;
    count[2] = 0;
    status = nc_inq_varid(l1aScanLineGrp, "delta_time", &varid);
    if (status != NC_NOERR) {
        fprintf(stderr, "-E- Error finding scan_delta_time_ms.\n");
        exit(EXIT_FAILURE);
    } 
    status = nc_get_vara_double(l1aScanLineGrp, varid, start, count, &scan_delta_time_ms);
    if (status != NC_NOERR) {
        fprintf(stderr, "-E- Error reading scan_delta_time_ms.\n");
        exit(EXIT_FAILURE);
    };
    // Set lastvalidtime, the start of this scan (in secs since 1/1/1970)
    if (scan_delta_time_ms > 0) {
        // starttime is the start of this granule (L1A "time_coverage_start"), converted to secs since 1/1/1970
        // scan_delta_time_ms is when this scan starts (number of milliseconds since the start of this granule)
        lastvalidtime = starttime + scan_delta_time_ms/1000;
    } else {
        lastvalidtime = lastvalidtime + (time_interval * (line - lastvalidscan));      
    }
    // Set scan_year, scan_day, scan_sec
    unix2yds(lastvalidtime, &scan_year, &scan_day, &scan_sec);
    // Store lastvalidtime in l1rec->scantime
    l1rec->scantime = lastvalidtime; 
    // Convert lastvalidtime to julian date
    yd2md(scan_year, scan_day, &scan_month, &scan_dom);
    current_julian_date = jday(scan_year, scan_month, scan_dom) + (scan_sec/(24*3600.0));
    if (scan == 0) {
        //printf("current_julian_date = %f, and scan_delta_time_ms = %f\n", current_julian_date, scan_delta_time_ms);
    }
    
    // Before calibrating this hawkeye scan, we need to interpolate the CCD_temperatures_cleansed (in time). 
    // The resulting interpolated array will be CCD_temperatures_this_scan[num_ccd_temps], which is needed by 
    // calibrate_hawkeye. Note that the interpolation routine requires the arrays to be of type double.
    if (line>=num_scans-1) {
        //printf("line %d has scan_delta_time_ms=%f\n", line, scan_delta_time_ms);
    }    
    interp_hawkeye_CCD_T(scan_delta_time_ms);
    
    // GROUP earth_view_data
    // Get dn[num_bands][num_pixels] from this GROUP just for THIS scan  
    start[0] = line;
    count[0] = 1;                                                               // 1 line at a time
    count[1] = num_pixels;
    count[2] = 0;
    for (band_num=0; band_num<num_bands; band_num++) {       
        status = nc_get_vara_double(l1aEarthViewGrp, dn_varid[band_num], start, count, &dn[band_num][0]);
        if (status != NC_NOERR) {
            fprintf(stderr, "-E- Error reading band %d data from the L1A file\n", band_num);
            exit(EXIT_FAILURE);
        };
    }    
    // Reset start[0] in case any other nc read is expecting (assuming) it to be the default of 0.
    start[0] = 0;

    // Apply instrument calibration to this scan
    // dn is BEFORE cal; Lt is AFTER cal.
    calibrate_hawkeye(current_julian_date);
    
    // Set l1rec->Lt and fudge lat/lon
    for (pixel_num = 0; pixel_num < num_pixels; pixel_num++) {
        //l1rec->lat[pixel_num] = 5*((float)scan/(float)num_scans);
        //l1rec->lon[pixel_num] = 5*((float)pixel_num/(float)num_pixels);
        for (band_num = 0; band_num < num_bands; band_num++) {
            //ipb = ip * nbands + ib
            l1rec->Lt[pixel_num*num_bands + band_num] = (float)Lt[band_num][pixel_num];
        }
    }

    

    
    //------------------------------------
    // if there is no geo file just return
    // This is used for l1info with only a L1A file and no GEO
    //-------------------------------------
    if (file->geofile[0] == '\0') {
        return 0;
    }

    // first check the scan quality flag
    // 1   SCE_side_A_B
    // 2   SCE_side_invalid
    // 4   Sector_rotation
    // 8   Encoder_degraded
    // 16  SAA
    // 32  Solar_eclipse
    // 64  Lunar_eclipse
    // 128 HAM_side
    //
    // Sector_rotation
    short scanQualityWarnMask = 2 | 8 | 128;
    short scanQualityFailMask = 4;
    static short scanQualityFlag = 0;

    // Skip because scanQualityId is not yet in the GEO file [SBL]
    /* 
    if (scan != prevScan) {
        start[0] = scan;
        status = nc_get_var1_short(geoScanLineGrp, scanQualityId, start, &scanQualityFlag);
        check_err(status, __LINE__, __FILE__);
    }
    if (scanQualityFlag & scanQualityFailMask) {
        for (ip = 0; ip < num_pixels; ip++)
            l1rec->flags[ip] |= NAVFAIL;
        return 0;
    }
    if (scanQualityFlag & scanQualityWarnMask) {
        for (ip = 0; ip < num_pixels; ip++)
            l1rec->flags[ip] |= NAVWARN;
    }
    */

    static unsigned char HAMSideVal = 0;
    // Commenting out HAM for now, as it's not in L1A yet
    //if (scan != prevScan) {
    //    start[0] = scan;
    //    status = nc_get_var1_uchar(l1aScanLineGrp, HAMSideId, start, &HAMSideVal);
    //    check_err(status, __LINE__, __FILE__);
    //}
    l1rec->mside = HAMSideVal;

    // set up to read all pixels of the line.
    start[0] = line;
    start[1] = 0;
    count[0] = 1;
    count[1] = num_pixels; // read all pixels

    status = nc_get_vara_float(geoGeolocationGrp, latId, start, count, l1rec->lat);
    check_err(status, __LINE__, __FILE__);
    for (i = 0; i < num_pixels; i++)
        if (l1rec->lat[i] == latGeoFillValue)
            l1rec->lat[i] = latL2FillValue;

    status = nc_get_vara_float(geoGeolocationGrp, lonId, start, count, l1rec->lon);
    check_err(status, __LINE__, __FILE__);
    for (i = 0; i < num_pixels; i++)
        if (l1rec->lon[i] == lonGeoFillValue)
            l1rec->lon[i] = lonL2FillValue;

    status = nc_get_vara_short(geoGeolocationGrp, solzId, start, count, tmpShort);
    check_err(status, __LINE__, __FILE__);
    for (i = 0; i < num_pixels; i++)
        if (tmpShort[i] == solzGeoFillValue)
            l1rec->solz[i] = solzL2FillValue;
        else
            l1rec->solz[i] = tmpShort[i] * 0.01;

    status = nc_get_vara_short(geoGeolocationGrp, solaId, start, count, tmpShort);
    check_err(status, __LINE__, __FILE__);
    for (i = 0; i < num_pixels; i++)
        if (tmpShort[i] == solaGeoFillValue)
            l1rec->sola[i] = solaL2FillValue;
        else
            l1rec->sola[i] = tmpShort[i] * 0.01;

    status = nc_get_vara_short(geoGeolocationGrp, senzId, start, count, tmpShort);
    check_err(status, __LINE__, __FILE__);
    for (i = 0; i < num_pixels; i++)
        if (tmpShort[i] == senzGeoFillValue)
            l1rec->senz[i] = senzL2FillValue;
        else
            l1rec->senz[i] = tmpShort[i] * 0.01;

    status = nc_get_vara_short(geoGeolocationGrp, senaId, start, count, tmpShort);
    check_err(status, __LINE__, __FILE__);
    for (i = 0; i < num_pixels; i++)
        if (tmpShort[i] == senaGeoFillValue)
            l1rec->sena[i] = senaL2FillValue;
        else
            l1rec->sena[i] = tmpShort[i] * 0.01;

    // Load Angles
    float ang[3]; // degrees
    float pos[3]; // km
    float vel[3]; // km/sec
    size_t s3[] = { scan, 0 };
    size_t c3[] = { 1, 3 };
    
    //status = nc_get_vara_float(geoNavigationGrp, angId, s3, c3, ang);
    //check_err(status, __LINE__, __FILE__);
    status = nc_get_vara_float(geoNavigationGrp, posId, s3, c3, pos);
    check_err(status, __LINE__, __FILE__);
    status = nc_get_vara_float(geoNavigationGrp, velId, s3, c3, vel);
    check_err(status, __LINE__, __FILE__);
    for (i = 0; i < 3; i++) {
        pos[i] /= 1000.; // m   -> km
        vel[i] /= 1000.; // m/s -> km/s
    }

    // Compute polarization rotation angles
    // Skipping for now since no ang [SBL]
    /* 
    float sen_mat[3][3], coeff[10];
    double mnorm[3];
    ocorient_(pos, vel, ang, sen_mat, coeff);
    for (i = 0; i < 3; i++)
        mnorm[i] = sen_mat[i][0];
    compute_alpha(l1rec->lon, l1rec->lat,
            l1rec->senz, l1rec->sena,
            mnorm, l1rec->npix, l1rec->alpha);
    */

    // Check pixel values 
    status = nc_get_vara_uchar(geoGeolocationGrp, pixelQualityId, start, count, tmpByte);
    check_err(status, __LINE__, __FILE__);
    // 1 Input_invalid
    // 2 Pointing_bad
    // 4 Terrain_bad
    unsigned char qualityFailMask = 1 || 2;
    unsigned char qualityWarnMask = 4;
    for (i = 0; i < num_pixels; i++) {
        if (tmpByte[i] & qualityFailMask)
            l1rec->flags[i] |= NAVFAIL;
        if (tmpByte[i] & qualityWarnMask)
            l1rec->flags[i] |= NAVWARN;
    }

    // Earth-sun distance correction for this scan
    static double esdist = -999.9;
    if (scan != prevScan) {
	int16_t year, day;
	double dsec;
	unix2yds(l1rec->scantime, &year, &day, &dsec);
	int32_t yr = (int32_t) year;
	int32_t dy = (int32_t) day;
	int32_t msec = (int32_t) (dsec * 1000.0);
        esdist = esdist_(&yr, &dy, &msec);
    }
    l1rec->fsol = pow(1.0 / esdist, 2);

    // Skipping this VIIRS-related stuff
    /* 
    int nbands = 16; //, nRSBbands = 10, nCIRbands = 1, nTEBbands = 5;

    // read in calibrated L1B data
    static float *l1bptrs[16] = {0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0,
        0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0};

    static char *bandType[16] = {"RSB", "RSB", "RSB", "RSB", "RSB", "RSB",
        "RSB", "RSB", "CIR", "RSB", "RSB", "TEB",
        "TEB", "TEB", "TEB", "TEB"};

    // Note: l1bptrs arrays are 3200 pixels wide
    int oldVerbose = want_verbose;
    want_verbose = 0;
    VcstViirsCal_calibrateMOD(line, nbands, l1bptrs);
    want_verbose = oldVerbose;
    */

    l1rec->detnum = line % file->ndets;

    // Skipping this VIIRS-related stuff
    /* 
    int irsb = 0, iteb = 0;
    int l1bptrs_scan_idx;
    for (ib = 0; ib < nbands; ib++) {

        // get specific f table cal correction  
        f_corr = (f_cal_corr == NULL) ? 1.0
                : f_cal_corr[ib][l1rec->detnum][l1rec->mside];

        if (strcmp(bandType[ib], "TEB") == 0) {

            for (ip = 0; ip < num_pixels; ip++) {
                ipb = ip * NBANDSIR + iteb;
                l1rec->Ltir[ipb] = 0;

                l1bptrs_scan_idx = l1rec->detnum * 3200 + ip + extract_pixel_start;

                if (l1bptrs[ib][l1bptrs_scan_idx] != -32767) {
                    l1rec->Ltir[ipb] = l1bptrs[ib][l1bptrs_scan_idx] / 10.0;

                    // Apply F-factor 
                    l1rec->Ltir[ipb] *= f_corr;
                }

            }
            iteb++;

        } else if (strcmp(bandType[ib], "CIR") == 0) {

            for (ip = 0; ip < num_pixels; ip++) {

                l1bptrs_scan_idx = l1rec->detnum * 3200 + ip + extract_pixel_start;

                if (l1bptrs[ib][l1bptrs_scan_idx] != -32767) {
                    l1rec->rho_cirrus[ip] = l1bptrs[ib][l1bptrs_scan_idx];

                    // Normalize reflectance by solar zenith angle 
                    l1rec->rho_cirrus[ip] /= cos(l1rec->solz[ip] / RADEG);

                    // Apply F-factor 
                    l1rec->rho_cirrus[ip] *= f_corr;
                }

            }

        } else if (strcmp(bandType[ib], "RSB") == 0) {

            l1rec->Fo[irsb] = Fobar[irsb] * l1rec->fsol;

            // copy to Lt record.
            for (ip = 0; ip < num_pixels; ip++) {
                ipb = ip * l1rec->l1file->nbands + irsb;

                l1bptrs_scan_idx = l1rec->detnum * 3200 + ip + extract_pixel_start;

                if (l1bptrs[ib][l1bptrs_scan_idx] != -32767) {
                    l1rec->Lt[ipb] = l1bptrs[ib][l1bptrs_scan_idx];

                    // convert from reflectance to radiance 
                    l1rec->Lt[ipb] *= l1rec->Fo[irsb] / PI;

                    // Apply F-factor 
                    l1rec->Lt[ipb] *= f_corr;
                }

            }
            irsb++;
        } // if RSB

    } // for ib
    */

    // Skipping for now to avoid "No brightness temperature conversion provided for this sensor" [SBL]
    //radiance2bt(l1rec, -1); // calculate brightness temperature

    // Bowtie stuff is specific to VIIRS, so skip it.
    //for (ip = 0; ip < num_pixels; ip++) {
    //    flag_bowtie_deleted(l1rec, ip, extract_pixel_start);
    //}

    prevScan = scan;
    
    return (LIFE_IS_GOOD);
}


/**
 * Close L1A file, GEO file, and free memory
 * @param file
 * @return 
 */
int closel1a_hawkeye(filehandle *file) {
    int status;

    printf("Closing hawkeye l1a file\n");
    status = nc_close(file->sd_id);
    check_err(status, __LINE__, __FILE__);

    if (file->geofile[0] != '\0') {
        printf("Closing hawkeye l1a GEO file\n");
        status = nc_close(geoFileId);
        check_err(status, __LINE__, __FILE__);

        // Free memory
        // From readl1a_hawkeye
        if (tmpShort) free(tmpShort);
        if (tmpByte) free(tmpByte);
        if (CCD_temperatures_this_scan) free(CCD_temperatures_this_scan);
        if (dn) free2d_double(dn);                         
        if (Lt) free2d_double(Lt);
        // From read_cal_hawkeye
        if (temperature_ref) free(temperature_ref);
        if (time_ref) free(time_ref);
        if (K1) free(K1); 
        if (K2) free3d_all_double(K2, num_bands, num_pixels, cal_num_dates);
        if (K3) free2d_double(K3);
        if (K4) free2d_float(K4);
        if (K5) free2d_double(K5);       
        // From openl1a_hawkeye
        if (CCD_temperatures) free2d_double(CCD_temperatures);
        if (CCD_temperatures_cleansed) free2d_double(CCD_temperatures_cleansed); 
        if (tlm_delta_time_ms) free(tlm_delta_time_ms);
        if (tlm_delta_time_ms_cleansed) free(tlm_delta_time_ms_cleansed);  
        if (dn_varid) free(dn_varid);
    
    }

    return (LIFE_IS_GOOD);
}


/**
 * Perform quality control on the CCD_temperatures[num_tlm_blocks][num_ccd_temps] 
 * and its time base tlm_delta_time_ms[num_tlm_blocks], creating cleansed versions of both of these arrays.
 * Cleansing includes the following rules:
 *    -If one or more CCD_temperatures for a given tlm block are FILL, replace them with 
 *     the mean value of the non-FILL temperatures.
 *    -Identify possibly bad (extreme) values for CCD temperatures e.g. if one 
       of the 4 CCD thermistors has a value much different from the others.
 *    -If ALL the CCD temperatures for a given tlm block are FILL, skip this tlm block. 
 *    -If there are no good tlm blocks, manufacture 2 good rows using a default value for CCD temperature.
 */
void qc_hawkeye_CCD_T() {
    
    // Misc declarations
    size_t tlm_block_num, tlm_block_cleansed_num, ccd_temp_num;
    bool *tlm_block_good;
    double *CCD_temperatures_Dim2, *nan_mean_weights;
    double CCD_temperature_max_deviation = 10.0;                    // deviation from median value
    double *CCD_temperature_sorted; 
    size_t *CCD_temperature_sort_index;
    double CCD_temperature_median;
    
    // Allocate memory to arrays
    tlm_block_good = (bool *) calloc(num_tlm_blocks, sizeof (bool));
    CCD_temperatures_Dim2 = (double *) calloc(num_ccd_temps, sizeof (double));
    nan_mean_weights = (double *) calloc(num_ccd_temps, sizeof (double));
    CCD_temperature_sorted = (double *) calloc(num_ccd_temps, sizeof (double));
    CCD_temperature_sort_index = (size_t *) calloc(num_ccd_temps, sizeof (size_t));
    
    // To test, set some or all CCD temperature values to FILL or bad values
    //for (tlm_block_num=7; tlm_block_num<8; tlm_block_num++) {
    //    for (ccd_temp_num=2; ccd_temp_num<3; ccd_temp_num++) {
    //        CCD_temperatures[tlm_block_num][ccd_temp_num] = fill_value_in_CCD_T;
    //        CCD_temperatures[tlm_block_num][ccd_temp_num] = 75.0;
    //    }
    //}
    
    // If a given row has all FILL values, skip that row.
    for (tlm_block_num=0; tlm_block_num<num_tlm_blocks; tlm_block_num++) {                          // Dim1
        // Apply rules. If at least one of the CCD_temperatures in this tlm_block is good, do not skip this tlm_block.
        tlm_block_good[tlm_block_num] = false;                                  //Assume worst case
        for (ccd_temp_num=0; ccd_temp_num<num_ccd_temps; ccd_temp_num++) {
            if (CCD_temperatures[tlm_block_num][ccd_temp_num] != fill_value_in_CCD_T) {
                tlm_block_good[tlm_block_num] = true;                           //At least one good value, so don't skip this row
                //printf("tlm_block_num %d is good\n", (int)tlm_block_num);
                break;
            }
        }
    }
    
    // If one or more CCD_temperatures for a given tlm block are FILL, replace them with the mean value of the 
    // non-FILL temperatures.
    // First, set nan_mean_weights
    for (ccd_temp_num=0; ccd_temp_num<num_ccd_temps; ccd_temp_num++) {
        nan_mean_weights[ccd_temp_num] = 1.0;
    }
    tlm_block_cleansed_num = 0;
    for (tlm_block_num=0; tlm_block_num<num_tlm_blocks; tlm_block_num++) {                          // Dim1
        if (tlm_block_good[tlm_block_num]) {
            // This is a good row, so populate CCD_temperatures_cleansed, tlm_delta_time_ms_cleansed and increment 
            // tlm_block_cleansed_num
            tlm_delta_time_ms_cleansed[tlm_block_cleansed_num] = tlm_delta_time_ms[tlm_block_num];
            // If any CCD_temperatures are FILL, replace them with the mean value of CCD_temperatures_Dim2
            // We have already removed any rows that have ALL values = FILL, so this should work.
            for (ccd_temp_num=0; ccd_temp_num<num_ccd_temps; ccd_temp_num++) {                      // Dim2
                // First pass to build an array just for this row i.e. CCD_temperatures_Dim2
                CCD_temperatures_Dim2[ccd_temp_num] = CCD_temperatures[tlm_block_num][ccd_temp_num]; 
            }
            for (ccd_temp_num=0; ccd_temp_num<num_ccd_temps; ccd_temp_num++) {                      // Dim2
                // Second pass to replace FILL values with the mean for this row
                if (CCD_temperatures_Dim2[ccd_temp_num] == fill_value_in_CCD_T) {
                    CCD_temperatures_cleansed[tlm_block_cleansed_num][ccd_temp_num] = nan_wmean(nan_mean_weights,CCD_temperatures_Dim2, 
                                                                                                num_ccd_temps, (double)fill_value_in_CCD_T);
                    //printf("CCD_temperatures_cleansed[%d][%d] was FILL, now set to the mean for this row = %f\n",
                    //        (int)tlm_block_cleansed_num,(int)ccd_temp_num,CCD_temperatures_cleansed[tlm_block_cleansed_num][ccd_temp_num]);
                } else {
                    CCD_temperatures_cleansed[tlm_block_cleansed_num][ccd_temp_num] = CCD_temperatures_Dim2[ccd_temp_num];
                }
            }
            tlm_block_cleansed_num = tlm_block_cleansed_num + 1;
        } else {
            // This is not a good row i.e. all FILL, so do nothing. This row is skipped.
        }
    }
    
    // At this point, FILL values have been replaced. Identify possibly bad (extreme) values for CCD temperatures e.g. if one 
    // of the 4 thermistors has a value much different from the others.
    for (tlm_block_num=0; tlm_block_num<num_tlm_blocks; tlm_block_num++) {                          // Dim1
        for (ccd_temp_num=0; ccd_temp_num<num_ccd_temps; ccd_temp_num++) {                          // Dim2
            // First pass to build an array just for this row i.e. CCD_temperatures_Dim2
            CCD_temperatures_Dim2[ccd_temp_num] = CCD_temperatures[tlm_block_num][ccd_temp_num];         
        }
        // Sort before calculating median value
        gsl_sort_index(CCD_temperature_sort_index, CCD_temperatures_Dim2, 1, num_ccd_temps);
        for (ccd_temp_num=0; ccd_temp_num<num_ccd_temps; ccd_temp_num++) {
            CCD_temperature_sorted[ccd_temp_num] = CCD_temperatures_Dim2[CCD_temperature_sort_index[ccd_temp_num]];
        }
        // Calculate the median value
        CCD_temperature_median = gsl_stats_median_from_sorted_data(CCD_temperature_sorted, 1,num_ccd_temps);
        // Loop through again to see if any values are too far from median value.
        for (ccd_temp_num=0; ccd_temp_num<num_ccd_temps; ccd_temp_num++) {                          // Dim2
            if (fabs(CCD_temperatures_Dim2[ccd_temp_num] - CCD_temperature_median) > CCD_temperature_max_deviation) {
                fprintf(stderr, "-W- In tlm block %d, CCD %d records a temperature of %f, more than %f degrees different than the median value of %f\n",
                        tlm_block_num, ccd_temp_num, CCD_temperatures_Dim2[ccd_temp_num], CCD_temperature_max_deviation, CCD_temperature_median);
            }
        }
    }
    
    
    // Handle special case where there are no good tlm blocks by setting CCD temperatures to a default value 
    if (tlm_block_cleansed_num < 2){
        // There were NO good tlm blocks
        for (tlm_block_num=0; tlm_block_num<num_tlm_blocks; tlm_block_num++) {                          // Dim1
            tlm_delta_time_ms_cleansed[tlm_block_num] = tlm_delta_time_ms[tlm_block_num];
            for (ccd_temp_num=0; ccd_temp_num<num_ccd_temps; ccd_temp_num++) {
                CCD_temperatures_cleansed[tlm_block_num][ccd_temp_num] = CCD_temperature_default;
            }
        }
        num_tlm_blocks_cleansed = num_tlm_blocks;
        fprintf(stderr, "-W- There were fewer than 2 good tlm_blocks of CCD_temperatures, so using the default value of %f\n", 
                CCD_temperature_default);
    } else {
        // There WAS one or more good tlm blocks
        num_tlm_blocks_cleansed = tlm_block_cleansed_num;
    }
    //printf("At delta_time = %f, CCD_temperatures_cleansed[0][0] = %f\n", 
    //        tlm_delta_time_ms_cleansed[0], CCD_temperatures_cleansed[0][0]);
    //printf("At delta_time = %f, CCD_temperatures_cleansed[num_tlm_blocks_cleansed-1][num_ccd_temps-1] = %f\n", 
    //        tlm_delta_time_ms_cleansed[num_tlm_blocks_cleansed-1], CCD_temperatures_cleansed[num_tlm_blocks_cleansed-1][num_ccd_temps-1]);
    
    // Clean up
    free(tlm_block_good);
    free(CCD_temperatures_Dim2);
    free(nan_mean_weights);
    free(CCD_temperature_sorted); 
    free(CCD_temperature_sort_index);
}


/**
 * Get calibration coefficients (for instrument calibration) from the specified calibration file
 * and store them in static arrays. These arrays K1-K5, time_ref, and temperature_ref have already been 
 * declared static. Here in read_cal_hawkeye, we determine their dimensions, allocate memory to them,
 * and populate them.
 * @param cal_path
 */
void read_cal_hawkeye(char *cal_path) {
    
    // Declarations for reading cal file
    int status, ncid_CAL, varid, dimid; 
    size_t start[3], count[3];
    
    // Misc declarations
    int band_num, pixel_num, cal_temperature_num, cal_date_num;
    size_t cal_num_bands = 0;
    size_t cal_num_pixels = 0;
    
    // Open LUT for calibrations
    status = nc_open(cal_path, NC_NOWRITE, &ncid_CAL);
    if (status == FAIL) {
        fprintf(stderr, "-E- %s line %d: nc_open(%s) failed.\n",
                __FILE__, __LINE__, cal_path);
        exit(EXIT_FAILURE);
    } 
    
    // Get dims from cal file
    // num_bands
    status = nc_inq_dimid(ncid_CAL, "phony_dim_0", &dimid);
    if (status != NC_NOERR) {
        fprintf(stderr, "-E- Error reading cal_num_bands.\n");
        exit(EXIT_FAILURE);
    };
    nc_inq_dimlen(ncid_CAL, dimid, &cal_num_bands);
    // num_pixels
    status = nc_inq_dimid(ncid_CAL, "phony_dim_1", &dimid);
    if (status != NC_NOERR) {
        fprintf(stderr, "-E- Error reading cal_num_pixels.\n");
        exit(EXIT_FAILURE);
    };
    nc_inq_dimlen(ncid_CAL, dimid, &cal_num_pixels);
    // cal_num_dates
    status = nc_inq_dimid(ncid_CAL, "phony_dim_2", &dimid);
    if (status != NC_NOERR) {
        fprintf(stderr, "-E- Error reading cal_num_dates.\n");
        exit(EXIT_FAILURE);
    };
    nc_inq_dimlen(ncid_CAL, dimid, &cal_num_dates);
    // cal_num_temperatures
    status = nc_inq_dimid(ncid_CAL, "phony_dim_3", &dimid);
    if (status != NC_NOERR) {
        fprintf(stderr, "-E- Error reading cal_num_temperatures.\n");
        exit(EXIT_FAILURE);
    };
    nc_inq_dimlen(ncid_CAL, dimid, &cal_num_temperatures);
    
    // Check to see if cal file dims = L1A file dims
    if (cal_num_bands != num_bands) {
        fprintf(stderr, "-E- num_bands is cal file not equal to num_bands in L1A file\n");
        exit(EXIT_FAILURE);
    }
    if (cal_num_pixels != num_pixels) {
        fprintf(stderr, "-E- num_pixels is cal file not equal to num_pixels in L1A file\n");
        exit(EXIT_FAILURE);
    }
    
    
    // Declare 3D array used for calibration coefficients
    // I tried declaring double ***K2, but call to nc_get_vara_double failed
    // So, the following declaration works for the call to nc_get_vara_double.
    double K2_temp[num_bands][num_pixels][cal_num_dates];                                    //Temporal correction, BEFORE interpolation 
    
    // Allocate space for arrays used for calibration coefficients 
    // 1D arrays
    time_ref = (double *) calloc(cal_num_dates, sizeof (double));
    temperature_ref = (double *) calloc(cal_num_temperatures, sizeof (double));
    K1  = (double *) calloc(num_bands, sizeof (double));        
    // 2D arrays
    // When calling alloc2d, switch order of dims i.e. it wants (w,h), where h is the first dim
    K3  = alloc2d_double(cal_num_temperatures,num_bands);
    K4  = alloc2d_float(num_pixels,num_bands);
    K5  = alloc2d_double(num_pixels,num_bands);
    // 3D arrays
    K2 = alloc3d_double(num_bands, num_pixels, cal_num_dates);
    
    // Initialize vars for reading arrays from cal file
    start[0] = 0;
    start[1] = 0;
    start[2] = 0;
    
    // Get variables from calibration LUT
    // Get K1
    status = nc_inq_varid (ncid_CAL, "K1", &varid);
    if (status != NC_NOERR) {
        fprintf(stderr, "-E- %s line %d: nc_inq_varid failed for K1.\n",
               __FILE__, __LINE__);
        exit(EXIT_FAILURE);
    } else {      
        count[0] = num_bands;
        count[1] = 0;
        count[2] = 0;
        status = nc_get_vara_double(ncid_CAL, varid, start, count, K1); 
    }
    // Get K2. 
    status = nc_inq_varid (ncid_CAL, "K2", &varid);
    if (status != NC_NOERR) {
        fprintf(stderr, "-E- %s line %d: nc_inq_varid failed for K2.\n",
               __FILE__, __LINE__);
        exit(EXIT_FAILURE);
    } else { 
        count[0] = num_bands;
        count[1] = num_pixels;
        count[2] = cal_num_dates;
        status = nc_get_vara_double(ncid_CAL, varid, start, count, &K2_temp[0][0][0]);         // Auto convert from float to double
    }
    // Stuff K2_temp into K2
    for (band_num=0; band_num<num_bands; band_num++) {
        for (pixel_num=0; pixel_num<num_pixels; pixel_num++) {
            for (cal_date_num=0; cal_date_num<cal_num_dates; cal_date_num++) {
                K2[band_num][pixel_num][cal_date_num] = K2_temp[band_num][pixel_num][cal_date_num];
            }
        }
    }
    
    // Get K2t
    status = nc_inq_varid (ncid_CAL, "K2t", &varid);
    if (status != NC_NOERR) {
        fprintf(stderr, "-E- %s line %d: nc_inq_varid failed for K2t.\n",
               __FILE__, __LINE__);
        exit(EXIT_FAILURE);
    } else {      
        count[0] = cal_num_dates;
        count[1] = 0;
        count[2] = 0;
        status = nc_get_vara_double(ncid_CAL, varid, start, count, time_ref); 
    }
    // Get K3
    status = nc_inq_varid (ncid_CAL, "K3", &varid);
    if (status != NC_NOERR) {
       fprintf(stderr, "-E- %s line %d: nc_inq_varid failed for K3.\n",
               __FILE__, __LINE__);
       exit(EXIT_FAILURE);
    } else { 
       count[0] = num_bands;
       count[1] = cal_num_temperatures;
       count[2] = 0;
       status = nc_get_vara_double(ncid_CAL, varid, start, count, &K3[0][0]);          // Auto convert from float to double
    }
    // Get K3T
    status = nc_inq_varid (ncid_CAL, "K3T", &varid);
    if (status != NC_NOERR) {
        fprintf(stderr, "-E- %s line %d: nc_inq_varid failed for K3T.\n",
               __FILE__, __LINE__);
        exit(EXIT_FAILURE);
    } else {      
        count[0] = cal_num_temperatures;
        count[1] = 0;
        count[2] = 0;
        status = nc_get_vara_double(ncid_CAL, varid, start, count, temperature_ref);    // Auto convert from float to double
    }
    // Get K4
    status = nc_inq_varid (ncid_CAL, "K4", &varid);
    if (status != NC_NOERR) {
       fprintf(stderr, "-E- %s line %d: nc_inq_varid failed for K4.\n",
               __FILE__, __LINE__);
       exit(EXIT_FAILURE);
    } else { 
       count[0] = num_bands;
       count[1] = num_pixels;
       count[2] = 0;
       status = nc_get_vara_float(ncid_CAL, varid, start, count, &K4[0][0]);         
    }
    // Get K5
    status = nc_inq_varid (ncid_CAL, "K5", &varid);
    if (status != NC_NOERR) {
       fprintf(stderr, "-E- %s line %d: nc_inq_varid failed for K5.\n",
               __FILE__, __LINE__);
       exit(EXIT_FAILURE);
    } else { 
       count[0] = num_bands;
       count[1] = num_pixels;
       count[2] = 0;
       status = nc_get_vara_double(ncid_CAL, varid, start, count, &K5[0][0]);          // Auto convert from float to double
    }
       
    // Cleanup
    // Close files
    status = nc_close(ncid_CAL);
    
       
}


/**
 * The array CCD_temperatures_cleansed has dimensions num_tlm_blocks_cleansed x num_ccd_temps. 
 * Given the input value of delta_time_ms_i, we interpolate on the first dimension, generating
 * the array CCD_temperatures_this_scan, having dimension 1 x num_ccd_temps
 * @param delta_time_ms_i
 */
void interp_hawkeye_CCD_T(double delta_time_ms_i) {
    
    // Misc declarations
    size_t tlm_block_num, ccd_temp_num; 
    double *CCD_temperatures_Dim1;
    
    // Allocate memory for arrays
    CCD_temperatures_Dim1 = (double *) calloc(num_tlm_blocks_cleansed, sizeof (double));
   
    gsl_interp_accel *acc_delta_time = gsl_interp_accel_alloc();
    gsl_spline *spline_delta_time = gsl_spline_alloc(gsl_interp_linear, num_tlm_blocks_cleansed);
    
    // First, before interpolating, make sure delta_time_ms_i is not outside the tlm_delta_time_ms_cleansed array
    //printf("Before prep_for_interp_double, delta_time_ms_i = %f\n", delta_time_ms_i);
    prep_for_interp_double(&delta_time_ms_i, tlm_delta_time_ms_cleansed, num_tlm_blocks_cleansed);
    //printf("After prep_for_interp_double, delta_time_ms_i = %f\n", delta_time_ms_i);
    for (ccd_temp_num=0; ccd_temp_num<num_ccd_temps; ccd_temp_num++) {                 // Dim2
        for (tlm_block_num=0; tlm_block_num<num_tlm_blocks_cleansed; tlm_block_num++) {         // Dim1
            // Temporarily store a column in CCD_temperatures_Dim1
            CCD_temperatures_Dim1[tlm_block_num] = CCD_temperatures_cleansed[tlm_block_num][ccd_temp_num];       
        }
        // Interpolate this column in time to get CCD_temperatures_i[ccd_temp_num]
        gsl_spline_init(spline_delta_time, tlm_delta_time_ms_cleansed, &CCD_temperatures_Dim1[0], num_tlm_blocks_cleansed);
        CCD_temperatures_this_scan[ccd_temp_num] = gsl_spline_eval(spline_delta_time, delta_time_ms_i, acc_delta_time);
    }
    
    // Clean up
    free(CCD_temperatures_Dim1);
    gsl_spline_free(spline_delta_time);
    gsl_interp_accel_free(acc_delta_time);
}


/**
 * Apply instrument calibration, converting dn into Lt.
 * @param current_julian_date
 * @return 
 */
int calibrate_hawkeye(double current_julian_date){
   
    // Declare demo vars   
    double this_dn;
    double current_temp_C;

    // Misc declarations
    int band_num, pixel_num, cal_temperature_num, cal_date_num, ccd_temp_num;
    double cal_coeff_1_4;
        
    // Declarations for interpolation
    double **K2i;                            //Temporal correction AFTER interpolation, [num_bands][num_pixels]
    double *K3i;                             //Temperature correction AFTER interpolation,  [num_bands] 
    gsl_interp_accel *acc_T = gsl_interp_accel_alloc();
    gsl_interp_accel *acc_date = gsl_interp_accel_alloc();
    gsl_spline *spline_T = gsl_spline_alloc(gsl_interp_linear, cal_num_temperatures);
    gsl_spline *spline_date = gsl_spline_alloc(gsl_interp_linear, cal_num_dates);
    
    // Allocate memory for arrays
    // When calling alloc2d, switch order of dims i.e. it wants (w,h), where h is the first dim
    K2i = alloc2d_double(num_pixels,num_bands);
    K3i = (double *) calloc(num_bands, sizeof (double));
       
    // Interpolate in time for K2, picking some time just for demo.
    //printf("Before prep_for_interp_double, current_julian_date = %f\n", current_julian_date);
    prep_for_interp_double(&current_julian_date, time_ref, cal_num_dates);
    //printf("After prep_for_interp_double, current_julian_date = %f\n", current_julian_date);
    for (band_num=0; band_num<num_bands; band_num++) {
        for (pixel_num=0; pixel_num<num_pixels; pixel_num++) {
            gsl_spline_init(spline_date, time_ref, &K2[band_num][pixel_num][0], cal_num_dates);
            K2i[band_num][pixel_num] = gsl_spline_eval(spline_date, current_julian_date, acc_date);
        }
    }

    // Interpolate in temperature for K3
    //printf("\n");
    for (band_num=0; band_num<num_bands; band_num++) {
        // Given the band_num, pick the correct CCD_temperature
        ccd_temp_num = (int)(band_num/2);
        //printf("For band_num = %d, getting temp from CCD_temperatures_this_scan[%d]\n", band_num, ccd_temp_num);
        current_temp_C = CCD_temperatures_this_scan[ccd_temp_num];
        // Make sure this temperature is in bounds
        //printf("Before prep_for_interp_double, current_temp_C = %f\n", current_temp_C);
        prep_for_interp_double(&current_temp_C, temperature_ref, cal_num_temperatures);
        //printf("After prep_for_interp_double, current_temp_C = %f\n", current_temp_C);
        // Interpolate K3, given this temperature
        gsl_spline_init(spline_T, temperature_ref, &K3[band_num][0], cal_num_temperatures);
        K3i[band_num] = gsl_spline_eval(spline_T, current_temp_C, acc_T);
    }

    // Display K values
    // time_ref
    //printf("\n");
    for (cal_date_num=0; cal_date_num<cal_num_dates; cal_date_num++) {
       //printf("time_ref[%d] = %f\n",cal_date_num, time_ref[cal_date_num]);
    }
    //printf("\n");
    // temperature_ref
    for (cal_temperature_num=0; cal_temperature_num<cal_num_temperatures; cal_temperature_num++) {
       //printf("temperature_ref[%d] = %f\n",cal_temperature_num, temperature_ref[cal_temperature_num]);
    }
    //printf("\n");
    // K1
    for (band_num=0; band_num<num_bands; band_num++) {
        //printf("K1[%d] = %f\n",band_num, K1[band_num]);
    }
    //printf("\n");
    // Subset of K2, K2i.
    band_num = num_bands-1;
    pixel_num = num_pixels-1;
    for (cal_date_num=0; cal_date_num<cal_num_dates; cal_date_num++) {
        //printf("K2[%d][%d][%d] = %f\n",band_num,pixel_num,cal_date_num, K2[band_num][pixel_num][cal_date_num]);
    }
    //printf("At t=%f, K2i[%d][%d] = %f\n",current_julian_date,band_num,pixel_num,K2i[band_num][pixel_num]);
    //printf("\n");
    // K3
    for (band_num=0; band_num<num_bands; band_num++) {
        for (cal_temperature_num=0; cal_temperature_num<cal_num_temperatures; cal_temperature_num++) {
            //printf("K3[%d][%d] = %f\n",band_num,cal_temperature_num, K3[band_num][cal_temperature_num]);
        }
        //printf("At T=%f, K3i[%d] = %f\n",current_temp_C,band_num,K3i[band_num]);
        //printf("\n");
    }
    // Subset of K4,K5
    pixel_num = num_pixels-1;
    for (band_num=0; band_num<num_bands; band_num++) {
        //printf("K4[%d][%d] = %f\n",band_num,pixel_num, K4[band_num][pixel_num]);
    }
    //printf("\n");
    pixel_num = num_pixels-1;
    for (band_num=0; band_num<num_bands; band_num++) {
        //printf("K5[%d][%d] = %f\n",band_num,pixel_num, K5[band_num][pixel_num]);
    }

    // Apply cal to sample scan
    for (band_num=0; band_num<num_bands; band_num++) {
        for (pixel_num=0; pixel_num<num_pixels; pixel_num++) {
            this_dn = dn[band_num][pixel_num];
            cal_coeff_1_4 = (K1[band_num])*(K2i[band_num][pixel_num])*(K3i[band_num])*(K4[band_num][pixel_num]);
            Lt[band_num][pixel_num] = this_dn*cal_coeff_1_4*(1 + (K5[band_num][pixel_num])*this_dn);
        }
    }
    
    // Clean up e.g. free memory
    free2d_double(K2i);
    free(K3i);
    gsl_spline_free(spline_T);
    gsl_interp_accel_free(acc_T);
    gsl_spline_free(spline_date);
    gsl_interp_accel_free(acc_date);
    
    // Done
    return(0);
}



double nan_wmean(double *weights, double *data, size_t n, double fill_value) {
    size_t i, fill_count;
    double wmean;
    //Set weight to zero if corresponding data element is FILL
    fill_count = 0;
    for (i = 0; i < n; i++) {
        if (data[i] == fill_value) {
            weights[i] = 0;
            fill_count++;
        }
    }
    //Call gsl routing for weighted mean
    if (fill_count < n) {
       wmean = gsl_stats_wmean(weights, 1, data, 1, n);
    } else {
       wmean = fill_value;
    }
    return wmean;

}


/**
 * If xi is outside the range of x, set it to the nearest neighbor.
 * @param xi
 * @param x
 * @param N
 */
void prep_for_interp_double(double *xi, double *x, int N) {
    // Adjust xi to be within range of x 
    
    double xmin, xmax;
    
    gsl_stats_minmax (&xmin, &xmax, x, 1, N);
    if (*xi < xmin) {
        *xi = xmin;
    }
    if (*xi > xmax) {
        *xi = xmax;
    }
}


double ***alloc3d_double(size_t m, size_t n, size_t o) {
    size_t i, j;
    double ***p;
    p = malloc(m * sizeof (double *));

    for (i = 0; i < m; i++) {
        p[ i ] = malloc(n * sizeof (double **));
        for (j = 0; j < n; j++) {
            p[ i ][ j ] = malloc(o * sizeof (double ***));
        }
    }

    return (p);
}


void free3d_all_double(double ***p, size_t m, size_t n, size_t o) {
    // Use the same dimensions as was passed to alloc3d_double
    size_t i, j;
    
    for (i = 0; i < m; i++) {
        for (j = 0; j < n; j++) {
            free(p[ i ][ j ]);
        }
        free(p[ i ]);
    }
    
    free(p);
}




