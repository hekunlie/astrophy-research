#ifndef GGLEN_FUNCTIONS_H
#define GGLEN_FUNCTIONS_H

#define MAX_JACK 2000
#define MY_FLOAT float
#define MAX_FORE_EXPO 100
#define MAX_BACK_EXPO 
struct GG_data_info
{
    int jack_num;
    int theta_bin_num;
    
    ////////////////////  for the SYM_PDF method of Fourier_Quad  ///////////////////////
    // bin num for G1,G2
    int mg_bin_num;
    // total source count in the PDF, maintained in the master thread
    // the first half belongs to g_t,
    // the rest belongs to g_x 
    // (from which the signals should be consistent with 0, systematic check)   
    double *total_chi_count[MAX_JACK];        
    // source count of each pair-calculation in each worker thread
    // it will be passed to the master thread when it finishes one calculation job
    double *chi_count;

    ///////////////////  the position informs of each foreground exposure //////////////
    MY_FLOAT *foreground_infoms;

    MY_FLOAT *foreground_expo_data;
    int foreground_expo_num;
    int foreground_data_col;
    int foreground_data_row;

    int foreground_ra_col;
    int foreground_dec_col;
    int foreground_z_col;



};


#endif
