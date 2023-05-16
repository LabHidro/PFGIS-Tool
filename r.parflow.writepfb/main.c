
/****************************************************************************
 *
 * MODULE:    r.parflow.writepfb
 * 
 * AUTHOR:    Tomas Carlotto ----------------- thomas.carl@hotmail.com
 *                               
 *               
 * PURPOSE:      Write PFB files for Parflow hydrologic model
 * 
 * COPYRIGHT:    (C) 2022 by the GRASS Development Team
 *
 *               This program is free software under the GNU General Public
 *   	    	 License (>=v2). Read the file COPYING that comes with GRASS
 *   	    	 for details.
 *
 ****************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stddef.h>
#include <stdint.h>

#include <grass/gis.h>
#include <grass/raster.h>
#include <grass/glocale.h>
#include <grass/gmath.h>

#include <grass/config.h>
#include <grass/spawn.h>

#include "functions.h"

typedef double ValueType;

//  Feature

typedef struct{
    
//-------------INPUTS FOR SIMULATIATION ('Simulation Inputs guisection')----------------
    
    struct Option *ValuesMap;        
    struct Option *pfbOutput;
    struct Option *PathOutput;
    

} 

paramType;
paramType param, outparam;		//Parameters 

void set_parameters(){
//---------------------------------------------------------------------------------------------

    // ============================================================
    // ================== Domain definition =======================
        
    param.ValuesMap = G_define_standard_option(G_OPT_R_INPUT);
    param.ValuesMap->key = "values";
    param.ValuesMap->label = _("Name of input raster map");
    param.ValuesMap->required = YES;
    param.ValuesMap->guisection = _("Required");
    
    param.pfbOutput= G_define_standard_option(G_OPT_R_OUTPUT);
    param.pfbOutput->key = "out_pfb";
    param.pfbOutput->label = _("Name for the PFB output file (.pfb)");
    param.pfbOutput->description = _(" ");
    param.pfbOutput->required = YES;
    param.pfbOutput->guisection = _("Required");
                            
    outparam.PathOutput = G_define_standard_option(G_OPT_M_DIR);//G_OPT_F_OUTPUT);
    outparam.PathOutput->key = "out_dir";
    outparam.PathOutput->type = TYPE_STRING;
    outparam.PathOutput->description = _("Name for output directory");
    outparam.PathOutput->required = YES;
    outparam.PathOutput->guisection = _("Required");    
        
  }

int main(int argc, char *argv[])
{
    struct Cell_head window; 
    int nrows, ncols;
   // int Temp_value;   
   
    struct History history;	// holds meta-data (title, comments,..) 
  
    struct GModule *module;	 // GRASS module for parsing arguments

    // Initialize GIS environment
    G_gisinit(argv[0]);		 // reads grass env, stores program name to G_program_name()

    // Header of window module
    module = G_define_module();
    G_add_keyword(_("ParFlow"));
    G_add_keyword(_("Hydrologic model"));
    module->description = _("Preprocessing tool to write PFB files for ParFlow hydrologic model");
    
    // set parameters for user (user window)
    set_parameters();

    // options and flags parser 
    if (G_parser(argc, argv))
	exit(EXIT_FAILURE);

//---------------------------------------------------------------------------------------------
// GET MAP INFORMATIONS ( nº of rows and cols and map resolution)------------------------------
//---------------------------------------------------------------------------------------------

// G_message(_("START"));
   
nrows = Rast_window_rows();
ncols = Rast_window_cols();    
    
G_get_window(&window);
double yres = (double)window.ns_res;
double xres = (double)window.ew_res;
        
double xcorner = (double)window.west;
double ycorner = (double)window.south;             

RASTER_MAP_TYPE data_type_input_map;

char *name_input_map = param.ValuesMap->answer;        // Name of input raster Map 
char *out_pfb_name = param.pfbOutput->answer;          // pfb file name

data_type_input_map = MAPdatatype(name_input_map);

int i,j,is,in_i;

double **pfvalues1 = (double**)malloc(nrows*sizeof(double*));
for (i = 0; i < nrows; i++){                      
    pfvalues1[i] = (double*)malloc(ncols*sizeof(double));
}

G_message(_("START"));
            
//========= LOAD DEM DATA ==========================
if(data_type_input_map == CELL_TYPE){
    G_message(_("LOAD INPUT MAP"));
    in_i = nrows-1;    

    CELL** pfvalues = (CELL**) G_malloc( nrows*sizeof(CELL*));    
    for (i = 0; i < nrows; i++){                    
         pfvalues[i] = (CELL*) G_malloc( ncols * sizeof(CELL));
    }
    //
    
    pfvalues = open_raster_C_variable(name_input_map);

    for (i=0;i<nrows;i++){
       for (j=0;j<ncols;j++){   
          pfvalues1[i][j] = pfvalues[in_i-i][j];   
          
       }
    }
      close_raster_C_variable(pfvalues);
}
else if(data_type_input_map == FCELL_TYPE){
    G_message(_("LOAD INPUT DATA"));
    in_i = nrows-1;
   
    FCELL** pfvalues = (FCELL**) G_malloc( nrows*sizeof(FCELL*));    
    for (i = 0; i < nrows; i++){                    
         pfvalues[i] = (FCELL*) G_malloc( ncols * sizeof(FCELL));
    }
    
    pfvalues = open_raster_F_variable(name_input_map);

//    printf("%s\n","DEM MAP IS FCELL TYPE");

    for (i=0;i<nrows;i++){
       for (j=0;j<ncols;j++){   
          pfvalues1[i][j] = pfvalues[in_i-i][j];
          
       }
    }
    
    close_raster_F_variable(pfvalues);
}
else if(data_type_input_map == DCELL_TYPE){
    G_message(_("LOAD INPUT DATA"));
    in_i = nrows-1;
    
    DCELL** pfvalues = (DCELL**) G_malloc( nrows*sizeof(DCELL*));    
    for (i = 0; i < nrows; i++){                    
        pfvalues[i] = (DCELL*) G_malloc( ncols * sizeof(DCELL));
    }    
    
    pfvalues = open_raster_D_variable(name_input_map);

    for (i=0;i<nrows;i++){
        for (j=0;j<ncols;j++){   
          pfvalues1[i][j] = pfvalues[in_i-i][j];
        }
    }
    close_raster_D_variable(pfvalues);
}

G_message(_("SAVING PFB file"));
char *out_dir = outparam.PathOutput->answer; // Save Path
char *PathOutPFB;//    
PathOutPFB = concat(out_dir, "/"); 
PathOutPFB = concat(PathOutPFB, out_pfb_name);   
PathOutPFB = concat(PathOutPFB, ".pfb");

pfb_maker(ncols,pfvalues1, nrows,1,xcorner,ycorner,0.0, xres, yres,0.0, PathOutPFB);
                     			
G_message(_("DONE"));
                                                               
G_free(out_dir); 
G_free(PathOutPFB);      
G_free(name_input_map); 
G_free(out_pfb_name);

exit(EXIT_SUCCESS);
  
}

