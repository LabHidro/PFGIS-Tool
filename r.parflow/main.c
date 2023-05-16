
/****************************************************************************
 *
 * MODULE:    r.parflow
 * 
 * AUTHOR:    Tomas Carlotto ----------------- thomas.carl@hotmail.com
 *                               
 *               
 * PURPOSE:      Create the input files for the ParFlow hydrological model
 * 
 * COPYRIGHT:    (C) 2017 by the GRASS Development Team
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
    
    struct Option *MaskMap;
    struct Option *ElevationMap;
    struct Option *BottomMap;
    struct Option *GridFormat;
        
    //struct Option *Slope_x_Map;
    //struct Option *Slope_y_Map;
    
    struct Option *Slope_x_Output;
    struct Option *Slope_y_Output;    
    struct Option *ElevOutput;
    struct Option *PathOutput;
    struct Option *Solid_pfsolOutput;
    struct Option *Solid_vtkOutput;
    struct Option *NZ_layers;    
    struct Option *dz_resolution;
    
    struct Flag *Out_vtk_flag; 
    struct Flag *SurfaceCond_flag;   
    
    struct Flag *SinkFilling_flag;
    struct Option *SinkFillingMethod;
    struct Option *SinkFillingIter;
    struct Option *SinkFillingIncrement;
    struct Option *SinkFillingOut;
    struct Option *SinkLakeMask;
    
    struct Option *SlopeMethod;
    struct Option *SlopeRiver;
    struct Option *SlopeFlowdir;
    struct Option *SlopeMinSlope; //minimum slope 
    struct Option *SlopeRiverSlope;
    struct Flag *Out_InitTCLfile_flag;
    
        
    struct Option *SurfaceCond;
    struct Option *SurfaceCond_OutName;
    
    struct Option *coordsOut;
    
} 

paramType;
paramType param, outparam;		//Parameters 

/* start points */
struct point
{
    int row;
    int col;
    double value;
    struct point *next;
};

void set_parameters(){
//---------------------------------------------------------------------------------------------
//---------------------------------------------------------------------------------------------
//----MODULE INPUTS -------------
//---------------------------------------------------------------------------------------------
  
    // ============================================================
    // ================== Grid definition =========================
    
    param.GridFormat = G_define_option();
    param.GridFormat->key = "grdf";
    param.GridFormat->type = TYPE_STRING;
    param.GridFormat->required = YES;
    param.GridFormat->description = _("Grid format");
    param.GridFormat->options = "Orthogonal Grid,Terrain Following Grid";
    param.GridFormat->answer = "Orthogonal Grid";
    param.GridFormat->guisection = _("Grid definition");
           
    // ============================================================
    // ================== Domain definition =======================
        
    param.MaskMap = G_define_standard_option(G_OPT_R_INPUT);
    param.MaskMap->key = "pfmsk";
    param.MaskMap->label = _("Raster mask name");
    param.MaskMap->description = _("In the mask the values equal to 0 indicate the inactive cells and values equal to 1 indicate the active cells. Values other than 0 can also be used to define patches at the domain boundaries. A value of -7 defines a surface patch (e.g., lakes or rivers)");
    param.MaskMap->required = YES;
    param.MaskMap->guisection = _("Domain definition");
    
    param.ElevationMap = G_define_standard_option(G_OPT_R_INPUT);
    param.ElevationMap->key = "pftop";
    param.ElevationMap->label = _("Upper boundary raster map name ");
    param.ElevationMap->description = _("The upper boundary, can be the digital elevation model.");
    param.ElevationMap->required = YES;
    param.ElevationMap->guisection = _("Domain definition");
    
    param.BottomMap = G_define_standard_option(G_OPT_R_INPUT);
    param.BottomMap->key = "pfbot";
    param.BottomMap->label = _("Bottom boundary raster map name");
    param.BottomMap->description = _("Raster that defines the bottom of the computational domain.");
    param.BottomMap->required = YES;
    param.BottomMap->guisection = _("Domain definition");
    
    param.Solid_pfsolOutput= G_define_standard_option(G_OPT_R_OUTPUT);
    param.Solid_pfsolOutput->key = "out_pfsol";
    param.Solid_pfsolOutput->label = _("Name for the solid output file (.pfsol)");
    param.Solid_pfsolOutput->description = _(" ");
    param.Solid_pfsolOutput->required = YES;
    param.Solid_pfsolOutput->guisection = _("Domain definition");
    
    param.NZ_layers = G_define_option();
    param.NZ_layers->key = "nzl";
    param.NZ_layers->label = _("Number of ParFlow vertical layers (NZ)");
    param.NZ_layers->required = YES;   
    param.NZ_layers->type = TYPE_INTEGER;
    param.NZ_layers->answer = "1";
    param.NZ_layers->guisection = _("Domain definition");
        
    // ============================================================
    // ================ Sink-Filling ==============================
    
    param.SinkFilling_flag = G_define_flag();
    param.SinkFilling_flag->key = 'f';
    param.SinkFilling_flag->label = _("Enable sink filling");
    param.SinkFilling_flag->description = _("When checking this option, the sink cells will be filled to prevent water retention.");
    param.SinkFilling_flag->guisection = _("Sink filling");
    
    param.SinkLakeMask = G_define_standard_option(G_OPT_R_INPUT);
    param.SinkLakeMask->key = "lake_msk";
    param.SinkLakeMask->label = _("Name of lake mask (Optional)");
    param.SinkLakeMask->description = _("The Lake mask will be used to avoid application over the lake area. If there is no lake this field can be neglected..");
    param.SinkLakeMask->required = NO;
    param.SinkLakeMask->guisection = _("Sink filling");
    
    param.SinkFillingMethod = G_define_option();
    param.SinkFillingMethod->key = "method";
    param.SinkFillingMethod->type = TYPE_STRING;
    param.SinkFillingMethod->required = NO;
    param.SinkFillingMethod->label = _("Methods");
    param.SinkFillingMethod->description = _("Pit-fill method: increases the elevation of each sink cell in the domain (adds a small increment) until there are no more sinks or the maximum number of iterations is reached. Moving average method: the elevation of the sink cells is replaced by the average elevation of the adjacent cells. A small increment is also considered.");
    param.SinkFillingMethod->options = "Pit‐fill,Moving average";
    param.SinkFillingMethod->answer = "Pit‐fill";
    param.SinkFillingMethod->guisection = _("Sink filling");
        
    param.SinkFillingIncrement = G_define_option();
    param.SinkFillingIncrement->key = "increment";
    param.SinkFillingIncrement->label = _("Increment for sink cells [meters]");
    param.SinkFillingIncrement->required = NO;   
    param.SinkFillingIncrement->type = TYPE_DOUBLE;
    param.SinkFillingIncrement->answer = "0.015";
    param.SinkFillingIncrement->guisection = _("Sink filling");    
    
    param.SinkFillingIter = G_define_option();
    param.SinkFillingIter->key = "max_iteration";
    param.SinkFillingIter->required = NO;   
    param.SinkFillingIter->label = _("Number of iterations");    
    param.SinkFillingIter->type = TYPE_INTEGER;
    param.SinkFillingIter->answer = "500";
    param.SinkFillingIter->guisection = _("Sink filling");    
    
    param.SinkFillingOut= G_define_standard_option(G_OPT_R_OUTPUT);
    param.SinkFillingOut->key = "dem_filled";
    param.SinkFillingOut->label = _("Name for the output raster map (Optional)");
    param.SinkFillingOut->description = _(" ");
    param.SinkFillingOut->required = NO;
    param.SinkFillingOut->guisection = _("Sink filling");
    
    // =====================================================================
    // ====================== Slope map ====================================
    
    param.SlopeMethod = G_define_option();
    param.SlopeMethod->key = "slope_method";
    param.SlopeMethod->type = TYPE_STRING;
    param.SlopeMethod->required = NO;
    param.SlopeMethod->label = _("Slope methods");
    param.SlopeMethod->description = _("Upwind scheme, Central difference and Global Slope Enforcement (GSE) method. The GSE method provides complete domain drainage, reducing water retention in flat and sink cells.");
//  param.SlopeMethod->options = "Upwind scheme, Central difference, Global slope enforcement";
    param.SlopeMethod->options = "Upwind scheme, Central difference, Global Slope Enforcement (GSE)";
    param.SlopeMethod->answer = "Central difference";
    param.SlopeMethod->guisection = _("Terrain slopes");    
    
    param.Slope_x_Output = G_define_standard_option(G_OPT_R_OUTPUT);
    param.Slope_x_Output->key = "xslope";
    param.Slope_x_Output->label = _("Name for the x-direction slope output file");
    param.Slope_x_Output->description = _(" ");
    param.Slope_x_Output->required = NO;
    param.Slope_x_Output->guisection = _("Terrain slopes");
    
    param.Slope_y_Output= G_define_standard_option(G_OPT_R_OUTPUT);
    param.Slope_y_Output->key = "yslope";
    param.Slope_y_Output->label = _("Name for the y-direction slope output file");
    param.Slope_y_Output->description = _(" ");
    param.Slope_y_Output->required = NO;
    param.Slope_y_Output->guisection = _("Terrain slopes"); 
    
    // =================================================================    
    // =============== Global slope enforcement method ====================
    
    param.SlopeRiver = G_define_standard_option(G_OPT_R_INPUT);
    param.SlopeRiver->key = "network";//stream
    param.SlopeRiver->label = _("Name of the input drainage network raster map (For GSE method only)");
    param.SlopeRiver->description = _("Drainage network raster map obtained with r.watershed module");
    param.SlopeRiver->required = NO;
    param.SlopeRiver->guisection = _("Terrain slopes");
    
    param.SlopeFlowdir = G_define_standard_option(G_OPT_R_INPUT);
    param.SlopeFlowdir->key = "drainage";//flow_dir
    param.SlopeFlowdir->label = _("Name of the input drainage directions raster map (For GSE method only)");
    param.SlopeFlowdir->description = _("Drainage direction raster map obtained with r.watershed module");
    param.SlopeFlowdir->required = NO;
    param.SlopeFlowdir->guisection = _("Terrain slopes");         
    
    param.SlopeRiverSlope = G_define_option();
    param.SlopeRiverSlope->key = "river_slope";
    param.SlopeRiverSlope->label = _("River slope (For GSE method only)");
    param.SlopeRiverSlope->description = _("Average slope of the drainage network");
    param.SlopeRiverSlope->required = NO;   
    param.SlopeRiverSlope->type = TYPE_DOUBLE;
//  param.SlopeRiverSlope->answer = "0.01";
    param.SlopeRiverSlope->guisection = _("Terrain slopes"); 
 
    param.SlopeMinSlope = G_define_option();
    param.SlopeMinSlope->key = "min_slope";
    param.SlopeMinSlope->label = _("Minimum slope (For GSE method only)");
    param.SlopeMinSlope->description = _("Minimum slope threshold");
    param.SlopeMinSlope->required = NO;   
    param.SlopeMinSlope->type = TYPE_DOUBLE;
//  param.SlopeMinSlope->answer = "0.01";
    param.SlopeMinSlope->guisection = _("Terrain slopes");
        
    // ============================================================
    // =============== Surface conditions ================
    
    param.SurfaceCond_flag = G_define_flag();
    param.SurfaceCond_flag->key = 'b';
    param.SurfaceCond_flag->label = _("Create surface boundary condition");
    param.SurfaceCond_flag->description = _("When checking this option, surface boundary condition will be created in .pfb format");
    param.SurfaceCond_flag->guisection = _("Surface condition");
    
    param.SurfaceCond = G_define_standard_option(G_OPT_R_INPUT);
    param.SurfaceCond->key = "surf_cond";
    param.SurfaceCond->label = _("Name of surface boundary condition map");
    param.SurfaceCond->description = _("This raster contains the values used as surface boundary conditions. The values used depend on the purpose of the user (values to be used with patches of type: PressureFile, FluxFile or OverlandFlowPFB).");
    param.SurfaceCond->required = NO;
    param.SurfaceCond->guisection = _("Surface condition");
    
    param.SurfaceCond_OutName = G_define_standard_option(G_OPT_R_INPUT);
    param.SurfaceCond_OutName->key = "cond_name";
    param.SurfaceCond_OutName->label = _("Name for the output surface boundary condition file");
    param.SurfaceCond_OutName->description = _(" ");
    param.SurfaceCond_OutName->required = NO;
    param.SurfaceCond_OutName->guisection = _("Surface condition");           
              
    // ==========================================================================
    // ================== Output definition =====================================
    param.Out_InitTCLfile_flag = G_define_flag();
    param.Out_InitTCLfile_flag->key = 'g';
    param.Out_InitTCLfile_flag->label = _("Create an initial setup file (tcl format) for running the ParFlow model.");
    param.Out_InitTCLfile_flag->description = _("When this option is checked, an initial setup file will be generated to run the ParFlow model using the generated files.");
    param.Out_InitTCLfile_flag->guisection = _("Output");
    
    param.Out_vtk_flag = G_define_flag();
    param.Out_vtk_flag->key = 's';
    param.Out_vtk_flag->label = _("Create output files in .vtk format.");
    param.Out_vtk_flag->description = _("When this option is checked, in addition to files in .pfb format, files in .vtk format will be generated.");
    param.Out_vtk_flag->guisection = _("Output");
         
    outparam.coordsOut = G_define_standard_option(G_OPT_M_COORDS);
    outparam.coordsOut->description = _("Overland flow observation points (Optional)");
    outparam.coordsOut->multiple = YES;
    outparam.coordsOut->required = NO;
    outparam.coordsOut->guisection = _("Output");     
         
    outparam.PathOutput = G_define_standard_option(G_OPT_M_DIR);//G_OPT_F_OUTPUT);
    outparam.PathOutput->key = "out_dir";
    outparam.PathOutput->type = TYPE_STRING;
    outparam.PathOutput->required = YES;
    outparam.PathOutput->description = _("Name for output directory");
    outparam.PathOutput->guisection = _("Output");    
    
  }

int main(int argc, char *argv[])
{
    struct Cell_head window; 
    int nrows, ncols;
    int i,j,ii,jj,k,is,in_i;
    // int Temp_value;
    
    // Observation points structures
    
    struct point *head_obs_pt = NULL;
    struct point *next_obs_pt;
    int obs_row, obs_col,have_points = 0;                 
    int npoints;
    double east, north;
   
    struct History history;	// holds meta-data (title, comments,..) 
  
    struct GModule *module;	 // GRASS module for parsing arguments

    // Initialize GIS environment
    G_gisinit(argv[0]);		 // reads grass env, stores program name to G_program_name()

    // Header of window module
    module = G_define_module();
    G_add_keyword(_("ParFlow"));
    G_add_keyword(_("Hydrological model"));
    module->description = _("Preprocessing tool for the ParFlow hydrological model");
    
    // set parameters for user (user window)
    set_parameters();

    // options and flags parser 
    if (G_parser(argc, argv))
	exit(EXIT_FAILURE);

//---------------------------------------------------------------------------------------------
//---------------------------------------------------------------------------------------------
// -----VERIFY THE INPUTS (MESSAGE ERRORS)-------------------------------------
//---------------------------------------------------------------------------------------------

char* slope_method = param.SlopeMethod->answer;

if(param.SinkFilling_flag->answer){

    if (!param.SinkFillingIncrement->answer){
        G_fatal_error(_("To perform the sink filling, it is necessary to inform the  (<%s>)."),param.SinkFillingIncrement->key);
    }
    if (!param.SinkFillingIter->answer){
        G_fatal_error(_("To perform the sink filling, it is necessary to inform the  (<%s>)."),param.SinkFillingIter->key);
    }

}

if (G_strcasestr(slope_method,"Global Slope Enforcement (GSE)")){

//if (param.GlobalSlopeEnforcement_flag->answer){
  
    if (!param.SlopeRiver->answer){
        G_fatal_error(_("To perform the global slope enforcement method, it is necessary to inform the  (<%s>)."),param.SlopeRiver->key);
    }
    if (!param.SlopeFlowdir->answer){
        G_fatal_error(_("To perform the global slope enforcement method, it is necessary to inform the  (<%s>)."),param.SlopeFlowdir->key);
    }
    
    if (!param.SlopeMinSlope->answer){
        G_fatal_error(_("To perform the global slope enforcement method, it is necessary to inform the  (<%s>)."),param.SlopeMinSlope->key);
    }
    if (!param.SlopeRiverSlope->answer){
        G_fatal_error(_("To perform the global slope enforcement method, it is necessary to inform the  (<%s>)."),param.SlopeRiverSlope->key);
    }
}

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


RASTER_MAP_TYPE data_type_dem_map;
RASTER_MAP_TYPE data_type_mask_map;
RASTER_MAP_TYPE data_type_bottom_map;

char *name_dem_map = param.ElevationMap->answer;        // Digital Elevation Map name 
char *name_mask_map = param.MaskMap->answer;            // mask map name
char *name_bottom_map = param.BottomMap->answer;        // Bottom limit map name 

char *out_dir = outparam.PathOutput->answer; // Save Path

char *out_pfsol_name = param.Solid_pfsolOutput->answer;  // Solid file name

char * solid_name = concat(out_pfsol_name,".pfsol");

char *grid_format = param.GridFormat->answer;            // Grid format

char *SlopeX_file;
char *SlopeY_file;
char *out_SlopeX_name = param.Slope_x_Output->answer;    // Output name (Slope x direction)
char *out_SlopeY_name = param.Slope_y_Output->answer;    // Output name (Slope y direction)

data_type_dem_map = MAPdatatype(name_dem_map);
data_type_mask_map = MAPdatatype(name_mask_map);
data_type_bottom_map = MAPdatatype(name_bottom_map);

char* filling_method = param.SinkFillingMethod->answer; // Sink filling method



double **pfdem1 = (double**)malloc(nrows*sizeof(double*));
double **pfmask1= (double**)malloc(nrows*sizeof(double*));
double **pfbottom1= (double**)malloc(nrows*sizeof(double*));
double **pfSlopex1= (double**)malloc(nrows*sizeof(double*));
double **pfSlopey1= (double**)malloc(nrows*sizeof(double*));
double **pfSlope= (double**)malloc(nrows*sizeof(double*));

for (i = 0; i < nrows; i++){                      
    pfdem1[i] = (double*)malloc(ncols*sizeof(double));
    pfmask1[i] = (double*)malloc(ncols*sizeof(double));
    pfbottom1[i] = (double*)malloc(ncols*sizeof(double));
    pfSlopex1[i] = (double*)malloc(ncols*sizeof(double));
    pfSlopey1[i] = (double*)malloc(ncols*sizeof(double));
    pfSlope[i] = (double*)malloc(ncols*sizeof(double));    
}

G_message(_("MODULE START"));
            
CELL** C_pfdem;
FCELL** F_pfdem;
DCELL** D_pfdem;
//========= LOAD DEM DATA ==========================
if(data_type_dem_map == CELL_TYPE){
//    G_message(_("LOAD DEM DATA"));
    in_i = nrows-1;    
    //CELL** pfdem;
    //
    C_pfdem = (CELL**) G_malloc( nrows*sizeof(CELL*));  
    CELL** pfdem = (CELL**) G_malloc( nrows*sizeof(CELL*));    
    for (i = 0; i < nrows; i++){                    
         C_pfdem[i] = (CELL*) G_malloc( ncols * sizeof(CELL));
         pfdem[i] = (CELL*) G_malloc( ncols * sizeof(CELL));
    }
    //
    
    pfdem = open_raster_C_variable(name_dem_map);
    if (param.SinkFillingOut->answer){
       C_pfdem = open_raster_C_variable(name_dem_map);
    }
//    printf("%s\n","DEM MAP IS CELL TYPE");

    for (i=0;i<nrows;i++){
       for (j=0;j<ncols;j++){   
          pfdem1[i][j] = pfdem[in_i-i][j];   
          
       }
    }
//    out_name = param.ElevOutput->answer;
//    save_raster_C_variable(pfdem,out_name);
      close_raster_C_variable(pfdem);
}
else if(data_type_dem_map == FCELL_TYPE){
//    G_message(_("LOAD DEM DATA"));
    in_i = nrows-1;
    //FCELL** pfdem;
    
    F_pfdem = (FCELL**) G_malloc( nrows*sizeof(FCELL*));  
    FCELL** pfdem = (FCELL**) G_malloc( nrows*sizeof(FCELL*));    
    for (i = 0; i < nrows; i++){                    
         F_pfdem[i] = (FCELL*) G_malloc( ncols * sizeof(FCELL));
         pfdem[i] = (FCELL*) G_malloc( ncols * sizeof(FCELL));
    }
    
    pfdem = open_raster_F_variable(name_dem_map);
    if (param.SinkFillingOut->answer){
       F_pfdem = open_raster_F_variable(name_dem_map);
    }
//    printf("%s\n","DEM MAP IS FCELL TYPE");

    for (i=0;i<nrows;i++){
       for (j=0;j<ncols;j++){   
          pfdem1[i][j] = pfdem[in_i-i][j];
          
       }
    }
    //out_name = param.ElevOutput->answer;
    //save_raster_F_variable(pfdem,out_name);
       
    
    close_raster_F_variable(pfdem);
}
else if(data_type_dem_map == DCELL_TYPE){
//    G_message(_("LOAD DEM DATA"));
    in_i = nrows-1;
    //DCELL** pfdem;
    D_pfdem = (DCELL**) G_malloc( nrows*sizeof(DCELL*));  
    DCELL** pfdem = (DCELL**) G_malloc( nrows*sizeof(DCELL*));    
    for (i = 0; i < nrows; i++){                    
         D_pfdem[i] = (DCELL*) G_malloc( ncols * sizeof(DCELL));
         pfdem[i] = (DCELL*) G_malloc( ncols * sizeof(DCELL));
    }    
    
    pfdem = open_raster_D_variable(name_dem_map);
    if (param.SinkFillingOut->answer){
      D_pfdem = open_raster_D_variable(name_dem_map);
    }  
//    printf("%s\n","DEM MAP IS DCELL TYPE");

    for (i=0;i<nrows;i++){
        for (j=0;j<ncols;j++){   
          pfdem1[i][j] = pfdem[in_i-i][j];
        //  printf("%lf\n", pfdem1[i][j]);
        }
    }
    //out_name = param.ElevOutput->answer;
    //save_raster_D_variable(pfdem,out_name);       
     close_raster_D_variable(pfdem);
}
// =====================================================
// ============= LOAD MASK DATA ========================
if(data_type_mask_map == CELL_TYPE){
//    G_message(_("LOAD MASK DATA"));
    in_i = nrows-1;    
    CELL** pfmask;
    pfmask = open_raster_C_variable(name_mask_map);
//    printf("%s\n","MASK MAP IS CELL TYPE");

    for (i=0;i<nrows;i++){
       for (j=0;j<ncols;j++){   
          pfmask1[i][j] = pfmask[in_i-i][j];
          // printf("%lf\n", pfdem1[i][j]);
       }
    }
    //out_name = param.ElevOutput->answer;
    //save_raster_C_variable(pfdem,out_name);
     close_raster_C_variable(pfmask);
}
else if(data_type_mask_map == FCELL_TYPE){
//    G_message(_("LOAD MASK DATA"));
    in_i = nrows-1;
    FCELL** pfmask;
    pfmask = open_raster_F_variable(name_mask_map);
//    printf("%s\n","MASK MAP IS FCELL TYPE");

    for (i=0;i<nrows;i++){
       for (j=0;j<ncols;j++){   
          pfmask1[i][j] = pfmask[in_i-i][j];
       }
    }
    //out_name = param.ElevOutput->answer;
    //save_raster_F_variable(pfdem,out_name);
     close_raster_F_variable(pfmask);
}
else if(data_type_mask_map == DCELL_TYPE){
//    G_message(_("LOAD MASK DATA"));
    in_i = nrows-1;
    DCELL** pfmask;
    pfmask = open_raster_D_variable(name_mask_map);
//    printf("%s\n","MASK MAP IS DCELL TYPE");

    for (i=0;i<nrows;i++){
        for (j=0;j<ncols;j++){   
          pfmask1[i][j] = pfmask[in_i-i][j];
        }
    }
    //out_name = param.ElevOutput->answer;
    //save_raster_D_variable(pfdem,out_name);
    close_raster_D_variable(pfmask);
}

// =====================================================
// ============== LOAD BOTTOM DATA =====================
if(data_type_bottom_map == CELL_TYPE){
//    G_message(_("LOAD BOTTOM DATA"));
    in_i = nrows-1;
    CELL** pfbottom;
    pfbottom = open_raster_C_variable(name_bottom_map);
//    printf("%s\n","BOTTOM MAP IS CELL TYPE");

    for (i=0;i<nrows;i++){
       for (j=0;j<ncols;j++){   
          pfbottom1[i][j] = pfbottom[in_i-i][j];
          // printf("%lf\n", pfdem1[i][j]);
       }
    }
    //out_name = param.ElevOutput->answer;
    //save_raster_C_variable(pfdem,out_name);
    close_raster_C_variable(pfbottom);
}
else if(data_type_bottom_map == FCELL_TYPE){
//    G_message(_("LOAD BOTTOM DATA"));
    in_i = nrows-1;
    FCELL** pfbottom;
    pfbottom = open_raster_F_variable(name_bottom_map);
//    printf("%s\n","BOTTOM MAP IS FCELL TYPE");

    for (i=0;i<nrows;i++){
       for (j=0;j<ncols;j++){   
          pfbottom1[i][j] = pfbottom[in_i-i][j];
       }
    }
    //out_name = param.ElevOutput->answer;
    //save_raster_F_variable(pfdem,out_name);
    close_raster_F_variable(pfbottom);
}
else if(data_type_bottom_map == DCELL_TYPE){
//    G_message(_("LOAD BOTTOM DATA"));
    in_i = nrows-1;
    DCELL** pfbottom;
    pfbottom = open_raster_D_variable(name_bottom_map);
//    printf("%s\n","BOTTOM MAP IS DCELL TYPE");

    for (i=0;i<nrows;i++){
        for (j=0;j<ncols;j++){   
          pfbottom1[i][j] = pfbottom[in_i-i][j];
        }
    }
    close_raster_D_variable(pfbottom);
    //out_name = param.ElevOutput->answer;
    //save_raster_D_variable(pfdem,out_name);
}

// ======================================================
// ======Calculation of slopes in the x and y direction
if ((param.SinkFilling_flag->answer)||(param.Slope_x_Output->answer)||(param.Slope_y_Output->answer)){  


int isl; // x direction
int jsl; // y direction

int itcount=0;
int maxiter;
int fill_flag=0;
int wsize=1;
double mav;
double counter, increment;
int li, hi, lj, hj, ii, jj;
 
 
  if((param.SinkFilling_flag->answer)){  
//  G_message(_("CALCULATION SLOPE AND SINK-FILLING"));
  sscanf(param.SinkFillingIter->answer, "%d", &maxiter);
  sscanf(param.SinkFillingIncrement->answer, "%lf", &increment);
  
    if(param.SinkLakeMask->answer){  
        RASTER_MAP_TYPE data_type_SinkLakeMask;
        char *name_sinklakemask = param.SinkLakeMask->answer;
        data_type_SinkLakeMask = MAPdatatype(name_sinklakemask);
        double **pfSinkLakeMask1= (double**)malloc(nrows*sizeof(double*));
        for (i = 0; i < nrows; i++){                      
            pfSinkLakeMask1[i] = (double*)malloc(ncols*sizeof(double));    
        }

    if(data_type_SinkLakeMask == CELL_TYPE){
        in_i = nrows-1;    
        CELL** pfSinkLakeMask;
        pfSinkLakeMask = open_raster_C_variable(name_sinklakemask);
        for (i=0;i<nrows;i++){
           for (j=0;j<ncols;j++){   
               pfSinkLakeMask1[i][j] = pfSinkLakeMask[in_i-i][j];
           }
        }
        close_raster_C_variable(pfSinkLakeMask);
    }
    else if(data_type_SinkLakeMask == FCELL_TYPE){
        in_i = nrows-1;
        FCELL** pfSinkLakeMask;
        pfSinkLakeMask = open_raster_F_variable(name_sinklakemask);
        for (i=0;i<nrows;i++){
            for (j=0;j<ncols;j++){   
               pfSinkLakeMask1[i][j] = pfSinkLakeMask[in_i-i][j];
            }
        }
        close_raster_F_variable(pfSinkLakeMask);
    }
    else if(data_type_SinkLakeMask == DCELL_TYPE){
        in_i = nrows-1;
        DCELL** pfSinkLakeMask;
        pfSinkLakeMask = open_raster_D_variable(name_sinklakemask);
        for (i=0;i<nrows;i++){
            for (j=0;j<ncols;j++){   
                pfSinkLakeMask1[i][j] = pfSinkLakeMask[in_i-i][j];
            }
        }
        close_raster_D_variable(pfSinkLakeMask);
    }   
  
    while (fill_flag == 0){   

    fill_flag=1;
    itcount = itcount + 1;  
    
       if (itcount == maxiter){        
          G_warning(_("(SINK-FILLING) Increase the maximum number of iterations\n"
                      "or increase the increment value!"));       
       }

       for (jsl = 0; jsl < nrows; jsl++){
         for (isl = 0; isl < ncols; isl++){
              
            if(pfSinkLakeMask1[jsl][isl] != 1){
// =================================================================
// Based in:
// https://github.com/parflow/parflow/blob/master/pftools/toposlopes.c            
// =================================================================
           
           
           // Slope x direction
           if (isl == 0){
             pfSlopex1[jsl][isl] = (pfdem1[jsl][isl+1]-pfdem1[jsl][isl])/xres;
           }
          else if (isl == (ncols-1)){
              pfSlopex1[jsl][isl] = (pfdem1[jsl][isl]-pfdem1[jsl][isl-1])/xres;       
           }
           else{               
               if (((pfdem1[jsl][isl]-pfdem1[jsl][isl-1])/xres >= 0.0) && ((pfdem1[jsl][isl+1]-pfdem1[jsl][isl])/xres <= 0.0)){
                  if (fabs((pfdem1[jsl][isl]-pfdem1[jsl][isl-1])/xres) > fabs((pfdem1[jsl][isl+1]-pfdem1[jsl][isl])/xres)){
//                  pfSlopex1[jsl][isl] = (pfdem1[jsl][isl+1]-pfdem1[jsl][isl])/xres;    // forward
                    pfSlopex1[jsl][isl] = (pfdem1[jsl][isl]-pfdem1[jsl][isl-1])/xres;    // backward
                  }
                  else{                  
                    pfSlopex1[jsl][isl] = (pfdem1[jsl][isl+1]-pfdem1[jsl][isl])/xres;
                  }
               }
               else if ((pfdem1[jsl][isl+1] > pfdem1[jsl][isl]) && (pfdem1[jsl][isl-1] > pfdem1[jsl][isl])){
                  pfSlopex1[jsl][isl] = 0.0;
               }  
               else if ((pfdem1[jsl][isl]-pfdem1[jsl][isl-1])/xres < 0.0 && (pfdem1[jsl][isl+1]-pfdem1[jsl][isl])/xres < 0.0){               
                 pfSlopex1[jsl][isl] = (pfdem1[jsl][isl]-pfdem1[jsl][isl-1])/xres;               
               }             
               else if ((pfdem1[jsl][isl]-pfdem1[jsl][isl-1])/xres > 0.0 && (pfdem1[jsl][isl+1]-pfdem1[jsl][isl])/xres > 0.0){
                 pfSlopex1[jsl][isl] = (pfdem1[jsl][isl+1]-pfdem1[jsl][isl])/xres;
               }
               else if (((pfdem1[jsl][isl+1]-pfdem1[jsl][isl])/xres) == 0.0 && ((pfdem1[jsl][isl]-pfdem1[jsl][isl-1])/xres) == 0.0){
                  pfSlopex1[jsl][isl] = 0.0;
               }               
           }  
           
           // Slope y direction
           if (jsl == 0){
              pfSlopey1[jsl][isl] = (pfdem1[jsl+1][isl]-pfdem1[jsl][isl])/yres;
           }
           else if (jsl == (nrows-1)){
              pfSlopey1[jsl][isl] = (pfdem1[jsl][isl]-pfdem1[jsl-1][isl])/yres;       
           }
           else{          
               if (((pfdem1[jsl][isl]-pfdem1[jsl-1][isl])/yres >= 0.0) && ((pfdem1[jsl+1][isl]-pfdem1[jsl][isl])/yres <= 0.0)){
              
                 if(fabs((pfdem1[jsl][isl]-pfdem1[jsl-1][isl])/yres) > fabs((pfdem1[jsl+1][isl]-pfdem1[jsl][isl])/yres)){
//               pfSlopey1[jsl][isl] = (pfdem1[jsl+1][isl]-pfdem1[jsl][isl])/yres;    // forward
                 pfSlopey1[jsl][isl] = (pfdem1[jsl][isl]-pfdem1[jsl-1][isl])/yres;    // backward
                 }
                 else{
                 pfSlopey1[jsl][isl] = (pfdem1[jsl+1][isl]-pfdem1[jsl][isl])/yres;                 
                 }
              }
              else if ((pfdem1[jsl+1][isl] > pfdem1[jsl][isl]) && (pfdem1[jsl-1][isl] > pfdem1[jsl][isl])){
                 pfSlopey1[jsl][isl] = 0.0;
              }              
              else if ((pfdem1[jsl][isl]-pfdem1[jsl-1][isl])/yres < 0.0 && (pfdem1[jsl+1][isl]-pfdem1[jsl][isl])/yres < 0.0){
                 pfSlopey1[jsl][isl] = (pfdem1[jsl][isl]-pfdem1[jsl-1][isl])/yres;
              }
              else if ((pfdem1[jsl][isl]-pfdem1[jsl-1][isl])/yres > 0.0 && (pfdem1[jsl+1][isl]-pfdem1[jsl][isl])/yres > 0.0){
                 pfSlopey1[jsl][isl] = (pfdem1[jsl+1][isl]-pfdem1[jsl][isl])/yres;
              }
              else if (((pfdem1[jsl+1][isl]-pfdem1[jsl][isl])/yres) == 0.0 && ((pfdem1[jsl][isl]-pfdem1[jsl-1][isl])/yres) == 0.0){
                 pfSlopey1[jsl][isl] = 0.0;
              }              
              
           }
                    
          pfSlope[jsl][isl] = sqrt((pfSlopex1[jsl][isl] * pfSlopex1[jsl][isl]) + (pfSlopey1[jsl][isl]*pfSlopey1[jsl][isl]));

// Sink‐filling routines            
            
          if ((pfSlope[jsl][isl] == 0.0) && (itcount<=maxiter)){
              
          fill_flag=0;

              // Pit‐fill algorithm (R. M. Maxwell) 
              if (G_strcasestr(filling_method,"Pit‐fill")){  
                                 
                
               pfdem1[jsl][isl] = pfdem1[jsl][isl] + increment;      
                
                if (param.SinkFillingOut->answer){
                 if(data_type_dem_map == CELL_TYPE){
                   C_pfdem[(nrows-1)-jsl][isl] = pfdem1[jsl][isl];
                  }                  
                 else if(data_type_dem_map == FCELL_TYPE){              
                  F_pfdem[(nrows-1)-jsl][isl] = pfdem1[jsl][isl];
                  }
                 else if(data_type_dem_map == DCELL_TYPE){
                  D_pfdem[(nrows-1)-jsl][isl] = pfdem1[jsl][isl];              
                 } 
                }
                
              
              }
             
              // Moving average algorithm (R. M. Maxwell)
              if (G_strcasestr(filling_method,"Moving average")){
                 mav = 0.0;
                 counter = 0.0;
                 li = isl - wsize;
                 hi = isl + wsize;
                 lj = jsl - wsize;
                 hj = jsl + wsize;
       
                 if (jsl <= wsize){
                    lj = 0;
                 }
                 else if (hj > ((nrows-1)-wsize)){
                    hj = nrows-1;
                 }
                 if (li <= wsize){
                    li = 0;
                 }
                 else if (hi > ((ncols-1)-wsize)){
                   hi = ncols-1;          
                 }    
           
                 for(jj = lj; jj<=hj; jj++){
                    for(ii = li; ii<=hi; ii++){              
                       if ((ii != isl)||(jj != jsl)){
                         mav = mav + pfdem1[jj][ii];
                         counter = counter + 1.0;
                       }                              
                    }           
                 } 
             
                 mav = mav + increment;
                 mav = mav/(counter);
                 pfdem1[jsl][isl] = mav;      
               
                                
               if (param.SinkFillingOut->answer){ 
                if(data_type_dem_map == CELL_TYPE){
                   C_pfdem[(nrows-1)-jsl][isl] = mav;
                  }                  
                 else if(data_type_dem_map == FCELL_TYPE){              
                  F_pfdem[(nrows-1)-jsl][isl] = mav;
                  }
                 else if(data_type_dem_map == DCELL_TYPE){
                  D_pfdem[(nrows-1)-jsl][isl] = mav;              
                }                 
               }
              
              }               

          }
// ====================================================================
          }
                   
        }
     }
  
   }   
   }
   else{
   
   while (fill_flag == 0){   

    fill_flag=1;
    itcount = itcount + 1;  
    
       if (itcount == maxiter){        
          G_warning(_("(SINK-FILLING) Increase the maximum number of iterations\n"
                      "or increase the increment value!"));       
       }

       for (jsl = 0; jsl < nrows; jsl++){
         for (isl = 0; isl < ncols; isl++){
              
// =================================================================
// Based in:
// https://github.com/parflow/parflow/blob/master/pftools/toposlopes.c            
// =================================================================
           
           
           // Slope x direction
           if (isl == 0){
             pfSlopex1[jsl][isl] = (pfdem1[jsl][isl+1]-pfdem1[jsl][isl])/xres;
           }
          else if (isl == (ncols-1)){
              pfSlopex1[jsl][isl] = (pfdem1[jsl][isl]-pfdem1[jsl][isl-1])/xres;       
           }
           else{               
               if (((pfdem1[jsl][isl]-pfdem1[jsl][isl-1])/xres >= 0.0) && ((pfdem1[jsl][isl+1]-pfdem1[jsl][isl])/xres <= 0.0)){
                  if (fabs((pfdem1[jsl][isl]-pfdem1[jsl][isl-1])/xres) > fabs((pfdem1[jsl][isl+1]-pfdem1[jsl][isl])/xres)){
//                  pfSlopex1[jsl][isl] = (pfdem1[jsl][isl+1]-pfdem1[jsl][isl])/xres;    // forward
                    pfSlopex1[jsl][isl] = (pfdem1[jsl][isl]-pfdem1[jsl][isl-1])/xres;    // backward
                  }
                  else{                  
                    pfSlopex1[jsl][isl] = (pfdem1[jsl][isl+1]-pfdem1[jsl][isl])/xres;
                  }
               }
               else if ((pfdem1[jsl][isl+1] > pfdem1[jsl][isl]) && (pfdem1[jsl][isl-1] > pfdem1[jsl][isl])){
                  pfSlopex1[jsl][isl] = 0.0;
               }  
               else if ((pfdem1[jsl][isl]-pfdem1[jsl][isl-1])/xres < 0.0 && (pfdem1[jsl][isl+1]-pfdem1[jsl][isl])/xres < 0.0){               
                 pfSlopex1[jsl][isl] = (pfdem1[jsl][isl]-pfdem1[jsl][isl-1])/xres;               
               }             
               else if ((pfdem1[jsl][isl]-pfdem1[jsl][isl-1])/xres > 0.0 && (pfdem1[jsl][isl+1]-pfdem1[jsl][isl])/xres > 0.0){
                 pfSlopex1[jsl][isl] = (pfdem1[jsl][isl+1]-pfdem1[jsl][isl])/xres;
               }
               else if (((pfdem1[jsl][isl+1]-pfdem1[jsl][isl])/xres) == 0.0 && ((pfdem1[jsl][isl]-pfdem1[jsl][isl-1])/xres) == 0.0){
                  pfSlopex1[jsl][isl] = 0.0;
               }               
           }  
           
           // Slope y direction
           if (jsl == 0){
              pfSlopey1[jsl][isl] = (pfdem1[jsl+1][isl]-pfdem1[jsl][isl])/yres;
           }
           else if (jsl == (nrows-1)){
              pfSlopey1[jsl][isl] = (pfdem1[jsl][isl]-pfdem1[jsl-1][isl])/yres;       
           }
           else{          
               if (((pfdem1[jsl][isl]-pfdem1[jsl-1][isl])/yres >= 0.0) && ((pfdem1[jsl+1][isl]-pfdem1[jsl][isl])/yres <= 0.0)){
              
                 if(fabs((pfdem1[jsl][isl]-pfdem1[jsl-1][isl])/yres) > fabs((pfdem1[jsl+1][isl]-pfdem1[jsl][isl])/yres)){
//               pfSlopey1[jsl][isl] = (pfdem1[jsl+1][isl]-pfdem1[jsl][isl])/yres;    // forward
                 pfSlopey1[jsl][isl] = (pfdem1[jsl][isl]-pfdem1[jsl-1][isl])/yres;    // backward
                 }
                 else{
                 pfSlopey1[jsl][isl] = (pfdem1[jsl+1][isl]-pfdem1[jsl][isl])/yres;                 
                 }
              }
              else if ((pfdem1[jsl+1][isl] > pfdem1[jsl][isl]) && (pfdem1[jsl-1][isl] > pfdem1[jsl][isl])){
                 pfSlopey1[jsl][isl] = 0.0;
              }              
              else if ((pfdem1[jsl][isl]-pfdem1[jsl-1][isl])/yres < 0.0 && (pfdem1[jsl+1][isl]-pfdem1[jsl][isl])/yres < 0.0){
                 pfSlopey1[jsl][isl] = (pfdem1[jsl][isl]-pfdem1[jsl-1][isl])/yres;
              }
              else if ((pfdem1[jsl][isl]-pfdem1[jsl-1][isl])/yres > 0.0 && (pfdem1[jsl+1][isl]-pfdem1[jsl][isl])/yres > 0.0){
                 pfSlopey1[jsl][isl] = (pfdem1[jsl+1][isl]-pfdem1[jsl][isl])/yres;
              }
              else if (((pfdem1[jsl+1][isl]-pfdem1[jsl][isl])/yres) == 0.0 && ((pfdem1[jsl][isl]-pfdem1[jsl-1][isl])/yres) == 0.0){
                 pfSlopey1[jsl][isl] = 0.0;
              }              
              
           }
                    
          pfSlope[jsl][isl] = sqrt((pfSlopex1[jsl][isl] * pfSlopex1[jsl][isl]) + (pfSlopey1[jsl][isl]*pfSlopey1[jsl][isl]));

// Sink‐filling routines            
            
          if ((pfSlope[jsl][isl] == 0.0) && (itcount<=maxiter)){
              
          fill_flag=0;

              // Pit‐fill algorithm (R. M. Maxwell) 
              if (G_strcasestr(filling_method,"Pit‐fill")){  
                                 
                
               pfdem1[jsl][isl] = pfdem1[jsl][isl] + increment;      
                
                if (param.SinkFillingOut->answer){
                 if(data_type_dem_map == CELL_TYPE){
                   C_pfdem[(nrows-1)-jsl][isl] = pfdem1[jsl][isl];
                  }                  
                 else if(data_type_dem_map == FCELL_TYPE){              
                  F_pfdem[(nrows-1)-jsl][isl] = pfdem1[jsl][isl];
                  }
                 else if(data_type_dem_map == DCELL_TYPE){
                  D_pfdem[(nrows-1)-jsl][isl] = pfdem1[jsl][isl];              
                 } 
                }
                
              
              }
             
              // Moving average algorithm (R. M. Maxwell)
              if (G_strcasestr(filling_method,"Moving average")){
                 mav = 0.0;
                 counter = 0.0;
                 li = isl - wsize;
                 hi = isl + wsize;
                 lj = jsl - wsize;
                 hj = jsl + wsize;
       
                 if (jsl <= wsize){
                    lj = 0;
                 }
                 else if (hj > ((nrows-1)-wsize)){
                    hj = nrows-1;
                 }
                 if (li <= wsize){
                    li = 0;
                 }
                 else if (hi > ((ncols-1)-wsize)){
                   hi = ncols-1;          
                 }    
           
                 for(jj = lj; jj<=hj; jj++){
                    for(ii = li; ii<=hi; ii++){              
                       if ((ii != isl)||(jj != jsl)){
                         mav = mav + pfdem1[jj][ii];
                         counter = counter + 1.0;
                       }                              
                    }           
                 } 
             
                 mav = mav + increment;
                 mav = mav/(counter);
                 pfdem1[jsl][isl] = mav;      
               
                                
               if (param.SinkFillingOut->answer){ 
                if(data_type_dem_map == CELL_TYPE){
                   C_pfdem[(nrows-1)-jsl][isl] = mav;
                  }                  
                 else if(data_type_dem_map == FCELL_TYPE){              
                  F_pfdem[(nrows-1)-jsl][isl] = mav;
                  }
                 else if(data_type_dem_map == DCELL_TYPE){
                  D_pfdem[(nrows-1)-jsl][isl] = mav;              
                }                 
               }
              
              }               

          }
// ====================================================================
                             
        }
     }
  
   } 
   
   
   }   
   
     if (param.SinkFillingOut->answer){
      if(data_type_dem_map == CELL_TYPE){
           char* out_filled_dem_name = param.SinkFillingOut->answer;
           save_raster_C_variable(C_pfdem,out_filled_dem_name);
           close_raster_C_variable(C_pfdem);
      }                  
      else if(data_type_dem_map == FCELL_TYPE){              
           char* out_filled_dem_name = param.SinkFillingOut->answer;
           save_raster_F_variable(F_pfdem,out_filled_dem_name);
           close_raster_F_variable(F_pfdem);
      }
      else if(data_type_dem_map == DCELL_TYPE){
           char* out_filled_dem_name = param.SinkFillingOut->answer;
           save_raster_D_variable(D_pfdem,out_filled_dem_name);        
           close_raster_D_variable(D_pfdem);      
      } 
     }
    
  }      

  if((param.Slope_x_Output->answer)||(param.Slope_y_Output->answer)){  
//   G_message(_("CALCULATION X AND Y SLOPE"));
  

   if (G_strcasestr(slope_method,"Global Slope Enforcement (GSE)")){
//   if (param.GlobalSlopeEnforcement_flag->answer){
   G_message(_("GLOBAL SLOPE ENFORCEMENT METHOD")); 
   
   double riverSlope, minSlope;
   
   sscanf(param.SlopeRiverSlope->answer, "%lf", &riverSlope);
   sscanf(param.SlopeMinSlope->answer, "%lf", &minSlope);
   
   char *name_flowdirection_map = param.SlopeFlowdir->answer;  // Flow direction map name
   char *name_stream_map = param.SlopeRiver->answer;           //  Stream segments map name

   CELL** flowDir;
//   int flowDir1[nrows][ncols];
   int **flowDir1 = (int**)malloc(nrows*sizeof(int*));
   in_i = nrows-1;
   flowDir = open_raster_C_variable(name_flowdirection_map);    
   for (i=0;i<nrows;i++){
      flowDir1[i] = (int*)malloc(ncols*sizeof(int));
      for (j=0;j<ncols;j++){   
       flowDir1[i][j] = flowDir[in_i-i][j];
      }
   }
   close_raster_C_variable(flowDir);
   CELL** streamSeg;
//   int streamSeg1[nrows][ncols];
   int **streamSeg1 = (int**)malloc(nrows*sizeof(int*));
   in_i = nrows-1;
   streamSeg = open_raster_C_variable(name_stream_map);    
   for (i=0;i<nrows;i++){
      streamSeg1[i] = (int*)malloc(ncols*sizeof(int));
      for (j=0;j<ncols;j++){   
       streamSeg1[i][j] = streamSeg[in_i-i][j];
      }
   }
   close_raster_C_variable(streamSeg);
   
   // Based in: 
   // Global Topographic Slope Enforcement to Ensure Connectivity and Drainage in an Urban Terrain
   // Michael L. Barnes; Claire Welty and Andrew J. Miller. 
   // Journal of Hydrologic Engineering
   // Volume 21 Issue 4 - April 2016    
   
   for (jsl=0;jsl<nrows;jsl++){
       for (isl=0;isl<ncols;isl++){
              
//     slopes magnitudes in x and y direction     
             
                   if (jsl < (nrows-1) && isl < (ncols-1)){                   
                        pfSlopey1[jsl][isl] =(pfdem1[jsl+1][isl]-pfdem1[jsl][isl])/(yres); 
                        pfSlopex1[jsl][isl] =(pfdem1[jsl][isl+1]-pfdem1[jsl][isl])/(xres);                        
                   }         
                             
                   if (jsl == (nrows-1) && isl < (ncols-1)){                   
                        pfSlopey1[jsl][isl] =pfSlopey1[jsl-1][isl];
                        pfSlopex1[jsl][isl] =(pfdem1[jsl][isl+1]-pfdem1[jsl][isl])/(xres);                         
                   }           
                             
                   if (jsl<(nrows-1) && isl == (ncols-1)){                   
                        pfSlopey1[jsl][isl] =(pfdem1[jsl+1][isl]-pfdem1[jsl][isl])/(yres);
                        pfSlopex1[jsl][isl] =pfSlopex1[jsl][isl-1];                        
                   }
                   
                   if (jsl == (nrows-1) && isl == (ncols-1)){
                        pfSlopey1[jsl][isl] =pfSlopey1[jsl-1][isl];
                        pfSlopex1[jsl][isl] =pfSlopex1[jsl][isl-1];
                   }
                   
//     impose minimum slope threshold

                   if (fabs(pfSlopex1[jsl][isl]) < minSlope){
                       pfSlopex1[jsl][isl] = minSlope;
                   }
                   if(fabs(pfSlopey1[jsl][isl]) < minSlope){
                       pfSlopey1[jsl][isl] = minSlope;
                   }
                   
//     give channel network a constant linear slope 

                   if (streamSeg1[jsl][isl] != 0){                   
                        pfSlopex1[jsl][isl] = riverSlope;
                        pfSlopey1[jsl][isl] = riverSlope;                        
                   }
                   
//     enforce single flow direction throughout domain
//     based on GRASS r.watershed derived flow direction
//     NORTH
                   if(fabs(flowDir1[jsl][isl]) == 2){ 
                        pfSlopey1[jsl][isl] = -fabs(pfSlopey1[jsl][isl]);
                        if(streamSeg1[jsl][isl] != 0){
                             pfSlopex1[jsl][isl] = 0.0;
                        }
                   //     continue;                        
                   }
//     WEST
                   else if(fabs(flowDir1[jsl][isl]) == 4){
                        pfSlopex1[jsl][isl] = fabs(pfSlopex1[jsl][isl]);
                        if(streamSeg1[jsl][isl] != 0){
                             pfSlopey1[jsl][isl] = 0.0;
                        }
                   //     continue;                        
                   }
//     SOUTH 
                   else if(fabs(flowDir1[jsl][isl]) == 6){
                        pfSlopey1[jsl][isl] = fabs(pfSlopey1[jsl][isl]);
                        if(streamSeg1[jsl][isl] != 0){
                             pfSlopex1[jsl][isl] = 0.0;
                        }
                  //      continue;                        
                   }
//     EAST         
                   else if(fabs(flowDir1[jsl][isl]) == 8){
                        pfSlopex1[jsl][isl] = -fabs(pfSlopex1[jsl][isl]);
                        if(streamSeg1[jsl][isl] != 0){
                             pfSlopey1[jsl][isl] = 0.0;
                        }
                 //       continue;
                   }     
                   else {                        
                        G_warning(_("GLOBAL SLOPE ENFORCEMENT METHOD: no conditions met"));                        
                   }
                   
             
          }
      }         
      
   }  
   //else if(G_strcasestr(slope_method,"Upwind scheme")){
   else if(G_strcasestr(slope_method,"Upwind scheme")){
   // Based in:
   // https://github.com/parflow/parflow/blob/master/pftools/toposlopes.c
   
   G_message(_("SLOPE: UPWIND SCHEME"));
   for (jsl = 0; jsl < nrows; jsl++){
         for (isl = 0; isl < ncols; isl++){                     
           // Slope x direction
           if (isl == 0){
             pfSlopex1[jsl][isl] = (pfdem1[jsl][isl+1]-pfdem1[jsl][isl])/xres;
           }
           else if (isl == (ncols-1)){
              pfSlopex1[jsl][isl] = (pfdem1[jsl][isl]-pfdem1[jsl][isl-1])/xres;       
           }
           else{
               if (((pfdem1[jsl][isl]-pfdem1[jsl][isl-1])/xres >= 0.0) && ((pfdem1[jsl][isl+1]-pfdem1[jsl][isl])/xres <= 0.0)){
                  if (fabs((pfdem1[jsl][isl]-pfdem1[jsl][isl-1])/xres) > fabs((pfdem1[jsl][isl+1]-pfdem1[jsl][isl])/xres)){
//                  pfSlopex1[jsl][isl] = (pfdem1[jsl][isl+1]-pfdem1[jsl][isl])/xres;    // forward
                    pfSlopex1[jsl][isl] = (pfdem1[jsl][isl]-pfdem1[jsl][isl-1])/xres;    // backward
                  }
                  else{                  
                    pfSlopex1[jsl][isl] = (pfdem1[jsl][isl+1]-pfdem1[jsl][isl])/xres;
                  }
               }
               else if ((pfdem1[jsl][isl+1] > pfdem1[jsl][isl]) && (pfdem1[jsl][isl-1] > pfdem1[jsl][isl])){
                  pfSlopex1[jsl][isl] = 0.0;
               }                 
               else if ((pfdem1[jsl][isl]-pfdem1[jsl][isl-1])/xres < 0.0 && (pfdem1[jsl][isl+1]-pfdem1[jsl][isl])/xres < 0.0){               
                 pfSlopex1[jsl][isl] = (pfdem1[jsl][isl]-pfdem1[jsl][isl-1])/xres;               
               }             
               else if ((pfdem1[jsl][isl]-pfdem1[jsl][isl-1])/xres > 0.0 && (pfdem1[jsl][isl+1]-pfdem1[jsl][isl])/xres > 0.0){
                 pfSlopex1[jsl][isl] = (pfdem1[jsl][isl+1]-pfdem1[jsl][isl])/xres;
               }      
               else if (((pfdem1[jsl][isl+1]-pfdem1[jsl][isl])/xres) == 0.0 && ((pfdem1[jsl][isl]-pfdem1[jsl][isl-1])/xres) == 0.0){
                  pfSlopex1[jsl][isl] = 0.0;
               }         
           }  
           
           // Slope y direction
           if (jsl == 0){
              pfSlopey1[jsl][isl] = (pfdem1[jsl+1][isl]-pfdem1[jsl][isl])/yres;
           }
           else if (jsl == (nrows-1)){
              pfSlopey1[jsl][isl] = (pfdem1[jsl][isl]-pfdem1[jsl-1][isl])/yres;       
           }
           else{            
               if (((pfdem1[jsl][isl]-pfdem1[jsl-1][isl])/yres >= 0.0) && ((pfdem1[jsl+1][isl]-pfdem1[jsl][isl])/yres <= 0.0)){
              
                 if(fabs((pfdem1[jsl][isl]-pfdem1[jsl-1][isl])/yres) > fabs((pfdem1[jsl+1][isl]-pfdem1[jsl][isl])/yres)){
//               pfSlopey1[jsl][isl] = (pfdem1[jsl+1][isl]-pfdem1[jsl][isl])/yres;    // forward
                 pfSlopey1[jsl][isl] = (pfdem1[jsl][isl]-pfdem1[jsl-1][isl])/yres;    // backward
                 }
                 else{
                 pfSlopey1[jsl][isl] = (pfdem1[jsl+1][isl]-pfdem1[jsl][isl])/yres;                 
                 }
              }
              else if ((pfdem1[jsl+1][isl] > pfdem1[jsl][isl]) && (pfdem1[jsl-1][isl] > pfdem1[jsl][isl])){
                 pfSlopey1[jsl][isl] = 0.0;
              }              
              else if ((pfdem1[jsl][isl]-pfdem1[jsl-1][isl])/yres < 0.0 && (pfdem1[jsl+1][isl]-pfdem1[jsl][isl])/yres < 0.0){
                 pfSlopey1[jsl][isl] = (pfdem1[jsl][isl]-pfdem1[jsl-1][isl])/yres;
              }
              else if ((pfdem1[jsl][isl]-pfdem1[jsl-1][isl])/yres > 0.0 && (pfdem1[jsl+1][isl]-pfdem1[jsl][isl])/yres > 0.0){
                 pfSlopey1[jsl][isl] = (pfdem1[jsl+1][isl]-pfdem1[jsl][isl])/yres;
              }
              else if (((pfdem1[jsl+1][isl]-pfdem1[jsl][isl])/yres) == 0.0 && ((pfdem1[jsl][isl]-pfdem1[jsl-1][isl])/yres) == 0.0){
                 pfSlopey1[jsl][isl] = 0.0;
              }              
              
           }
      }
   }
   }  
   else if (G_strcasestr(slope_method,"Central difference")){
   G_message(_("SLOPE: CENTRAL DIFFERENCE METHOD"));
   
   for (jsl = 0; jsl < nrows; jsl++){
     for (isl = 0; isl < ncols; isl++){       
       // Slope x direction
       if (isl == 0){
          pfSlopex1[jsl][isl] = (pfdem1[jsl][isl+1]-pfdem1[jsl][isl])/xres;
       }
       else if (isl == (ncols-1)){
          pfSlopex1[jsl][isl] = (pfdem1[jsl][isl]-pfdem1[jsl][isl-1])/xres;       
       }
       else{       
           if ((pfdem1[jsl][isl+1] > pfdem1[jsl][isl]) && (pfdem1[jsl][isl-1] > pfdem1[jsl][isl])){
              pfSlopex1[jsl][isl] = 0.0;
           }
           else{
               pfSlopex1[jsl][isl] = (pfdem1[jsl][isl+1]-pfdem1[jsl][isl-1])/(2*xres); // central
           }
       }       
       // Slope y direction
       if (jsl == 0){
          pfSlopey1[jsl][isl] = (pfdem1[jsl+1][isl]-pfdem1[jsl][isl])/yres;
       }
       else if (jsl == (nrows-1)){
          pfSlopey1[jsl][isl] = (pfdem1[jsl][isl]-pfdem1[jsl-1][isl])/yres;       
       }
       else{ 
           if ((pfdem1[jsl+1][isl] > pfdem1[jsl][isl]) && (pfdem1[jsl-1][isl] > pfdem1[jsl][isl])){
              pfSlopey1[jsl][isl] = 0.0;
           }
           else{
               pfSlopey1[jsl][isl] = (pfdem1[jsl+1][isl]-pfdem1[jsl-1][isl])/(2*yres); // central
           }
       }
     }
   }
   
   }
  }  
}

//G_message(_("SAVING DIGITAL ELEVATION MODEL"));

char *PathOutElev;// Digital Elevation Model   
PathOutElev = concat(out_dir, "/elevations.pfb");     
pfb_maker(ncols,pfdem1, nrows,1,xcorner,ycorner,0.0, xres, yres,0.0, PathOutElev);


// Domain Max and Min values

double max_dem_value=-10000000000, min_bottom_value= 10000000000, max_subsurfacedepth=-10000000000;//pfdem1[0][0]  pfbottom1[0][0];
 //   int i,j;

    for (i=0; i<nrows; i++){
	for (j=0; j<ncols; j++){	

            if((pfdem1[i][j]> max_dem_value)){
	      max_dem_value= pfdem1[i][j];
	    }

	    if((pfbottom1[i][j]<min_bottom_value)){
              min_bottom_value= pfbottom1[i][j];
            }			
	
	}
	
    }         
      


int NZ;
double zres;
double zcorner;

sscanf(param.NZ_layers->answer, "%d", &NZ);
//sscanf(param.dz_resolution->answer, "%lf", &zres);
zres = 0.00;
  
// ======================================================
// ============== Make Terrain Following Grid ===========

if(G_strcasestr(grid_format,"Terrain Following Grid")){

G_message(_("GRID FORMAT: TERRAIN FOLLOWING GRID"));
                       
    for (i=0;i<nrows;i++){
        for (j=0;j<ncols;j++){  
            
            //if((zres==0.00) && (pfmask1[i][j]==1)){
            //    zres = (pfdem1[i][j] - pfbottom1[i][j])/NZ;            
            //}
            
            if((pfmask1[i][j]==1)&&((pfdem1[i][j] - pfbottom1[i][j])>max_subsurfacedepth)){
                max_subsurfacedepth = (pfdem1[i][j] - pfbottom1[i][j]);            
            }
            
            zres = max_subsurfacedepth/NZ;
                             
          pfbottom1[i][j] = (max_dem_value) + (pfbottom1[i][j]-pfdem1[i][j]); 
          pfdem1[i][j] = max_dem_value;
          
        }
    }   
    // 
    zcorner = max_dem_value - max_subsurfacedepth;    
}
else{
   G_message(_("GRID FORMAT: ORTHOGONAL GRID"));

   zcorner=min_bottom_value;
   
   zres = (max_dem_value-min_bottom_value)/NZ;
	
}

// ============================================================================
// ============== SURFACE BOUNDARY CONDITION DATA TO .pfb FORMAT =====================

if(param.SurfaceCond_flag->answer){

	RASTER_MAP_TYPE data_type_BC_P_surf;
	char *name_BC_P_surf = param.SurfaceCond->answer;        // Surface Condition map name 
	char *out_BC_P_name = param.SurfaceCond_OutName->answer;
	data_type_BC_P_surf = MAPdatatype(name_BC_P_surf);
	
	double **pfBC_P_surf1 = (double**)malloc(nrows*sizeof(double*));
	for (i = 0; i < nrows; i++){ 
	    pfBC_P_surf1[i] = (double*)malloc(ncols*sizeof(double));
	}


	if(data_type_BC_P_surf == CELL_TYPE){
//	    G_message(_("LOAD SURFACE CONDITION DATA"));
	    in_i = nrows-1;
	    CELL** pfBC_P_surf;
	    pfBC_P_surf = open_raster_C_variable(name_BC_P_surf);


	    for (i=0;i<nrows;i++){
	       for (j=0;j<ncols;j++){   
	          pfBC_P_surf1[i][j] = pfBC_P_surf[in_i-i][j];

	       }
	    }

	    close_raster_C_variable(pfBC_P_surf);
	}
	else if(data_type_BC_P_surf == FCELL_TYPE){
//	    G_message(_("LOAD SURFACE CONDITION DATA"));
	    in_i = nrows-1;
	    FCELL** pfBC_P_surf;
	    pfBC_P_surf = open_raster_F_variable(name_BC_P_surf);


	    for (i=0;i<nrows;i++){
	       for (j=0;j<ncols;j++){   
	          pfBC_P_surf1[i][j] = pfBC_P_surf[in_i-i][j];
	       }
	    }

	    close_raster_F_variable(pfBC_P_surf);
	}
	else if(data_type_BC_P_surf == DCELL_TYPE){
//	    G_message(_("LOAD SURFACE CONDITION DATA"));
	    in_i = nrows-1;
	    DCELL** pfBC_P_surf;
	    pfBC_P_surf = open_raster_D_variable(name_BC_P_surf);


	    for (i=0;i<nrows;i++){
	        for (j=0;j<ncols;j++){   
	          pfBC_P_surf1[i][j] = pfBC_P_surf[in_i-i][j];
	        }
	    }
	    close_raster_D_variable(pfBC_P_surf);

	}
	
	// Write Surface condition
//        G_message(_("SAVING SURFACE CONDITION MAP"));
        char *PathOutSPBC;        
        PathOutSPBC = concat(out_dir, "/");
        PathOutSPBC = concat(PathOutSPBC, out_BC_P_name);
        PathOutSPBC = concat(PathOutSPBC, ".pfb");
        
        if(G_strcasestr(grid_format,"Terrain Following Grid")){
           TFG_Surf_condition_write(pfBC_P_surf1,ncols, nrows,NZ,xcorner,ycorner,zcorner, xres, yres,zres, PathOutSPBC);
        } else if(G_strcasestr(grid_format,"Orthogonal Grid")){
           OG_Surf_condition_write(pfdem1,pfBC_P_surf1,ncols, nrows,NZ,xcorner,ycorner,zcorner, xres, yres,zres, PathOutSPBC);        
        }
                           
        if (param.Out_vtk_flag-> answer){           
//        	G_message(_("CREATING SURFACE CONDITION (VTK FILE)"));  
        	char *SurfBCid = "SurfaceCondition";
               char *PathOutSurfBCVTK;
               PathOutSurfBCVTK = concat(out_dir, "/"); 
               PathOutSurfBCVTK = concat(PathOutSurfBCVTK, out_BC_P_name);
               PathOutSurfBCVTK = concat(PathOutSurfBCVTK, ".vtk");
               FILE *file_SurfBCVTK;
               file_SurfBCVTK = fopen(PathOutSurfBCVTK, "wb");            
	       PrintVTK_2D(file_SurfBCVTK,pfBC_P_surf1,ncols,nrows,xcorner,ycorner,0.0,xres,yres,SurfBCid);  
	       fclose(file_SurfBCVTK);  
        }
        	
}
        
// ********************************************************
//              Start output files .pfb
// ********************************************************	
G_message(_("..."));
			
char *out_pfsol_file = param.Solid_pfsolOutput->answer;
		
char *PathOutMask;// Mask (Active domain and patches)	
PathOutMask = concat(out_dir, "/mask.pfb");	
pfb_maker(ncols,pfmask1, nrows,1,xcorner,ycorner,0.0,xres,yres,0.0,PathOutMask); // 
                       
char *PathOutTop;// Topography (Digital Elevation Model)     
PathOutTop = concat(out_dir, "/topography.pfb");     
pfb_maker(ncols,pfdem1, nrows,1,xcorner,ycorner,0.0, xres,yres,0.0,PathOutTop); //           
                    
char *PathOutBot;// Bottom domain;   
PathOutBot = concat(out_dir, "/bottom.pfb");
pfb_maker(ncols,pfbottom1, nrows,1,xcorner,ycorner,0.0,xres,yres,0.0,PathOutBot); //  
       
if (param.Slope_x_Output->answer){
//	G_message(_("CREATING SLOPE X-DIRECTION (PFB FILE)"));
        SlopeX_file = concat(out_SlopeX_name,".pfb");
	char *PathOutSlopeX;// Digital Elevation Model   
	PathOutSlopeX = concat(out_dir, "/");
	PathOutSlopeX = concat(PathOutSlopeX, SlopeX_file);
	pfb_maker(ncols,pfSlopex1, nrows,1,xcorner,ycorner,0.0,xres, yres,0.0, PathOutSlopeX);
	G_free(PathOutSlopeX);
                                                                     
	if (param.Out_vtk_flag-> answer){
//               G_message(_("CREATING SLOPE X-DIRECTION (VTK FILE)"));
               char *slopeidx = "x_slope";
               char *PathOutSlopeXVTK;
               PathOutSlopeXVTK = concat(out_dir, "/"); 
               PathOutSlopeXVTK = concat(PathOutSlopeXVTK, out_SlopeX_name);
               PathOutSlopeXVTK = concat(PathOutSlopeXVTK, ".vtk");
               FILE *file_slopexVTK;
               file_slopexVTK = fopen(PathOutSlopeXVTK, "wb");            
	       PrintVTK_2D(file_slopexVTK,pfSlopex1,ncols,nrows,xcorner,ycorner,0.0,xres,yres,slopeidx);  
	       fclose(file_slopexVTK);  
	          
	}
}

if (param.Slope_y_Output->answer){
//	G_message(_("CREATING SLOPE Y-DIRECTION (PFB FILE)"));
        SlopeY_file = concat(out_SlopeY_name,".pfb");
	char *PathOutSlopeY;// Digital Elevation Model   
	PathOutSlopeY = concat(out_dir, "/");
	PathOutSlopeY = concat(PathOutSlopeY, SlopeY_file);
	pfb_maker(ncols,pfSlopey1, nrows,1,xcorner,ycorner,0.0, xres, yres,0.0, PathOutSlopeY);
	G_free(PathOutSlopeY);
	          
	if (param.Out_vtk_flag-> answer){
//		G_message(_("CREATING SLOPE Y-DIRECTION (VTK FILE)"));
		char *slopeidy = "y_slope";     
		char *PathOutSlopeYVTK;
               PathOutSlopeYVTK = concat(out_dir, "/"); 
               PathOutSlopeYVTK = concat(PathOutSlopeYVTK, out_SlopeY_name);
               PathOutSlopeYVTK = concat(PathOutSlopeYVTK, ".vtk");
               FILE *file_slopeyVTK;
               file_slopeyVTK = fopen(PathOutSlopeYVTK, "wb");
               PrintVTK_2D(file_slopeyVTK,pfSlopey1,ncols,nrows,xcorner,ycorner,0.0,xres,yres,slopeidy);  
	       fclose(file_slopeyVTK);       
	}           

} 
         
// ********************************************************
//          Working directory changed to: out_dir   
// ********************************************************        
        
int st_dir = chdir(out_dir); 
       
// ********************************************************
//         CALL Parflow function (pfpatchysolid)
// ********************************************************
	
	
//G_message(_("CREATING SOLID FILE"));
        	
char *PathSolidFile;
PathSolidFile = concat(out_dir,"/");
PathSolidFile = concat(PathSolidFile,solid_name);
//PathSolidFile = concat(PathSolidFile,".pfsol");
	
FILE *file_solidfile;
file_solidfile = fopen(PathSolidFile, "w");	
	
char *PathSolidVTKFile;	
PathSolidVTKFile = concat(out_dir,"/");
PathSolidVTKFile = concat(PathSolidVTKFile,out_pfsol_name);
PathSolidVTKFile = concat(PathSolidVTKFile,".vtk");
	
FILE *file_solidfileVTK;
file_solidfileVTK = fopen(PathSolidVTKFile, "w");		

char * PatchesNames;
char * PatchesNames2;	

PatchesNames = MakePatchySolid(ncols,pfmask1,pfdem1,pfbottom1,nrows,1,xcorner,ycorner, xres, yres,0,0,file_solidfile,file_solidfileVTK);

PatchesNames2 = PatchesNames;

// *******************************************************
//               Create initial .tcl file
// *******************************************************    
if(param.Out_InitTCLfile_flag->answer){

   char *PathRunTclFile;
   PathRunTclFile = concat(out_dir, "/InitialRunFile.tcl");
   FILE *run_tcl_file;
   run_tcl_file = fopen(PathRunTclFile, "w");   
   fprintf(run_tcl_file, "%s\n","#Initial setup file for running the Parflow model.");
   fprintf(run_tcl_file, "%s\n","#");  
   fprintf(run_tcl_file, "%s\n","set tcl_precision 17");
   fprintf(run_tcl_file, "%s\n","#");
   fprintf(run_tcl_file, "%s\n","# Import the ParFlow TCL package");
   fprintf(run_tcl_file, "%s\n","#");
   fprintf(run_tcl_file, "%s\n","lappend auto_path $env(PARFLOW_DIR)/bin"); 
   fprintf(run_tcl_file, "%s\n","package require parflow");
   fprintf(run_tcl_file, "%s\n","namespace import Parflow::*");
   fprintf(run_tcl_file, "%s\n","pfset FileVersion 4");

   fprintf(run_tcl_file, "%s\n","set runname InitialRunFile");

   fprintf(run_tcl_file, "%s\n","pfset Process.Topology.P        [lindex $argv 0]");
   fprintf(run_tcl_file, "%s\n","pfset Process.Topology.Q        [lindex $argv 1]");
   fprintf(run_tcl_file, "%s\n","pfset Process.Topology.R        [lindex $argv 2]");

   fprintf(run_tcl_file, "%s%lf\n","set dx ",xres);
   fprintf(run_tcl_file, "%s\n","set dy $dx");
   fprintf(run_tcl_file, "%s%lf\n","set dz ", zres);

   fprintf(run_tcl_file, "%s%d\n","set nx ", ncols);
   fprintf(run_tcl_file, "%s%d\n","set ny ", nrows);
   fprintf(run_tcl_file, "%s%d\n","set nz ", NZ);
        
   fprintf(run_tcl_file, "%s%lf\n","set x0 ", xcorner);
   fprintf(run_tcl_file, "%s%lf\n","set y0 ", ycorner);
   fprintf(run_tcl_file, "%s%lf\n","set z0 ", zcorner);

   fprintf(run_tcl_file, "%s\n","set xmax [expr $x0 + ($nx * $dx)]");
   fprintf(run_tcl_file, "%s\n","set ymax [expr $y0 + ($ny * $dy)]");
   fprintf(run_tcl_file, "%s\n","set zmax [expr $z0 + ($nz * $dz)]");
   
   fprintf(run_tcl_file, "%s\n","#---------------------------------------------------------");
   fprintf(run_tcl_file, "%s\n","# Computational Grid");
   fprintf(run_tcl_file, "%s\n","#---------------------------------------------------------");
   fprintf(run_tcl_file, "%s\n","pfset ComputationalGrid.Lower.X                 $x0");
   fprintf(run_tcl_file, "%s\n","pfset ComputationalGrid.Lower.Y                 $y0");
   fprintf(run_tcl_file, "%s\n","pfset ComputationalGrid.Lower.Z                 $z0");
   fprintf(run_tcl_file, "%s\n","pfset ComputationalGrid.DX	                 $dx");
   fprintf(run_tcl_file, "%s\n","pfset ComputationalGrid.DY                      $dy");
   fprintf(run_tcl_file, "%s\n","pfset ComputationalGrid.DZ	                 $dz");
   fprintf(run_tcl_file, "%s\n","pfset ComputationalGrid.NX                      $nx");
   fprintf(run_tcl_file, "%s\n","pfset ComputationalGrid.NY                      $ny");
   fprintf(run_tcl_file, "%s\n","pfset ComputationalGrid.NZ                      $nz");
   fprintf(run_tcl_file, "%s\n","#---------------------------------------------------------");
   fprintf(run_tcl_file, "%s\n","# The Names of the GeomInputs");
   fprintf(run_tcl_file, "%s\n","#---------------------------------------------------------");
   
   fprintf(run_tcl_file, "%s\n","pfset GeomInput.Names                 \"domainbox solidfile\"");
   fprintf(run_tcl_file, "%s\n","pfset GeomInput.domainbox.InputType    Box");
   fprintf(run_tcl_file, "%s\n","pfset GeomInput.domainbox.GeomName     background");
   fprintf(run_tcl_file, "%s\n","pfset GeomInput.solidfile.InputType    SolidFile");
   fprintf(run_tcl_file, "%s\n","pfset GeomInput.solidfile.GeomNames    domain");
   fprintf(run_tcl_file, "%s%s\n","pfset GeomInput.solidfile.FileName     ",solid_name);
   fprintf(run_tcl_file, "%s%s%s\n","pfset Geom.domain.Patches              \"",PatchesNames,"\"");
   
   fprintf(run_tcl_file, "%s\n","#---------------------------------------------------------");
   fprintf(run_tcl_file, "%s\n","# Domain ");
   fprintf(run_tcl_file, "%s\n","#---------------------------------------------------------");
   fprintf(run_tcl_file, "%s\n","pfset Domain.GeomName domain");
   fprintf(run_tcl_file, "%s\n","#---------------------------------------------------------");
   fprintf(run_tcl_file, "%s\n","# Domain Geometry ");
   fprintf(run_tcl_file, "%s\n","#---------------------------------------------------------");
   fprintf(run_tcl_file, "%s\n","pfset Geom.background.Lower.X                        $x0");
   fprintf(run_tcl_file, "%s\n","pfset Geom.background.Lower.Y                        $y0");
   fprintf(run_tcl_file, "%s\n","pfset Geom.background.Lower.Z                        $z0");
   fprintf(run_tcl_file, "%s\n","pfset Geom.background.Upper.X                        $xmax");
   fprintf(run_tcl_file, "%s\n","pfset Geom.background.Upper.Y                        $ymax");
   fprintf(run_tcl_file, "%s\n","pfset Geom.background.Upper.Z                        $zmax");
   fprintf(run_tcl_file, "%s\n","pfset Geom.background.Patches             \"x-lower x-upper y-lower y-upper z-lower z-upper\"");
   
   fprintf(run_tcl_file, "%s\n","#-----------------------------------------------------------------------------");
   fprintf(run_tcl_file, "%s\n","# Perm");
   fprintf(run_tcl_file, "%s\n","#-----------------------------------------------------------------------------");
   fprintf(run_tcl_file, "%s\n","pfset Geom.Perm.Names                 \"domain\"");
   fprintf(run_tcl_file, "%s\n","pfset Geom.domain.Perm.Type            Constant");
   fprintf(run_tcl_file, "%s\n","pfset Geom.domain.Perm.Value           0.001");

   fprintf(run_tcl_file, "%s\n","pfset Perm.TensorType               TensorByGeom");
   fprintf(run_tcl_file, "%s\n","pfset Geom.Perm.TensorByGeom.Names  \"domain\"");
   fprintf(run_tcl_file, "%s\n","pfset Geom.domain.Perm.TensorValX  1.0");
   fprintf(run_tcl_file, "%s\n","pfset Geom.domain.Perm.TensorValY  1.0");
   fprintf(run_tcl_file, "%s\n","pfset Geom.domain.Perm.TensorValZ  1.0");
   fprintf(run_tcl_file, "%s\n","#-----------------------------------------------------------------------------");
   fprintf(run_tcl_file, "%s\n","# Specific Storage");
   fprintf(run_tcl_file, "%s\n","#-----------------------------------------------------------------------------");
   fprintf(run_tcl_file, "%s\n","pfset SpecificStorage.Type            Constant");
   fprintf(run_tcl_file, "%s\n","pfset SpecificStorage.GeomNames       \"domain \"");
   fprintf(run_tcl_file, "%s\n","pfset Geom.domain.SpecificStorage.Value 1.0e-6");
   fprintf(run_tcl_file, "%s\n","#-----------------------------------------------------------------------------");
   fprintf(run_tcl_file, "%s\n","# Phases");
   fprintf(run_tcl_file, "%s\n","#-----------------------------------------------------------------------------");
   fprintf(run_tcl_file, "%s\n","pfset Phase.Names \"water\"");
   fprintf(run_tcl_file, "%s\n","pfset Phase.water.Density.Type	Constant");
   fprintf(run_tcl_file, "%s\n","pfset Phase.water.Density.Value	1.0");
   fprintf(run_tcl_file, "%s\n","pfset Phase.water.Viscosity.Type	Constant");
   fprintf(run_tcl_file, "%s\n","pfset Phase.water.Viscosity.Value	1.0");
   fprintf(run_tcl_file, "%s\n","#-----------------------------------------------------------------------------");
   fprintf(run_tcl_file, "%s\n","# Contaminants");
   fprintf(run_tcl_file, "%s\n","#-----------------------------------------------------------------------------");
   fprintf(run_tcl_file, "%s\n","pfset Contaminants.Names			\"\"");
   fprintf(run_tcl_file, "%s\n","#-----------------------------------------------------------------------------");
   fprintf(run_tcl_file, "%s\n","# Retardation");
   fprintf(run_tcl_file, "%s\n","#-----------------------------------------------------------------------------");
   fprintf(run_tcl_file, "%s\n","pfset Geom.Retardation.GeomNames           \"\"");
   fprintf(run_tcl_file, "%s\n","#-----------------------------------------------------------------------------");
   fprintf(run_tcl_file, "%s\n","# Gravity");
   fprintf(run_tcl_file, "%s\n","#-----------------------------------------------------------------------------");
   fprintf(run_tcl_file, "%s\n","pfset Gravity				1.0");
   fprintf(run_tcl_file, "%s\n","#-----------------------------------------------------------------------------");
   fprintf(run_tcl_file, "%s\n","# Porosity");
   fprintf(run_tcl_file, "%s\n","#-----------------------------------------------------------------------------");
   fprintf(run_tcl_file, "%s\n","pfset Geom.Porosity.GeomNames           \"domain\"");
   fprintf(run_tcl_file, "%s\n","pfset Geom.domain.Porosity.Type         Constant");
   fprintf(run_tcl_file, "%s\n","pfset Geom.domain.Porosity.Value        0.25");
   fprintf(run_tcl_file, "%s\n","#-----------------------------------------------------------------------------");
   fprintf(run_tcl_file, "%s\n","# Relative Permeability");
   fprintf(run_tcl_file, "%s\n","#-----------------------------------------------------------------------------");
   fprintf(run_tcl_file, "%s\n","pfset Phase.RelPerm.Type           VanGenuchten");
   fprintf(run_tcl_file, "%s\n","pfset Phase.RelPerm.GeomNames      \"domain\"");
   fprintf(run_tcl_file, "%s\n","pfset Geom.domain.RelPerm.Alpha    3.47");
   fprintf(run_tcl_file, "%s\n","pfset Geom.domain.RelPerm.N        7.24");
   fprintf(run_tcl_file, "%s\n","#-----------------------------------------------------------------------------");
   fprintf(run_tcl_file, "%s\n","# Saturation");
   fprintf(run_tcl_file, "%s\n","#-----------------------------------------------------------------------------");
   fprintf(run_tcl_file, "%s\n","pfset Phase.Saturation.Type              VanGenuchten");
   fprintf(run_tcl_file, "%s\n","pfset Phase.Saturation.GeomNames         \"domain\"");
   fprintf(run_tcl_file, "%s\n","pfset Geom.domain.Saturation.Alpha        3.47");
   fprintf(run_tcl_file, "%s\n","pfset Geom.domain.Saturation.N            7.24");
   fprintf(run_tcl_file, "%s\n","pfset Geom.domain.Saturation.SRes         0.12");
   fprintf(run_tcl_file, "%s\n","pfset Geom.domain.Saturation.SSat         1.0");
   fprintf(run_tcl_file, "%s\n","#-----------------------------------------------------------------------------");
   fprintf(run_tcl_file, "%s\n","# Mobility");
   fprintf(run_tcl_file, "%s\n","#-----------------------------------------------------------------------------");
   fprintf(run_tcl_file, "%s\n","pfset Phase.water.Mobility.Type        Constant");
   fprintf(run_tcl_file, "%s\n","pfset Phase.water.Mobility.Value       1.0");
   fprintf(run_tcl_file, "%s\n","#-----------------------------------------------------------------------------");
   fprintf(run_tcl_file, "%s\n","# Wells");
   fprintf(run_tcl_file, "%s\n","#-----------------------------------------------------------------------------");
   fprintf(run_tcl_file, "%s\n","pfset Wells.Names                           \"\"");
   fprintf(run_tcl_file, "%s\n","#-----------------------------------------------------------------------------");
   fprintf(run_tcl_file, "%s\n","# Setup timing info");
   fprintf(run_tcl_file, "%s\n","#-----------------------------------------------------------------------------");
   fprintf(run_tcl_file, "%s\n","# The UNITS on this simulation are HOURS");
   fprintf(run_tcl_file, "%s\n","pfset TimingInfo.BaseUnit        0.01");
   fprintf(run_tcl_file, "%s\n","pfset TimingInfo.StartCount      0");
   fprintf(run_tcl_file, "%s\n","pfset TimingInfo.StartTime       0.0");
   fprintf(run_tcl_file, "%s\n","pfset TimingInfo.StopTime        2.0");
   fprintf(run_tcl_file, "%s\n","pfset TimingInfo.DumpInterval    0.1");
   fprintf(run_tcl_file, "%s\n","pfset TimeStep.Type              Constant");
   fprintf(run_tcl_file, "%s\n","pfset TimeStep.Value             0.01");
   fprintf(run_tcl_file, "%s\n","#-----------------------------------------------------------------------------");
   fprintf(run_tcl_file, "%s\n","# Time Cycles");
   fprintf(run_tcl_file, "%s\n","#-----------------------------------------------------------------------------");
   fprintf(run_tcl_file, "%s\n","pfset Cycle.Names                       \"constant rainrec\"");
   fprintf(run_tcl_file, "%s\n","pfset Cycle.constant.Names              \"alltime\"");
   fprintf(run_tcl_file, "%s\n","pfset Cycle.constant.alltime.Length      200");
   fprintf(run_tcl_file, "%s\n","pfset Cycle.constant.Repeat             -1");
   
   fprintf(run_tcl_file, "%s\n","# rainfall and recession time periods are defined here");
   fprintf(run_tcl_file, "%s\n","pfset Cycle.rainrec.Names                  \"rain rec\"");
   fprintf(run_tcl_file, "%s\n","# nrain/BaseUnit");
   fprintf(run_tcl_file, "%s\n","pfset Cycle.rainrec.rain.Length            100");
   fprintf(run_tcl_file, "%s\n","# nrec/BaseUnit");
   fprintf(run_tcl_file, "%s\n","pfset Cycle.rainrec.rec.Length             100");
   fprintf(run_tcl_file, "%s\n","pfset Cycle.rainrec.Repeat                 -1");   
   fprintf(run_tcl_file, "%s\n","#-----------------------------------------------------------------------------");
   fprintf(run_tcl_file, "%s\n","# Boundary Conditions: Pressure");
   fprintf(run_tcl_file, "%s\n","#-----------------------------------------------------------------------------");

   fprintf(run_tcl_file, "%s%s%s\n","pfset BCPressure.PatchNames              \"",PatchesNames,"\"");
   int ct =0;
   char * pch;
   pch = strtok (PatchesNames2," ,.-");
   while (pch != NULL){
       ct++;
       if(ct == 2){
          fprintf(run_tcl_file, "%s%s%s\n","pfset Patch.",pch,".BCPressure.Type             OverlandFlow");
          fprintf(run_tcl_file, "%s%s%s\n","pfset Patch.",pch,".BCPressure.Cycle	      \"rainrec\"");
          fprintf(run_tcl_file, "%s%s%s\n","pfset Patch.",pch,".BCPressure.rain.Value	      -0.005");
          fprintf(run_tcl_file, "%s%s%s\n","pfset Patch.",pch,".BCPressure.rec.Value	      0.000");
       } 
       else {                 
          fprintf(run_tcl_file, "%s%s%s\n","pfset Patch.",pch,".BCPressure.Type		     FluxConst");
          fprintf(run_tcl_file, "%s%s%s\n","pfset Patch.",pch,".BCPressure.Cycle		      \"constant\"");
          fprintf(run_tcl_file, "%s%s%s\n","pfset Patch.",pch,".BCPressure.alltime.Value	      0.0");
       }       
       pch = strtok (NULL, " ,.-");       
   }

   
   fprintf(run_tcl_file, "%s\n","#-----------------------------------------------------------------------------");
   fprintf(run_tcl_file, "%s\n","# Topo slopes in x-direction");
   fprintf(run_tcl_file, "%s\n","#-----------------------------------------------------------------------------");
   if (param.Slope_x_Output->answer){
      fprintf(run_tcl_file, "%s%s%s\n","file copy -force \"",SlopeX_file,"\" slope_x.pfb");
      fprintf(run_tcl_file, "%s\n","pfset TopoSlopesX.Type                 \"PFBFile\"");
      fprintf(run_tcl_file, "%s\n","pfset TopoSlopesX.GeomNames            \"domain\"");
      fprintf(run_tcl_file, "%s\n","pfset TopoSlopesX.FileName              slope_x.pfb");
      fprintf(run_tcl_file, "%s\n","pfdist -nz 1 slope_x.pfb");
   }
   fprintf(run_tcl_file, "%s\n","#-----------------------------------------------------------------------------");
   fprintf(run_tcl_file, "%s\n","# Topo slopes in y-direction");
   fprintf(run_tcl_file, "%s\n","#-----------------------------------------------------------------------------");
   if (param.Slope_y_Output->answer){
      fprintf(run_tcl_file, "%s%s%s\n","file copy -force \"",SlopeY_file,"\" slope_y.pfb");
      fprintf(run_tcl_file, "%s\n","pfset TopoSlopesY.Type                 \"PFBFile\"");
      fprintf(run_tcl_file, "%s\n","pfset TopoSlopesY.GeomNames            \"domain\"");
      fprintf(run_tcl_file, "%s\n","pfset TopoSlopesY.FileName              slope_y.pfb");
      fprintf(run_tcl_file, "%s\n","pfdist -nz 1 slope_y.pfb");
   }
   if(G_strcasestr(grid_format,"Terrain Following Grid")){
      fprintf(run_tcl_file, "%s\n","# This define the DEM for the deflect terrain option (TFG only)");
      fprintf(run_tcl_file, "%s\n","file copy -force \"elevations.pfb\" \"Top_Surf.pfb\"");
      fprintf(run_tcl_file, "%s\n","pfset TopoSlopes.Elevation.FileName  \"Top_Surf.pfb\"");
      fprintf(run_tcl_file, "%s\n","pfdist -nz 1 \"Top_Surf.pfb\"");
   }
   
   fprintf(run_tcl_file, "%s\n","#-----------------------------------------------------------------------------");
   fprintf(run_tcl_file, "%s\n","# Mannings coefficient ");
   fprintf(run_tcl_file, "%s\n","#-----------------------------------------------------------------------------");
   
   double manning =3.3e-5;
   
   fprintf(run_tcl_file, "%s%lf\n","set mn ",manning);
   fprintf(run_tcl_file, "%s\n","pfset Mannings.Type \"Constant\"");   
   fprintf(run_tcl_file, "%s\n","pfset Mannings.GeomNames \"domain\"");
   fprintf(run_tcl_file, "%s\n","pfset Mannings.Geom.domain.Value $mn");
   fprintf(run_tcl_file, "%s\n","#-----------------------------------------------------------------------------");
   fprintf(run_tcl_file, "%s\n","# Phase sources:");
   fprintf(run_tcl_file, "%s\n","#-----------------------------------------------------------------------------");
   fprintf(run_tcl_file, "%s\n","pfset PhaseSources.water.Type                         Constant");
   fprintf(run_tcl_file, "%s\n","pfset PhaseSources.water.GeomNames                    domain");
   fprintf(run_tcl_file, "%s\n","pfset PhaseSources.water.Geom.domain.Value            0.0");
   fprintf(run_tcl_file, "%s\n","#-----------------------------------------------------------------------------");
   fprintf(run_tcl_file, "%s\n","# Exact solution specification for error calculations");
   fprintf(run_tcl_file, "%s\n","#-----------------------------------------------------------------------------");
   fprintf(run_tcl_file, "%s\n","pfset KnownSolution                                    NoKnownSolution");
   fprintf(run_tcl_file, "%s\n","#-----------------------------------------------------------------------------");
   fprintf(run_tcl_file, "%s\n","# Initial conditions: water pressure");
   fprintf(run_tcl_file, "%s\n","#-----------------------------------------------------------------------------");
   fprintf(run_tcl_file, "%s\n","pfset ICPressure.Type                                   HydroStaticPatch");
   fprintf(run_tcl_file, "%s\n","pfset ICPressure.GeomNames                              domain");
   fprintf(run_tcl_file, "%s\n","pfset Geom.domain.ICPressure.Value                      -1.0");
   fprintf(run_tcl_file, "%s\n","pfset Geom.domain.ICPressure.RefGeom                    domain");
   fprintf(run_tcl_file, "%s\n","pfset Geom.domain.ICPressure.RefPatch                   Top");
   fprintf(run_tcl_file, "%s\n","#-----------------------------------------------------------------------------");
   fprintf(run_tcl_file, "%s\n","# Set solver parameters");
   fprintf(run_tcl_file, "%s\n","#-----------------------------------------------------------------------------");
   if(G_strcasestr(grid_format,"Terrain Following Grid")){
      fprintf(run_tcl_file, "%s\n","pfset Solver.TerrainFollowingGrid 		             True");
   }
   fprintf(run_tcl_file, "%s\n","pfset Solver                                             Richards");
   fprintf(run_tcl_file, "%s\n","pfset Solver.MaxIter                                     2500000");
   fprintf(run_tcl_file, "%s\n","pfset OverlandFlowDiffusive                              0");
   fprintf(run_tcl_file, "%s\n","pfset Solver.Nonlinear.MaxIter                           1000");
   fprintf(run_tcl_file, "%s\n","pfset Solver.Nonlinear.ResidualTol                       1e-10");
   fprintf(run_tcl_file, "%s\n","pfset Solver.Nonlinear.EtaChoice                         Walker1");
   fprintf(run_tcl_file, "%s\n","pfset Solver.Nonlinear.EtaChoice                         EtaConstant");
   fprintf(run_tcl_file, "%s\n","pfset Solver.Nonlinear.EtaValue                          0.001");
   fprintf(run_tcl_file, "%s\n","pfset Solver.Nonlinear.UseJacobian                       False");
   fprintf(run_tcl_file, "%s\n","pfset Solver.Nonlinear.DerivativeEpsilon                 1e-16");
   fprintf(run_tcl_file, "%s\n","pfset Solver.Nonlinear.StepTol		             1e-30");
   fprintf(run_tcl_file, "%s\n","pfset Solver.Nonlinear.Globalization                     LineSearch");
   fprintf(run_tcl_file, "%s\n","pfset Solver.Linear.KrylovDimension                      50");
   fprintf(run_tcl_file, "%s\n","pfset Solver.Linear.MaxRestart                           3");
   fprintf(run_tcl_file, "%s\n","pfset Solver.Linear.Preconditioner                       PFMG");
   fprintf(run_tcl_file, "%s\n","pfset Solver.Linear.Preconditioner.PFMG.MaxIter           5");
   fprintf(run_tcl_file, "%s\n","pfset Solver.Linear.Preconditioner.PFMG.Smoother          RBGaussSeidelNonSymmetric");
   fprintf(run_tcl_file, "%s\n","pfset Solver.Linear.Preconditioner.PFMG.NumPreRelax       1");
   fprintf(run_tcl_file, "%s\n","pfset Solver.Linear.Preconditioner.PFMG.NumPostRelax      1");
   fprintf(run_tcl_file, "%s\n","pfset Solver.Drop                                       1E-20");
   fprintf(run_tcl_file, "%s\n","pfset Solver.AbsTol                                      1E-9");
   fprintf(run_tcl_file, "%s\n","pfset Solver.WriteSiloSubsurfData True");
   fprintf(run_tcl_file, "%s\n","pfset Solver.WriteSiloPressure True");
   fprintf(run_tcl_file, "%s\n","pfset Solver.WriteSiloSaturation True");
   fprintf(run_tcl_file, "%s\n","pfset Solver.WriteSiloConcentration True");

   fprintf(run_tcl_file, "%s\n","#-----------------------------------------------------------------------------");
   fprintf(run_tcl_file, "%s\n","# Run and Unload the ParFlow output files");
   fprintf(run_tcl_file, "%s\n","#-----------------------------------------------------------------------------");
   fprintf(run_tcl_file, "%s\n","pfrun $runname");
   fprintf(run_tcl_file, "%s\n","pfundist $runname");
   if(G_strcasestr(grid_format,"Terrain Following Grid")){
      fprintf(run_tcl_file, "%s\n","pfundist \"Top_Surf.pfb\"");
   }
   
   
   // ********************************************************
   //           Observation points definition
   // ********************************************************

    /* TODO: use r.cost method to create a list of observation points */
    if (outparam.coordsOut->answer) {
    
        fprintf(run_tcl_file, "%s\n","#-----------------------------------------------------------------------------");
        fprintf(run_tcl_file, "%s\n","# Overland flow observation points definition");
        fprintf(run_tcl_file, "%s\n","#-----------------------------------------------------------------------------");
        fprintf(run_tcl_file, "%s\n","#-----------------------------------------------------------------------------");
        fprintf(run_tcl_file, "%s\n","# Calculate overland flow at a point using Manning’s equation");
        fprintf(run_tcl_file, "%s\n","#-----------------------------------------------------------------------------");    
        fprintf(run_tcl_file, "%s\n","#Set the locations");
    
    	   // code based in r.path module 
        npoints = 0;
                
        int NfileResults = 20;
        int zlocp = 0;
        
    	fprintf(run_tcl_file, "%s\n","#Get the slope at the point");    	
	fprintf(run_tcl_file, "%s\n","set slopex   [pfload slope_x.pfb]");
	fprintf(run_tcl_file, "%s\n","set slopey   [pfload slope_y.pfb]");
    
	for (i = 0; outparam.coordsOut->answers[i] != NULL; i += 2) {
	    G_scan_easting(outparam.coordsOut->answers[i], &east, G_projection());
	    G_scan_northing(outparam.coordsOut->answers[i + 1], &north, G_projection());
	    obs_col = (int)Rast_easting_to_col(east, &window);
	    obs_row = nrows - (int)Rast_northing_to_row(north, &window);
	    

	    if (obs_row < 0 || obs_row > nrows ||
		obs_col < 0 || obs_col > ncols) {
		G_warning(_("Observation point %d is outside the current region"),
			  i + 1);
		continue;
	    }
	    npoints++;
	    
	    fprintf(run_tcl_file, "%s%d%s%d\n","set Xloc_",npoints," ", obs_col);
	    fprintf(run_tcl_file, "%s%d%s%d\n","set Yloc_",npoints," ", obs_row);
	    
	    if(G_strcasestr(grid_format,"Terrain Following Grid")){	    
	       fprintf(run_tcl_file, "%s%d%s%d\n","set Zloc_",npoints," ", NZ-1);	       
	    }
	    else{
	       zlocp = 0;
	       for(k=0;k<NZ;k++){
                   for (ii=0;ii<nrows;ii++){
	              for(jj=0;jj<ncols;jj++){	   
//	                count=count+1;
	                if(((k*zres)<=pfdem1[ii][jj]) && ((k*zres)>=(pfdem1[ii][jj]-zres))){         
	                   if((ii == obs_row) && (jj == obs_col)){
	                         zlocp = k;
	                   }
	                }	            
	              }                     	        
	           }
	        }
	        
	        fprintf(run_tcl_file, "%s%d%s%d\n","set Zloc_",npoints," ", zlocp);	    
	    }	    
	    fprintf(run_tcl_file, "%s\n","#This should be a z location on the surface of your domain");
	    fprintf(run_tcl_file, "%s\n","#Set the Mannings roughness coefficient");	    
	    fprintf(run_tcl_file, "%s%d%s\n","set n_",npoints," $mn");
	    
	    fprintf(run_tcl_file, "%s%d%s%d%s\n","set QT",npoints," [open \"OverladFlowP",npoints,".txt\" w+]");
	    
	    fprintf(run_tcl_file, "%s%d%s%d%s%d%s\n","set sx",npoints," [pfgetelt $slopex $Xloc_",npoints," $Yloc_",npoints," 0]");
	    fprintf(run_tcl_file, "%s%d%s%d%s%d%s\n","set sy",npoints," [pfgetelt $slopey $Xloc_",npoints," $Yloc_",npoints," 0]");
	    fprintf(run_tcl_file, "%s%d%s%d%s%d%s\n","set S",npoints," [expr ($sx",npoints,"**2+$sy",npoints,"**2)**0.5]");
	    
	    //have_points = 1;
	    //next_obs_pt = (struct point *)(G_malloc(sizeof(struct point)));

	    //next_obs_pt->row = obs_row;
	    //next_obs_pt->col = obs_col;
	    //next_obs_pt->value = npoints;
	    //next_obs_pt->next = head_obs_pt;
	    //head_obs_pt = next_obs_pt;
	}
	npoints=0;
	fprintf(run_tcl_file, "%s%d%s\n","for {set i 0} {$i <= ",NfileResults,"} {incr i} {");
	fprintf(run_tcl_file, "%s\n","#Get the pressure at the point");
	fprintf(run_tcl_file, "%s\n","set filename [format \"%s.out.press.%05d.pfb\" $runname $i]");
	fprintf(run_tcl_file, "%s\n","set press [pfload $filename]");		
	for (i = 0; outparam.coordsOut->answers[i] != NULL; i += 2) {
		npoints++;
		fprintf(run_tcl_file, "%s%d%s%d%s%d%s%d%s\n","set P",npoints," [pfgetelt $press $Xloc_",npoints," $Yloc_",npoints," $Zloc_",npoints,"]");
		fprintf(run_tcl_file, "%s\n","#If the pressure is less than zero set to zero");
		fprintf(run_tcl_file, "%s%d%s%d%s\n","if {$P",npoints," < 0} { set P",npoints," 0 }");
		fprintf(run_tcl_file, "%s%d%s%d%s%d%s%d%s\n","puts $QT",npoints," [expr ($dx/$n_",npoints,")*($S",npoints,"**0.5)*($P",npoints,"**(5./3.))]");
				
	}

	fprintf(run_tcl_file, "%s\n","}");
	npoints=0;
	for (i = 0; outparam.coordsOut->answers[i] != NULL; i += 2) {
		npoints++;
		fprintf(run_tcl_file, "%s%d\n","close $QT",npoints);
	}
		
    }      
   
   fclose(run_tcl_file); 
}

printf("X0: %lf\n",xcorner);
printf("Y0: %lf\n",ycorner);
printf("Z0: %lf\n",zcorner);

printf("Maximum value of the DEM: %lf\n",max_dem_value);
//printf("min_bottom_value: %lf\n",min_bottom_value);
                			
G_message(_("DONE"));
                                                        
//G_free(funcCall);       
G_free(out_dir); 
G_free(PathOutElev);
G_free(PathOutBot); 
G_free(PathOutTop);
G_free(PathOutMask);
//G_free(PathTclFile);        
G_free(name_dem_map); 
G_free(name_mask_map);
G_free(name_bottom_map);
G_free(PathSolidFile);
G_free(PathSolidVTKFile);
//G_free(name_slope_x_map);
//G_free(name_slope_y_map);

exit(EXIT_SUCCESS);
  
}

