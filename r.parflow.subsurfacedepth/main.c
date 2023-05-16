
/****************************************************************************
 *
 * MODULE:    r.parflow.subsurfacedepth
 * 
 * AUTHOR:    Tomas Carlotto ----------------- thomas.carl@hotmail.com
 *                               
 *               
 * PURPOSE:      Estimate subsurface depths and bottom surface from topography to be used in the Parflow hydrological model.
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
    
    struct Option *SoilDepthMethods;
    struct Option *ElevationMap;
    struct Option *DminParam;    
    struct Option *DmaxParam;           
    struct Option *SlopeMethod;        
    struct Option *SoilDepthOut;    
    struct Option *BottomSurfaceOut;
    struct Option *CatchmentMask;
    
} 

paramType;
paramType param, outparam;		//Parameters 

void set_parameters(){
   
    // ============================================================
    // ================== Soil depth methods ======================
    
    param.SoilDepthMethods = G_define_option();
    param.SoilDepthMethods->key = "method";
    param.SoilDepthMethods->type = TYPE_STRING;
    param.SoilDepthMethods->required = YES;
    param.SoilDepthMethods->label = _("Methods for estimating subsurface depths");
    param.SoilDepthMethods->description = _("Method for estimating subsurface depths based on topography.");
    param.SoilDepthMethods->options = "Elevation based,Slope based";
    param.SoilDepthMethods->answer = "Elevation based";
    param.SoilDepthMethods->guisection = _("Required");
    
    param.ElevationMap = G_define_standard_option(G_OPT_R_INPUT);
    param.ElevationMap->key = "pftop";
    param.ElevationMap->label = _("Digital elevation model");
    param.ElevationMap->required = YES;
    param.ElevationMap->guisection = _("Required");
    
    param.CatchmentMask = G_define_standard_option(G_OPT_R_INPUT);
    param.CatchmentMask->key = "basin_msk";
    param.CatchmentMask->label = _("Name of catchment mask");
    param.CatchmentMask->description = _("The catchment mask ensures that the area considered for calculations does not exceed the limits of the catchment. The mask is composed of values 1 in the catchment and 0 in the external areas.");
    param.CatchmentMask->required = YES;
    param.CatchmentMask->guisection = _("Required");
    
    param.DminParam = G_define_option();
    param.DminParam->key = "d_min";
    param.DminParam->label = _("Minimum depth [meters]");
    param.DminParam->required = YES;   
    param.DminParam->type = TYPE_DOUBLE;
    param.DminParam->guisection = _("Required"); 
    
    param.DmaxParam = G_define_option();
    param.DmaxParam->key = "d_max";
    param.DmaxParam->label = _("Maximum depth [meters]");
    param.DmaxParam->required = YES;   
    param.DmaxParam->type = TYPE_DOUBLE;
    param.DmaxParam->guisection = _("Required"); 
    
    param.SlopeMethod = G_define_option();
    param.SlopeMethod->key = "slope_method";
    param.SlopeMethod->type = TYPE_STRING;
    param.SlopeMethod->required = NO;
    param.SlopeMethod->label = _("Slope methods [Only for slope based method]");
    param.SlopeMethod->description = _(" ");
    param.SlopeMethod->options = "Upwind scheme, Central difference";
    param.SlopeMethod->answer = "Upwind scheme";
    param.SlopeMethod->guisection = _("Required"); 
    
    param.SoilDepthOut= G_define_standard_option(G_OPT_R_OUTPUT);
    param.SoilDepthOut->key = "depth_out";
    param.SoilDepthOut->label = _("Name for the resulting subsurface depth map");
    param.SoilDepthOut->required = YES;
    param.SoilDepthOut->guisection = _("Required");
    
    param.BottomSurfaceOut= G_define_standard_option(G_OPT_R_OUTPUT);
    param.BottomSurfaceOut->key = "bottom_out";
    param.BottomSurfaceOut->label = _("Name for the resulting bottom surface map");
    param.BottomSurfaceOut->required = YES;
    param.BottomSurfaceOut->guisection = _("Required");

  }

int main(int argc, char *argv[]){
    
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
module->description = _("Tool to estimate subsurface depths and bottom surface from topography to be used in the Parflow hydrological model");
    
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
           
RASTER_MAP_TYPE data_type_dem_map;
//RASTER_MAP_TYPE data_type_bottom_map;

char *name_dem_map = param.ElevationMap->answer;        // Digital Elevation Map name ]
char *soil_depth_method = param.SoilDepthMethods->answer;
char *slope_method = param.SlopeMethod->answer;

data_type_dem_map = MAPdatatype(name_dem_map);
//data_type_bottom_map = MAPdatatype(name_bottom_map);

double Dmax, Dmin;

sscanf(param.DminParam->answer, "%lf", &Dmin);
sscanf(param.DmaxParam->answer, "%lf", &Dmax);

RASTER_MAP_TYPE data_type_CatchmentMask;
char *name_CatchmentMask = param.CatchmentMask->answer;
data_type_CatchmentMask = MAPdatatype(name_CatchmentMask);

int i,j,is,in_i;

double **pfCatchmentMask1= (double**)malloc(nrows*sizeof(double*));
double **pfdem1 = (double**)malloc(nrows*sizeof(double*));
//double **pfbottom1= (double**)malloc(nrows*sizeof(double*));
double **pfSlopex1= (double**)malloc(nrows*sizeof(double*));
double **pfSlopey1= (double**)malloc(nrows*sizeof(double*));
double **pfSlope= (double**)malloc(nrows*sizeof(double*));

for (i = 0; i < nrows; i++){                      
    pfCatchmentMask1[i] = (double*)malloc(ncols*sizeof(double));
    pfdem1[i] = (double*)malloc(ncols*sizeof(double));
//    pfbottom1[i] = (double*)malloc(ncols*sizeof(double));
    pfSlopex1[i] = (double*)malloc(ncols*sizeof(double));
    pfSlopey1[i] = (double*)malloc(ncols*sizeof(double));
    pfSlope[i] = (double*)malloc(ncols*sizeof(double));    
}

G_message(_("START"));
            
//========= LOAD DEM DATA ==========================
if(data_type_dem_map == CELL_TYPE){
    G_message(_("LOAD DEM DATA"));
    in_i = nrows-1;    

    CELL** pfdem = (CELL**) G_malloc( nrows*sizeof(CELL*));    
    for (i = 0; i < nrows; i++){                    
         pfdem[i] = (CELL*) G_malloc( ncols * sizeof(CELL));
    }
    //    
    pfdem = open_raster_C_variable(name_dem_map);

    for (i=0;i<nrows;i++){
       for (j=0;j<ncols;j++){   
          pfdem1[i][j] = pfdem[in_i-i][j];   
          
       }
    }
    close_raster_C_variable(pfdem);
}
else if(data_type_dem_map == FCELL_TYPE){
    G_message(_("LOAD DEM DATA"));
    in_i = nrows-1;

    FCELL** pfdem = (FCELL**) G_malloc( nrows*sizeof(FCELL*));    
    for (i = 0; i < nrows; i++){                    
         pfdem[i] = (FCELL*) G_malloc( ncols * sizeof(FCELL));
    }
    
    pfdem = open_raster_F_variable(name_dem_map);

    for (i=0;i<nrows;i++){
       for (j=0;j<ncols;j++){   
          pfdem1[i][j] = pfdem[in_i-i][j];
          
       }
    }
    close_raster_F_variable(pfdem);
}
else if(data_type_dem_map == DCELL_TYPE){
    G_message(_("LOAD DEM DATA"));
    in_i = nrows-1; 
    DCELL** pfdem = (DCELL**) G_malloc( nrows*sizeof(DCELL*));    
    for (i = 0; i < nrows; i++){                    
         pfdem[i] = (DCELL*) G_malloc( ncols * sizeof(DCELL));
    }    
    
    pfdem = open_raster_D_variable(name_dem_map);
 
    for (i=0;i<nrows;i++){
        for (j=0;j<ncols;j++){   
          pfdem1[i][j] = pfdem[in_i-i][j];

        }
    }
           
     close_raster_D_variable(pfdem);
}

//========= LOAD CATCHMENT MASK DATA ==========================
if(data_type_CatchmentMask == CELL_TYPE){
    in_i = nrows-1;    
    CELL** pfCatchmentMask;
    pfCatchmentMask = open_raster_C_variable(name_CatchmentMask);
    for (i=0;i<nrows;i++){
        for (j=0;j<ncols;j++){   
            pfCatchmentMask1[i][j] = pfCatchmentMask[in_i-i][j];
        }
    }
    
    // Add mask contour

    for (i=0;i<nrows;i++){
        for (j=0;j<ncols;j++){   
          if(pfCatchmentMask[(in_i-i)][j]==1){
            if(pfCatchmentMask[(in_i-i)+1][j]==0){  
                 pfCatchmentMask1[i-1][j] = 1;
            }
            if(pfCatchmentMask[(in_i-i)-1][j]==0){  
                 pfCatchmentMask1[i+1][j] = 1;
            }
            if(pfCatchmentMask[(in_i-i)][j+1]==0){  
                 pfCatchmentMask1[i][j+1] = 1;
            }
            if(pfCatchmentMask[(in_i-i)][j-1]==0){  
                 pfCatchmentMask1[i][j-1] = 1;
            }
            // corners
            if(pfCatchmentMask[(in_i-i)+1][j+1]==0){  
                 pfCatchmentMask1[i-1][j+1] = 1;
            }
            if(pfCatchmentMask[(in_i-i)-1][j+1]==0){  
                 pfCatchmentMask1[i+1][j+1] = 1;
            }
            if(pfCatchmentMask[(in_i-i)+1][j-1]==0){  
                 pfCatchmentMask1[i-1][j-1] = 1;
            }
            if(pfCatchmentMask[(in_i-i)-1][j-1]==0){  
                 pfCatchmentMask1[i+1][j-1] = 1;
            }
          }
        }
    }
    
    close_raster_C_variable(pfCatchmentMask);
}
else if(data_type_CatchmentMask == FCELL_TYPE){
    in_i = nrows-1;
    FCELL** pfCatchmentMask;
    pfCatchmentMask = open_raster_F_variable(name_CatchmentMask);
    for (i=0;i<nrows;i++){
        for (j=0;j<ncols;j++){   
            pfCatchmentMask1[i][j] = pfCatchmentMask[in_i-i][j];
        }
    }
    
    // Add mask contour

    for (i=0;i<nrows;i++){
        for (j=0;j<ncols;j++){   
          if(pfCatchmentMask[(in_i-i)][j]==1){
            if(pfCatchmentMask[(in_i-i)+1][j]==0){  
                 pfCatchmentMask1[i-1][j] = 1;
            }
            if(pfCatchmentMask[(in_i-i)-1][j]==0){  
                 pfCatchmentMask1[i+1][j] = 1;
            }
            if(pfCatchmentMask[(in_i-i)][j+1]==0){  
                 pfCatchmentMask1[i][j+1] = 1;
            }
            if(pfCatchmentMask[(in_i-i)][j-1]==0){  
                 pfCatchmentMask1[i][j-1] = 1;
            }
            // corners
            if(pfCatchmentMask[(in_i-i)+1][j+1]==0){  
                 pfCatchmentMask1[i-1][j+1] = 1;
            }
            if(pfCatchmentMask[(in_i-i)-1][j+1]==0){  
                 pfCatchmentMask1[i+1][j+1] = 1;
            }
            if(pfCatchmentMask[(in_i-i)+1][j-1]==0){  
                 pfCatchmentMask1[i-1][j-1] = 1;
            }
            if(pfCatchmentMask[(in_i-i)-1][j-1]==0){  
                 pfCatchmentMask1[i+1][j-1] = 1;
            }
          }
        }
    }
    
    close_raster_F_variable(pfCatchmentMask);
}
else if(data_type_CatchmentMask == DCELL_TYPE){
    in_i = nrows-1;
    DCELL** pfCatchmentMask;
    pfCatchmentMask = open_raster_D_variable(name_CatchmentMask);
    for (i=0;i<nrows;i++){
        for (j=0;j<ncols;j++){   
            pfCatchmentMask1[i][j] = pfCatchmentMask[in_i-i][j];
        }
    }
    
// Add mask contour

    for (i=0;i<nrows;i++){
        for (j=0;j<ncols;j++){   
          if(pfCatchmentMask[(in_i-i)][j]==1){
            if(pfCatchmentMask[(in_i-i)+1][j]==0){  
                 pfCatchmentMask1[i-1][j] = 1;
            }
            if(pfCatchmentMask[(in_i-i)-1][j]==0){  
                 pfCatchmentMask1[i+1][j] = 1;
            }
            if(pfCatchmentMask[(in_i-i)][j+1]==0){  
                 pfCatchmentMask1[i][j+1] = 1;
            }
            if(pfCatchmentMask[(in_i-i)][j-1]==0){  
                 pfCatchmentMask1[i][j-1] = 1;
            }
            // corners
            if(pfCatchmentMask[(in_i-i)+1][j+1]==0){  
                 pfCatchmentMask1[i-1][j+1] = 1;
            }
            if(pfCatchmentMask[(in_i-i)-1][j+1]==0){  
                 pfCatchmentMask1[i+1][j+1] = 1;
            }
            if(pfCatchmentMask[(in_i-i)+1][j-1]==0){  
                 pfCatchmentMask1[i-1][j-1] = 1;
            }
            if(pfCatchmentMask[(in_i-i)-1][j-1]==0){  
                 pfCatchmentMask1[i+1][j-1] = 1;
            }
          }
        }
    }
    
    close_raster_D_variable(pfCatchmentMask);
}   



   
// ======================================================
// ======Calculation of slopes in the x and y direction 

int isl; // x direction
int jsl; // y direction
   
  if(G_strcasestr(soil_depth_method,"Slope based")){ 
   //else if(G_strcasestr(slope_method,"Upwind scheme")){
   if(G_strcasestr(slope_method,"Upwind scheme")){
   // Based in:
   // https://github.com/parflow/parflow/blob/master/pftools/toposlopes.c
   
   G_message(_("UPWIND SCHEME"));
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
           
           pfSlope[jsl][isl] = sqrt((pfSlopex1[jsl][isl] * pfSlopex1[jsl][isl]) + (pfSlopey1[jsl][isl]*pfSlopey1[jsl][isl]));
           
      }
   }
   }  
   else if (G_strcasestr(slope_method,"Central difference")){
   G_message(_("CENTRAL DIFFERENCE METHOD"));
   
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
       
       pfSlope[jsl][isl] = sqrt((pfSlopex1[jsl][isl] * pfSlopex1[jsl][isl]) + (pfSlopey1[jsl][isl]*pfSlopey1[jsl][isl]));
       
     }
   }
   
   }
  }  

// DEM Max and Min values

double ELmax=-10000000000, ELmin= 10000000000;//pfdem1[0][0]  pfbottom1[0][0];
double SLmax = -100, SLmin = 100;
 //   int i,j;

    for (i=0; i<nrows; i++){
	for (j=0; j<ncols; j++){	
          
          if(pfCatchmentMask1[i][j] == 1){
          //if(G_strcasestr(soil_depth_method,"Elevation based")){  
          
            if((pfdem1[i][j]> ELmax)){
	      ELmax= pfdem1[i][j];
	    }

	    if((pfdem1[i][j]<ELmin)){
              ELmin= pfdem1[i][j];
            }
            			
          //} 
          //else {
          
            if((pfSlope[i][j]> SLmax)){
	      SLmax= pfSlope[i][j];
	    }

	    if((pfSlope[i][j]<SLmin)){
              SLmin= pfSlope[i][j];
            }
            
          //}
          
	  }
	}
	
    }         
      
// =====================================
// ============== Soil depth ===========

DCELL** pfbottom = (DCELL**) G_malloc( nrows*sizeof(DCELL*));
DCELL** SoilDepth = (DCELL**) G_malloc( nrows*sizeof(DCELL*));    
for (i = 0; i < nrows; i++){                    
    pfbottom[i] = (DCELL*) G_malloc( ncols * sizeof(DCELL));
    SoilDepth[i] = (DCELL*) G_malloc( ncols * sizeof(DCELL));
} 

if(G_strcasestr(soil_depth_method,"Elevation based")){
   G_message(_("METHOD: ELEVATION BASED"));
 
// Dmax é a profundidade máxima na bacia (em metros)
// Dmin é a profundidade minima na bacia (em metros)
// pfdem1(i,j) é a elevação da célula; 
// ELmax é a elevação máxima na bacia;
// ELmin é a elevação mínima na bacia;

//   Dmax = ELmin - Dmax;
//   Dmin = ELmax - Dmin;

   for (i=0; i<nrows; i++){
      for (j=0; j<ncols; j++){
         if(pfCatchmentMask1[i][j] == 1){
            SoilDepth[(nrows-1)-i][j] = (Dmax-((Dmax-Dmin)/(ELmax-ELmin))*(pfdem1[i][j]-ELmin));
	    pfbottom[(nrows-1)-i][j] = pfdem1[i][j] - (Dmax-((Dmax-Dmin)/(ELmax-ELmin))*(pfdem1[i][j]-ELmin));
	 }
	 else{
	    SoilDepth[(nrows-1)-i][j] = Dmin;
	    pfbottom[(nrows-1)-i][j] = pfdem1[i][j]-Dmin;	 
	 }    
      }	
   } 
	
}
else{
   G_message(_("METHOD: SLOPE BASED"));

   //Dmax = ELmin - Dmax; 
   //Dmin = ELmax - Dmin; 

   for (i=0; i<nrows; i++){
      for (j=0; j<ncols; j++){
         if(pfCatchmentMask1[i][j] == 1){
	    SoilDepth[(nrows-1)-i][j] = (Dmax-((Dmax-Dmin)/(SLmax-SLmin))*(pfSlope[i][j]-SLmin)); 
	    pfbottom[(nrows-1)-i][j] = pfdem1[i][j]-(Dmax-((Dmax-Dmin)/(SLmax-SLmin))*(pfSlope[i][j]-SLmin));   
	 }
	 else{
	    SoilDepth[(nrows-1)-i][j] = Dmin;
	    pfbottom[(nrows-1)-i][j] = pfdem1[i][j]-Dmin;
	 }
      }	
   } 
}

char* out_soil_depth_name = param.SoilDepthOut->answer;
char* out_bottom_name = param.BottomSurfaceOut->answer;
save_raster_D_variable(SoilDepth,out_soil_depth_name);        
save_raster_D_variable(pfbottom,out_bottom_name); 
close_raster_D_variable(SoilDepth);
close_raster_D_variable(pfbottom);
                			
G_message(_("DONE"));
                                                              
G_free(name_dem_map);
G_free(out_soil_depth_name); 

exit(EXIT_SUCCESS);
  
}

