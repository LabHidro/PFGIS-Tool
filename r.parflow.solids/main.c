
/****************************************************************************
 *
 * MODULE:    r.parflow.solids
 * 
 * AUTHOR:    Tomas Carlotto ----------------- thomas.carl@hotmail.com
 *                               
 *               
 * PURPOSE:      Create solid files (.pfsol) to define the computational domain or to specify regions with different properties in the Parflow hydrologic model
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
    
    struct Option *MaskMap;
    struct Option *ElevationMap;
    struct Option *BottomMap;
    struct Option *GridFormat;
        
    struct Option *ElevOutput;
    struct Option *PathOutput;
    struct Option *Solid_pfsolOutput;
    struct Option *Solid_vtkOutput;
    
    struct Flag *Out_vtk_flag; 

} 

paramType;
paramType param, outparam;		//Parameters 

void set_parameters(){
//---------------------------------------------------------------------------------------------
//---------------------------------------------------------------------------------------------
//----INPUTS OF SIMULATION MODULE -------------
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
    param.MaskMap->guisection = _("Build solid file");
    
    param.ElevationMap = G_define_standard_option(G_OPT_R_INPUT);
    param.ElevationMap->key = "pftop";
    param.ElevationMap->label = _("Upper boundary raster map name ");
    param.ElevationMap->description = _("The upper boundary can be the digital elevation model.");
    param.ElevationMap->required = YES;
    param.ElevationMap->guisection = _("Build solid file");
    
    param.BottomMap = G_define_standard_option(G_OPT_R_INPUT);
    param.BottomMap->key = "pfbot";
    param.BottomMap->label = _("Bottom boundary raster map name");
    param.BottomMap->description = _("Raster that defines the bottom of the computational domain.");
    param.BottomMap->required = YES;
    param.BottomMap->guisection = _("Build solid file");
    
    param.Solid_pfsolOutput= G_define_standard_option(G_OPT_R_OUTPUT);
    param.Solid_pfsolOutput->key = "out_pfsol";
    param.Solid_pfsolOutput->label = _("Name for the output solid file (.pfsol)");
    param.Solid_pfsolOutput->description = _(" ");
    param.Solid_pfsolOutput->required = YES;
    param.Solid_pfsolOutput->guisection = _("Build solid file");
                       
    // ==========================================================================
    // ================== Output definition =====================================
    /*
    param.Out_vtk_flag = G_define_flag();
    param.Out_vtk_flag->key = 's';
    param.Out_vtk_flag->label = _("Create output files in .vtk format.");
    param.Out_vtk_flag->description = _("When this option is checked, in addition to files in .pfb format, files in .vtk format will also be generated.");
    param.Out_vtk_flag->guisection = _("Output");
    */     
    outparam.PathOutput = G_define_standard_option(G_OPT_M_DIR);//G_OPT_F_OUTPUT);
    outparam.PathOutput->key = "out_dir";
    outparam.PathOutput->type = TYPE_STRING;
    outparam.PathOutput->description = _("Name for output directory");
    outparam.PathOutput->required = YES;
    outparam.PathOutput->guisection = _("Build solid file");    
        
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
    module->description = _("Preprocessing tool to create solid files (.pfsol) for the ParFlow hydrologic model");
    
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
RASTER_MAP_TYPE data_type_mask_map;
RASTER_MAP_TYPE data_type_bottom_map;

char *out_dir = outparam.PathOutput->answer; // Save Path

char *name_dem_map = param.ElevationMap->answer;        // Digital Elevation Map name 
char *name_mask_map = param.MaskMap->answer;            // mask map name
char *name_bottom_map = param.BottomMap->answer;        // Bottom limit map name
char *out_pfsol_name = param.Solid_pfsolOutput->answer;  // Solid file name
char *grid_format = param.GridFormat->answer;            // Grid format

data_type_dem_map = MAPdatatype(name_dem_map);
data_type_mask_map = MAPdatatype(name_mask_map);
data_type_bottom_map = MAPdatatype(name_bottom_map);

int i,j,is,in_i;

double **pfdem1 = (double**)malloc(nrows*sizeof(double*));
double **pfmask1= (double**)malloc(nrows*sizeof(double*));
double **pfbottom1= (double**)malloc(nrows*sizeof(double*));

for (i = 0; i < nrows; i++){                      
    pfdem1[i] = (double*)malloc(ncols*sizeof(double));
    pfmask1[i] = (double*)malloc(ncols*sizeof(double));
    pfbottom1[i] = (double*)malloc(ncols*sizeof(double));
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
     
  
// ======================================================
// ============== Make Terrain Following Grid ===========
double Z0;
if(G_strcasestr(grid_format,"Terrain Following Grid")){
//double s_depth=0.00;
G_message(_("GRID FORMAT - TERRAIN FOLLOWING GRID"));
                       
    for (i=0;i<nrows;i++){
        for (j=0;j<ncols;j++){  
            
            //if((s_depth==0.00) && (pfmask1[i][j]==1)){
            //    s_depth = (pfdem1[i][j] - pfbottom1[i][j]);            
            //}
            
            if((pfmask1[i][j]==1)&&((pfdem1[i][j] - pfbottom1[i][j])>max_subsurfacedepth)){
                max_subsurfacedepth = (pfdem1[i][j] - pfbottom1[i][j]);            
            }
                             
          pfbottom1[i][j] = (max_dem_value) + (pfbottom1[i][j]-pfdem1[i][j]); 
          pfdem1[i][j] = max_dem_value;
          
        }
    }       
    Z0 = max_dem_value-max_subsurfacedepth;//s_depth;      
}
else{
   G_message(_("GRID FORMAT - ORTHOGONAL GRID"));
   Z0 = min_bottom_value;
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
PathSolidFile = concat(PathSolidFile,out_pfsol_name);
PathSolidFile = concat(PathSolidFile,".pfsol");
	
FILE *file_solidfile;
file_solidfile = fopen(PathSolidFile, "w");	
	
char *PathSolidVTKFile;	
PathSolidVTKFile = concat(out_dir,"/");
PathSolidVTKFile = concat(PathSolidVTKFile,out_pfsol_name);
PathSolidVTKFile = concat(PathSolidVTKFile,".vtk");
	
FILE *file_solidfileVTK;
file_solidfileVTK = fopen(PathSolidVTKFile, "w");			

int status;
		
status = MakePatchySolid(ncols,pfmask1,pfdem1,pfbottom1,nrows,1,xcorner,ycorner, xres, yres,0,0,file_solidfile,file_solidfileVTK);

printf("X0 = %lf\n",xcorner);
printf("Y0 = %lf\n",ycorner);
printf("Z0 = %lf\n",Z0);
printf("top maximum value = %lf\n",max_dem_value);
//printf("bottom minimum value = %lf\n",min_bottom_value); 
                			
G_message(_("DONE"));
                                                               
G_free(out_dir); 
     
G_free(name_dem_map); 
G_free(name_mask_map);
G_free(name_bottom_map);
G_free(PathSolidFile);
G_free(PathSolidVTKFile);

exit(EXIT_SUCCESS);
  
}

