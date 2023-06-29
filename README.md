# PFGIS-Tool
Preprocessing tool for the ParFlow hydrological model

Authors: Tomas Carlotto, Julian Klaus and Pedro Luiz Borges Chaffe

Contact email: thomas.carl@hotmail.com

ParFlow is an open-source high performance parallel hydrological model designed for the integrated modeling of the hydrological cycle from the hillslope to continental scales. Due to its complexity, ParFlow requires intrincate data preprocessing steps that might hinder its application by the wider modeling community. Here, we present the PFGIS-Tool module that provides agility and efficiency in the preprocessing of ParFlow model input files with a friendly graphical user interface. The new module performs: (i) topographic processing; (ii) creation of solid objects that delimit the computational domain and assign face identifiers in an orthogonal or terrain following grid; (iii) creation of files with the terrain slopes in the x and y directions; and (iv) creation of initial and boundary condition files. We provide enhancements that allow modeling of lake catchments with the insertion of multiple boundary conditions into the top face of the computational domain. Other new features allow the estimation of distributed subsurface depths in the catchment and the definition of the bottom boundary of the computational domain based on terrain elevation and slope.

A step by step for installing the PFGIS-Tool are available in the file Install-PFGIS-Tool.txt or in the link: https://pfgis-tool.readthedocs.io/en/latest/start.html

The PFGIS-Tool user manual is available at: https://pfgis-tool.readthedocs.io/en/latest/index.html (Under development).
