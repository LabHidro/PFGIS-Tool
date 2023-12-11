.. PFGIS-Tool documentation master file, created by
   sphinx-quickstart on Monday, June 26, 2023.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to PFGIS-Tool User's Manual!
===================================

.. image:: Fig_PFGIS.png

**Tomas Carlotto** [1]_, **Julian Klaus** [2]_, **Pedro Luiz Borges Chaffe** [3]_

**PFGIS-Tool** provides agility and efficiency in the preprocessing of ParFlow model input files with a friendly graphical user interface. The new module performs: (i) topographic processing; (ii) creation of solid objects that delimit the computational domain and assign face identifiers in an orthogonal or terrain following grid; (iii) creation of files with the terrain slopes in the x and y directions; and (iv) creation of initial and boundary condition files. We provide enhancements that allow modeling of lake catchments with the insertion of multiple boundary conditions into the top face of the computational domain. Other new features allow the estimation of distributed subsurface depths in the catchment and the definition of the bottom boundary of the computational domain based on terrain elevation and slope.

A detailed description of PFGIS-Tool functionalities is available at: https://doi.org/10.1016/j.envsoft.2023.105824

.. note::
   This project is under active development.

.. toctree::
   :maxdepth: 2
   :caption: Contents:
   :numbered:

   start
   rparflow
   solids
   subsurfacedepth
   writepfb
   
.. [1]
   *Graduate Program in Environmental Engineering, Federal University of Santa Catarina, Florianópolis, Santa Catarina, Brazil.* thomas.carl@hotmail.com

.. [2]
   *Department of Geography, University of Bonn, Bonn, North Rhine-Westphalia, Germany.* jklaus@uni-bonn.de

.. [3]
   *Department of Sanitary and Environmental Engineering, Federal University of Santa Catarina, Florianópolis, Santa Catarina, Brazil* pedro.chaffe@ufsc.br

This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or any later version.
