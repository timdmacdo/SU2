/*!
 * \file wall_model.hpp
 * \brief Headers for the wall model functions for hom large eddy simulations.
 * \author E. van der Weide, T. Economon, P. Urbanczyk
 * \version 5.0.0 "Raven"
 *
 * SU2 Lead Developers: Dr. Francisco Palacios (Francisco.D.Palacios@boeing.com).
 *                      Dr. Thomas D. Economon (economon@stanford.edu).
 *
 * SU2 Developers: Prof. Juan J. Alonso's group at Stanford University.
 *                 Prof. Piero Colonna's group at Delft University of Technology.
 *                 Prof. Nicolas R. Gauger's group at Kaiserslautern University of Technology.
 *                 Prof. Alberto Guardone's group at Polytechnic University of Milan.
 *                 Prof. Rafael Palacios' group at Imperial College London.
 *                 Prof. Edwin van der Weide's group at the University of Twente.
 *                 Prof. Vincent Terrapon's group at the University of Liege.
 *
 * Copyright (C) 2012-2017 SU2, the open-source CFD code.
 *
 * SU2 is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 *
 * SU2 is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with SU2. If not, see <http://www.gnu.org/licenses/>.
 */

#pragma once

#include "../../Common/include/mpi_structure.hpp"

#include <iostream>
#include <cmath>

/* Classes from other parts of the code need to be forward declared */
class CSolver;
class CConfig;
class CGeometry;

using namespace std;

/*!
 * \class CSGSModel
 * \brief Base class for defining the LES subgrid scale model.
 * \author: E. van der Weide, T. Economon, P. Urbanczyk
 * \version 5.0.0 "Raven"
 */
class CWallModel {

public:
  /*!
   * \brief Constructor of the class.
   */
  CWallModel(void);

  /*!
   * \brief Destructor of the class.
   */
  virtual ~CWallModel(void);

  virtual void ComputeWallShear();

  virtual void ComputeWallHeatFlux();

  virtual void Initialize(CGeometry *geometry, CConfig *config, CSolver *solver, unsigned short int thisMarker);

  virtual void SetUpExchange();

  virtual void SetPoints();

};

class CWallModel1DEQ : public CWallModel {

public:

  /*!
   * \brief Constructor of the class.
   */
  CWallModel1DEQ(void);

  /*!
   * \brief Destructor of the class.
   */
  ~CWallModel1DEQ(void);

  void ComputeWallShear();

  void ComputeWallHeatFlux();

  void Initialize(CGeometry *geometry, CConfig *config, CSolver *solver, unsigned short int thisMarker);

  void SetUpExchange();

  void SetPoints();

};


