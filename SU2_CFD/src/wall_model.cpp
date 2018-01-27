/*!
 * \file wall_model.cpp
 * \brief Function for the wall model functions for hom large eddy simulations.
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

#include "../include/wall_model.hpp"

CWallModel::CWallModel(void){}

CWallModel::~CWallModel(void){}

void CWallModel::ComputeWallShear(void){}

void CWallModel::ComputeWallHeatFlux(void){}

void CWallModel::Initialize(CGeometry *geometry, CConfig *config, CSolver *solver, unsigned short int thisMarker){}

void CWallModel::SetUpExchange(void){}

void CWallModel::SetPoints(void){}

CWallModel1DEQ::CWallModel1DEQ(void) : CWallModel(){}

CWallModel1DEQ::~CWallModel1DEQ(void){}

void CWallModel1DEQ::ComputeWallShear(void){}

void CWallModel1DEQ::ComputeWallHeatFlux(void){}

void CWallModel1DEQ::Initialize(CGeometry *geometry, CConfig *config, CSolver *solver, unsigned short int thisMarker){}

void CWallModel1DEQ::SetUpExchange(void){}

void CWallModel1DEQ::SetPoints(void){}
