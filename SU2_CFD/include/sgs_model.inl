/*!
 * \file sgs_model.inl
 * \brief In-Line subroutines of the <i>sgs_model.hpp</i> file.
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

inline CSGSModel::CSGSModel(void){}
inline CSGSModel::~CSGSModel(void){}

inline su2double CSGSModel::ComputeEddyViscosity_2D(const su2double rho,
                                                   const su2double dudx,
                                                   const su2double dudy,
                                                   const su2double dvdx,
                                                   const su2double dvdy,
                                                   const su2double lenScale,
                                                   const su2double distToWall) {
  return 0.0;
}

inline su2double CSGSModel::ComputeEddyViscosity_3D(const su2double rho,
                                                    const su2double dudx,
                                                    const su2double dudy,
                                                    const su2double dudz,
                                                    const su2double dvdx,
                                                    const su2double dvdy,
                                                    const su2double dvdz,
                                                    const su2double dwdx,
                                                    const su2double dwdy,
                                                    const su2double dwdz,
                                                    const su2double lenScale,
                                                    const su2double distToWall) {
  return 0.0;
}

inline void CSGSModel::ComputeGradEddyViscosity_2D(const su2double rho,
                                                   const su2double drhodx,
                                                   const su2double drhody,
                                                   const su2double dudx,
                                                   const su2double dudy,
                                                   const su2double dvdx,
                                                   const su2double dvdy,
                                                   const su2double d2udx2,
                                                   const su2double d2udy2,
                                                   const su2double d2udxdy,
                                                   const su2double d2vdx2,
                                                   const su2double d2vdy2,
                                                   const su2double d2vdxdy,
                                                   const su2double lenScale,
                                                   const su2double distToWall,
                                                         su2double &dMuTdx,
                                                         su2double &dMuTdy) {
  dMuTdx = dMuTdy = 0.0;
}

inline void CSGSModel::ComputeGradEddyViscosity_3D(const su2double rho,
                                                   const su2double drhodx,
                                                   const su2double drhody,
                                                   const su2double drhodz,
                                                   const su2double dudx,
                                                   const su2double dudy,
                                                   const su2double dudz,
                                                   const su2double dvdx,
                                                   const su2double dvdy,
                                                   const su2double dvdz,
                                                   const su2double dwdx,
                                                   const su2double dwdy,
                                                   const su2double dwdz,
                                                   const su2double d2udx2,
                                                   const su2double d2udy2,
                                                   const su2double d2udz2,
                                                   const su2double d2udxdy,
                                                   const su2double d2udxdz,
                                                   const su2double d2udydz,
                                                   const su2double d2vdx2,
                                                   const su2double d2vdy2,
                                                   const su2double d2vdz2,
                                                   const su2double d2vdxdy,
                                                   const su2double d2vdxdz,
                                                   const su2double d2vdydz,
                                                   const su2double d2wdx2,
                                                   const su2double d2wdy2,
                                                   const su2double d2wdz2,
                                                   const su2double d2wdxdy,
                                                   const su2double d2wdxdz,
                                                   const su2double d2wdydz,
                                                   const su2double lenScale,
                                                   const su2double distToWall,
                                                         su2double &dMuTdx,
                                                         su2double &dMuTdy,
                                                         su2double &dMuTdz) {
  dMuTdx = dMuTdy = dMuTdz = 0.0;
}

inline CSmagorinskyModel::CSmagorinskyModel(void) : CSGSModel() {
  const_smag  = 0.1;
  filter_mult = 2.0;
}

inline CSmagorinskyModel::~CSmagorinskyModel(void){}

inline su2double CSmagorinskyModel::ComputeEddyViscosity_2D(const su2double rho,
                                                            const su2double dudx,
                                                            const su2double dudy,
                                                            const su2double dvdx,
                                                            const su2double dvdy,
                                                            const su2double lenScale,
                                                            const su2double distToWall) {
  /* Constant coefficient Smagorinsky SGS is calculated:
   * ( C_s * L_c )^2 * |S(x,t)|
   * C_s = Smagorinsky constant
   * L_c = Filter width
   * S(x,t) = Rate of Strain Tensor ( 1/2 [ du_i/dx_j + du_j/dx_i] )
   */
  const su2double C_s_filter_width = const_smag*filter_mult*lenScale;

  const su2double S12          = 0.5*(dudy + dvdx);
  const su2double strain_rate2 = 2.0*(dudx*dudx + dvdy*dvdy + 2.0*S12*S12);

  /* Return the SGS dynamic viscosity. */
  return rho*C_s_filter_width*C_s_filter_width*sqrt(strain_rate2);
}

inline su2double CSmagorinskyModel::ComputeEddyViscosity_3D(const su2double rho,
                                                            const su2double dudx,
                                                            const su2double dudy,
                                                            const su2double dudz,
                                                            const su2double dvdx,
                                                            const su2double dvdy,
                                                            const su2double dvdz,
                                                            const su2double dwdx,
                                                            const su2double dwdy,
                                                            const su2double dwdz,
                                                            const su2double lenScale,
                                                            const su2double distToWall) {
  /* Constant coefficient Smagorinsky SGS is calculated:
   * ( C_s * L_c )^2 * |S(x,t)|
   * C_s = Smagorinsky constant
   * L_c = Filter width
   * S(x,t) = Rate of Strain Tensor ( 1/2 [ du_i/dx_j + du_j/dx_i] )
   */
  const su2double C_s_filter_width = const_smag*filter_mult*lenScale;

  const su2double S12 = 0.5*(dudy + dvdx);
  const su2double S13 = 0.5*(dudz + dwdx);
  const su2double S23 = 0.5*(dvdz + dwdy);

  const su2double strain_rate2 = 2.0*(dudx*dudx + dvdy*dvdy + dwdz*dwdz
                               +      2.0*(S12*S12 + S13*S13 + S23*S23));

  /* Return the SGS dynamic viscosity. */
  return rho*C_s_filter_width*C_s_filter_width*sqrt(strain_rate2);
}

inline void CSmagorinskyModel::ComputeGradEddyViscosity_2D(const su2double rho,
                                                           const su2double drhodx,
                                                           const su2double drhody,
                                                           const su2double dudx,
                                                           const su2double dudy,
                                                           const su2double dvdx,
                                                           const su2double dvdy,
                                                           const su2double d2udx2,
                                                           const su2double d2udy2,
                                                           const su2double d2udxdy,
                                                           const su2double d2vdx2,
                                                           const su2double d2vdy2,
                                                           const su2double d2vdxdy,
                                                           const su2double lenScale,
                                                           const su2double distToWall,
                                                                 su2double &dMuTdx,
                                                                 su2double &dMuTdy) {
  cout << "CSmagorinskyModel::ComputeGradEddyViscosity_2D: Not implemented yet" << endl;
  exit(1);
}

inline void CSmagorinskyModel::ComputeGradEddyViscosity_3D(const su2double rho,
                                                           const su2double drhodx,
                                                           const su2double drhody,
                                                           const su2double drhodz,
                                                           const su2double dudx,
                                                           const su2double dudy,
                                                           const su2double dudz,
                                                           const su2double dvdx,
                                                           const su2double dvdy,
                                                           const su2double dvdz,
                                                           const su2double dwdx,
                                                           const su2double dwdy,
                                                           const su2double dwdz,
                                                           const su2double d2udx2,
                                                           const su2double d2udy2,
                                                           const su2double d2udz2,
                                                           const su2double d2udxdy,
                                                           const su2double d2udxdz,
                                                           const su2double d2udydz,
                                                           const su2double d2vdx2,
                                                           const su2double d2vdy2,
                                                           const su2double d2vdz2,
                                                           const su2double d2vdxdy,
                                                           const su2double d2vdxdz,
                                                           const su2double d2vdydz,
                                                           const su2double d2wdx2,
                                                           const su2double d2wdy2,
                                                           const su2double d2wdz2,
                                                           const su2double d2wdxdy,
                                                           const su2double d2wdxdz,
                                                           const su2double d2wdydz,
                                                           const su2double lenScale,
                                                           const su2double distToWall,
                                                                 su2double &dMuTdx,
                                                                 su2double &dMuTdy,
                                                                 su2double &dMuTdz) {
  cout << "CSmagorinskyModel::ComputeGradEddyViscosity_3D: Not implemented yet" << endl;
  exit(1);
}

inline CWALEModel::CWALEModel(void) : CSGSModel() {
  const_smag = 0.1;
  const_mult = 2.5;
}

inline CWALEModel::~CWALEModel(void){}

inline su2double CWALEModel::ComputeEddyViscosity_2D(const su2double rho,
                                                     const su2double dudx,
                                                     const su2double dudy,
                                                     const su2double dvdx,
                                                     const su2double dvdy,
                                                     const su2double lenScale,
                                                     const su2double distToWall) {
  /*!
   * The WALE sgs viscosity is calculated:
   * C_wale * sqrt( B_beta / (alpha_ij * alpha_ij) )
   * alpha_ij = du_j / dx_i
   * B_beta = beta_11*beta_22 - beta_12^2 + beta_11*beta_33 - beta_13^2 + beta_22*beta_33 - beta_23^2
   * beta_ij = delta_m^2 * alpha_mi * alpha_mi
   * C_wale = 2.5 * C_smagorinsky^2;
   * delta = min(lenScale, distToWall)
   */

  /* Determine the length scale to use. This is the lesser of the characteristic
   * element length scale or the distance to the wall
   */
  su2double delta = std::min(lenScale,distToWall);
  su2double delta_2 = delta*delta;

  /* Calculate two-dimensional beta values */
  su2double beta_11 = delta_2 * ( dudx*dudx + dudy*dudy );
  su2double beta_22 = delta_2 * ( dvdx*dvdx + dvdy*dvdy );
  su2double beta_12 = delta_2 * ( dudx*dvdx + dudy*dvdy );

  /* Calculate B_beta */
  su2double B_beta = beta_11*beta_22 - beta_12*beta_12;

  /* Calculate denominator (sum of squares of velocity gradients) */
  su2double sum_veloc_gradients_2 = ( dudx*dudx + dudy*dudy + dvdx*dvdx + dvdy*dvdy );

  /* Calculate C_wale */
  su2double C_wale = const_mult * const_smag * const_smag;

  /* Return the WALE SGS dynamic viscosity */

  return rho * C_wale * sqrt(B_beta / sum_veloc_gradients_2);
}

inline su2double CWALEModel::ComputeEddyViscosity_3D(const su2double rho,
                                                     const su2double dudx,
                                                     const su2double dudy,
                                                     const su2double dudz,
                                                     const su2double dvdx,
                                                     const su2double dvdy,
                                                     const su2double dvdz,
                                                     const su2double dwdx,
                                                     const su2double dwdy,
                                                     const su2double dwdz,
                                                     const su2double lenScale,
                                                     const su2double distToWall) {
  /*!
   * The WALE sgs viscosity is calculated:
   * C_wale * sqrt( B_beta / (alpha_ij * alpha_ij) )
   * alpha_ij = du_j / dx_i
   * B_beta = beta_11*beta_22 - beta_12^2 + beta_11*beta_33 - beta_13^2 + beta_22*beta_33 - beta_23^2
   * beta_ij = delta_m^2 * alpha_mi * alpha_mi
   * C_wale = 2.5 * C_smagorinsky^2;
   * delta = min(lenScale, distToWall)
   */

  /* Determine the length scale to use. This is the lesser of the characteristic
   * element length scale or the distance to the wall
   */
  su2double delta = std::min(lenScale,distToWall);
  su2double delta_2 = delta*delta;

  /* Calculate three-dimensional beta values */
  su2double beta_11 = delta_2 * ( dudx*dudx + dudy*dudy + dudz*dudz );
  su2double beta_22 = delta_2 * ( dvdx*dvdx + dvdy*dvdy + dvdz*dvdz );
  su2double beta_33 = delta_2 * ( dwdx*dwdx + dwdy*dwdy + dwdz*dwdz );

  su2double beta_12 = delta_2 * ( dudx*dvdx + dudy*dvdy + dudz*dvdz );
  su2double beta_13 = delta_2 * ( dudx*dwdx + dudy*dwdy + dudz*dwdz );
  su2double beta_23 = delta_2 * ( dvdx*dwdx + dvdy*dwdy + dvdz*dwdz );

  /* Calculate B_beta */
  su2double B_beta = beta_11*beta_22 - beta_12*beta_12 + beta_11*beta_33 - beta_13*beta_13 + beta_22*beta_33 - beta_23*beta_23;

  /* Calculate denominator (sum of squares of velocity gradients) */
  su2double sum_veloc_gradients_2 = ( dudx*dudx + dudy*dudy + dudz*dudz + dvdx*dvdx + dvdy*dvdy + dvdz*dvdz + dwdx*dwdx + dwdy*dwdy + dwdz*dwdz );

  /* Calculate C_wale */
  su2double C_wale = const_mult * const_smag * const_smag;

  /* Return the WALE SGS dynamic viscosity */

  return rho * C_wale * sqrt(B_beta / sum_veloc_gradients_2);
}

inline void CWALEModel::ComputeGradEddyViscosity_2D(const su2double rho,
                                                    const su2double drhodx,
                                                    const su2double drhody,
                                                    const su2double dudx,
                                                    const su2double dudy,
                                                    const su2double dvdx,
                                                    const su2double dvdy,
                                                    const su2double d2udx2,
                                                    const su2double d2udy2,
                                                    const su2double d2udxdy,
                                                    const su2double d2vdx2,
                                                    const su2double d2vdy2,
                                                    const su2double d2vdxdy,
                                                    const su2double lenScale,
                                                    const su2double distToWall,
                                                          su2double &dMuTdx,
                                                          su2double &dMuTdy) {
  cout << "CWALEModel::ComputeGradEddyViscosity_2D: Not implemented yet" << endl;
  exit(1);
}

inline void CWALEModel::ComputeGradEddyViscosity_3D(const su2double rho,
                                                    const su2double drhodx,
                                                    const su2double drhody,
                                                    const su2double drhodz,
                                                    const su2double dudx,
                                                    const su2double dudy,
                                                    const su2double dudz,
                                                    const su2double dvdx,
                                                    const su2double dvdy,
                                                    const su2double dvdz,
                                                    const su2double dwdx,
                                                    const su2double dwdy,
                                                    const su2double dwdz,
                                                    const su2double d2udx2,
                                                    const su2double d2udy2,
                                                    const su2double d2udz2,
                                                    const su2double d2udxdy,
                                                    const su2double d2udxdz,
                                                    const su2double d2udydz,
                                                    const su2double d2vdx2,
                                                    const su2double d2vdy2,
                                                    const su2double d2vdz2,
                                                    const su2double d2vdxdy,
                                                    const su2double d2vdxdz,
                                                    const su2double d2vdydz,
                                                    const su2double d2wdx2,
                                                    const su2double d2wdy2,
                                                    const su2double d2wdz2,
                                                    const su2double d2wdxdy,
                                                    const su2double d2wdxdz,
                                                    const su2double d2wdydz,
                                                    const su2double lenScale,
                                                    const su2double distToWall,
                                                          su2double &dMuTdx,
                                                          su2double &dMuTdy,
                                                          su2double &dMuTdz) {
  cout << "CWALEModel::ComputeGradEddyViscosity_3D: Not implemented yet" << endl;
  exit(1);
}
