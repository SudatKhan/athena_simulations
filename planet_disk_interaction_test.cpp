#include <algorithm>  // min
#include <cmath>      // sqrt
#include <cstdlib>    // srand
#include <cstring>    // strcmp()
#include <fstream>
#include <iostream>   // endl
#include <limits>
#include <sstream>    // stringstream
#include <stdexcept>  // runtime_error
#include <string>     // c_str()

// Athena++ headers
#include "../athena.hpp"
#include "../athena_arrays.hpp"
#include "../bvals/bvals.hpp"
#include "../coordinates/coordinates.hpp"
#include "../eos/eos.hpp"
#include "../field/field.hpp"
#include "../globals.hpp"
#include "../hydro/hydro.hpp"
#include "../mesh/mesh.hpp"
#include "../orbital_advection/orbital_advection.hpp"
#include "../parameter_input.hpp"

using namespace std;

namespace {
void GetCylCoord(Coordinates *pco,Real &rad,Real &phi,Real &z,int i,int j,int k);
Real PoverR(const Real rad, const Real phi, const Real z);
// problem parameters which are useful to make global to this file
Real gm_star, rstar, r0, rho0, dslope, r_in, r_out, p0_over_r0, pslope, gamma_gas, gm_planet, gm_planet2, alpha, nu_iso, scale, z, phi, r, rp, rp2, d, dfloor, Omega0, cosine_term, sine_term, epsilon;
} // namespace

// User-defined boundary conditions for disk simulations
void OutflowInner(MeshBlock *pmb, Coordinates *pco,
                  AthenaArray<Real> &prim, FaceField &b,
                  Real time, Real dt,
                  int il, int iu, int jl, int ju, int kl, int ku, int ngh);
void Steady_State_Outer(MeshBlock *pmb, Coordinates *pco,
                  AthenaArray<Real> &prim, FaceField &b,
                  Real time, Real dt,
                  int il, int iu, int jl, int ju, int kl, int ku, int ngh);

//========================================================================================
//! \fn void Mesh::InitUserMeshData(ParameterInput *pin)
//! \brief Function to initialize problem-specific data in mesh class.  Can also be used
//! to initialize variables which are global to (and therefore can be passed to) other
//! functions in this file.  Called in Mesh constructor.
//========================================================================================

void Mesh::InitUserMeshData(ParameterInput *pin) {
  Real x1, x2, x3;
  // Get parameters for gravitatonal potential of central star mass
  r0 = pin->GetOrAddReal("problem","r0",1.0);

  // Get parameters for initial density and velocity
  rho0 = pin->GetReal("problem","rho0");
  dslope = pin->GetOrAddReal("problem","dslope",0.0);

  // Parameters for the radial bounds of the protoplanetary disk
  r_in = pin->GetReal("mesh", "x1min");
  r_out = pin->GetReal("mesh", "x1max");

  // Get parameters for gravitational potential of orbiting protoplanet
  gm_star = pin->GetOrAddReal("problem","starmass",0.0);
  gm_planet = pin -> GetOrAddReal("problem", "planetgm", 0.0);
  gm_planet2 = pin -> GetOrAddReal("problem", "planetgm2", 0.0);
  rp = pin -> GetOrAddReal("problem", "ptosr", 1.0);
  rp2 = pin -> GetOrAddReal("problem", "ptosr2", 0.0);

  // Get viscosity parameters and scale ratio
  alpha = pin -> GetOrAddReal("problem", "alpha", 0.0);
  nu_iso = pin -> GetOrAddReal("problem", "nu_iso", 0.0);
  scale = pin -> GetOrAddReal("hydro", "iso_sound_speed", 0.05);

  // Get parameters of initial pressure and cooling parameters
  if (NON_BAROTROPIC_EOS) {
    p0_over_r0 = pin->GetOrAddReal("problem","p0_over_r0",0.0025);
    pslope = pin->GetOrAddReal("problem","pslope",0.0);
    gamma_gas = pin->GetReal("hydro","gamma");
  } else {
    p0_over_r0=SQR(pin->GetReal("hydro","iso_sound_speed"));
  }
  Real float_min = std::numeric_limits<float>::min();
  dfloor=pin->GetOrAddReal("hydro","dfloor",(1024*(float_min)));

  Omega0 = pin->GetOrAddReal("orbital_advection","Omega0",0.0);

  // enroll user-defined boundary condition
  if (mesh_bcs[BoundaryFace::inner_x1] == GetBoundaryFlag("user")) {
    EnrollUserBoundaryFunction(BoundaryFace::inner_x1, OutflowInner);
  }
  if (mesh_bcs[BoundaryFace::outer_x1] == GetBoundaryFlag("user")) {
    EnrollUserBoundaryFunction(BoundaryFace::outer_x1, Steady_State_Outer);
  }

  void StarandPlanet(MeshBlock *pmb, const Real time, const Real dt, const AthenaArray<Real> &prim,
              const AthenaArray<Real> &prim_scalar, const AthenaArray<Real> &bcc,
              AthenaArray<Real> &cons, AthenaArray<Real> &cons_scalar);
  EnrollUserExplicitSourceFunction(StarandPlanet);

  void Viscosity(HydroDiffusion *phdif, MeshBlock *pmb, const AthenaArray<Real> &prim, const AthenaArray<Real> &bcc, 
            int is, int ie, int js, int je, int ks, int ke);
  EnrollViscosityCoefficient(Viscosity);

  /*Real Torque(MeshBlock *pmb, int iout);
  Real Torque2(MeshBlock *pmb, int iout);
  AllocateUserHistoryOutput(2);
  EnrollUserHistoryOutput(0, Torque, "first planet torque");
  EnrollUserHistoryOutput(1, Torque2, "second planet torque");*/

  Real Inner_Lindblad_Torque1(MeshBlock *pmb, int iout);
  Real Outer_Lindblad_Torque1(MeshBlock *pmb, int iout);
  AllocateUserHistoryOutput(2);
  EnrollUserHistoryOutput(0, Inner_Lindblad_Torque1, "First Planet Inner Lindblad Torque");
  EnrollUserHistoryOutput(1, Outer_Lindblad_Torque1, "First Planet Outer Lindblad Torque");
  return;
}
/*void MeshBlock::InitUserMeshBlockData(ParameterInput *pin) {
    AllocateUserOutputVariables(2);
    return;
}*/
  
//========================================================================================
//! \fn void MeshBlock::ProblemGenerator(ParameterInput *pin)
//! \brief Initializes Keplerian accretion disk.
//========================================================================================

void MeshBlock::ProblemGenerator(ParameterInput *pin) {
  Real rad(0.0), phi(0.0), z(0.0);
  Real den, vel;
  Real x1,x2,x3;

  OrbitalVelocityFunc &vK = porb->OrbitalVelocity;
  for (int k=ks; k<=ke; ++k) {
    z = pcoord->x3v(k);
    for (int j=js; j<=je; ++j) {
      phi = pcoord->x2v(j);
      for (int i=is; i<=ie; ++i) {
        r = pcoord->x1v(i);
        GetCylCoord(pcoord,rad,phi,z,i,j,k); // convert to cylindrical coordinates
        Real surface_density = rho0 / sqrt(r);
        Real v_r = -3.0/2.0 * alpha * pow(scale,2) * sqrt((gm_star+gm_planet)/r);
        Real v_phi = r * sqrt(1-0.5*pow(scale,2)) * sqrt(gm_star+gm_planet)* sqrt(1 / pow(r,3));
        phydro->u(IDN,k,j,i) = surface_density;
        phydro->u(IM1,k,j,i) = surface_density * v_r;
        if (std::strcmp(COORDINATE_SYSTEM, "cylindrical") == 0) {
          phydro->u(IM2,k,j,i) = surface_density * v_phi;
          phydro->u(IM3,k,j,i) = 0.0;
        }
        if (NON_BAROTROPIC_EOS) {
          Real p_over_r = PoverR(rad,phi,z);
          phydro->u(IEN,k,j,i) = p_over_r*phydro->u(IDN,k,j,i)/(gamma_gas - 1.0);
          phydro->u(IEN,k,j,i) += 0.5*(SQR(phydro->u(IM1,k,j,i))+SQR(phydro->u(IM2,k,j,i))
                                       + SQR(phydro->u(IM3,k,j,i)))/phydro->u(IDN,k,j,i);
        }
      }
    }
  }
  return;
}

void AddGOneObject(Real x, Real y, Real xp, Real yp, Real & gx, Real & gy, Real d_smooth, Real gm) {
  // add one object's acceleration to gx, gy
  Real dx = x-xp;
  Real dy = y-yp;
  Real d = sqrt(dx*dx+dy*dy);
  // g = -gm * d / |d| / (|d^2| + d_smooth^2)
  gx += -gm * dx / d / (d*d + d_smooth*d_smooth);
  gy += -gm * dy / d / (d*d + d_smooth*d_smooth);
  return;
}

void StarandPlanet(MeshBlock *pmb, const Real time, const Real dt, const AthenaArray<Real> &prim,
            const AthenaArray<Real> &prim_scalar, const AthenaArray<Real> &bbc,
            AthenaArray<Real> &cons, AthenaArray<Real> &cons_scalar) {
  // star and planet positions (xy)
  Real R_star[2] = {0.0, 0.0};    // star
  Real R_planet1[2] = {0.0, 0.0}; // planet 1
  Real R_planet2[2] = {0.0, 0.0}; // planet 2
  // compute the positions
  Real period = 2.*M_PI*sqrt(pow(rp,3)/(gm_star + gm_planet));
  Real phip1 = 2.*(M_PI / period)*time;
  period = 2.*M_PI*sqrt(pow(rp2,3)/(gm_star + gm_planet + gm_planet2));
  Real phip2 = 2.*(M_PI / period)*time;
  R_planet1[0] = rp*cos(phip1);
  R_planet1[1] = rp*sin(phip1);
  R_planet2[0] = rp2*cos(phip2);
  R_planet2[1] = rp2*sin(phip2);
  R_star[0] = - R_planet1[0]*gm_planet/gm_star - R_planet2[0]*gm_planet2/gm_star;
  R_star[1] = - R_planet1[1]*gm_planet/gm_star - R_planet2[1]*gm_planet2/gm_star;
  // apply forces + enforce isothermal
  for (int k = pmb->ks; k <= pmb->ke; ++k) {
    z = pmb->pcoord->x3v(k);
    for (int j = pmb->js; j <= pmb->je; ++j) {
      phi = pmb->pcoord->x2v(j);
      for (int i = pmb->is; i <= pmb->ie; ++i) {
        r = pmb->pcoord->x1v(i);
        //density initialization
        Real dens = prim(IDN,k,j,i);
        //compute acceleration
        Real x, y, gx, gy, gr, gphi;
        x = r*cos(phi);
        y = r*sin(phi);
        gx=0.; gy=0.;
        AddGOneObject(x,y,R_star[0],R_star[1],gx,gy,0.,gm_star);
        Real R_H = rp*cbrt(gm_planet/(3*gm_star));
        AddGOneObject(x,y,R_planet1[0],R_planet1[1],gx,gy,0.3*R_H,gm_planet);
        R_H = rp2*cbrt(gm_planet2/(3*gm_star));
        AddGOneObject(x,y,R_planet2[0],R_planet2[1],gx,gy,0.3*R_H,gm_planet2);
        // convert gx gy to gr gphi
        gr = (x*gx+y*gy)/r;
        gphi = (x*gy-y*gx)/r;
        // apply force & energy source term
        cons(IM1, k,j,i) += prim(IDN,k,j,i)*gr*dt;
        cons(IM2, k,j,i) += prim(IDN,k,j,i)*gphi*dt;
        if (NON_BAROTROPIC_EOS) cons(IEN,k,j,i) += (prim(IDN,k,j,i)*prim(IVX,k,j,i)*gr + prim(IDN,k,j,i)*prim(IVY,k,j,i)*gphi) * dt;
        // update temperature
        //Real gamma = (rho0*p0_over_r0) / (pow(r0, dslope));
        //Real beta = rho0/(pow(r0, dslope));
        //Real pressure_0 = gamma * pow(r,pslope+dslope);
        //Real surface_density_0 = beta * pow(r, dslope);
        //Real pressure = dens * (pressure_0/surface_density_0); //definition of isothermal eos
        //if (NON_BAROTROPIC_EOS) cons(IEN,k,j,i) += 3.0/2.0 * (pressure-prim(IPR,k,j,i));    
        if (NON_BAROTROPIC_EOS) {
          Real p_over_r = PoverR(r,phi,z); //temperature profile which scales with radius
          cons(IEN,k,j,i) = p_over_r*cons(IDN,k,j,i)/(gamma_gas - 1.0);
          cons(IEN,k,j,i) += 0.5*(SQR(cons(IM1,k,j,i))+SQR(cons(IM2,k,j,i))
                                       + SQR(cons(IM3,k,j,i)))/cons(IDN,k,j,i);
        }
      }
    }
  }
  return;
}

void Viscosity(HydroDiffusion *phdif, MeshBlock *pmb, const AthenaArray<Real> &prim, const AthenaArray<Real> &bcc, 
               int is, int ie, int js, int je, int ks, int ke) {
    if (phdif->nu_iso > 0.0) {
      for (int k = ks; k <= ke; ++k) {
        z = pmb->pcoord->x3v(k);
        for (int j = js; j <= je; ++j) {
          phi = pmb->pcoord->x2v(j);
          for (int i = is; i <= ie; ++i) {
            r = pmb->pcoord->x1v(i);
            Real omega = sqrt((gm_star + gm_planet)/(pow(r,3)));
            Real sound_speed = scale * omega*r;
            Real kinematic_viscosity = alpha * sound_speed * (sound_speed/omega); 
            phdif->nu(HydroDiffusion::DiffProcess::iso,k,j,i) = kinematic_viscosity;
        }
      }
    }
  }
}

/*void MeshBlock::UserWorkBeforeOutput(ParameterInput *pin) {
  Real time1 = pmy_mesh -> time;
  for (int k = ks; k <= ke; ++k) {
    z = pcoord->x3v(k);
    for (int j = js; j <= je; ++j) {
      phi = pcoord->x2v(j);
      for (int i = is; i <= ie; ++i) {
        r = pcoord->x1v(i);
        Real period = 2*M_PI*sqrt(pow(rp,3)/gm_star);
        Real phip = 2*(M_PI / period)*time1;
        d = sqrt(pow(rp,2) + pow(r,2) - 2*rp*r*cos(phi - phip));
        epsilon = 0.3;
        Real R_H = rp*cbrt(gm_planet/(3*gm_star));
        Real g_mag = -1*((gm_planet*d) / (sqrt(pow(pow(d,2) + pow(epsilon,2)*pow(R_H,2), 3))));
        cosine_term = (pow(r,2)*(pow(cos(phi),2)) - r*rp*cos(phi)*cos(phip) + pow(r,2)*(pow(sin(phi),2)) - r*rp*sin(phi)*sin(phip)) / (r*d);
        sine_term = (r*rp*cos(phi)*sin(phip) - r*rp*sin(phi)*cos(phip)) / (r*d);
        user_out_var(0,k,j,i) = g_mag*cosine_term;
        user_out_var(1,k,j,i) = -g_mag*sine_term;
      }
    }
  }
}*/

/*Real Torque(MeshBlock *pmb, int iout) { //planet one torque
  int is=pmb->is, ie=pmb->ie, js=pmb->js, je=pmb->je, ks=pmb->ks, ke=pmb->ke;
  Real sum_torque = 0;
  Real time2 = pmb->pmy_mesh->time;
  for(int k=ks; k<=ke; k++) {
    z = pmb->pcoord->x3v(k);
    for(int j=js; j<=je; j++) {
      phi = pmb->pcoord->x2v(j);
      for(int i=is; i<=ie; i++) {
        r = pmb->pcoord->x1v(i);
        Real period = 2 * M_PI * sqrt(pow(rp, 3) / (gm_star + gm_planet));
        Real phip = 2 * (M_PI / period) * time2;
        Real d = sqrt(pow(rp,2) + pow(r,2) - 2*rp*r*cos(phi - phip));
        Real R_H = rp*cbrt(gm_planet/(3*gm_star));
        Real g_mag = -1*((gm_planet*d) / (sqrt(pow(pow(d,2) + pow(epsilon,2)*pow(R_H,2), 3))));
        Real dens = pmb->phydro->u(IDN,k,j,i);
        Real volume = pmb ->pcoord->GetCellVolume(k,j,i);
        Real sine_term = (r*rp*cos(phi)*sin(phip) - r*rp*sin(phi)*cos(phip)) / (r*d);
        sum_torque +=  dens * volume *r * g_mag * sine_term;
      }
    }
  }
  return sum_torque;
}

Real Torque2 (MeshBlock *pmb, int iout) { //planet two torque
  int is=pmb->is, ie=pmb->ie, js=pmb->js, je=pmb->je, ks=pmb->ks, ke=pmb->ke;
  Real sum_torque2 = 0;
  Real time3 = pmb->pmy_mesh->time;
  for(int k=ks; k<=ke; k++) {
    z = pmb->pcoord->x3v(k); 
    for(int j=js; j<=je; j++) {
      phi = pmb->pcoord->x2v(j);
      for(int i=is; i<=ie; i++) {
        r = pmb->pcoord->x1v(i);
        Real period = 2 * M_PI * sqrt(pow(rp2, 3) / (gm_star +gm_planet + gm_planet2));
        Real phip = 2 * (M_PI / period) * time3;
        Real d = sqrt(pow(rp2,2) + pow(r,2) - 2*rp2*r*cos(phi - phip));
        Real R_H = rp2*cbrt(gm_planet2/(3*gm_star));
        Real g_mag = -1*((gm_planet2*d) / (sqrt(pow(pow(d,2) + pow(epsilon,2)*pow(R_H,2), 3))));
        Real dens = pmb->phydro->u(IDN,k,j,i);
        Real volume = pmb ->pcoord->GetCellVolume(k,j,i);
        Real sine_term = (r*rp2*cos(phi)*sin(phip) - r*rp2*sin(phi)*cos(phip)) / (r*d);
        sum_torque2 +=  dens * volume *r * g_mag * sine_term;       
      }
    }
  }
  return sum_torque2;
} */

/*Real Lindblad_Torque1 (MeshBlock *pmb, int iout) {
  int is=pmb->is, ie=pmb->ie, js=pmb->js, je=pmb->je, ks=pmb->ks, ke=pmb->ke;
  Real sum_lindblad_torque_inner = 0;
  Real sum_lindblad_torque_outer = 0;
  Real time4= pmb->pmy_mesh->time;
  for (int k=ks; k<ke; k++) {
    z = pmb->pcoord->x3v(k);
    for (int j=js; j<je; j++) {
      phi = pmb->pcoord->x2v(j);
      for (int i=is; i<ie; i++) {
        r = pmb->pcoord->x1v(i);        
        Real mass_ratio = (gm_planet/(pow(scale,3.0)));
        Real horseshoe = scale * (rp) * ((1.05 * pow(mass_ratio, 0.5) + 3.4 * pow(mass_ratio, 7.0/3.0)) / (1.0 + 2.0*pow(mass_ratio, 2.0)));
        Real inner_horseshoe = rp - horseshoe;
        Real outer_horseshoe = rp + horseshoe;
          Real period = 2 * M_PI * sqrt(pow(rp, 3) / (gm_star + gm_planet));
          Real phip = 2 * (M_PI / period) * time4;
          Real d = sqrt(pow(rp,2) + pow(r,2) - 2*rp*r*cos(phi - phip));
          Real R_H = rp*cbrt(gm_planet/(3*gm_star));
          Real g_mag = -1*((gm_planet*d) / (sqrt(pow(pow(d,2) + pow(epsilon,2)*pow(R_H,2), 3))));
          Real dens = pmb->phydro->u(IDN,k,j,i);
          Real area = pmb ->pcoord->GetCellVolume(k,j,i);
          Real sine_term = (r*rp*cos(phi)*sin(phip) - r*rp*sin(phi)*cos(phip)) / (r*d);
          sum_lindblad_torque_inner += dens * r * g_mag * sine_term * area;

        if (r >= outer_horseshoe && r <= r_out) {
          Real period = 2 * M_PI * sqrt(pow(rp, 3) / (gm_star + gm_planet));
          Real phip = 2 * (M_PI / period) * time4;
          Real d = sqrt(pow(rp,2) + pow(r,2) - 2*rp*r*cos(phi - phip));
          Real R_H = rp*cbrt(gm_planet/(3*gm_star));
          Real g_mag = -1*((gm_planet*d) / (sqrt(pow(pow(d,2) + pow(epsilon,2)*pow(R_H,2), 3))));
          Real dens = pmb->phydro->u(IDN,k,j,i);
          Real area = pmb ->pcoord->GetCellVolume(k,j,i);
          Real sine_term = (r*rp*cos(phi)*sin(phip) - r*rp*sin(phi)*cos(phip)) / (r*d);
          sum_lindblad_torque_outer += dens * r * g_mag * sine_term * area;
        }
      }
    }
  }
  return sum_lindblad_torque_inner;
  // return sum_lindblad_torque_outer;
}*/

Real Inner_Lindblad_Torque1 (MeshBlock *pmb, int iout) {
  int is=pmb->is, ie=pmb->ie, js=pmb->js, je=pmb->je, ks=pmb->ks, ke=pmb->ke;
  Real sum_lindblad_torque_inner = 0;
  Real time4= pmb->pmy_mesh->time;
  for(int k=ks; k<=ke; k++) {
    z = pmb->pcoord->x3v(k);
    for(int j=js; j<=je; j++) {
      phi = pmb->pcoord->x2v(j);
      for(int i=is; i<=ie; i++) {
        r = pmb->pcoord->x1v(i);
        Real mass_ratio = (gm_planet/(pow(scale,3.0)));
        Real horseshoe = scale * (rp) * ((1.05 * pow(mass_ratio, 0.5) + 3.4 * pow(mass_ratio, 7.0/3.0)) / (1.0 + 2.0*pow(mass_ratio, 2.0)));
        Real inner_horseshoe = rp - horseshoe;
        Real outer_horseshoe = rp + horseshoe;
        if (r <= inner_horseshoe) {
          Real period = 2 * M_PI * sqrt(pow(rp, 3) / (gm_star + gm_planet));
          Real phip = 2 * (M_PI / period) * time4;
          Real d = sqrt(pow(rp,2) + pow(r,2) - 2*rp*r*cos(phi - phip));
          Real R_H = rp*cbrt(gm_planet/(3*gm_star));
          Real g_mag = -1*((gm_planet*d) / (sqrt(pow(pow(d,2) + pow(epsilon,2)*pow(R_H,2), 3))));
          Real dens = pmb->phydro->u(IDN,k,j,i);
          Real area = pmb ->pcoord->GetCellVolume(k,j,i);
          Real sine_term = (r*rp*cos(phi)*sin(phip) - r*rp*sin(phi)*cos(phip)) / (r*d);
          sum_lindblad_torque_inner +=  dens * r * g_mag * sine_term * area;
        }
      }
    }
  }
  return sum_lindblad_torque_inner;
}

Real Outer_Lindblad_Torque1 (MeshBlock *pmb, int iout) {
  int is=pmb->is, ie=pmb->ie, js=pmb->js, je=pmb->je, ks=pmb->ks, ke=pmb->ke;
  Real sum_lindblad_torque_outer = 0;
  Real time5 = pmb->pmy_mesh->time;
  for(int k=ks; k<=ke; k++) {
    z = pmb->pcoord->x3v(k);
    for(int j=js; j<=je; j++) {
      phi = pmb->pcoord->x2v(j);
      for(int i=is; i<=ie; i++) {
        r = pmb->pcoord->x1v(i);
        Real mass_ratio = (gm_planet/(pow(scale,3.0)));
        Real horseshoe = scale * (rp) * ((1.05 * pow(mass_ratio, 0.5) + 3.4 * pow(mass_ratio, 7.0/3.0)) / (1.0 + 2.0*pow(mass_ratio, 2.0)));
        Real inner_horseshoe = rp - horseshoe;
        Real outer_horseshoe = rp + horseshoe;
        if (r >= outer_horseshoe) {
          Real period = 2 * M_PI * sqrt(pow(rp, 3) / (gm_star + gm_planet));
          Real phip = 2 * (M_PI / period) * time5;
          Real d = sqrt(pow(rp,2) + pow(r,2) - 2*rp*r*cos(phi - phip));
          Real R_H = rp*cbrt(gm_planet/(3*gm_star));
          Real g_mag = -1*((gm_planet*d) / (sqrt(pow(pow(d,2) + pow(epsilon,2)*pow(R_H,2), 3))));
          Real dens = pmb->phydro->u(IDN,k,j,i);
          Real area = pmb ->pcoord->GetCellVolume(k,j,i);
          Real sine_term = (r*rp*cos(phi)*sin(phip) - r*rp*sin(phi)*cos(phip)) / (r*d);
          sum_lindblad_torque_outer +=  dens * r * g_mag * sine_term * area;
        }
      }
    }
  }
  return sum_lindblad_torque_outer;
}

namespace {
//----------------------------------------------------------------------------------------
//! transform to cylindrical coordinate

void GetCylCoord(Coordinates *pco,Real &rad,Real &phi,Real &z,int i,int j,int k) {
  if (std::strcmp(COORDINATE_SYSTEM, "cylindrical") == 0) {
    rad=pco->x1v(i);
    phi=pco->x2v(j);
    z=pco->x3v(k);
  } else if (std::strcmp(COORDINATE_SYSTEM, "spherical_polar") == 0) {
    rad=std::abs(pco->x1v(i)*std::sin(pco->x2v(j)));
    phi=pco->x3v(k);
    z=pco->x1v(i)*std::cos(pco->x2v(j));
  }
  return;
}

//----------------------------------------------------------------------------------------
//! computes pressure/density in cylindrical coordinates

Real PoverR(const Real rad, const Real phi, const Real z) {
  Real poverr;
  poverr = p0_over_r0*std::pow(rad/r0, pslope);
  return poverr;
}
} // namespace

//----------------------------------------------------------------------------------------
//! User-defined Boundary and Initial Conditions
void Steady_State_Inner(MeshBlock *pmb, Coordinates *pco,
                    AthenaArray<Real> &prim, FaceField &b,
                    Real time, Real dt,
                    int il, int iu, int jl, int ju, int kl, int ku, int ngh) {
  for (int k=kl; k<=ku; ++k) {
    z = pmb->pcoord->x3v(k);
    for (int j=jl; j<=ju; ++j) {
      phi = pmb->pcoord->x2v(j);
      for (int i=1; i<=ngh; ++i) {
        r = pmb->pcoord->x1v(il-i);
        Real gamma = (rho0*p0_over_r0) / (pow(r0, dslope));
        Real beta = rho0/(pow(r0, dslope));
        Real pressure_0 = gamma * pow(r, pslope+dslope);
        Real surface_density_0 = beta * pow(r, dslope);
        Real surface_density = rho0 / sqrt(r);
        Real pressure = surface_density * (pressure_0/surface_density_0);
        Real v_r = -3.0/2.0 * alpha * pow(scale,2) * sqrt((gm_star)/r);
        Real v_phi = r * sqrt(1-0.5*pow(scale,2)) * sqrt(gm_star)* sqrt(1 / pow(r,3));
        prim(IDN,k,j,il-i) = surface_density;
        prim(IPR,k,j,il-i) = pressure;
        prim(IVX,k,j,il-i) = v_r;
        prim(IVY,k,j,il-i) = v_phi;
      }
    }
  }
}

//----------------------------------------------------------------------------------------
//! User-defined Boundary and Initial Conditions
void OutflowInner(MeshBlock *pmb,Coordinates *pco, AthenaArray<Real> &prim, FaceField &b,
                 Real time, Real dt,
                 int il, int iu, int jl, int ju, int kl, int ku, int ngh) {
  if (std::strcmp(COORDINATE_SYSTEM, "cylindrical") == 0) {
    for (int k=kl; k<=ku; ++k) {
      z = pmb->pcoord->x3v(k);
      for (int j=jl; j<=ju; ++j) {
        phi = pmb->pcoord->x2v(j);
        for (int i=1; i<=ngh; ++i) {
          Real r_active = pmb->pcoord->x1v(il);
          Real r_ghost = pmb->pcoord->x1v(il-i);
          Real omega = sqrt((gm_star + gm_planet)/(pow(r_active,3)));
          Real sound_speed = scale * omega*r_active;
          Real kinematic_viscosity = alpha * sound_speed * (sound_speed/omega); 
          prim(IDN,k,j,il-i) = prim(IDN,k,j,il)* 1.0/sqrt(r_ghost/r_active);
          prim(IVX,k,j,il-i) = prim(IVX,k,j,il)* 1.0/sqrt(r_ghost/r_active);
          //if (abs(prim(IVX,k,j,il-i)) > 3.0/2.0 * kinematic_viscosity/r_active)
            //prim(IVX,k,j,i,il-i) = -3.0/2.0 * kinematic_viscosity/r_active;
          prim(IVY,k,j,il-i) = prim(IVY,k,j,il) * 1.0/sqrt(r_ghost/r_active);
          prim(IVZ,k,j,il-i) = prim(IVZ,k,j,il);
          if (NON_BAROTROPIC_EOS) 
            prim(IPR,k,j,il-i) = prim(IPR,k,j,il) * pow((r_ghost/r_active), -3.0/2.0); 
        }
      }
    }
  }
}

void Steady_State_Outer(MeshBlock *pmb, Coordinates *pco,
                    AthenaArray<Real> &prim, FaceField &b,
                    Real time, Real dt,
                    int il, int iu, int jl, int ju, int kl, int ku, int ngh) {
  for (int k=kl; k<=ku; ++k) {
    z = pmb->pcoord->x3v(k);
    for (int j=jl; j<=ju; ++j) {
      phi = pmb->pcoord->x2v(j);
      for (int i=1; i<=ngh; ++i) {
        r = pmb->pcoord->x1v(iu+i);
        Real gamma = (rho0*p0_over_r0) / (pow(r0, dslope));
        Real beta = rho0/(pow(r0, dslope));
        Real pressure_0 = gamma * pow(r, pslope+dslope);
        Real surface_density_0 = beta * pow(r, dslope);
        Real surface_density = rho0 / sqrt(r);
        Real pressure = surface_density * (pressure_0/surface_density_0);
        Real v_r = -3.0/2.0 * alpha * pow(scale,2) * sqrt((gm_star)/r);
        Real v_phi = r * sqrt(1-0.5*pow(scale,2)) * sqrt(gm_star)* sqrt(1 / pow(r,3));
        prim(IDN,k,j,iu+i) = surface_density;
        prim(IPR,k,j,iu+i) = pressure;
        prim(IVX,k,j,iu+i) = v_r;
        prim(IVY,k,j,iu+i) = v_phi;
      }
    }
  }
}