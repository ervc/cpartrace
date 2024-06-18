#include <stdlib.h>
#include <stdio.h>
#include <math.h>

// #include "newmesh.h"
// #include "domain.h"

//TODO: put these into model struct?

// disk values
#define ALPHA    ( 1.0e-3 )   // alpha viscosity
#define ASPECT     0.05       // aspect ratio
#define FLARING    0.25       // flaring angle

// Planet position
#define PLANETMASS ( 0.0 * MEARTH )
#define PLANX    ( R0 )
#define PLANY      0.0
#define PLANZ      0.0

// Sun's position
#define SUNX     ( -PLANX * PLANETMASS/MSUN )
#define SUNY       0.0
#define SUNZ       0.0


typedef struct Model {
    // a struct to hold info for the Many fields in a fargo model
    char fargodir[100];
    int nout;
    size_t nx;
    size_t ny;
    size_t nz;

    double alpha;
    double aspect;
    double flaring;
    double planetmass;

    // data
    MeshField *gasdens;
    MeshField *gasvphi;
    MeshField *gasvr;
    MeshField *gasvtheta;

    // domain
    Domain *domain;

    // derived fields
    MeshField *gasvx;
    MeshField *gasvy;
    MeshField *gasvz;

    MeshField *drhodx;
    MeshField *drhody;
    MeshField *drhodz;

    // MeshField *cartXgrid;
    // MeshField *cartYgrid;
    // MeshField *cartZgrid;
    

} Model;

void make_cartvels(Model *model);
void init_gradrho(Model *model);

Model *init_Model(char* fargodir, int nout, size_t nx, size_t ny, size_t nz) {
    // allocate memory
    Model *model = (Model*)malloc(sizeof(*model));
    if (!model) {
        perror("Malloc Failed on Model Creation");
        exit(1);
    }
    // string copy the fargo directory
    snprintf(model->fargodir,100,"%s",fargodir);
    // store model params
    model->nout = nout;
    model->nx = nx;
    model->ny = ny;
    model->nz = nz;
    // read in the files
    char rhofile[100];
    char vphifile[100];
    char vrfile[100];
    char vthetafile[100];
    int cx;
    cx = snprintf(rhofile,100,"%s/gasdens%d.dat",fargodir,nout);
    cx = snprintf(vphifile,100,"%s/gasvx%d.dat",fargodir,nout);
    cx = snprintf(vrfile,100,"%s/gasvy%d.dat",fargodir,nout);
    cx = snprintf(vthetafile,100,"%s/gasvz%d.dat",fargodir,nout);
    if (cx>100) {
        perror("Fargodir is too long!");
        exit(1);
    }
    // make the MeshFields
    model->gasdens   = init_MeshField_fromFile(rhofile,nx,ny,nz,RHO);
    model->gasvphi   = init_MeshField_fromFile(vphifile,nx,ny,nz,VPHI);
    model->gasvr     = init_MeshField_fromFile(vrfile,nx,ny,nz,VR);
    model->gasvtheta = init_MeshField_fromFile(vthetafile,nx,ny,nz,VTHETA);

    // initialize the domain
    model->domain = init_Domain(fargodir,nx,ny,nz);

    // get the cartesian gas velocities
    printf("Making cartvels...\n");
    make_cartvels(model);
    
    // make the density gradients in cartesian
    printf("Getting Gradrho...\n");
    init_gradrho(model);

    return model;
}

void free_Model(Model *model) {
    free_MeshField(model->gasdens);
    free_MeshField(model->gasvphi);
    free_MeshField(model->gasvr);
    free_MeshField(model->gasvtheta);

    free_Domain(model->domain);

    free_MeshField(model->gasvx);
    free_MeshField(model->gasvy);
    free_MeshField(model->gasvz);
    free_MeshField(model->drhodx);
    free_MeshField(model->drhody);
    free_MeshField(model->drhodz);
    // free_MeshField(model->cartXgrid);
    // free_MeshField(model->cartYgrid);
    // free_MeshField(model->cartZgrid);

    free (model);
}

void make_cartvels(Model *model) {
    // printf("Making cartesian velocity grids\n");
    // pull out some useful values
    size_t nx = model->nx;
    size_t ny = model->ny;
    size_t nz = model->nz;
    Domain *domain = model->domain;

    // init empty mesh fields
    model->gasvx = init_MeshField(nx, ny, nz);
    model->gasvy = init_MeshField(nx, ny, nz);
    model->gasvz = init_MeshField(nx, ny, nz);

    double phi, theta;
    // double x,y,z;
    double vphi, vr, vtheta;
    double vx, vy, vz;
    size_t idx = 0;
    for (size_t k=0; k<nz; k++) {
        theta = domain->thetaCenters[k];
        for (size_t j=0; j<ny; j++) {
            // r = domain->rCenters[j];
            for (size_t i=0; i<nx; i++) {
                phi = domain->phiCenters[i];

                // x = r * cos(phi) * sin(theta);
                // y = r * sin(phi) * sin(theta);
                // z = r * cos(theta);
                // model->cartXgrid->data[idx] = x;
                // model->cartYgrid->data[idx] = y;
                // model->cartZgrid->data[idx] = z;

                vphi   = model->gasvphi->data[idx];
                vr     = model->gasvr->data[idx];
                vtheta = model->gasvtheta->data[idx];

                // NOTE: don't multiply by r in derivs because 
                // vphi, vtheta = r*dot(phi), r*dot(theta)
                vx = vr*cos(phi)*sin(theta) + -vphi*sin(phi)*sin(theta) + vtheta*cos(phi)*cos(theta);
                vy = vr*sin(phi)*sin(theta) +  vphi*cos(phi)*sin(theta) + vtheta*sin(phi)*cos(theta);
                vz = vr*cos(theta) + -vtheta*sin(theta);

                model->gasvx->data[idx] = vx;
                model->gasvy->data[idx] = vy;
                model->gasvz->data[idx] = vz;

                idx++;
            }
        }
    }

    
}

MeshField *get_mesh(Model *model, int which) {
    
    MeshField *data_values;

    switch (which)
    {
    case RHO:
        data_values = model->gasdens;
        break;
    case VPHI:
        data_values = model->gasvphi;
        break;
    case VR:
        data_values = model->gasvr;
        break;
    case VTHETA:
        data_values = model->gasvtheta;
        break;
    case VX:
        data_values = model->gasvx;
        break;
    case VY:
        data_values = model->gasvy;
        break;
    case VZ:
        data_values = model->gasvz;
        break;
    case DRHODX:
        data_values = model->drhodx;
        break;
    case DRHODY:
        data_values = model->drhody;
        break;
    case DRHODZ:
        data_values = model->drhodz;
        break;
    default:
        perror("Invalid choice in get_array!");
        exit(1);
        break;
    }
    return data_values;
}

double *get_array(Model *model, int which) {
    return get_mesh(model, which)->data;
}

// analytic values
double get_Omega(double r) {
    return sqrt(G*MSUN/r/r/r);
}

double get_scaleheight(double r) {
    return r*ASPECT*pow(r/R0,FLARING);
}

double get_soundspeed(Model *model, double r) {
    double H = get_scaleheight(r);
    double OM = get_Omega(r);
    return H*OM;
}

void init_gradrho(Model *model) {
    // init empty mesh fields
    model->drhodx = init_MeshField(model->nx, model->ny, model->nz);
    model->drhody = init_MeshField(model->nx, model->ny, model->nz);
    model->drhodz = init_MeshField(model->nx, model->ny, model->nz);

    // double *rho = get_array(model, RHO);
    MeshField *rho = get_mesh(model, RHO);
    double *PHI = model->domain->phiCenters;
    double *R = model->domain->rCenters;
    double *THETA = model->domain->thetaCenters;

    double drhodphi, drhodr, drhodtheta;
    double theta, r, phi;
    size_t idx = 0;
    for (size_t k=0; k<model->nz; k++) {
        theta = THETA[k];
        for (size_t j=0; j<model->ny; j++) {
            r = R[j];
            for (size_t i=0; i<model->nx; i++) {
                phi = PHI[i];

    if (k==0) {
        drhodtheta = (get_data(rho,k+1,j,i)-get_data(rho,k  ,j,i))/(THETA[k+1] - THETA[k  ]);
    } else if (k==model->nz-1) {
        drhodtheta = (get_data(rho,k  ,j,i)-get_data(rho,k-1,j,i))/(THETA[k  ] - THETA[k-1]);
    } else {
        drhodtheta = (get_data(rho,k+1,j,i)-get_data(rho,k-1,j,i))/(THETA[k+1] - THETA[k-1]);
    }

    if (j==0) {
        drhodr = (get_data(rho,k,j+1,i)-get_data(rho,k,j  ,i))/(R[j+1] - R[j  ]);
    } else if (j==model->ny-1) {
        drhodr = (get_data(rho,k,j  ,i)-get_data(rho,k,j-1,i))/(R[j  ] - R[j-1]);
    } else {
        drhodr = (get_data(rho,k,j+1,i)-get_data(rho,k,j-1,i))/(R[j+1] - R[j-1]);
    }

    if (i==0) {
        drhodphi = (get_data(rho,k,j,i+1)-get_data(rho,k,j,i  ))/(PHI[i+1] - PHI[i  ]);
    } else if (i==model->nx-1) {
        drhodphi = (get_data(rho,k,j,i  )-get_data(rho,k,j,i-1))/(PHI[i  ] - PHI[i-1]);
    } else {
        drhodphi = (get_data(rho,k,j,i+1)-get_data(rho,k,j,i-1))/(PHI[i+1] - PHI[i-1]);
    }

    // chain rules
    double dphidx = -sin(phi)/r/sin(theta);
    double dphidy =  cos(phi)/r/sin(theta);
    double dphidz =  0;
    double drdx = cos(phi)*sin(theta);
    double drdy = sin(phi)*sin(theta);
    double drdz = cos(theta);
    double dthetadx =  cos(phi)*cos(theta)/r;
    double dthetady =  sin(phi)*cos(theta)/r;
    double dthetadz = -sin(theta)/r;

    model->drhodx->data[idx] = drhodphi*dphidx + drhodr*drdx + drhodtheta*dthetadx;
    model->drhody->data[idx] = drhodphi*dphidy + drhodr*drdy + drhodtheta*dthetady;
    model->drhodz->data[idx] = drhodphi*dphidz + drhodr*drdz + drhodtheta*dthetadz;

    idx++;
            }
        }
    }
}