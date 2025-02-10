#include <stdlib.h>
#include <stdio.h>
#include <math.h>

// #include "newmesh.h"
// #include "domain.h"


// TODO: read in Jupiter outputs
// Alternative: change Jupiter to fake fargo output
// Not worrying about temp for now

typedef struct Model {
    // a struct to hold info for the Many fields in a fargo model
    char fargodir[100];
    char nout[5]; // this is usually an int, but save as char* so it can be 'avg'
    size_t nx;
    size_t ny;
    size_t nz;

    double alpha;
    double aspect;
    double flaring;
    double planetmass;
    double omegaframe;
    double planetpos[3];
    double sunpos[3];
    double planetEnvelope;

    // data
    MeshField *gasdens;
    // MeshField *gasvphi;
    // MeshField *gasvr;
    // MeshField *gasvtheta;

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

void make_cartvels(Model *model, MeshField *gasvphi, MeshField *gasvr, MeshField *gasvtheta);
void init_gradrho(Model *model);
void get_planetVars(Variables *var, Model *model);
double get_soundspeed(Model *model, double r);
Model *init_fargo_Model(char* fargodir, char* nout, size_t nx, size_t ny, size_t nz);
Model *init_Jupiter_Model(char* fargodir, char* nout, size_t nx, size_t ny, size_t nz);
Model *init_RADMC_Model(char* fargodir, char*nout, size_t nx, size_t ny, size_t nz);

Model *init_Model(int which, char* fargodir, char* nout, size_t nx, size_t ny, size_t nz) {
    switch (which) {
        case FARGO_MODEL:
            return init_fargo_Model(fargodir, nout, nx, ny, nz);
        case JUPITER_MODEL:
            return init_Jupiter_Model(fargodir, nout, nx, ny, nz);
        case RADMC_MODEL:
            return init_RADMC_Model(fargodir, nout, nx, ny, nz);
        default:
            printf("Assuming model input is from Fargo\n");
            return init_fargo_Model(fargodir, nout, nx, ny, nz);
    }
}

Model *init_fargo_Model(char* fargodir, char* nout, size_t nx, size_t ny, size_t nz) {
    // allocate memory
    Model *model = (Model*)malloc(sizeof(*model));
    if (!model) {
        perror("Malloc Failed on Model Creation");
        exit(1);
    }
    // string copy the fargo directory
    snprintf(model->fargodir,100,"%s",fargodir);
    // store model params
    strcpy(model->nout, nout);
    model->nx = nx;
    model->ny = ny;
    model->nz = nz;
    // read in the files
    char rhofile[100];
    char vphifile[100];
    char vrfile[100];
    char vthetafile[100];
    char varfile[100];
    int cx;
    cx = snprintf(rhofile,100,"%s/gasdens%s.dat",fargodir,nout);
    cx = snprintf(vphifile,100,"%s/gasvx%s.dat",fargodir,nout);
    cx = snprintf(vrfile,100,"%s/gasvy%s.dat",fargodir,nout);
    cx = snprintf(vthetafile,100,"%s/gasvz%s.dat",fargodir,nout);
    cx = snprintf(varfile,100,"%s/variables.par",fargodir);
    if (cx>100) {
        perror("Fargodir is too long!");
        exit(1);
    }
    // make the MeshFields
    model->gasdens   = init_MeshField_fromFile(rhofile,nx,ny,nz,RHO);
    MeshField *gasvphi   = init_MeshField_fromFile(vphifile,nx,ny,nz,VPHI);
    MeshField *gasvr     = init_MeshField_fromFile(vrfile,nx,ny,nz,VR);
    MeshField *gasvtheta = init_MeshField_fromFile(vthetafile,nx,ny,nz,VTHETA);

    // initialize the domain
    model->domain = init_Domain(fargodir,nx,ny,nz);

    // get the cartesian gas velocities
    printf("Making cartvels...\n");
    make_cartvels(model, gasvphi, gasvr, gasvtheta);
    
    // make the density gradients in cartesian
    printf("Getting Gradrho...\n");
    init_gradrho(model);

    printf("Getting variables...\n");
    Variables *var = init_Variables_fromFile(varfile);
    printf("Getting planet variables...\n");
    get_planetVars(var,model);
    // read the variables and rescale where necessary
    model->alpha = get_value(var,"ALPHA");              // unitless
    model->aspect = get_value(var,"ASPECTRATIO");       // unitless
    model->flaring = get_value(var,"FLARINGINDEX");     // unitless
    model->planetmass = get_value(var,"PLANETMASS")     * MSUN;
    model->omegaframe = get_value(var,"PLANETROTFRAME") * 1/TIME;

    // get the planet(s) position(s)
    model->planetpos[0] = get_value(var,"PLANX") * R0;
    model->planetpos[1] = get_value(var,"PLANY") * R0;
    model->planetpos[2] = get_value(var,"PLANZ") * R0;

    // CoM is at origin;
    model->sunpos[0] = ( -model->planetpos[0] * model->planetmass/MSUN );
    model->sunpos[1] = ( -model->planetpos[1] * model->planetmass/MSUN );
    model->sunpos[2] = ( -model->planetpos[2] * model->planetmass/MSUN );

    // semi major axis of planet
    double sma = sqrt(model->planetpos[0]*model->planetpos[0] + model->planetpos[1]*model->planetpos[1]);
    double hillRadius = sma * pow(model->planetmass/3/MSUN,1.0/3.0);
    double soundspeed = get_soundspeed(model,sma);
    double bondiRadius = 2*G*model->planetmass/soundspeed/soundspeed;
    // planet enevlope is min(hillRadius/4, bondiRadius)
    model->planetEnvelope = (hillRadius/4. < bondiRadius) ? hillRadius/4 : bondiRadius;


    printf("Read in variables: \n");
    print_variables(var);

    free_Variables(var);

    free_MeshField(gasvphi);
    free_MeshField(gasvr);
    free_MeshField(gasvtheta);

    return model;
}

Model *init_Jupiter_Model(char* fargodir, char* nout, size_t nx, size_t ny, size_t nz) {
    Model *model = (Model*)malloc(sizeof(*model));
    if (!model) {
        perror("Malloc Failed on Model Creation");
        exit(1);
    }
    // string copy the fargo directory
    snprintf(model->fargodir,100,"%s",fargodir);
    // store model params
    strcpy(model->nout, nout);
    model->nx = nx;
    model->ny = ny;
    model->nz = nz;
    // read in the files
    char rhofile[100];
    char vxfile[100];
    char vyfile[100];
    char vzfile[100];
    int cx;
    cx = snprintf(rhofile,100,"%s/gasdens%s.dat",fargodir,nout);
    cx = snprintf(vxfile,100,"%s/gasvx%s.dat",fargodir,nout);
    cx = snprintf(vyfile,100,"%s/gasvy%s.dat",fargodir,nout);
    cx = snprintf(vzfile,100,"%s/gasvz%s.dat",fargodir,nout);
    if (cx>100) {
        perror("Fargodir is too long!");
        exit(1);
    }
    // make the MeshFields
    // do not rescale any of the values
    model->gasdens = init_MeshField_fromFile(rhofile,nx,ny,nz,0);
    // read the vx, vy, and vz directly
    model->gasvx   = init_MeshField_fromFile(vxfile,nx,ny,nz,0);
    model->gasvy   = init_MeshField_fromFile(vyfile,nx,ny,nz,0);
    model->gasvz   = init_MeshField_fromFile(vzfile,nx,ny,nz,0);

    // initialize the domain
    model->domain = init_Jupiter_Domain(fargodir,nx,ny,nz);

    // get the gradients
    init_gradrho(model);

    // lets just fill these in for now
    // TODO: read these in
    model->alpha = 1.0e-3;
    model->aspect = 0.05;
    model->flaring = 0.25;
    model->planetmass = 1.e-3 * MSUN;
    model->omegaframe = 1/TIME;

     // get the planet(s) position(s)
    model->planetpos[0] = R0;
    model->planetpos[1] = 0.0;
    model->planetpos[2] = 0.0;

    // CoM is at origin;
    model->sunpos[0] = ( -model->planetpos[0] * model->planetmass/MSUN );
    model->sunpos[1] = ( -model->planetpos[1] * model->planetmass/MSUN );
    model->sunpos[2] = ( -model->planetpos[2] * model->planetmass/MSUN );

    // semi major axis of planet
    double sma = sqrt(model->planetpos[0]*model->planetpos[0] + model->planetpos[1]*model->planetpos[1]);
    double hillRadius = sma * pow(model->planetmass/3/MSUN,1.0/3.0);
    double soundspeed = get_soundspeed(model,sma);
    double bondiRadius = 2*G*model->planetmass/soundspeed/soundspeed;
    // planet enevlope is min(hillRadius/4, bondiRadius)
    model->planetEnvelope = (hillRadius/4. < bondiRadius) ? hillRadius/4 : bondiRadius;

    return model;
}

Model *init_RADMC_Model(char* fargodir, char* nout, size_t nx, size_t ny, size_t nz) {
    Model *model = (Model*)malloc(sizeof(*model));
    if (!model) {
        perror("Malloc Failed on Model Creation");
        exit(1);
    }
    // string copy the fargo directory
    snprintf(model->fargodir,100,"%s",fargodir);
    // store model params
    strcpy(model->nout, nout);
    model->nx = nx;
    model->ny = ny;
    model->nz = nz;
    // read in the files
    char rhofile[100];
    char vphifile[100];
    char vrfile[100];
    char vthetafile[100];
    char varfile[100];
    int cx;
    cx = snprintf(rhofile,100,"%s/gasdens%s.dat",fargodir,nout);
    cx = snprintf(vphifile,100,"%s/gasvx%s.dat",fargodir,nout);
    cx = snprintf(vrfile,100,"%s/gasvy%s.dat",fargodir,nout);
    cx = snprintf(vthetafile,100,"%s/gasvz%s.dat",fargodir,nout);
    cx = snprintf(varfile,100,"%s/variables.par",fargodir);
    if (cx>100) {
        perror("Fargodir is too long!");
        exit(1);
    }
    // make the MeshFields
    // do not rescale any of the values
    model->gasdens = init_MeshField_fromFile(rhofile,nx,ny,nz,0);
    MeshField *gasvphi   = init_MeshField_fromFile(vphifile,nx,ny,nz,0);
    MeshField *gasvr     = init_MeshField_fromFile(vrfile,nx,ny,nz,0);
    MeshField *gasvtheta = init_MeshField_fromFile(vthetafile,nx,ny,nz,0);

    // initialize the domain
    // reuse Jupiter domain because we have no ghost cells and no rescaling
    model->domain = init_Jupiter_Domain(fargodir,nx,ny,nz);

    // get the cartesian gas velocities
    printf("Making cartvels...\n");
    make_cartvels(model, gasvphi, gasvr, gasvtheta);

    // get the gradients
    init_gradrho(model);

    printf("Getting variables...\n");
    Variables *var = init_Variables_fromFile(varfile);
    printf("Getting planet variables...\n");
    get_planetVars(var,model);
    // read the variables and rescale where necessary
    model->alpha = get_value(var,"ALPHA");               // unitless
    model->aspect = get_value(var,"ASPECTRATIO");        // unitless
    model->flaring = get_value(var,"FLARINGINDEX");      // unitless
    model->planetmass = get_value(var,"PLANETMASS");     // do not rescale
    model->omegaframe = get_value(var,"PLANETROTFRAME"); // do not rescale

    // get the planet(s) position(s)
    model->planetpos[0] = get_value(var,"PLANX");
    model->planetpos[1] = get_value(var,"PLANY");
    model->planetpos[2] = get_value(var,"PLANZ");

    // CoM is at origin;
    model->sunpos[0] = ( -model->planetpos[0] * model->planetmass/MSUN );
    model->sunpos[1] = ( -model->planetpos[1] * model->planetmass/MSUN );
    model->sunpos[2] = ( -model->planetpos[2] * model->planetmass/MSUN );

    // semi major axis of planet
    double sma = sqrt(model->planetpos[0]*model->planetpos[0] + model->planetpos[1]*model->planetpos[1]);
    double hillRadius = sma * pow(model->planetmass/3/MSUN,1.0/3.0);
    double soundspeed = get_soundspeed(model,sma);
    double bondiRadius = 2*G*model->planetmass/soundspeed/soundspeed;
    // planet enevlope is min(hillRadius/4, bondiRadius)
    model->planetEnvelope = (hillRadius/4. < bondiRadius) ? hillRadius/4 : bondiRadius;

    return model;
}

void free_Model(Model *model) {
    free_MeshField(model->gasdens);
    // free_MeshField(model->gasvphi);
    // free_MeshField(model->gasvr);
    // free_MeshField(model->gasvtheta);

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

void free_Models(Model **models, int nlvl) {
    for (int i=0; i<nlvl; i++) {
        free_Model(models[i]);
    }
    // free (models);
}

void make_cartvels(Model *model, MeshField *gasvphi, MeshField *gasvr, MeshField *gasvtheta) {
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

                vphi   = gasvphi->data[idx];
                vr     = gasvr->data[idx];
                vtheta = gasvtheta->data[idx];

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
    // case VPHI:
    //     data_values = model->gasvphi;
    //     break;
    // case VR:
    //     data_values = model->gasvr;
    //     break;
    // case VTHETA:
    //     data_values = model->gasvtheta;
    //     break;
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

double get_scaleheight(Model *model, double r) {
    double ASPECT = model->aspect;
    double FLARING = model->flaring;
    return r*ASPECT*pow(r/R0,FLARING);
}

double get_soundspeed(Model *model, double r) {
    double H = get_scaleheight(model, r);
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

void get_planetVars(Variables *var, Model *model) {
    char planetfile[100];
    int cx;
    cx = snprintf(planetfile,100,"%s/planet0.dat",model->fargodir);
    if (cx>100) {
        perror("Fargo dir too long!");
        exit(1);
    }
    FILE *file;
    file = fopen(planetfile,"r");
    if (file == NULL) {
        perror("Cannot open planet file");
        exit(1);
    }
    char* line = NULL;
    size_t len=0;
    ssize_t read=0;
    while ((read = getline(&line, &len, file)) != -1) {
        char *split_str;
        char *nout, *x, *y, *z, *vx, *vy, *vz, *mass, *time, *rotframe;
        split_str = strtok(line, " \t\n");
        size_t idx = 0;
        while (split_str != NULL) {
            switch (idx)
            {
            case 0:
                nout = split_str;
                break;
            case 1:
                x = split_str;
                break;
            case 2:
                y = split_str;
                break;
            case 3:
                z = split_str;
                break;
            case 4:
                vx = split_str;
                break;
            case 5:
                vy = split_str;
                break;
            case 6:
                vz = split_str;
                break;
            case 7:
                mass = split_str;
                break;
            case 8:
                time = split_str;
                break;
            case 9:
                rotframe = split_str;
                break;
            default:
                break;
            }
            idx++;
            split_str = strtok(NULL," \t\n");
        }
        if (strcmp(nout,model->nout)==0) {
            add_variable(var, "PLANX", atof(x));
            add_variable(var, "PLANY", atof(y));
            add_variable(var, "PLANZ", atof(z));
            add_variable(var, "PLANVX", atof(vx));
            add_variable(var, "PLANVY", atof(vy));
            add_variable(var, "PLANVZ", atof(vz));
            add_variable(var, "PLANETMASS", atof(mass));
            add_variable(var, "PLANETTIME", atof(time));
            add_variable(var, "PLANETROTFRAME", atof(rotframe));
            break;
        }
    }
}