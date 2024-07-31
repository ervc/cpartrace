#include "partrace.h"

// defaults
#define NX 2048
#define NY 256
#define NZ 32

int solve_grid(Model *model, char outfile[100], double s);

int main(int argc, char **argv) {
    double s;
    char fargodir[100];
    char nout[5];
    char outfile[100];
    if (argc==1) {
        // DEFAULTS
        s = 0.01;
        strcpy(fargodir,"/Users/ericvc/fargo/outputs/alpha3_mplan300");
        strcpy(nout,"avg");
        strcpy(outfile,"outputs/equilibrium_velocities.dat");
    } else if (argc>=4) {
        s = atof(argv[1]);
        strcpy(fargodir,argv[2]);
        strcpy(nout,argv[3]);
        if (argc==5) {
            strcpy(outfile,argv[4]);
        } else {
            strcpy(outfile,"outputs/equilibrium_velocities.dat");
        }
    } else {
        printf(
"Please provide grainsize, fargodirectory, nout [and optional outfile] or nothing for defaults\n\
    USAGE: ./find_grain_velocities [grainsize fargodir nout [outfile]]\n\
Note: this is all or nothing cannot provide only some values\n");
        exit(EXIT_FAILURE);
    }
    size_t nx = NX;
    size_t ny = NY;
    size_t nz = NZ;

    Model *model = init_Model(fargodir,nout,nx,ny,nz);

    if ( solve_grid(model,outfile,s) != EXIT_SUCCESS) {
        printf("Not successful :(\n");
    }
}

/**
 * @brief Get the equilibrium velocity at a given location
 * 
 * @param model IN - partrace gas model
 * @param x IN - X position (cm)
 * @param y IN - y position (cm)
 * @param z IN - z position (cm)
 * @param s IN - grain size (cm)
 * @param vxout OUT - equilibrium x velocity (cm/s)
 * @param vyout OUT - equilibrium y velocity (cm/s)
 * @param vzout OUT - equilibrium z velocity (cm/s)
 * @return int 0 success; 1 failure
 */
int get_equilibrium_velocity(Model *model,
        double x, double y, double z, double s,
        double *vxout, double *vyout, double *vzout) {
    double omegaframe = model->omegaframe;
    double rho_s = PARTDENSITY;

    double GMSUN = G*MSUN;
    double GMPLAN = G*model->planetmass;

    double r = sqrt(x*x + y*y + z*z);
    double omega = sqrt(GMSUN/r)/r;
    double scaleheight = model->aspect*r*pow(r/R0,model->flaring);
    double cs = scaleheight*omega;
    double rho_g = cart_trilinterp_one(model,RHO,x,y,z);
    double tstop = (rho_s*s)/(rho_g*cs);
    double invts = 1/tstop;
    double d = sqrt( (x-model->sunpos[0])*(x-model->sunpos[0])
        + (y-model->sunpos[1])*(y-model->sunpos[1])
        + (z-model->sunpos[2])*(z-model->sunpos[2]) );
    double Gsun = GMSUN/d/d/d;
    double Gx = Gsun*(x-model->sunpos[0]);
    double Gy = Gsun*(y-model->sunpos[1]);
    double Gz = Gsun*(z-model->sunpos[2]);
    double dp = sqrt( (x-model->planetpos[0])*(x-model->planetpos[0])
        + (y-model->planetpos[1])*(y-model->planetpos[1])
        + (z-model->planetpos[2])*(z-model->planetpos[2]) );
    if (dp < model->planetEnvelope) {
        *vxout = 0.0;
        *vyout = 0.0;
        *vzout = 0.0;
        return EXIT_SUCCESS;
    }
    double Gplan = GMPLAN/dp/dp/dp;
    double Gpx = Gplan*(x-model->planetpos[0]);
    double Gpy = Gplan*(y-model->planetpos[1]);
    double Gpz = Gplan*(z-model->planetpos[2]);

    double gasvx = cart_trilinterp_one(model,VX,x,y,z);
    double gasvy = cart_trilinterp_one(model,VY,x,y,z);
    double gasvz = cart_trilinterp_one(model,VZ,x,y,z);

    // initial guesses
    double vx = gasvx;
    double vy = gasvy;
    double vz = gasvz;

    //solver
    int status = -1;
    int step = 0;
    int maxN = 100;
    double tol = 1.e-10;
    double ax,ay,az;
    while (status==-1) {
        step++;
        double vphi = (x*vy - y*vx)/r;
        ax = (invts*gasvx - invts*vx
            - Gx - Gpx
            + 2*vy*omegaframe + x*omegaframe*omegaframe
            + vphi*vphi/r/r*x
        );
        ay = (invts*gasvy - invts*vy
            - Gy - Gpy
            + -2*vx*omegaframe + y*omegaframe*omegaframe
            + vphi*vphi/r/r*y
        );
        az = (invts*gasvz - invts*vz
            -Gz - Gpz
        );
        if ((fabs(ax)<tol) && (fabs(ay)<tol) && (fabs(az)<tol)) {
            status = EXIT_SUCCESS;
            break;
        }
        vx += ax*tstop;
        vy += ay*tstop;
        vz += az*tstop;

        if (step >= maxN) {
            status=EXIT_FAILURE;
        }
        
    }
    *vxout = vx;
    *vyout = vy;
    *vzout = vz;

    return status;
}

int py_get_equilibrium_velocity(
        char fargodir[100], char nout[5],
        double x, double y, double z, double s,
        double *vxout, double *vyout, double *vzout)
        {
    printf("fargodir: %s\n",fargodir);
    Model *model = init_Model(fargodir,nout,NX,NY,NZ);

    return get_equilibrium_velocity(model,x,y,z,s,vxout,vyout,vzout);
}

int solve_grid(Model *model, char outfile[100], double s) {
    Domain *domain = model->domain;
    double omegaframe = model->omegaframe;
    double rho_s = PARTDENSITY;

    double GMSUN = G*MSUN;
    double GMPLAN = G*model->planetmass;


    size_t size = model->nx*model->ny*model->nz*3;
    // equilibrium velocity arrays
    double *eqvels = calloc(size,sizeof(double));

    size_t k0 = 0;
    size_t j0 = 0;
    size_t i0 = 0;
    size_t idx = get_idx(model->gasdens,k0,j0,i0);
    for (size_t k=k0; k<model->nz; k++) {
        double theta = domain->thetaCenters[k];
        for (size_t j=j0; j<model->ny; j++) {
            double r = domain->rCenters[j];
            double r3 = r*r*r;
            double omega = sqrt(GMSUN/r3);
            double vkep_norot = omega*r;
            double scaleheight = 0.05*r*pow(r/R0,0.25);
            double cs = scaleheight*omega;
            printf("Starting %zu, %zu\n",k,j);
            for (size_t i=i0; i<model->nx; i++) {
                double phi = domain->phiCenters[i];
                double x = r*cos(phi)*sin(theta);
                double y = r*sin(phi)*sin(theta);
                double z = r*cos(theta);
                double rho_g = model->gasdens->data[idx];
                double tstop = (rho_s*s)/(rho_g*cs);
                double Stokes = tstop*omega;
                double invts = 1/tstop;
                double d = sqrt( (x-model->sunpos[0])*(x-model->sunpos[0])
                    + (y-model->sunpos[1])*(y-model->sunpos[1])
                    + (z-model->sunpos[2])*(z-model->sunpos[2]) );
                double Gsun = GMSUN/d/d/d;
                double Gx = Gsun*(x-model->sunpos[0]);
                double Gy = Gsun*(y-model->sunpos[1]);
                double Gz = Gsun*(z-model->sunpos[2]);
                double dp = sqrt( (x-model->planetpos[0])*(x-model->planetpos[0])
                    + (y-model->planetpos[1])*(y-model->planetpos[1])
                    + (z-model->planetpos[2])*(z-model->planetpos[2]) );
                if (dp < model->planetEnvelope) {
                    eqvels[3*idx] = 0.0;
                    eqvels[3*idx+1] = 0.0;
                    eqvels[3*idx+2] = 0.0;
                    idx++;
                    break;
                }
                double Gplan = GMPLAN/dp/dp/dp;
                double Gpx = Gplan*(x-model->planetpos[0]);
                double Gpy = Gplan*(y-model->planetpos[1]);
                double Gpz = Gplan*(z-model->planetpos[2]);

                // get the gas velocities
                double gasvx = model->gasvx->data[idx];
                double gasvy = model->gasvy->data[idx];
                double gasvz = model->gasvz->data[idx];
                double gasvp = (x*gasvy - y*gasvx)/r + r*omegaframe;
                double gasvr = (x*gasvx + y*gasvy)/r;
                double eta = 1-(gasvp*gasvp/vkep_norot/vkep_norot);
                
                // armitage approximations for no planet, smooth disk
                double armvr = (gasvr - vkep_norot*eta*Stokes)/(1+Stokes*Stokes);
                double armvphi = gasvp - 0.5*Stokes*armvr - r*omegaframe;
                double armvx = armvr*cos(phi) - armvphi*sin(phi);
                double armvy = armvr*sin(phi) + armvphi*cos(phi);
                double acent = armvphi*armvphi/r/r;

                // initial guesses
                double vx = gasvx + tstop*(
                    -Gx - Gpx 
                    + 2*armvy*omegaframe + x*omegaframe*omegaframe
                    + x*acent);
                double vy = gasvy + tstop*(
                    -Gy - Gpy
                    + -2*armvx*omegaframe + y*omegaframe*omegaframe
                    + y*acent
                );
                double vz = gasvz + tstop*(
                    -Gz - Gpz
                );

                // printf("phi/pi, r/AU, z/r = %e, %e, %e\n",phi/M_PI, r/AU, z/r);
                // printf("tstop, rho, cs = %e, %e, %e\n",tstop,rho_g,cs);
                // printf("Initial guess:\n  %e\n  %e\n  %e\n",vx,vy,vz);
                // return EXIT_FAILURE;

                // write a solver
                int status = 0;
                int step = 0;
                int maxN = 100000;
                double tol = 1.e-10;
                double ax,ay,az;
                while (status==0) {
                    step++;
                    double vphi = (x*vy - y*vx)/r;
                    ax = (invts*gasvx - invts*vx
                        - Gx - Gpx
                        + 2*vy*omegaframe + x*omegaframe*omegaframe
                        + vphi*vphi/r/r*x
                    );
                    ay = (invts*gasvy - invts*vy
                        - Gy - Gpy
                        + -2*vx*omegaframe + y*omegaframe*omegaframe
                        + vphi*vphi/r/r*y
                    );
                    az = (invts*gasvz - invts*vz
                        -Gz - Gpz
                    );
                    if ((fabs(ax)<tol) && (fabs(ay)<tol) && (fabs(az)<tol)) {
                        status = 1;
                        break;
                    }
                    vx += ax*tstop*0.01;
                    vy += ay*tstop*0.01;
                    vz += az*tstop*0.01;

                    if (step >= maxN) {
                        printf("Reached maxN: %zu, %zu, %zu\n",k,j,i);
                        status=-1;
                    }
                }
                // if (status!=1){
                //     return EXIT_FAILURE;
                // }
                // printf("Final Guess:\n  %e\n  %e\n  %e\n",vx,vy,vz);
                eqvels[3*idx] = vx;
                eqvels[3*idx+1] = vy;
                eqvels[3*idx+2] = vz;

                idx++;
            }
        }
    }

    FILE *eqvelfile;
    eqvelfile = fopen(outfile,"wb");
    fwrite(eqvels,sizeof(double),size,eqvelfile);
    return EXIT_SUCCESS;
}

int py_solve_grid(char fargodir[100], char nout[5], char outfile[100],
    double s) {
        printf("fargodir: %s\n",fargodir);
        Model *model = init_Model(fargodir,nout,NX,NY,NZ);

        return solve_grid(model,outfile,s);
    }