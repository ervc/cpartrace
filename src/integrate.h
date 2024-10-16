// integrator status values
#define COMPLETE    1
#define OOB         2
#define ACCRETED    3
#define CROSSED     4

// if Particle crosses this radius then it is crossed from outer to inner
#define CROSS_R  ( R0 )
// Define the stokes number where particles are fully coupled with the gas
#define COUPLING  ( 1.0e-3 )
// diffusion boolean
// #define DIFFUSION   1

/**
 * @brief Returns a random double between 0.0 and 1.0
 * 
 * @return double 
 */
double random_double() {
    return ((double)rand()) / ((double)(RAND_MAX));
}

/**
 * @brief Return a random double between a and b
 * 
 * @param a 
 * @param b 
 * @return double 
 */
double random_range(double a, double b) {
    double R = random_double();
    return (b-a)*R + a;
}

void get_gradrho(Particle *particle, double Y[6], double gradrho[3]) {
    double x,y,z;
    double phi,r,theta;
    Model *model = particle->model;

    // get things we'll need from the particle
    x = Y[0];
    y = Y[1];
    z = Y[2];

    // find the "lower left" corner for cube containing our value
    // deal with negative z values here
    // deal with phi between inside and outside edge in the interpolation
    int zflag = 0;
    if (z<0) {
        zflag=1;
        z=-z;
    }
    // get velocity and acceleration as though we are in the top of the disk
    phi = atan2(y,x);
    r = sqrt(x*x + y*y + z*z);
    theta = acos(z/r);

    double drhodx, drhody, drhodz;
    // we've already checked for negative z, so only way this is the case is if
    // we're right near the midplane. Due to symmetry, values are the same above
    // and below midplane, so no need to interpolate, just take the midplane values

    // interpolate for values
    if (theta > model->domain->thetaCenters[model->nz-1]) {
        drhodx = midplane_interp(model, DRHODX, phi, r);
        drhody = midplane_interp(model, DRHODY, phi, r);
        drhodz = midplane_interp(model, DRHODZ, phi, r);
    } else {
        size_t corner[3];
        size_t i,j,k;
        // these values don't change so let's just find them once
        get_corner(model,phi,r,theta,corner);
        i=corner[0]; j=corner[1]; k=corner[2];
        drhodx = fast_linterp(model, DRHODX, i,j,k, phi, r, theta);
        drhody = fast_linterp(model, DRHODY, i,j,k, phi, r, theta);
        drhodz = fast_linterp(model, DRHODZ, i,j,k, phi, r, theta);
    }
    if (zflag) {drhodz = -drhodz;}
    gradrho[0]=drhodx; gradrho[1]=drhody; gradrho[2]=drhodz;
}

void get_vrho(Particle *particle, double Y[6], double vrho[3]) {
    double x,y,z;
    double phi,r,theta;
    Model *model = particle->model;

    // get things we'll need from the particle
    x = Y[0];
    y = Y[1];
    z = Y[2];

    // find the "lower left" corner for cube containing our value
    // deal with negative z values here
    // deal with phi between inside and outside edge in the interpolation
    int zflag = 0;
    if (z<0) {
        zflag=1;
        z=-z;
    }
    phi = atan2(y,x);
    r = sqrt(x*x + y*y + z*z);
    theta = acos(z/r);
    double r3 = r*r*r;

    double rho_g;
     // interpolate for values
    if (theta > model->domain->thetaCenters[model->nz-1]) {
        rho_g  = midplane_interp(model, RHO, phi, r);
    } else {
        size_t corner[3];
        size_t i,j,k;
        // these values don't change so let's just find them once
        get_corner(model,phi,r,theta,corner);
        i=corner[0]; j=corner[1]; k=corner[2];
        rho_g  = fast_linterp(model, RHO, i,j,k, phi, r, theta);
    }

    // scalars
    double ASPECT = model->aspect;
    double FLARING = model->flaring;
    double H = r*ASPECT*pow(r/R0,FLARING);
    double OM = sqrt(G*MSUN/r3);
    double ALPHA = model->alpha;
    double gasdiff = ALPHA*H*H*OM;
    double St = particle->size*particle->density/H/rho_g;
    double invSt2 = 1/(1+St*St);
    double partdiff = gasdiff*invSt2;

    double gradrho[3];
    double drhodx, drhody, drhodz;
    get_gradrho(particle, Y, gradrho);
    drhodx=gradrho[0]; drhody=gradrho[1]; drhodz=gradrho[2];

    vrho[0] = partdiff/rho_g*drhodx;
    vrho[1] = partdiff/rho_g*drhody;
    double vrhoz = partdiff/rho_g*drhodz;
    if (zflag) {vrhoz = -vrhoz;};
    vrho[2] = vrhoz;
}

double get_3mag(double V[3]) {
    double vx,vy,vz;
    vx=V[0]; vy=V[1]; vz=V[2];
    return sqrt(vx*vx + vy*vy + vz*vz);
}

double get_vrhomag(Particle *particle, double Y[6]) {
    double vrho[3];
    get_vrho(particle, Y, vrho);
    return get_3mag(vrho);
}

void get_vdiff(Particle *particle, double Y[6], double vdiff[3]) {
    double x,y,z;
    double phi,r,theta;
    Model *model = particle->model;

    // get things we'll need from the particle
    x = Y[0];
    y = Y[1];
    z = Y[2];

    // find the "lower left" corner for cube containing our value
    // deal with negative z values here
    // deal with phi between inside and outside edge in the interpolation
    int zflag = 0;
    if (z<0) {
        zflag=1;
        z=-z;
    }
    // get velocity and acceleration as though we are in the top of the disk
    phi = atan2(y,x);
    r = sqrt(x*x + y*y + z*z);
    theta = acos(z/r);
    double r3 = r*r*r;
    double invr = 1./r;

    double rho_g, drhodx, drhody, drhodz;
    // we've already checked for negative z, so only way this is the case is if
    // we're right near the midplane. Due to symmetry, values are the same above
    // and below midplane, so no need to interpolate, just take the midplane values

    // we have some useful function defined already, but I'm not calling them here
    // to save some time by not recalculating things over and over again

    // interpolate for values
    if (theta > model->domain->thetaCenters[model->nz-1]) {
        rho_g  = midplane_interp(model, RHO, phi, r);
        drhodx = midplane_interp(model, DRHODX, phi, r);
        drhody = midplane_interp(model, DRHODY, phi, r);
        drhodz = midplane_interp(model, DRHODZ, phi, r);
    } else {
        size_t corner[3];
        size_t i,j,k;
        // these values don't change so let's just find them once
        get_corner(model,phi,r,theta,corner);
        i=corner[0]; j=corner[1]; k=corner[2];
        rho_g  = fast_linterp(model, RHO, i,j,k, phi, r, theta);
        drhodx = fast_linterp(model, DRHODX, i,j,k, phi, r, theta);
        drhody = fast_linterp(model, DRHODY, i,j,k, phi, r, theta);
        drhodz = fast_linterp(model, DRHODZ, i,j,k, phi, r, theta);
    }
    

    // scalars
    double ASPECT = model->aspect;
    double FLARING = model->flaring;
    double H = r*ASPECT*pow(r/R0,FLARING);
    double OM = sqrt(G*MSUN/r3);
    double ALPHA = model->alpha;
    double gasdiff = ALPHA*H*H*OM;
    double St = particle->size*particle->density/H/rho_g;
    double invSt2 = 1.0/(1.0+St*St);

    // derivatives

    // chain rules
    // double dphidx = -sin(phi)/r/sin(theta);
    // double dphidy =  cos(phi)/r/sin(theta);
    // double dphidz =  0;
    double drdx = cos(phi)*sin(theta);
    double drdy = sin(phi)*sin(theta);
    double drdz = cos(theta);
    // double dthetadx =  cos(phi)*cos(theta)/r;
    // double dthetady =  sin(phi)*cos(theta)/r;
    // double dthetadz = -sin(theta)/r;

    // dHdphi = 0
    double dHdr = (FLARING+1.0)*ASPECT*pow(r/R0,FLARING);
    // dHdtheta = 0
    double dHdx = dHdr*drdx;
    double dHdy = dHdr*drdy;
    double dHdz = dHdr*drdz;
    // dOMdphi = 0;
    double dOMdr = -3.0/2.0*OM*invr;
    // dOMdtheta = 0;
    double dOMdx = dOMdr*drdx;
    double dOMdy = dOMdr*drdy;
    double dOMdz = dOMdr*drdz;
    double dStdx = -St/H*dHdx - St/rho_g*drhodx;
    double dStdy = -St/H*dHdy - St/rho_g*drhody;
    double dStdz = -St/H*dHdz - St/rho_g*drhodz;
    double dgasdiffdx = 2.*ALPHA*H*OM*dHdx + ALPHA*H*H*dOMdx;
    double dgasdiffdy = 2.*ALPHA*H*OM*dHdy + ALPHA*H*H*dOMdy;
    double dgasdiffdz = 2.*ALPHA*H*OM*dHdz + ALPHA*H*H*dOMdz;
    vdiff[0] = invSt2*dgasdiffdx - 2.*gasdiff*invSt2*invSt2*dStdx;
    vdiff[1] = invSt2*dgasdiffdy - 2.*gasdiff*invSt2*invSt2*dStdy;
    double dpartdiffdz = invSt2*dgasdiffdz - 2.*gasdiff*invSt2*invSt2*dStdz;
    if (zflag) {dpartdiffdz = -dpartdiffdz;}
    vdiff[2] = dpartdiffdz;
}

double get_vdiffmag(Particle *particle, double Y[6]) {
    double vdiff[3];
    get_vdiff(particle, Y, vdiff);
    return get_3mag(vdiff);
}

/**
 * @brief Give the effective velocity and acceleration for a particle with
 * a given position-velocity vector (not necessarily at the particle's location)
 * 
 * @param particle 
 * @param Y 
 * @param derivative 
 */
void dYdt(Particle *particle, double Y[6], double derivative[6]) {
    double x,y,z,vx,vy,vz;
    double phi,r,theta;
    // TODO: Here
    // if x,y,z are inside level 4: model = particle->model4
    // elif x,y,z are inside level 3: model = particle->model3
    // ...
    Model *model = particle->model;

    // get things we'll need from the particle
    x = Y[0];
    y = Y[1];
    z = Y[2];
    vx = Y[3];
    vy = Y[4];
    vz = Y[5];

    // find the "lower left" corner for cube containing our value
    // deal with negative z values here
    // deal with phi between inside and outside edge in the interpolation
    int zflag = 0;
    if (z<0) {
        zflag=1;
        z=-z;
        vz=-vz;
    }
    // get velocity and acceleration as though we are in the top of the disk
    phi = atan2(y,x);
    r = sqrt(x*x + y*y + z*z);
    theta = acos(z/r);
    double r3 = r*r*r;
    double invr = 1./r;

    double rho_g, gasvx, gasvy, gasvz, drhodx, drhody, drhodz;
    // we've already checked for negative z, so only way this is the case is if
    // we're right near the midplane. Due to symmetry, values are the same above
    // and below midplane, so no need to interpolate, just take the midplane values

    // we have some useful function defined already, but I'm not calling them here
    // to save some time by not recalculating things over and over again

    // interpolate for values
    if (theta > model->domain->thetaCenters[model->nz-1]) {
        rho_g  = midplane_interp(model, RHO, phi, r);
        gasvx  = midplane_interp(model, VX , phi, r);
        gasvy  = midplane_interp(model, VY , phi, r);
        gasvz  = midplane_interp(model, VZ , phi, r);
        drhodx = midplane_interp(model, DRHODX, phi, r);
        drhody = midplane_interp(model, DRHODY, phi, r);
        drhodz = midplane_interp(model, DRHODZ, phi, r);
    } else {
        size_t corner[3];
        size_t i,j,k;
        // these values don't change so let's just find them once
        get_corner(model,phi,r,theta,corner);
        i=corner[0]; j=corner[1]; k=corner[2];
        rho_g  = fast_linterp(model, RHO, i,j,k, phi, r, theta);
        gasvx  = fast_linterp(model, VX , i,j,k, phi, r, theta);
        gasvy  = fast_linterp(model, VY , i,j,k, phi, r, theta);
        gasvz  = fast_linterp(model, VZ , i,j,k, phi, r, theta);
        drhodx = fast_linterp(model, DRHODX, i,j,k, phi, r, theta);
        drhody = fast_linterp(model, DRHODY, i,j,k, phi, r, theta);
        drhodz = fast_linterp(model, DRHODZ, i,j,k, phi, r, theta);
    }
    

    // scalars
    double ASPECT = model->aspect;
    double FLARING = model->flaring;
    double H = r*ASPECT*pow(r/R0,FLARING);
    double OM = sqrt(G*MSUN/r3);
    double ALPHA = model->alpha;
    double gasdiff = ALPHA*H*H*OM;
    double St = particle->size*particle->density/H/rho_g;
    double invSt2 = 1./(1.+St*St);
    double partdiff = gasdiff*invSt2;

    if (St < COUPLING) {
        // if stokes number is low then the particle is perfectly coupled with the gas
        // set the v_eff as gas velocity
        // We will also change the particle velocity in the update step
        derivative[0] = gasvx;
        derivative[1] = gasvy;
        if (zflag) {gasvz = -gasvz;}
        derivative[2] = gasvz;
        derivative[3] = 0;
        derivative[4] = 0;
        derivative[5] = 0;
        return;
    }

    // derivatives

    // chain rules
    // double dphidx = -sin(phi)/r/sin(theta);
    // double dphidy =  cos(phi)/r/sin(theta);
    // double dphidz =  0;
    double drdx = cos(phi)*sin(theta);
    double drdy = sin(phi)*sin(theta);
    double drdz = cos(theta);
    // double dthetadx =  cos(phi)*cos(theta)/r;
    // double dthetady =  sin(phi)*cos(theta)/r;
    // double dthetadz = -sin(theta)/r;

    // dHdphi = 0
    double dHdr = (FLARING+1.)*ASPECT*pow(r/R0,FLARING);
    // dHdtheta = 0
    double dHdx = dHdr*drdx;
    double dHdy = dHdr*drdy;
    double dHdz = dHdr*drdz;
    // dOMdphi = 0;
    double dOMdr = -3./2.*OM*invr;
    // dOMdtheta = 0;
    double dOMdx = dOMdr*drdx;
    double dOMdy = dOMdr*drdy;
    double dOMdz = dOMdr*drdz;
    double dStdx = -St/H*dHdx - St/rho_g*drhodx;
    double dStdy = -St/H*dHdy - St/rho_g*drhody;
    double dStdz = -St/H*dHdz - St/rho_g*drhodz;
    double dgasdiffdx = 2.*ALPHA*H*OM*dHdx + ALPHA*H*H*dOMdx;
    double dgasdiffdy = 2.*ALPHA*H*OM*dHdy + ALPHA*H*H*dOMdy;
    double dgasdiffdz = 2.*ALPHA*H*OM*dHdz + ALPHA*H*H*dOMdz;
    double dpartdiffdx = invSt2*dgasdiffdx - 2.*gasdiff*invSt2*invSt2*dStdx;
    double dpartdiffdy = invSt2*dgasdiffdy - 2.*gasdiff*invSt2*invSt2*dStdy;
    double dpartdiffdz = invSt2*dgasdiffdz - 2.*gasdiff*invSt2*invSt2*dStdz;

    // Veff = V + D/rho_g grad(rho) + grad(D)
    double veffx = vx + partdiff/rho_g*drhodx + dpartdiffdx;
    double veffy = vy + partdiff/rho_g*drhody + dpartdiffdy;
    double veffz = vz + partdiff/rho_g*drhodz + dpartdiffdz;

    double dsun = sqrt((x-model->sunpos[0])*(x-model->sunpos[0])
                     + (y-model->sunpos[1])*(y-model->sunpos[1]) 
                     + (z-model->sunpos[2])*(z-model->sunpos[2]));
    double dplan = sqrt((x-model->planetpos[0])*(x-model->planetpos[0])
                      + (y-model->planetpos[1])*(y-model->planetpos[1])
                      + (z-model->planetpos[2])*(z-model->planetpos[2]));
    double dsun3 = dsun*dsun*dsun;
    double dplan3 = dplan*dplan*dplan;

    double agravmag = -G*MSUN/dsun3;
    double agravx = agravmag * (x-model->sunpos[0]);
    double agravy = agravmag * (y-model->sunpos[1]);
    double agravz = agravmag * (z-model->sunpos[2]);
    double aplanmag = -G*model->planetmass/dplan3;
    double aplanx = aplanmag * (x-model->planetpos[0]);
    double aplany = aplanmag * (y-model->planetpos[1]);
    double aplanz = aplanmag * (z-model->planetpos[2]);
    double adragmag = -OM/St;
    double adragx = adragmag * (vx-gasvx);
    double adragy = adragmag * (vy-gasvy);
    double adragz = adragmag * (vz-gasvz);
    double OMEGAFRAME = model->omegaframe;
    double arotx =  2.*vy*OMEGAFRAME + x*OMEGAFRAME*OMEGAFRAME;
    double aroty = -2.*vx*OMEGAFRAME + y*OMEGAFRAME*OMEGAFRAME;
    double arotz =  0.;

    // total acceleration
    double atotx = agravx + aplanx + adragx + arotx;
    double atoty = agravy + aplany + adragy + aroty;
    double atotz = agravz + aplanz + adragz + arotz;

    // if z was negative, flip the z velocity and acceleration
    if (zflag) {
        veffz = -veffz;
        atotz = -atotz;
    }

    derivative[0] = veffx;
    derivative[1] = veffy;
    derivative[2] = veffz;
    derivative[3] = atotx;
    derivative[4] = atoty;
    derivative[5] = atotz;
}

void add_6vectors(double Y[6], double dt, double deriv[6], double result[6]) {
    for (int i=0; i<6; i++) {
        result[i] = Y[i] + dt*deriv[i];
    }
}

/**
 * @brief Runge-Kutta 4 integration, gives new position and velocity vector
 * for a particle at a given location
 * 
 * @param particle 
 * @param h 
 * @param result 
 */
void rk4(Particle *particle, double h, double result[6]) {
    // derivative we will be using
    double k1[6];
    double k2[6];
    double k3[6];
    double k4[6];
    // 6d posvel vectors
    double Y0[6];
    double Y1[6];
    double Y2[6];
    double Y3[6];
    // initial position and velocity vector of the particle
    get_posvelVector(particle, Y0);
    // first step to get the derivative
    dYdt(particle, Y0, k1);
    // get the new location to evalutate from
    add_6vectors(Y0, h/2., k1, Y1);
    // second step
    dYdt(particle, Y1, k2);
    // etc...
    add_6vectors(Y0, h/2., k2, Y2);
    dYdt(particle, Y2, k3);
    
    add_6vectors(Y0, h, k3, Y3);
    dYdt(particle, Y3, k4);

    for (int i=0; i<6; i++) {
        result[i] = Y0[i] + 1.0/6.0*(k1[i]+2.*k2[i]+2.*k3[i]+k4[i])*h;
    }
}

/**
 * @brief Step the particle using the Runge-Kutta 4 integration. Add diffusion
 * if turned on.
 * 
 * @param particle 
 * @param dt 
 */
void rkstep_particle(Particle *particle, double dt, int DIFFUSION) {
    double result[6];
    rk4(particle, dt, result);
    double St = get_Stokes(particle);
    if (St < COUPLING) {
        double gasvx, gasvy, gasvz;
        double x,y,z,phi,r,theta;
        x = result[0];
        y = result[1];
        z = result[2];
        // interp for the gas velocities to give the particle
        // trilinterp_one function checks for -z so we don't have to
        phi = atan2(y,x);
        r = sqrt(x*x + y*y + z*z);
        theta = acos(z/r);
        gasvx = trilinterp_one(particle->model, VX, phi,r,theta);
        gasvy = trilinterp_one(particle->model, VY, phi,r,theta);
        gasvz = trilinterp_one(particle->model, VZ, phi,r,theta);
        result[3] = gasvx;
        result[4] = gasvy;
        result[5] = gasvz;
    }
    if (DIFFUSION) {
        double gradpartdiff[3];
        double Y[6];
        get_posvelVector(particle, Y);
        // note that vdiff = grad D
        get_vdiff(particle, Y, gradpartdiff);
        double x,y,z;
        x = particle->x;
        y = particle->y;
        z = particle->z;
        double xprime = x+0.5*gradpartdiff[0]*dt;
        double yprime = y+0.5*gradpartdiff[1]*dt;
        double zprime = z+0.5*gradpartdiff[2]*dt;
        // grad partdif at Yprime
        // we really only care about the position, so velocity can be 0
        double Yprime[6] = {xprime, yprime, zprime, 0, 0, 0};
        double Dxprime = get_partdiff(particle, Yprime);
        // 3 random numbers between -1. and 1.
        double Rx = random_range(-1.0,1.0);
        double Ry = random_range(-1.0,1.0);
        double Rz = random_range(-1.0,1.0);
        // 6.0 = 2.0/(1/3) = 2/std of distribution
        double randkick = sqrt(6.0*Dxprime*dt);
        result[0] += Rx*randkick;
        result[1] += Ry*randkick;
        result[2] += Rz*randkick;
    }
    update_Particle(particle, result);
}

double inv_sigmoid(double x) {
    double scale = 10.0;
    double center = log10(COUPLING);
    double xp = (x-center)*scale;
    return 1.0+exp(-xp);
}

/**
 * @brief Physical limitations on the timestep
 * 
 * @param particle 
 * @return double 
 */
double get_dt(Particle *particle) {
    double x,y,z;
    x=particle->x;
    y=particle->y;
    z=particle->z;
    double r = sqrt(x*x + y*y + z*z);
    double St = get_Stokes(particle);
    double OM = get_Omega(r);
    double tstop = St/OM * inv_sigmoid(log10(St));
    double torbit = 1.0/50.0 * 2.*M_PI/OM;
    double xi = 1.0e-6;
    double zeta = 1.0e-5;
    double Y[6];
    get_posvelVector(particle, Y);
    double tdiff = xi * r/( get_vdiffmag(particle, Y) );
    double trho = zeta * r/( get_vrhomag(particle, Y) );

    double dt = torbit;
    if (tstop < dt) {dt=tstop;}
    if (tdiff < dt) {dt=tdiff;}
    if (trho  < dt) {dt=trho; }
    return dt;
}

/**
 * @brief maximum possible timestep, including computational limits
 * 
 * @param particle 
 * @param time 
 * @param tf 
 * @param tout 
 * @return double 
 */
double max_dt(Particle *particle, double time, double tf, double tout) {
    double dt1 = get_dt(particle);
    double dt2 = 0.1 * YR;
    double dt3 = tf-time;
    double dt4 = tout-time;

    double dt;
    dt = dt1;
    if (dt2<dt) {dt=dt2;}
    if (dt3<dt) {dt=dt3;}
    if (dt4<dt) {dt=dt4;}
    double mindt = YR*1.e-12;
    if (dt<mindt) {
        dt = mindt;
        printf("Warning, dt smaller than mindt. setting to mindt\n");
    }
    return dt;

}

void write_output(FILE* file, Particle *part, double time) {
    if (file==NULL) return;
    fprintf(file, "%e\t%e\t%e\t%e\t%e\t%e\t%e\n",
    time,part->x,part->y,part->z,part->vx,part->vy,part->vz);
}

void write_ending(FILE* file, int status) {
    if (file==NULL) return;
    printf("Finished with status: %d\n",status);
}

void write_crossing(char* crossFilename, Particle *part, double time) {
    FILE *crossfile = fopen(crossFilename, "a");
    write_output(crossfile,part,time);
    fclose(crossfile);
}

void heartbeat(double time, double tf) {
    printf("%4.0f %% \r",time/tf*100.);
    fflush(stdout);
}

int integrate(Particle *particle, double t0, double tf, double dtout, int DIFFUSION,
        char* filename, char* resFilename, char* velFilename, char* crossFilename) {
    Model* model = particle->model;
    Domain* domain = model->domain;

    FILE *file;
    int track_trajectory = 0;
    if ( strcmp(filename, "NULL") == 0 ) {
        file = NULL;
    } else {
        printf("Tracking Particle\n");
        track_trajectory=1;
        file = fopen(filename,"w+");
    }

    double *resTimes;
    // top and bottom of disk minus edges
    size_t bigSize = 2*model->nz*model->ny*model->nx;
    int track_resTimes = 0;
    if ( strcmp(resFilename,"NULL") != 0 ) {
        printf("Tracking restime\n");
        track_resTimes = 1;
        resTimes = calloc(bigSize,sizeof(double));
    }

    double *velocities;
    size_t velSize = 2*model->nz*model->ny*model->nx*3;
    int track_velocities = 0;
    if ( strcmp(velFilename,"NULL") != 0 ) {
        printf("Tracking velocities\n");
        track_velocities = 1;
        velocities = calloc(velSize,sizeof(double));
    }

    int track_crossings = 0;
    if ( strcmp(crossFilename,"NULL") != 0 ) {
        printf("Tracking Crossings\n");
        track_crossings = 1;
    }

    double time = t0;
    int status = 0; // status of integration
    double dt = 0.0;
    // output every 1 year
    double next_tout = time + dtout;
    double next_heartbeat = 0.0;
    write_output(file, particle, time);
    while (status == 0) {
        dt = max_dt(particle, time, tf, next_tout);
        rkstep_particle(particle, dt, DIFFUSION);
        time += dt;
        
        double x,y,z;
        x=particle->x;
        y=particle->y;
        z=particle->z;
        double r = sqrt(x*x + y*y + z*z);
        double theta = acos(z/r);
        // check inbounds
        if ((theta < domain->thetaCenters[1]) || (theta > M_PI-domain->thetaCenters[1])) {
            status = OOB;
        }
        if ((r < domain->rCenters[1]) || (r > domain->rCenters[domain->ny-2])) {
            status = OOB;
        }
        if (status == OOB) {
            write_output(file, particle, time);
        }
        // check distance from planet
        double dx = x-model->planetpos[0];
        double dy = y-model->planetpos[1];
        double dz = z-model->planetpos[2];
        double dr = sqrt(dx*dx + dy*dy + dz*dz);
        if (dr < model->planetEnvelope) {
            status = ACCRETED;
        }

        // track restimes or velocities as needed
        if ((track_resTimes) || (track_velocities)) {
            size_t corner[3] = {0,0,0};
            get_edge_from_cart(model,x,y,z,corner);
            if (z<0) {
                corner[2] = 2*model->nz-1-corner[2];
            }
            if (track_resTimes) {
                size_t idx = corner[2]*model->nx*model->ny + corner[1]*model->nx + corner[0];
                resTimes[idx] += dt;
            }
            if (track_velocities) {
                size_t idx = corner[2]*model->nx*model->ny*3 + corner[1]*model->nx*3 + corner[0]*3;
                velocities[idx]   += particle->vx;
                velocities[idx+1] += particle->vy;
                velocities[idx+2] += particle->vz;
            }
        }

        if ((track_crossings) & (r < CROSS_R)) {
            status = CROSSED;
            write_output(file,particle,time);
            write_crossing(crossFilename,particle,time);
        }

        if (time>=tf) {
            // simulation end
            status = COMPLETE;
            write_output(file, particle, time);
        }
        if ((time>=next_tout) && (status == 0)) {
            next_tout += dtout;
            write_output(file, particle, time);
        }
        if (time/tf > next_heartbeat) {
            heartbeat(time,tf);
            next_heartbeat += 0.01;
        }
        
    }
    write_ending(file, status);

    if (track_resTimes) {
        FILE *resFile;
        // if the file exists, read the previous residence times then update
        // our tracker, then re-open the file in write mode to write.
        // format, array[0] = Npart, array[i>0] = restime[i-1]
        resFile = fopen(resFilename,"rb");
        double nparts = 1.0;
        if (resFile != NULL) {
            printf("Reading residence times\n");
            double *old_resTimes = calloc(bigSize,sizeof(double));

            printf("Allocated space for old restime\n");
            fread(&nparts, sizeof(double), 1, resFile);
            printf("Read in number of particles\n");
            fread(old_resTimes, sizeof(double), bigSize, resFile);
            nparts++;
            printf("Add old restimes to new restimes...\n");
            for (size_t i=0; i<bigSize; i++) {
                resTimes[i] += old_resTimes[i];
            }
            free(old_resTimes);
            fclose(resFile);
        }
        resFile = fopen(resFilename,"wb");
        fwrite(&nparts, sizeof(double), 1, resFile);
        fwrite(resTimes, sizeof(double), bigSize, resFile);

        fclose(resFile);
        free(resTimes);
    }

    if (track_velocities) {
        FILE *velFile;
        // if the file exists, read the previous velocities then update
        // our tracker, then re-open the file in write mode to write.
        // format, array[0] = Npart, array[i>0] = velocities[i-1], shape=nz,ny,nx,3
        velFile = fopen(velFilename,"rb");
        double nparts = 1.0;
        if (velFile != NULL) {
            double *old_velocities = calloc(velSize,sizeof(double));
            fread(&nparts, sizeof(double), 1, velFile);
            fread(old_velocities, sizeof(double), velSize, velFile);
            nparts++;
            for (size_t i=0; i<velSize; i++) {
                velocities[i] += old_velocities[i];
            }
            free(old_velocities);
            fclose(velFile);
        }
        velFile = fopen(velFilename,"wb");
        fwrite(&nparts, sizeof(double), 1, velFile);
        fwrite(velocities, sizeof(double), velSize, velFile);

        fclose(velFile);
        free(velocities);
    }

    if (track_trajectory) {
        fclose(file);
    }
    
    return status;
}