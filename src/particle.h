
// #include <stdlib.h>
// #include "model.h"

#define PARTDENSITY 2.0 // g cm-3

typedef struct Particle {
    Model *model;

    double size;
    double density;

    double x;
    double y;
    double z;

    double vx;
    double vy;
    double vz;
} Particle;

typedef struct ParticleArray {
    size_t np; //number of particles
    Particle *partArray;
} ParticleArray;

void set_initialVelocity(Particle *particle);
void get_posvelVector(Particle *part, double Y[6]);
void print_vector(char* prefix, double Y[6]);

Particle *init_Particle(Model *model, double size, double x, double y, double z) {
    Particle *particle = (Particle*)malloc(sizeof(*particle));

    particle->model = model;
    particle->size = size;
    particle->density = PARTDENSITY;
    particle->x = x;
    particle->y = y;
    particle->z = z;

    set_initialVelocity(particle);
    return particle;
}

void free_Particle(Particle *part) {
    free(part);
}

ParticleArray *init_ParticleArray(Model *model, double np, double *sizes, double *xs, double *ys, double *zs) {
    Particle *particles = (Particle*)malloc(sizeof(*particles) * np);

    for (int i=0; i<np; i++) {
        particles[i].model = model;
        particles[i].size = sizes[i];
        particles[i].density = PARTDENSITY;
        particles[i].x = xs[i];
        particles[i].y = ys[i];
        particles[i].z = zs[i];
        set_initialVelocity(&particles[i]);
    }

    ParticleArray *PA = (ParticleArray*)malloc(sizeof(*PA));
    PA->np = np;
    PA->partArray = particles;

    return PA;
}

void free_ParticleArray(ParticleArray *PA) {
    free(PA);
}

double get_Stokes(Particle *particle) {
    Model *model = particle->model;
    double x = particle->x;
    double y = particle->y;
    double z = particle->z;
    double r,phi,theta;
    phi = atan2(y,x);
    r = sqrt(x*x + y*y + z*z);
    theta = acos(z/r);
    double rho_g = trilinterp_one(model, RHO, phi, r, theta);
    double H = get_scaleheight(r);
    return particle->size*particle->density/H/rho_g;
}

double get_partdiff(Particle *particle, double Y[6]) {
    double x,y,z;
    x=Y[0]; y=Y[1]; z=Y[2];
    double phi,r,theta;
    phi = atan2(y,x);
    r = sqrt(x*x + y*y + z*z);
    theta = acos(z/r);
    double rho_g = trilinterp_one(particle->model, RHO, phi, r, theta);
    double H = get_scaleheight(r);
    double OM = get_Omega(r);
    double St = particle->size*particle->density/H/rho_g;
    double invSt2 = 1/(1+St*St);
    double gasdiff = ALPHA*H*H*OM;
    double partdiff = gasdiff*invSt2;
    return partdiff;

}

void set_initialVelocity(Particle *particle) {
    double x,y,z;
    x = particle->x;
    y = particle->y;
    z = particle->z;
    double r, phi, theta;
    r = sqrt(x*x + y*y + z*z);
    phi = atan2(y,x);
    theta = acos(z/r);

    Model* model = particle->model;
    double vkep = sqrt(G*MSUN/r);
    double St = get_Stokes(particle);
    double Om = get_Omega(r);
    double tstop = St/Om;
    // get non-rotating gas vphi at particle location
    double vphi_gas = trilinterp_one(model, RHO, phi, r, M_PI_2) + r*sin(theta)*OMEGAFRAME;
    double eta = 1.0-pow(vphi_gas/vkep,2.0);
    // from armitage notes eq (140)
    double vr = -eta/(St + 1.0/St) * vkep;
    // armitage notes eq (136)
    // in the rotating frame
    double vphi = vphi_gas - 0.5*tstop*vr*vkep/r - r*sin(theta)*OMEGAFRAME;

    particle->vx = vr*cos(phi)*sin(theta) - vphi*sin(phi)*sin(theta);
    particle->vy = vr*sin(phi)*sin(theta) + vphi*cos(phi)*sin(theta);
    particle->vz = 0.0;
}

void get_posvelVector(Particle *part, double Y[6]) {
    Y[0] = part->x;
    Y[1] = part->y;
    Y[2] = part->z;
    Y[3] = part->vx;
    Y[4] = part->vy;
    Y[5] = part->vz;
}

void update_Particle(Particle *part, double Y[6]) {
    part->x = Y[0];
    part->y = Y[1];
    part->z = Y[2];
    part->vx = Y[3];
    part->vy = Y[4];
    part->vz = Y[5];
}

void print_vector(char* prefix, double Y[6]) {
    printf("%s[",prefix);
    for (int i=0; i<6; i++) {
        printf("%e, ",Y[i]);
    }
    printf("]\n");
}

void print_3vector(char* prefix, double Y[3]) {
    printf("%s[",prefix);
    for (int i=0; i<3; i++) {
        printf("%e, ",Y[i]);
    }
    printf("]\n");
}