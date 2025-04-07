#include "partrace.h"
// this technically makes this a c++ file
// #include "mlinterp.hpp"

// defaults
#define NX 512
#define NY 256
#define NZ 32
#define NLVL 5

void init_random_particles(Inputs *inputs, double *sizes, double *xs, double *ys, double *zs);
void read_partifle(Inputs *inputs, double *sizes, double *xs, double *ys, double *zs);

int main(int argc, char **argv) {
    printf("*** CPARTRACE VERSION %s ***\n",VERSION);
    
    // read inputs
    char infile[100];
    Inputs *inputs = init_Inputs();
    if (argc <=1 ) {
        printf("No input supplied, using defaults\n");
    } else {
        strcpy(infile, argv[1]);
        printf("Reading input file: %s\n", infile);
        inputs = read_inputs(infile);
    }
    if ( makedir(inputs->outputdir) < 0 ) { exit(1); }
    FILE *fin;
    char inputout[100];
    sprintf(inputout,"%s/inputs.in",inputs->outputdir);
    fin = fopen(inputout,"w");
    if (fin==NULL) {
        fprintf(stderr,"Cannot create input file in output directory!");
        exit(1);
    }
    fprintf_Inputs(fin,inputs);
    fclose(fin);
    
    const size_t nx = NX;
    const size_t ny = NY;
    const size_t nz = NZ;
    
    // make the model
    printf("making Model...\n");
    Model *models[NLVL];
    int nlvl = NLVL;
    if (inputs->modeltype==JUPITER_MODEL) {
        nlvl = NLVL;
        size_t nxs[NLVL] = {680, 120, 120, 120, 120};
        size_t nys[NLVL] = {215, 120, 120, 120, 120};
        size_t nzs[NLVL] = {20, 34, 62, 86, 86};
        for (int i=0; i<NLVL; i++) {
            char leveldir[100];
            sprintf(leveldir,"%s/fargolev%d/",inputs->fargodir,i);
            models[i] = init_Model(inputs->modeltype,leveldir,"0",nxs[i],nys[i],nzs[i]);
        }
        // default model level0
        printf("Models initialized\n");
    } else {
        nlvl = 1;
        models[0] = init_Model(inputs->modeltype,inputs->fargodir,inputs->nout,nx,ny,nz);
        printf("Model initialized!\n");
    }

    // seed the random number generator
    srand(time(NULL));

    
    double size0 = inputs->partsize;
    int np = inputs->nparts;
    double sizes[np];
    double xs[np];
    double ys[np];
    double zs[np];
    if ( strcmp(inputs->partfile, "NULL") == 0 ) {
        init_random_particles(inputs, sizes, xs, ys, zs);
    } else {
        read_partfile(inputs, sizes, xs, ys, zs);
    }

    double t0 = inputs->t0;
    double tf = inputs->tf;
    double dtout = inputs->dtout;
    Intout result;
    result.status = 0;
    result.tf = 0.0;
    int all_final[np];

    char resFilename[100];
    if (inputs->residenceTimes) {
        if (nlvl > 1) {
            printf("Cannot currently track residence times with multilevel model\n");
            return 1;
        }
        Model* model = models[0];
        sprintf(resFilename, "%s/residenceTimes.dat",inputs->outputdir);
        if (inputs->reset) {
            printf("!!! Resetting Residence Times !!!\n");
            FILE *resFile;
            resFile = fopen(resFilename,"wb");
            size_t bigSize = 2*model->nz*model->ny*model->nx;
            double zero = 0;
            // bigSize+1 because we also include the number of particles as the first double
            for (int i=0; i<bigSize+1; i++) {
                fwrite(&zero, sizeof(double), 1, resFile);
            }
            fclose(resFile);
        }
    } else {
        strcpy(resFilename,"NULL");
    }

    char velFilename[100];
    if (inputs->velocities) {
        sprintf(velFilename, "%s/velocities.dat",inputs->outputdir);
    } else {
        strcpy(velFilename,"NULL");
    }

    char crossFilename[100];
    if (inputs->crossings) {
        sprintf(crossFilename, "%s/partCrossings.txt",inputs->outputdir);
        // if reset option then write a blank file
        if (inputs->reset) {
            FILE *crossFile;
            crossFile = fopen(crossFilename,"w");
            fclose(crossFile);
        }
    } else {
        strcpy(crossFilename,"NULL");
    }

    char allpartsFilename[100];
    sprintf(allpartsFilename, "%s/allparts.txt",inputs->outputdir);
    FILE *allpartsf;
    // if the file doesn't exist yet or reset is picked, create it and write the header
    if(!fileExists(allpartsFilename) || inputs->reset) {
        allpartsf = fopen(allpartsFilename,"w");
        fprintf(allpartsf,"tf\tx0\ty0\tz0\txf\tyf\tzf\tstatus\n");
        fclose(allpartsf);
    }

    // TODO: Parallelize this loop
    for (int i=0; i<np; i++) {
        printf("Starting loop\n");
        char filename[100];
        // save every dsave-th output
        if ((i%inputs->dsave) == 0) {
            sprintf(filename, "%s/particle%d.txt",inputs->outputdir,i+inputs->nstart);
        } else {
            strcpy(filename,"NULL");
        }
        printf("Starting number: %d\n",i);
        if (strcmp(filename,"NULL") != 0) {
            printf("Saving output to %s\n",filename);
        }
        Particle *p = init_Particle(models, nlvl, sizes[i], xs[i], ys[i], zs[i]);
        printf("Integrating...\n");
        result = integrate(p, t0, tf, dtout, inputs->diffusion,
                                 filename, resFilename, velFilename, crossFilename);
        // save to the allparts file
        allpartsf = fopen(allpartsFilename,"a");
        fprintf(allpartsf, "%f\t%f\t%f\t%f\t%f\t%f\t%f\t%d\n",result.tf,xs[i],ys[i],zs[i],p->x,p->y,p->z,result.status);
        fclose(allpartsf);
        all_final[i] = result.status;
        free_Particle(p);
    }

    printf("All statuses: ");
    for (int i=0; i<np; i++) {
        printf("%d, ",all_final[i]);
    }
    printf("\n");

    free_Inputs(inputs);
    free_Models(models, nlvl);
    return 0;
}

void init_random_particles(Inputs *inputs, double *sizes, double *xs, double *ys, double *zs) {
    double rmin = inputs->rmin;
    double rmax = inputs->rmax;
    double phimin = inputs->phimin;
    double phimax = inputs->phimax;
    double thetamin = inputs->thetamin;
    double thetamax = inputs->thetamax;
    for (int i=0; i<inputs->nparts; i++) {
        sizes[i] = inputs->partsize;
        double phi = random_range(phimin,phimax);
        double r = random_range(rmin,rmax);
        double theta = random_range(thetamin,thetamax);
        xs[i] = r*cos(phi)*sin(theta);
        ys[i] = r*sin(phi)*sin(theta);
        zs[i] = r*cos(theta);
    }
}

void read_partfile(Inputs *inputs, double *sizes, double *xs, double *ys, double *zs) {
    FILE *file;
    file = fopen(inputs->partfile, "r");
    if (file==NULL) {
        exit(EXIT_FAILURE);
    }
    double s, x, y, z;
    int nline = 0;
    while ( fscanf(file, "%lf %lf %lf %lf", &s, &x, &y, &z) == 4 ) {
        if ( nline>inputs->nparts ) {
            fprintf(stdout, "TOO MANY LINES IN PARTFILE!\n");
            exit(EXIT_FAILURE);
        }
        sizes[nline] = s;
        xs[nline] = x;
        ys[nline] = y;
        zs[nline] = z;
        nline++;
    }
}

int run_partrace(char *inputfile) {
    char *argv[2];
    strcpy(argv[0],"./partrace");
    strcpy(argv[1],inputfile);
    return main(2,argv);
}