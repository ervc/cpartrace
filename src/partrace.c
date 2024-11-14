#include "partrace.h"
// this technically makes this a c++ file
// #include "mlinterp.hpp"

// defaults
#define NX 2048
#define NY 256
#define NZ 32

int main(int argc, char **argv) {
    printf("*** CPARTRACE VERSION %s ***\n",VERSION);
    
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
    
    printf("making Model...\n");
    Model *model = init_Model(inputs->fargodir,inputs->nout,nx,ny,nz);
    printf("Model initialized!\n");

    // seed the random number generator
    srand(time(NULL));

    
    double size0 = inputs->partsize;
    int np = inputs->nparts;
    double sizes[np];
    double xs[np];
    double ys[np];
    double zs[np];
    double rmin = inputs->rmin;
    double rmax = inputs->rmax;
    double phimin = inputs->phimin;
    double phimax = inputs->phimax;
    double thetamin = inputs->thetamin;
    double thetamax = inputs->thetamax;
    for (int i=0; i<np; i++) {
        sizes[i] = size0;
        double phi = random_range(phimin,phimax);
        double r = random_range(rmin,rmax);
        double theta = random_range(thetamin,thetamax);
        xs[i] = r*cos(phi)*sin(theta);
        ys[i] = r*sin(phi)*sin(theta);
        zs[i] = r*cos(theta);
    }

    double t0 = inputs->t0;
    double tf = inputs->tf;
    double dtout = inputs->dtout;
    int final_status=0;
    int all_final[np];

    char resFilename[100];
    if (inputs->residenceTimes) {
        sprintf(resFilename, "%s/residenceTimes.dat",inputs->outputdir);
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
    } else {
        strcpy(crossFilename,"NULL");
    }

    char allpartsFilename[100];
    sprintf(allpartsFilename, "%s/allparts.txt",inputs->outputdir);
    FILE *allpartsf;
    // if the file doesn't exist yet, create it and write the header
    if(!fileExists(allpartsFilename)) {
        allpartsf = fopen(allpartsFilename,"w");
        fprintf(allpartsf,"x0\ty0\tz0\txf\tyf\tzf\n");
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
        Particle *p = init_Particle(model, sizes[i], xs[i], ys[i], zs[i]);
        printf("Integrating...\n");
        final_status = integrate(p, t0, tf, dtout, inputs->diffusion,
                                 filename, resFilename, velFilename, crossFilename);
        // save to the allparts file
        allpartsf = fopen(allpartsFilename,"a");
        fprintf(allpartsf, "%f\t%f\t%f\t%f\t%f\t%f\n",xs[i],ys[i],zs[i],p->x,p->y,p->z);
        fclose(allpartsf);
        all_final[i] = final_status;
        free_Particle(p);
    }

    printf("All statuses: ");
    for (int i=0; i<np; i++) {
        printf("%d, ",all_final[i]);
    }
    printf("\n");

    free_Inputs(inputs);
    free_Model(model);
    return 0;
}

int run_partrace(char *inputfile) {
    char *argv[2];
    strcpy(argv[0],"./partrace");
    strcpy(argv[1],inputfile);
    return main(2,argv);
}