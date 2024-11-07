#include "partrace.h"
// this technically makes this a c++ file
// #include "mlinterp.hpp"

// defaults
#define NX 680
#define NY 215
#define NZ 20
#define NLVL 5

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
    double rmin = 4.2*AU;
    double rmax = 6.2*AU;
    double phimin = -1.*M_PI/4.;
    double phimax =  1.*M_PI/4.;
    double zmin = -1.0; //scaleheight
    double zmax =  1.0; //scaleheight
    for (int i=0; i<np; i++) {
        sizes[i] = size0; // /( (double)pow(10.0,i) );
        double phi = random_range(phimin,phimax);
        double r = random_range(rmin,rmax);
        // if (r > 4.2*AU) {r+= 2.0*AU;}
        // double r = rmin+ i*(rmax-rmin)/np;
        double z = r*0.05*random_range(zmin,zmax);
        xs[i] = r*cos(phi);
        ys[i] = r*sin(phi);
        zs[i] = z;
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
    allpartsf = fopen(allpartsFilename,"w+");
    fprintf(allpartsf,"size\tx0\ty0\tz0\txf\tyf\tzf\n");
    fclose(allpartsf);

    // TODO: Parallelize this loop
    for (int i=0; i<np; i++) {
        printf("Starting loop\n");
        char filename[100];
        // save every 10th output
        if ((i%1) == 0) {
            sprintf(filename, "%s/particle%d.txt",inputs->outputdir,i);
        } else {
            strcpy(filename,"NULL");
        }
        printf("Starting number: %d\n",i);
        if (strcmp(filename,"NULL") != 0) {
            printf("Saving output to %s\n",filename);
        }
        Particle *p = init_Particle(models, nlvl, sizes[i], xs[i], ys[i], zs[i]);
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
    free_Models(models, nlvl);
    return 0;
}

int run_partrace(char *inputfile) {
    char *argv[2];
    strcpy(argv[0],"./partrace");
    strcpy(argv[1],inputfile);
    return main(2,argv);
}