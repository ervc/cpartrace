# include "src/partrace.h"

#define NX 512
#define NY 256
#define NZ 32

double flinterp_mesh(Domain *domain, MeshField *meshdata, int i, size_t j, size_t k,
                        double phi, double r, double theta, int VERBOSE) {
    if (phi<domain->phiCenters[0]) {
        phi += 2*M_PI;
    }
    if (k==domain->nz-1) {
        // midplane interpolation in 2D
        double c00,c10,c01,c11;
        double x0,x1,y0,y1;
        if (i==domain->nx-1) {
        c00 = get_data(meshdata,k  ,j  ,i  );
        c10 = get_data(meshdata,k  ,j  ,0  );
        c01 = get_data(meshdata,k  ,j+1,i  );
        c11 = get_data(meshdata,k  ,j+1,0  );

        x0 = domain->phiCenters[i];
        x1 = domain->phiCenters[0]+(2*M_PI);
        y0 = domain->rCenters[j];
        y1 = domain->rCenters[j+1];
        } else {
        c00 = get_data(meshdata,k  ,j  ,i  );
        c10 = get_data(meshdata,k  ,j  ,i+1);
        c01 = get_data(meshdata,k  ,j+1,i  );
        c11 = get_data(meshdata,k  ,j+1,i+1);

        x0 = domain->phiCenters[i];
        x1 = domain->phiCenters[i+1];
        y0 = domain->rCenters[j];
        y1 = domain->rCenters[j+1];
        }

        double c0 = (x1-phi)/(x1-x0)*c00 + (phi-x0)/(x1-x0)*c10;
        double c1 = (x1-phi)/(x1-x0)*c01 + (phi-x0)/(x1-x0)*c11;

        return (y1-r)/(y1-y0)*c0 + (r-y0)/(y1-y0)*c1;
    } else {
    // this is from wikipedia page on trilinear interpolation
    // values at cube of points around phi,r,theta
    // some special catches for when phi is between nx-1 and 0
    double c000,c100,c010,c110,c001,c101,c011,c111;
    double x0,x1,y0,y1,z0,z1;
    if (i==domain->nx-1) { // phi near pi
    c000 = get_data(meshdata,k  ,j  ,i  );
    c100 = get_data(meshdata,k  ,j  ,0  );
    c010 = get_data(meshdata,k  ,j+1,i  );
    c110 = get_data(meshdata,k  ,j+1,0  );
    c001 = get_data(meshdata,k+1,j  ,i  );
    c101 = get_data(meshdata,k+1,j  ,0  );
    c011 = get_data(meshdata,k+1,j+1,i  );
    c111 = get_data(meshdata,k+1,j+1,0  );

    x0 = domain->phiCenters[i];
    x1 = domain->phiCenters[0]+(2*M_PI);
    y0 = domain->rCenters[j];
    y1 = domain->rCenters[j+1];
    z0 = domain->thetaCenters[k];
    z1 = domain->thetaCenters[k+1];
    if (VERBOSE==1){
        fprintf(stdout, "x0, x1 = %e, %e\n", x0, x1);
        fprintf(stdout, "y0, y1 = %e, %e\n", y0, y1);
        fprintf(stdout, "z0, z1 = %e, %e\n", z0, z1);
    }
    } else {
    c000 = get_data(meshdata,k  ,j  ,i  );
    c100 = get_data(meshdata,k  ,j  ,i+1);
    c010 = get_data(meshdata,k  ,j+1,i  );
    c110 = get_data(meshdata,k  ,j+1,i+1);
    c001 = get_data(meshdata,k+1,j  ,i  );
    c101 = get_data(meshdata,k+1,j  ,i+1);
    c011 = get_data(meshdata,k+1,j+1,i  );
    c111 = get_data(meshdata,k+1,j+1,i+1);
    if (VERBOSE) {
        fprintf(stdout, "c000, c100, c010, c110 = %e, %e, %e, %e\n", c000, c100, c010, c110);
        fprintf(stdout, "c001, c101, c011, c111 = %e, %e, %e, %e\n", c001, c101, c011, c111);
    }

    x0 = domain->phiCenters[i];
    x1 = domain->phiCenters[i+1];
    y0 = domain->rCenters[j];
    y1 = domain->rCenters[j+1];
    z0 = domain->thetaCenters[k];
    z1 = domain->thetaCenters[k+1];
    }

    double xd = (phi-x0)/(x1-x0);
    double yd = (r-y0)/(y1-y0);
    double zd = (theta-z0)/(z1-z0);

    double c00 = c000*(1-xd)+c100*xd;
    double c01 = c001*(1-xd)+c101*xd;
    double c10 = c010*(1-xd)+c110*xd;
    double c11 = c011*(1-xd)+c111*xd;

    double c0 = c00*(1-yd)+c10*yd;
    double c1 = c01*(1-yd)+c11*yd;

    return c0*(1-zd)+c1*zd;
    }
}

int main(int argc, char **argv) {
    if (argc!=2) {
        fputs("PLEASE PROVIDE OUTPUT DIRECTORY\n", stdout);
        return 1;
    }
    // char *outdir = "outputs/Jup_beyondCO_a1E-05/";
    char outdir[50];
    int jsn;
    jsn = snprintf(outdir, sizeof(outdir), argv[1]);
    if (jsn>=sizeof(outdir)) {
        fputs("OUTPUT DIRECTORY TOO LONG!\n", stdout);
        return 1;
    }
    int npart = 1000;

    char infile[50];
    jsn = snprintf(infile, sizeof(infile), "%s/inputs.in", outdir);
    if (jsn>=sizeof(infile)) {
        fputs("INPUT FILE NAME TOO LONG\n", stdout);
    }
    Inputs *inputs = read_inputs(infile);
    Domain *domain = init_Jupiter_Domain(inputs->fargodir, NX, NY, NZ);

    char tempfile[100];
    jsn = snprintf(tempfile, sizeof(tempfile), "%s/dusttemp0.dat", inputs->fargodir);
    if (jsn>=sizeof(tempfile)) {
        fputs("TEMPERATURE FILE NAME TOO LONG\n", stdout);
    }
    MeshField *tempfield = init_MeshField_fromFile(tempfile, NX, NY, NZ, 0);

    double time;
    double x,y,z;
    double vx, vy, vz;
    int status;

    int nt = (int)((inputs->tf - inputs->t0)/inputs->dtout)+1;

    double *dust_temps;
    dust_temps = (double*)malloc(sizeof(double)*npart*nt);

    char partfile[50];
    int nline;
    fprintf(stdout, "Interpolating [%d] particles from output dir: %s\n", npart, outdir);
    for (int n=0; n<npart; n++) {
        jsn = snprintf(partfile, sizeof(partfile), "%s/particle%d.txt", outdir, n);
        if (jsn>=sizeof(partfile)) {
            fputs("PARTICLE FILE NAME TOO LONG\n", stdout);
        }
        FILE *file = fopen(partfile, "r");
        if (file == NULL) {
            printf("Problem opening file: %s", partfile);
            return 1;
        }
        nline = 0;
        while (fscanf(file, "%lf %lf %lf %lf %lf %lf %lf %d", &time, &x, &y, &z, &vx, &vy, &vz, &status) == 8) {
            if (z<0) {
                z=-z;
            }
            double phi = atan2(y, x);
            double r = sqrt(x*x + y*y + z*z);
            double theta = acos(z/r);

            int i, j, k;
            // get the corner
            if (phi<domain->phiCenters[0]) {
                phi += 2*M_PI;
                i = domain->nx-1;
            } else {
                for (i=0; i<domain->nx; i++) {
                    if(domain->phiCenters[i] > phi) {
                        break;
                    }
                }
                i--;
            }
            for (j=0;j<domain->ny;j++) {
                if(domain->rCenters[j] > r) {
                    break;
                }
            }
            j--;
            for (k=0;k<domain->nz;k++) {
                if(domain->thetaCenters[k] > theta) {
                    break;
                }
            }
            k--;

            int index = n*nt + nline;
            if (nline > nt) {
                fprintf(stdout, "TOO MANY LINES!\n");
                return 1;
            }
            int VERBOSE = 0;
            dust_temps[index] = flinterp_mesh(domain, tempfield, i, j, k, phi, r, theta, VERBOSE);
            
            nline++;
        }
        fclose(file);
        fprintf(stdout, "%d/%d\r", n+1, npart);
        fflush(stdout);
    }
    fprintf(stdout, "\n");
    char outfile[50];
    jsn = snprintf(outfile, sizeof(outfile), "%s/ctemp.bdat", outdir);
    if (jsn>sizeof(outfile)) {
        fputs("OUTPUT FILE NAME TOO LONG\n", stdout);
    }
    fprintf(stdout, "Saving to %s ...\n", outfile);
    FILE *file = fopen(outfile, "wb");
    // if (file) {
    //     for (int i=0; i<npart*nt; i++) {
    //         fprintf(file, "%e\n", dust_temps[i]);
    //     }
    // }
    if (file==NULL) {
        fputs("ERROR OPENING OUTPUT FILE\n", stdout);
        return 1;
    }
    /// I should be error checking here but I'm lazy
    fwrite(&npart, sizeof(npart), 1, file);
    fwrite(&nt, sizeof(nt), 1, file);
    if (fwrite(dust_temps, sizeof(double), npart*nt, file) != npart*nt) {
        fputs("Error writing to binary file!\n", stdout);
        return 1;
    }
    fclose(file);

    free(dust_temps);

    fprintf(stdout, "Done!\n");
    return 0;
}