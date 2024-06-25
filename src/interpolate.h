// #include "newmesh.h"
// #include "domain.h"
// #include "model.h"

double midplane_interp(Model *model, int which, double phi, double r);

/**
 * @brief fast trilinar interpolation assuming you already have corner of cube.
 * does not check for negative z!
 * 
 * @param model 
 * @param which 
 * @param i 
 * @param j 
 * @param k 
 * @param phi 
 * @param r 
 * @param theta 
 * @return double 
 */
double fast_linterp(Model *model, int which, int i, size_t j, size_t k, 
                    double phi, double r, double theta) {
    MeshField *meshdata = get_mesh(model, which);
    Domain *domain = model->domain;
    if (phi<domain->phiCenters[0]) {
        phi += 2*M_PI;
    }
    // this is from wikipedia page on trilinear interpolation
    // values at cube of points around phi,r,theta
    // some special catches for when phi is between nx-1 and 0
    double c000,c100,c010,c110,c001,c101,c011,c111;
    double x0,x1,y0,y1,z0,z1;
    if (i==model->nx-1) { // phi near pi
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
    } else {
    c000 = get_data(meshdata,k  ,j  ,i  );
    c100 = get_data(meshdata,k  ,j  ,i+1);
    c010 = get_data(meshdata,k  ,j+1,i  );
    c110 = get_data(meshdata,k  ,j+1,i+1);
    c001 = get_data(meshdata,k+1,j  ,i  );
    c101 = get_data(meshdata,k+1,j  ,i+1);
    c011 = get_data(meshdata,k+1,j+1,i  );
    c111 = get_data(meshdata,k+1,j+1,i+1);

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

/**
 * @brief Finds the lower corner for a point within the 3d mesh. If i=nx-1, then
 * you are at the edge and need to think about boundary conditions.
 * 
 * @param model 
 * @param phi 
 * @param r 
 * @param theta 
 * @param corner 
 */
void get_corner(Model *model, double phi, double r, double theta, size_t *corner) {
    Domain *domain = model->domain;
    size_t i=0;
    size_t j=0;
    size_t k=0;

    if ((r < domain->rCenters[0]) || (r > domain->rCenters[model->ny-1])) {
        perror("r outside of range!");
        exit(1);
    }
    if ((theta < domain->thetaCenters[0]) || (theta > domain->thetaCenters[model->nz-1])) {
        printf("Theta = %f, theta bounds = (%f, %f)\n",theta,domain->thetaCenters[0],domain->thetaCenters[model->nz-1]);
        perror("theta outside of range!");
        exit(1);
    }

    // if phi is less than phiCenters[0] then it is between phi[-1] and phi[0]
    // this is easier to deal with assuming phi > pi for interpolation purposes.
    if (phi<domain->phiCenters[0]) {
        phi += 2*M_PI;
        i = model->nx-1;
    } else {
    for (i=0;i<model->nx;i++) {
        if(domain->phiCenters[i] > phi) {
            break;
        }
    }
    i--;
    // Note if phi>phiCenters[nx-1] or phi<phiCenters[0] 
    // then we are between end points
    }
    for (j=0;j<model->ny;j++) {
        if(domain->rCenters[j] > r) {
            break;
        }
    }
    j--;
    for (k=0;k<model->nz;k++) {
        if(domain->thetaCenters[k] > theta) {
            break;
        }
    }
    k--;

    corner[0]=i; corner[1]=j; corner[2]=k;
}

/**
 * @brief Returns the lower corner of a cartesian location **As though z were positive**
 * 
 * @param model 
 * @param x 
 * @param y 
 * @param z 
 * @param corner 
 */
void get_corner_from_cart(Model *model, double x, double y, double z, size_t *corner) {
    Domain *domain = model->domain;
    size_t i=0;
    size_t j=0;
    size_t k=0;
    if (z<0) {
        z=-z;
    }
    double phi = atan2(y,x);
    double r = sqrt(x*x + y*y + z*z);
    double theta = acos(z/r);

    if ((r < domain->rCenters[0]) || (r > domain->rCenters[domain->ny-1])) {
        perror("r outside of range!");
        exit(1);
    }
    if (theta < domain->thetaCenters[0]) {
        printf("Theta = %f, theta bounds = (%f, %f)\n",theta,domain->thetaCenters[0],domain->thetaCenters[domain->nz-1]);
        perror("theta outside of range!");
        exit(1);
    }
    if (phi<domain->phiCenters[0]) {
        i = domain->nx-1;
    } else {
        for (i=0;i<domain->nx;i++) {
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
    if (theta > domain->thetaCenters[domain->nz-1]) {
        k = domain->nz-1;
    } else {
        for (k=0;k<domain->nz;k++) {
            if(domain->thetaCenters[k] > theta) {
                break;
            }
        }
        k--;
    }
    corner[0]=i; corner[1]=j; corner[2]=k;

}

/**
 * @brief Returns the lower corner of the cell edge at a cartesian location **As though z were positive**
 * 
 * @param model 
 * @param x 
 * @param y 
 * @param z 
 * @param corner 
 */
void get_edge_from_cart(Model *model, double x, double y, double z, size_t *corner) {
    Domain *domain = model->domain;
    size_t i=0;
    size_t j=0;
    size_t k=0;
    if (z<0) {
        z=-z;
    }
    double phi = atan2(y,x);
    double r = sqrt(x*x + y*y + z*z);
    double theta = acos(z/r);

    if ((r < domain->rEdges[0]) || (r > domain->rEdges[domain->ny])) {
        printf("R = %f, r bounds = (%f, %f)\n",r,domain->rEdges[0],domain->rEdges[domain->ny]);
        perror("r outside of range!");
        exit(1);
    }
    if (theta < domain->thetaEdges[0]) {
        printf("Theta = %f, theta bounds = (%f, %f)\n",theta,domain->thetaEdges[0],domain->thetaEdges[domain->nz-1]);
        perror("theta outside of range!");
        exit(1);
    }
    if (phi<domain->phiEdges[0]) {
        i = domain->nx-1;
    } else {
        for (i=0;i<domain->nx;i++) {
            if(domain->phiEdges[i] > phi) {
                break;
            }
        }
        i--;
    }
    for (j=0;j<domain->ny;j++) {
        if(domain->rEdges[j] > r) {
            break;
        }
    }
    j--;
    if (theta > domain->thetaEdges[domain->nz-1]) {
        k = domain->nz-1;
    } else {
        for (k=0;k<domain->nz;k++) {
            if(domain->thetaEdges[k] > theta) {
                break;
            }
        }
        k--;
    }
    corner[0]=i; corner[1]=j; corner[2]=k;

}

/**
 * @brief trilinear interpolation of data, checks for negative z when called.
 * 
 * @param model 
 * @param which 
 * @param phi 
 * @param r 
 * @param theta 
 * @return double 
 */
double trilinterp_one(Model *model, int which, double phi, double r, double theta) {
    size_t i=0;
    size_t j=0;
    size_t k=0;

    double negz = 1.0;
    // this is a negative z
    if (theta > M_PI_2) {
        theta = M_PI-theta;
        negz = -1.0;
    }
    double result;
    // this is near the midplane (we've already accounted for negative z)
    if (theta > model->domain->thetaCenters[model->nz-1]) {
        result = midplane_interp(model, which, phi, r);
    } else {
        size_t corner[3];
        get_corner(model,phi,r,theta,corner);
        i=corner[0]; j=corner[1]; k=corner[2];
        result = fast_linterp(model, which, i,j,k, phi,r,theta);
    }
    switch (which)
    {
        // cases for which the value should be negative if z is negative.
        case VTHETA:
        case VZ:
        case DRHODZ:
            result = negz*result;
            break;
        
        default:
            break;
    }
    return result;
}

double fast_midplane_interp(Model *model, int which, size_t i, size_t j, 
                            double phi, double r) {
    MeshField *meshdata = get_mesh(model, which);
    Domain *domain = model->domain;
    // this is from wikipedia page on bilinear interpolation
    // values at square of points around phi,r
    size_t k = model->nz-1;
    double c00,c10,c01,c11;
    double x0,x1,y0,y1;
    if (i==model->nx-1) {
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
}

/**
 * @brief Interpolation for when the particle is at (or near) the midplane
 * 
 * @param model 
 * @param which 
 * @param phi 
 * @param r 
 * @return double 
 */
double midplane_interp(Model *model, int which, double phi, double r) {
    size_t i=0;
    size_t j=0;
    size_t k=model->nz-1;
    Domain *domain = model->domain;

    if (phi<domain->phiCenters[0]) {
        phi += 2*M_PI;
        i = model->nx-1;
    } else if (phi>domain->phiCenters[model->nx-1]) {
        i = model->nx-1;
    } else {
        for (i=0;i<model->nx;i++) {
            if(domain->phiCenters[i] > phi) {
                break;
            }
        }
        i--;
    }
    for (j=0;j<model->ny;j++) {
        if(domain->rCenters[j] > r) {
            break;
        }
    }
    j--;

    MeshField *meshdata = get_mesh(model, which);
    // this is from wikipedia page on bilinear interpolation
    // values at square of points around phi,r

    double c00,c10,c01,c11;
    double x0,x1,y0,y1;
    if (i==model->nx-1) {
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
}

double cart_trilinterp_one(Model *model, int which, double x, double y, double z) {
    double phi,r,theta;
    phi = atan2(y,x);
    r = sqrt(x*x + y*y +z*z);
    theta = acos(z/r);
    return trilinterp_one(model, which, phi, r, theta);
}


void cart_trilinterp(Model *model, int which, int ni, double *results,
                    double *xi, double *yi, double *zi) {
    double x,y,z;
    double phi,r,theta;
    for (int i=0; i<ni; i++) {
        x = xi[i];
        y = yi[i];
        z = zi[i];
        phi = atan2(y,x);
        r = sqrt(x*x + y*y + z*z);
        theta = acos(z/r);
        results[i] = trilinterp_one(model,which,phi,r,theta);
    }
}

void trilinterp(Model *model, int which, int ni, double *results, 
                double *phii, double *ri, double *thetai) {
    for (int i=0; i<ni; i++) {
        results[i] = trilinterp_one(model,which,phii[i],ri[i],thetai[i]);
    }    
}

