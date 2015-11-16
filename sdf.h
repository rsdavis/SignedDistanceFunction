
#include "math.h"

class SDF {

    private:
    
        inline double pos_sq(double x);
        inline double neg_sq(double x);
        inline double max(double a, double b);
        inline int index(int i, int j, int k, int * dims);
        double norm_grad_2d(double * lsf, double * sgn, int i, int j, int * dims);
        double norm_grad_3d(double * lsf, double * sgn, int i, int j, int k, int * dims);

    public:
        
        int construct(double * lsf, int * dims, int ndims, double band, double tolerance);
};

inline double SDF :: pos_sq(double x)
{
    if (x>0) return x*x;
    else return 0;
}

inline double SDF :: neg_sq(double x)
{
    if (x<0) return x*x;
    else return 0;
}

inline double SDF :: max(double a, double b)
{
    if (a >= b) return a;
    else return b;
}

inline int SDF :: index(int i, int j, int k, int * dims)
{
    if (i == dims[0]) i = 0;
    if (i == -1) i = dims[0] - 1;

    if (j == dims[1]) j = 0;
    if (j == -1) j = dims[1] - 1;

    if (k == dims[2]) k = 0;
    if (k == -1) k = dims[2] - 1;

    return i*dims[1]*dims[2] + j*dims[2] + k;
}

double SDF :: norm_grad_3d(double * lsf, double * sgn, int i, int j, int k, int * dims)
{
    unsigned p = index(i,j,k,dims);

    double a = lsf[p] - lsf[index(i-1,j,k,dims)];
    double b = lsf[index(i+1,j,k,dims)] - lsf[p];

    double c = lsf[p] - lsf[index(i,j-1,k,dims)];
    double d = lsf[index(i,j+1,k,dims)] - lsf[p];

    double e = lsf[p] - lsf[index(i,j,k-1,dims)];
    double f = lsf[index(i,j,k+1,dims)] - lsf[p];

    double grad_x, grad_y, grad_z;

    if (sgn[p]>0) {
        grad_x = max(pos_sq(a), neg_sq(b));
        grad_y = max(pos_sq(c), neg_sq(d));
        grad_z = max(pos_sq(e), neg_sq(f));
    } else {
        grad_x = max(pos_sq(b), neg_sq(a));
        grad_y = max(pos_sq(d), neg_sq(c));
        grad_z = max(pos_sq(f), neg_sq(e));
    }

    return sqrt(grad_x + grad_y + grad_z);
}

double SDF :: norm_grad_2d(double * lsf, double * sgn, int i, int j, int * dims)
{
    int k = 0;
    unsigned p = index(i,j,k,dims);

    double a = lsf[p] - lsf[index(i-1,j,k,dims)];
    double b = lsf[index(i+1,j,k,dims)] - lsf[p];

    double c = lsf[p] - lsf[index(i,j-1,k,dims)];
    double d = lsf[index(i,j+1,k,dims)] - lsf[p];

    double grad_x, grad_y;

    if (sgn[p] > 0) {
        grad_x = max( pos_sq(a), neg_sq(b) );
        grad_y = max( pos_sq(c), neg_sq(d) );
    } else {
        grad_x = max( pos_sq(b), neg_sq(a) );
        grad_y = max( pos_sq(d), neg_sq(c) );
    }

    return sqrt(grad_x + grad_y);
}

int SDF :: construct(double * lsf, int * dims, int ndims, double band, double tolerance)
{
    double dt = 0.2;
    double max_dt = tolerance + 1;
    double * sign, * sdf_dt;
    int sdf_dims[3];
    
    sdf_dims[0] = dims[0];
    sdf_dims[1] = dims[1];
    if (ndims == 2) sdf_dims[2] = 1;
    else if (ndims == 3) sdf_dims[2] = dims[2];
    
    sign = new double [dims[0]*dims[1]*dims[2]];
    sdf_dt = new double [dims[0]*dims[1]*dims[2]];

    for (int i=0; i<dims[0]; i++)
    for (int j=0; j<dims[1]; j++)
    for (int k=0; k<dims[2]; k++)
    {
        int ndx = i*dims[1]*dims[2] + j*dims[2] + k;

        if (lsf[ndx] > 0) sign[ndx] = 1;
        else if (lsf[ndx] < 0) sign[ndx] = -1;
        else              sign[ndx] = 0;
    }

    while (max_dt > tolerance)
    {
        max_dt = 0.0;

        for (int i=0; i<dims[0]; i++)
        for (int j=0; j<dims[1]; j++)
        for (int k=0; k<dims[2]; k++)
        {
            int ndx = i*dims[1]*dims[2] + j*dims[2] + k;

            if (fabs(lsf[ndx]) >= band) {
                sdf_dt[ndx] = 0;
            } else {
                double s, grad;

                s = lsf[ndx] / sqrt(lsf[ndx]*lsf[ndx] + 1.0);

                if (ndims == 2) grad = norm_grad_2d(lsf, sign, i, j, sdf_dims);
                else if (ndims == 3) grad = norm_grad_3d(lsf, sign, i, j, k, sdf_dims);
                else return -1;

                sdf_dt[ndx] = s * (1.0 - grad);

                if (fabs(sdf_dt[ndx]) > max_dt) max_dt = fabs(sdf_dt[ndx]);
            }
        }

        for (int i=0; i<dims[0]; i++)
        for (int j=0; j<dims[1]; j++)
        for (int k=0; k<dims[2]; k++)
        {
            int ndx = i*dims[1]*dims[2] + j*dims[2] + k;
            lsf[ndx] += dt * sdf_dt[ndx];
        }
    }

    delete [] sign;
    delete [] sdf_dt;
    return 0;
}
