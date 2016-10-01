#include "qksolver_boltzmannbgk_1d2v_shocktube.h"

// STL includes
#include <cmath>
#include <fstream>

#define PI (3.14159265359)

#define EXPORT_MOMENTS
//#define EXPORT_DISTRIBUTION

//#define NO_MOMENTS
//#define NO_BGK_COLLISIONS

#define ASSUME_ALIGNED_DOUBLE(p)  (double *) __builtin_assume_aligned(p,32)
#define ASSUME_ALIGNED_DOUBLE_CONST(p)  (const double *) __builtin_assume_aligned(p,32)

namespace boltzmann_solver
{

    void apply_boundary_conditions(
            QKExtendedDatachunk & n_array
            )
    {
        // Here we do copy out boundary conditions on all boundaries

        // There are six chunks to apply BCs to

        const int Nx = n_array.getLength(0);
        const int Mx = n_array.getLength(1);
        const int Mr = n_array.getLength(2);

        const int xStride = Mx * Mr;
        const int vxStride = Mr;
        const int vrStride = 1;

        const int numGhostLayers = n_array.getInternalRange().getLower(0) - n_array.getLower(0);

        double * __restrict__ n = n_array.getData();

        if(numGhostLayers == 0){
            return;
        }

        // x - copy out
        {
            // Left Wall
            for(int i = 0; i < numGhostLayers; i++){
                for(int j = 0; j < Mx; j++){
                    for(int k = 0; k < Mr; k++){
                        const int toIndex = i * xStride + j * vxStride + k * vrStride;
                        const int fromIndex = (2 * numGhostLayers - 1 - i) * xStride + j * vxStride + k * vrStride;
                        n[toIndex] = n[fromIndex];
                        //printf("xl: Pulling index (%2i, %2i, %2i) and putting in (%2i, %2i, %2i)\n", (2 * numGhostLayers - 1 - i),j,k, i,j,k);
                    }
                }
            }

            // Right Wall
            for(int i = 0; i < numGhostLayers; i++){
                for(int j = 0; j < Mx; j++){
                    for(int k = 0; k < Mr; k++){
                        const int toIndex = (Nx - 1 - i) * xStride + j * vxStride + k * vrStride;
                        const int fromIndex = (Nx - 2 * numGhostLayers + i) * xStride + j * vxStride + k * vrStride;
                        n[toIndex] = n[fromIndex];
                        //printf("xr: Pulling index (%2i, %2i, %2i) and putting in (%2i, %2i, %2i)\n", (Nx - 2 * numGhostLayers + i),j,k, (Nx - 1 - i),j,k);

                    }
                }
            }
        }

        // vx - copy out (helps with integration)
        {
            // Left Wall
            for(int i = 0; i < numGhostLayers; i++){
                for(int j = 0; j < Nx; j++){
                    for(int k = 0; k < Mr; k++){
                        const int toIndex = i * vxStride + j * xStride + k * vrStride;
                        const int fromIndex = (2 * numGhostLayers - 1 - i) * vxStride + j * xStride + k * vrStride;
                        n[toIndex] = n[fromIndex];
                        //printf("vxl: Pulling index (%2i, %2i, %2i) and putting in (%2i, %2i, %2i)\n", j,(2 * numGhostLayers - 1 - i),k, j,i,k);
                    }
                }
            }

            // Right Wall
            for(int i = 0; i < numGhostLayers; i++){
                for(int j = 0; j < Nx; j++){
                    for(int k = 0; k < Mr; k++){
                        const int toIndex = (Mx - 1 - i) * vxStride + j * xStride + k * vrStride;
                        const int fromIndex = (Mx - 2 * numGhostLayers + i) * vxStride + j * xStride + k * vrStride;
                        n[toIndex] = n[fromIndex];
                        //printf("vxr: Pulling index (%2i, %2i, %2i) and putting in (%2i, %2i, %2i)\n", j,(Mx - 2 * numGhostLayers + i),k, j,(Mx - 1 - i),k);
                    }
                }
            }
        }

        // vr - inverse (helps with integration)
        {
            // Left Wall
            for(int i = 0; i < numGhostLayers; i++){
                for(int j = 0; j < Nx; j++){
                    for(int k = 0; k < Mx; k++){
                        const int toIndex = i * vrStride + j * xStride + k * vxStride;
                        const int fromIndex = (2 * numGhostLayers - 1 - i) * vrStride + j * xStride + k * vxStride;
                        n[toIndex] = -n[fromIndex];
                        //printf("vrl: Pulling index (%2i, %2i, %2i) and putting in (%2i, %2i, %2i)\n", j,k,(2 * numGhostLayers - 1 - i), j,k,i);
                    }
                }
            }

            // Right Wall
            for(int i = 0; i < numGhostLayers; i++){
                for(int j = 0; j < Nx; j++){
                    for(int k = 0; k < Mx; k++){
                        const int toIndex = (Mr - 1 - i) * vrStride + j * xStride + k * vxStride;
                        const int fromIndex = (Mr - 2 * numGhostLayers + i) * vrStride + j * xStride + k * vxStride;
                        n[toIndex] = -n[fromIndex];
                        //printf("vrr: Pulling index (%2i, %2i, %2i) and putting in (%2i, %2i, %2i)\n", j,k,(Mr - 2 * numGhostLayers + i), j,k,(Mr - 1 - i));
                    }
                }
            }
        }
    }


}


// 3rd order
    const int num_quad_points = 3;
    const double quad_weights[num_quad_points] = {5./18.,8./18.,5./18.};
    const double quad_points[num_quad_points] = {-0.5*std::sqrt(3./5.),0.,0.5*std::sqrt(3./5.)};

// 4th order
//const int num_quad_points = 4;
//const double quad_weights[num_quad_points] = {
//                      0.5*(18. + std::sqrt(30.))/36.,
//                      0.5*(18. + std::sqrt(30.))/36.,
//                      0.5*(18. - std::sqrt(30.))/36.,
//                      0.5*(18. - std::sqrt(30.))/36.
//                     };
//const double quad_points[num_quad_points] = {
//                      -0.5*std::sqrt(3./7. - 2./7. * std::sqrt(6./5.)),
//                       0.5*std::sqrt(3./7. - 2./7. * std::sqrt(6./5.)),
//                      -0.5*std::sqrt(3./7. + 2./7. * std::sqrt(6./5.)),
//                       0.5*std::sqrt(3./7. + 2./7. * std::sqrt(6./5.)),
//                     };

// 5th order
//    const int num_quad_points = 5;
//    const double quad_weights[num_quad_points] = {0.5*128./255.,
//                          0.5*(322.+13*std::sqrt(70.))/900.,
//                          0.5*(322.+13*std::sqrt(70.))/900.,
//                          0.5*(322.-13*std::sqrt(70.))/900.,
//                          0.5*(322.-13*std::sqrt(70.))/900.
//                         };
//    const double quad_points[num_quad_points] = {0.,
//                          -0.5*std::sqrt(5. - 2. * std::sqrt(10./7.)) / 3.,
//                           0.5*std::sqrt(5. - 2. * std::sqrt(10./7.)) / 3.,
//                          -0.5*std::sqrt(5. + 2. * std::sqrt(10./7.)) / 3.,
//                           0.5*std::sqrt(5. + 2. * std::sqrt(10./7.)) / 3.
//                         };



QKSolver_BoltzmannBGK_1D2V_ShockTube::QKSolver_BoltzmannBGK_1D2V_ShockTube()
{

}

QKSolver_BoltzmannBGK_1D2V_ShockTube::~QKSolver_BoltzmannBGK_1D2V_ShockTube()
{

}

void
QKSolver_BoltzmannBGK_1D2V_ShockTube::setup(const QKRectilinearGrid & grid, const QKSolver_BoltzmannBGK_1D2V_ShockTube_Configuration & config)
{
    _config = config;
    _grid = grid;

    _f.resize(_grid);
    _fs.resize(_grid);
    _fp.resize(_grid);

    {
        int lower=0;
        int upper=_grid.getLength(0);
        QKRange momentsRange(1,&lower,&upper);
        _n.resize(momentsRange);
        _ux.resize(momentsRange);
        _uth.resize(momentsRange);
        _qx.resize(momentsRange);

    }

    // Now we initialize a gaussian in f

    // Get dx for grid
    const double dx[3] = {_grid.getDx(0), _grid.getDx(1), _grid.getDx(2)};

    // Shift xmin to get cell centers
    const double xmin[3] = {_grid.getStart(0) + 0.5*_grid.getDx(0), _grid.getStart(1) + 0.5*_grid.getDx(1), _grid.getStart(2) + 0.5*_grid.getDx(2)};

    // BGK collision terms
    {
        QKRange vxtermsRange = _grid;
        vxtermsRange.set(2,0,1);
        _vxterms.resize(vxtermsRange);

        QKRange vrtermsRange = _grid;
        vrtermsRange.set(0,0,_grid.getLength(0));
        vrtermsRange.set(1,0,1);
        vrtermsRange.set(2,0,num_quad_points*_grid.getLength(2));
        _vrterms.resize(vrtermsRange);


    }

    // Setup integration coefficients
    _i_vx.setup(dx[1],xmin[1],grid.getLength(1));
    _i_vr.setup(dx[2],xmin[2],grid.getLength(2));

    // Iterate through global mesh
    for(QKIndexer chunkIndexer = _f.getIndexer(); chunkIndexer.exists(); chunkIndexer.next()){
        QKExtendedDatachunk & f = _f[chunkIndexer];
        const QKRange & iRange = f.getInternalRange();

        generateInitialConditions(f);
        boltzmann_solver::apply_boundary_conditions(f);
        for(int i = iRange.getLower(0); i < iRange.getUpper(0); i++){
            generateMoments(f,i);
        }
    }


}

void
QKSolver_BoltzmannBGK_1D2V_ShockTube::advance(const double time, const double dt)
{

    for(QKIndexer indexer = _f.indexer(); indexer.exists(); indexer.next()){

        QKExtendedDatachunk & f = _f[indexer];
        QKExtendedDatachunk & fs = _fs[indexer];
        QKExtendedDatachunk & fp = _fp[indexer];

        // Initialize fs with f
        memcpy(fs.data(), f.data(), f.volume() * sizeof(double));

        // Step and apply BCs
        step(dt, f, fs);
        boltzmann_solver::apply_boundary_conditions(fs);

        // Initialize fp with (f + fs)/2
        const double * __restrict__ f_data = f.getData();
        double * __restrict__ fs_data = fs.getData();
        double * __restrict__ fp_data = fp.getData();
        for(int i = 0; i < f.getVolume(); i++){
            fp_data[i] = 0.5*(f_data[i]+fs_data[i]);
        }

        // Step and apply BCs
        step(0.5*dt, fs, fp);
        boltzmann_solver::apply_boundary_conditions(fp);

        // Swap the datasets in preparation for the next time step
        _f[indexer].swap(_fp[indexer]);
    }
}

void
QKSolver_BoltzmannBGK_1D2V_ShockTube::step(
        const double dt,
        const QKExtendedDatachunk & f_edc,
        QKExtendedDatachunk & fp_edc
        )
{

    // The range we wish to iterate over
    const QKRange & iRange = f_edc.getInternalRange();

    // Get dx for grid
    const double dx[3] = {_grid.getDx(0), _grid.getDx(1), _grid.getDx(2)};

    // Shift xmin to get cell centers
    const double xmin[3] = {_grid.getStart(0) + 0.5*dx[0], _grid.getStart(1) + 0.5*dx[1], _grid.getStart(2) + 0.5*dx[2]};

    // Some coefficients
    const double dtodx = dt / dx[0];
    const double cdx2=dx[1]*dx[1]/12.;


    // NOTE : all iteration is in the global ijk coordinates

#pragma omp parallel for
    for(int i = iRange.getLower(0)-1; i < iRange.getUpper(0)+1; i++){
        int index[3] = {i,0,0};
        double x[3] = {0};
        x[0] = xmin[0] + i*dx[0];

        double int_vx0vr0[4] = {0};
        double int_vx1vr0[4] = {0};
        double int_vx2vr0[4] = {0};
        double int_vx3vr0[4] = {0};
        double int_vx0vr2[4] = {0};
        double int_vx1vr2[4] = {0};

        // Calcuate advection terms
        {
            for(int j = iRange.getLower(1); j < iRange.getUpper(1); j++){
                index[1] = j;
                index[2] = iRange.getLower(2);

                x[1] = xmin[1] + j*dx[1];
                const int pos = (x[1]>0) ? 1 : -1;

//                const double C0 = dtodx*fabs(x[1]);


                const double vxsgn = (0. < x[1]) - (x[1] < 0.);

                // We must convert between globak ijk and local ijk
                const double * __restrict__ fc = ASSUME_ALIGNED_DOUBLE_CONST(f_edc.getData(index));
                const double * __restrict__ fl = ASSUME_ALIGNED_DOUBLE_CONST(fc - pos*f_edc.getStride(0));
                const double * __restrict__ fr = ASSUME_ALIGNED_DOUBLE_CONST(fc + pos*f_edc.getStride(0));

                const double * __restrict__ fb = ASSUME_ALIGNED_DOUBLE_CONST(fc-f_edc.getStride(1));
                const double * __restrict__ ft = ASSUME_ALIGNED_DOUBLE_CONST(fc+f_edc.getStride(1));

                const double * __restrict__ fbl = ASSUME_ALIGNED_DOUBLE_CONST(fb-f_edc.getStride(0));
                const double * __restrict__ fbr = ASSUME_ALIGNED_DOUBLE_CONST(fb+f_edc.getStride(0));
                const double * __restrict__ ftl = ASSUME_ALIGNED_DOUBLE_CONST(ft-f_edc.getStride(0));
                const double * __restrict__ ftr = ASSUME_ALIGNED_DOUBLE_CONST(ft+f_edc.getStride(0));

                double * __restrict__ fpc = ASSUME_ALIGNED_DOUBLE(fp_edc.getData(index));
                double * __restrict__ fpr = ASSUME_ALIGNED_DOUBLE(fpc + pos*f_edc.getStride(0));

                int idx = 0;
                for(int k = iRange.getLower(2); k < iRange.getUpper(2); k++){

                    const double fvx_l = (ftl[idx] - fbl[idx]) / (2.*dx[1]);
                    const double fvx_c = (ft [idx] - fb [idx]) / (2.*dx[1]);
                    const double fvx_r = (ftr[idx] - fbr[idx]) / (2.*dx[1]);

                    // Calculate the cell center flux integrated across vx and vr for c,l,r
//                    const double Fl = dx[1]*dx[1]*fvx_l/12. + x[1]*fl[idx];
//                    const double Fc = dx[1]*dx[1]*fvx_c/12. + x[1]*fc[idx];
//                    const double Fr = dx[1]*dx[1]*fvx_r/12. + x[1]*fr[idx];

                    const double Fl = vxsgn*(cdx2*fvx_l + x[1]*fl[idx]);
                    const double Fc = vxsgn*(cdx2*fvx_c + x[1]*fc[idx]);
                    const double Fr = vxsgn*(cdx2*fvx_r + x[1]*fr[idx]);

                    // Gradient ratio
                    //const double r = (fl[idx] != fc[idx]) ? (fr[idx] - fc[idx]) / (fc[idx] - fl[idx]) : 1.0; // 2 add, 1 div : 3 FLO
                    //const double r = (Fl != Fc) ? (Fr - Fc) / (Fc - Fl) : 1.0;
                    const double r = (Fr - Fc) / (Fc - Fl + 1.e-8);

                    // van Albada 1 (leaves spiidxe at front of discontinuity)
                    const double phi = (r*r + r) / (r*r + 1.0);

                    // van Albada 2 (highly diffusive)
                    //const double phi = 2. * r / (r*r+1);

                    // Superbee (sharp for leading wave, but slow)
                    //const double phi = std::max(0.,std::max(std::min(2.*r,1.),std::min(r,2.)));

                    // Flux
//                    const double fluxr = C0 * (fc[idx] + 0.25 * phi * (fr[idx] - fl[idx])); // 2 add, 3 mul : 5 FLO
                    const double fluxr = dtodx*(Fc + 0.25 * phi * (Fr-Fl));

                    // Add flux to downstream element
                    fpr[idx] += fluxr;

                    // Subtract flux from center element
                    fpc[idx] -= fluxr;

#ifndef NO_MOMENTS

                    const double f0 = (44.*fc[idx] + 17.*(fl[idx]+fr[idx]+fb[idx]+ft[idx]) - 10.*(ftl[idx]+ftr[idx]+fbl[idx]+fbr[idx])) / 72.;
                    const double fvx = (fr[idx]+fbr[idx]+ftr[idx] - (fl[idx]+fbl[idx]+ftl[idx])) / (6.*dx[1]);
                    const double fvr = (ft[idx]+ftl[idx]+ftr[idx] - (fb[idx]+fbl[idx]+fbr[idx])) / (6.*dx[2]);
                    const double fvxvx = (fl[idx]+fbl[idx]+ftl[idx] - 2.*(fc[idx]+ft[idx]+fb[idx]) + fr[idx]+fbr[idx]+ftr[idx]) / (6.*dx[1]*dx[1]);
                    const double fvxvr = (ftr[idx]+fbl[idx] - (fbr[idx]+ftl[idx])) / (4.*dx[1]*dx[2]);
                    const double fvrvr = (fb[idx]+fbl[idx]+fbr[idx] - 2.*(fc[idx]+fl[idx]+fr[idx]) + ft[idx]+ftl[idx]+ftr[idx]) / (6.*dx[2]*dx[2]);


#define KAHAN_ADD(a,b) {a[2] = b - a[1]; a[3] = a[0] + a[2]; a[1] = (a[3] - a[0]) - a[2]; a[0] = a[3];}
//#define KAHAN_ADD(a,b) {a[0] += b;}
                    // Add element's contribution to the conserved moments
                    KAHAN_ADD(int_vx0vr0, f0*_i_vx._w0v0   *_i_vr._w0v0     /*+ fvx*_i_vx._w1v0   *_i_vr._w0v0   *//*+ fvr*_i_vx._w0v0   *_i_vr._w1v0   */  + fvxvx*_i_vx._w2v0   *_i_vr._w0v0     /*+ fvxvr*_i_vx._w1v0   *_i_vr._w1v0   */  + fvrvr*_i_vx._w0v0   *_i_vr._w2v0     );
                    KAHAN_ADD(int_vx1vr0, f0*_i_vx._w0v1[j]*_i_vr._w0v0       + fvx*_i_vx._w1v1   *_i_vr._w0v0     /*+ fvr*_i_vx._w0v1[j]*_i_vr._w1v0   */  + fvxvx*_i_vx._w2v1[j]*_i_vr._w0v0     /*+ fvxvr*_i_vx._w1v1   *_i_vr._w1v0   */  + fvrvr*_i_vx._w0v1[j]*_i_vr._w2v0     );
                    KAHAN_ADD(int_vx2vr0, f0*_i_vx._w0v2[j]*_i_vr._w0v0       + fvx*_i_vx._w1v2[j]*_i_vr._w0v0     /*+ fvr*_i_vx._w0v2[j]*_i_vr._w1v0   */  + fvxvx*_i_vx._w2v2[j]*_i_vr._w0v0     /*+ fvxvr*_i_vx._w1v2[j]*_i_vr._w1v0   */  + fvrvr*_i_vx._w0v2[j]*_i_vr._w2v0     );
                    KAHAN_ADD(int_vx3vr0, f0*_i_vx._w0v3[j]*_i_vr._w0v0       + fvx*_i_vx._w1v3[j]*_i_vr._w0v0     /*+ fvr*_i_vx._w0v3[j]*_i_vr._w1v0   */  + fvxvx*_i_vx._w2v3[j]*_i_vr._w0v0     /*+ fvxvr*_i_vx._w1v3[j]*_i_vr._w1v0   */  + fvrvr*_i_vx._w0v3[j]*_i_vr._w2v0     );

                    KAHAN_ADD(int_vx0vr2, f0*_i_vx._w0v0   *_i_vr._w0v2[k]  /*+ fvx*_i_vx._w1v0   *_i_vr._w0v2[k]*/  + fvr*_i_vx._w0v0   *_i_vr._w1v2[k]    + fvxvx*_i_vx._w2v0   *_i_vr._w0v2[k]  /*+ fvxvr*_i_vx._w1v0   *_i_vr._w1v2[k]*/  + fvrvr*_i_vx._w0v0   *_i_vr._w2v2[k]  );
                    KAHAN_ADD(int_vx1vr2, f0*_i_vx._w0v1[j]*_i_vr._w0v2[k]    + fvx*_i_vx._w1v1   *_i_vr._w0v2[k]    + fvr*_i_vx._w0v1[j]*_i_vr._w1v2[k]    + fvxvx*_i_vx._w2v1[j]*_i_vr._w0v2[k]    + fvxvr*_i_vx._w1v1   *_i_vr._w1v2[k]    + fvrvr*_i_vx._w0v1[j]*_i_vr._w2v2[k]  );

#undef KAHAN_ADD

#endif

                    idx++;
                }
            }
        }


#ifndef NO_MOMENTS

        {
            double & n = _n.getData(index)[0];
            double & ux = _ux.getData(index)[0];
            double & uth = _uth.getData(index)[0];
            double & qx = _qx.getData(index)[0];

            // Set n
            n = int_vx0vr0[0];

            // Calculate ux from conserved moments
            ux = int_vx1vr0[0] / n;

            // Calculate pressures
            const double Pxx = int_vx2vr0[0] - n*ux*ux;
            const double Prr = int_vx0vr2[0];

            // Calculate vth from pressures
            uth = std::sqrt((Pxx+Prr) / (3.*n));

            // Caclulate conserved heat flux
            const double Hx = int_vx3vr0[0] + int_vx1vr2[0];

            // Calculate qx from other moments
            qx = 0.5*(Hx - ux*(Pxx+Prr) - 2.*ux*Pxx - n*ux*ux*ux);
        }

        // The moments have been calculated. Now apply them in the collision operator
#ifndef NO_BGK_COLLISIONS
        {

            const double n = _n.getData(index)[0];
            const double ux = _ux.getData(index)[0];
            const double uth = 0.9998*_uth.getData(index)[0];
            //const double uth = 0.9988*_uth.getData(index)[0];
            //const double uth = 0.9997*_uth.getData(index)[0];
            //const double uth = 0.99968*_uth.getData(index)[0];

            const double T = uth*uth;
            const double nu = _config.nu_tau * n / (std::sqrt(2*T) * T);

            const double d = 2. * PI * uth * uth;
            const double C0 = - dt * nu;
            const double C1 = 2. * PI * dt * nu * n / (d * std::sqrt(d));

            index[1]=iRange.getLower(1);
            index[2]=iRange.getLower(2);
            int idx;



            double * __restrict__ vrterms = _vrterms.getData(index);
            idx=0;
            for(int k = iRange.getLower(2); k < iRange.getUpper(2); k++){
                vrterms[idx]=0.;
                for(int n=0;n<num_quad_points;n++){
                    const double wr = (xmin[2] + (k + quad_points[n])*dx[2])/uth;
                    vrterms[idx] += quad_weights[n]*std::exp(-0.5*wr*wr)*wr*uth;
                }
                idx++;
            }

            for(int j = iRange.getLower(1); j < iRange.getUpper(1); j++){

                double fvx=0.;
                for(int n=0;n<num_quad_points;n++){
                    const double wx = (xmin[1] + (j + quad_points[n])*dx[1] - ux)/uth;
                    fvx += quad_weights[n]*std::exp(-0.5*wx*wx);
                }
                fvx*=C1;

                index[1]=j;
                idx=0;
                const double * __restrict__ fc  = f_edc.getData(index);
                double * __restrict__ fpc = fp_edc.getData(index);

                for(int k = iRange.getLower(2); k < iRange.getUpper(2); k++){
                    fpc[idx] += C0*fc[idx] + vrterms[idx]*fvx;
                    idx++;
                }

            }

//            double * __restrict__ vrterms = _vrterms.getData(index);
//            idx=0;
//            for(int k = iRange.getLower(2); k < iRange.getUpper(2); k++){
//                const double vr = xmin[2] + k * dx[2];
//                const double wrc = vr / uth;
//                vrterms[idx++] = std::exp(-0.5*wrc*wrc)*vr;
//            }

//            // Apply the distribution
//            index[2]=iRange.getLower(2);
//            for(int j = iRange.getLower(1); j < iRange.getUpper(1); j++){
//                index[1]=j;
//                const double * __restrict__ fc  = f_edc.getData(index);
//                double * __restrict__ fpc = fp_edc.getData(index);
//                idx = 0;
//                for(int k = iRange.getLower(2); k < iRange.getUpper(2); k++){
//                    fpc[idx] += C0 * fc[idx] + vxterms[j]*vrterms[k];
//                    idx++;
//                }
//            }
        }

#endif // NO_BGK_COLLISIONS

#endif // NO_MOMENTS

        // Mission complete!

    }
}

void
QKSolver_BoltzmannBGK_1D2V_ShockTube::generateInitialConditions(QKExtendedDatachunk & f_edc)
{
    // Range to iterate over
    const QKRange & iRange = f_edc.getInternalRange();

    // Get left and right conditions
    const double uth_l = std::sqrt(_config.T_l / _config.m);
    const double uth_r = std::sqrt(_config.T_r / _config.m);
//    const double c_l = 2. * PI * _config.n_l / std::pow(2. * PI * uth_l * uth_l,1.5);
//    const double c_r = 2. * PI * _config.n_r / std::pow(2. * PI * uth_r * uth_r,1.5);
    //const double C = 2. * PI / std::pow(2. * PI * uth_l * uth_l,1.5);

    // Get grid info
    const double dx[3] = {_grid.getDx(0), _grid.getDx(1), _grid.getDx(2)};
    const double xmin[3] = {_grid.getStart(0) + 0.5*dx[0], _grid.getStart(1) + 0.5*dx[1], _grid.getStart(2) + 0.5*dx[2]};

    // For the quadrature points
    double vx[num_quad_points] = {0};
    double vr[num_quad_points] = {0};
    double fvx[num_quad_points] = {0};
    double fvr[num_quad_points] = {0};

    const double dn = 0.5 * (_config.n_r - _config.n_l);
    const double n0 = 0.5 * (_config.n_r + _config.n_l);
    const double dux = 0.5 * (_config.vx_r - _config.vx_l);
    const double ux0 = 0.5 * (_config.vx_r + _config.vx_l);
    const double duth = 0.5 * (uth_r - uth_l);
    const double uth0 = 0.5 * (uth_r + uth_l);

    // For indexing into f
    int index[3] = {0,0,iRange.getLower(2)};

    for(int i = iRange.getLower(0); i < iRange.getUpper(0); i++){
        index[0] = i;

        const double x  = (xmin[0] + dx[0] * i);

        const double amp = std::tanh(2. * x / _config.transitionWidth);

        const double n = n0 + dn * amp;
        const double uth = uth0 + duth * amp;
        const double ux = ux0 + dux * amp;

        const double c = 2. * PI * n / std::pow(2. * PI * uth * uth,1.5);

//        const double ux  = (x<0.) ? _config.vx_l : _config.vx_r;
//        //const double c   = (x<0.) ? c_l          : c_r;
//        const double c   = C*((_config.n_r - _config.n_l) * amp + _config.n_r + _config.n_l);
//        const double uth = (x<0.) ? uth_l        : uth_r;

        for(int j = iRange.getLower(1); j < iRange.getUpper(1); j++){
            index[1] = j;

            for(int n=0;n<num_quad_points;n++){
                vx[n] = (xmin[1] + (j + quad_points[n])*dx[1] - ux)/uth;
                fvx[n] = c*quad_weights[n]*std::exp(-0.5*vx[n]*vx[n]);
            }

            int idx = 0;
            double * __restrict__ f = f_edc.getData(index);
            for(int k = iRange.getLower(2); k < iRange.getUpper(2); k++){

                for(int n=0;n<num_quad_points;n++){
                    vr[n] = (xmin[2] + (k + quad_points[n])*dx[2])/uth;
                    fvr[n] = quad_weights[n]*std::exp(-0.5*vr[n]*vr[n])*vr[n]*uth;
                }

                f[idx] = 0.;
                for(int n=0; n<num_quad_points;n++){
                    for(int m=0;m<num_quad_points;m++){
                        f[idx] += fvx[n]*fvr[m];
                    }
                }
                idx++;
            }
        }

        generateMoments(f_edc,i);
    }

}

void
QKSolver_BoltzmannBGK_1D2V_ShockTube::generateMoments(const QKExtendedDatachunk & f, const int i)
{

    const QKRange & iRange = f.getInternalRange();
//    const double gradC[3] = {0.5/_grid.getDx(0),0.5/_grid.getDx(1),0.5/_grid.getDx(2)};
//    const double gradC2[3] = {gradC[0]/_grid.getDx(0),gradC[1]/_grid.getDx(1),gradC[2]/_grid.getDx(2)};
    const double dx[3] = {_grid.getDx(0),_grid.getDx(1),_grid.getDx(2)};

    int index[3]= {i,0,iRange.getLower(2)};

    double int_vx0vr0[4] = {0};
    double int_vx1vr0[4] = {0};
    double int_vx2vr0[4] = {0};
    double int_vx3vr0[4] = {0};
    double int_vx0vr2[4] = {0};
    double int_vx1vr2[4] = {0};

    // For each element in vx space
    for(int j = iRange.getLower(1); j < iRange.getUpper(1); j++){
        index[1]=j;
        const double * __restrict__ fc = f.getData(index);
        const double * __restrict__ fl = fc-f.getStride(0);
        const double * __restrict__ fr = fc+f.getStride(0);
        const double * __restrict__ fb = fc-f.getStride(1);
        const double * __restrict__ ft = fc+f.getStride(1);

        const double * __restrict__ fbl = fb-f.getStride(0);
        const double * __restrict__ fbr = fb+f.getStride(0);
        const double * __restrict__ ftl = ft-f.getStride(0);
        const double * __restrict__ ftr = ft+f.getStride(0);

        int idx = 0;
        // For each element in vr space
        for(int k = iRange.getLower(2); k < iRange.getUpper(2); k++){

            // Define cell center f0, and the gradients in vx and vr
//            const double f0  = (26.*fc[idx]-fr[idx]-fl[idx])/24.;
//            const double fvx = gradC[1]*(fr[idx]-fl[idx]);
//            const double fvr = gradC[2]*(ft[idx]-fb[idx]);

//            const double fvxvx=gradC2[1]*(fl[idx]-2.*fc[idx]+fr[idx]);
//            //const double fvxvr;
//            const double fvrvr=gradC2[2]*(fb[idx]-2.*fc[idx]+ft[idx]);

            const double f0 = (44.*fc[idx] + 17.*(fl[idx]+fr[idx]+fb[idx]+ft[idx]) - 10.*(ftl[idx]+ftr[idx]+fbl[idx]+fbr[idx])) / 72.;
            const double fvx = (fr[idx]+fbr[idx]+ftr[idx] - (fl[idx]+fbl[idx]+ftl[idx])) / (6.*dx[1]);
            const double fvr = (ft[idx]+ftl[idx]+ftr[idx] - (fb[idx]+fbl[idx]+fbr[idx])) / (6.*dx[2]);
            const double fvxvx = (fl[idx]+fbl[idx]+ftl[idx] - 2.*(fc[idx]+ft[idx]+fb[idx]) + fr[idx]+fbr[idx]+ftr[idx]) / (6.*dx[1]*dx[1]);
            const double fvxvr = (ftr[idx]+fbl[idx] - (fbr[idx]+ftl[idx])) / (4.*dx[1]*dx[2]);
            const double fvrvr = (fb[idx]+fbl[idx]+fbr[idx] - 2.*(fc[idx]+fl[idx]+fr[idx]) + ft[idx]+ftl[idx]+ftr[idx]) / (6.*dx[2]*dx[2]);



            idx++;

#define KAHAN_ADD(a,b) {a[2] = b - a[1]; a[3] = a[0] + a[2]; a[1] = (a[3] - a[0]) - a[2]; a[0] = a[3];}
//#define KAHAN_ADD(a,b) {a[0] += b;}
            // Add element's contribution to the conserved moments
            KAHAN_ADD(int_vx0vr0, f0*_i_vx._w0v0   *_i_vr._w0v0     /*+ fvx*_i_vx._w1v0   *_i_vr._w0v0   *//*+ fvr*_i_vx._w0v0   *_i_vr._w1v0   */  + fvxvx*_i_vx._w2v0   *_i_vr._w0v0     /*+ fvxvr*_i_vx._w1v0   *_i_vr._w1v0   */  + fvrvr*_i_vx._w0v0   *_i_vr._w2v0     );
            KAHAN_ADD(int_vx1vr0, f0*_i_vx._w0v1[j]*_i_vr._w0v0       + fvx*_i_vx._w1v1   *_i_vr._w0v0     /*+ fvr*_i_vx._w0v1[j]*_i_vr._w1v0   */  + fvxvx*_i_vx._w2v1[j]*_i_vr._w0v0     /*+ fvxvr*_i_vx._w1v1   *_i_vr._w1v0   */  + fvrvr*_i_vx._w0v1[j]*_i_vr._w2v0     );
            KAHAN_ADD(int_vx2vr0, f0*_i_vx._w0v2[j]*_i_vr._w0v0       + fvx*_i_vx._w1v2[j]*_i_vr._w0v0     /*+ fvr*_i_vx._w0v2[j]*_i_vr._w1v0   */  + fvxvx*_i_vx._w2v2[j]*_i_vr._w0v0     /*+ fvxvr*_i_vx._w1v2[j]*_i_vr._w1v0   */  + fvrvr*_i_vx._w0v2[j]*_i_vr._w2v0     );
            KAHAN_ADD(int_vx3vr0, f0*_i_vx._w0v3[j]*_i_vr._w0v0       + fvx*_i_vx._w1v3[j]*_i_vr._w0v0     /*+ fvr*_i_vx._w0v3[j]*_i_vr._w1v0   */  + fvxvx*_i_vx._w2v3[j]*_i_vr._w0v0     /*+ fvxvr*_i_vx._w1v3[j]*_i_vr._w1v0   */  + fvrvr*_i_vx._w0v3[j]*_i_vr._w2v0     );

            KAHAN_ADD(int_vx0vr2, f0*_i_vx._w0v0   *_i_vr._w0v2[k]  /*+ fvx*_i_vx._w1v0   *_i_vr._w0v2[k]*/  + fvr*_i_vx._w0v0   *_i_vr._w1v2[k]    + fvxvx*_i_vx._w2v0   *_i_vr._w0v2[k]  /*+ fvxvr*_i_vx._w1v0   *_i_vr._w1v2[k]*/  + fvrvr*_i_vx._w0v0   *_i_vr._w2v2[k]  );
            KAHAN_ADD(int_vx1vr2, f0*_i_vx._w0v1[j]*_i_vr._w0v2[k]    + fvx*_i_vx._w1v1   *_i_vr._w0v2[k]    + fvr*_i_vx._w0v1[j]*_i_vr._w1v2[k]    + fvxvx*_i_vx._w2v1[j]*_i_vr._w0v2[k]    + fvxvr*_i_vx._w1v1   *_i_vr._w1v2[k]    + fvrvr*_i_vx._w0v1[j]*_i_vr._w2v2[k]  );

//            KAHAN_ADD(n_k,   f0*_wx0vx0   *_wr0vr0      /*+ fvx*_wx1vx0   *_wr0vr0   *//*+ fvr*_wx0vx0   *_wr1vr0   */);
//            KAHAN_ADD(px_k,  f0*_wx0vx1[j]*_wr0vr0        + fvx*_wx1vx1   *_wr0vr0     /*+ fvr*_wx0vx1[j]*_wr1vr0   */);
//            KAHAN_ADD(exx_k, f0*_wx0vx2[j]*_wr0vr0        + fvx*_wx1vx2[j]*_wr0vr0     /*+ fvr*_wx0vx2[j]*_wr1vr0   */);
//            KAHAN_ADD(err_k, f0*_wx0vx0   *_wr0vr2[k]   /*+ fvx*_wx1vx0   *_wr0vr2[k]*/  + fvr*_wx0vx0   *_wr1vr2[k]  );
#undef KAHAN_ADD
        }
    }

    // Reference the proper place in n,vx,vth to place data - offset since moments are in global frame
    double & n = _n.getData(index)[0];
    double & ux = _ux.getData(index)[0];
    double & uth = _uth.getData(index)[0];
    double & qx = _qx.getData(index)[0];

    // Set n
    n = int_vx0vr0[0];

    // Calculate ux from conserved moments
    ux = int_vx1vr0[0] / n;

    // Calculate pressures
    const double Pxx = int_vx2vr0[0] - n*ux*ux;
    const double Prr = int_vx0vr2[0];

    // Calculate vth from pressures
    uth = std::sqrt((Pxx+Prr) / (3.*n));

    // Caclulate conserved heat flux
    const double Hx = int_vx3vr0[0] + int_vx1vr2[0];

    // Calculate qx from other moments
    qx = 0.5*(Hx - ux*(Pxx+Prr) - 2.*ux*Pxx - n*ux*ux*ux);


}


void
QKSolver_BoltzmannBGK_1D2V_ShockTube::writeVTK(const std::string & prefix, const std::string & suffix) const
{

#ifdef EXPORT_DISTRIBUTION
    // Iterate through the domains in n
    for(QKIndexer indexer = _f.getIndexer(); indexer.exists(); indexer.next()){
        const QKExtendedDatachunk & f = _f[indexer];
        const QKRange & iRange = f.getInternalRange();

        // Generate Filename

        // TODO: Add range information
        std::string rangeInfo = "";

        std::string filename = prefix + rangeInfo + suffix;

        std::ofstream file(filename.c_str());
        if(!file){
            std::cout << "QKDatachunk::writeCSV : File " << filename << " could not be created\n";
            exit(EXIT_FAILURE);
        }

        file << "# vtk DataFile Version 2.0\n";
        file << "Stuff for dataset\n";
        file << "ASCII\n";

        file << "DATASET STRUCTURED_POINTS\n";
        file << "DIMENSIONS " << iRange.getLength(0)+1 << " " << iRange.getLength(1)+1 << " " << iRange.getLength(2)+1 << "\n";
        file << "ORIGIN " << _grid.getStart(0) + 0.5*_grid.getWidth(0) << " " << _grid.getStart(1) + 0.5*_grid.getWidth(1) << " " << _grid.getStart(2) + 0.5*_grid.getWidth(2) << "\n";
        file << "SPACING " << _grid.getDx(0) << " " << _grid.getDx(1) << " " << _grid.getDx(2) << "\n\n";

        file << "CELL_DATA " << iRange.getVolume() << std::endl;

        // Write data
        file << "SCALARS f float" << std::endl;
        file << "LOOKUP_TABLE default" << std::endl;

        const double * d = f.getData();
        for(int k = iRange.getLower(2)-f.getLower(2); k < iRange.getUpper(2)-f.getLower(2); k++){
            for(int j = iRange.getLower(1)-f.getLower(1); j < iRange.getUpper(1)-f.getLower(1); j++){
                for(int i = iRange.getLower(0)-f.getLower(0); i < iRange.getUpper(0)-f.getLower(0); i++){
                    file << d[i*f.getStride(0)+j*f.getStride(1)+k*f.getStride(2)] << ' ';
                }
            }
        }
        file << "\n\n";

        file.close();

    }

#endif //EXPORT_DISTRIBUTION

#ifndef NO_MOMENTS

#ifdef EXPORT_MOMENTS

    // Write out the moments

    // Iterate through the domains in n
    {
        // Generate Filename

        // TODO: Add range information
        std::string rangeInfo = "_moments_";

        const QKRange iRange = _n.getInternalRange();

        std::string filename = prefix + rangeInfo + suffix;

        std::ofstream file(filename.c_str());
        if(!file){
            std::cout << "QKDatachunk::writeCSV : File " << filename << " could not be created\n";
            exit(EXIT_FAILURE);
        }

        file << "# vtk DataFile Version 2.0\n";
        file << "Stuff for dataset\n";
        file << "ASCII\n";

        file << "DATASET STRUCTURED_POINTS\n";
        file << "DIMENSIONS " << iRange.getLength(0)+1 << " 1 1\n";
        file << "ORIGIN " << _grid.getStart(0) + 0.5*_grid.getWidth(0) << " 0 0\n";
        file << "SPACING " << _grid.getDx(0) << " " << _grid.getDx(1) << " " << _grid.getDx(2) << "\n\n";

        file << "CELL_DATA " << iRange.getVolume() << std::endl;

        // Write data
        int index[1] = {iRange.getLower(0)};
        int idx;

        const double * n = _n.getData(index);
        const double * ux = _ux.getData(index);
        const double * uth = _uth.getData(index);
        const double * qx = _qx.getData(index);

        file << "SCALARS n float" << std::endl;
        file << "LOOKUP_TABLE default" << std::endl;
        idx = 0;
        for(int i = iRange.getLower(0); i < iRange.getUpper(0); i++){
            file << n[idx++] << ' ';
        }
        file << "\n\n";

        file << "SCALARS u_x float" << std::endl;
        file << "LOOKUP_TABLE default" << std::endl;
        idx = 0;
        for(int i = iRange.getLower(0); i < iRange.getUpper(0); i++){
            file << ux[idx++] << ' ';
        }
        file << "\n\n";

        file << "SCALARS P float" << std::endl;
        file << "LOOKUP_TABLE default" << std::endl;
        idx = 0;
        for(int i = iRange.getLower(0); i < iRange.getUpper(0); i++){
            file << n[idx]*uth[idx]*uth[idx] << ' ';
            idx++;
        }
        file << "\n\n";

        file << "SCALARS q_x float" << std::endl;
        file << "LOOKUP_TABLE default" << std::endl;
        idx = 0;
        for(int i = iRange.getLower(0); i < iRange.getUpper(0); i++){
            file << qx[idx] << ' ';
            idx++;
        }
        file << "\n\n";

        file.close();

    }

#endif //EXPORT_MOMENTS

#endif //NO_MOMENTS

}

void
QKSolver_BoltzmannBGK_1D2V_ShockTube::integral_set_t::setup(double dv, double vmin, double num)
{
    _w0v0 = dv;
    _w0v1.resize(num, 0.);
    _w0v2.resize(num, 0.);
    _w0v3.resize(num, 0.);
    _w1v0 = 0.;
    _w1v1 = dv*dv*dv/12.;
    _w1v2.resize(num, 0.);
    _w1v3.resize(num, 0.);
    _w2v0 = _w1v1;
    _w2v1.resize(num, 0.);
    _w2v2.resize(num, 0.);
    _w2v3.resize(num, 0.);
    const double dv5 = dv*dv*dv*dv*dv;
    for(int i = 0; i < num; i++){
        const double v = vmin + i * dv;
        _w0v1[i] = dv*v;
        _w0v2[i] = v*v*dv+_w1v1;
        _w0v3[i] = v*v*v*dv+3.*v*_w1v1;

        _w1v2[i] = 2.*v*_w1v1;
        _w1v3[i] = 3.*v*v*_w1v1 + dv5/80.;

        _w2v1[i] = v*_w1v1;
        _w2v2[i] = v*v*_w1v1 + dv5/80.;
        _w2v3[i] = v*v*v*_w1v1+3.*v*dv5/80.;
    }

}
