#include "qksolver_boltzmannbgk_1d2v_shocktube_var.h"

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

namespace boltzmann_solver_var
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



QKSolver_BoltzmannBGK_1D2V_ShockTube_Var::QKSolver_BoltzmannBGK_1D2V_ShockTube_Var():
    _offsetIsSet(false),
    _offset_density(1.),
    _offset_pressurexx(1.),
    _offset_pressurerr(1.)
{

}

QKSolver_BoltzmannBGK_1D2V_ShockTube_Var::~QKSolver_BoltzmannBGK_1D2V_ShockTube_Var()
{

}

void
QKSolver_BoltzmannBGK_1D2V_ShockTube_Var::setup(const QKStructuredGrid & grid, const QKSolver_BoltzmannBGK_1D2V_ShockTube_Configuration & config)
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
//    _i_vx.setup(dx[1],xmin[1],grid.getLength(1));
//    _i_vr.setup(dx[2],xmin[2],grid.getLength(2));

    // Now we initialize a gaussian in f

    // Iterate through global mesh
    for(QKIndexer chunkIndexer = _f.getIndexer(); chunkIndexer.exists(); chunkIndexer.next()){
        QKExtendedDatachunk & f = _f[chunkIndexer];
        const QKRange & iRange = f.getInternalRange();

        generateInitialConditions(f);
        boltzmann_solver_var::apply_boundary_conditions(f);
#ifndef NO_MOMENTS
        for(int i = iRange.getLower(0); i < iRange.getUpper(0); i++){
            generateMoments(f,i);
        }
#endif
    }


}

void
QKSolver_BoltzmannBGK_1D2V_ShockTube_Var::advance(const double time, const double dt)
{

    for(QKIndexer indexer = _f.getIndexer(); indexer.exists(); indexer.next()){

        QKExtendedDatachunk & f = _f[indexer];
        QKExtendedDatachunk & fs = _fs[indexer];
        QKExtendedDatachunk & fp = _fp[indexer];

        // Initialize fs with f
        memcpy(fs.getData(), f.getData(), f.getVolume() * sizeof(double));

        // Step and apply BCs
        step(dt, f, fs);
        boltzmann_solver_var::apply_boundary_conditions(fs);

        // Initialize fp with (f + fs)/2
        const double * __restrict__ f_data = f.getData();
        double * __restrict__ fs_data = fs.getData();
        double * __restrict__ fp_data = fp.getData();
        for(int i = 0; i < f.getVolume(); i++){
            fp_data[i] = 0.5*(f_data[i]+fs_data[i]);
        }

        // Step and apply BCs
        step(0.5*dt, fs, fp);
        boltzmann_solver_var::apply_boundary_conditions(fp);

        // Swap the datasets in preparation for the next time step
        _f[indexer].swap(_fp[indexer]);
    }
}

void
QKSolver_BoltzmannBGK_1D2V_ShockTube_Var::step(
        const double dt,
        const QKExtendedDatachunk & f_edc,
        QKExtendedDatachunk & fp_edc
        )
{

    // The range we wish to iterate over
    const QKRange & iRange = f_edc.getInternalRange();

    // Assume dx is constant

    // NOTE : all iteration is in the global ijk coordinates

#pragma omp parallel for
    for(int i = iRange.getLower(0)-1; i < iRange.getUpper(0)+1; i++){
        int index[3] = {i,0,0};
        double dx = 0.;

        // Grid doesn't exist in ghost cells so we have to fake it by reflecting the grid across the wall
        if(i < iRange.getLower(0)){
            dx = _grid.getDx(0,iRange.getLower(0));
        } else if(i >= iRange.getUpper(0)){
            dx = _grid.getDx(0,iRange.getUpper(0)-1);
        } else {
            dx = _grid.getDx(0,i);
        }

        // Some coefficients
        const double dtodx = dt / dx;

        // Calcuate advection terms
        {
            for(int j = iRange.getLower(1); j < iRange.getUpper(1); j++){
                index[1] = j;
                index[2] = iRange.getLower(2);

                const double dvx = _grid.getDx(1,j);
                const double vxc = _grid.getCentroid(1,j);

                const int pos = (vxc>0) ? 1 : -1;
                const double absvx = std::fabs(vxc);

                // We must convert between globak ijk and local ijk
                const double * __restrict__ fc = ASSUME_ALIGNED_DOUBLE_CONST(f_edc.getData(index));
                const double * __restrict__ fl = ASSUME_ALIGNED_DOUBLE_CONST(fc - pos*f_edc.getStride(0));
                const double * __restrict__ fr = ASSUME_ALIGNED_DOUBLE_CONST(fc + pos*f_edc.getStride(0));

                double * __restrict__ fpc = ASSUME_ALIGNED_DOUBLE(fp_edc.getData(index));
                double * __restrict__ fpr = ASSUME_ALIGNED_DOUBLE(fpc + pos*f_edc.getStride(0));

                int idx = 0;
                for(int k = iRange.getLower(2); k < iRange.getUpper(2); k++){

                    const double Fl = absvx*fl[idx];
                    const double Fc = absvx*fc[idx];
                    const double Fr = absvx*fr[idx];

                    // Gradient ratio
                    //const double r = (fl[idx] != fc[idx]) ? (fr[idx] - fc[idx]) / (fc[idx] - fl[idx]) : 1.0; // 2 add, 1 div : 3 FLO
                    const double r = (Fl != Fc) ? (Fr - Fc) / (Fc - Fl) : 1.0;
                    //const double r = (Fr - Fc) / (Fc - Fl + 1.e-8);

                    // van Albada 1 (leaves spiidxe at front of discontinuity)
                    //const double phi = (r*r + r) / (r*r + 1.0);

                    // van Albada 2 (highly diffusive)
                    const double phi = 2. * r / (r*r+1);

                    // Superbee (sharp for leading wave, but slow)
                    //const double phi = std::max(0.,std::max(std::min(2.*r,1.),std::min(r,2.)));

                    // Flux
                    const double fluxr = dtodx*(Fc + 0.25 * phi * (Fr-Fl));

                    // Add flux to downstream element
                    fpr[idx] += fluxr;

                    // Subtract flux from center element
                    fpc[idx] -= fluxr;

                    idx++;
                }
            }
        }


#ifndef NO_MOMENTS
        const bool in_domain = i >= iRange.getLower(0) && i < iRange.getUpper(0);
        if(in_domain){
            generateMoments(f_edc,i);
        }

        // The moments have been calculated. Now apply them in the collision operator
#ifndef NO_BGK_COLLISIONS
        {

            const double n = _n.getData(index)[0];
            const double ux = _ux.getData(index)[0];
            const double uth = _uth.getData(index)[0];
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
                const double vrc = _grid.getCentroid(2,k);
                const double dvr = _grid.getDx(2,k);
                for(int n=0;n<num_quad_points;n++){
                    const double wr = (vrc + quad_points[n]*dvr)/uth;
                    vrterms[idx] += quad_weights[n]*std::exp(-0.5*wr*wr)*wr*uth;
                }
                idx++;
            }

            for(int j = iRange.getLower(1); j < iRange.getUpper(1); j++){

                double fvx=0.;
                const double vxc = _grid.getCentroid(1,j);
                const double dvx = _grid.getDx(1,j);
                for(int n=0;n<num_quad_points;n++){
                    const double wx = (vxc + quad_points[n]*dvx - ux)/uth;
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

        }

#endif // NO_BGK_COLLISIONS

#endif // NO_MOMENTS

        // Mission complete!

    }
}

void
QKSolver_BoltzmannBGK_1D2V_ShockTube_Var::generateInitialConditions(QKExtendedDatachunk & f_edc)
{
    // Range to iterate over
    const QKRange & iRange = f_edc.getInternalRange();

    // Get left and right conditions
    const double uth_l = std::sqrt(_config.T_l / _config.m);
    const double uth_r = std::sqrt(_config.T_r / _config.m);

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

        const double x = _grid.getCentroid(0,i);

        const double amp = std::tanh(2. * x / _config.transitionWidth);

        const double n = n0 + dn * amp;
        const double uth = uth0 + duth * amp;
        const double ux = ux0 + dux * amp;

        const double c = 2. * PI * n / std::pow(2. * PI * uth * uth,1.5);

        for(int j = iRange.getLower(1); j < iRange.getUpper(1); j++){
            index[1] = j;

            const double dvx = _grid.getDx(1,j);
            const double vxc = _grid.getCentroid(1,j);
            double fvx = 0.;
            for(int n=0;n<num_quad_points;n++){
                const double vx = (vxc - ux + quad_points[n]*dvx)/uth;
                fvx += c*quad_weights[n]*std::exp(-0.5*vx*vx);
            }

            int idx = 0;
            double * __restrict__ f = f_edc.getData(index);
            for(int k = iRange.getLower(2); k < iRange.getUpper(2); k++){
                const double dvr = _grid.getDx(2,k);
                const double vrc = _grid.getCentroid(2,k);
                double fvr = 0.;
                for(int n=0;n<num_quad_points;n++){
                    const double vr = (vrc + quad_points[n]*dvr)/uth;
                    fvr += quad_weights[n]*std::exp(-0.5*vr*vr)*vr*uth;
                }

                f[idx] = fvr*fvx;
                idx++;
            }
        }
#ifndef NO_MOMENTS
        generateMoments(f_edc,i);
#endif
    }

}

void
QKSolver_BoltzmannBGK_1D2V_ShockTube_Var::generateMoments(const QKExtendedDatachunk & f, const int i)
{

    const QKRange & iRange = f.getInternalRange();
//    const double gradC[3] = {0.5/_grid.getDx(0),0.5/_grid.getDx(1),0.5/_grid.getDx(2)};
//    const double gradC2[3] = {gradC[0]/_grid.getDx(0),gradC[1]/_grid.getDx(1),gradC[2]/_grid.getDx(2)};
    double dx[3] = {0};
    double x[3] = {0};

    x[0] = _grid.getCentroid(0,i);
    dx[0] = _grid.getDx(0,i);

    int index[3]= {i,0,iRange.getLower(2)};

    double int_vx0vr0[4] = {0};
    double int_vx1vr0[4] = {0};
    double int_vx2vr0[4] = {0};
    double int_vx3vr0[4] = {0};
    double int_vx0vr2[4] = {0};
    double int_vx1vr2[4] = {0};

    double c_vx0vr0[3] = {1.,0.,0.};
    double c_vx1vr0[3] = {0};
    double c_vx2vr0[3] = {0};
    double c_vx3vr0[3] = {0};
    double c_vx0vr2[3] = {0};
    double c_vx1vr2[3] = {0};

    double c_vx[3] = {0};
    double c_vr[3] = {0};

    // For each element in vx space
    for(int j = iRange.getLower(1); j < iRange.getUpper(1); j++){
        index[1]=j;
        x[1] = _grid.getCentroid(1,j);
        dx[1] = _grid.getDx(1,j);

        c_vx1vr0[0] = x[1];
        c_vx1vr0[1] = dx[1]*dx[1]/12.;

        c_vx2vr0[0] = dx[1]*dx[1]/12. + x[1]*x[1];
        c_vx2vr0[1] = dx[1]*dx[1]*x[1]/6.;

        c_vx3vr0[0] = x[1]*(dx[1]*dx[1]/4.+x[1]*x[1]);
        c_vx3vr0[1] = dx[1]*dx[1]*(dx[1]*dx[1]+20.*x[1]*x[1])/80.;

        const int jp = (j==iRange.getUpper(1)-1)? j : j+1;
        const int jn = (j==iRange.getLower(1))? j : j-1;

        c_vx[0] = -1./(dx[1]+_grid.getDx(1,jn));
        c_vx[2] =  1./(dx[1]+_grid.getDx(1,jp));
        c_vx[1] = -1./(c_vx[0] + c_vx[2]);

        const double * __restrict__ fc = f.getData(index);
        const double * __restrict__ fl = fc-f.getStride(0);
        const double * __restrict__ fr = fc+f.getStride(0);
        const double * __restrict__ fb = fc-f.getStride(1);
        const double * __restrict__ ft = fc+f.getStride(1);

        int idx = 0;
        // For each element in vr space
        for(int k = iRange.getLower(2); k < iRange.getUpper(2); k++){

            x[2] = _grid.getCentroid(2,k);
            dx[2] = _grid.getDx(2,k);
            const double dV = dx[2]*dx[1];

            c_vx0vr2[0] = dx[2]*dx[2]/12. + x[2]*x[2];
            c_vx0vr2[2] = dx[2]*dx[2]*x[2]/6.;

            c_vx1vr2[0] = x[1]*(dx[2]*dx[2]/12. + x[2]*x[2]);
            c_vx1vr2[1] = dx[1]*dx[1]*(dx[2]*dx[2] + 12.*x[2]*x[2])/144.;
            c_vx1vr2[2] = dx[2]*dx[2]*x[2]*x[1]/6.;

            const int kp = (k==iRange.getUpper(2)-1)? k : k+1;
            const int kn = (k==iRange.getLower(2))? k : k-1;

            c_vr[0] = -1./(dx[2]+_grid.getDx(2,kn));
            c_vr[2] =  1./(dx[2]+_grid.getDx(2,kp));
            c_vr[1] = -1./(c_vr[0] + c_vr[2]);

            // Define cell center f0, and the gradients in vx and vr
            const double f0 = fc[idx];
            const double fvx = 0.;//c_vx[0]*fl[idx]+c_vx[1]*fc[idx]+c_vx[2]*fr[idx];
            const double fvr = 0.;//c_vr[0]*fb[idx]+c_vr[1]*fc[idx]+c_vr[2]*ft[idx];
            idx++;

#define KAHAN_ADD(a,b) {a[2] = b - a[1]; a[3] = a[0] + a[2]; a[1] = (a[3] - a[0]) - a[2]; a[0] = a[3];}
//#define KAHAN_ADD(a,b) {a[0] += b;}
            // Add element's contribution to the conserved moments
            KAHAN_ADD(int_vx0vr0, dV*(c_vx0vr0[0]*f0+c_vx0vr0[1]*fvx+c_vx0vr0[2]*fvr));
            KAHAN_ADD(int_vx1vr0, dV*(c_vx1vr0[0]*f0+c_vx1vr0[1]*fvx+c_vx1vr0[2]*fvr));
            KAHAN_ADD(int_vx2vr0, dV*(c_vx2vr0[0]*f0+c_vx2vr0[1]*fvx+c_vx2vr0[2]*fvr));
            KAHAN_ADD(int_vx3vr0, dV*(c_vx3vr0[0]*f0+c_vx3vr0[1]*fvx+c_vx3vr0[2]*fvr));

            KAHAN_ADD(int_vx0vr2, dV*(c_vx0vr2[0]*f0+c_vx0vr2[1]*fvx+c_vx0vr2[2]*fvr));
            KAHAN_ADD(int_vx1vr2, dV*(c_vx1vr2[0]*f0+c_vx1vr2[1]*fvx+c_vx1vr2[2]*fvr));
#undef KAHAN_ADD
        }
    }

    // Reference the proper place in n,vx,vth to place data - offset since moments are in global frame
    double & n = _n.getData(index)[0];
    double & ux = _ux.getData(index)[0];
    double & uth = _uth.getData(index)[0];
    double & qx = _qx.getData(index)[0];

//    printf("vx0vr0 % e, vx1vr0 % e, vx2vr0 = % e, vx3vr0 = % e, vx0vr2 = % e, vx1vr2 = % e \n",int_vx0vr0[0], int_vx1vr0[0], int_vx2vr0[0], int_vx3vr0[0], int_vx0vr2[0], int_vx1vr2[0]);

    if(i == 0 && !_offsetIsSet){

        _offset_density = _config.n_l / int_vx0vr0[0];
        _offset_pressurexx = (_config.n_l * (_config.T_l + _config.m*_config.vx_l*_config.vx_l)) / int_vx2vr0[0];
        _offset_pressurerr = (2.*_config.n_l*_config.T_l) / int_vx0vr2[0];
        _offsetIsSet=true;
    }

    // Set n
    n = _offset_density*int_vx0vr0[0];

    // Calculate ux from conserved moments
    ux = int_vx1vr0[0] / n;

    // Calculate pressures
    const double Pxx = _offset_pressurexx*(int_vx2vr0[0] - n*ux*ux);
    const double Prr = _offset_pressurerr*(int_vx0vr2[0]);

//    printf("Pxx %e, Prr/2 %e\n",Pxx,Prr/2.);

    // Calculate vth from pressures
    uth = std::sqrt((Pxx+Prr) / (3.*n));

    // Caclulate conserved heat flux
    const double Hx = int_vx3vr0[0] + int_vx1vr2[0];

    // Calculate qx from other moments
    qx = 0.5*(Hx - ux*(Pxx+Prr) - 2.*ux*Pxx - n*ux*ux*ux);

}


void
QKSolver_BoltzmannBGK_1D2V_ShockTube_Var::writeVTK(const std::string & prefix, const std::string & suffix) const
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

        file << "DATASET RECTILINEAR_GRID\n";
        file << "DIMENSIONS " << iRange.getLength(0)+1 << " " << iRange.getLength(1)+1 << " " << iRange.getLength(2)+1 << "\n";
        //file << "ORIGIN " << _grid.getStart(0) + 0.5*_grid.getWidth(0) << " " << _grid.getStart(1) + 0.5*_grid.getWidth(1) << " " << _grid.getStart(2) + 0.5*_grid.getWidth(2) << "\n";
        //file << "SPACING " << _grid.getDx(0) << " " << _grid.getDx(1) << " " << _grid.getDx(2) << "\n\n";
        {
            const int dim = 0;
            file << "X_COORDINATES "<<_grid.getLength(dim)+1<<" float\n";
            double x = _grid.getStart(dim);
            file << x << " ";
            for(int i = 0; i < iRange.getLength(dim);i++){
                x += _grid.getDx(dim,i);
                file << x << " ";
            }
        }
        file << "\n\n";

        {
            const int dim = 1;
            file << "Y_COORDINATES "<<_grid.getLength(dim)+1<<" float\n";
            double x = _grid.getStart(dim);
            file << x << " ";
            for(int i = 0; i < iRange.getLength(dim);i++){
                x += _grid.getDx(dim,i);
                file << x << " ";
            }
        }
        file << "\n\n";

        {
            const int dim = 2;
            file << "Z_COORDINATES "<<_grid.getLength(dim)+1<<" float\n";
            double x = _grid.getStart(dim);
            file << x << " ";
            for(int i = 0; i < iRange.getLength(dim);i++){
                x += _grid.getDx(dim,i);
                file << x << " ";
            }
        }
        file << "\n\n";

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
        file << "SPACING " << _grid.getDx(0,0) << " 0 0\n\n";

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

