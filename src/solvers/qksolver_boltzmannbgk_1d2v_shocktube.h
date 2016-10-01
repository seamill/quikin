#ifndef QKSOLVER_BOLTZMANNBGK_1D2V_SHOCKTUBE_H
#define QKSOLVER_BOLTZMANNBGK_1D2V_SHOCKTUBE_H

// STL includes
#include <cmath>

// QK includes
#include "solvers/qksolver.h"

struct QKSolver_BoltzmannBGK_1D2V_ShockTube_Configuration
{

    QKSolver_BoltzmannBGK_1D2V_ShockTube_Configuration():
        T_l(1.0),
        T_r(1.0),
        n_l(1.0),
        n_r(0.125),
        vx_l(0.0),
        vx_r(0.0),
        m(1),
        nu_tau(100.0),
        transitionWidth(0.0)
    {

    }

    double vth_max(){
        return std::sqrt(std::max(T_l, T_r) / m);
    }

    double T_l;
    double T_r;
    double n_l;
    double n_r;
    double vx_l;
    double vx_r;

    double m;

    double nu_tau;
    double transitionWidth;


};

class QKSolver_BoltzmannBGK_1D2V_ShockTube:
        public QKSolver
{
public:
    QKSolver_BoltzmannBGK_1D2V_ShockTube();
    ~QKSolver_BoltzmannBGK_1D2V_ShockTube();

    void setup(const QKRectilinearGrid & grid, const QKSolver_BoltzmannBGK_1D2V_ShockTube_Configuration & config);

    void advance(const double time, const double dt);

    void writeVTK(const std::string & prefix, const std::string & suffix) const;

protected:
    struct integral_set_t{

        void setup(double dv, double vmin, double num);

        double _w0v0;
        std::vector<double> _w0v1;
        std::vector<double> _w0v2;
        std::vector<double> _w0v3;
        double _w1v0;
        double _w1v1;
        std::vector<double> _w1v2;
        std::vector<double> _w1v3;
        double _w2v0;
        std::vector<double> _w2v1;
        std::vector<double> _w2v2;
        std::vector<double> _w2v3;
    };

    void step(const double dt, const QKExtendedDatachunk & f_edc, QKExtendedDatachunk & fp_edc);

    void generateInitialConditions(QKExtendedDatachunk & f_edc);
    void generateMoments(const QKExtendedDatachunk & f, const int i);

    // Time step variables
    QKDataset _f;
    QKDataset _fs;
    QKDataset _fp;

    // Moments - locally stored
    QKExtendedDatachunk _n;
    QKExtendedDatachunk _ux;
    QKExtendedDatachunk _uth;
    QKExtendedDatachunk _qx;

    // Configuration
    QKSolver_BoltzmannBGK_1D2V_ShockTube_Configuration _config;

    // Grid
    QKRectilinearGrid _grid;

    // Integration stuff
    integral_set_t _i_vx;
    integral_set_t _i_vr;

    // Storage for BGK collision operator
    QKExtendedDatachunk _vxterms;
    QKExtendedDatachunk _vrterms;

};

#endif // QKSOLVER_BOLTZMANNBGK_1D2V_SHOCKTUBE_H
