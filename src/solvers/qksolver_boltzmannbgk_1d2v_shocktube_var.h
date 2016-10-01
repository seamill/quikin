#ifndef QKSOLVER_BOLTZMANNBGK_1D2V_SHOCKTUBE_VAR_H
#define QKSOLVER_BOLTZMANNBGK_1D2V_SHOCKTUBE_VAR_H

// STL includes
#include <cmath>

// QK includes
#include "solvers/qksolver.h"
#include "qksolver_boltzmannbgk_1d2v_shocktube.h"


class QKSolver_BoltzmannBGK_1D2V_ShockTube_Var:
        public QKSolver
{
public:
    QKSolver_BoltzmannBGK_1D2V_ShockTube_Var();
    ~QKSolver_BoltzmannBGK_1D2V_ShockTube_Var();

    void setup(const QKStructuredGrid & grid, const QKSolver_BoltzmannBGK_1D2V_ShockTube_Configuration & config);

    void advance(const double time, const double dt);

    void writeVTK(const std::string & prefix, const std::string & suffix) const;

protected:

    void step(const double dt, const QKExtendedDatachunk & f_edc, QKExtendedDatachunk & fp_edc);

    void generateInitialConditions(QKExtendedDatachunk & f_edc);
    void generateMoments(const QKExtendedDatachunk & f, const int i);

    bool _offsetIsSet;
    double _offset_density;
    double _offset_pressurexx;
    double _offset_pressurerr;

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
    QKStructuredGrid _grid;

    // Storage for BGK collision operator
    QKExtendedDatachunk _vxterms;
    QKExtendedDatachunk _vrterms;

};

#endif // QKSOLVER_BOLTZMANNBGK_1D2V_SHOCKTUBE_VAR_H
