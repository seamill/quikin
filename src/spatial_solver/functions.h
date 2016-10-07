


#define _QK_FV_UPWIND_RECON_SOLVE_START(dim, rhs_size, q_data, rhs_data)      \
    const double d = grid.dx(dim);                              \
    const int stride = q_data.stride(dim);                      \
    qk::range itr_range(q_data.internal_range());               \
    itr_range.expand(dim, -1, 0);                               \
    itr_range.expand(itr_range.num_dims()-1, 0, 1);             \
    const double * __restrict__ q_REST = q_data.data();         \
    double * __restrict__ rhs_REST = rhs_data.data();           \
    double ql[rhs_size],qr[rhs_size],F[rhs_size];               \
    qk::indexer idx = q_data.indexer(itr_range);                \
    const int imin = idx.linear_index();                        \
    const int imax = idx.final_linear_index()+1;                \
    for (int i = imin; i < imax; i+=rhs_size) {                 \
        const int in = i - stride;                              \
        const int ip = i + stride;                              \
        const int ipp = ip + stride;                            \
        for (int j = 0; j < rhs_size; ++j) {                    \
            const double & qn = q_REST[in + j];                 \
            const double & qc = q_REST[i + j];                  \
            const double & qp = q_REST[ip + j];                 \
            const double & qpp = q_REST[ipp + j];               \
            const double rl = (qp - qc) / (qc - qn + 1.e-12);   \
            const double rr = (qp - qc) / (qpp - qp + 1.e-12);  \
            const double phil = (rl * rl + rl) / (rl * rl + 1); \
            const double phir = (rr * rr + rr) / (rr * rr + 1); \
            ql[j] = qc + 0.5 * phil * (qc - qn);                \
            qr[j] = qp + 0.5 * phir * (qp - qpp);               \
        }

#define _QK_FV_UPWIND_RECON_SOLVE_END(rhs_size)                 \
    for (int j = 0; j < rhs_size; ++j) {                        \
        rhs_REST[i + j] -= F[j] / d;                            \
        rhs_REST[ip + j] += F[j] / d;                           \
    }                                                           \
}
