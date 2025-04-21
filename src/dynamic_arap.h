//
// Created by dylanmac on 4/21/25.
//

#ifndef REDUCED_ARAP_DEFORMATION_DYNAMIC_ARAP_H
#define REDUCED_ARAP_DEFORMATION_DYNAMIC_ARAP_H
#include "ARAP.h"

class DynamicARAP: public ARAP_energy {
    DynamicARAP(vector<HalfEdge>& half_edges,VectorXd& weights,VectorXd& area,int method,double lambda, const MatrixXd& verts):
    ARAP_energy(half_edges, weights, area, method, lambda),prev_res(verts){}
    double operator()(MatrixXd& res) override {
        double E = ARAP_energy::operator()(res);

        const double h = 0.1;
        const double g = 9.8;
        const double m = 1.0;
        const double C = m / (2 * h * h);

        MatrixXd pred = 2*res - prev_res;
        pred.col(2).array() -= h*h * g;
        MatrixXd diff = res - pred;
        double Ek = diff.squaredNorm();
        E += C * Ek;
        prev_res = res;
        return E;
    }

private:
    MatrixXd prev_res;
};

inline void getRHSDynamics(MatrixXd* R, const MatrixXd& prev_res, const MatrixXd& res, const vector<HalfEdge>& half_edges,const VectorXd& weights,const vector<int>& fix, const vector<VectorXd>& fix_vec,MatrixXd& RHS,int func,bool is_first){
    getRHS(R,half_edges,weights,fix,fix_vec,RHS,func,is_first);
    const double h = 0.01;
    const double g = 1.0;
    const double m = 0.005;
    //const double C = m / (2 * h * h);

    RHS.col(1).array() -= g*m;
}
inline void global_phase_dynamics(MatrixXd* R, MatrixXd& prev_res, const vector<HalfEdge>& half_edges, const VectorXd& weights, const vector<int>& fix,const vector<VectorXd>& fix_vec,MatrixXd& RHS,int func,bool first, MatrixXd& res,SimplicialLDLT<SparseMatrix<double> >& dir_solver){
    getRHSDynamics(R, prev_res, res, half_edges,weights,fix,fix_vec,RHS,func,first);
    prev_res = res;
    res=dir_solver.solve(RHS);
}



#endif
