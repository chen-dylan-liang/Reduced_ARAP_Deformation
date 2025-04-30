//
// Created by dylanmac on 4/29/25.
//

#ifndef REDUCED_ARAP_DEFORMATION_REDUCED_ENERGY_H
#define REDUCED_ARAP_DEFORMATION_REDUCED_ENERGY_H
#include "energy.h"

class ReducedLocalGlobalEnergy: public LocalGlobalEnergy{
    public:
        void local_global_solve() override{
            local_phase();
            compute_subspace();
            res = subspace*global_phase();
        }
        void solver_compute() override{
            reduced_solver.compute(subspace.transpose()*laplacian*subspace);
        }
    private:
        MatrixXd global_phase() override{
            return reduced_solver.solve(subspace.transpose()*compute_rhs());
        }
        void compute_subspace(){
            // farthest sampling
            S.setZero();
            for(int i=0;i<subspace_indices.size();i++) S(i,subspace_indices[i])=1;
            T.setZero();
            for(int i=0; i<complement_indices.size();i++) T(i, complement_indices[i])=1;
            MatrixXd Qtt = T * Q * T.transpose();
            MatrixXd Qts = T * Q * S.transpose();
            subspace = S.transpose() - T.transpose() * Qtt.inverse()*Qts;
        }
        // n x h
        MatrixXd subspace;
        // related to subspace computations
        MatrixXd Q;
        // h x n
        MatrixXd S;
        // (n-h) x n
        MatrixXd T;
        // of size h
        vector<int> subspace_indices;
        // of size n-h
        vector<int> complement_indices;
        //MatrixXd J=Identity for anchors only;
        LDLT<MatrixXd> reduced_solver;
};
#endif
