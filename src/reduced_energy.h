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
            // sample subspace point indices
            subspace_index_set.clear();
            for(int i=0; i<anchors.size(); i++) subspace_index_set.insert(anchors[i]);
            while(subspace_index_set.size()<subspace_dim-anchors.size()){
                int idx = std::rand()%verts.cols();
                if(subspace_index_set.find(idx)==subspace_index_set.end()) subspace_index_set.insert(idx);
            }
            int subspace_cnt = 0, complement_cnt=0;
            S.setZero();
            T.setZero();
            for(int i=0; i<verts.cols();i++){
                if(subspace_index_set.find(i)==subspace_index_set.end()){
                    T(complement_cnt, i)=1;
                    complement_cnt++;
                }
                else{
                    S(subspace_cnt, i)=1;
                    subspace_cnt++;
                }
            }
            // compute subspace
            MatrixXd Qtt = T * Q * T.transpose();
            MatrixXd Qts = T * Q * S.transpose();
            subspace = S.transpose() - T.transpose() * Qtt.inverse()*Qts;
        }
        // h
        int subspace_dim;
        // n x h
        MatrixXd subspace;
        // related to subspace computations
        MatrixXd Q;
        // h x n
        MatrixXd S;
        // (n-h) x n
        MatrixXd T;
        // of size h
        std::set<int> subspace_index_set;
        // of size n-h
        //MatrixXd J=Identity for anchors only;
        Eigen::LDLT<MatrixXd> reduced_solver;
};
#endif
