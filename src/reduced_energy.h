//
// Created by dylanmac on 4/29/25.
//

#ifndef REDUCED_ARAP_DEFORMATION_REDUCED_ENERGY_H
#define REDUCED_ARAP_DEFORMATION_REDUCED_ENERGY_H
#include "energy.h"
#include <map>

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
            subspace_indices.clear();
            subspace_index_vec.clear();
            complement_index_vec.clear();
            for(int i=0; i<anchors.size(); i++) subspace_indices[i]=anchors[i];
            while(subspace_indices.size()<subspace_dim-anchors.size()){
                int idx = std::rand()%verts.size();
                if(subspace_indices.find(idx)==subspace_indices.end()){
                    int now_size = subspace_indices.size();
                    subspace_indices[now_size] = idx;
                }
            }
            for(int i=0; i<verts.size();i++){
                if(subspace_indices.find(i)==subspace_indices.end()) complement_index_vec.push_back(i);
                else subspace_index_vec.push_back(i);
            }
            // compute subspace
            S.setZero();
            for(int i=0;i<subspace_dim;i++) S(i,subspace_index_vec[i])=1;
            T.setZero();
            for(int i=0; i<complement_index_vec.size();i++) T(i, complement_index_vec[i])=1;
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
        std::map<int,int> subspace_indices;
        std::vector<int> subspace_index_vec;
        std::vector<int> complement_index_vec;
        // of size n-h
        //MatrixXd J=Identity for anchors only;
        LDLT<MatrixXd> reduced_solver;
};
#endif
