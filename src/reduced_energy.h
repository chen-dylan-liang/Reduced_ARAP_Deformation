//
// Created by dylanmac on 4/29/25.
//

#ifndef REDUCED_ARAP_DEFORMATION_REDUCED_ENERGY_H
#define REDUCED_ARAP_DEFORMATION_REDUCED_ENERGY_H
#include "energy.h"

class ReducedLocalGlobalEnergy: public LocalGlobalEnergy{
    public:
        ReducedLocalGlobalEnergy(string input_mesh, int method, double lambda, double mass_d, double g, int subspace_dim,
                                 bool restore_rest_pose=false):
                                 LocalGlobalEnergy(input_mesh, method, lambda, mass_d, g), linearly_precise(restore_rest_pose),
                                 specified_subspace_dim(subspace_dim){

            S.resize(subspace_dim, verts.cols());
            T.resize((verts.cols()-subspace_dim),verts.cols());
            Q.resize(verts.cols(), verts.cols());
            vector< Triplet<double> > tripletList;
            for(int i=0; i<half_edges.size(); i++)
                if(half_edges[i].InverseIdx==-1){
                int a =half_edges[i].Endpoints(0);
                int b = half_edges[i].Endpoints(1);
                int c =half_edges[i].OppositePoint;
                double w_bc = weights(b*laplacian.rows()+c);
                double w_ca = weights(c*laplacian.rows()+a);
                tripletList.push_back(Triplet<double>(a,a,w_ca));
                tripletList.push_back(Triplet<double>(b,b,w_bc));
                tripletList.push_back(Triplet<double>(a,b,w_bc));
                tripletList.push_back(Triplet<double>(b,a,w_ca));
                tripletList.push_back(Triplet<double>(a,c,-w_ca-w_bc));
            }
            normal_derivates.setFromTriplets(tripletList.begin(), tripletList.end());
        }

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
            sample_subspace();
            int subspace_cnt = 0, complement_cnt=0;
            int subspace_dim = subspace_index_set.size();
            if(subspace_dim>specified_subspace_dim){
                S.resize(subspace_dim, verts.cols());
                T.resize((verts.cols()-subspace_dim),verts.cols());
            }
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
            compute_q();
            MatrixXd Qtt = T * Q * T.transpose();
            MatrixXd Qts = T * Q * S.transpose();
            subspace = S.transpose() - T.transpose() * Qtt.inverse()*Qts;
        }
        void sample_subspace(){
            // sample subspace point indices
            subspace_index_set.clear();
            for(int i=0; i<anchors.size(); i++) subspace_index_set.insert(anchors[i]);
            while(subspace_index_set.size()<specified_subspace_dim-anchors.size()){
                int idx = std::rand()%verts.cols();
                if(subspace_index_set.find(idx)==subspace_index_set.end()) subspace_index_set.insert(idx);
            }
        }
        void compute_q(){
            if(!linearly_precise){
                Q = laplacian.transpose()*laplacian;
            }
            else{
                Q = (laplacian+normal_derivates).transpose()*(laplacian+normal_derivates);
            }
        }
        // h
        int specified_subspace_dim;
        // n x h
        MatrixXd subspace;
        // related to subspace computations
        MatrixXd Q;
        SparseMatrix<double> normal_derivates;
        // h x n
        MatrixXd S;
        // (n-h) x n
        MatrixXd T;
        // of size h
        std::set<int> subspace_index_set;
        // of size n-h
        //MatrixXd J=Identity for anchors only;
        Eigen::LDLT<MatrixXd> reduced_solver;
        // whether being able to restore the rest pose
        bool linearly_precise;
};
#endif
