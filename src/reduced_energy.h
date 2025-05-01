//
// Created by dylanmac on 4/29/25.
//

#ifndef REDUCED_ARAP_DEFORMATION_REDUCED_ENERGY_H
#define REDUCED_ARAP_DEFORMATION_REDUCED_ENERGY_H
#include "energy.h"

class ReducedLocalGlobalEnergy: public LocalGlobalEnergy{
    public:
        ReducedLocalGlobalEnergy(string input_mesh, int method, double lambda, double mass, double g, double dt, VectorXd offset,int subspace_dim,
                                 bool restore_rest_pose=true):
                                 LocalGlobalEnergy(input_mesh, method, lambda, mass, g,dt, offset), linearly_precise(restore_rest_pose),
                                 subspace_dim(subspace_dim){

            S.resize(subspace_dim, verts.cols());
            T.resize((verts.cols()-subspace_dim),verts.cols());
            Q.resize(verts.cols(), verts.cols());
            compute_q();
        }

        void local_global_solve() override{
            local_phase();
            res = subspace*global_phase();
            for(int i=0;i<anchors.size();i++)
                res.row(anchors[i]) = anchor_points[i].transpose();  // re-pin the anchors
        }

        void solver_compute() override{
            compute_subspace();
            reduced_solver.compute(subspace.transpose()*laplacian*subspace);
        }

    private:
        MatrixXd global_phase() override{
            return reduced_solver.solve(subspace.transpose()*compute_rhs());
        }
        void compute_subspace(){
            sample_subspace();
            vector< Triplet<double> > tripletListS, tripletListT;
            int si = 0, ti = 0;
            for (int i = 0; i < verts.cols(); ++i) {
                if (subspace_index_set.find(i)!=subspace_index_set.end()) {
                   tripletListS.push_back(Triplet<double>(si++, i, 1.0));
                }
                else {
                    tripletListT.push_back(Triplet<double>(ti++, i, 1.0));
                }
            }
            S.setFromTriplets(tripletListS.begin(), tripletListS.end());
            T.setFromTriplets(tripletListT.begin(), tripletListT.end());
            SparseMatrix<double> Qtt = T * Q * T.transpose();
            SparseMatrix<double> Qts = T * Q * S.transpose();
            solver.compute(Qtt);
            MatrixXd X = solver.solve(Qts);
            subspace = S.transpose() - T.transpose() * X;
            for (int idx : anchors) {
                subspace.row(idx).setZero();
            }
        }
        void sample_subspace(){
            // sample subspace point indices
            subspace_index_set.clear();
            int max_size = std::min(subspace_dim, (int)(verts.cols()-anchors.size()));
            while(subspace_index_set.size()<max_size){
                int idx = std::rand()%verts.cols();
                if(subspace_index_set.find(idx)==subspace_index_set.end()&&anchor_set.find(idx)==anchor_set.end())
                    subspace_index_set.insert(idx);
            }
        }

        void compute_q(){
            SparseMatrix<double> L;
            vector< Triplet<double> > tripletList;
            // compute L. Graph Laplacian of mesh
            L.resize(verts.cols(), verts.cols());
            VectorXd row_sum;
            row_sum.resize(verts.cols());
            row_sum.fill(0);
            for(int i=0;i<half_edges.size();i++){
                long long unsigned int a=half_edges[i].Endpoints(0);
                long long unsigned int b=half_edges[i].Endpoints(1);
                int inv_idx=half_edges[i].InverseIdx;
                double w=weights(a*laplacian.rows()+b);
                row_sum(a)+=w;
                tripletList.push_back(Triplet<double>(a,b,-w));
                // ab on the boundary
                if(inv_idx==-1){
                    row_sum(b)+=w;
                    tripletList.push_back(Triplet<double>(b,a,-w));
                }
            }
            for(int i=0;i<row_sum.rows();i++) tripletList.push_back(Triplet<double>(i,i,row_sum(i)));
            L.setFromTriplets(tripletList.begin(),tripletList.end());

            if(!linearly_precise){
                Q = L.transpose()*L;
            }
            else{
                // compute N, the normal derivatives
                SparseMatrix<double> N;
                tripletList.clear();
                N.resize(verts.cols(), verts.cols());
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
                N.setFromTriplets(tripletList.begin(), tripletList.end());
                Q = (L+N).transpose()*(L+N);
            }
        }

        // h
        int subspace_dim;
        // n x h
        MatrixXd subspace;
        // related to subspace computations
        SparseMatrix<double> Q;
        // h x n
        SparseMatrix<double> S;
        // (n-h) x n
        SparseMatrix<double> T;
        // of size h
        std::set<int> subspace_index_set;
        // of size n-h
        //MatrixXd J=Identity for anchors only;
        Eigen::LDLT<MatrixXd> reduced_solver;
        // whether being able to restore the rest pose
        bool linearly_precise;
};
#endif
