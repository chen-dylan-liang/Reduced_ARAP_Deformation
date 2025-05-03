//
//  ARAP.h
//  arap
//
//  Created by Liang Chen on 2021/6/14.
//

#ifndef ARAP_h
#define ARAP_h

#include "helper_algebra.h"
#include "helper_geometry.h"
#include "helper_init.h"
#include "Eigen/Dense"
#include "Eigen/Sparse"
#include <Eigen/CholmodSupport>
#include <Eigen/SparseCholesky>
#include <set>

class LocalGlobalEnergy{
    public:
        double operator()(const MatrixXd& res){
            double E=0;
            MatrixXd L(2,2);
            for(int i=0;i<half_edges.size()/3;i++){
                //get cross-covariance matrix or jacobi matrix from (u,v) to (x,y)
                MatrixXd J=getJacobian2x2(half_edges,res,i);
        //            assert(J.determinant()<0);
                // J=U*Sigma*V^T
               JacobiSVD<MatrixXd> SVD_solver;
               if(method==ARAP||method==ASAP) SVD_solver.compute(J,ComputeThinU | ComputeThinV);
               else if(method==Hybrid) SVD_solver.compute(J);
                // R=U*V^T
                Vector2d singular_value=SVD_solver.singularValues();
                MatrixXd SI(2,2);
                SI.fill(0);
                SI(0,0)=SI(1,1)=0.5*(singular_value(0)+singular_value(1));
                if(method==ARAP||method==ASAP){
                if(method==ARAP) L=SVD_solver.matrixU()*SVD_solver.matrixV().transpose();
                else if(method==ASAP) L=SVD_solver.matrixU()*SI*SVD_solver.matrixV().transpose();
              if( L.determinant()<0){
                  //cout<<R[i]<<endl;
                  Matrix2d newV;
                  newV <<SVD_solver.matrixV().transpose()(0, 0), SVD_solver.matrixV().transpose()(0, 1),
                     -(SVD_solver.matrixV().transpose()(1, 0)), -(SVD_solver.matrixV().transpose()(1, 1));
                  SI(0,0)=SI(1,1)=0.5*(singular_value(0)-singular_value(1));
                  if(method==ARAP) L=SVD_solver.matrixU()*newV;
                  else if(method==ASAP) L=SVD_solver.matrixU()*SI*newV;
                }
                }
                if(method==Hybrid){
                    double C1=0;
                    double C2=0;
                    double C3=0;
                    for(int j=0;j<3;j++){
                        Vector2d v=half_edges[i*3+j].EdgeVec;
                        int a=half_edges[i*3+j].Endpoints(0);
                        int b=half_edges[i*3+j].Endpoints(1);
                        Vector2d u=(res.row(a)-res.row(b)).transpose();
                        C1+=weights(i*3+j)*(v(0)*v(0)+v(1)*v(1));
                        C2+=weights(i*3+j)*(u(0)*v(0)+u(1)*v(1));
                        C3+=weights(i*3+j)*(u(0)*v(1)-u(1)*v(0));
                    }
                    /*
                    VectorXd init(2);
                    init(0)=init(1)=1;
                    R_Jacobian R_J(lamda,C1,C2,C3);
                    R_Hessian R_H(lamda,C1);
                    newton_optimizerNto1(init,R_J,R_H);*/
                    double entry_1=0;
                    double coe3=2*lambda*(1+C3*C3/C2/C2),coe1=C1-2*lambda,coe0=-C2;
                    std::function<double(double)> cubic_fun=[&](double x)->double
                    {
                        return coe3*x*x*x+coe1*x+coe0;
                    };
                    binary_find_root(entry_1,cubic_fun);
                    L(0,0)=L(1,1)=entry_1;
                    L(0,1)=entry_1*C3/C2;
                    L(1,0)=-entry_1*C3/C2;
                }
                E+=area(i)* ( ( (J-L).transpose() )*(J-L) ).trace() ;
            }
            // gravity
            if(g!=0.0){
                MatrixXd pred_res = 2 * res - prev_res;
                pred_res.col(1).array() -= g*dt*dt;
                E += ((res-pred_res).transpose()*(res-pred_res)).trace()/dt/dt/2;
            }
            return E;
        }

        int vert_size() const{
            return verts.cols();
        }
        virtual void local_global_solve(){
            local_phase();
            res=global_phase();
        }

        void compute_laplacian(){
            vector< Triplet<double> > tripletList;
            VectorXd row_sum;
            row_sum.resize(laplacian.rows());
            row_sum.fill(0);
            for(int i=0;i<half_edges.size();i++){
                long long unsigned int a=half_edges[i].Endpoints(0);
                long long unsigned int b=half_edges[i].Endpoints(1);
                int inv_idx=half_edges[i].InverseIdx;
                bool a_is_free = anchor_set.find(a)==anchor_set.end();
                bool b_is_free = anchor_set.find(b)==anchor_set.end();
                double w=weights(a*laplacian.rows()+b);
                if(a_is_free) row_sum(a)+=w;
                // knock out L(a,b) when b is anchor
                if(a_is_free && b_is_free) tripletList.push_back(Triplet<double>(a,b,-w));
                    // ab on the boundary
                    if(inv_idx==-1){
                        if(b_is_free) row_sum(b)+=w;
                        // knock out L(b,a) when a is anchor
                        if(a_is_free && b_is_free) tripletList.push_back(Triplet<double>(b,a,-w));
                     }
            }
            if(g!=0.0){
            // gravity
            for(int i=0; i<vert_size(); i++) row_sum(i)+=mass/dt/dt;
            }
            // set the diagonals of anchors to 1
            for(int i=0; i<anchors.size();i++) row_sum(anchors[i])=1;
            for(int i=0;i<laplacian.rows();i++) tripletList.push_back(Triplet<double>(i,i,row_sum(i)));
            laplacian.setFromTriplets(tripletList.begin(),tripletList.end());
        }


        virtual void solver_compute(){
            solver.compute(laplacian);
        }

        LocalGlobalEnergy(string input_mesh, int method, double lambda, double mass, double g, double dt, VectorXd offset):method(method),lambda(lambda), mass(mass),g(g),
        dt(dt){
            int vert_num=readObj(input_mesh, faces, edges, verts, half_edges, area);
            for(int i=0; i<vert_num; i++) verts.col(i) += offset;
            laplacian.resize(vert_num,vert_num);
            local_rotations.reserve(vert_num);
            res.resize(vert_num,3);
            prev_res.resize(vert_num,3);
            res = verts.transpose();
            prev_res = verts.transpose();
            for(int i=0;i<vert_num;i++) {
                local_rotations.push_back(MatrixXd::Zero(3,3));
                local_rotations[i](0,0)= local_rotations[i](1,1)= local_rotations[i](2,2)=1;
            }
            weights.resize(vert_num*vert_num);
            getWeights(half_edges,verts,weights,Cotangent_2);
            sort(half_edges.begin(),half_edges.end(),compare);
            neighbors.reserve(vert_num);
            for(int i=0;i<vert_num;i++) {
                neighbors.push_back(vector<int>());
            }
            getNeighbors(half_edges, neighbors);
        }
        // getters
        const MatrixXi& get_faces() const {return faces;}
        MatrixXi get_faces() {return faces;}
        const MatrixXd& get_res()const {return res;}
        const vector<int>& get_anchors() const {return anchors;}
        const vector<VectorXd>& get_anchor_points()  const {return anchor_points;}
        vector<int>& get_anchors() {return anchors;}
        vector<VectorXd>& get_anchor_points()  {return anchor_points;}
        bool anchor_in_region(int i) {return (bool)anchor_in_batch[i];}
        // setters
        void add_anchor(int idx, Eigen::Vector3d point, bool batch_op=false){
            if(anchor_set.find(idx)==anchor_set.end()) {
              anchor_set.insert(idx);
              anchors.push_back(idx);
              anchor_points.push_back(point);
              anchor_in_batch.push_back((int)batch_op);
            }
        }
        void clear_anchors(){
            anchors.clear();
            anchor_set.clear();
            anchor_points.clear();
            anchor_in_batch.clear();
        }

        void restore_rest_pose(){
            for(int i=0 ;i<anchor_points.size();i++) anchor_points[i]=verts.col(anchors[i]).transpose();
            //res = verts.transpose();
        }

    protected:
        // user-specified values or pre-computed values from mesh
        MatrixXd verts;
        MatrixXi edges;
        MatrixXi faces;
        VectorXd weights;
        VectorXd area;
        vector<HalfEdge> half_edges;
        vector<vector<int>> neighbors;
        double g;
        double mass;
        double lambda;
        double dt;
        int method;
        // anchors specifying during interaction
        vector<int> anchors;
        set<int> anchor_set;
        vector<VectorXd> anchor_points;
        vector<int> anchor_in_batch;
        // global phase
        MatrixXd res;
        MatrixXd prev_res;
        SparseMatrix<double> laplacian;

        virtual MatrixXd global_phase(){
            return solver.solve(compute_rhs());
        }

        MatrixXd compute_rhs(){
            MatrixXd rhs;
            rhs.resize(verts.cols(),3);
            rhs.fill(0);
            for(int i=0;i<anchors.size();i++){
                rhs.row(anchors[i])=anchor_points[i].transpose();
            }
            for(int i=0;i<half_edges.size();i++){
                long long unsigned int a=half_edges[i].Endpoints(0);
                long long unsigned  int b=half_edges[i].Endpoints(1);
                    int inv_idx=half_edges[i].InverseIdx;
                    MatrixXd coeR1;
                    MatrixXd coeR2;
                bool a_is_free = anchor_set.find(a)==anchor_set.end();
                bool b_is_free = anchor_set.find(b)==anchor_set.end();
                double w=weights(a*laplacian.rows()+b);
                if(a_is_free){
                    rhs.row(a)+=(0.5*w*(local_rotations[a]+local_rotations[b])*half_edges[i].EdgeVec).transpose();
                    // a is free and connected to an anchor
                    if(!b_is_free) rhs.row(a)+= w*rhs.row(b);
                }
                // ab on the boundary
                if(inv_idx==-1){
                    if(b_is_free){
                        rhs.row(b)-=(0.5*w*(local_rotations[a]+local_rotations[b])*half_edges[i].EdgeVec).transpose();
                    // b is free and connected to an anchor
                    if(!a_is_free) rhs.row(b)+= w*rhs.row(a);
                    }
                }
            }
            if(g!=0.0){
                // gravity m/dt/dt * (2Xn-Xn-1+dt*dt*g)
                MatrixXd pred_res = (2 * res - prev_res)*mass/dt/dt;
                pred_res.col(1).array() -= mass*g;
                // zero out anchor entries for pred_res
                for(int i=0;i<anchors.size();i++){
                    pred_res.row(anchors[i])=RowVector3d::Zero();
                }
                rhs = rhs + pred_res;
            }
            // update prev_res
            prev_res = res;
            return rhs;
        }
    SimplicialLDLT<SparseMatrix<double>> solver;
        // local phase
        vector<MatrixXd> local_rotations;
        void local_phase(){
            int vert_num=verts.cols();
            for(int i=0;i<vert_num;i++){
                MatrixXd S=getCovariance3x3(neighbors[i],verts,res,weights,i);
                JacobiSVD<MatrixXd> SVD_solver;
                SVD_solver.compute(S,ComputeThinU | ComputeThinV);
                // similarity matrix for asap
                VectorXd singular_value=SVD_solver.singularValues();
                Matrix3d similarity;
                similarity.fill(0);
                similarity(0,0)=similarity(1,1)=similarity(2,2)=singular_value.sum()/3.0;
                if(method==ARAP) local_rotations[i]=SVD_solver.matrixU()*SVD_solver.matrixV().transpose();
                else if(method==ASAP) local_rotations[i]=SVD_solver.matrixU()*similarity*SVD_solver.matrixV().transpose();
                if(local_rotations[i].determinant()<0){
                    double smallest_sv=1e16;
                    double smallest_sv_idx=-1;
                    for(int j=0;j<3;j++)
                        if(singular_value(j)<smallest_sv){
                            smallest_sv=singular_value(j);
                            smallest_sv_idx=j;
                        }
                    MatrixXd newV=SVD_solver.matrixV().transpose();
                    newV.row(smallest_sv_idx)=-1.0*newV.row(smallest_sv_idx);
                    if(method==ARAP) local_rotations[i]=SVD_solver.matrixU()*newV;
                    else if(method==ASAP){
                        double val=0.0;
                        for(int j=0;j<3;j++){
                            if(j==smallest_sv_idx) val-=singular_value(j);
                            else val+=singular_value(j);
                        }
                        similarity(0,0)=similarity(1,1)=similarity(2,2)=val;
                        local_rotations[i]=SVD_solver.matrixU()*similarity*newV;
                    }
                }
                assert(local_rotations[i].determinant()>0);
            }
        }
};

#endif /* ARAP_h */
