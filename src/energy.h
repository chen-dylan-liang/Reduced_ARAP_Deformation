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
        E += g*mass.transpose()*res.col(1);
        return E;
    }

    void compute_laplacian(){
        vector< Triplet<double> > tripletList;
        VectorXd row_sum;
        row_sum.resize(laplacian.rows());
        row_sum.fill(0);
        for(int i=0;i<anchors.size();i++)
            row_sum(anchors[i])=2*beta;
        for(int i=0;i<half_edges.size();i++){
            long long unsigned int a=half_edges[i].Endpoints(0);
            long long unsigned int b=half_edges[i].Endpoints(1);
            int inv_idx=half_edges[i].InverseIdx;
            double w=0;
            w=weights(a*laplacian.rows()+b);
            if(inv_idx!=-1){
                row_sum(a)+=w;
                tripletList.push_back(Triplet<double>(a,b,-w));
            }
            else{
                row_sum(a)+=w;
                row_sum(b)+=w;
                tripletList.push_back(Triplet<double>(a,b,-w));
                tripletList.push_back(Triplet<double>(b,a,-w));
            }
        }
        for(int i=0;i<laplacian.rows();i++) tripletList.push_back(Triplet<double>(i,i,row_sum(i)));
        laplacian.setFromTriplets(tripletList.begin(),tripletList.end());

    }

    void local_phase(const MatrixXd& res){
        int vert_num=verts.cols();
        for(int i=0;i<vert_num;i++){
            MatrixXd S=getCovariance3x3(neighbors[i],verts,res,weights,i);
            JacobiSVD<MatrixXd> SVD_solver;
            SVD_solver.compute(S,ComputeThinU | ComputeThinV);
            local_rotations[i]=SVD_solver.matrixU()*SVD_solver.matrixV().transpose();
            VectorXd singular_value=SVD_solver.singularValues();
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
                local_rotations[i]=SVD_solver.matrixU()*newV;
            }
            assert(local_rotations[i].determinant()>0);
        }
    }

    MatrixXd global_phase(){
        compute_rhs();
        return solver.solve(rhs);
    }

    LocalGlobalEnergy(string input_mesh, int method, double lambda, double mass_d, double g):method(method),lambda(lambda), g(g){
        int vert_num=readObj(input_mesh, faces, edges, verts, half_edges, area);
        laplacian.resize(vert_num,vert_num);
        local_rotations.reserve(vert_num);
        for(int i=0;i<vert_num;i++) {
            local_rotations[i].resize(3,3);
            local_rotations[i].setZero();
            local_rotations[i](0,0)= local_rotations[i](1,1)= local_rotations[i](2,2)=1;
        }
        weights.resize(vert_num*vert_num);
        getWeights(half_edges,verts,weights,Cotangent_2);
        sort(half_edges.begin(),half_edges.end(),compare);
        neighbors.reserve(vert_num);
        getNeighbors(half_edges,neighbors);
        rhs.resize(vert_num,3);
        mass = mass_d * VectorXd::Ones(vert_num);
    }
    const MatrixXi& faces(){return faces;}
private:
    // user-specified values or pre-computed values from mesh
    MatrixXd verts;
    MatrixXi edges;
    MatrixXi faces;
    VectorXd weights;
    VectorXd area;
    vector<HalfEdge> half_edges;
    vector<int> neighbors;
    double g;
    VectorXd mass;
    double lambda;
    int method;
    // anchors specifying during interaction
    vector<int>& anchors;
    vector<VectorXd>& anchor_vec;
    // global phase
    MatrixXd rhs;
    SparseMatrix<double> laplacian;
    void compute_rhs(){
        rhs.fill(0);
        for(int i=0;i<anchors.size();i++)
            rhs.row(anchors[i])=2*beta*anchor_vec[i].transpose();
        for(int i=0;i<half_edges.size();i++){
            long long unsigned int a=half_edges[i].Endpoints(0);
            long long unsigned  int b=half_edges[i].Endpoints(1);
            int inv_idx=half_edges[i].InverseIdx;
            MatrixXd coeR1;
            MatrixXd coeR2;
            if(inv_idx!=-1){
                rhs.row(a)+=(0.5*weights(a*rhs.rows()+b)*(local_rotations[a]+local_rotations[b])*half_edges[i].EdgeVec).transpose();
            }
            else{
                rhs.row(a)+=(0.5*weights(a*rhs.rows()+b)*(local_rotations[a]+local_rotations[b])*half_edges[i].EdgeVec).transpose();
                rhs.row(b)-=(0.5*weights(a*rhs.rows()+b)*(local_rotations[a]+local_rotations[b])*half_edges[i].EdgeVec).transpose();
            }
        }
        // gravity
        rhs.col(1) -= g*mass;
    }
    SimplicialLDLT<SparseMatrix<double>>& solver;
    // local phase
    vector<MatrixXd> local_rotations;
};

#endif /* ARAP_h */
