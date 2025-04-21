//
//  call_back.h
//  arap
//
//  Created by Liang Chen on 2021/6/14.
//

#ifndef call_back_h
#define call_back_h

#include "global_var.h"
#include "helper_print.h"
#include "ARAP.h"
#include "dynamic_arap.h"

MatrixXd prev_res;
inline void updateViewer(igl::opengl::glfw::Viewer& viewer, int& itr, int mode, MatrixXd* R, vector<int>* neighbors,const MatrixXd& verts,const vector<HalfEdge>& half_edges, const VectorXd& weights, vector<int>& fix,vector<VectorXd>& fix_vec,MatrixXd& RHS,bool& first,bool moved, MatrixXd& res,SimplicialLDLT<SparseMatrix<double> >& dir_solver)
 {
   // predefined colors
   const Eigen::RowVector3d orange(1.0,0.7,0.2);
   const Eigen::RowVector3d yellow(1.0,0.9,0.2);
   const Eigen::RowVector3d blue(0.2,0.3,0.8);
   const Eigen::RowVector3d green(0.2,0.6,0.3);
   if(mode==1)
   {
      viewer.data().set_vertices(res);
      viewer.data().set_colors(yellow);
      MatrixXd CV(fix.size(),3);
      for(int i=0;i<fix.size();i++)
      CV.row(i)=fix_vec[i].transpose();
      viewer.data().set_points(CV,green);
   }
   else if(mode==2){
       viewer.data().set_vertices(res);
       viewer.data().set_colors(orange);
       MatrixXd CV(fix.size(),3);
       for(int i=0;i<fix.size();i++)
           CV.row(i)=fix_vec[i].transpose();
       viewer.data().set_points(CV,green);
   }
   else{
       if(moved){
           //local and global optimizations
    if(!first) local_phase_deform(R,neighbors,verts,res,weights);
       global_phase_dynamics(R, prev_res, half_edges,weights,fix,fix_vec,RHS,DEFORM,first,res,dir_solver);
           if(first) first=false;
       for(int i=0;i<fix.size();i++) fix_vec[i]=(res.row(fix[i]).transpose());
           itr++;
       }
     viewer.data().set_vertices(res);
     viewer.data().set_colors(blue);
    MatrixXd CU(fix.size(),3);
    for(int i=0;i<fix.size();i++)
    CU.row(i)=fix_vec[i].transpose();
    viewer.data().set_points(CU,orange);
   }
   viewer.data().compute_normals();
 };


class callbackMouseDown{
public:
    bool operator() (igl::opengl::glfw::Viewer& viewer, int, int)
  {
    *last_mouse = Eigen::RowVector3f(viewer.current_mouse_x,viewer.core().viewport(3)-viewer.current_mouse_y,0);
    if(*mode==1)
    {
      // Find the closest point on mesh to mouse position
      int fid;
      Eigen::Vector3f bary;
      if(igl::unproject_onto_mesh((*last_mouse).head(2),viewer.core().view,viewer.core().proj,viewer.core().viewport,*res,*F,fid, bary))
      {
        long c;
        bary.maxCoeff(&c);
        Eigen::RowVector3d new_c = (*res).row((*F)(fid,c));
        MatrixXd CV((*fix).size(),3);
        for(int i=0;i<(*fix).size();i++)
          CV.row(i)=((*fix_vec)[i]).transpose();
        if((*fix).size()==0 || (CV.rowwise()-new_c).rowwise().norm().minCoeff() > 0)
        {
          (*fix).push_back((*F)(fid,c));
            (*fix_vec).push_back(new_c.transpose());
            *changed=true;
          updateViewer(viewer,*itr,*mode,R, neighbors,*verts,*half_edges, *weights, *fix,*fix_vec,*RHS,*first,*moved,*res,*dir_solver);
          return true;
        }
      }
    }
    else if(*mode==0)
    {
      // Move closest control point
      Eigen::MatrixXf CP;
        MatrixXd CU((*fix).size(),3);
        for(int i=0;i<(*fix).size();i++)
        CU.row(i)= (*fix_vec)[i].transpose();
      igl::project(Eigen::MatrixXf(CU.cast<float>()),viewer.core().view,viewer.core().proj, viewer.core().viewport, CP);
      Eigen::VectorXf D = (CP.rowwise()-(*last_mouse)).rowwise().norm();
      *sel = (D.minCoeff(sel) < 30)? (*sel):-1;
      if(*sel != -1)
      {
        (*last_mouse)(2) = CP(*sel,2);
        return true;
      }
    }
    return false;
  }
  callbackMouseDown(int* mode,
                    int* itr,
                    RowVector3f* last_mouse,
                    MatrixXd* res,
                    vector<int>* fix,
                    vector<VectorXd>* fix_vec,
                    MatrixXd* R,
                    vector<int>* neighbors,
                    MatrixXd* verts,
                    vector<HalfEdge>* half_edges,
                    VectorXd* weights,
                    MatrixXd* RHS,
                    MatrixXi* F,
                    bool* first,
                    bool* moved,
                    bool* changed,
                    long* sel,
                    SimplicialLDLT<SparseMatrix<double>>* dir_solver
                    ):mode(mode),itr(itr),res(res),fix(fix),fix_vec(fix_vec),R(R),neighbors(neighbors),verts(verts),half_edges(half_edges),weights(weights),RHS(RHS),F(F),first(first),sel(sel),dir_solver(dir_solver),last_mouse(last_mouse),moved(moved),changed(changed){}
private:
    int* mode;
    int* itr;
    MatrixXd* res;
    vector<int>* fix;
    vector<VectorXd>* fix_vec;
    MatrixXd* R;
    vector<int>* neighbors;
    MatrixXd* verts;
    vector<HalfEdge>* half_edges;
    VectorXd* weights;
    MatrixXd* RHS;
    MatrixXi* F;
    bool* first;
    bool* moved;
    bool* changed;
    long* sel;
    SimplicialLDLT<SparseMatrix<double>>* dir_solver;
    RowVector3f* last_mouse;
};

class callbackMouseMove{
public:
    bool operator() (igl::opengl::glfw::Viewer & viewer, int,int)
  {
    if(*sel!=-1)
    {
      Eigen::RowVector3f drag_mouse(viewer.current_mouse_x,viewer.core().viewport(3) - viewer.current_mouse_y,(*last_mouse)(2));
      Eigen::RowVector3f drag_scene,last_scene;
      igl::unproject(drag_mouse,viewer.core().view,viewer.core().proj,viewer.core().viewport,drag_scene);
      igl::unproject(*last_mouse,viewer.core().view,viewer.core().proj,viewer.core().viewport,last_scene);
      (*fix_vec)[*sel] +=(drag_scene-last_scene).cast<double>();
      *last_mouse = drag_mouse;
      *itr=0;
      *moved=true;
      return true;
    }
    if(*mode==2){
        Eigen::RowVector3f drag_mouse(viewer.current_mouse_x, viewer.core().viewport(3) - viewer.current_mouse_y, (*last_mouse)(2));
        Eigen::RowVector3f start(std::min(drag_mouse.x(), last_mouse->x()),std::min(drag_mouse.y(), last_mouse->y()), 0.0);
        Eigen::RowVector3f end(std::max(drag_mouse.x(), last_mouse->x()),std::max(drag_mouse.y(), last_mouse->y()), 0.0);
        fix_vec->clear();
        for(int i=0; i<res->rows();i++){
            Eigen::RowVector3f projected=igl::project(Eigen::Matrix<float, 3,1>(res->row(i).transpose().cast<float>()), viewer.core().view, viewer.core().proj,viewer.core().viewport);
            if(projected.x()>start.x()&&projected.x()<end.x()&&projected.y()>start.y()&&projected.y()<end.y()){
                fix_vec->push_back(res->row(i).transpose());
            }
        }
        MatrixXd CV(fix_vec->size(),3);
        const Eigen::RowVector3d green(0.2,0.6,0.3);
        for(int i=0;i<fix_vec->size();i++)
            CV.row(i)=(*fix_vec)[i].transpose();
        viewer.data().set_points(CV,green);

    }
    return false;
  }
    callbackMouseMove(int* mode, MatrixXd* res, vector<VectorXd>* fix_vec,RowVector3f* last_mouse,bool* moved, long* sel,int *itr):
    fix_vec(fix_vec),moved(moved),sel(sel),itr(itr),last_mouse(last_mouse),mode(mode),res(res){}
private:
    vector<VectorXd>* fix_vec;
    int* mode;
    int *itr;
    bool *moved;
    long *sel;
    MatrixXd* res;
   RowVector3f* last_mouse;
};

class callbackKeyPressed{
public:
bool operator()(igl::opengl::glfw::Viewer & viewer, unsigned int key, int mod)
 {
   switch(key)
   {
     case 'U':
     case 'u':
     {
         updateViewer(viewer,*itr,*mode, R, neighbors,*verts,*half_edges,*weights, *fix,*fix_vec,*RHS,*first,*moved,*res,*dir_solver);
       break;
     }
     case ' ':
        *mode = (*mode + 1)%3;
       if((*mode==0) && (*fix).size()>0)
       {
         // Switching to deformation mode
         //compute Laplace. Only analyze pattern when fixed points changed
        getLaplace(*half_edges,*weights,*fix,*Laplace,DEFORM);
        double h = 0.01;
        double m = 0.1;
        double mu = 1.0/(h*h);
        VectorXd mass_per_vertex = VectorXd::Ones(Laplace->rows()) * m;
        //SparseMatrix<double> Aaug = *Laplace;
        //Aaug.diagonal().array() += mu * mass_per_vertex.array();
           if(*changed){
               (*dir_solver).analyzePattern(*Laplace);
               *changed=false;
           }
           (*dir_solver).factorize(*Laplace);
       }
        updateViewer(viewer,*itr,*mode, R, neighbors,*verts,*half_edges, *weights, *fix,*fix_vec,*RHS,*first,*moved,*res,*dir_solver);
       break;
       case 'G':
       case 'g':
           printObjModel(*res,*F,output_name,*num);
           *num=(*num)+1;
           break;
     default:
       return false;
   }
   return true;
 }
    callbackKeyPressed(int* mode,
                       int* itr,
                       int* num,
                      MatrixXd* res,
                      vector<int>* fix,
                      vector<VectorXd>* fix_vec,
                      MatrixXd* R,
                      vector<int>* neighbors,
                      MatrixXd* verts,
                      vector<HalfEdge>* half_edges,
                      VectorXd* weights,
                      MatrixXd* RHS,
                      MatrixXi* F,
                      bool* first,
                       bool* moved,
                       bool* changed,
                      SparseMatrix<double>* Laplace,
                      SimplicialLDLT<SparseMatrix<double>>* dir_solver,
                       string output_name
                      ):mode(mode),itr(itr),res(res),fix(fix),fix_vec(fix_vec),R(R),neighbors(neighbors),verts(verts),half_edges(half_edges),weights(weights),RHS(RHS),F(F),first(first),Laplace(Laplace),dir_solver(dir_solver),moved(moved),changed(changed),output_name(output_name),num(num){}
private:
    int* mode;
    int *itr;
    int *num;
    MatrixXd* res;
    vector<int>* fix;
    vector<VectorXd>* fix_vec;
    MatrixXd* R;
    vector<int>* neighbors;
    MatrixXd* verts;
    vector<HalfEdge>* half_edges;
    VectorXd* weights;
    MatrixXd* RHS;
    MatrixXi* F;
    bool* first;
    bool* moved;
    bool* changed;
    SparseMatrix<double> *Laplace;
    SimplicialLDLT<SparseMatrix<double>>* dir_solver;
    string output_name;
};

class callbackPreDraw{
public:
bool operator()(igl::opengl::glfw::Viewer &viewer)
{
    if(viewer.core().is_animating&&(*mode==0)){
        if(*pause){
            if(!(*moved))  updateViewer(viewer,*itr,*mode, R, neighbors,*verts,*half_edges, *weights, *fix,*fix_vec,*RHS,*first,*moved,*res,*dir_solver);
        }
        else if(!(*inf_itr)){
            if(*itr<4)  updateViewer(viewer,*itr,*mode, R, neighbors,*verts,*half_edges, *weights, *fix,*fix_vec,*RHS,*first,*moved,*res,*dir_solver);
        }
        else  updateViewer(viewer,*itr,*mode, R, neighbors,*verts,*half_edges, *weights, *fix,*fix_vec,*RHS,*first,*moved,*res,*dir_solver);
    }
    return false;
  }
    callbackPreDraw(int* mode,
                    int* itr,
                    MatrixXd* res,
                    vector<int>* fix,
                    vector<VectorXd>* fix_vec,
                    MatrixXd* R,
                    vector<int>* neighbors,
                    MatrixXd* verts,
                    vector<HalfEdge>* half_edges,
                    VectorXd* weights,
                    MatrixXd* RHS,
                    MatrixXi* F,
                    bool* first,
                    bool* moved,
                    bool* pause,
                    bool* inf_itr,
                    SimplicialLDLT<SparseMatrix<double>>* dir_solver
                    ):mode(mode),itr(itr),res(res),fix(fix),fix_vec(fix_vec),R(R),neighbors(neighbors),verts(verts),half_edges(half_edges),weights(weights),RHS(RHS),F(F),first(first),moved(moved),pause(pause),inf_itr(inf_itr),dir_solver(dir_solver){}
private:
    int* mode;
    int* itr;
    MatrixXd* res;
    vector<int>* fix;
    vector<VectorXd>* fix_vec;
    MatrixXd* R;
    vector<int>* neighbors;
    MatrixXd* verts;
    vector<HalfEdge>* half_edges;
    VectorXd* weights;
    MatrixXd* RHS;
    MatrixXi* F;
    bool* first;
    bool* moved;
    bool* pause;
    bool* inf_itr;
    SimplicialLDLT<SparseMatrix<double>>* dir_solver;
};

#endif /* call_back_h */
