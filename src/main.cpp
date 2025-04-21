//
//  main.cpp
//  arap
//
//  Created by Liang Chen on 2021/5/20.
//

#include "global_var.h"
#include "helper_init.h"
#include "helper_print.h"
#include "helper_geometry.h"
#include "ARAP.h"
#include "call_back.h"
#include <Eigen/Dense>
#include <Eigen/Core>
#include <Eigen/SparseCore>
#include <Eigen/SparseCholesky>
template class Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>>;

int main(int argc, const char * argv[]) {
    /*initialization*/
    string input_name;
    int itrs=4;
    int method=ARAP;
    bool flip_avoid=false;
    bool print_pic=false;
    bool print_vtkfile=false;
    bool print_txtfile=true;
    bool print_each_frame=false;
    bool pause=false;
    bool inf_itr=false;
    bool show_texture=false;
    string texture_name;
    double lamda=0;
    string slamda="0";
    processArgv(argc,argv,input_name,itrs,method,flip_avoid,print_txtfile,print_vtkfile,print_pic,print_each_frame,pause,inf_itr,show_texture,texture_name,lamda,slamda);
    MatrixXd verts;
    MatrixXi edges;
    vector<HalfEdge> half_edges;
    VectorXd area;
    MatrixXi F;
    int vert_num=readObj(input_name,F,edges,verts,half_edges,area);
    
    MatrixXd res;
    MatrixXd RHS;
    VectorXd weights;
    SparseMatrix<double> Laplace;
    Laplace.resize(vert_num,vert_num);
    Eigen::SimplicialLDLT<SparseMatrix<double> > dir_solver;
    MatrixXd* R=NULL;
    vector<int> fix;
    vector<VectorXd> fix_vec;
    string output_name=genOutputName(input_name,method,slamda,flip_avoid);
    igl::opengl::glfw::Viewer viewer;

    /*DEFORMATION*/
    /*Preparation for deformation*/
    int mode=1;
    bool first=true;
    long sel=-1;
    int now_itr=0;
    int print_num=0;
    bool moved=false;
    bool changed=false;
    RowVector3f last_mouse(0,0,0);
    edges.resize(0,0);
    R=new MatrixXd[vert_num];
    for(int i=0;i<vert_num;i++) {
            R[i].resize(3,3);
            R[i].setZero();
            R[i](0,0)= R[i](1,1)= R[i](2,2)=1;
    }
    weights.resize(vert_num*vert_num);
    getWeights(half_edges,verts,weights,Cotangent_2);
    sort(half_edges.begin(),half_edges.end(),compare);
    vector<int>* neighbors=new vector<int>[vert_num];
    getNeighbors(half_edges,neighbors);
    res.resize(vert_num,3);
    RHS.resize(vert_num,3);
    for(int i=0;i<vert_num;i++) res.row(i)=verts.col(i).transpose();
    prev_res = res;
    //set viewer. The code of initial guess and ARAP optimizations for deformation is in the file call_back.h
    callbackMouseDown mouseDown(&mode,&now_itr,&last_mouse,&res,&fix,&fix_vec,R,neighbors,&verts,&half_edges,&weights,&RHS,&F,&first,&moved,&changed,&sel,&dir_solver);
    callbackMouseMove mouseMove(&mode, &res, &fix_vec,&last_mouse,&moved,&sel,&now_itr);
    callbackKeyPressed keyPressed(&mode,&now_itr,&print_num,&res,&fix,&fix_vec,R,neighbors,&verts,&half_edges,&weights,&RHS,&F,&first,&moved,&changed,&Laplace,&dir_solver,output_name);
    callbackPreDraw preDraw(&mode,&now_itr,&res,&fix,&fix_vec,R,neighbors,&verts,&half_edges,&weights,&RHS,&F,&first,&moved,&pause,&inf_itr,&dir_solver);
    viewer.callback_mouse_down=mouseDown;
    viewer.callback_mouse_move=mouseMove;
    viewer.callback_key_pressed=keyPressed;
    viewer.callback_pre_draw=preDraw;
    viewer.callback_mouse_up = [&](igl::opengl::glfw::Viewer&, int, int)->bool{ sel = -1; return false;};
    viewer.data().set_mesh(res,F);
    viewer.data().show_lines = false;
    viewer.core().is_animating = true;
    viewer.data().face_based = true;
    viewer.data().set_vertices(res);
    viewer.data().set_colors(RowVector3d(1.0,0.9,0.2));
    viewer.launch(false,"deformation",256,256);
    if(R) delete[] R;
    if(neighbors) delete[] neighbors;
    return EXIT_SUCCESS;
}

 
