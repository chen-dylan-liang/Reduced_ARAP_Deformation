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
#include "energy.h"
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
    double lambda=0;
    double g=1.0;
    double mass=0.001;
    string slamda="0";
    processArgv(argc,argv,input_name,itrs,method,flip_avoid,print_txtfile,print_vtkfile,print_pic,print_each_frame,pause,inf_itr,show_texture,texture_name,lamda,slamda);

    MatrixXd res;

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

    for(int i=0;i<vert_num;i++) res.row(i)=verts.col(i).transpose();
    LocalGlobalEnergy e(input_name, mass, method, lambda, g);
    //set viewer. The code of initial guess and ARAP optimizations for deformation is in the file call_back.h
    callbackMouseDown mouseDown(&e, &mode,&now_itr, &first, &moved,&changed,&sel, &last_mouse, &res);
    callbackMouseMove mouseMove(&e, &mode, &now_itr, &last_mouse, &res, &last_mouse,&moved,&sel);
    callbackKeyPressed keyPressed(&e, &mode,&now_itr,&print_num,&res,&first,&moved,&changed, output_name);
    callbackPreDraw preDraw(&e, &mode,&now_itr,&res,&first,&moved,&pause,&inf_itr);
    viewer.callback_mouse_down=mouseDown;
    viewer.callback_mouse_move=mouseMove;
    viewer.callback_key_pressed=keyPressed;
    viewer.callback_pre_draw=preDraw;
    viewer.callback_mouse_up = [&](igl::opengl::glfw::Viewer&, int, int)->bool{ sel = -1; return false;};
    viewer.data().set_mesh(res,e.faces());
    viewer.data().show_lines = false;
    viewer.core().is_animating = true;
    viewer.data().face_based = true;
    viewer.data().set_vertices(res);
    viewer.data().set_colors(RowVector3d(1.0,0.9,0.2));
    viewer.launch(false,"deformation",256,256);
    return EXIT_SUCCESS;
}

 
