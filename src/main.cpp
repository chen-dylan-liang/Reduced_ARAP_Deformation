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
    processArgv(argc,argv,input_name,itrs,method,flip_avoid,print_txtfile,print_vtkfile,print_pic,print_each_frame,pause,inf_itr,show_texture,texture_name,lambda,slamda);
    string output_name=genOutputName(input_name,method,slamda,flip_avoid);
    igl::opengl::glfw::Viewer viewer;
    LocalGlobalEnergy e(input_name, mass, method, lambda, g);
    InteractiveHelper helper;
    helper.pause=pause;
    helper.inf_itr=inf_itr;
    callbackMouseDown mouseDown(&e, &helper);
    callbackMouseMove mouseMove(&e, &helper);
    callbackKeyPressed keyPressed(&e, &helper);
    callbackPreDraw preDraw(&e, &helper);
    viewer.callback_mouse_down=mouseDown;
    viewer.callback_mouse_move=mouseMove;
    viewer.callback_key_pressed=keyPressed;
    viewer.callback_pre_draw=preDraw;
    viewer.callback_mouse_up = [&](igl::opengl::glfw::Viewer&, int, int)->bool{ sel = -1; return false;};
    viewer.data().set_mesh(e.get_res(),e.get_faces());
    viewer.data().show_lines = false;
    viewer.core().is_animating = true;
    viewer.data().face_based = true;
    viewer.data().set_vertices(e.get_res());
    viewer.data().set_colors(RowVector3d(1.0,0.9,0.2));
    viewer.launch(false,"deformation",256,256);
    return EXIT_SUCCESS;
}

 
