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
#include "reduced_energy.h"
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
    int method=ASAP;
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
    double g = 5.0;
    double mass = 1e-4;
    double dt = 0.1;
    string slamda="0";
    processArgv(argc,argv,input_name,itrs,method,flip_avoid,print_txtfile,print_vtkfile,print_pic,print_each_frame,pause,inf_itr,show_texture,texture_name,lambda,slamda);
    string output_name=genOutputName(input_name,method,slamda,flip_avoid);
    igl::opengl::glfw::Viewer viewer;
    vector<LocalGlobalEnergy*> energy_vec;
    double resolutions[7] = { 0.1, 0.2, 0.4, 0.6, 0.8};
    double x_step = 1;
    double y_step = -1;
    Vector3d offsets[7] = {
            Vector3d(x_step,0,0),  Vector3d(2*x_step,0,0),
            Vector3d(0,y_step,0), Vector3d(x_step,y_step,0), Vector3d(2*x_step,y_step,0),
    };
    method=ARAP;
   energy_vec.push_back(new LocalGlobalEnergy(input_name,  method, lambda, mass, g, dt, Vector3d(0,0,0)));
    for(int i=0;i<5;i++){
        energy_vec.push_back(new ReducedLocalGlobalEnergy(input_name,  method, lambda, mass, g, dt, offsets[i],resolutions[i],false));
    }
    printf("vert num=%d",energy_vec[0]->vert_size());
    global_res.resize(6*energy_vec[0]->vert_size(),3);
    global_face.resize(6*energy_vec[0]->get_faces().rows(),3);
    int verts_num = energy_vec[0]->vert_size();
    int faces_num = energy_vec[0]->get_faces().rows();
    for(int i=0;i<energy_vec.size();i++){
        global_res.middleRows(i*verts_num, verts_num)=energy_vec[i]->get_res();
        MatrixXi faces = energy_vec[i]->get_faces();
        for(int fi=0;fi<faces.rows(); fi++)
            for(int fj=0; fj<faces.cols(); fj++) faces(fi,fj) += i*verts_num;
        global_face.middleRows(i*faces_num, faces_num)=faces;
    }
    InteractiveHelper helper;
    helper.pause=pause;
    helper.inf_itr=inf_itr;
    callbackMouseDown mouseDown(energy_vec, &helper);
    callbackMouseMove mouseMove(energy_vec, &helper);
    callbackKeyPressed keyPressed(energy_vec, &helper);
    callbackPreDraw preDraw(energy_vec, &helper);
    viewer.callback_mouse_down=mouseDown;
    viewer.callback_mouse_move=mouseMove;
    viewer.callback_key_pressed=keyPressed;
    viewer.callback_pre_draw=preDraw;
    viewer.callback_mouse_up = [&](igl::opengl::glfw::Viewer&, int, int)->bool{ helper.selected_anchor = -1; return false;};
    viewer.data().set_mesh(global_res, global_face);
    viewer.data().show_lines = false;
    viewer.core().is_animating = true;
    viewer.data().face_based = true;
    viewer.data().set_vertices(global_res);
    viewer.data().set_colors(RowVector3d(1.0,0.9,0.2));
    viewer.launch(false,"deformation",256,256);
    return EXIT_SUCCESS;
}

 
