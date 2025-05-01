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
#include "energy.h"

inline static MatrixXd global_res;
inline static MatrixXi global_face;

struct InteractiveHelper{
    // mode == 0: anchor selection
    // mode == 1: move anchors + simulation
    int mode=0;
    // the position of mouse last time
    RowVector3f last_mouse=RowVector3f::Zero();
    bool anchor_changed=false;
    bool anchor_moved=false;
    long selected_anchor=-1;
    int itr=0;
    bool inf_itr=false;
    bool pause=false;
    string output_name;
    int output_num=0;

};

inline void updateViewer(igl::opengl::glfw::Viewer& viewer, vector<LocalGlobalEnergy*>& energy, InteractiveHelper* helper)
 {
   // predefined colors
   const Eigen::RowVector3d orange(1.0,0.7,0.2);
   const Eigen::RowVector3d yellow(1.0,0.9,0.2);
   const Eigen::RowVector3d blue(0.2,0.3,0.8);
   const Eigen::RowVector3d green(0.2,0.6,0.3);
   int verts_num = energy[0]->vert_size();
   // selection mode
   if(helper->mode==0 || helper->mode==2){
      for(int i=0;i<energy.size();i++){
          global_res.middleRows(i*verts_num, verts_num)=energy[i]->get_res();
      }
      viewer.data().set_vertices(global_res);
      if (helper->mode ==0) viewer.data().set_colors(yellow);
      else  viewer.data().set_colors(orange);
      MatrixXd CV(energy[0]->get_anchors().size(),3);
      for(int i=0;i<energy[0]->get_anchors().size();i++) CV.row(i)=energy[0]->get_anchor_points()[i].transpose();
      viewer.data().set_points(CV,green);
   }
   // solve
   else{

           //local and global optimizations
           for(int i=0;i<energy.size();i++){
               energy[i]->local_global_solve();
               global_res.middleRows(i*verts_num, verts_num)=energy[i]->get_res();
           }
           helper->itr++;
     viewer.data().set_vertices(global_res);
     viewer.data().set_colors(blue);
     // render anchor points
     MatrixXd CU(energy[0]->get_anchors().size(),3);
    for(int i=0;i<energy[0]->get_anchors().size();i++)
    CU.row(i)=energy[0]->get_anchor_points()[i].transpose();
    viewer.data().set_points(CU,orange);
   }
   viewer.data().compute_normals();
 }


class callbackMouseDown{
    public:
        bool operator() (igl::opengl::glfw::Viewer& viewer, int, int)
    {
        helper->last_mouse = Eigen::RowVector3f(viewer.current_mouse_x,viewer.core().viewport(3)-viewer.current_mouse_y,0);
        if(helper->mode==0){
          // Find the closest point on mesh to mouse position
          int fid;
          Eigen::Vector3f bary;
          if(igl::unproject_onto_mesh((helper->last_mouse).head(2),
                                      viewer.core().view,
                                      viewer.core().proj,
                                      viewer.core().viewport,
                                      energy[0]->get_res(),
                                      energy[0]->get_faces(),
                                      fid,
                                      bary)){
            long c;
            bary.maxCoeff(&c);
            Eigen::RowVector3d new_c = (energy[0]->get_res()).row((energy[0]->get_faces())(fid,c));
            MatrixXd CV((energy[0]->get_anchors()).size(),3);
            for(int i=0;i<(energy[0]->get_anchors()).size();i++) CV.row(i)=((energy[0]->get_anchor_points())[i]).transpose();
            if(energy[0]->get_anchors().empty() || (CV.rowwise()-new_c).rowwise().norm().minCoeff() > 0)
            {
                for(int i=0;i<energy.size();i++){
                    Eigen::RowVector3d new_c_i = (energy[i]->get_res()).row((energy[i]->get_faces())(fid,c));
                    energy[i]->add_anchor((energy[i]->get_faces())(fid,c), new_c_i.transpose());
                }
                helper->anchor_changed=true;
                updateViewer(viewer,energy, helper);
              return true;
            }
          }
        }
        else if(helper->mode==1){
          // Move the closest control point
          Eigen::MatrixXf CP;
          MatrixXd CU((energy[0]->get_anchors()).size(),3);
          for(int i=0;i<(energy[0]->get_anchors()).size();i++) CU.row(i)= (energy[0]->get_anchor_points())[i].transpose();
          igl::project(Eigen::MatrixXf(CU.cast<float>()),
                       viewer.core().view,
                       viewer.core().proj,
                       viewer.core().viewport,
                       CP);
          Eigen::VectorXf D = (CP.rowwise()-(helper->last_mouse)).rowwise().norm();
          helper->selected_anchor = (D.minCoeff(&helper->selected_anchor) < 30)? (helper->selected_anchor):-1;
          if(helper->selected_anchor != -1)
          {
            (helper->last_mouse)(2) = CP(helper->selected_anchor,2);
            return true;
          }
        }
        return false;
    }

    callbackMouseDown(vector<LocalGlobalEnergy*> energy,InteractiveHelper* helper):
                                        energy(energy),helper(helper){}
    private:
        vector<LocalGlobalEnergy*> energy;
        InteractiveHelper* helper;
};

class callbackMouseMove{
    public:
        bool operator() (igl::opengl::glfw::Viewer & viewer, int, int)
      {
            // simulation mode
        if(helper->mode==1){
            if(helper->selected_anchor!=-1)
            {
              Eigen::RowVector3f drag_mouse(viewer.current_mouse_x,viewer.core().viewport(3) - viewer.current_mouse_y,(helper->last_mouse)(2));
              Eigen::RowVector3f drag_scene,last_scene;
              igl::unproject(drag_mouse,viewer.core().view,viewer.core().proj,viewer.core().viewport,drag_scene);
              igl::unproject(helper->last_mouse,viewer.core().view,viewer.core().proj,viewer.core().viewport,last_scene);
              if(!energy[0]->anchor_in_region(helper->selected_anchor)){
                for(int i=0; i<energy.size(); i++){
                    (energy[i]->get_anchor_points())[helper->selected_anchor] +=(drag_scene-last_scene).cast<double>();
                }
              }
              else{
                  for(int i=0; i<energy.size(); i++){
                      for(int j=0; j<energy[i]->get_anchors().size();j++)
                        if(energy[i]->anchor_in_region(j)){
                          (energy[i]->get_anchor_points())[j] += (drag_scene-last_scene).cast<double>();
                      }
                  }
              }
              helper->last_mouse = drag_mouse;
              helper->itr=0;
              helper->anchor_moved=true;
              return true;
            }
        }
        // region selection
        else if (helper->mode==2){
            Eigen::RowVector3f drag_mouse(viewer.current_mouse_x, viewer.core().viewport(3) - viewer.current_mouse_y, (helper->last_mouse)(2));
            Eigen::RowVector3f start(std::min(drag_mouse.x(), helper->last_mouse.x()),std::min(drag_mouse.y(), helper->last_mouse.y()), 0.0);
            Eigen::RowVector3f end(std::max(drag_mouse.x(),helper->last_mouse.x()),std::max(drag_mouse.y(), helper->last_mouse.y()), 0.0);
            for(int i=0; i<energy[0]->vert_size();i++){
                Eigen::RowVector3f projected=igl::project(Eigen::Matrix<float, 3,1>(energy[0]->get_res().row(i).transpose().cast<float>()),
                        viewer.core().view,
                        viewer.core().proj,
                        viewer.core().viewport);
                if(projected.x()>start.x()&&projected.x()<end.x()&&projected.y()>start.y()&&projected.y()<end.y()){
                    for(int j=0; j<energy.size();j++) energy[j]->add_anchor(i, energy[j]->get_res().row(i).transpose(),true);
                }
            }
            updateViewer(viewer, energy,helper);
        }
        return false;
      }
        callbackMouseMove(vector<LocalGlobalEnergy*> energy,InteractiveHelper* helper):
                energy(energy),helper(helper){}
    private:
    vector<LocalGlobalEnergy*> energy;
        InteractiveHelper* helper;
};

class callbackKeyPressed{
    public:
    bool operator()(igl::opengl::glfw::Viewer & viewer, unsigned int key, int mod)
     {
       switch(key)
       {
           // anchor selection mode
           case '1':
               helper->mode =0;
               updateViewer(viewer,energy, helper);
               break;
           // region selection mode
           case '2':
               helper->mode =2;
               updateViewer(viewer,energy, helper);
               break;
           // remove anchors
           case 'R':
           case 'r':
               if(helper->mode!=1){
               for(int i=0;i<energy.size();i++) energy[i]->clear_anchors();
               updateViewer(viewer,energy, helper);}

               break;
         // simulation mode, when already in simulation mode, this alternates subspace
         case ' ':
            helper->mode = 1;
           if(energy[0]->get_anchors().size()>0)
           {
             //compute Laplace. Only analyze pattern when fixed points changed
             for(int i=0;i<energy.size();i++) energy[i]->compute_laplacian();
             if(helper->anchor_changed){
                   helper->anchor_changed=false;
               }
               for(int i=0;i<energy.size();i++) energy[i]->solver_compute();
           }
            updateViewer(viewer,energy, helper);
           break;
           case 'G':
           case 'g':
                printObjModel(energy[0]->get_res(),energy[0]->get_faces(),helper->output_name,helper->output_num);
                helper->output_num++;
               break;
         default:
           return false;
       }
       return true;
     }
    callbackKeyPressed(vector<LocalGlobalEnergy*> energy,InteractiveHelper* helper):
            energy(energy),helper(helper){}
    private:
        vector<LocalGlobalEnergy*> energy;
        InteractiveHelper* helper;
};

class callbackPreDraw{
    public:
    bool operator()(igl::opengl::glfw::Viewer &viewer)
    {
        if(viewer.core().is_animating&&(helper->mode==1)){
            if(helper->pause){
                if(!(helper->anchor_moved))  updateViewer(viewer,energy, helper);
            }
            else if(!(helper->inf_itr)){
                if(helper->itr<4)  updateViewer(viewer,energy, helper);
            }
            else  updateViewer(viewer,energy, helper);
        }
        return false;
      }
        callbackPreDraw(vector<LocalGlobalEnergy*> energy,InteractiveHelper* helper):
                energy(energy),helper(helper){}
    private:
        vector<LocalGlobalEnergy*> energy;
        InteractiveHelper* helper;
};

#endif /* call_back_h */
