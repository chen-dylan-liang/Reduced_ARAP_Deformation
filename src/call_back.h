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

inline void updateViewer(igl::opengl::glfw::Viewer& viewer, LocalGlobalEnergy& energy, InteractiveHelper& helper)
 {
   // predefined colors
   const Eigen::RowVector3d orange(1.0,0.7,0.2);
   const Eigen::RowVector3d yellow(1.0,0.9,0.2);
   const Eigen::RowVector3d blue(0.2,0.3,0.8);
   const Eigen::RowVector3d green(0.2,0.6,0.3);
   // anchor selection
   if(helper.mode==0){
      viewer.data().set_vertices(energy.get_res());
      viewer.data().set_colors(yellow);
      MatrixXd CV(energy.get_anchors().size(),3);
      for(int i=0;i<energy.get_anchors().size();i++) CV.row(i)=energy.get_anchor_points()[i].transpose();
      viewer.data().set_points(CV,green);
   }
   // solve
   else{
       if(helper.anchor_moved){
           //local and global optimizations
           energy.local_global_solve();
           printf("solved!\n");
           for(int i=0;i<energy.get_anchors().size();i++) energy.get_anchor_points()[i]=(energy.get_res().row(energy.get_anchors()[i]).transpose());
           helper.itr++;
       }
     viewer.data().set_vertices(energy.get_res());
     viewer.data().set_colors(blue);
     // render anchor points
     MatrixXd CU(energy.get_anchors().size(),3);
    for(int i=0;i<energy.get_anchors().size();i++)
    CU.row(i)=energy.get_anchor_points()[i].transpose();
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
                                      energy->get_res(),
                                      energy->get_faces(),
                                      fid,
                                      bary)){
            long c;
            bary.maxCoeff(&c);
            Eigen::RowVector3d new_c = (energy->get_res()).row((energy->get_faces())(fid,c));
            MatrixXd CV((energy->get_anchors()).size(),3);
            for(int i=0;i<(energy->get_anchors()).size();i++) CV.row(i)=((energy->get_anchor_points())[i]).transpose();
            if(energy->get_anchors().empty() || (CV.rowwise()-new_c).rowwise().norm().minCoeff() > 0)
            {
                energy->add_anchor((energy->get_faces())(fid,c), new_c.transpose());
                helper->anchor_changed=true;
                updateViewer(viewer,*energy, *helper);
              return true;
            }
          }
        }
        else if(helper->mode==1){
          // Move the closest control point
          Eigen::MatrixXf CP;
          MatrixXd CU((energy->get_anchors()).size(),3);
          for(int i=0;i<(energy->get_anchors()).size();i++) CU.row(i)= (energy->get_anchor_points())[i].transpose();
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

    callbackMouseDown(LocalGlobalEnergy* energy,InteractiveHelper* helper):
                                        energy(energy),helper(helper){}
    private:
        LocalGlobalEnergy* energy;
        InteractiveHelper* helper;
};

class callbackMouseMove{
    public:
        bool operator() (igl::opengl::glfw::Viewer & viewer, int, int)
      {
        if(helper->selected_anchor!=-1)
        {
          Eigen::RowVector3f drag_mouse(viewer.current_mouse_x,viewer.core().viewport(3) - viewer.current_mouse_y,(helper->last_mouse)(2));
          Eigen::RowVector3f drag_scene,last_scene;
          igl::unproject(drag_mouse,viewer.core().view,viewer.core().proj,viewer.core().viewport,drag_scene);
          igl::unproject(helper->last_mouse,viewer.core().view,viewer.core().proj,viewer.core().viewport,last_scene);
          (energy->get_anchor_points())[helper->selected_anchor] +=(drag_scene-last_scene).cast<double>();
          helper->last_mouse = drag_mouse;
          helper->itr=0;
          helper->anchor_moved=true;
          return true;
        }
        return false;
      }
        callbackMouseMove(LocalGlobalEnergy* energy,InteractiveHelper* helper):
                energy(energy),helper(helper){}
    private:
        LocalGlobalEnergy* energy;
        InteractiveHelper* helper;
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
             updateViewer(viewer,*energy, *helper);
           break;
         }
         // switch mode
         case ' ':
            helper->mode = (helper->mode + 1)%2;
           if((helper->mode==1) && (energy->get_anchors()).size()>0)
           {
             //compute Laplace. Only analyze pattern when fixed points changed
             energy->compute_laplacian();
             if(helper->anchor_changed){
                   helper->anchor_changed=false;
               }
               energy->solver_compute();
           }
            updateViewer(viewer,*energy, *helper);
           break;
           case 'G':
           case 'g':
                printObjModel(energy->get_res(),energy->get_faces(),helper->output_name,helper->output_num);
                helper->output_num++;
               break;
         default:
           return false;
       }
       return true;
     }
    callbackKeyPressed(LocalGlobalEnergy* energy,InteractiveHelper* helper):
            energy(energy),helper(helper){}
    private:
        LocalGlobalEnergy* energy;
        InteractiveHelper* helper;
};

class callbackPreDraw{
    public:
    bool operator()(igl::opengl::glfw::Viewer &viewer)
    {
        if(viewer.core().is_animating&&(helper->mode==1)){
            if(helper->pause){
                if(!(helper->anchor_moved))  updateViewer(viewer,*energy, *helper);
            }
            else if(!(helper->inf_itr)){
                if(helper->itr<4)  updateViewer(viewer,*energy, *helper);
            }
            else  updateViewer(viewer,*energy, *helper);
        }
        return false;
      }
        callbackPreDraw(LocalGlobalEnergy* energy,InteractiveHelper* helper):
                energy(energy),helper(helper){}
    private:
        LocalGlobalEnergy* energy;
        InteractiveHelper* helper;
};

#endif /* call_back_h */
