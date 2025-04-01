#ifndef CLASS_MESH_HPP
  #define CLASS_MESH_HPP
  #include <vector>

  class class_mesh {
    private:
      
    public:
      std::vector<double> x; // x1, x2, x3, ...
      std::vector<double> y; // y1, y2, y3, ...
          
      std::vector<int> face_point1; // f1p1, f2p1, ...
      std::vector<int> face_point2; // f1p2, f2p2, ...

      std::vector<int> face_cell1; // f1c1, f2c1, ...
      std::vector<int> face_cell2; // f1c2, f2c2, ...

      std::vector<double> face_area; // f1a, f2a, f3a, ...
      
      std::vector<double> face_nx; // f1nx, f2nx, f3nx, ...
      std::vector<double> face_ny; // f1ny, f2ny, f3ny, ...
      std::vector<int> face_cell1_n_out; // if normal vector is point out of cell -> 1 else -> -1
      std::vector<int> face_cell2_n_out; // if normal vector is point out of cell -> 1 else -> -1

      std::vector<double> face_cell1_rx; //f1rx, f2rx, ...
      std::vector<double> face_cell2_rx; //f1rx, f2rx, ...
      std::vector<double> face_cell1_ry; //f1ry, f2ry, ...
      std::vector<double> face_cell2_ry; //f1ry, f2ry, ...

      std::vector<double> face_cell1_rx_norm; //f1rx, f2rx, ...
      std::vector<double> face_cell1_ry_norm; //f1ry, f2ry, ...
      std::vector<double> face_cell2_rx_norm; //f1rx, f2rx, ...
      std::vector<double> face_cell2_ry_norm; //f1ry, f2ry, ...

      std::vector<double> face_centerx; // f1x, f2x, f3x, ...
      std::vector<double> face_centery; // f1y, f2y, f3y, ...

      std::vector<double> cell_vol; // c1v, c2v, c3v, ...
      std::vector<double> cell_centerx; // c1x, c2x, c3x, ...
      std::vector<double> cell_centery; // c1y, c2y, c3y, ...


      std::vector<int> cell_faces; // c1f1, c1f2, c1f3, c1f4, ..., c2f1, c2f2, c2f3, c2f4, ... if no face 4, 5, or 6 = -2 
      std::vector<int> cell_face_count; // c1fc, c2fc, c3fc, ...

      std::vector<double> cell_face_rx;
      std::vector<double> cell_face_ry;
      std::vector<double> cell_face_rx_norm;
      std::vector<double> cell_face_ry_norm;
      std::vector<double> cell_face_n_out;

      void init(int size) {
        face_cell1_n_out.assign(size,-101);
        face_cell2_n_out.assign(size,-101);

        face_cell1_rx.assign(size,-101);
        face_cell2_rx.assign(size,-101);
        face_cell1_ry.assign(size,-101);
        face_cell2_ry.assign(size,-101);

        face_cell1_rx_norm.assign(size,-101);
        face_cell1_ry_norm.assign(size,-101);
        face_cell2_rx_norm.assign(size,-101);
        face_cell2_ry_norm.assign(size,-101);
      }

      int find_cell_face(int cell_num, int cell_face) {
        return cell_faces[cell_num*6 + cell_face];
      }

      double get_cell_face_rx(int cell_num, int cell_face) {
        return cell_face_rx[cell_num*6 + cell_face];
      }

      // std::vector<int> face_cell1_Qface_num; // Face number of Cell Qi that corresponds to 
      // std::vector<int> face_cell2_Qface_num;

      // std::vector<double> cell_face_rx; //c1f1r, c1f2r, c1f3r, c1f4r, c1f5r, c1f6r, c2f1r, ...
      // std::vector<double> cell_face_ry; //c1f1r, c1f2r, c1f3r, c1f4r, c1f5r, c1f6r, c2f1r, ...

      // std::vector<double> cell_face_rx_norm; //c1f1r, c1f2r, c1f3r, c1f4r, c1f5r, c1f6r, c2f1r, ...
      // std::vector<double> cell_face_ry_norm; //c1f1r, c1f2r, c1f3r, c1f4r, c1f5r, c1f6r, c2f1r, ...

      // int find(int cell_num, int cell_face) {
      //   return cell_num*6 + cell_face;
      // }

      // double find_cell_face_pointx(int cell_num, int cell_face, int point_num) {
      //   return x[face_points[2*cell_faces[cell_num*6 + cell_face]+point_num]];
      // }

      // double find_cell_face_pointy(int cell_num, int cell_face, int point_num) {
      //   return y[face_points[2*cell_faces[cell_num*6 + cell_face]+point_num]];
      // }

      // double find_cell_face_area(int cell_num, int cell_face) {
      //   return face_area[cell_faces[cell_num*6 + cell_face]];
      // }

      // double find_cell_face_center_x(int cell_num, int cell_face) {
      //   return face_centerx[cell_faces[cell_num*6 + cell_face]];
      // }

      // double find_cell_face_center_y(int cell_num, int cell_face) {
      //   return face_centery[cell_faces[cell_num*6 + cell_face]];
      // }

      // double find_cell_face_nx(int cell_num, int cell_face) {
      //   return face_nx[cell_faces[cell_num*6 + cell_face]];
      // }

      // double find_cell_face_ny(int cell_num, int cell_face) {
      //   return face_ny[cell_faces[cell_num*6 + cell_face]];
      // }

      // double get_cell_face_n_out(int cell_num, int cell_face){
      //   return cell_face_n_out[cell_num*6 + cell_face];
      // }
  };
#endif