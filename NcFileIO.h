#ifndef NCFILEIO
#define NCFILEIO
#include "Pooma/Arrays.h"
#include <map>
#include "netcdf.hh"

class NcFileIO{
    typedef double Element_t;
    typedef Array<2, double> Array_t;
protected:
    std::string filename;
    std::map<std::string, Array_t> ncvar;
    Array_t buffer_in;
    long current;
    int nx, ny, gl;
    double xlen, ylen;
    double dx, dy;
    double start, end, dt;
    int frame;
public:
    NcFileIO(){}
    NcFileIO(std::string _fname){
        filename = _fname;
        ncread(_fname);
    }
    friend std::ostream& operator<< (std::ostream &os, const NcFileIO &other){
        os << "========== Experiment file name : " << other.filename << " ==========" << std::endl
            << "Number of Grids in X: " << other.nx << std::endl
            << "Number of Grids in Y: " << other.ny << std::endl
            << "Total length in X: " << other.xlen / 1000. << " km" << std::endl
            << "Total length in Y: " << other.ylen / 1000. << " km" << std::endl
            << "Grid size in X: " << other.dx / 1000. << " km" << std::endl
            << "Grid size in Y: " << other.dy / 1000. << " km" <<std::endl
            << "Time start: " << other.start << " s" << std::endl
            << "Time end: " << other.end << " s" << std::endl
            << "Time step: " << other.dt<< " s" << std::endl
            << "Times per frame: " << other.frame << std::endl;
        return os;
    }
    void ncread(std::string _fname){
        NcFile dataFile(_fname.c_str(), NcFile::ReadOnly);
        current = dataFile.get_dim("time")->size() - 1;
        if (!dataFile.is_valid()){ ASSERT_FILE_NOT_FOUND(_fname); }
        Array_t buffer_in;
        for (int i = 0; i < dataFile.num_vars(); i++){
            NcVar *data = dataFile.get_var(i);
            long *edges = data->edges();
            switch (data->num_dims()){
                case 1:
                    buffer_in.initialize(1, edges[0]);
                    data->get(&buffer_in(0, 0), edges[0]);
                    break;
                case 2:
                    buffer_in.initialize(edges[1], edges[0]);
                    data->get(&buffer_in(0, 0), edges[0], edges[1]);
                    break;
                case 3:
                    buffer_in.initialize(edges[2], edges[1]);
                    data->set_cur(current, 0, 0);
                    data->get(&buffer_in(0, 0), 1, edges[1], edges[2]);
                    break;
            }
            ncvar[data->name()].initialize(buffer_in.domain());
            ncvar[data->name()] = buffer_in;
        }
        nx      = dataFile.get_att("nx")->as_int(0);
        ny      = dataFile.get_att("ny")->as_int(0);
        gl      = dataFile.get_att("gl")->as_int(0);
        xlen    = dataFile.get_att("xlen")->as_double(0);
        ylen    = dataFile.get_att("ylen")->as_double(0);
        start   = ncvar["time"](0, current);
        end     = dataFile.get_att("end")->as_double(0);
        dt      = dataFile.get_att("step")->as_double(0);
        frame   = dataFile.get_att("frame")->as_int(0);
        dx      = xlen / (nx - 1);
        dy      = ylen / (ny - 1);
    }
    /*
    template<class A, class B>
    void ncwrite(const A& var, const B& vattr, double time){
        current++;
        NcFile dataFile(filename.c_str(), NcFile::Write);
        for (size_t i = 0; i < var.size(); i++){
            buffer_out = var[i](cij);
            dataFile.get_var(vattr[i].name.c_str())->put_rec(&buffer_out(0, 0), current);
        }
        dataFile.get_var("time")->put_rec(&time, current);
    }*/
};

#endif
