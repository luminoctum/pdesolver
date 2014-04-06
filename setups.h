#ifndef SETUPS
#define SETUPS
#include "Pooma/Arrays.h"
#include "netcdf.hh"
#include "utils.h"
#include "jansson.h"
#include <map>

namespace setups{
    struct StrCmp{
        bool operator()(const char *a, const char *b){
            printf("%s %s\n", a, b);
            printf("%d\n", strcmp(a, b));
            return strcmp(a, b) < 0;
        }
    };
    int nx, ny, frame;
    long current;
    double xlen, ylen, dx, dy, start, end, dt;
    std::string ncfile = "dynamics.nc";
    std::map<std::string, Array<2, double> > ncvar;
    Interval<2> cij, sij, sicj, cisj;

    void ncConfigure(){
        NcFile dataFile(ncfile.c_str(), NcFile::ReadOnly);
        current = dataFile.get_dim("time")->size() - 1;
        if (!dataFile.is_valid()){ ASSERT_FILE_NOT_FOUND(ncfile); }
        Array<2, double> buffer;
        for (int i = 0; i < dataFile.num_vars(); i++){
            NcVar *data = dataFile.get_var(i);
            long *edges = data->edges();
            switch (data->num_dims()){
                case 1:
                    buffer.initialize(1, edges[0]);
                    data->get(&buffer(0, 0), edges[0]);
                    break;
                case 2:
                    buffer.initialize(edges[1], edges[0]);
                    data->get(&buffer(0, 0), edges[0], edges[1]);
                    break;
                case 3:
                    buffer.initialize(edges[2], edges[1]);
                    data->set_cur(current, 0, 0);
                    data->get(&buffer(0, 0), 1, edges[1], edges[2]);
                    break;
            }
            ncvar[data->name()].initialize(buffer.domain());
            ncvar[data->name()] = buffer;
        }
        nx      = dataFile.get_att("num_grids_in_x")->as_int(0);
        ny      = dataFile.get_att("num_grids_in_y")->as_int(0);
        xlen    = dataFile.get_att("length_x")->as_double(0);
        ylen    = dataFile.get_att("length_y")->as_double(0);
        start   = ncvar["time_start"](0, current);
        dt      = dataFile.get_att("time_step")->as_double(0);
        end     = dataFile.get_att("time_end")->as_double(0);
        frame   = dataFile.get_att("steps_per_frame")->as_int(0);
        dx      = xlen / nx;
        dy      = ylen / ny;
        cij     = Interval<2>(Interval<1>(0, nx - 1), Interval<1>(0, ny - 1));
        sij[0]  = Interval<1>(cij[0].first() - 1, cij[0].last() + 1);
        sij[1]  = Interval<1>(cij[1].first() - 1, cij[1].last() + 1);
        sicj    = Interval<2>(sij[0], cij[1]);
        cisj    = Interval<2>(cij[0], sij[1]);
    }
}

#endif
