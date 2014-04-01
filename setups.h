#ifndef SETUPS
#define SETUPS
#include "Pooma/Arrays.h"
#include "netcdf.hh"
#include "utils.h"
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
    Interval<2> cij;

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
        nx      = dataFile.get_att("nx")->as_int(0);
        ny      = dataFile.get_att("ny")->as_int(0);
        xlen    = dataFile.get_att("xlen")->as_double(0);
        ylen    = dataFile.get_att("ylen")->as_double(0);
        start   = ncvar["time"](0, current);
        end     = dataFile.get_att("end")->as_double(0);
        dt      = dataFile.get_att("step")->as_double(0);
        frame   = dataFile.get_att("frame")->as_int(0);
        dx      = xlen / (nx - 1);
        dy      = ylen / (ny - 1);
        cij     = Interval<2>(Interval<1>(0, nx - 1), Interval<1>(0, ny - 1));
    }
}

#endif
