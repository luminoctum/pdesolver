#ifndef VARIABLE
#define VARIABLE
#include "Pooma/Arrays.h"
#include "Boundary.h"

class VariableBase{
public:
    virtual ~VariableBase(){}
    virtual void printInfo() = 0;
    virtual void fixBoundary() = 0;
};

template<int N, class T, class B = Periodic>
class Variable : public VariableBase, public B
{
    typedef T Element_t;
    typedef B Boundary_t;
public:
    Variable(){}
    Variable(char *_name, Interval<2> _cij) :
        name(_name), cij(_cij){
            Interval<1> cx(cij[0].first() - N, cij[0].last() + N);
            Interval<1> cy(cij[1].first() - N, cij[1].last() + N);
            cxy = Interval<2>(cx, cy);
            cell.initialize(cxy); cell = 0;
            cell_t.initialize(cij); cell_t = 0;
            wallx.initialize(cxy); wallx = 0;
            wally.initialize(cxy); wally = 0;
        }
    void printInfo(){
        printf("Variable: %s\n", name);
        std::cout << "Domain: " << cij << " , " << cxy << std::endl;
        printf("Cell Value: \n");
        std::cout << cell;
        printf("Wall Value: \n");
        std::cout << wallx;
        std::cout << "Maximum Value: " << max(cell) << ", Minimum Value: " << min(cell) << std::endl;
        std::cout << std::endl;
    }
    void fixBoundary(){ fix(cell, cij); }
    char *name;
    Interval<2> cxy, cij;
    Array<2, Element_t> cell, cell_t;
    Array<2, Vector<2, Element_t> > wallx, wally;
};

// default template arguments may not be used in partial specializations
template<int N, int S, class T, class B>
class Variable<N, Vector<S, T>, B> : public VariableBase, public B
{
    typedef Vector<S, T> Element_t;
    typedef B Boundary_t;
public:
    Variable(){}
    Variable(char *_name[], Interval<2> _cij){
        cij = _cij;
        for (int i = 0; i < S; i++){
            name[i] = new char[strlen(_name[i])];
            strcpy(name[i], _name[i]);
        };
        Interval<1> cx(cij[0].first() - N, cij[0].last() + N);
        Interval<1> cy(cij[1].first() - N, cij[1].last() + N);
        cxy = Interval<2>(cx, cy);
        cell.initialize(cxy); cell = 0;
        cell_t.initialize(cij); cell_t = 0;
        wallx.initialize(cxy); wallx = 0;
        wally.initialize(cxy); wally = 0;
    }
    ~Variable(){ for (int i = 0; i < S; i++) delete name[i]; }
    void printInfo(){
        printf("Variable: (");
        for (int i = 0; i < S - 1; i++) printf("%s, ", name[i]);
        printf("%s)\n", name[S - 1]);
        std::cout << "Domain: " << cij << " , " << cxy << std::endl;
        printf("Cell Value: \n");
        std::cout << cell;
        Element_t vmax, vmin;
        for (int i = 0; i < S; i++){
            vmax(i) = max(cell.comp(i));
            vmin(i) = min(cell.comp(i));
        }
        std::cout << "Maximum Value: " << vmax << ", Minimum Value: " << vmin << std::endl;
        std::cout << std::endl;
    }
    void fixBoundary(){ fix(cell, cij); }
    char *name[S];
    Interval<2> cxy, cij;
    Array<2, Element_t> cell, cell_t;
    Array<2, Vector<2, Element_t> > wallx, wally;
};

#endif
