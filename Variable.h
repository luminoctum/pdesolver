#ifndef VARIABLE
#define VARIABLE
#include "Pooma/Arrays.h"
#include "Boundary.h"
#include "setups.h"

class VariableBase{
public:
    VariableBase(){ cij = setups::cij; }
    virtual ~VariableBase(){}
    virtual void printInfo(int level = 0) = 0;
    virtual void fixBoundary() = 0;
    Interval<2> cxy, cij;
};

template<int N, class T, class B = Periodic>
class Variable : public VariableBase, public B
{
    typedef T Element_t;
    typedef B Boundary_t;
public:
    Variable(){}
    Variable(char* _name) : name(_name){
        Interval<1> cx(cij[0].first() - N, cij[0].last() + N);
        Interval<1> cy(cij[1].first() - N, cij[1].last() + N);
        cxy = Interval<2>(cx, cy);
        cell.initialize(cxy); cell = 0;
        cell_t.initialize(cij); cell_t = 0;
        wallx.initialize(cxy); wallx = 0;
        wally.initialize(cxy); wally = 0;
        cell(cij) = setups::ncvar[name];
    }
    void printInfo(int level = 0){
        std::cout << "Variable: " << name << std::endl;
        std::cout << "Domain: " << cij << " , " << cxy << std::endl;
        std::cout << "Boundary: " << B() << std::endl;
        if (level > 0){
            std::cout << "Maximum Value: " << max(cell) << ", Minimum Value: " << min(cell) << std::endl;
        }
        Interval<2> pij; pij = cij;
        while (pij[0].length() > 10 && pij[1].length() > 10){
            pij[0] = Interval<1>(0.75 * pij[0].first() + 0.25 * pij[0].last(),
                    0.25 * pij[0].first() + 0.75 * pij[0].last());
            pij[1] = Interval<1>(0.75 * pij[1].first() + 0.25 * pij[1].last(),
                    0.25 * pij[1].first() + 0.75 * pij[1].last());
        }
        if (level > 1){
            std::cout << "Sample Cell Value: " << pij << std::endl;
            std::cout << cell(pij);
        }
        if (level > 2){
            std::cout << "Sample Wall Value: " << pij << std::endl;
            std::cout << wallx(pij);
        }
        std::cout << std::endl;
    }
    void fixBoundary(){ fix(cell, cij); }
    //char *name;
    std::string name;
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
    Variable(){for (int i = 0; i < S; i++) name[i] = new char[1];}
    Variable(char* _name[]){
        for (int i = 0; i < S; i++){ name[i] = std::string(_name[i]); };
        Interval<1> cx(cij[0].first() - N, cij[0].last() + N);
        Interval<1> cy(cij[1].first() - N, cij[1].last() + N);
        cxy = Interval<2>(cx, cy);
        cell.initialize(cxy); cell = 0;
        cell_t.initialize(cij); cell_t = 0;
        wallx.initialize(cxy); wallx = 0;
        wally.initialize(cxy); wally = 0;
        for (int i = 0; i < S; i++) cell(cij).comp(i) = setups::ncvar[name[i]];
    }
    //~Variable(){ for (int i = 0; i < S; i++) delete name[i]; }
    void printInfo(int level = 0){
        std::cout << "Variable: (";
        for (int i = 0; i < S - 1; i++) std::cout << name[i] << ", ";
        std::cout << name[S -1] << ")" << std::endl;
        std::cout << "Domain: " << cij << " , " << cxy << std::endl;
        std::cout << "Boundary: " << B() << std::endl;
        if (level > 0){
            Element_t vmax, vmin;
            for (int i = 0; i < S; i++){
                vmax(i) = max(cell.comp(i));
                vmin(i) = min(cell.comp(i));
            }
            std::cout << "Maximum Value: " << vmax << ", Minimum Value: " << vmin << std::endl;
        }
        Interval<2> pij; pij = cij;
        while (pij[0].length() > 10 && pij[1].length() > 10){
            pij[0] = Interval<1>(0.75 * pij[0].first() + 0.25 * pij[0].last(),
                    0.25 * pij[0].first() + 0.75 * pij[0].last());
            pij[1] = Interval<1>(0.75 * pij[1].first() + 0.25 * pij[1].last(),
                    0.25 * pij[1].first() + 0.75 * pij[1].last());
        }
        if (level > 1){
            std::cout << "Sample Cell Value: " << pij << std::endl;
            std::cout << cell(pij);
        }
        if (level > 2){
            std::cout << "Sample Wall Value: " << pij << std::endl;
            std::cout << wallx(pij);
        }
        std::cout << std::endl;
    }
    void fixBoundary(){ fix(cell, cij); }
    //char *name[S];
    std::string name[S];
    Array<2, Element_t> cell, cell_t;
    Array<2, Vector<2, Element_t> > wallx, wally;
};

#endif
