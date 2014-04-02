#ifndef VARIABLE
#define VARIABLE
#include "Pooma/Arrays.h"
#include "Boundary.h"
#include "setups.h"

class VariableBase{
public:
    VariableBase(){ cij = setups::cij; cell.initialize(cij); }
    virtual ~VariableBase(){}
    virtual void printInfo(int level = 0) = 0;
    virtual void fixBoundary() = 0;
    virtual const char* getName(int s = 0) = 0;
    //virtual Array<2, T> getData(int s = 0) = 0;
    Interval<2> cxy, cij;
    std::string vtype;
    int size;
    Array<2> cell;
};

template<int N, class T, class B = Periodic>
class Variable : public VariableBase, public B
{
    typedef T Element_t;
    typedef B Boundary_t;
public:
    Variable(){}
    Variable(char* _name) : name(_name) {
        vtype = "scalar"; size = 1;
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
    const char* getName(int s = 0) { return name.c_str(); }
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
        vtype = "vector"; size = S;
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
    const char* getName(int s = 0) { return name[s].c_str(); }
    std::string name[S];
    Array<2, Element_t> cell, cell_t;
    Array<2, Vector<2, Element_t> > wallx, wally;
};

class VariableList{
public:
    Array<2, double> buffer;
    Interval<2> cij;
    VariableBase* mass;
    std::vector<VariableBase*> mvar, uvar, pvar;
    VariableList(){ 
        cij = setups::cij;
        buffer.initialize(cij); 
    }
    void printInfo(int level = 0){
        std::cout << "****  mass coupled variable  ****" << std::endl;
        for (int i = 0; i < mvar.size(); i++) mvar[i]->printInfo(level);
        std::cout << "**** mass uncoupled variable ****" << std::endl;
        for (int i = 0; i < uvar.size(); i++) uvar[i]->printInfo(level);
        std::cout << "*********************************" << std::endl;
    }
    void fixBoundary(){
        for (int i = 0; i < mvar.size(); i++) mvar[i]->fixBoundary();
        for (int i = 0; i < uvar.size(); i++) uvar[i]->fixBoundary();
    }
    void ncwrite(double time){
        setups::current++;
        NcFile dataFile(setups::ncfile.c_str(), NcFile::Write);
        std::cout << "0" << std::endl;
        for (int i = 0; i < mvar.size(); i++){
            std::cout << "1" << std::endl;
            if (mvar[i]->vtype == "scalar"){
                std::cout << mvar[i]->getName() << std::endl;
                std::cout << mvar[i]->cell(cij) << std::endl;
                std::cout << mass->cell(cij) << std::endl;
                buffer = mvar[i]->cell(cij) / mass->cell(cij);
                std::cout << "3" << std::endl;
                dataFile.get_var(mvar[i]->getName())->put_rec(&buffer(0, 0), setups::current);
                std::cout << "4" << std::endl;
            } else if (mvar[i]->vtype == "vector"){
                std::cout << "1" << std::endl;
                for (int s = 0; s < mvar[i]->size; s++){
                    buffer = mvar[i]->cell(cij).comp(s) / mass->cell(cij);
                    dataFile.get_var(mvar[i]->getName(s))->put_rec(&buffer(0, 0), setups::current);
                }
            }
        }
        std::cout << "2" << std::endl;
        for (int i = 0; i < uvar.size(); i++){
            if (uvar[i]->vtype == "scalar"){
                buffer = uvar[i]->cell(cij);
                dataFile.get_var(uvar[i]->getName())->put_rec(&buffer(0, 0), setups::current);
            } else if (uvar[i]->vtype == "vector"){
                for (int s = 0; s < uvar[i]->size; s++){
                    buffer = uvar[i]->cell(cij).comp(s);
                    dataFile.get_var(uvar[i]->getName(s))->put_rec(&buffer(0, 0), setups::current);
                }
            }
        }
        std::cout << "3" << std::endl;
        dataFile.get_var("time")->put_rec(&time, setups::current);
    }
};

#endif
