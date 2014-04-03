#ifndef VARIABLE
#define VARIABLE
#include "Pooma/Arrays.h"
#include "Boundary.h"
#include "EnoWeno.h"
#include "setups.h"

template<class T>
class VariableBase{
public:
    VariableBase(){}
    virtual ~VariableBase(){}
    virtual void printInfo(int level = 0) = 0;
    virtual void fixBoundary() = 0;
    virtual const char* getName(int s = 0) = 0;
    virtual Array<2, T> getData(int s = 0) = 0;
    virtual void updateTendency() = 0;
    virtual int getSize() = 0;
    virtual void wallConstruct() = 0;
};

template<int N, class T, class B>
class Variable : public VariableBase<T>, public B
{
    typedef T Element_t;
    typedef B Boundary_t;
public:
    Variable(){}
    Variable(char* _name) : name(_name) {
        cij = setups::cij;
        sicj = setups::sicj;
        cisj = setups::cisj;

        Interval<1> cx(cij[0].first() - N, cij[0].last() + N);
        Interval<1> cy(cij[1].first() - N, cij[1].last() + N);
        cxy = Interval<2>(cx, cy);

        cell.initialize(cxy); cell = 0;
        cell_t.initialize(cij); cell_t = 0;
        wallx.initialize(sicj); wallx = 0;
        wally.initialize(cisj); wally = 0;

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

    void updateTendency(){ cell(cij) += cell_t * setups::dt; cell_t = 0; };

    const char* getName(int) { return name.c_str(); }

    Array<2, T> getData(int){ 
        Array<2, T> result(cij);
        result = cell(cij);
        return result; 
    }

    int getSize(){ return 1; };

    void wallConstruct(){ eno.construct(cell, wallx, wally, cij); }

    std::string name;
    Interval<2> cxy, cij, sicj, cisj;
    Array<2, Element_t> cell, cell_t;
    Array<2, Vector<2, Element_t> > wallx, wally;
    EnoScheme<N> eno;
};

// default template arguments may not be used in partial specializations
template<int N, int S, class T, class B>
class Variable<N, Vector<S, T>, B> : public VariableBase<T>, public B
{
    typedef Vector<S, T> Element_t;
    typedef B Boundary_t;
public:
    Variable(){for (int i = 0; i < S; i++) name[i] = new char[1];}
    Variable(char* _name[]){
        for (int i = 0; i < S; i++){ name[i] = std::string(_name[i]); };
        cij = setups::cij;
        sicj = setups::sicj;
        cisj = setups::cisj;

        Interval<1> cx(cij[0].first() - N, cij[0].last() + N);
        Interval<1> cy(cij[1].first() - N, cij[1].last() + N);
        cxy = Interval<2>(cx, cy);

        cell.initialize(cxy); cell = 0;
        cell_t.initialize(cij); cell_t = 0;
        // this has to be cxy because eno(cell, wallx)
        wallx.initialize(sicj); wallx = 0;
        wally.initialize(cisj); wally = 0;

        for (int i = 0; i < S; i++) cell(cij).comp(i) = setups::ncvar[name[i]];
    }

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
            std::cout << "Sample WallY Value: " << pij << std::endl;
            std::cout << wally(pij);
        }
        if (level > 2){
            std::cout << "Sample WallX Value: " << pij << std::endl;
            std::cout << wallx(pij);
            std::cout << "Sample WallY Value: " << pij << std::endl;
            std::cout << wally(pij);
        }
        std::cout << std::endl << std::endl;
    }

    void fixBoundary(){ fix(cell, cij); }

    void updateTendency(){ cell(cij) += cell_t * setups::dt; cell_t = 0;};

    const char* getName(int s) { return name[s].c_str(); }

    Array<2, T> getData(int s){ 
        Array<2, T> result(cij);
        result = cell(cij).comp(s);
        return result; 
    }

    int getSize(){ return S; };

    void wallConstruct(){ eno.construct(cell, wallx, wally, cij); }

    std::string name[S];
    Interval<2> cxy, cij, sicj, cisj;
    Array<2, Element_t> cell, cell_t;
    Array<2, Vector<2, Element_t> > wallx, wally;
    EnoScheme<N> eno;
};

template<class T>
class VariableList{
public:
    Array<2, T> buffer;
    Interval<2> cij;
    VariableBase<T>* mass;
    std::vector<VariableBase<T>*> mvar, uvar, pvar;
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
    void updateTendency(){
        for (int i = 0; i < pvar.size(); i++) pvar[i]->updateTendency();
    }
    void wallConstruct(){
        for (int i = 0; i < pvar.size(); i++) pvar[i]->wallConstruct();
    }
    void ncwrite(double time){
        setups::current++;
        NcFile dataFile(setups::ncfile.c_str(), NcFile::Write);
        for (int i = 0; i < mvar.size(); i++){
            for (int s = 0; s < mvar[i]->getSize(); s++){
                buffer = mvar[i]->getData(s) / mass->getData();
                dataFile.get_var(mvar[i]->getName(s))->put_rec(&buffer(0, 0), setups::current);
            }
        }
        for (int i = 0; i < uvar.size(); i++){
            for (int s = 0; s < uvar[i]->getSize(); s++){
                buffer = uvar[i]->getData(s);
                dataFile.get_var(uvar[i]->getName(s))->put_rec(&buffer(0, 0), setups::current);
            }
        }
        dataFile.get_var("time")->put_rec(&time, setups::current);
    }
};

#endif

