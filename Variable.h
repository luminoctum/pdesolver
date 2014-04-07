#ifndef VARIABLE
#define VARIABLE
#include "Pooma/Arrays.h"
#include "Boundary.h"
#include "setups.h"
#include "NumericalScheme.h"
enum{ SnapshotView, MassView, AverageView, MassAverageView};

class VariableBase{
public:
    VariableBase(){}
    virtual ~VariableBase(){}
    virtual void printInfo(int level = 0) = 0;
    virtual void fixBoundary() = 0;
    virtual const char* getName(int s = 0) = 0;
    virtual void updateTendency() = 0;
    virtual int getSize() = 0;
    virtual void wallConstruct() = 0;
    virtual void updateView() = 0;
    virtual Array<2, _AtomicType_>& getView(int s = 0) = 0;
    _Constructor_<_SpatialOrder_> constructor;
    static long step_count;
    static Array<2, _AtomicType_> *mass;
};
long VariableBase::step_count = 0;
Array<2, _AtomicType_> *VariableBase::mass = 0;

template<class T = double, class B = Dependent, int ViewType = SnapshotView>
class Variable : public VariableBase, public B
{
    typedef T Element_t;
    typedef B Boundary_t;
public:
    Variable(){}
    Variable(std::string _name) : name(_name) {
        cij = setups::cij;
        sicj = setups::sicj;
        cisj = setups::cisj;
        cxy = setups::cxy;

        cell.initialize(cxy); cell = 0;
        tendency.initialize(cij); tendency = 0;
        view.initialize(cij); view = 0;
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

    void updateTendency(){ cell(cij) += tendency * setups::dt; tendency = 0; };

    void wallConstruct(){ constructor(cell, wallx, wally, cij); }
    
    void updateView(){
        switch (ViewType){
            case SnapshotView:
                if (step_count % setups::frame == 0) view = cell(cij); 
                break;
            case MassView:
                if (step_count % setups::frame == 0) view = cell(cij) / mass->read(cij); 
                break;
            case AverageView:
                view += cell(cij);
                if (step_count % setups::frame == 0) view = view / setups::frame;
                break;
            case MassAverageView:
                view += cell(cij) / mass->read(cij); 
                if (step_count % setups::frame == 0) view = view / setups::frame;
                break;
        }
    }

    const char* getName(int) { return name.c_str(); }

    Array<2, T>& getView(int){ return view; }

    int getSize(){ return 1; };


    std::string name;
    Interval<2> cxy, cij, sicj, cisj;
    Array<2, Element_t> cell, tendency, view;
    Array<2, Vector<2, Element_t> > wallx, wally;
};

// default template arguments may not be used in partial specializations
template<int S, class T, class B, int ViewType>
class Variable<Vector<S, T>, B, ViewType> : public VariableBase, public B
{
    typedef Vector<S, T> Element_t;
    typedef B Boundary_t;
public:
    Variable(){for (int i = 0; i < S; i++) name[i] = new char[1];}
    Variable(std::initializer_list<std::string> _name){
        int j = 0;
        for (auto i = _name.begin(), j = 0; i != _name.end(); i++, j++) { name[j] = *i; }
        cij = setups::cij;
        sicj = setups::sicj;
        cisj = setups::cisj;
        cxy = setups::cxy;

        cell.initialize(cxy); cell = 0;
        tendency.initialize(cij); tendency = 0;
        for (int i = 0; i < S; i++) {
            view[i].initialize(cij); view[i] = 0;
        }
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

    void updateTendency(){ cell(cij) += tendency * setups::dt; tendency = 0;};

    void wallConstruct(){ constructor(cell, wallx, wally, cij); }

    void updateView(){
        switch (ViewType){
            case SnapshotView:
                if (step_count % setups::frame == 0) 
                    for (int i = 0; i < S; i++) view[i] = cell(cij).comp(i);
                break;
            case MassView:
                if (step_count % setups::frame == 0) 
                    for (int i = 0; i < S; i++) view[i] = cell(cij).comp(i) / mass->read(cij); 
                break;
            case AverageView:
                for (int i = 0; i < S; i++) view[i] += cell(cij).comp(i);
                if (step_count % setups::frame == 0)
                    for (int i = 0; i < S; i++) view[i] = view[i] / setups::frame;
                break;
            case MassAverageView:
                for (int i = 0; i < S; i++) view[i] += cell(cij).comp(i) / mass->read(cij); 
                if (step_count % setups::frame == 0) 
                    for (int i = 0; i < S; i++) view[i] = view[i] / setups::frame;
                break;
        }
        return;
    }

    const char* getName(int s) { return name[s].c_str(); }

    Array<2, T>& getView(int s){ return view[s]; }

    int getSize(){ return S; };

    std::string name[S];
    Interval<2> cxy, cij, sicj, cisj;
    Array<2, Element_t> cell, tendency;
    Array<2, Vector<2, Element_t> > wallx, wally;
    Array<2, T> view[S];
};

class VariableList{
public:
    std::vector<VariableBase*> pvar, var;
    Array<2, _AtomicType_>*& mass;
    long& step_count;
    VariableList() : 
        step_count(VariableBase::step_count),
        mass(VariableBase::mass){}
    void printInfo(int level = 0){
        std::cout << "******  variable infomation *****" << std::endl;
        for (int i = 0; i < var.size(); i++) var[i]->printInfo(level);
        std::cout << "*********************************" << std::endl;
    }
    void fixBoundary(){
        for (int i = 0; i < var.size(); i++) var[i]->fixBoundary();
    }
    void updateTendency(){
        for (int i = 0; i < pvar.size(); i++) pvar[i]->updateTendency();
    }
    void wallConstruct(){
        for (int i = 0; i < pvar.size(); i++) pvar[i]->wallConstruct();
    }
    void updateView(){
        for (int i = 0; i < var.size(); i++) var[i]->updateView();
    }
    void ncwrite(double time){
        setups::current++;
        NcFile dataFile(setups::ncfile.c_str(), NcFile::Write);
        for (int i = 0; i < var.size(); i++){
            for (int s = 0; s < var[i]->getSize(); s++){
                dataFile.get_var(var[i]->getName(s))->put_rec(&var[i]->getView(s)(0, 0), setups::current);
                var[i]->getView(s) = 0;
            }
        }
        dataFile.get_var("time")->put_rec(&time, setups::current);
    }
};

#endif

