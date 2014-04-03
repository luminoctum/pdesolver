#ifndef PHYSICS
#define PHYSICS
template<int O, class T>
class Advection{
    typedef T Element_t;
private:
    Array<2, Element_t> x_velocity, y_velocity;
    Advection(){};
    template<class A, class B>
    void set_velocity(const A& a, const B& b){
        x_velocity.initialize(a.cxy);
        y_velocity.initialize(b.cxy);
    }
}

#endif
