#ifndef ENOSCHEME
#define ENOSCHEME
#include "Godunov.h"
#include "Pooma/Arrays.h"

template<int ORDER>
class ENOScheme{
protected:
    static Array<2> enoc;
    static void setEnoc(){
        enoc.initialize(Interval<1>(-1, ORDER - 1), Interval<1>(0, ORDER - 1));
        double coeff1, coeff2, coeff3, result;
        for(int r = -1; r < ORDER; r++) for (int j = 0; j < ORDER; j++){
            result = 0.;
            for (int m = j + 1; m < ORDER + 1; m++){
                coeff2 = 0.; coeff3 = 1.;
                for (int l = 0; l < ORDER + 1; l++){
                    if (l == m) continue;
                    coeff1 = 1.;
                    for (int q = 0; q < ORDER + 1; q++){
                        if (q == m || q == l) continue;
                        coeff1 *= r - q + 1;
                    }
                    coeff2 += coeff1;
                    coeff3 *= m - l;
                }
                result += coeff2 / coeff3;
            }
            enoc(r, j) = result;
        }
    }
    template<class ET, class EG> Array<1, ET>
    static inline polyCoeff(const Array<1, ET, EG>& arr){
        Array<1, ET> result(arr.domain());
        result = arr;
        int n = arr.domain().length();
        for (int i = 1; i < n; i++)
            for (int j = n - 1; j >= i; j--)
                result(j) = (result(j) - result(j - 1)) / i;
        return result;
    }
    template<class A> void
    static inline adaptiveStencil(const A& ua, int& jmin, int& jmax){
        for (int k = 1; k < ORDER; k++)
            if ( fabs(polyCoeff(ua(Interval<1>(jmin - 1, jmax)))(k)) <
                 fabs(polyCoeff(ua(Interval<1>(jmin, jmax + 1)))(k))
                 )
                jmin = jmin - 1;
            else jmax = jmax + 1;
    }
public:
    ENOScheme(){ setEnoc(); }
    /* scalar reconstruction */
    template<class A, class I> void 
    static construct(A, I ci, Int2Type<0>){
        I si(ci.first(), ci.last() + 1);
        int jmin, jmax;
        A::wallm = 0; A::wallp = 0;
        for (int i = ci.first() - 1; i <= ci.last() + 1; i++){
            jmin = i; jmax = i;
            adaptiveStencil(A::cell, jmin, jmax);
            for (int j = 0; j < ORDER; j++){
                A::wallp(i) += enoc(i - jmin - 1, j) * A::cell(jmin + j);
                A::wallm(i + 1) += enoc(i - jmin, j) * A::cell(jmin + j);
            }
        }
    }
    /* componentwise reconstruction */
    template<class A, class I> void 
    static construct(A, I ci, Int2Type<1>){
        I si(ci.first(), ci.last() + 1);
        int jmin, jmax;
        enum {ncomp = A::Element_t::d1};
        A::wallm = 0; A::wallp = 0;
        for (int i = ci.first() - 1; i <= ci.last() + 1; i++){
            for (int s = 0; s < ncomp; s++){
                jmin = i; jmax = i;
                adaptiveStencil(A::cell.comp(s), jmin, jmax);
                for (int j = 0; j < ORDER; j++){
                    A::wallp(i)(s) += enoc(i - jmin - 1, j) * A::cell(jmin + j)(s);
                    A::wallm(i + 1)(s) += enoc(i - jmin, j) * A::cell(jmin + j)(s);
                }
            }
        }
    }
    /* characteristic reconstruction */
    template<class A, class I> void 
    static construct(A, const I& ci, Int2Type<2>){
        I si(ci.first(), ci.last() + 1);
        int jmin, jmax;
        enum {ncomp = A::Element_t::d1};
        // assuming double here
        Tensor<ncomp> lp, ep, rp;
        Vector<ncomp> vleft, vright;
        for (int i = ci.first() - 1; i <= ci.last() + 1; i++){
            vleft = 0.; vright = 0.;
            typename A::Array_t va(I(i - ORDER + 1, i + ORDER - 1));
            // want to get the number of components
            // For characteristic reconstruction, A must have dfluxdu method
            EigenSystem<ncomp>::solve(A::dfluxdu(i), lp, ep, rp);
            for (int j = i - ORDER + 1; j <= i + ORDER - 1; j++){
                va(j) = dot(transpose(lp), A::cell(j));
            }
            for (int s = 0; s < ncomp; s++){
                jmin = i; jmax = i;
                adaptiveStencil(va.comp(s), jmin, jmax);
                for (int j = 0; j < ORDER; j++){
                    vleft(s) += enoc(i - jmin - 1, j) * va(jmin + j)(s);
                    vright(s) += enoc(i - jmin, j) * va(jmin + j)(s);
                }
            }
            A::wallp(i) = dot(vleft, rp);
            A::wallm(i + 1) = dot(vright, rp);
        }
    }
};

// note: put Array behind ENOScheme
template<int ORDER> Array<2> ENOScheme<ORDER>::enoc;

//Array<2> template<int ORDER, class SYS> ENOScheme<ORDER, SYS>::enoc;

/* disalbe WENOScheme for the time being
template<int ORDER, int DIM, class SYS>
class WENOScheme : public ENOScheme<ORDER, DIM, SYS>
{
private:
    double eps;
    Array<1> delta, wenoc;
    template<class A>
    inline Array<1> setBeta(const A& ua, int i){
        Array<1> beta(ORDER);
        beta(0) = 13./12. * pow(ua(i - 2) - 2 * ua(i - 1) + ua(i), 2)
            + 0.25 * pow(ua(i - 2) - 4 * ua(i - 1) + 3 * ua(i), 2);
        beta(1) = 13./12. * pow(ua(i - 1) - 2 * ua(i) + ua(i + 1), 2)
            + 0.25 * pow(ua(i - 1) - ua(i + 1), 2);
        beta(2) = 13./12. * pow(ua(i) - 2 * ua(i + 1) + ua(i + 2), 2)
            + 0.25 * pow(3 * ua(i) - 4 * ua(i + 1) + ua(i + 2), 2);
        return beta;
    }
    void setDelta(){
        delta(0) = 1./10.;
        delta(1) = 3./5.;
        delta(2) = 3./10.;
    }
    template<class A>
    void setWenoc(const A& ua, int i){
        double sum = 0;
        Array<1> beta(ORDER);
        beta = setBeta(ua, i);
        for (int j = 0; j < ORDER; j++)
            sum += delta(j) / pow(eps + beta(j), 2);
        for (int j = 0; j < ORDER; j++)
            wenoc(j) = delta(j) / pow(eps + beta(j), 2) / sum;
    }
public:
    WENOScheme(){
        delta.initialize(ORDER);
        wenoc.initialize(ORDER);
        this -> setEnoc();
        this -> setDelta();
        eps = 1.E-6;
    }
    // wallFlux for vector 
    template<int SIZE, class ET, class EG> Array<DIM, Vector<SIZE, ET> > 
    wallFlux(const Array<DIM, Vector<SIZE, ET>, EG>& ua, const Interval<DIM>& ci){
        Interval<DIM> si(ci.first(), ci.last() + 1);
        int jmin, jmax;
        Tensor<SIZE, ET> lp, ep, rp;
        Vector<SIZE, ET> vleft, vright, vleft_weno, vright_weno;

        for (int i = ci.first() - 1; i <= ci.last() + 1; i++){
            vleft_weno = 0; vright_weno = 0;
            EigenSystem<SIZE>::solve(SYS::dFluxdu(ua(i)), lp, ep, rp);
            for (int j = i - ORDER + 1; j <= i + ORDER - 1; j++){
                _va_(j) = dot(transpose(lp), ua(j));
            }
            for (int jmin = - ORDER + 1; jmin <= 0; jmin++){
                vleft = 0.; vright = 0;
                for (int s = 0; s < SIZE; s++){
                    setWenoc(_va_.comp(s), i);
                    for (int j = 0; j < ORDER; j++){
                        vleft(s) += this->enoc(- jmin - 1, j) * _va_(jmin + i + j)(s);
                        vright(s) += this->enoc(- jmin, j) * _va_(jmin + i + j)(s);
                    }
                    vleft_weno(s) += wenoc(jmin + ORDER - 1) * vleft(s);
                    vright_weno(s) += wenoc(jmin + ORDER - 1) * vright(s);
                }
            }
            _uleftn_(i) = dot(vleft_weno, rp);
            _urightn_(i) = dot(vright_weno, rp);
        }
        return Godunov<DIM, SYS>::wallFlux(_urightn_(si - 1), _uleftn_(si));
    }
};
*/

#endif
