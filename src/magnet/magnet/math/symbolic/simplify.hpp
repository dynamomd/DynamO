/*  dynamo:- Event driven molecular dynamics simulator 
    http://www.dynamomd.org
    Copyright (C) 2011  Marcus N Campbell Bannerman <m.bannerman@gmail.com>

    This program is free software: you can redistribute it and/or
    modify it under the terms of the GNU General Public License
    version 3 as published by the Free Software Foundation.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/
#pragma once

namespace magnet {
  namespace math {

    template<class T> struct is_ratio { static const bool value = false; };
    template<std::intmax_t N, std::intmax_t D> struct is_ratio<ratio<N,D> > { static const bool value = true; };

    //The implementations below perform basic simplification of expressions
    //
    //These simplify expressions dramatically, so they have the highest priority
    template<class RHS>
    typename std::enable_if<!is_ratio<RHS>::value, NullSymbol>::type 
    multiply(const NullSymbol&, const RHS&, detail::choice<0>) { return NullSymbol(); }
    
    template<class LHS>
    typename std::enable_if<!is_ratio<LHS>::value, NullSymbol>::type 
    multiply(const LHS&, const NullSymbol&, detail::choice<0>) { return NullSymbol(); }

    template<class LHS>
    typename std::enable_if<!is_ratio<LHS>::value, LHS>::type 
    multiply(const LHS& l, const UnitySymbol& r, detail::choice<0>) { return l; }
    template<class RHS> 
    typename std::enable_if<!is_ratio<RHS>::value, RHS>::type 
    multiply(const UnitySymbol& l, const RHS& r, detail::choice<0>) { return r; }
    

    template<class LHS> 
    typename std::enable_if<!is_ratio<LHS>::value, LHS>::type 
    add(const LHS& l, const NullSymbol& r, detail::choice<0>) { return l; }
    template<class RHS> 
    typename std::enable_if<!is_ratio<RHS>::value, RHS>::type
    add(const NullSymbol& l, const RHS& r, detail::choice<0>) { return r; }

    template<class LHS> 
    typename std::enable_if<!is_ratio<LHS>::value, LHS>::type 
    subtract(const LHS& l, const NullSymbol&, detail::choice<0>) { return l; }
    template<class RHS> 
    auto subtract(const NullSymbol&, const RHS& r, detail::choice<0>) -> typename std::enable_if<!is_ratio<RHS>::value, decltype(-r)>::type { return -r; }
    
    template<class LHS> LHS divide(const LHS& l, const UnitySymbol&, detail::choice<0>) { return l; }
    template<char Letter> UnitySymbol divide(const Variable<Letter>& l, const Variable<Letter>&, detail::choice<0>) { return UnitySymbol(); }

    template<char Letter>
    PowerOp<Variable<Letter>, 2> multiply(const Variable<Letter>&, const Variable<Letter>&, detail::choice<0>)
    { return PowerOp<Variable<Letter>, 2>(Variable<Letter>()); }
    
    template<char Letter, size_t Order>
    PowerOp<Variable<Letter>, Order+1> multiply(const PowerOp<Variable<Letter>, Order>&, const Variable<Letter>&, detail::choice<0>)
    { return PowerOp<Variable<Letter>, Order+1>(Variable<Letter>()); }
    
    template<char Letter, size_t Order>
    PowerOp<Variable<Letter>, Order+1> multiply(const Variable<Letter>&, const PowerOp<Variable<Letter>, Order>&, detail::choice<0>)
    { return PowerOp<Variable<Letter>, Order+1>(Variable<Letter>()); }


    template<class stdratio>
    using ratio_wrap = ratio<stdratio::num, stdratio::den>;

    //Ratio operators (these are lower priority than above
    template<std::intmax_t Num1, std::intmax_t Denom1, std::intmax_t Num2, std::intmax_t Denom2>
    ratio_wrap<std::ratio_multiply<std::ratio<Num1, Denom1>, std::ratio<Num2, Denom2> > >
    multiply(const ratio<Num1, Denom1>&, const ratio<Num2, Denom2>&, detail::choice<1>)
    { return {};}

    template<std::intmax_t Num1, std::intmax_t Denom1, std::intmax_t Num2, std::intmax_t Denom2>
    ratio_wrap<std::ratio_add<std::ratio<Num1, Denom1>, std::ratio<Num2, Denom2> > >
    add(const ratio<Num1, Denom1>&, const ratio<Num2, Denom2>&, detail::choice<1>)
    { return {};}

    template<std::intmax_t Num1, std::intmax_t Denom1, std::intmax_t Num2, std::intmax_t Denom2>
    ratio_wrap<std::ratio_divide<std::ratio<Num1, Denom1>, std::ratio<Num2, Denom2> > >
    divide(const ratio<Num1, Denom1>&, const ratio<Num2, Denom2>&, detail::choice<1>)
    { return {};}

    template<std::intmax_t Num1, std::intmax_t Denom1, std::intmax_t Num2, std::intmax_t Denom2>
    ratio_wrap<std::ratio_subtract<std::ratio<Num1, Denom1>, std::ratio<Num2, Denom2> > >
    subtract(const ratio<Num1, Denom1>&, const ratio<Num2, Denom2>&, detail::choice<1>)
    { return {};}

    template<std::intmax_t Num1, std::intmax_t Denom1, std::intmax_t Num2, std::intmax_t Denom2>
    constexpr bool operator==(const ratio<Num1, Denom1>&, const ratio<Num2, Denom2>&)
    { return std::ratio_equal<ratio<Num1, Denom1>, ratio<Num2, Denom2> >::value; }

    template<class ratio_arg, class factor, class offset = std::ratio<0> >
    struct is_whole_factor {
      static const bool value = (std::ratio_divide<std::ratio_subtract<ratio_arg, offset>, factor>::den == 1);
    };

    //Specialisations of sine cosine for whole multiples of pi/2
    template<std::intmax_t num, std::intmax_t den,
	     typename = typename std::enable_if<is_whole_factor<std::ratio<num, den>, pi>::value>::type>
    constexpr NullSymbol sin(const ratio<num, den>&) { return NullSymbol(); }

    template<std::intmax_t num, std::intmax_t den,
	     typename = typename std::enable_if<is_whole_factor<std::ratio<num, den>, pi, decltype(pi()/ratio<2>())>::value>::type>
    constexpr UnitySymbol sin(const ratio<num, den>&) { return UnitySymbol(); }

    template<std::intmax_t num, std::intmax_t den,
	     typename = typename std::enable_if<is_whole_factor<std::ratio<num, den>, pi, decltype(pi()/ratio<2>())>::value>::type>
    constexpr NullSymbol cos(const ratio<num, den>&) { return NullSymbol(); }

    template<std::intmax_t num, std::intmax_t den,
	     typename = typename std::enable_if<is_whole_factor<std::ratio<num, den>, pi >::value>::type>
    constexpr UnitySymbol cos(const ratio<num, den>&) { return UnitySymbol(); }

    //Removal of sign via abs on compile-time constants!
    template<std::intmax_t num, std::intmax_t den>
    constexpr ratio<(num >= 0) ? num : -num, den> abs(const ratio<num, den>&) { return ratio<(num >= 0) ? num : -num, den>(); }

//    /*! \brief Simplify multiplication and addition operations. */
//    template<class LHS1, class RHS1, class RHS>
//    auto simplify(const MultiplyOp<AddOp<LHS1, RHS1>, RHS>& f)
//      -> decltype(simplify(f._l._l * f._r + f._l._r * f._r)) {
//      return simplify(f._l._l * f._r + f._l._r * f._r);
//    }
//
//    /*! \brief Simplify multiplication and addition operations. */
//    template<class LHS1, class RHS1, class RHS>
//    auto simplify(const MultiplyOp<RHS, AddOp<LHS1, RHS1> >& f)
//      -> decltype(simplify(f._l * f._r._l + f._l * f._r._r)) {
//      return simplify(f._l * f._r._l + f._l * f._r._r);
//    }
//
//    /*! \brief Simplify multiplication and addition operations. */
//    template<class LHS1, class RHS1, class LHS2, class RHS2>
//    auto simplify(const MultiplyOp<AddOp<LHS1, RHS1>,  AddOp<LHS2, RHS2> >& f)
//      -> decltype(simplify(f._l._l * f._r._l + f._l._l * f._r._r + f._l._r * f._r._l + f._l._r * f._r._r)) {
//      return simplify(f._l._l * f._r._l + f._l._l * f._r._r + f._l._r * f._r._l + f._l._r * f._r._r);
//    }
    
    namespace detail {
      template<class T>
      auto try_simplify_imp(const T& a, int) -> decltype(simplify(a)) {
	return simplify(a);
      }
      
      template<class T>
      const T& try_simplify_imp(const T& a, long) {
	return a;
      }
    }

    template<class T>
    auto try_simplify(const T& a) -> decltype(detail::try_simplify_imp(a,0)) {
      return detail::try_simplify_imp(a, 0);
    }



    template<class Arg, size_t Power>
    auto simplify_powerop_impl(const PowerOp<Arg, Power>& f, detail::choice<0>) -> decltype(simplify(PowerOpSubstitution<Power>::eval(simplify(f._arg))))
    { return PowerOpSubstitution<Power>::eval(simplify(f._arg)); }
    
    template<class Arg, size_t Power>
    auto simplify_powerop_impl(const PowerOp<Arg, Power>& f, detail::choice<1>) -> decltype(PowerOpSubstitution<Power>::eval(simplify(f._arg)))
    { return PowerOpSubstitution<Power>::eval(simplify(f._arg)); }
    
    template<class Arg, size_t Power>
    auto simplify_powerop_impl(const PowerOp<Arg, Power>& f, detail::choice<2>) -> decltype(simplify(PowerOpSubstitution<Power>::eval(f._arg)))
    { return simplify(PowerOpSubstitution<Power>::eval(f._arg)); }

    template<class Arg, size_t Power>
    auto simplify_powerop_impl(const PowerOp<Arg, Power>& f, detail::choice<3>) -> decltype(PowerOpSubstitution<Power>::eval(f._arg))
    { return PowerOpSubstitution<Power>::eval(f._arg); }

    template<class T> struct PowerOpEnableExpansion { static const bool value = true; };
    template<char Letter> struct PowerOpEnableExpansion<Variable<Letter> > { static const bool value = false; };

    /*! \brief Expansion operator for PowerOp types. 
    
      This implementation only works if the argument has an simplify
      function defined for its argument.
     */
    template<class Arg, size_t Power,
	     typename = typename std::enable_if<PowerOpEnableExpansion<Arg>::value>::type>
    auto simplify(const PowerOp<Arg, Power>& f) -> decltype(simplify_powerop_impl(f, detail::select_overload{}))
    { return simplify_powerop_impl(f, detail::select_overload{}); }

    /*! \brief Optimisation of a Polynomial LHS multiplied by a
      Variable. */
    template<char Letter, size_t Order, class Real>
    Polynomial<Order+1, Real, Letter> simplify(const MultiplyOp<Variable<Letter>, Polynomial<Order, Real, Letter> >& f)
    { 
      Polynomial<Order+1, Real, Letter> retval;
      retval[0] = 0;
      std::copy(f._r.begin(), f._r.end(), retval.begin() + 1);
      return retval;
    }

    /*! \brief Optimisation of a Polynomial RHS multiplied by a
      Variable. */
    template<char Letter, size_t Order, class Real>
    Polynomial<Order+1, Real, Letter> simplify(const MultiplyOp<Polynomial<Order, Real, Letter>, Variable<Letter> >& f)
    {
      Polynomial<Order+1, Real, Letter> retval;
      retval[0] = 0;
      std::copy(f._l.begin(), f._l.end(), retval.begin() + 1);
      return retval;
    }

    /*! \brief Optimisation of a Polynomial LHS added to a
      Variable. */
    template<char Letter, size_t Order, class Real>
    Polynomial<Order+1, Real, Letter> simplify(const AddOp<Polynomial<Order, Real, Letter>, Variable<Letter> >& f)
    { 
      Polynomial<(Order > 0) ? Order : 1, Real, Letter> retval(f._l);
      retval[1] += 1;
      return retval;
    }

    /*! \brief Optimisation of a Polynomial RHS added to a
      Variable. */
    template<char Letter, size_t Order, class Real>
    Polynomial<Order+1, Real, Letter> simplify(const AddOp<Variable<Letter>, Polynomial<Order, Real, Letter> >& f)
    { 
      Polynomial<(Order > 0) ? Order : 1, Real, Letter> retval(f._r);
      retval[1] += 1;
      return retval;
    }

    /*! \brief Optimisation of a Polynomial LHS subtracted by a
      Variable. */
    template<char Letter, size_t Order, class Real>
    Polynomial<Order+1, Real, Letter> simplify(const SubtractOp<Polynomial<Order, Real, Letter>, Variable<Letter> >& f)
    {
      Polynomial<(Order > 0) ? Order : 1, Real, Letter> retval(f._l);
      retval[1] -= 1;
      return retval;
    }

    /*! \brief Optimisation of a Polynomial RHS subtracted by a
      Variable. */
    template<char Letter, size_t Order, class Real>
    Polynomial<Order+1, Real, Letter> simplify(const SubtractOp<Variable<Letter>, Polynomial<Order, Real, Letter> >& f)
    { 
      Polynomial<(Order > 0) ? Order : 1, Real, Letter> retval(-f._r);
      retval[1] += 1;
      return retval;
    }

    /*! \brief Optimisation of a Polynomial LHS added to a
      PowerOp of the Variable. */
    template<char Letter, size_t Order, class Real, size_t POrder>
    Polynomial<(Order > POrder) ? Order : POrder, Real, Letter> simplify(const AddOp<Polynomial<Order, Real, Letter>, PowerOp<Variable<Letter>, POrder> > & f)
    {
      Polynomial<(Order > POrder) ? Order : POrder, Real, Letter> retval(f._l);
      retval[POrder] += 1;
      return retval;
    }

    /*! \brief Optimisation of a Polynomial RHS added to a
      PowerOp of the Variable. */
    template<char Letter, size_t Order, class Real, size_t POrder>
    Polynomial<(Order > POrder) ? Order : POrder, Real, Letter> simplify(const AddOp<PowerOp<Variable<Letter>, POrder>, Polynomial<Order, Real, Letter> > & f)
    {
      Polynomial<(Order > POrder) ? Order : POrder, Real, Letter> retval(f._r);
      retval[POrder] += 1;
      return retval;
    }

    /*! \brief Optimisation of a Polynomial LHS subtracted from a
      PowerOp of the Variable. */
    template<char Letter, size_t Order, class Real, size_t POrder>
    Polynomial<(Order > POrder) ? Order : POrder, Real, Letter> simplify(const SubtractOp<Polynomial<Order, Real, Letter>, PowerOp<Variable<Letter>, POrder> >& f)
    {
      Polynomial<(Order > POrder) ? Order : POrder, Real, Letter> retval(f._r);
      retval[POrder] -= 1;
      return retval;
    }

    /*! \brief Optimisation of a Polynomial RHS subtracted from a
      PowerOp of the Variable. */
    template<char Letter, size_t Order, class Real, size_t POrder>
    Polynomial<(Order > POrder) ? Order : POrder, Real, Letter> simplify(const SubtractOp<PowerOp<Variable<Letter>, POrder>, Polynomial<Order, Real, Letter> >& f)
    {
      Polynomial<(Order > POrder) ? Order : POrder, Real, Letter> retval(-f._r);
      retval[POrder] += Real(1);
      return retval;
    }

    /*! \brief Optimisation of a Polynomial LHS multiplied by a
      PowerOp of a Variable. */
    template<char Letter, size_t Order, class Real, size_t POrder>
    Polynomial<Order+POrder, Real, Letter> simplify(const MultiplyOp<PowerOp<Variable<Letter>, POrder>, Polynomial<Order, Real, Letter> >& f)
    {
      Polynomial<Order+POrder, Real, Letter> retval;
      std::copy(f._r.begin(), f._r.end(), retval.begin() + POrder);
      return retval;
    }

    /*! \brief Optimisation of a Polynomial RHS multiplied by a
      PowerOp of a Variable. */
    template<char Letter, size_t Order, class Real, size_t POrder>
    Polynomial<Order+POrder, Real, Letter> simplify(const MultiplyOp<Polynomial<Order, Real, Letter>, PowerOp<Variable<Letter>, POrder> >& f)
    {
      Polynomial<Order+POrder, Real, Letter> retval;
      std::copy(f._l.begin(), f._l.end(), retval.begin() + POrder);
      return retval;
    }

    /*! \brief Optimisation of a Polynomial LHS added to a
      PowerOp of the Variable. */
    template<char Letter, size_t Order, class Real, size_t POrder>
    Polynomial<(Order > POrder) ? Order : POrder, Real, Letter> simplify(const AddOp<Polynomial<Order, Real, Letter>, UnitySymbol> & f)
    {
      Polynomial<(Order > POrder) ? Order : POrder, Real, Letter> retval(f._l);
      retval[0] += 1;
      return retval;
    }

    /*! \brief Optimisation of a Polynomial RHS added to a
      PowerOp of the Variable. */
    template<char Letter, size_t Order, class Real>
    Polynomial<Order, Real, Letter> simplify(const AddOp<UnitySymbol, Polynomial<Order, Real, Letter> > & f)
    {
      Polynomial<Order, Real, Letter> retval(f._r);
      retval[0] += 1;
      return retval;
    }

    /*! \brief Optimisation of a Polynomial LHS subtracted by a
      PowerOp of the Variable. */
    template<char Letter, size_t Order, class Real>
    Polynomial<Order, Real, Letter> simplify(const SubtractOp<Polynomial<Order, Real, Letter>, UnitySymbol> & f)
    {
      Polynomial<Order, Real, Letter> retval(f._l);
      retval[0] -= 1;
      return retval;
    }

    /*! \brief Optimisation of a Polynomial RHS subtracted by a
      PowerOp of the Variable. */
    template<char Letter, size_t Order, class Real>
    Polynomial<Order, Real, Letter> simplify(const SubtractOp<UnitySymbol, Polynomial<Order, Real, Letter> > & f)
    {
      Polynomial<Order, Real, Letter> retval(-f._r);
      retval[0] += 1;
      return retval;
    }

    /*! \brief Conversion of PowerOp RHS multiplied by a constant to a
      Polynomial. */
    template<char Letter, size_t Order, class Real>
    typename std::enable_if<std::is_arithmetic<Real>::value, Polynomial<Order, Real, Letter> >::type 
    simplify(const MultiplyOp<PowerOp<Variable<Letter>, Order>, Real>& f)
    {
      Polynomial<Order, Real, Letter> retval;
      retval[Order] = f._r;
      return retval;
    }

    /*! \brief Conversion of PowerOp LHS multiplied by a constant to a
      Polynomial. */
    template<char Letter, size_t Order, class Real>
    typename std::enable_if<std::is_arithmetic<Real>::value, Polynomial<Order, Real, Letter> >::type 
    simplify(const MultiplyOp<Real, PowerOp<Variable<Letter>, Order> >& f)
    {
      Polynomial<Order, Real, Letter> retval;
      retval[Order] = f._l;
      return retval;
    }

    /*! \brief Conversion of Variable RHS multiplied by a constant to a
      Polynomial. */
    template<char Letter, class Real>
    typename std::enable_if<std::is_arithmetic<Real>::value, Polynomial<1, Real, Letter> >::type 
    simplify(const MultiplyOp<Variable<Letter>, Real> & f)
    { return Polynomial<1, Real, Letter>{0, f._r}; }

    /*! \brief Conversion of Variable LHS multiplied by a constant to a
      Polynomial. */
    template<char Letter, class Real>
    typename std::enable_if<std::is_arithmetic<Real>::value, Polynomial<1, Real, Letter> >::type 
    simplify(const MultiplyOp<Real, Variable<Letter> > & f)
    { return Polynomial<1, Real, Letter>{0, f._l}; }

    /*! \brief Conversion of PowerOp LHS added with a constant to a
      Polynomial. */
    template<char Letter, size_t Order, class Real>
    typename std::enable_if<std::is_arithmetic<Real>::value, Polynomial<Order, Real, Letter> >::type 
    simplify(const AddOp<PowerOp<Variable<Letter>, Order>, Real>& f)
    { 
      Polynomial<Order, Real, Letter> retval;
      retval[Order] = 1;
      retval[0] = f._r;
      return retval;
    }

    /*! \brief Conversion of PowerOp LHS added with a constant to a
      Polynomial. */
    template<char Letter, size_t Order, class Real>
    typename std::enable_if<std::is_arithmetic<Real>::value, Polynomial<Order, Real, Letter> >::type 
    simplify(const AddOp<Real, PowerOp<Variable<Letter>, Order> >& f)
    { 
      Polynomial<Order, Real, Letter> retval;
      retval[Order] = 1;
      retval[0] = f._r;
      return retval;
    }

    /*! \brief Conversion of PowerOp LHS subtracted with a constant to a
      Polynomial. */
    template<char Letter, size_t Order, class Real>
    typename std::enable_if<std::is_arithmetic<Real>::value, Polynomial<Order, Real, Letter> >::type 
    simplify(const SubtractOp<PowerOp<Variable<Letter>, Order>, Real>& f)
    {
      Polynomial<Order, Real, Letter> retval;
      retval[Order] = 1;
      retval[0] = -f._r;
      return retval;
    }

    /*! \brief Conversion of PowerOp LHS subtracted with a constant to a
      Polynomial. */
    template<char Letter, size_t Order, class Real>
    typename std::enable_if<std::is_arithmetic<Real>::value, Polynomial<Order, Real, Letter> >::type 
    simplify(const SubtractOp<Real, PowerOp<Variable<Letter>, Order> >& f)
    { 
      Polynomial<Order, Real, Letter> retval;
      retval[Order] = -1;
      retval[0] = f._l;
      return retval;
    }

    /*! \brief Conversion of a Variable RHS added with a constant. */
    template<char Letter, class Real>
    typename std::enable_if<std::is_arithmetic<Real>::value, Polynomial<1, Real, Letter> >::type 
    simplify(const AddOp<Variable<Letter>, Real>& f)
    { return Polynomial<1, Real, Letter>{f._r, Real(1)}; }

    /*! \brief Conversion of a Variable LHS added with a constant. */
    template<char Letter, class Real>
    typename std::enable_if<std::is_arithmetic<Real>::value, Polynomial<1, Real, Letter> >::type 
    simplify(const AddOp<Real, Variable<Letter> >& f)
    { return Polynomial<1, Real, Letter>{f._r, Real(1)}; }


    /*! \brief Conversion of a Variable added with a Variable. */
    template<char Letter>
    Polynomial<1, int, Letter>
    simplify(const AddOp<Variable<Letter>, Variable<Letter> >& f)
    { return Polynomial<1, int, Letter>{0, 2}; }


    /*! \brief Conversion of a Variable RHS subtracted with a constant. */
    template<char Letter, class Real>
    typename std::enable_if<std::is_arithmetic<Real>::value, Polynomial<1, Real, Letter> >::type 
    simplify(const SubtractOp<Variable<Letter>, Real>& f)
    { return Polynomial<1, Real, Letter>{-f._r, Real(1)}; }

    /*! \brief Conversion of a Variable LHS subtracted with a constant. */
    template<char Letter, class Real>
    typename std::enable_if<std::is_arithmetic<Real>::value, Polynomial<1, Real, Letter> >::type
    simplify(const SubtractOp<Real, Variable<Letter> >& f)
    { return Polynomial<1, Real, Letter>{f._l, Real(-1)}; }


    //                      FUNCTION SIMPLIFICATION
    template<class Arg, size_t FuncID>
    auto simplify(const Function<Arg, FuncID>& f) -> decltype(Function<decltype(simplify(f._arg)), FuncID>(simplify(f._arg)))
    { return Function<decltype(simplify(f._arg)), FuncID>(simplify(f._arg)); }

    template<class LHS, class Arg>
    auto simplify(const MultiplyOp<LHS, arbsignF<Arg> >& f) 
      -> decltype(arbsign(try_simplify(f._l * f._r._arg)))
    { return arbsign(try_simplify(f._l * f._r._arg)); }

    template<class RHS, class Arg>
    auto simplify(const MultiplyOp<arbsignF<Arg>, RHS>& f)
      -> decltype(arbsign(try_simplify(f._l._arg * f._r)))
    { return arbsign(try_simplify(f._l._arg * f._r)); }

    template<class Arg1, class Arg2>
    auto simplify(const MultiplyOp<arbsignF<Arg1>, arbsignF<Arg2> >& f)
      -> decltype(arbsign(try_simplify(f._l._arg * f._r._arg)))
    { return arbsign(try_simplify(f._l._arg * f._r._arg)); }

    template<class LHS, class Arg>
    auto simplify(const DivideOp<LHS, arbsignF<Arg> >& f) 
      -> decltype(arbsign(try_simplify(f._l / f._r._arg)))
    { return arbsign(try_simplify(f._l / f._r._arg)); }

    template<class RHS, class Arg>
    auto simplify(const DivideOp<arbsignF<Arg>, RHS>& f)
      -> decltype(arbsign(try_simplify(f._l._arg / f._r)))
    { return arbsign(try_simplify(f._l._arg / f._r)); }

    template<class Arg1, class Arg2>
    auto simplify(const DivideOp<arbsignF<Arg1>, arbsignF<Arg2> >& f)
      -> decltype(arbsign(try_simplify(f._l._arg / f._r._arg)))
    { return arbsign(try_simplify(f._l._arg / f._r._arg)); }

    template<class Arg>
    auto simplify(const arbsignF<arbsignF<Arg> >& f)
      -> decltype(arbsign(try_simplify(f._arg._arg)))
    { return arbsign(try_simplify(f._arg._arg)); }

    //For even powers, remove the sign term
    template<class Arg, size_t Power>
    auto simplify(const PowerOp<arbsignF<Arg>,Power>& f)
      -> typename std::enable_if<!(Power % 2), decltype(pow<Power>(f._arg._arg))>::type
    { return pow<Power>(f._arg._arg); }

    //For odd powers, move the sign term outside
    template<class Arg, size_t Power>
    auto simplify(const PowerOp<arbsignF<Arg>,Power>& f)
      -> typename std::enable_if<Power % 2, decltype(arbsign(pow<Power>(f._arg._arg)))>::type
    { return arbsign(pow<Power>(f._arg._arg)); }    
  }
}
