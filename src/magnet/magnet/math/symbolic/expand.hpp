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
//    /*! \brief Expand multiplication and addition operations. */
//    template<class LHS1, class RHS1, class RHS>
//    auto expand(const MultiplyOp<AddOp<LHS1, RHS1>, RHS>& f)
//      -> decltype(expand(f._l._l * f._r + f._l._r * f._r)) {
//      return expand(f._l._l * f._r + f._l._r * f._r);
//    }
//
//    /*! \brief Expand multiplication and addition operations. */
//    template<class LHS1, class RHS1, class RHS>
//    auto expand(const MultiplyOp<RHS, AddOp<LHS1, RHS1> >& f)
//      -> decltype(expand(f._l * f._r._l + f._l * f._r._r)) {
//      return expand(f._l * f._r._l + f._l * f._r._r);
//    }
//
//    /*! \brief Expand multiplication and addition operations. */
//    template<class LHS1, class RHS1, class LHS2, class RHS2>
//    auto expand(const MultiplyOp<AddOp<LHS1, RHS1>,  AddOp<LHS2, RHS2> >& f)
//      -> decltype(expand(f._l._l * f._r._l + f._l._l * f._r._r + f._l._r * f._r._l + f._l._r * f._r._r)) {
//      return expand(f._l._l * f._r._l + f._l._l * f._r._r + f._l._r * f._r._l + f._l._r * f._r._r);
//    }
    
    namespace detail {
      template<class T>
      auto try_expand_imp(const T& a, int) -> decltype(expand(a)) {
	return expand(a);
      }
      
      template<class T>
      const T& try_expand_imp(const T& a, long) {
	return a;
      }
    }

    template<class T>
    auto try_expand(const T& a) -> decltype(detail::try_expand_imp(a,0)) {
      return detail::try_expand_imp(a, 0);
    }



    template<class Arg, size_t Power>
    auto expand_powerop_impl(const PowerOp<Arg, Power>& f, detail::choice<0>) -> decltype(expand(PowerOpSubstitution<Power>::eval(expand(f._arg))))
    { return PowerOpSubstitution<Power>::eval(expand(f._arg)); }
    
    template<class Arg, size_t Power>
    auto expand_powerop_impl(const PowerOp<Arg, Power>& f, detail::choice<1>) -> decltype(PowerOpSubstitution<Power>::eval(expand(f._arg)))
    { return PowerOpSubstitution<Power>::eval(expand(f._arg)); }
    
    template<class Arg, size_t Power>
    auto expand_powerop_impl(const PowerOp<Arg, Power>& f, detail::choice<2>) -> decltype(expand(PowerOpSubstitution<Power>::eval(f._arg)))
    { return expand(PowerOpSubstitution<Power>::eval(f._arg)); }

    template<class Arg, size_t Power>
    auto expand_powerop_impl(const PowerOp<Arg, Power>& f, detail::choice<3>) -> decltype(PowerOpSubstitution<Power>::eval(f._arg))
    { return PowerOpSubstitution<Power>::eval(f._arg); }

    template<class T> struct PowerOpEnableExpansion { static const bool value = true; };
    template<char Letter> struct PowerOpEnableExpansion<Variable<Letter> > { static const bool value = false; };

    /*! \brief Expansion operator for PowerOp types. 
    
      This implementation only works if the argument has an expand
      function defined for its argument.
     */
    template<class Arg, size_t Power,
	     typename = typename std::enable_if<PowerOpEnableExpansion<Arg>::value>::type>
    auto expand(const PowerOp<Arg, Power>& f) -> decltype(expand_powerop_impl(f, detail::select_overload{}))
    { return expand_powerop_impl(f, detail::select_overload{}); }

    /*! \brief Optimisation of a Polynomial LHS multiplied by a
      Variable. */
    template<char Letter, size_t Order, class Real>
    Polynomial<Order+1, Real, Letter> expand(const MultiplyOp<Variable<Letter>, Polynomial<Order, Real, Letter> >& f)
    { 
      Polynomial<Order+1, Real, Letter> retval;
      retval[0] = 0;
      std::copy(f._r.begin(), f._r.end(), retval.begin() + 1);
      return retval;
    }

    /*! \brief Optimisation of a Polynomial RHS multiplied by a
      Variable. */
    template<char Letter, size_t Order, class Real>
    Polynomial<Order+1, Real, Letter> expand(const MultiplyOp<Polynomial<Order, Real, Letter>, Variable<Letter> >& f)
    {
      Polynomial<Order+1, Real, Letter> retval;
      retval[0] = 0;
      std::copy(f._l.begin(), f._l.end(), retval.begin() + 1);
      return retval;
    }

    /*! \brief Optimisation of a Polynomial LHS added to a
      Variable. */
    template<char Letter, size_t Order, class Real>
    Polynomial<Order+1, Real, Letter> expand(const AddOp<Polynomial<Order, Real, Letter>, Variable<Letter> >& f)
    { 
      Polynomial<(Order > 0) ? Order : 1, Real, Letter> retval(f._l);
      retval[1] += 1;
      return retval;
    }

    /*! \brief Optimisation of a Polynomial RHS added to a
      Variable. */
    template<char Letter, size_t Order, class Real>
    Polynomial<Order+1, Real, Letter> expand(const AddOp<Variable<Letter>, Polynomial<Order, Real, Letter> >& f)
    { 
      Polynomial<(Order > 0) ? Order : 1, Real, Letter> retval(f._r);
      retval[1] += 1;
      return retval;
    }

    /*! \brief Optimisation of a Polynomial LHS subtracted by a
      Variable. */
    template<char Letter, size_t Order, class Real>
    Polynomial<Order+1, Real, Letter> expand(const SubtractOp<Polynomial<Order, Real, Letter>, Variable<Letter> >& f)
    {
      Polynomial<(Order > 0) ? Order : 1, Real, Letter> retval(f._l);
      retval[1] -= 1;
      return retval;
    }

    /*! \brief Optimisation of a Polynomial RHS subtracted by a
      Variable. */
    template<char Letter, size_t Order, class Real>
    Polynomial<Order+1, Real, Letter> expand(const SubtractOp<Variable<Letter>, Polynomial<Order, Real, Letter> >& f)
    { 
      Polynomial<(Order > 0) ? Order : 1, Real, Letter> retval(-f._r);
      retval[1] += 1;
      return retval;
    }

    /*! \brief Optimisation of a Polynomial LHS added to a
      PowerOp of the Variable. */
    template<char Letter, size_t Order, class Real, size_t POrder>
    Polynomial<(Order > POrder) ? Order : POrder, Real, Letter> expand(const AddOp<Polynomial<Order, Real, Letter>, PowerOp<Variable<Letter>, POrder> > & f)
    {
      Polynomial<(Order > POrder) ? Order : POrder, Real, Letter> retval(f._l);
      retval[POrder] += 1;
      return retval;
    }

    /*! \brief Optimisation of a Polynomial RHS added to a
      PowerOp of the Variable. */
    template<char Letter, size_t Order, class Real, size_t POrder>
    Polynomial<(Order > POrder) ? Order : POrder, Real, Letter> expand(const AddOp<PowerOp<Variable<Letter>, POrder>, Polynomial<Order, Real, Letter> > & f)
    {
      Polynomial<(Order > POrder) ? Order : POrder, Real, Letter> retval(f._r);
      retval[POrder] += 1;
      return retval;
    }

    /*! \brief Optimisation of a Polynomial LHS subtracted from a
      PowerOp of the Variable. */
    template<char Letter, size_t Order, class Real, size_t POrder>
    Polynomial<(Order > POrder) ? Order : POrder, Real, Letter> expand(const SubtractOp<Polynomial<Order, Real, Letter>, PowerOp<Variable<Letter>, POrder> >& f)
    {
      Polynomial<(Order > POrder) ? Order : POrder, Real, Letter> retval(f._r);
      retval[POrder] -= 1;
      return retval;
    }

    /*! \brief Optimisation of a Polynomial RHS subtracted from a
      PowerOp of the Variable. */
    template<char Letter, size_t Order, class Real, size_t POrder>
    Polynomial<(Order > POrder) ? Order : POrder, Real, Letter> expand(const SubtractOp<PowerOp<Variable<Letter>, POrder>, Polynomial<Order, Real, Letter> >& f)
    {
      Polynomial<(Order > POrder) ? Order : POrder, Real, Letter> retval(-f._r);
      retval[POrder] += Real(1);
      return retval;
    }

    /*! \brief Optimisation of a Polynomial LHS multiplied by a
      PowerOp of a Variable. */
    template<char Letter, size_t Order, class Real, size_t POrder>
    Polynomial<Order+POrder, Real, Letter> expand(const MultiplyOp<PowerOp<Variable<Letter>, POrder>, Polynomial<Order, Real, Letter> >& f)
    {
      Polynomial<Order+POrder, Real, Letter> retval;
      std::copy(f._r.begin(), f._r.end(), retval.begin() + POrder);
      return retval;
    }

    /*! \brief Optimisation of a Polynomial RHS multiplied by a
      PowerOp of a Variable. */
    template<char Letter, size_t Order, class Real, size_t POrder>
    Polynomial<Order+POrder, Real, Letter> expand(const MultiplyOp<Polynomial<Order, Real, Letter>, PowerOp<Variable<Letter>, POrder> >& f)
    {
      Polynomial<Order+POrder, Real, Letter> retval;
      std::copy(f._l.begin(), f._l.end(), retval.begin() + POrder);
      return retval;
    }

    /*! \brief Optimisation of a Polynomial LHS added to a
      PowerOp of the Variable. */
    template<char Letter, size_t Order, class Real, size_t POrder>
    Polynomial<(Order > POrder) ? Order : POrder, Real, Letter> expand(const AddOp<Polynomial<Order, Real, Letter>, UnitySymbol> & f)
    {
      Polynomial<(Order > POrder) ? Order : POrder, Real, Letter> retval(f._l);
      retval[0] += 1;
      return retval;
    }

    /*! \brief Optimisation of a Polynomial RHS added to a
      PowerOp of the Variable. */
    template<char Letter, size_t Order, class Real>
    Polynomial<Order, Real, Letter> expand(const AddOp<UnitySymbol, Polynomial<Order, Real, Letter> > & f)
    {
      Polynomial<Order, Real, Letter> retval(f._r);
      retval[0] += 1;
      return retval;
    }

    /*! \brief Optimisation of a Polynomial LHS subtracted by a
      PowerOp of the Variable. */
    template<char Letter, size_t Order, class Real>
    Polynomial<Order, Real, Letter> expand(const SubtractOp<Polynomial<Order, Real, Letter>, UnitySymbol> & f)
    {
      Polynomial<Order, Real, Letter> retval(f._l);
      retval[0] -= 1;
      return retval;
    }

    /*! \brief Optimisation of a Polynomial RHS subtracted by a
      PowerOp of the Variable. */
    template<char Letter, size_t Order, class Real>
    Polynomial<Order, Real, Letter> expand(const SubtractOp<UnitySymbol, Polynomial<Order, Real, Letter> > & f)
    {
      Polynomial<Order, Real, Letter> retval(-f._r);
      retval[0] += 1;
      return retval;
    }

    /*! \brief Conversion of PowerOp RHS multiplied by a constant to a
      Polynomial. */
    template<char Letter, size_t Order, class Real>
    typename std::enable_if<std::is_arithmetic<Real>::value, Polynomial<Order, Real, Letter> >::type 
    expand(const MultiplyOp<PowerOp<Variable<Letter>, Order>, Real>& f)
    {
      Polynomial<Order, Real, Letter> retval;
      retval[Order] = f._r;
      return retval;
    }

    /*! \brief Conversion of PowerOp LHS multiplied by a constant to a
      Polynomial. */
    template<char Letter, size_t Order, class Real>
    typename std::enable_if<std::is_arithmetic<Real>::value, Polynomial<Order, Real, Letter> >::type 
    expand(const MultiplyOp<Real, PowerOp<Variable<Letter>, Order> >& f)
    {
      Polynomial<Order, Real, Letter> retval;
      retval[Order] = f._l;
      return retval;
    }

    /*! \brief Conversion of Variable RHS multiplied by a constant to a
      Polynomial. */
    template<char Letter, class Real>
    typename std::enable_if<std::is_arithmetic<Real>::value, Polynomial<1, Real, Letter> >::type 
    expand(const MultiplyOp<Variable<Letter>, Real> & f)
    { return Polynomial<1, Real, Letter>{0, f._r}; }

    /*! \brief Conversion of Variable LHS multiplied by a constant to a
      Polynomial. */
    template<char Letter, class Real>
    typename std::enable_if<std::is_arithmetic<Real>::value, Polynomial<1, Real, Letter> >::type 
    expand(const MultiplyOp<Real, Variable<Letter> > & f)
    { return Polynomial<1, Real, Letter>{0, f._l}; }

    /*! \brief Conversion of PowerOp LHS added with a constant to a
      Polynomial. */
    template<char Letter, size_t Order, class Real>
    typename std::enable_if<std::is_arithmetic<Real>::value, Polynomial<Order, Real, Letter> >::type 
    expand(const AddOp<PowerOp<Variable<Letter>, Order>, Real>& f)
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
    expand(const AddOp<Real, PowerOp<Variable<Letter>, Order> >& f)
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
    expand(const SubtractOp<PowerOp<Variable<Letter>, Order>, Real>& f)
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
    expand(const SubtractOp<Real, PowerOp<Variable<Letter>, Order> >& f)
    { 
      Polynomial<Order, Real, Letter> retval;
      retval[Order] = -1;
      retval[0] = f._l;
      return retval;
    }

    /*! \brief Conversion of a Variable RHS added with a constant. */
    template<char Letter, class Real>
    typename std::enable_if<std::is_arithmetic<Real>::value, Polynomial<1, Real, Letter> >::type 
    expand(const AddOp<Variable<Letter>, Real>& f)
    { return Polynomial<1, Real, Letter>{f._r, Real(1)}; }

    /*! \brief Conversion of a Variable LHS added with a constant. */
    template<char Letter, class Real>
    typename std::enable_if<std::is_arithmetic<Real>::value, Polynomial<1, Real, Letter> >::type 
    expand(const AddOp<Real, Variable<Letter> >& f)
    { return Polynomial<1, Real, Letter>{f._r, Real(1)}; }


    /*! \brief Conversion of a Variable added with a Variable. */
    template<char Letter>
    Polynomial<1, int, Letter>
    expand(const AddOp<Variable<Letter>, Variable<Letter> >& f)
    { return Polynomial<1, int, Letter>{0, 2}; }


    /*! \brief Conversion of a Variable RHS subtracted with a constant. */
    template<char Letter, class Real>
    typename std::enable_if<std::is_arithmetic<Real>::value, Polynomial<1, Real, Letter> >::type 
    expand(const SubtractOp<Variable<Letter>, Real>& f)
    { return Polynomial<1, Real, Letter>{-f._r, Real(1)}; }

    /*! \brief Conversion of a Variable LHS subtracted with a constant. */
    template<char Letter, class Real>
    typename std::enable_if<std::is_arithmetic<Real>::value, Polynomial<1, Real, Letter> >::type 
    expand(const SubtractOp<Real, Variable<Letter> >& f)
    { return Polynomial<1, Real, Letter>{f._r, Real(-1)}; }
  }
}
