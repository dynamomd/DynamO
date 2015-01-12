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
    /*! \brief Function types are symbolic representation of
        non-polynomial functions.
    */
    template<class Arg, size_t FUNC_ID>
    class Function {
    public:
      /*! \brief The symbolic expression which is the argument of the
          Function. 
      */
      Arg _arg;
      
      /*! \brief Construct a Function with its argument. */
      Function(const Arg& a): _arg(a) {}
    };

    /*! \brief Enabling of symbolic operators for Function types 
     */
    template<class Arg, size_t FuncID>
    struct SymbolicOperators<Function<Arg, FuncID> > {
      static const bool value = true;
    };

#define CREATE_FUNCTION(HELPERNAME, NUMERIC_IMPL, ARG_DERIVATIVE, ID)	\
    template<class A,							\
	     typename = typename std::enable_if<!std::is_arithmetic<A>::value>::type> \
    Function<A, ID> HELPERNAME(const A& a) { return Function<A, ID>(a); } \
									\
    template<class A,							\
	     typename = typename std::enable_if<std::is_arithmetic<A>::value>::type> \
    constexpr auto HELPERNAME(const A& x) -> decltype(NUMERIC_IMPL)		\
    { return NUMERIC_IMPL; }						\
									\
    template<char Letter, class Arg1, class Arg2>			\
    auto substitution(const Function<Arg1, ID>& f, const VariableSubstitution<Letter, Arg2>& x) \
      -> decltype(HELPERNAME(substitution(f._arg, x)))			\
    { return HELPERNAME(substitution(f._arg, x)); }			\
									\
    template<class A>							\
    inline std::ostream& operator<<(std::ostream& os, const Function<A, ID>& s) \
    { return os << #HELPERNAME "(" << s._arg << ")"; }			\
									\
    template<char dVariable, class A>					\
    auto derivative(const Function<A, ID>& f, Variable<dVariable>)	\
      -> decltype(derivative(f._arg, Variable<dVariable>()) * (ARG_DERIVATIVE)) \
    { return derivative(f._arg, Variable<dVariable>()) * (ARG_DERIVATIVE); }
    
    CREATE_FUNCTION(sin, std::sin(x), cos(f._arg), 0)
    CREATE_FUNCTION(cos, std::cos(x), -sin(f._arg), 1)

    constexpr NullSymbol sin(const NullSymbol&) { return NullSymbol(); }
    constexpr UnitySymbol cos(const NullSymbol&) { return UnitySymbol(); }
  }
}
