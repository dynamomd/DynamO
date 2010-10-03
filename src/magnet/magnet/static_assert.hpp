#pragma once

#define M_STATIC_ASSERT(expr, msg)			\
  do {							\
  extern char ASSERT_FAIL__##msg[(expr)? 1: -1];	\
  if (!expr) ASSERT_FAIL__##msg[0]++;			\
  } while (false)
