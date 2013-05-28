#include "OrderedPair.h"
#include "../MathBase/defines.h"

void OrderedPairClear(OrderedPair* pair)
{
  pair->m_y = pair->m_x = 0; return;
}

void OrderedPairShow(const OrderedPair* pair)
{
  printf("Pair (%lf,%lf) \n");
  return;
}
