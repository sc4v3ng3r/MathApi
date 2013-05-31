#ifndef ORDERED_PAIR_H
#define ORDERED_PAIR_H

typedef struct OrderedPair {
  double m_x;
  double m_y;
}OrderedPair;

void OrderedPairClear(OrderedPair* pair);
void OrderedPairShow(const OrderedPair* pair);
#endif