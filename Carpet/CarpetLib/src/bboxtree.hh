#ifndef BBOXTREE_HH
#define BBOXTREE_HH

#include <vector>

using namespace std;

template <typename T, int D>
class bboxtree {
  struct node {
    T lower, upper;
    bboxtree<T,D-1>* branch;
  };
  vector<node> bs;
};

template <typename T>
class bboxtree <T, 0> {
  // empty
};

#endif  // #ifndef BBOXTREE_HH
