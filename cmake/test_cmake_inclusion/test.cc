#include "AABBtree.hh"

int main()
{
  AABBtree::Box<double, 2> box_d;
  AABBtree::Box<float, 2>  box_f;
  AABBtree::Ray<double, 2> ray_d;
  AABBtree::Ray<float, 2>  ray_f;
  AABBtree::AABBtree<double, 2> tree_d;
  AABBtree::AABBtree<float, 2>  tree_f;
  return 0;
}