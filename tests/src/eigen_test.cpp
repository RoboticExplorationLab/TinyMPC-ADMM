#include <Eigen.h>
#include <iostream>

#ifdef __cplusplus
extern "C" {
#endif

int main()
{
  Eigen::Matrix2f m(2,2);
  m(0,0) = 3;
  m(1,0) = 2.5;
  m(0,1) = -1;
  m(1,1) = m(1,0) + m(0,1);
  std::cout << m << std::endl;

  Eigen::Vector2f vec, ub, lb;
  vec << -3, 3;
  ub << 1, 2;
  lb << -2, -1;
  vec = vec.cwiseMin(ub).cwiseMax(lb); 
  std::cout << vec << std::endl;
}

#ifdef __cplusplus
}
#endif