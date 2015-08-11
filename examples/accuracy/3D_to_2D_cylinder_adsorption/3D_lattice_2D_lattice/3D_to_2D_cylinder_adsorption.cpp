#include "../../GreensFunction2DAbsSym.hpp"
#include "../../GreensFunction2DRadAbs.hpp"
#include <iostream>

int main()
{
  using namespace greens_functions;

  const Real D(1e-12); //Diffusion coefficient
  const Real a(8e-6); //Distance from middle to the outer boundary
  const Real tau(a * a / (4 * D));

	GreensFunction2DAbsSym gf(D, a);
  std::cout << "log interval=0.01,world length_x=25,world length_y=800,world length_z=800,voxel radius=1e-08,[/:A][0]=1e-08" << std::endl;
  for (unsigned int i = 0; i < 625; ++i)
  {
    Real const T = (i + 1) * 0.001 * tau;
    std::cout << T << "," << gf.p_survival(T) << std::endl;
  }
  return 0;
}
