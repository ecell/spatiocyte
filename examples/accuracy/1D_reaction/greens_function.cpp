#include "../../GreensFunction1DRadAbs.hpp"
#include <iostream>

int main()
{
  using namespace greens_functions;

  const Real D(1e-12); //Diffusion coefficient
  const Real a(8e-6); //Distance from middle to the outer boundary

  const Real kf(4.4e+6); //Reaction rate
  const Real r0(0.8e-6); //Distance from the middle to the molecule
  const Real sigma(10e-9); //Distance from the middle to the surface of filament

  const Real tau(a * a / (2 * D));

	GreensFunction1DRadAbs gf(D, kf, r0, sigma, a);
  std::cout << "log interval=0.01,world length_x=25,world length_y=800,world length_z=800,voxel radius=1e-08,[/:A][0]=1e-08" << std::endl;
  for (unsigned int i = 0; i < 313; ++i)
  {
    Real const T = (i + 1) * 0.001 * tau;
    std::cout << T << "," << gf.p_survival(T) << std::endl;
  }
  return 0;
}
