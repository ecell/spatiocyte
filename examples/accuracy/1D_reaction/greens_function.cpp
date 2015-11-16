#include "../../GreensFunction1DRadAbs.hpp"
#include <iostream>

int main()
{
  using namespace greens_functions;
  //From spatiocyte:
  const Real p(0.001);
  const Real r(10e-9);

  //Greens function
  const Real D(1e-12); //Diffusion coefficient
  const Real a(8e-6); //Distance from middle to the outer boundary

  const Real kf(D*p/(2*r)); //Spatiocyte rate is k = (DA+DB)*p/(2*rv) (unit m/s)
  const Real r0(0.8e-6); //Distance from the middle to the molecule
  const Real sigma(10e-9); //Distance from the middle to the surface of center molecule 

  const Real tau(a * a / (2 * D));

	GreensFunction1DRadAbs gf(D, kf, r0, sigma, a);
  std::cout << "log interval=0.01,world length_x=25,world length_y=800,world length_z=800,voxel radius=1e-08,[/:A][0]=1e-08" << std::endl;
  for (unsigned int i = 0; i < 3124; ++i)
  {
    Real const T = (i + 1) * 0.001 * tau;
    std::cout << T << "," << gf.p_survival(T) << std::endl;
  }
  return 0;
}
