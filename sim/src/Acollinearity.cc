#include "Acollinearity.hh"

#include "G4RandomTools.hh"
#include "G4SystemOfUnits.hh"

#include <cmath>

void applyAcollinearityDeviation(G4ThreeVector & v)
{
    G4ThreeVector vCopy(v);
    constexpr double Sigma = 0.5 * deg / 2.35482;
    const double angle = G4RandGauss::shoot(0, Sigma);
    v.rotate(angle, v.orthogonal());
    v.rotate(G4UniformRand() * 2.0 * M_PI, vCopy);
}
