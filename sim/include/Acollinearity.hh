#ifndef acollinearity_h
#define acollinearity_h

#include "G4ThreeVector.hh"

// Applies the same angular deviation used for annihilation gamma acollinearity
// (second gamma of a back-to-back pair). Called at primary generation so logged
// directions match the simulated tracks.
void applyAcollinearityDeviation(G4ThreeVector & direction);

#endif // acollinearity_h
