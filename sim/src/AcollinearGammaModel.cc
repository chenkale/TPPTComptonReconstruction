#include "AcollinearGammaModel.hh"
#include "G4Gamma.hh"

G4bool AcollinearGammaModel::IsApplicable(const G4ParticleDefinition & particle)
{
    return (&particle == G4Gamma::GammaDefinition());
}

G4bool AcollinearGammaModel::ModelTrigger(const G4FastTrack & fastTrack)
{
    // Acollinearity is applied at primary generation (see generateSecondGamma).
    // Keep the fast-sim model registered for compatibility but do not modify tracks here.
    (void)fastTrack;
    return false;
}

void AcollinearGammaModel::DoIt(const G4FastTrack & fastTrack, G4FastStep & step)
{
    (void)fastTrack;
    (void)step;
}
