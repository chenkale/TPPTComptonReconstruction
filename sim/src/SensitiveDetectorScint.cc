#include "SensitiveDetectorScint.hh"
#include "SessionManager.hh"
#include "SimMode.hh"
#include "out.hh"

#include <sstream>
#include <cstdint>

#include "G4Step.hh"
#include "G4ThreeVector.hh"
#include "G4SystemOfUnits.hh"
#include "G4EventManager.hh"
#include "G4Event.hh"

#include <G4VProcess.hh>

G4bool SensitiveDetectorScint_MultipleEvents::ProcessHits(G4Step *step, G4TouchableHistory *)
{
    const double edep = step->GetTotalEnergyDeposit();
    if (edep == 0) return true;

    SessionManager & SM = SessionManager::getInstance();
    ModeDepositionScint * Mode = static_cast<ModeDepositionScint*>(SM.SimMode);

    const G4StepPoint * postP  = step->GetPostStepPoint();
    const G4StepPoint * preP   = step->GetPreStepPoint();

    bool bTransport = false;
    const G4VProcess  * proc = postP->GetProcessDefinedStep();
    if (proc) bTransport = ( (proc->GetProcessType() == fTransportation) );

    //There could be energy deposition on Transportation (charged particles)!
    //if (!bTransport && postP->GetPhysicalVolume()->GetName() != "Scint") out("AAAAAAAAAAAAAAAAAAAAAAAAA");

    const int iScint = (bTransport ? preP ->GetPhysicalVolume()->GetCopyNo()
                                   : postP->GetPhysicalVolume()->GetCopyNo() );

    int64_t eventId = -1;
    const G4Event * event = G4EventManager::GetEventManager()->GetConstCurrentEvent();
    if (event) eventId = SM.makeGlobalEventId(event->GetEventID());

    DepositionNodeRecord newNode(postP->GetGlobalTime(), edep, eventId);
    std::vector<DepositionNodeRecord> & Nodes = Mode->DepositionData[iScint];
    if (!Mode->bDoCluster || Nodes.empty() || !Nodes.back().isCluster(newNode, Mode->MaxTimeDif))
    {
        if (Nodes.size() == Mode->MaxCapacity) Mode->saveData();
        Nodes.push_back(newNode);
    }
    else
        Nodes.back().merge(newNode);

    return true;
}
