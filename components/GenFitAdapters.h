#ifndef GENFIT_ADAPTERS_H
#define GENFIT_ADAPTERS_H

#include "GenFit/AbsBField.h"
#include "GenFit/AbsMaterialInterface.h"
#include "DD4hep/Fields.h"
#include "DDRec/MaterialManager.h"
#include "DDRec/Vector3D.h"
#include "TVector3.h"

/**
 * Adapter class to make DD4hep magnetic field usable in GenFit
 */
class DD4hepFieldAdapter : public genfit::AbsBField {
public:
    DD4hepFieldAdapter(const dd4hep::OverlayedField& field) : m_field(field) {}
    
    /**
     * Get magnetic field at position (in cm)
     */
    TVector3 get(const TVector3& pos) const override {
        // Convert TVector3 to dd4hep::Position (pos is in cm as required by DD4hep)
        dd4hep::Position ddPos(pos.X(), pos.Y(), pos.Z());
        
        // Get field at position (DD4hep returns field in tesla)
        dd4hep::Direction fieldVector = m_field.magneticField(ddPos);
        
        // Convert to TVector3 (units are kGauss for GenFit, 1 Tesla = 10 kGauss)
        return TVector3(fieldVector.x() * 10.0, 
                        fieldVector.y() * 10.0, 
                        fieldVector.z() * 10.0);
    }
    
private:
    dd4hep::OverlayedField m_field;
};

/**
 * Adapter class to make DD4hep material information usable in GenFit
 */
class DD4hepMaterialAdapter : public genfit::AbsMaterialInterface {
public:
    DD4hepMaterialAdapter(dd4hep::rec::MaterialManager* manager) : m_manager(manager) {}
    
    void getMaterialParameters(double& density, 
                               double& Z, 
                               double& A, 
                               double& radiationLength, 
                               double& meanExcitationEnergy,
                               TVector3& position) override {
        // Convert TVector3 to dd4hep::Vector3D (position is in cm)
        dd4hep::rec::Vector3D pos(position.X(), position.Y(), position.Z());
        
        // Get material at this position
        dd4hep::Material material = m_manager->materialAt(pos);
        
        // Set the requested parameters
        density = material.density() / (dd4hep::g/dd4hep::cm3); // Convert to g/cm³
        Z = material.Z();
        A = material.A() / (dd4hep::g/dd4hep::mole); // Convert to g/mole
        radiationLength = material.radLength() / dd4hep::cm; // Convert to cm
        meanExcitationEnergy = 16.0e-9 * pow(Z, 0.9); // Approximation in GeV
        
        // No need to modify position
    }
    
    bool initTrack(double& mom, 
                   int& pdg, 
                   TVector3& position) override {
        // This method is called at the beginning of a track extrapolation
        // We don't need to modify anything here for DD4hep
        return true;
    }
    
private:
    dd4hep::rec::MaterialManager* m_manager;
};

#endif // GENFIT_ADAPTERS_H