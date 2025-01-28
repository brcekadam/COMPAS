#ifndef __MS_gt_07_h__
#define __MS_gt_07_h__

#include "constants.h"
#include "typedefs.h"
#include "profiling.h"
#include "utils.h"

#include "BaseStar.h"
#include "MainSequence.h"

class BaseStar;
class MainSequence;

class MS_gt_07: virtual public BaseStar, public MainSequence {

public:

    MS_gt_07() { m_StellarType = STELLAR_TYPE::MS_GT_07; };
    
    MS_gt_07(const BaseStar &p_BaseStar, const bool p_Initialise = true) : BaseStar(p_BaseStar), MainSequence(p_BaseStar) {
        m_StellarType = STELLAR_TYPE::MS_GT_07;                                                                                                         // Set stellar type
        if (p_Initialise) Initialise();                                                                                                                 // Initialise if required
    }

    MS_gt_07* Clone(const OBJECT_PERSISTENCE p_Persistence, const bool p_Initialise = true) {
        MS_gt_07* clone = new MS_gt_07(*this, p_Initialise); 
        clone->SetPersistence(p_Persistence); 
        return clone; 
    }

    static MS_gt_07* Clone(MS_gt_07& p_Star, const OBJECT_PERSISTENCE p_Persistence, const bool p_Initialise = true) {
        MS_gt_07* clone = new MS_gt_07(p_Star, p_Initialise); 
        clone->SetPersistence(p_Persistence); 
        return clone; 
    }


protected:

    void Initialise() {
        CalculateTimescales();                                                                                                                          // Initialise timescales
        // Age for MS_GT_07 is carried over from CH stars switching to MS after spinning down, so not set to 0.0 here
        
        // Initialise core mass, luminosity, radius, and temperature if Shikauchi core mass prescription is used
        // Only do this once - this should not be done if a CH star spins down and becomes a MS star (when using CHE_MODE::PESSIMISTIC)
        if (OPTIONS->MainSequenceCoreMassPrescription() == CORE_MASS_PRESCRIPTION::SHIKAUCHI &&                                                         // Shikauchi core mass prescription?
            utils::Compare(m_MZAMS, SHIKAUCHI_LOWER_MASS_LIMIT) >= 0                         &&                                                         // ZAMS mass >= SHIKAUCHI_LOWER_MASS_LIMIT?
            m_Time <= 0.0) {                                                                                                                            // star not yet aged past creation?
                                                                                                                                                        // yes - initialise
            m_InitialMainSequenceCoreMass = MainSequence::CalculateInitialMainSequenceCoreMass(m_MZAMS);
            m_MainSequenceCoreMass        = m_InitialMainSequenceCoreMass;
            m_Luminosity                  = MainSequence::CalculateLuminosityShikauchi(m_MainSequenceCoreMass, m_InitialHeliumAbundance, m_Age);
            m_Radius                      = MainSequence::CalculateRadiusOnPhase(m_Mass, m_Age, m_RZAMS0);
            m_Temperature                 = BaseStar::CalculateTemperatureOnPhase_Static(m_Luminosity, m_Radius);
        }
    }


    // member functions - alphabetically

    double      CalculateCriticalMassRatioClaeys14(const bool p_AccretorIsDegenerate) const ;
    double      CalculateCriticalMassRatioHurleyHjellmingWebbink() const { return 0.33; }                                                               // As coded in BSE. Using the inverse owing to how qCrit is defined in COMPAS. See Hurley et al. 2002 sect. 2.6.1 for additional details.
    double      CalculateMassLossRateHurley();
    double      CalculateMassTransferRejuvenationFactor();

    ENVELOPE    DetermineEnvelopeType() const;
};

#endif // __MS_gt_07_h__
