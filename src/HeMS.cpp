#include "HeMS.h"
#include "HeWD.h"


///////////////////////////////////////////////////////////////////////////////////////
//                                                                                   //
//             PARAMETERS, MISCELLANEOUS CALCULATIONS AND FUNCTIONS ETC.             //
//                                                                                   //
///////////////////////////////////////////////////////////////////////////////////////

/*
 * Calculate the helium abundance in the core of the star
 * 
 * Currently just a simple linear model. Should be updated to match detailed models.
 *
 * double CalculateHeliumAbundanceCore(const double p_Tau)
 * 
 * @param   [IN]    p_Tau                       Fraction of main sequence lifetime
 * 
 * @return                                      Helium abundance in the core (Y_c)
 */
double HeMS::CalculateHeliumAbundanceCoreOnPhase(const double p_Tau) const {
    double heliumAbundanceCoreMax = 1.0 - m_Metallicity;
    return heliumAbundanceCoreMax * (1.0 - p_Tau);
}


/*
 * Calculate timescales in units of Myr
 *
 * Timescales depend on a star's mass, so this needs to be called at least each timestep
 *
 * Vectors are passed by reference here for performance - preference would be to pass const& and
 * pass modified value back by functional return, but this way is faster - and this function is
 * called many, many times.
 *
 *
 * void CalculateTimescales(const double p_Mass, DBL_VECTOR &p_Timescales)
 *
 * @param   [IN]        p_Mass                  Mass in Msol
 * @param   [IN/OUT]    p_Timescales            Timescales
 */
void HeMS::CalculateTimescales(const double p_Mass, DBL_VECTOR &p_Timescales) {
#define timescales(x) p_Timescales[static_cast<int>(TIMESCALE::x)]  // for convenience and readability - undefined at end of function

    TPAGB::CalculateTimescales(p_Mass, p_Timescales);               // calculate common values

    timescales(tHeMS) = CalculateLifetimeOnPhase_Static(p_Mass);    // recalculate tHeMS

#undef timescales
}


/*
 * Calculate Giant Branch (GB) parameters per Hurley et al. 2000
 *
 * Giant Branch Parameters depend on a star's mass, so this needs to be called at least each timestep
 *
 * Vectors are passed by reference here for performance - preference would be to pass const& and
 * pass modified value back by functional return, but this way is faster - and this function is
 * called many, many times.
 *
 * void CalculateGBParams(const double p_Mass, DBL_VECTOR &p_GBParams)
 *
 * @param   [IN]        p_Mass                  Mass in Msol
 * @param   [IN/OUT]    p_GBParams              Giant Branch Parameters - calculated here
 */
void HeMS::CalculateGBParams(const double p_Mass, DBL_VECTOR &p_GBParams) {
#define gbParams(x) p_GBParams[static_cast<int>(GBP::x)]    // for convenience and readability - undefined at end of function

    GiantBranch::CalculateGBParams(p_Mass, p_GBParams);                                 // calculate common values (actually, all)

    // recalculate HeMS specific values

	gbParams(B) = CalculateCoreMass_Luminosity_B_Static();
	gbParams(D) = CalculateCoreMass_Luminosity_D_Static(p_Mass);

    gbParams(p) = CalculateCoreMass_Luminosity_p_Static(p_Mass, m_MassCutoffs);
    gbParams(q) = CalculateCoreMass_Luminosity_q_Static(p_Mass, m_MassCutoffs);
    
    gbParams(Mx) = GiantBranch::CalculateCoreMass_Luminosity_Mx_Static(p_GBParams);      // depends on B, D, p & q - recalculate if any of those are changed
    gbParams(Lx) = GiantBranch::CalculateCoreMass_Luminosity_Lx_Static(p_GBParams);      // depends on B, D, p, q & Mx - recalculate if any of those are changed

#undef gbParams
}


///////////////////////////////////////////////////////////////////////////////////////
//                                                                                   //
//                              LUMINOSITY CALCULATIONS                              //
//                                                                                   //
///////////////////////////////////////////////////////////////////////////////////////


/*
 * Calculate luminosity at ZAMS for a Helium Main Sequence star
 *
 * Hurley et al. 2000, eq 77
 *
 *
 * double CalculateLuminosityAtZAMS_Static(const double p_Mass)
 *
 * @param   [IN]    p_Mass                      Mass in Msol
 * @return                                      Luminosity at ZAMS for a Helium Main Sequence star in Lsol
 */
double HeMS::CalculateLuminosityAtZAMS_Static(const double p_Mass) {

    // pow() is slow - use multiplication (sqrt() is much faster than pow())
    double m_0_5   = std::sqrt(p_Mass);
    double m_3     = p_Mass * p_Mass * p_Mass;
    double m_6     = m_3 * m_3;
    double m_7_5   = m_6 * p_Mass * m_0_5;
    double m_9     = m_6 * m_3;
    double m_10_25 = m_9 * p_Mass * std::sqrt(m_0_5);

    return (15262.0 * m_10_25) / (m_9 + (29.54 * m_7_5) + (31.18 * m_6) + 0.0469);
}


/*
 * Calculate luminosity at the end of the Helium Main Sequence
 *
 * Hurley et al. 2000, eq 80 at tau = 1
 *
 *
 * double CalculateLuminosityAtPhaseEnd_Static(const double p_Mass)
 *
 * @param   [IN]    p_Mass                      Mass in Msol
 * @return                                      Luminosity at the end of the Helium Main Sequence in Lsol
 */
double HeMS::CalculateLuminosityAtPhaseEnd_Static(const double p_Mass) {
    return CalculateLuminosityAtZAMS_Static(p_Mass) * (1.0 + 0.45 + std::max(0.0, 0.85 - (0.08 * p_Mass)));
}


/*
 * Calculate luminosity for a Helium Main Sequence star (during central He burning)
 *
 * Hurley et al. 2000, eqs 80 & 82
 *
 *
 * double CalculateLuminosityOnPhase_Static(const double Mass, const double p_Tau)
 *
 * @param   [IN]    p_Mass                      Mass in Msol
 * @param   [IN]    p_Tau                       tHeMS relative age
 * @return                                      Luminosity for a Helium Main Sequence star in Lsol
 */
double HeMS::CalculateLuminosityOnPhase_Static(const double p_Mass, const double p_Tau) {

    double alpha = std::max(0.0, 0.85 - 0.08 * p_Mass);

    return CalculateLuminosityAtZAMS_Static(p_Mass) * (1.0 + (p_Tau * (0.45 + (alpha * p_Tau))));
}


///////////////////////////////////////////////////////////////////////////////////////
//                                                                                   //
//                                RADIUS CALCULATIONS                                //
//                                                                                   //
///////////////////////////////////////////////////////////////////////////////////////

/*
 * Calculate convective core radius
 *
 * Assume equal to total radius at start (for continuity with stripped CHeB or HG star), continuity with HeHG at end of phase, linear growth
 *
 *
 * double CalculateConvectiveCoreRadius()
 *
 * @return                                      Convective core radius (solar radii)
 */
double HeMS::CalculateConvectiveCoreRadius() const {

    // We need core radius at end of phase, which is just the core radius at the start of the HeHG phase.
    // Since we are on the He main sequence here, we can clone this object as an HeHG object
    // and, as long as it is initialised (to correctly set Tau to 0.0 on the HeHG phase),
    // we can query the cloned object for its core mass.
    //
    // The clone should not evolve, and so should not log anything, but to be sure the
    // clone does not participate in logging, we set its persistence to EPHEMERAL.
      
    HeHG *clone = HeHG::Clone(static_cast<HeHG&>(const_cast<HeMS&>(*this)), OBJECT_PERSISTENCE::EPHEMERAL);
    double finalConvectiveCoreRadius = clone->CalculateConvectiveCoreRadius();                  // get core radius from clone
    delete clone; clone = nullptr;                                                              // return the memory allocated for the clone

    double initialConvectiveCoreRadius = m_Radius;
    return (initialConvectiveCoreRadius - m_Tau * (initialConvectiveCoreRadius - finalConvectiveCoreRadius));
}


/*
 * Calculate radius at ZAMS for a Helium Main Sequence star
 *
 * Hurley et al. 2000, eq 78
 *
 *
 * double CalculateRadiusAtZAMS_Static(const double p_Mass)
 *
 * @param   [IN]    p_Mass                      Mass in Msol
 * @return                                      Radius at ZAMS for a Helium Main Sequence star in Rsol
 */
double HeMS::CalculateRadiusAtZAMS_Static(const double p_Mass) {
    // pow() is slow - use multiplication
    double m_3 = p_Mass * p_Mass * p_Mass;
    double m_4 = m_3 * p_Mass;

    return (0.2391 * PPOW(p_Mass, 4.6)) / (m_4 + (0.162 * m_3) + 0.0065);
}


/*
 * Calculate radius for a Helium Main Sequence star (during central He burning)
 *
 * Hurley et al. 2000, eq 81
 *
 *
 * double CalculateRadiusOnPhase_Static(const double p_Mass, const double p_Tau)
 *
 * @param   [IN]    p_Mass                      Mass in Msol
 * @param   [IN]    p_Tau                       tHeMS relative age
 * @return                                      Radius for a Helium Main Sequence star in Rsol
 */
double HeMS::CalculateRadiusOnPhase_Static(const double p_Mass, const double p_Tau) {

    // sanity check for mass - just return 0.0 if mass <= 0
    if (utils::Compare(p_Mass, 0.0) <= 0) return 0.0;

    double tau_6 = p_Tau * p_Tau * p_Tau * p_Tau * p_Tau * p_Tau;   // pow() is slow - use multiplication
    double beta  = std::max(0.0, 0.4 - 0.22 * log10(p_Mass));

    return CalculateRadiusAtZAMS_Static(p_Mass) * (1.0 + (beta * (p_Tau - tau_6)));
}


/*
 * Calculate the radius at the end of the helium main sequence
 *
 * Hurley et al. 2000, eq 81 at tau = 1
 *
 *
 * double CalculateRadiusAtPhaseEnd_Static(const double p_Mass)
 *
 * @param   [IN]    p_Mass                      Mass in Msol
 * @return                                      Radius at the end of the helium main sequence (RTHe)
 */
double HeMS::CalculateRadiusAtPhaseEnd_Static(const double p_Mass) {
    return CalculateRadiusOnPhase_Static(p_Mass, 1.0);
}


///////////////////////////////////////////////////////////////////////////////////////
//                                                                                   //
//                                 MASS CALCULATIONS                                 //
//                                                                                   //
///////////////////////////////////////////////////////////////////////////////////////


/*
 * Calculate convective core mass
 *
 * Assume equal to total mass at start (for continuity with stripped CHeB or HG star), continuity with HeHG at end of phase, linear growth
 *
 *
 * double CalculateConvectiveCoreMass()
 *
 * @return                                      Convective core mass (solar masses)
 */
double HeMS::CalculateConvectiveCoreMass() const {

    // We need core mass at end of phase, which is just the core mass at the start of the HeHG phase.
    // Since we are on the He main sequence here, we can clone this object as an HeHG object
    // and, as long as it is initialised (to correctly set Tau to 0.0 on the HeHG phase),
    // we can query the cloned object for its core mass.
    //
    // The clone should not evolve, and so should not log anything, but to be sure the
    // clone does not participate in logging, we set its persistence to EPHEMERAL.
      
    HeHG *clone = HeHG::Clone(static_cast<HeHG&>(const_cast<HeMS&>(*this)), OBJECT_PERSISTENCE::EPHEMERAL, true);
    double finalConvectiveCoreMass = clone->CoreMass();                                         // get core mass from clone
    delete clone; clone = nullptr;                                                              // return the memory allocated for the clone

    double initialConvectiveCoreMass = m_Mass;
    return (initialConvectiveCoreMass - m_Tau * (initialConvectiveCoreMass - finalConvectiveCoreMass));
}


/*
 * Calculate rejuvenation factor for stellar age based on mass lost/gained during mass transfer
 *
 * Description?
 *
 *
 * double CalculateMassTransferRejuvenationFactor()
 *
 * @return                                      Rejuvenation factor
 */
double HeMS::CalculateMassTransferRejuvenationFactor() {

    double fRej = 1.0;                                                                          // default value

    switch (OPTIONS->MassTransferRejuvenationPrescription()) {

        case MT_REJUVENATION_PRESCRIPTION::HURLEY:                                              // use default Hurley et al. 2000 prescription = 1.0
            break;

        case MT_REJUVENATION_PRESCRIPTION::STARTRACK:                                           // StarTrack 2008 prescription - section 5.6 of http://arxiv.org/pdf/astro-ph/0511811v3.pdf
            fRej = utils::Compare(m_Mass, m_MassPrev) <= 0 ? 1.0 : m_MassPrev / m_Mass;         // rejuvenation factor is unity for mass losing stars
            break;

        default:                                                                                // unknown prescription
            // the only way this can happen is if someone added a MT_REJUVENATION_PRESCRIPTION
            // and it isn't accounted for in this code.  We should not default here, with or without a warning.
            // We are here because the user chose a prescription this code doesn't account for, and that should
            // be flagged as an error and result in termination of the evolution of the star or binary.
            // The correct fix for this is to add code for the missing prescription or, if the missing
            // prescription is superfluous, remove it from the option.

            THROW_ERROR(ERROR::UNKNOWN_MT_REJUVENATION_PRESCRIPTION);                           // throw error
    }

    return fRej;
}


/*
 * Calculate the dominant mass loss mechanism and associated rate for the star
 * at the current evolutionary phase.
 *
 * According to Hurley et al. 2000
 *
 * double CalculateMassLossRateHurley()
 *
 * @return                                      Mass loss rate in Msol per year
 */
double HeMS::CalculateMassLossRateHurley() {
    double rateNJ = CalculateMassLossRateNieuwenhuijzenDeJager();
    double rateKR = CalculateMassLossRateKudritzkiReimers();
    double rateWR = OPTIONS->WolfRayetFactor()  * CalculateMassLossRateWolfRayet(0.0);        // use mu = 0.0 for Helium stars 

    m_DominantMassLossRate = MASS_LOSS_TYPE::GB;
    double dominantRate    = std::max(rateNJ, rateKR);

    if (utils::Compare(rateWR, dominantRate) > 0) {
        dominantRate           = rateWR;
        m_DominantMassLossRate = MASS_LOSS_TYPE::WR;
    }

    return dominantRate;
}


/*
 * Calculate the dominant mass loss mechanism and associated rate for the star at the current evolutionary phase
 *
 * According to Vink - based on implementation in StarTrack
 *
 * double CalculateMassLossRateBelczynski2010()
 *
 * @return                                      Mass loss rate in Msol per year
 */
double HeMS::CalculateMassLossRateBelczynski2010() {
    m_DominantMassLossRate = MASS_LOSS_TYPE::WR;
    return OPTIONS->WolfRayetFactor() * CalculateMassLossRateWolfRayetZDependent(0.0);  
}


/*
 * Calculate the mass-loss rate for Wolf--Rayet stars according to the
 * prescription of Shenar et al. 2019 (https://ui.adsabs.harvard.edu/abs/2019A%26A...627A.151S/abstract)
 * 
 * See their Eq. 6 and Table 5
 * 
 * We use the fitting coefficients for hydrogen poor WR stars 
 * The C4 (X_He) term is = 0 and is omitted
 * 
 * double CalculateMassLossRateWolfRayetShenar2019()
 *
 *
 * @return                                      Mass loss rate (in Msol yr^{-1})
 */
double HeMS::CalculateMassLossRateWolfRayetShenar2019() const {

    double logMdot = 0.0;
    double Teff    = m_Temperature * TSOL;

    // For H-poor WR stars (X_H < 0.05)
    const double C1 = -7.99;
    const double C2 =  0.97;
    const double C3 = -0.07;
    const double C5 =  0.89;

    logMdot = C1 + (C2 * log10(m_Luminosity)) + (C3 * log10(Teff)) + (C5 * m_Log10Metallicity); 

    return PPOW(10.0, logMdot); // Mdot 
}


/*
 * Calculate the mass loss rate for helium stars in the updated prescription
 * Uses Sander & Vink 2020 for Wolf--Rayet stars
 * 
 * double CalculateMassLossRateMerritt2024()
 *
 * @return                                      Mass loss rate in Msol per year
 */
double HeMS::CalculateMassLossRateMerritt2024() {

    double MdotWR = 0.0;

    m_DominantMassLossRate = MASS_LOSS_TYPE::WR;                                                                // set dominant mass loss rate

    switch (OPTIONS->WRMassLossPrescription()) {                                                                // which WR mass loss prescription?

        case WR_MASS_LOSS_PRESCRIPTION::SANDERVINK2023: {
            // calculate Sander & Vink 2020 mass-loss rate
            double MdotSanderVink2020 = CalculateMassLossRateWolfRayetSanderVink2020(0.0);

            // apply the Sander et al. 2023 temperature correction to the Sander & Vink 2020 rate
            double MdotSander2023     = CalculateMassLossRateWolfRayetTemperatureCorrectionSander2023(MdotSanderVink2020);

            // calculate Vink 2017 mass-loss rate
            double MdotVink2017       = CalculateMassLossRateHeliumStarVink2017();

            // use whichever gives the highest mass loss rate -- will typically be Vink 2017 for
            // low Mass or Luminosity, and Sander & Vink 2020 for high Mass or Luminosity

            MdotWR = OPTIONS->WolfRayetFactor() * std::max(MdotSander2023, MdotVink2017);

        } break;

        case WR_MASS_LOSS_PRESCRIPTION::SHENAR2019: {
            // mass loss rate for WR stars from Shenar+ 2019
            double MdotShenar2019 = CalculateMassLossRateWolfRayetShenar2019();                                 // OPTIONS->WolfRayetFactor()  is applied in Shenar2019 function

            // calculate Vink 2017 mass-loss rate
            double MdotVink2017   = CalculateMassLossRateHeliumStarVink2017();

            // apply a minimum of Vink 2017 mass-loss rate to avoid extrapolating to low luminosity
            MdotWR = OPTIONS->WolfRayetFactor() * std::max(MdotShenar2019, MdotVink2017);

        } break;

        case WR_MASS_LOSS_PRESCRIPTION::BELCZYNSKI2010:
            MdotWR = CalculateMassLossRateBelczynski2010(); // OPTIONS->WolfRayetFactor() is applied in Belczynski2010 function
            break;

        default:                                                                                                // unknown prescription
            // the only way this can happen is if someone added a WR_MASS_LOSS_PRESCRIPTION
            // and it isn't accounted for in this code.  We should not default here, with or without a warning.
            // We are here because the user chose a prescription this code doesn't account for, and that should
            // be flagged as an error and result in termination of the evolution of the star or binary.
            // The correct fix for this is to add code for the missing prescription or, if the missing
            // prescription is superfluous, remove it from the option.

            THROW_ERROR(ERROR::UNKNOWN_WR_MASS_LOSS_PRESCRIPTION);                                              // throw error
    }

    return MdotWR;
}


/*
 * Determines if mass transfer is unstable according to the critical mass ratio.
 *
 * See e.g de Mink et al. 2013, Claeys et al. 2014, and Ge et al. 2010, 2015, 2020 for discussions.
 *
 * Assumes this star is the donor; relevant accretor details are passed as parameters.
 * Critical mass ratio is defined as qCrit = mAccretor/mDonor.
 *
 * double HeMS::CalculateCriticalMassRatioClaeys14(const bool p_AccretorIsDegenerate) 
 *
 * @param   [IN]    p_AccretorIsDegenerate      Boolean indicating if accretor in degenerate (true = degenerate)
 * @return                                      Critical mass ratio for unstable MT 
 */
double HeMS::CalculateCriticalMassRatioClaeys14(const bool p_AccretorIsDegenerate) const {

    double qCrit;
                                                                                                                            
    qCrit = p_AccretorIsDegenerate
                ? OPTIONS->MassTransferCriticalMassRatioHeliumMSDegenerateAccretor()        // degenerate accretor
                : OPTIONS->MassTransferCriticalMassRatioHeliumMSNonDegenerateAccretor();    // non-degenerate accretor
                                                                                                                        
    return qCrit;
}


///////////////////////////////////////////////////////////////////////////////////////
//                                                                                   //
//                            LIFETIME / AGE CALCULATIONS                            //
//                                                                                   //
///////////////////////////////////////////////////////////////////////////////////////


/*
 * Calculate lifetime for a Helium Main Sequence star
 *
 * Hurley et al. 2000, eq 79
 *
 *
 * double CalculateLifetimeOnPhase_Static(const double p_Mass)
 *
 * @param   [IN]    p_Mass                      Mass in Msol
 * @return                                      Lifetime for a Helium Main Sequence star in Myr
 */
double HeMS::CalculateLifetimeOnPhase_Static(const double p_Mass) {

    // pow() is slow - use multiplication (sqrt() is much faster than pow())
    double m_4   = p_Mass * p_Mass * p_Mass * p_Mass;
    double m_6   = m_4 * p_Mass * p_Mass;
    double m_6_5 = m_6 * std::sqrt(p_Mass);

    return (0.4129 + (18.81 * m_4) + (1.853 * m_6)) / m_6_5;
}


/*
 * Recalculates the star's age after mass loss
 *
 * Hurley et al. 2000, section 7.1
 *
 * Modifies attribute m_Age
 *
 *
 * UpdateAgeAfterMassLoss()
 *
 */
void HeMS::UpdateAgeAfterMassLoss() {

    double tHeMS      = m_Timescales[static_cast<int>(TIMESCALE::tHeMS)];
    double tHeMSprime = CalculateLifetimeOnPhase_Static(m_Mass);

    m_Age *= tHeMSprime / tHeMS;
    
    CalculateTimescales(m_Mass, m_Timescales);
}


///////////////////////////////////////////////////////////////////////////////////////
//                                                                                   //
//                    MISCELLANEOUS FUNCTIONS / CONTROL FUNCTIONS                    //
//                                                                                   //
///////////////////////////////////////////////////////////////////////////////////////


/*
 * Choose timestep for evolution
 *
 * Given in the discussion in Hurley et al. 2000
 *
 *
 * ChooseTimestep(const double p_Time)
 *
 * @param   [IN]    p_Time                      Current age of star in Myr
 * @return                                      Suggested timestep (dt)
 */
double HeMS::ChooseTimestep(const double p_Time) const {
#define timescales(x) m_Timescales[static_cast<int>(TIMESCALE::x)]  // for convenience and readability - undefined at end of function

    double dtk = 0.05 * timescales(tHeMS);
    double dte = timescales(tHeMS) - p_Time;

    return std::max(std::min(dtk, dte), NUCLEAR_MINIMUM_TIMESTEP);

#undef timescales
}


/*
 * Modify the star after it loses its envelope
 *
 * Where necessary updates attributes of star (depending upon stellar type):
 *
 *     - m_StellarType
 *     - m_Timescales
 *     - m_GBParams
 *     - m_Luminosity
 *     - m_Radius
 *     - m_Mass
 *     - m_Mass0
 *     - m_CoreMass
 *     - m_HeCoreMass
 *     - m_COCoreMass
 *     - m_Age
 *
 *
 * STELLAR_TYPE ResolveEnvelopeLoss(bool p_Force)
 *
 * @param   [IN]    p_Force                     Boolean to indicate whether the resolution of the loss of the envelope should be performed
 *                                              without checking the precondition(s).
 *                                              Default is false.
 *
 * @return                                      Stellar Type to which star should evolve after losing envelope
 */
STELLAR_TYPE HeMS::ResolveEnvelopeLoss(bool p_Force) {

    STELLAR_TYPE stellarType = m_StellarType;

    if (p_Force || utils::Compare(m_Mass, 0.0) <= 0) {
        stellarType = STELLAR_TYPE::MASSLESS_REMNANT;
        m_Radius    = 0.0;
        m_Mass      = 0.0;
    }

    return stellarType;
}


/*
 * Set parameters for evolution to next phase and return Stellar Type for next phase
 *
 *
 * STELLAR_TYPE EvolveToNextPhase()
 *
 * @return                                      Stellar Type for next phase
 */
STELLAR_TYPE HeMS::EvolveToNextPhase() {
#define timescales(x) m_Timescales[static_cast<int>(TIMESCALE::x)]  // for convenience and readability - undefined at end of function

    m_Age = timescales(tHeMS);

    return STELLAR_TYPE::NAKED_HELIUM_STAR_HERTZSPRUNG_GAP;

#undef timescales
}



/* 
 * Interpolate Ge+ Critical Mass Ratios, for H-poor stars
 * 
 * Function to interpolate in mass and radius to calculate the stellar response of a He star to mass loss.
 * Functionally works the same as the interpolator for H-rich stars, except that there is only one variation
 * for the H-poor stars. 
 *
 * Function takes no input (unlike in the H-rich case) because the existing table only applies for fully conservative 
 * mass transfer and the GE fully adiabatic response, not the artificially isentropic one. Also only for Z=Zsol.
 *
 * Interpolation is done linearly in logM and logR. 
 * 
 * double HeMS::InterpolateGeEtAlQCrit()
 * 
 * @return                                       Interpolated value of either the critical mass ratio or zeta for given stellar mass / radius
 */ 
double HeMS::InterpolateGeEtAlQCrit() {

    // Get vector of masses from qCritTable
    GE_QCRIT_TABLE_HE qCritTable = QCRIT_GE_HE_STAR;
    DBL_VECTOR massesFromQCritTable = std::get<0>(qCritTable);
    GE_QCRIT_RADII_QCRIT_VECTOR_HE radiiQCritsFromQCritTable = std::get<1>(qCritTable);

    INT_VECTOR indices = utils::BinarySearch(massesFromQCritTable, m_Mass);
    int lowerMassIndex = indices[0];
    int upperMassIndex = indices[1];
    
    if (lowerMassIndex == -1) {                                                   // if masses are out of range, set to endpoints
        lowerMassIndex = 0; 
        upperMassIndex = 1;
    } 
    else if (upperMassIndex == -1) { 
        lowerMassIndex = massesFromQCritTable.size() - 2; 
        upperMassIndex = massesFromQCritTable.size() - 1;
    } 
    
    // Get vector of radii from qCritTable for the lower and upper mass indices
    std::vector<double> logRadiusVectorLowerMass = std::get<0>(radiiQCritsFromQCritTable[lowerMassIndex]);
    std::vector<double> logRadiusVectorUpperMass = std::get<0>(radiiQCritsFromQCritTable[upperMassIndex]);

    // Get the qCrit vector for the lower and upper mass bounds 
    std::vector<double> qCritVectorLowerMass = std::get<1>(radiiQCritsFromQCritTable[lowerMassIndex]);
    std::vector<double> qCritVectorUpperMass = std::get<1>(radiiQCritsFromQCritTable[upperMassIndex]);

    
    // Get vector of radii from qCritTable for both lower and upper masses
    INT_VECTOR indicesR0          = utils::BinarySearch(logRadiusVectorLowerMass, log10(m_Radius));
    int lowerRadiusLowerMassIndex = indicesR0[0];
    int upperRadiusLowerMassIndex = indicesR0[1];
    
    if (lowerRadiusLowerMassIndex == -1) {                                        // if radii are out of range, set to endpoints
        lowerRadiusLowerMassIndex = 0; 
        upperRadiusLowerMassIndex = 1; 
    }
    else if (upperRadiusLowerMassIndex == -1) {                                                   
        lowerRadiusLowerMassIndex = logRadiusVectorLowerMass.size() - 2; 
        upperRadiusLowerMassIndex = logRadiusVectorLowerMass.size() - 1; 
    }
    
    INT_VECTOR indicesR1          = utils::BinarySearch(logRadiusVectorUpperMass, log10(m_Radius));
    int lowerRadiusUpperMassIndex = indicesR1[0];
    int upperRadiusUpperMassIndex = indicesR1[1];
    
    if (lowerRadiusUpperMassIndex == -1) {                                        // if radii are out of range, set to endpoints
        lowerRadiusUpperMassIndex = 0; 
        upperRadiusUpperMassIndex = 1; 
    }
    else if (upperRadiusUpperMassIndex == -1) {                                                   
        lowerRadiusUpperMassIndex = logRadiusVectorUpperMass.size() - 2; 
        upperRadiusUpperMassIndex = logRadiusVectorUpperMass.size() - 1; 
    }
    
    // Set the 4 boundary points for the 2D interpolation
    double qLowLow = qCritVectorLowerMass[lowerRadiusLowerMassIndex];
    double qLowUpp = qCritVectorLowerMass[upperRadiusLowerMassIndex];
    double qUppLow = qCritVectorUpperMass[lowerRadiusUpperMassIndex];
    double qUppUpp = qCritVectorUpperMass[upperRadiusUpperMassIndex];
    
    double lowerLogRadiusLowerMass = logRadiusVectorLowerMass[lowerRadiusLowerMassIndex];
    double upperLogRadiusLowerMass = logRadiusVectorLowerMass[upperRadiusLowerMassIndex];
    double lowerLogRadiusUpperMass = logRadiusVectorUpperMass[lowerRadiusUpperMassIndex];
    double upperLogRadiusUpperMass = logRadiusVectorUpperMass[upperRadiusUpperMassIndex];
    
    double logLowerMass = log10(massesFromQCritTable[lowerMassIndex]);
    double logUpperMass = log10(massesFromQCritTable[upperMassIndex]);
    
    // Interpolate on logR first, then logM, using nearest neighbor for extrapolation
    double logRadius = log10(m_Radius);
    double qCritLowerMass = (logRadius < lowerLogRadiusLowerMass) ? qLowLow
                          : (logRadius > upperLogRadiusLowerMass) ? qLowUpp
                          : qLowLow + (upperLogRadiusLowerMass - logRadius) / (upperLogRadiusLowerMass - lowerLogRadiusLowerMass) * (qLowUpp - qLowLow);
    double qCritUpperMass = (logRadius < lowerLogRadiusUpperMass) ? qUppLow
                          : (logRadius > upperLogRadiusUpperMass) ? qUppUpp 
                          : qUppLow + (upperLogRadiusUpperMass - logRadius) / (upperLogRadiusUpperMass - lowerLogRadiusUpperMass) * (qUppUpp - qUppLow);

    double logMass = log10(m_Mass);
    return   (logMass < logLowerMass) ? qCritLowerMass
           : (logMass > logUpperMass) ? qCritUpperMass
           : qCritLowerMass + (logUpperMass - logMass) / (logUpperMass - logLowerMass) * (qCritUpperMass - qCritLowerMass);
}
