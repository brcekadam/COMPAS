#include "GiantBranch.h"
#include "HeMS.h"
#include "WhiteDwarfs.h"
#include "ONeWD.h"
#include "NS.h"
#include "BH.h"

#include "stellarUtils.h"


///////////////////////////////////////////////////////////////////////////////////////
//                                                                                   //
//                     COEFFICIENT AND CONSTANT CALCULATIONS ETC.                    //
//                                                                                   //
///////////////////////////////////////////////////////////////////////////////////////


/*
 * Calculate the Hydrogen rate constant AH'
 *
 * Hurley et al. 2000, just after eq 43 (before eq 44)
 *
 *
 * double CalculateHRateConstant_Static(const double p_Mass)
 *
 * @param   [IN]    p_Mass                      Mass in Msol
 * @return                                      Hydrogen rate constant (Msol Lsol^-1 Myr^-1)
 */
double GiantBranch::CalculateHRateConstant_Static(const double p_Mass) {
    return PPOW(10.0, std::max(-4.8, std::min((-5.7 + (0.8 * p_Mass)), (-4.1 + (0.14 * p_Mass)))));
}


///////////////////////////////////////////////////////////////////////////////////////
//                                                                                   //
//             PARAMETERS, MISCELLANEOUS CALCULATIONS AND FUNCTIONS ETC.             //
//                                                                                   //
///////////////////////////////////////////////////////////////////////////////////////


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
void GiantBranch::CalculateTimescales(const double p_Mass, DBL_VECTOR &p_Timescales) {
#define timescales(x) p_Timescales[static_cast<int>(TIMESCALE::x)]      // for convenience and readability - undefined at end of function
#define gbParams(x) m_GBParams[static_cast<int>(GBP::x)]                // for convenience and readability - undefined at end of function

    double p1   = gbParams(p) - 1.0;
    double q1   = gbParams(q) - 1.0;
    double p1_p = p1 / gbParams(p);
    double q1_q = q1 / gbParams(q);

    double LBGB = CalculateLuminosityAtPhaseBase_Static(p_Mass, m_AnCoefficients);

    MainSequence::CalculateTimescales(p_Mass, p_Timescales);   // calculate common values

    timescales(tinf1_FGB) = timescales(tBGB) + ((1.0 / (p1 * gbParams(AH) * gbParams(D))) * PPOW((gbParams(D) / LBGB), p1_p));
    timescales(tMx_FGB)   = timescales(tinf1_FGB) - ((timescales(tinf1_FGB) - timescales(tBGB)) * PPOW((LBGB / gbParams(Lx)), p1_p));
    timescales(tinf2_FGB) = timescales(tMx_FGB) + ((1.0 / (q1 * gbParams(AH) * gbParams(B))) * PPOW((gbParams(B) / gbParams(Lx)), q1_q));

    timescales(tHeI)      = CalculateLifetimeToHeIgnition(p_Mass, timescales(tinf1_FGB), timescales(tinf2_FGB));
    timescales(tHeMS)     = HeMS::CalculateLifetimeOnPhase_Static(p_Mass);

#undef gbParams
#undef timescales
}


/*
 * Calculate the Core mass - Luminosity relation parameter B
 *
 * Hurley et al. 2000, eqs 31 - 38
 *
 *
 * double CalculateCoreMass_Luminosity_B_Static(const double p_Mass)
 *
 * @param   [IN]    p_Mass                      Mass in Msol
 * @return                                      Core mass - Luminosity relation parameter B
 */
double GiantBranch::CalculateCoreMass_Luminosity_B_Static(const double p_Mass) {
    return std::max(3.0E4, (500.0 + (1.75E4 * PPOW(p_Mass, 0.6))));
}


/*
 * Calculate the Core mass - Luminosity relation parameter D
 *
 * Hurley et al. 2000, eqs 31 - 38
 *
 *
 * double CalculateCoreMass_Luminosity_D_Static(const double p_Mass, const double p_LogMetallicityXi, DBL_VECTOR &p_MassCutoffs)
 *
 * @param   [IN]    p_Mass                      Mass in Msol
 * @param   [IN]    p_LogMetallicityXi          log10(Metallicity / Zsol) - called xi in Hurley et al 2000
 * @param   [IN]    p_MassCutoffs               Mass cutoffs vector
 * @return                                      Core mass - Luminosity relation parameter D
 */
double GiantBranch::CalculateCoreMass_Luminosity_D_Static(const double p_Mass, const double p_LogMetallicityXi, const DBL_VECTOR &p_MassCutoffs) {
#define massCutoffs(x) p_MassCutoffs[static_cast<int>(MASS_CUTOFF::x)]  // for convenience and readability - undefined at end of function

    double D0   = 5.37 + (0.135 * p_LogMetallicityXi);
    double D1   = (0.975 * D0) - (0.18 * p_Mass);
    double logD = D0;

    if (utils::Compare(p_Mass, massCutoffs(MHeF)) > 0) {
        if (utils::Compare(p_Mass, 2.5) >= 0) {
            double D2 = (0.5 * D0) - (0.06 * p_Mass);
            logD      = std::max(std::max(-1.0, D1), D2);
        }
        else {  // Linear interpolation between end points
            double gradient  = (D0 - D1) / (massCutoffs(MHeF) - 2.5);
            double intercept = D0 - (massCutoffs(MHeF) * gradient);
            logD = (gradient * p_Mass) + intercept;
        }
    }

    return PPOW(10.0, logD);

#undef massCutoffs
}


/*
 * Calculate the Core mass - Luminosity relation parameter p
 *
 * Hurley et al. 2000, eqs 31 - 38
 *
 *
 * double CalculateCoreMass_Luminosity_p_Static(const double p_Mass, DBL_VECTOR &p_MassCutoffs)
 *
 * @param   [IN]    p_Mass                      Mass in Msol
 * @param   [IN]    p_MassCutoffs               Mass cutoffs vector
 * @return                                      Core mass - Luminosity relation parameter p
 */
double GiantBranch::CalculateCoreMass_Luminosity_p_Static(const double p_Mass, const DBL_VECTOR &p_MassCutoffs) {
#define massCutoffs(x) p_MassCutoffs[static_cast<int>(MASS_CUTOFF::x)]  // for convenience and readability - undefined at end of function

    double p = 6.0;

    if (utils::Compare(p_Mass, massCutoffs(MHeF)) > 0) {
        if (utils::Compare(p_Mass, 2.5) >= 0) {
            p = 5.0;
        }
        else {                                                          // linear interpolation between end points
            double gradient  = 1.0 / (massCutoffs(MHeF) - 2.5);         // will be negative
            double intercept = 5.0 - (2.5 * gradient);
            p = (gradient * p_Mass) + intercept;
        }
    }

    return p;

#undef massCutoffs
}


/*
 * Calculate the Core mass - Luminosity relation parameter q
 *
 * Hurley et al. 2000, eqs 31 - 38
 *
 *
 * double CalculateCoreMass_Luminosity_q_Static(const double p_Mass, DBL_VECTOR &p_MassCutoffs)
 *
 * @param   [IN]    p_Mass                      Mass in Msol
 * @param   [IN]    p_MassCutoffs               Mass cutoffs vector
 * @return                                      Core mass - Luminosity relation parameter q
 */
double GiantBranch::CalculateCoreMass_Luminosity_q_Static(const double p_Mass, const DBL_VECTOR &p_MassCutoffs) {
#define massCutoffs(x) p_MassCutoffs[static_cast<int>(MASS_CUTOFF::x)]  // for convenience and readability - undefined at end of function

    double q = 3.0;

    if (utils::Compare(p_Mass, massCutoffs(MHeF)) > 0) {
        if (utils::Compare(p_Mass, 2.5) >= 0) {
            q = 2.0;
        }
        else {                                                          // linear interpolation between end points
            double gradient  = 1.0 / (massCutoffs(MHeF) - 2.5);
            double intercept = 2.0 - (2.5 * gradient);
            q = (gradient * p_Mass) + intercept;
        }
    }

    return q;

#undef massCutoffs
}


/*
 * Calculate the Core mass - Luminosity relation parameter Mx
 *
 * Mx is the point at which the low- and high-luminosity approximations cross
 *
 * Hurley et al. 2000, eq 38
 *
 *
 * double CalculateCoreMass_Luminosity_Mx_Static(const DBL_VECTOR &p_GBParams)
 *
 * @param   [IN]    p_GBParams                  Giant Branch Parameters
 * @return                                      Core mass - Luminosity relation parameter Mx
 */
double GiantBranch::CalculateCoreMass_Luminosity_Mx_Static(const DBL_VECTOR &p_GBParams) {
#define gbParams(x) p_GBParams[static_cast<int>(GBP::x)]    // for convenience and readability - undefined at end of function

    return PPOW((gbParams(B) / gbParams(D)), (1.0 / (gbParams(p) - gbParams(q))));

#undef gbParams
}


/*
 * Calculate First Giant Branch (FGB) Core mass - Luminosity relation parameter Lx
 *
 * Hurley et al. 2000, eq 37
 *
 *
 * double CalculateCoreMass_Luminosity_Lx_Static(const DBL_VECTOR &p_GBParams)
 *
 * @param   [IN]    p_GBParams                  Giant Branch Parameters
 * @return                                      Core mass - Luminosity relation parameter Lx
 */
double GiantBranch::CalculateCoreMass_Luminosity_Lx_Static(const DBL_VECTOR &p_GBParams) {
#define gbParams(x) p_GBParams[static_cast<int>(GBP::x)]    // for convenience and readability - undefined at end of function
    // since the mass used here is the mass at crossover (Mx), these
    // should give the same answer - but we'll take the minimum anyway
    return std::min((gbParams(B) * PPOW(gbParams(Mx), gbParams(q))), (gbParams(D) * PPOW(gbParams(Mx), gbParams(p))));

#undef gbParams
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
void GiantBranch::CalculateGBParams(const double p_Mass, DBL_VECTOR &p_GBParams) {
#define gbParams(x) p_GBParams[static_cast<int>(GBP::x)]    // for convenience and readability - undefined at end of function

    gbParams(AH)     = CalculateHRateConstant_Static(p_Mass);
    gbParams(AHHe)   = CalculateHHeRateConstant_Static();
    gbParams(AHe)    = CalculateHeRateConstant_Static();

    gbParams(B)      = CalculateCoreMass_Luminosity_B_Static(p_Mass);
    gbParams(D)      = CalculateCoreMass_Luminosity_D_Static(p_Mass, LogMetallicityXi(), m_MassCutoffs);

    gbParams(p)      = CalculateCoreMass_Luminosity_p_Static(p_Mass, m_MassCutoffs);
    gbParams(q)      = CalculateCoreMass_Luminosity_q_Static(p_Mass, m_MassCutoffs);

    gbParams(Mx)     = CalculateCoreMass_Luminosity_Mx_Static(p_GBParams);      // depends on B, D, p & q - recalculate if any of those are changed
    gbParams(Lx)     = CalculateCoreMass_Luminosity_Lx_Static(p_GBParams);      // depends on B, D, p, q & Mx - recalculate if any of those are changed

    gbParams(McBAGB) = CalculateCoreMassAtBAGB(p_Mass);
    gbParams(McDU)   = CalculateCoreMassAt2ndDredgeUp_Static(gbParams(McBAGB));
    gbParams(McBGB)  = CalculateCoreMassAtBGB(p_Mass, p_GBParams);

    gbParams(McSN)   = CalculateCoreMassAtSupernova_Static(gbParams(McBAGB));

#undef gbParams
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
 *
 * This function exists to facilitate the calculation of gbParams in EAGB::ResolveEnvelopeLoss() for the
 * stellar type to which the star will evolve.  The calculations of some of the stellar attributes there
 * depend on new gbParams.  
 * JR: This really needs to be revisited one day - these calculations should really be performed after
 *     switching to the new stellar type, but other calculations are done (in the legacy code) before the 
 *     switch (see evolveOneTimestep() in star.cpp for EAGB stars in the legacy code).
 *
 *
 * void CalculateGBParams_Static(const double      p_Mass, 
 *                               const double      p_LogMetallicityXi, 
 *                               const DBL_VECTOR &p_MassCutoffs, 
 *                               const DBL_VECTOR &p_AnCoefficients, 
 *                               const DBL_VECTOR &p_BnCoefficients,* 
 *                                     DBL_VECTOR &p_GBParams)
 *
 * @param   [IN]        p_Mass                  Mass in Msol
 * @param   [IN]        p_LogMetallicityXi      log10(Metallicity / Zsol) - called xi in Hurley et al 2000
 * @param   [IN]        p_MassCutoffs           Mass cutoffs
 * @param   [IN]        p_AnCoefficients        a(n) coefficients
 * @param   [IN]        p_BnCoefficients        b(n) coefficients
 * @param   [IN/OUT]    p_GBParams              Giant Branch Parameters - calculated here
 */
void GiantBranch::CalculateGBParams_Static(const double      p_Mass, 
                                           const double      p_LogMetallicityXi, 
                                           const DBL_VECTOR &p_MassCutoffs, 
                                           const DBL_VECTOR &p_AnCoefficients, 
                                           const DBL_VECTOR &p_BnCoefficients, 
                                                 DBL_VECTOR &p_GBParams) {

#define gbParams(x) p_GBParams[static_cast<int>(GBP::x)]    // for convenience and readability - undefined at end of function

    gbParams(AH)     = CalculateHRateConstant_Static(p_Mass);
    gbParams(AHHe)   = CalculateHHeRateConstant_Static();
    gbParams(AHe)    = CalculateHeRateConstant_Static();

    gbParams(B)      = CalculateCoreMass_Luminosity_B_Static(p_Mass);
    gbParams(D)      = CalculateCoreMass_Luminosity_D_Static(p_Mass, p_LogMetallicityXi, p_MassCutoffs);

    gbParams(p)      = CalculateCoreMass_Luminosity_p_Static(p_Mass, p_MassCutoffs);
    gbParams(q)      = CalculateCoreMass_Luminosity_q_Static(p_Mass, p_MassCutoffs);

    gbParams(Mx)     = CalculateCoreMass_Luminosity_Mx_Static(p_GBParams);      // depends on B, D, p & q - recalculate if any of those are changed
    gbParams(Lx)     = CalculateCoreMass_Luminosity_Lx_Static(p_GBParams);      // depends on B, D, p, q & Mx - recalculate if any of those are changed

    gbParams(McDU)   = CalculateCoreMassAt2ndDredgeUp_Static(gbParams(McBAGB));
    gbParams(McBAGB) = CalculateCoreMassAtBAGB_Static(p_Mass, p_BnCoefficients);
    gbParams(McBGB)  = CalculateCoreMassAtBGB_Static(p_Mass, p_MassCutoffs, p_AnCoefficients, p_GBParams);

    gbParams(McSN)   = CalculateCoreMassAtSupernova_Static(gbParams(McBAGB));

#undef gbParams
}


/*
 * Calculate the perturbation parameter mu
 *
 * Hurley et al. 2000, eqs 97 & 98
 *
 *
 * double CalculatePerturbationMu()
 *
 * @return                                      Perturbation parameter mu
 */
double GiantBranch::CalculatePerturbationMu() const {
    double kappa = -0.5;
    double L0    = 7.0E4;

    return ((m_Mass - m_CoreMass) / m_Mass) * (std::min(5.0, std::max(1.2, PPOW((m_Luminosity / L0), kappa))));
}


/*
 * Perturb Luminosity and Radius
 *
 * Perturbs Luminosity and Radius per Hurley et al. 2000, section 6.3
 * The attributes of the star are updated.
 *
 * Perturbation is disabled by default when DEBUG is enabled - except when
 * DEBUG_PERTURB is defined (see below).  The stellar class perturbation
 * functions are defined away if DEBUG is defined - so the generic Star
 * function is called (and does nothing). (So far FGB is the only class that
 * defines this function where it actually does anything)
 *
 * If DEBUG_PERTURB is defined then perturbation is not disabled while debugging.
 * To enable perturbation while DEBUG is enabled, define DEBUG_PERTURB.
 *
 *
 * void PerturbLuminosityAndRadius()
 */
#if !defined(DEBUG) || defined(DEBUG_PERTURB)   // don't perturb if debugging and DEBUG_PERTURB not defined
void GiantBranch::PerturbLuminosityAndRadius() {

    if (utils::Compare(m_Mu, 1.0) < 0) {   // perturb only if mu < 1.0

        double Lc = CalculateRemnantLuminosity();
        double Rc = CalculateRemnantRadius();

        double s = CalculatePerturbationS(m_Mu, m_Mass);
        double r = CalculatePerturbationR(m_Mu, m_Mass, m_Radius, Rc);

        m_Luminosity = Lc * PPOW((m_Luminosity / Lc), s);        
        m_Radius     = Rc * PPOW((m_Radius / Rc), r);
    }
}
#else
void GiantBranch::PerturbLuminosityAndRadius() { }
#endif



///////////////////////////////////////////////////////////////////////////////////////
//                                                                                   //
//                              LUMINOSITY CALCULATIONS                              //
//                                                                                   //
///////////////////////////////////////////////////////////////////////////////////////


/*
 * Calculate luminosity at the Base of the Giant Branch
 *
 * Hurley et al. 2000, eq 10
 *
 *
 * double CalculateLuminosityAtPhaseBase_Static(const double p_Mass, const DBL_VECTOR &p_AnCoefficients)
 *
 * @param   [IN]    p_Mass                      Mass in Msol
 * @param   [IN]    p_AnCoefficients            a(n) coefficients
 * @return                                      Luminosity at the Base of the Giant Branch in Lsol
 *
 * p_AnCoefficients passed as parameter so function can be declared static
 */
double GiantBranch::CalculateLuminosityAtPhaseBase_Static(const double p_Mass, const DBL_VECTOR &p_AnCoefficients) {
#define a p_AnCoefficients  // for convenience and readability - undefined at end of function

    double top    = (a[27] * PPOW(p_Mass, a[31])) + (a[28] * PPOW(p_Mass, C_COEFF.at(2)));
    double bottom = a[29] + (a[30] * PPOW(p_Mass, C_COEFF.at(3))) + PPOW(p_Mass, a[32]);

    return top / bottom;

#undef a
}


/*
 * Calculate luminosity on the Zero Age Horizontal Branch
 * (for Low Mass stars)
 *
 * Hurley et al. 2000, eq 53
 *
 *
 * double CalculateLuminosityOnZAHB_Static(const double      p_Mass,
 *                                         const double      p_CoreMass,
 *                                         const double      p_Alpha1,
 *                                         const double      p_MHeF,
 *                                         const double      p_MFGB,
 *                                         const double      p_MinimumLuminosityOnPhase,
 *                                         const DBL_VECTOR &p_BnCoefficients)
 *
 * @param   [IN]    p_Mass                      Mass in Msol
 * @param   [IN]    p_CoreMass                  Core Mass in Msol
 * @param   [IN]    p_Alpha1                    alpha1 in Hurly et al. 2000 (just after eq 49)
 * @param   [IN]    p_MHeF                      Maximum initial mass for which helium ignites degenerately in a Helium Flash
 * @param   [IN]    p_MFGB                      Maximum initial mass for which helium ignites on the First Giant Branch
 * @param   [IN]    p_BnCoefficients            b(n) coefficients
 * @return                                      Luminosity on the Zero Age Horizontal Branch in Lsol
 *
 * p_Alpha1, p_MHeF, p_MFGB and p_BnCoefficients passed as parameters so function can be declared static
 */
double GiantBranch::CalculateLuminosityOnZAHB_Static(const double      p_Mass,
                                                     const double      p_CoreMass,
                                                     const double      p_Alpha1,
                                                     const double      p_MHeF,
                                                     const double      p_MFGB,
                                                     const double      p_MinimumLuminosityOnPhase,
                                                     const DBL_VECTOR &p_BnCoefficients) {
#define b p_BnCoefficients  // for convenience and readability - undefined at end of function

    double LZHe   = HeMS::CalculateLuminosityAtZAMS_Static(p_CoreMass);
    double LminHe = p_MinimumLuminosityOnPhase;

    double mu     = (p_Mass - p_CoreMass) / (p_MHeF - p_CoreMass);
    double alpha2 = (b[18] + LZHe - LminHe) / (LminHe - LZHe);

    double first  = (1.0 + b[20]) / (1.0 + (b[20] * PPOW(mu, 1.6479)));
    double second = (b[18] * PPOW(mu, b[19])) / (1.0 + alpha2 * exp(15.0 * (p_Mass - p_MHeF)));

    return LZHe + (first * second);

#undef b
}


/*
 * Calculate luminosity at Helium Ignition
 *
 * Hurley et al. 2000, eq 49
 *
 *
 * double CalculateLuminosityAtHeIgnition_Static(const double      p_Mass,
 *                                               const double      p_Alpha1,
 *                                               const double      p_MHeF,
 *                                               const DBL_VECTOR &p_BnCoefficients)
 *
 * @param   [IN]    p_Mass                      Mass in Msol
 * @param   [IN]    p_Alpha1                    alpha1 in Hurley et al. 2000 (just after eq 49)
 * @param   [IN]    p_MHeF                      Maximum initial mass for which helium ignites degenerately in a Helium Flash
 * @param   [IN]    p_BnCoefficients            b(n) coefficients
 * @return                                      Luminosity at Helium Ignition in Lsol
 *
 * p_Alpha1, p_MHeF and p_BnCoefficients passed as parameter so function can be declared static
 */
double GiantBranch::CalculateLuminosityAtHeIgnition_Static(const double      p_Mass,
                                                           const double      p_Alpha1,
                                                           const double      p_MHeF,
                                                           const DBL_VECTOR &p_BnCoefficients) {
#define b p_BnCoefficients  // for convenience and readability - undefined at end of function

    return (utils::Compare(p_Mass, p_MHeF) < 0)
            ? (b[9] * PPOW(p_Mass, b[10])) / (1.0 + (p_Alpha1 * exp(15.0 * (p_Mass - p_MHeF))))
            : (b[11] + (b[12] * PPOW(p_Mass, 3.8))) / (b[13] + (p_Mass * p_Mass));

#undef b
}


/*
 * Calculate luminosity of the remnant the star would become if it lost all of its
 * envelope immediately (i.e. M = Mc, coreMass)
 *
 * Hurley et al. 2000, just after eq 105
 *
 *
 * double CalculateRemnantLuminosity()
 *
 * @return                                      Luminosity of remnant core in Lsol
 */
double GiantBranch::CalculateRemnantLuminosity() const {
#define massCutoffs(x) m_MassCutoffs[static_cast<int>(MASS_CUTOFF::x)]  // for convenience and readability - undefined at end of function

    return (utils::Compare(m_Mass0, massCutoffs(MHeF)) > 0)
            ? HeMS::CalculateLuminosityAtZAMS_Static(m_CoreMass)
            : WhiteDwarfs::CalculateLuminosityOnPhase_Static(m_CoreMass, 0.0, m_Metallicity, WD_Baryon_Number.at(STELLAR_TYPE::HELIUM_WHITE_DWARF));

#undef massCutoffs
}


///////////////////////////////////////////////////////////////////////////////////////
//                                                                                   //
//                                RADIUS CALCULATIONS                                //
//                                                                                   //
///////////////////////////////////////////////////////////////////////////////////////


/*
 * Calculate radius on the Giant Branch
 *
 * Hurley et al. 2000, eq 46
 *
 *
 * double CalculateRadiusOnPhase_Static(const double p_Mass, const double p_Luminosity, const DBL_VECTOR &p_BnCoefficients)
 *
 * @param   [IN]    p_Mass                      Mass in Msol
 * @param   [IN]    p_Luminosity                Luminosity in Lsol
 * @param   [IN]    p_BnCoefficients            b(n) coefficients
 * @return                                      Radius on the First Giant Branch in Rsol
 *
 * p_BnCoefficients passed as parameter so function can be declared static
 */
double GiantBranch::CalculateRadiusOnPhase_Static(const double p_Mass, const double p_Luminosity, const DBL_VECTOR &p_BnCoefficients) {
#define b p_BnCoefficients  // for convenience and readability - undefined at end of function

    double A = std::min((b[4] * PPOW(p_Mass, -b[5])), (b[6] * PPOW(p_Mass, -b[7])));  // Hurley et al. 2000, just before eq 47

    return A * (PPOW(p_Luminosity, b[1]) + (b[2] * PPOW(p_Luminosity, b[3])));        // Hurley et al. 2000, eq 46

#undef b
}


/*
 * Calculate radius on the Zero Age Horizontal Branch
 *
 * Hurley et al. 2000, eq 54
 *
 *
 * double CalculateRadiusOnZAHB_Static(const double      p_Mass,
 *                                     const double      p_CoreMass,
 *                                     const double      p_Alpha1,
 *                                     const double      p_MHeF,
 *                                     const double      p_MFGB,
 *                                     const double      p_MinimumLuminosityOnPhase,
 *                                     const DBL_VECTOR &p_BnCoefficients)
 *
 * @param   [IN]    p_Mass                      Mass in Msol
 * @param   [IN]    p_CoreMass                  Core Mass in Msol
 * @param   [IN]    p_Alpha1                    alpha1 in Hurly et al. 2000 (just after eq 49)
 * @param   [IN]    p_MHeF                      Maximum initial mass for which helium ignites degenerately in a Helium Flash
 * @param   [IN]    p_MFGB                      Maximum initial mass for which helium ignites don the First Giant Branch
 * @param   [IN]    p_MinimumLuminosityOnPhase  Minimum luminosity on phase (only required for CHeB stars) - calculated once per star
 * @param   [IN]    p_BnCoefficients            b(n) coefficients
 * @return                                      Radius on the Zero Age Horizontal Branch in Rsol
 *
 * p_Alpha1, p_MHeF and p_BnCoefficients passed as parameters so function can be declared static
 */
double GiantBranch::CalculateRadiusOnZAHB_Static(const double      p_Mass,
                                                 const double      p_CoreMass,
                                                 const double      p_Alpha1,
                                                 const double      p_MHeF,
                                                 const double      p_MFGB,
                                                 const double      p_MinimumLuminosityOnPhase,
                                                 const DBL_VECTOR &p_BnCoefficients) {
#define b p_BnCoefficients  // for convenience and readability - undefined at end of function

    double RZHe  = HeMS::CalculateRadiusAtZAMS_Static(p_CoreMass);
    double LZAHB = CalculateLuminosityOnZAHB_Static(p_Mass, p_CoreMass, p_Alpha1, p_MHeF, p_MFGB, p_MinimumLuminosityOnPhase, p_BnCoefficients);
    double RGB   = GiantBranch::CalculateRadiusOnPhase_Static(p_Mass, LZAHB, p_BnCoefficients);

    double mu    = (p_Mass - p_CoreMass) / (p_MHeF - p_CoreMass);
    double f     = ((1.0 + b[21]) * PPOW(mu, b[22])) / (1.0 + b[21] * PPOW(mu, b[23]));

    return ((1.0 - f)) * RZHe + (f * RGB);

#undef b
}


/*
 * Calculate radius at Helium Ignition
 *
 * Hurley et al. 2000, eq 50
 *
 *
 * double CalculateRadiusAtHeIgnition(const double p_Mass)
 *
 * @param   [IN]    p_Mass                      Mass in Msol
 * @return                                      Radius at Helium Ignition in Rsol
 */
double GiantBranch::CalculateRadiusAtHeIgnition(const double p_Mass) const {
#define massCutoffs(x) m_MassCutoffs[static_cast<int>(MASS_CUTOFF::x)]  // for convenience and readability - undefined at end of function

    double RHeI = 0.0;                                                      // Radius at Helium Ignition

    double LHeI      = CalculateLuminosityAtHeIgnition_Static(p_Mass, m_Alpha1, massCutoffs(MHeF), m_BnCoefficients);
    double RmHe      = CHeB::CalculateMinimumRadiusOnPhase_Static(p_Mass, m_CoreMass, m_Alpha1, massCutoffs(MHeF), massCutoffs(MFGB), m_MinimumLuminosityOnPhase, m_BnCoefficients);
    double RGB_LHeI  = CalculateRadiusOnPhase(p_Mass, LHeI);

    if (utils::Compare(p_Mass, massCutoffs(MFGB)) <= 0) {
        RHeI = RGB_LHeI;
    }
    else if (utils::Compare(p_Mass, std::max(massCutoffs(MFGB), HIGH_MASS_THRESHOLD)) >= 0) {
        double RAGB_LHeI = EAGB::CalculateRadiusOnPhase_Static(p_Mass, LHeI, massCutoffs(MHeF), m_BnCoefficients);
        RHeI             = std::min(RmHe, RAGB_LHeI);                        // Hurley et al. 2000, eq 55
    }
    else {
        double mu = log10(p_Mass / HIGH_MASS_THRESHOLD) / log10(massCutoffs(MFGB) / HIGH_MASS_THRESHOLD);
        RHeI      = RmHe * PPOW((RGB_LHeI / RmHe), mu);
    }

    return RHeI;

#undef massCutoffs
}


/*
 * Calculate radius of the remnant the star would become if it lost all of its
 * envelope immediately (i.e. M = Mc, coreMass)
 *
 * Hurley et al. 2000, just after eq 105
 *
 *
 * double CalculateRemnantRadius()
 *
 * @return                                      Radius of remnant core in Rsol
 */
double GiantBranch::CalculateRemnantRadius() const {
#define massCutoffs(x) m_MassCutoffs[static_cast<int>(MASS_CUTOFF::x)]  // for convenience and readability - undefined at end of function
    return (utils::Compare(m_Mass0, massCutoffs(MHeF)) > 0)
            ? HeMS::CalculateRadiusAtZAMS_Static(m_CoreMass)
            : WhiteDwarfs::CalculateRadiusOnPhase_Static(m_CoreMass);

#undef massCutoffs
}

/*
 * Calculate the radial extent of the convective outer envelope
 *
 * Combination of Hurley et al. 2000, end of sec. 7.2, and Hurley et al. 2002, sec. 2.3, particularly subsec. 2.3.1, eqs 39-40
 *
 *
 * double CalculateRadialExtentConvectiveEnvelope()
 *
 * @return                                      Radial extent of the convective outer envelope
 */
double GiantBranch::CalculateRadialExtentConvectiveEnvelope() const{
    double convectiveEnvelopeMass, convectiveEnvelopeMassMax;
    std::tie(convectiveEnvelopeMass, convectiveEnvelopeMassMax) = CalculateConvectiveEnvelopeMass();
    if (utils::Compare(convectiveEnvelopeMass, 0.0) <= 0 || utils::Compare(convectiveEnvelopeMassMax, 0.0) <= 0 ) return 0.0;   // massless convective envelope has zero radial extent
    
    double convectiveCoreMass   = CalculateConvectiveCoreMass();
    double convectiveCoreRadius = CalculateConvectiveCoreRadius();
    // assume that the final radiative intershell (if any) would have a density that is a geometric mean of the core density and total density
    double radiativeIntershellMass           = m_Mass - convectiveCoreMass - convectiveEnvelopeMassMax;
    double convectiveCoreRadiusCubed         = convectiveCoreRadius * convectiveCoreRadius * convectiveCoreRadius;
    double radiativeIntershellDensity        = 1.0 / (4.0 /3.0 * M_PI) * std::sqrt(convectiveCoreMass / convectiveCoreRadiusCubed * m_Mass / m_Radius / m_Radius / m_Radius);
    double outerEdgeRadiativeIntershellCubed = radiativeIntershellMass / (4.0 / 3.0 * M_PI * radiativeIntershellDensity) + convectiveCoreRadiusCubed;
    
    return std::sqrt(convectiveEnvelopeMass/convectiveEnvelopeMassMax) * (m_Radius - std::cbrt(outerEdgeRadiativeIntershellCubed));
}


///////////////////////////////////////////////////////////////////////////////////////
//                                                                                   //
//                                 MASS CALCULATIONS                                 //
//                                                                                   //
///////////////////////////////////////////////////////////////////////////////////////


/*
 * Calculate core mass at the Base of the Asymptotic Giant Branch
 *
 * Hurley et al. 2000, eq 66
 *
 *
 * double CalculateCoreMassAtBAGB(const double p_Mass)
 *
 * @param   [IN]    p_Mass                      Mass in Msol
 * @return                                      Core mass at the Base of the Asymptotic Giant Branch in Msol
 */
double GiantBranch::CalculateCoreMassAtBAGB(const double p_Mass) const {
#define b m_BnCoefficients  // for convenience and readability - undefined at end of function
    return std::sqrt(std::sqrt((b[36] * PPOW(p_Mass, b[37])) + b[38]));   // sqrt() is much faster than PPOW()
#undef b
}


/*
 * Calculate core mass at the Base of the Asymptotic Giant Branch
 *
 * Hurley et al. 2000, eq 66
 *
 * Static version required by CalculateGBParams_Static()
 * 
 *
 * double CalculateCoreMassAtBAGB_Static(const double p_Mass, const DBL_VECTOR &p_BnCoefficients)
 *
 * @param   [IN]    p_Mass                      Mass in Msol
 * @param   [IN]    p_BnCoefficients            b(n) coefficients
 * @return                                      Core mass at the Base of the Asymptotic Giant Branch in Msol
 */
double GiantBranch::CalculateCoreMassAtBAGB_Static(const double p_Mass, const DBL_VECTOR &p_BnCoefficients) {
#define b p_BnCoefficients  // for convenience and readability - undefined at end of function
    return std::sqrt(std::sqrt((b[36] * PPOW(p_Mass, b[37])) + b[38]));   // sqrt() is much faster than PPOW()
#undef b
}


/*
 * Calculate core mass at the Base of the Giant Branch
 *
 * Hurley et al. 2000, eq 44
 *
 * For large enough M, we have McBGB ~ 0.098*Mass**(1.35)
 *
 *
 * double CalculateCoreMassAtBGB(const double p_Mass, const DBL_VECTOR &p_GBParams)
 *
 * @param   [IN]    p_Mass                      Mass in Msol
 * @param   [IN]    p_GBParams                  Giant Branch parameters
 * @return                                      Core mass at the Base of the Giant Branch in Msol
 */
double GiantBranch::CalculateCoreMassAtBGB(const double p_Mass, const DBL_VECTOR &p_GBParams) {
#define gbParams(x) p_GBParams[static_cast<int>(GBP::x)]                // for convenience and readability - undefined at end of function
#define massCutoffs(x) m_MassCutoffs[static_cast<int>(MASS_CUTOFF::x)]  // for convenience and readability - undefined at end of function

    if (utils::Compare(p_Mass, massCutoffs(MHeF)) <= 0)     return 0.0;                                                 // No McBGB for stars with mass below the helium flash threshold, see text above Eq. (44) of Hurley+ (2000)
    
    double luminosity = GiantBranch::CalculateLuminosityAtPhaseBase_Static(massCutoffs(MHeF), m_AnCoefficients);
    double Mc_MHeF    = BaseStar::CalculateCoreMassGivenLuminosity_Static(luminosity, p_GBParams);
    double c          = (Mc_MHeF * Mc_MHeF * Mc_MHeF * Mc_MHeF) - (MC_L_C1 * PPOW(massCutoffs(MHeF), MC_L_C2));         // pow() is slow - use multiplication
    
    return std::min((0.95 * gbParams(McBAGB)), std::sqrt(std::sqrt(c + (MC_L_C1 * PPOW(p_Mass, MC_L_C2)))));            // sqrt is much faster than PPOW()

#undef massCutoffs
#undef gbParams
}


/*
 * Calculate core mass at the Base of the Giant Branch
 *
 * Hurley et al. 2000, eq 44
 *
 * For large enough M, we have McBGB ~ 0.098*Mass**(1.35)
 *
 * Static version required by CalculateGBParams_Static()
 *
 *
 * double CalculateCoreMassAtBGB_Static(const double p_Mass, const DBL_VECTOR &p_MassCutoffs, const DBL_VECTOR &p_AnCoefficients, const DBL_VECTOR &p_GBParams)
 *
 * @param   [IN]    p_Mass                      Mass in Msol
 * @param   [IN]    p_MassCutoffs               Mass cutoffs
 * @param   [IN]    p_AnCoefficients            a(n) coefficients
 * @param   [IN]    p_GBParams                  Giant Branch parameters
 * @return                                      Core mass at the Base of the Giant Branch in Msol
 */
double GiantBranch::CalculateCoreMassAtBGB_Static(const double      p_Mass, 
                                                  const DBL_VECTOR &p_MassCutoffs, 
                                                  const DBL_VECTOR &p_AnCoefficients, 
                                                  const DBL_VECTOR &p_GBParams) {
#define gbParams(x) p_GBParams[static_cast<int>(GBP::x)]                // for convenience and readability - undefined at end of function
#define massCutoffs(x) p_MassCutoffs[static_cast<int>(MASS_CUTOFF::x)]  // for convenience and readability - undefined at end of function

    double luminosity = GiantBranch::CalculateLuminosityAtPhaseBase_Static(massCutoffs(MHeF), p_AnCoefficients);
    double Mc_MHeF    = BaseStar::CalculateCoreMassGivenLuminosity_Static(luminosity, p_GBParams);
    double c          = (Mc_MHeF * Mc_MHeF * Mc_MHeF * Mc_MHeF) - (MC_L_C1 * PPOW(massCutoffs(MHeF), MC_L_C2));     // pow() is slow - use multiplication

    return std::min((0.95 * gbParams(McBAGB)), std::sqrt(std::sqrt(c + (MC_L_C1 * PPOW(p_Mass, MC_L_C2)))));        // sqrt is much faster than PPOW()

#undef massCutoffs
#undef gbParams
}


/*
 * Calculate the core mass at which the Asymptotic Giant Branch phase is terminated in a SN/loss of envelope
 *
 * Hurley et al. 2000, eq 75 -- but note we use MECS rather than MCH
 *
 *
 * double CalculateCoreMassAtSupernova_Static(const double p_McBAGB)
 *
 * @param   [IN]    p_McBAGB                    Core mass at the Base of the Asymptotic Giant Branch in Msol
 * @return                                      Maximum core mass before supernova on the Asymptotic Giant Branch
 */
double GiantBranch::CalculateCoreMassAtSupernova_Static(const double p_McBAGB) {
    return std::max(MECS, (0.773 * p_McBAGB) - 0.35);
}


/*
 * Calculate core mass at Helium Ignition
 *
 * Hurley et al. 2000, eq 37 for low mass stars (coreMass - Luminosity relation)
 *
 * For M >= X use Hurley et al. 2000, eq 44, replacing Mc(LBGB(MHeF)) with Mc(LHeI(MHeF))
 * (see end of section 5.3 for this note)
 *
 *
 * double CalculateCoreMassAtHeIgnition(const double p_Mass)
 *
 * @param   [IN]    p_Mass                      Mass in Msol
 * @return                                      Core mass at Helium Ignition in Msol
 */
double GiantBranch::CalculateCoreMassAtHeIgnition(const double p_Mass) const {
#define massCutoffs(x) m_MassCutoffs[static_cast<int>(MASS_CUTOFF::x)]  // for convenience and readability - undefined at end of function

    double coreMass;

    if (utils::Compare(p_Mass, massCutoffs(MHeF)) < 0) {
        double luminosity = CalculateLuminosityAtHeIgnition_Static(p_Mass, m_Alpha1, massCutoffs(MHeF), m_BnCoefficients);

        coreMass          = BaseStar::CalculateCoreMassGivenLuminosity_Static(luminosity, m_GBParams);
    }
    else {
        double luminosity_MHeF = CalculateLuminosityAtHeIgnition_Static(massCutoffs(MHeF), m_Alpha1, massCutoffs(MHeF), m_BnCoefficients);
        double Mc_MHeF         = BaseStar::CalculateCoreMassGivenLuminosity_Static(luminosity_MHeF, m_GBParams);
        double McBAGB          = CalculateCoreMassAtBAGB(p_Mass);
        double c               = (Mc_MHeF * Mc_MHeF * Mc_MHeF * Mc_MHeF) - (MC_L_C1 * PPOW(massCutoffs(MHeF), MC_L_C2));    // pow() is slow - use multiplication

        coreMass               = std::min((0.95 * McBAGB), std::sqrt(std::sqrt(c + (MC_L_C1 * PPOW(p_Mass, MC_L_C2)))));    // sqrt() is much faster than PPOW()
    }

    return coreMass;

#undef massCutoffs
}


/*
 * Calculate the mass of the core during second dredge up
 *
 * Hurley et al. 2000, eq 69 (not numbered in the paper, but between eqs 68 & 70)
 *
 *
 * double CalculateCoreMassAt2ndDredgeUp(const double p_McBAGB)
 *
 * @param       [IN]    p_McBAGB                Core mass at the Base of the Asymptotic Giant Branch
 * @return                                      Core Mass at second dredge up (McDU)
 */
double GiantBranch::CalculateCoreMassAt2ndDredgeUp_Static(const double p_McBAGB) {
    return utils::Compare(p_McBAGB, 0.8) >= 0 ? ((0.44 * p_McBAGB) + 0.448) : p_McBAGB;
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
double GiantBranch::CalculateMassLossRateHurley() {

    double rateNJ = CalculateMassLossRateNieuwenhuijzenDeJager();
    double rateKR = CalculateMassLossRateKudritzkiReimers();
    double rateWR = CalculateMassLossRateWolfRayet(m_Mu);

    m_DominantMassLossRate = MASS_LOSS_TYPE::GB;
    double dominantRate    = std::max(rateNJ, rateKR);

    if (utils::Compare(rateWR, dominantRate) > 0) {
        dominantRate           = rateWR;
        m_DominantMassLossRate = MASS_LOSS_TYPE::WR;
    }

    return dominantRate;
}


/*
 * Determines if mass transfer is unstable according to the critical mass ratio.
 *
 * See e.g de Mink et al. 2013, Claeys et al. 2014, and Ge et al. 2010, 2015, 2020 for discussions.
 *
 * Assumes this star is the donor; relevant accretor details are passed as parameters.
 * Critical mass ratio is defined as qCrit = mAccretor/mDonor.
 *
 * double GiantBranch::CalculateCriticalMassRatioClaeys14(const bool p_AccretorIsDegenerate) 
 *
 * @param   [IN]    p_AccretorIsDegenerate      Boolean indicating if accretor in degenerate (true = degenerate)
 * @return                                      Critical mass ratio for unstable MT 
 */
double GiantBranch::CalculateCriticalMassRatioClaeys14(const bool p_AccretorIsDegenerate) const {

    double qCrit;
                                                                                                                            
    if (p_AccretorIsDegenerate) {                                                                                       // degenerate accretor
        qCrit = OPTIONS->MassTransferCriticalMassRatioGiantDegenerateAccretor();
    }
    else {                                                                                                              // non-degenerate accretor 
        qCrit = OPTIONS->MassTransferCriticalMassRatioGiantNonDegenerateAccretor();
        if (qCrit < 0.0) {                                                                                              // default value of -1 recalculates qCrit with the following function 
            double coreMassRatio   = m_HeCoreMass / m_Mass;
            double coreMassRatio_4 = coreMassRatio * coreMassRatio * coreMassRatio * coreMassRatio * coreMassRatio;
            double x               = BaseStar::CalculateGBRadiusXExponent();                                            // x from Hurley et al 2000, Eq. 47 - Depends on composition
            qCrit                  = 2.13 / ( 1.67 - x + 2.0 * coreMassRatio_4);                                        // Claeys+ 2014, Table 2
        }
    }

    return qCrit;
}


/*
 * Determines if mass transfer is unstable according to the critical mass ratio.
 *
 * See Hurley et al. 2002 sect. 2.6.1.
 *
 * Assumes this star is the donor.
 * Critical mass ratio is defined as qCrit = mAccretor/mDonor.
 *
 * double GiantBranch::CalculateCriticalMassRatioHurleyHjellmingWebbink() 
 *
 * @return                                      Critical mass ratio for unstable MT 
 */
double GiantBranch::CalculateCriticalMassRatioHurleyHjellmingWebbink() const {
    
    double qCrit = 0.362 + 1.0 / (3.0 * (1.0 - CoreMass() / Mass()));   // Defined as mDonor/mAccretor in Hurley et al. 2002, equation found after eq. 57 (no label)
    
    return 1.0 / qCrit;
}


/*
 * Calculate the Adiabatic Exponent (for convective-envelope giant-like stars)
 *
 *
 * double CalculateZetaConvectiveEnvelopeGiant(ZETA_PRESCRIPTION p_ZetaPrescription) const
 *
 * @param   [IN]    p_ZetaPrescription          Prescription for computing ZetaStar
 * @return                                      Adiabatic exponent
 */
double GiantBranch::CalculateZetaConvectiveEnvelopeGiant(ZETA_PRESCRIPTION p_ZetaPrescription) {
    
    double zeta = 0.0;                                                                                          // default value

    switch (p_ZetaPrescription) {                                                                               // which prescription?
            
        case ZETA_PRESCRIPTION::HURLEY:                                                                         // HURLEY: Hurley, Tout, and Pols, 2002, eq 56
            zeta = CalculateZetaAdiabaticHurley2002(m_CoreMass);
            break;
            
        case ZETA_PRESCRIPTION::ARBITRARY:                                                                      // ARBITRARY: user program options thermal zeta value
            zeta = OPTIONS->ZetaAdiabaticArbitrary();
            break;

        case ZETA_PRESCRIPTION::SOBERMAN:                                                                       // SOBERMAN: Soberman, Phinney, and van den Heuvel, 1997, eq 61
            zeta = CalculateZetaAdiabaticSPH(m_CoreMass);
            break;

        default:                                                                                                // unknown prescription
            // the only ways this can happen are if someone added a ZETA_PRESCRIPTION
            // and it isn't accounted for in this code, or if there is a defect in the code that causes
            // this function to be called with a bad parameter.  We should not default here, with or without a
            // warning - this is a code defect, so we flag it as an error and that will result in termination of
            // the evolution of the star or binary.
            // The correct fix for this is to add code for the missing prescription or, if the missing
            // prescription is superfluous, remove it from ZETA_PRESCRIPTION - or find and fix the code defect.

            THROW_ERROR(ERROR::UNKNOWN_ZETA_PRESCRIPTION);                                                      // throw error             
    }
    
    return zeta;
}



/*
 * Calculate the mass-radius response exponent Zeta
 *
 * Hurley et al. 2000, eqs 97 & 98
 *
 *
 * double CalculateZetaConstantsByEnvelope(ZETA_PRESCRIPTION p_ZetaPrescription)
 *
 * @param   [IN]    p_ZetaPrescription          Prescription for computing ZetaStar
 * @return                                      mass-radius response exponent Zeta
 */
double GiantBranch::CalculateZetaConstantsByEnvelope(ZETA_PRESCRIPTION p_ZetaPrescription) {
    
    double zeta = 0.0;                                                                              // default value
    
    switch (DetermineEnvelopeType()) {                                                              // which envelope type?
        case ENVELOPE::RADIATIVE:
            zeta = OPTIONS->ZetaRadiativeEnvelopeGiant();
            break;
            
        case ENVELOPE::CONVECTIVE:
            zeta = CalculateZetaConvectiveEnvelopeGiant(p_ZetaPrescription);
            break;

        case ENVELOPE::REMNANT:
            zeta = 0.0;
            break;

        default:                                                                                    // unknown envelope type
            // the only way this can happen is if someone added a new envelope type (to ENVELOPE)
            // and it isn't accounted for in this code.  We should not default here, with or without a warning.
            // We are here because DetermineEnvelopeType() returned an envelope type this code doesn't account
            // for, and that should be flagged as an error and result in termination of the evolution of the
            // star or binary.
            // The correct fix for this is to add code for the missing envelope type or, if the missing envelope
            // type is incorrect/superfluous, remove it from ENVELOPE.

            THROW_ERROR(ERROR::UNKNOWN_ENVELOPE_TYPE);                                              // throw error
    }
    
    return zeta;
}


/*
 * Approximates the mass of the outer convective envelope.
 *
 * This is needed for the Hirai & Mandel (2022) two-stage CE formalism.
 * Follows the fits of Picker, Hirai, Mandel (2024), arXiv:2402.13180
 *
 *
 * std::tuple<double, double> GiantBranch::CalculateConvectiveEnvelopeMass()
 *
 * @return                                      Tuple containing the mass of the outer convective envelope and its maximum value
 */
DBL_DBL GiantBranch::CalculateConvectiveEnvelopeMass() const {
    
    double MinterfMcoref = -0.023 * m_Log10Metallicity - 0.0023;                                                            // eq. (8) of Picker+ 2024

    // We need the temperature of the star just after BAGB, which is the temperature at the
    // start of the EAGB phase.  Since we are on the giant branch here, we can clone this
    // object as an EAGB object and, as long as it is initialised (to the start of the phase),
    // we can query the cloned object for its temperature.
    //
    // To ensure the clone does not participate in logging, we set its persistence to EPHEMERAL.
    //
    // Furthermore, 'this' is const in this function, so we first remove its const-ness (required
    // to call Clone()) via the use of const_cast<>().  Since we don't know what class the
    // underlying object is, we cast it to EAGB&.

    EAGB *clone = EAGB::Clone(static_cast<EAGB&>(const_cast<GiantBranch&>(*this)), OBJECT_PERSISTENCE::EPHEMERAL);
    clone->UpdateAttributesAndAgeOneTimestep(0.0, 0.0, 0.0, true);                                                          // otherwise, temperature not updated
    double Tmin = clone->Temperature();                                                                                     // get temperature of clone
    delete clone; clone = nullptr;                                                                                          // return the memory allocated for the clone
    
    // Use Eq. 6 of Mandel, Hirai, Picker (2024) rather than Eq. 6 of Picker+ 2024 for Tonset to avoid issues caused by
    // differences between temperatures in MESA models (used in Picker+ fits) and Pols models (used in Hurley+ SSE tracks)
    double Tonset     = Tmin / std::min(0.0695 - 0.057 * m_Log10Metallicity, 0.95);                                         // eq. (6) of Mandel, Hirai, Picker, 2024
    
    double mCoreFinal = CalculateCoreMassAtBAGB(m_Mass0);
    double mConvMax   = std::max(m_Mass - mCoreFinal * (1.0 + MinterfMcoref), 0.0);                                         // eq. (9) of Picker+ 2024

    // Picker+ 2024 fits were only made for stars above 8.0 solar masses, with runs down to 5.0 solar masses, 
    // so using the final core mass as an approximate threshold of validity
    if(utils::Compare(mCoreFinal, 1.5) < 0) mConvMax = std::max(m_Mass - mCoreFinal, 0.0);                                  // unlike massive stars, intermediate-mass stars have almost no radiative intershell at maximum convective envelope extent
    
    double convectiveEnvelopeMass = mConvMax / (1.0 + exp(4.6 * (Tmin + Tonset - 2.0 * m_Temperature) / (Tmin - Tonset)));  // eq. (7) of Picker+ 2024
    
    return std::tuple<double, double> (convectiveEnvelopeMass, mConvMax);
}


///////////////////////////////////////////////////////////////////////////////////////
//                                                                                   //
//                            LIFETIME / AGE CALCULATIONS                            //
//                                                                                   //
///////////////////////////////////////////////////////////////////////////////////////


/*
 * Calculate the time to Helium ignition tHeI
 *
 * Hurley et al. 2000, eq 43
 *
 *
 * double CalculateLifetimeToHeIgnition(const double p_Mass, const double p_Tinf1_FGB, const double p_Tinf2_FGB)
 *
 *
 * @param   [IN]    p_Mass                      Mass in Msol
 * @param   [IN]    p_Tinf1_FGB                 tinf1_FGB (Timescales[TIMESCALE::tinf1_FGB]) - First Giant Branch tinf1
 * @param   [IN]    p_Tinf2_FGB                 tinf2_FGB (Timescales[TIMESCALE::tinf2_FGB]) - First Giant Branch tinf2
 * @return                                      Lifetime to He ignition (tHeI)
 */
double GiantBranch::CalculateLifetimeToHeIgnition(const double p_Mass, const double p_Tinf1_FGB, const double p_Tinf2_FGB) {
#define gbParams(x) m_GBParams[static_cast<int>(GBP::x)]                // for convenience and readability - undefined at end of function
#define massCutoffs(x) m_MassCutoffs[static_cast<int>(MASS_CUTOFF::x)]  // for convenience and readability - undefined at end of function

    double LHeI = CalculateLuminosityAtHeIgnition_Static(p_Mass, m_Alpha1, massCutoffs(MHeF), m_BnCoefficients);
    double p1   = gbParams(p) - 1.0;
    double q1   = gbParams(q) - 1.0;

    return utils::Compare(LHeI, gbParams(Lx)) <= 0
            ? p_Tinf1_FGB - (1.0 / (p1 * gbParams(AH) * gbParams(D))) * PPOW((gbParams(D) / LHeI), (p1 / gbParams(p)))
            : p_Tinf2_FGB - (1.0 / (q1 * gbParams(AH) * gbParams(B))) * PPOW((gbParams(B) / LHeI), (q1 / gbParams(q)));

#undef massCutoffs
#undef gbParams
}


///////////////////////////////////////////////////////////////////////////////////////
//                                                                                   //
//                               ROTATION CALCULATIONS                               //
//                                                                                   //
///////////////////////////////////////////////////////////////////////////////////////

/*
 * Calculate moment of inertia
 *
 * Hurley et al., 2000, paragraph immediately following eq 109 
 *
 * 
 * double GiantBranch::CalculateMomentOfInertia()
 * 
 * @return                                      Moment of inertia (Msol AU^2)
 */
double GiantBranch::CalculateMomentOfInertia() const {
    double Rc = CalculateRemnantRadius();
    
    return (0.1 * (m_Mass - m_CoreMass) * m_Radius * m_Radius) + (0.21 * m_CoreMass * Rc * Rc);
}



///////////////////////////////////////////////////////////////////////////////////////
//                                                                                   //
//                              SUPERNOVA CALCULATIONS                               //
//                                                                                   //
///////////////////////////////////////////////////////////////////////////////////////


/*
 * Calculate remnant type given COCoreMass
 *
 * Muller et al. 2016 as presented in appendix B of Vigna-Gomez et al. 2018 (arXiv:1805.07974)
 *
 *
 * STELLAR_TYPE CalculateRemnantTypeByMuller2016(const double p_COCoreMass)
 *
 * @param   [IN]    p_COCoreMass                COCoreMass in Msol
 * @return                                      Remnant type (stellar type)
 */
STELLAR_TYPE GiantBranch::CalculateRemnantTypeByMuller2016(const double p_COCoreMass) {

    STELLAR_TYPE stellarType;

         if (utils::Compare(p_COCoreMass, 3.6 ) < 0) { stellarType = STELLAR_TYPE::NEUTRON_STAR; }
    else if (utils::Compare(p_COCoreMass, 4.05) < 0) { stellarType = STELLAR_TYPE::BLACK_HOLE; }
    else if (utils::Compare(p_COCoreMass, 4.6 ) < 0) { stellarType = STELLAR_TYPE::NEUTRON_STAR; }
    else if (utils::Compare(p_COCoreMass, 5.7 ) < 0) { stellarType = STELLAR_TYPE::BLACK_HOLE; }
    else if (utils::Compare(p_COCoreMass, 6.0 ) < 0) { stellarType = STELLAR_TYPE::NEUTRON_STAR; }
    else                                             { stellarType = STELLAR_TYPE::BLACK_HOLE; }

    return stellarType;
}


/*
 * Calculate remnant mass given COCoreMass according to the Schneider et al. 2020 prescription (arxiv:2008.08599)
 *
 * Note that Schneider only prescribes remnant masses for the simple cases of single episode
 * Mass Transfer, so some of the double episode cases here are a bit uncertain, and may
 * need to be refined at a later date.
 *
 * double CalculateRemnantMassBySchneider2020(const double p_COCoreMass, const bool p_UseSchneiderAlt)
 *
 * @param   [IN]    p_COCoreMass                COCoreMass in Msol
 * @param   [IN]    p_UseSchneiderAlt           Whether to use the Schneider alt prescription 
 * @return                                      Remnant mass in Msol
 */
double GiantBranch::CalculateRemnantMassBySchneider2020(const double p_COCoreMass, const bool p_UseSchneiderAlt) {

    double logRemnantMass;

    ST_VECTOR mtHist                    = MassTransferDonorHistory();                                                   // mass transfer history vector
    MT_CASE   massTransferCase = MT_CASE::OTHER;

    // determine MT history - this will tell us which Schneider MT case prescription should be used
    if (mtHist.size() == 0) {                                                                                           // no history of MT - effectively single star
        massTransferCase = MT_CASE::NONE;
    }
    else {                                                                                                              // star was MT donor at least once
        // determine MT_CASE of most recent MT event
        STELLAR_TYPE mostRecentDonorType = mtHist[mtHist.size() - 1];                                                   // stellar type at most recent MT event (as donor)
        BaseStar* newStar                = stellarUtils::NewStar(mostRecentDonorType);                                  // create new (empty) star of correct stellar type
        massTransferCase        = newStar->DetermineMassTransferTypeAsDonor();                                          // get MT type as donor
        delete newStar; newStar = nullptr;                                                                              // return the memory allocated for the new star
    }

    // apply the appropriate remnant mass prescription for the chosen MT case
    switch (massTransferCase) {                                                                                         // which MT_CASE?

        case MT_CASE::NONE:                                                                                             // no history of MT
        case MT_CASE::OTHER:                                                                                            // if MT happens from naked He stars, WDs, etc., assume that the core properties are not affected
            if (!p_UseSchneiderAlt) {                                                                                   // use standard or alternative remnant mass prescription for effectively single stars?
                     // standard prescription
                     if (utils::Compare(p_COCoreMass, 6.357)  < 0) { logRemnantMass = log10(0.03357 * p_COCoreMass + 1.31780); }
                else if (utils::Compare(p_COCoreMass, 7.311)  < 0) { logRemnantMass = -0.02466 * p_COCoreMass + 1.28070; }
                else if (utils::Compare(p_COCoreMass, 12.925) < 0) { logRemnantMass = log10(0.03357 * p_COCoreMass + 1.31780); }
                else                                               { logRemnantMass = 0.01940 * p_COCoreMass + 0.98462; }
            }
            else {  
                     // alternative prescription
                     if (utils::Compare(p_COCoreMass, 6.357)  < 0) { logRemnantMass = log10(0.04199 * p_COCoreMass + 1.28128); }
                else if (utils::Compare(p_COCoreMass, 7.311)  < 0) { logRemnantMass = -0.02466 * p_COCoreMass + 1.28070; }
                else if (utils::Compare(p_COCoreMass, 12.925) < 0) { logRemnantMass = log10(0.04701 * (p_COCoreMass * p_COCoreMass) - 0.91403 * p_COCoreMass + 5.93380); }
                else                                               { logRemnantMass = 0.01940*p_COCoreMass + 0.98462; }
            }
            break;

        case MT_CASE::A:                                                                                                // case A MT

                 if (utils::Compare(p_COCoreMass, 7.064)  < 0) { logRemnantMass = log10(0.02128 * p_COCoreMass + 1.35349); }
            else if (utils::Compare(p_COCoreMass, 8.615)  < 0) { logRemnantMass = 0.03866 * p_COCoreMass + 0.64417; }
            else if (utils::Compare(p_COCoreMass, 15.187) < 0) { logRemnantMass = log10(0.02128 * p_COCoreMass + 1.35349); }
            else                                               { logRemnantMass = 0.02573 * p_COCoreMass + 0.79027; }
            break;

        case MT_CASE::B:                                                                                                // case B MT

                 if (utils::Compare(p_COCoreMass, 7.548)  < 0) { logRemnantMass = log10(0.01909 * p_COCoreMass + 1.34529); }
            else if (utils::Compare(p_COCoreMass, 8.491)  < 0) { logRemnantMass = 0.03306 * p_COCoreMass + 0.68978; }
            else if (utils::Compare(p_COCoreMass, 15.144) < 0) { logRemnantMass = log10(0.01909 * p_COCoreMass + 1.34529); }
            else                                               { logRemnantMass = 0.02477 * p_COCoreMass + 0.80614; }
            break;

        case MT_CASE::C:                                                                                                // case C MT

                 if (utils::Compare(p_COCoreMass, 6.357)  < 0) { logRemnantMass = log10(0.03781 * p_COCoreMass + 1.36363); }
            else if (utils::Compare(p_COCoreMass, 7.311)  < 0) { logRemnantMass = 0.05264 * p_COCoreMass + 0.58531; }
            else if (utils::Compare(p_COCoreMass, 14.008) < 0) { logRemnantMass = log10(0.03781 * p_COCoreMass + 1.36363); }
            else                                               { logRemnantMass = 0.01795 * p_COCoreMass + 0.98797; }
            break;

        default:                                                                                                        // unknown MT_CASE
            // the only way this can happen is if someone added an MT_CASE
            // and it isn't accounted for in this code.  We should not default here, with or without a warning.
            // We are here because DetermineMassTransferTypeAsDonor() returned an MT_CASE this code doesn't
            // account for, and that should be flagged as an error and result in termination of the evolution
            // of the star or binary.
            // The correct fix for this is to add code for the missing MT_CASE or, if the missing MT_CASE is
            // incorrect/superfluous, remove it from the possible MT_CASE values.

            THROW_ERROR(ERROR::UNKNOWN_MT_CASE);                                                                        // throw error   
    }
    
    // convert to linear value, and limit to the pre-SN He Core mass
    return std::min(PPOW(10.0, logRemnantMass), m_SupernovaDetails.HeCoreMassAtCOFormation);
}


/*
 * Calculate remnant mass according to the Maltsev et al. 2024 prescription
 *
 * Note that Maltsev2024 prescription is based on the donor type at the first MT event
 *
 * double CalculateRemnantMassByMaltsev2024(const double p_COCoreMass)
 *
 * @param   [IN]    p_COCoreMass                COCoreMass in Msol
 * @param   [IN]    p_HeCoreMass                HeCoreMass in Msol
 * @return                                      Remnant mass in Msol
 */
double GiantBranch::CalculateRemnantMassByMaltsev2024(const double p_COCoreMass, const double p_HeCoreMass) {

    ST_VECTOR mtHist           = MassTransferDonorHistory();                                                            // mass transfer history vector
    MT_CASE   massTransferCase = MT_CASE::OTHER;
    double    log10Z           = m_Log10Metallicity - LOG10_ZSOL;                                                       // log_{10} (Z/Zsol), for convenience
    double    M1, M2, M3;

    if (utils::Compare(p_COCoreMass, MALTSEV2024_MMIN) < 0)                                                             // NS formation regardless of metallicity and MT history
        return CalculateRemnantNSMassMullerMandel(p_COCoreMass, p_HeCoreMass);
    
    if (utils::Compare(p_COCoreMass, MALTSEV2024_MMAX) > 0)                                                             // BH formation regardless of metallicity and MT history
        return p_HeCoreMass;
    
    // determine MT history - this will tell us which Schneider MT case prescription should be used
    if (mtHist.size() == 0) {                                                                                           // no history of MT - effectively single star
        massTransferCase = MT_CASE::NONE;
    }
    else {                                                                                                              // star was MT donor at least once
        // determine MT_CASE of first MT event
        STELLAR_TYPE mostRecentDonorType = mtHist[0];                                                                   // stellar type at first MT event (as donor)
        BaseStar* newStar                = stellarUtils::NewStar(mostRecentDonorType);                                  // create new (empty) star of correct stellar type
        massTransferCase        = newStar->DetermineMassTransferTypeAsDonor();                                          // get MT type as donor
        delete newStar; newStar = nullptr;                                                                              // return the memory allocated for the new star
    }

    // apply the appropriate remnant mass prescription for the chosen MT case
    switch (massTransferCase) {                                                                                         // which MT_CASE?

        case MT_CASE::NONE:                                                                                             // no history of MT
        case MT_CASE::OTHER:                                                                                            // if MT happens from naked He stars, WDs, etc., assume that the core properties are not affected
            M1 = MALTSEV2024_M1S + (MALTSEV2024_M1S - MALTSEV2024_M1SZ01) * log10Z;
            M2 = MALTSEV2024_M2S + (MALTSEV2024_M2S - MALTSEV2024_M2SZ01) * log10Z;
            M3 = MALTSEV2024_M3S + (MALTSEV2024_M3S - MALTSEV2024_M3SZ01) * log10Z;
            break;

        case MT_CASE::A:                                                                                                // case A MT
            M1 = MALTSEV2024_M1A + (MALTSEV2024_M1A - MALTSEV2024_M1AZ01) * log10Z;
            M2 = MALTSEV2024_M2A + (MALTSEV2024_M2A - MALTSEV2024_M2AZ01) * log10Z;
            M3 = MALTSEV2024_M3A + (MALTSEV2024_M3A - MALTSEV2024_M3AZ01) * log10Z;
            break;

        case MT_CASE::B:                                                                                                // case B MT
            M1 = MALTSEV2024_M1B + (MALTSEV2024_M1B - MALTSEV2024_M1BZ01) * log10Z;
            M2 = MALTSEV2024_M2B + (MALTSEV2024_M2B - MALTSEV2024_M2BZ01) * log10Z;
            M3 = MALTSEV2024_M3B + (MALTSEV2024_M3B - MALTSEV2024_M3BZ01) * log10Z;
            break;

        case MT_CASE::C:                                                                                                // case C MT
            M1 = MALTSEV2024_M1C + (MALTSEV2024_M1C - MALTSEV2024_M1CZ01) * log10Z;
            M2 = MALTSEV2024_M2C + (MALTSEV2024_M2C - MALTSEV2024_M2CZ01) * log10Z;
            M3 = MALTSEV2024_M3C + (MALTSEV2024_M3C - MALTSEV2024_M3CZ01) * log10Z;
            break;

        default:                                                                                                        // unknown MT_CASE
            // the only way this can happen is if someone added an MT_CASE
            // and it isn't accounted for in this code.  We should not default here, with or without a warning.
            // We are here because DetermineMassTransferTypeAsDonor() returned an MT_CASE this code doesn't
            // account for, and that should be flagged as an error and result in termination of the evolution
            // of the star or binary.
            // The correct fix for this is to add code for the missing MT_CASE or, if the missing MT_CASE is
            // incorrect/superfluous, remove it from the possible MT_CASE values.

            THROW_ERROR(ERROR::UNKNOWN_MT_CASE);                                                                        // throw error
    }
    
    if( utils::Compare(p_COCoreMass, M3) >=0 || (utils::Compare(p_COCoreMass, M1) >= 0 && utils::Compare(p_COCoreMass, M2) <= 0) )              // Complete fallback into BH
        return p_HeCoreMass;
    else if ( utils::Compare(p_COCoreMass, M2) > 0 && utils::Compare(p_COCoreMass, M3) < 0 && utils::Compare(RAND->Random(0, 1), 0.1) <= 0 )    // Partial fallback BH formation
        return CalculateFallbackBHMassMullerMandel(p_COCoreMass, p_HeCoreMass);
    return CalculateRemnantNSMassMullerMandel(p_COCoreMass, p_HeCoreMass);
}



/*
 * Calculate remnant mass given COCoreMass and HeCoreMass in the Mandel & Mueller (2020) prescription
 *
 * Mandel & Mueller, 2020
 *
 *
 * double CalculateRemnantMassByMullerMandel (const double p_COCoreMass, const double p_HeCoreMass)
 *
 * @param   [IN]    p_COCoreMass                COCoreMass in Msol
 * @param   [IN]    p_HeCoreMass                HeCoreMass in Msol
 * @return                                      Remnant mass in Msol
 */
double GiantBranch::CalculateRemnantMassByMullerMandel(const double p_COCoreMass, const double p_HeCoreMass) {

    double remnantMass       = 0.0;   
    double pBH               = 0.0;
    double pCompleteCollapse = 0.0;
   
    if (utils::Compare(p_COCoreMass, MULLERMANDEL_M1) < 0 || utils::Compare(p_HeCoreMass, OPTIONS->MaximumNeutronStarMass()) <= 0 )
	    pBH = 0.0;
    else if (utils::Compare(p_COCoreMass, MULLERMANDEL_M3) < 0)
    	pBH = 1.0 / (MULLERMANDEL_M3-MULLERMANDEL_M1) * (p_COCoreMass-MULLERMANDEL_M1);
    else
	    pBH = 1.0;

    if (utils::Compare(RAND->Random(0, 1), pBH) < 0) {  // this is a BH
        if (utils::Compare(p_COCoreMass, MULLERMANDEL_M4) < 0)
		    pCompleteCollapse = 1.0 / (MULLERMANDEL_M4 - MULLERMANDEL_M1) * (p_COCoreMass - MULLERMANDEL_M1);
        else
		    pCompleteCollapse = 1.0;

	    if (utils::Compare(RAND->Random(0, 1), pCompleteCollapse) < 0)
		    remnantMass = p_HeCoreMass;
	    else
            remnantMass = CalculateFallbackBHMassMullerMandel(p_COCoreMass, p_HeCoreMass);
    }
    else                                               // this is an NS
        remnantMass = CalculateRemnantNSMassMullerMandel(p_COCoreMass, p_HeCoreMass);

    return remnantMass;
}

/*
 * Calculate NS remnant mass given COCoreMass and HeCoreMass
 *
 * Mandel & Mueller, 2020; also used by Maltsev2024, so separated out into helper function
 *
 *
 * double CalculateRemnantNSMassByMullerMandel (const double p_COCoreMass, const double p_HeCoreMass)
 *
 * @param   [IN]    p_COCoreMass                COCoreMass in Msol
 * @param   [IN]    p_HeCoreMass                HeCoreMass in Msol
 * @return                                      Remnant mass in Msol
 */
double GiantBranch::CalculateRemnantNSMassMullerMandel(const double p_COCoreMass, const double p_HeCoreMass) {

    double remnantMass      = 0.0;
    std::size_t iterations  = 0;

    if (utils::Compare(p_COCoreMass, MULLERMANDEL_M1) < 0) {
        while (iterations++ < MULLERMANDEL_REMNANT_MASS_MAX_ITERATIONS            &&
               (utils::Compare(remnantMass, MULLERMANDEL_MINNS) < 0                ||
                utils::Compare(remnantMass, OPTIONS->MaximumNeutronStarMass()) > 0 ||
                utils::Compare(remnantMass, p_HeCoreMass) > 0)) {
            remnantMass = MULLERMANDEL_MU1 + RAND->RandomGaussian(MULLERMANDEL_SIGMA1);
        }
        if (iterations >= MULLERMANDEL_REMNANT_MASS_MAX_ITERATIONS) // failure to find a solution implies a narrow range; just pick a midpoint in this case
            remnantMass = (std::min(OPTIONS->MaximumNeutronStarMass(), p_HeCoreMass) + MULLERMANDEL_MINNS) / 2.0;
    }
    else if (utils::Compare(p_COCoreMass, MULLERMANDEL_M2) < 0) {
        while (iterations++ < MULLERMANDEL_REMNANT_MASS_MAX_ITERATIONS            &&
               (utils::Compare(remnantMass, MULLERMANDEL_MINNS) < 0                ||
                utils::Compare(remnantMass, OPTIONS->MaximumNeutronStarMass()) > 0 ||
                utils::Compare(remnantMass, p_HeCoreMass) > 0)) {
            remnantMass = MULLERMANDEL_MU2A + MULLERMANDEL_MU2B / (MULLERMANDEL_M2 - MULLERMANDEL_M1) * (p_COCoreMass - MULLERMANDEL_M1) + RAND->RandomGaussian(MULLERMANDEL_SIGMA2);
        }
        if (iterations >= MULLERMANDEL_REMNANT_MASS_MAX_ITERATIONS) // failure to find a solution implies a narrow range; just pick a midpoint in this case
            remnantMass = (std::min(OPTIONS->MaximumNeutronStarMass(), p_HeCoreMass) + MULLERMANDEL_MINNS) / 2.0;
    }
    else {
        while (iterations++ < MULLERMANDEL_REMNANT_MASS_MAX_ITERATIONS            &&
               (utils::Compare(remnantMass, MULLERMANDEL_MINNS) < 0                ||
                utils::Compare(remnantMass, OPTIONS->MaximumNeutronStarMass()) > 0 ||
                utils::Compare(remnantMass, p_HeCoreMass) > 0)) {
            remnantMass = MULLERMANDEL_MU3A + MULLERMANDEL_MU3B / (MULLERMANDEL_M3 - MULLERMANDEL_M2) * (p_COCoreMass - MULLERMANDEL_M2) + RAND->RandomGaussian(MULLERMANDEL_SIGMA3);
        }
        if (iterations >= MULLERMANDEL_REMNANT_MASS_MAX_ITERATIONS) // failure to find a solution implies a narrow range; just pick a midpoint in this case
            remnantMass = (std::min(OPTIONS->MaximumNeutronStarMass(), p_HeCoreMass) + MULLERMANDEL_MINNS) / 2.0;
    }
    return remnantMass;
}


/*
 * Calculate fallback BH mass given COCoreMass, HeCoreMass
 *
 * Mandel & Mueller, 2020; also used by Maltsev2024, so separated out into helper function
 *
 *
 * double CalculateFallbackBHMassByMullerMandel (const double p_COCoreMass, const double p_HeCoreMass)
 *
 * @param   [IN]    p_COCoreMass                COCoreMass in Msol
 * @param   [IN]    p_HeCoreMass                HeCoreMass in Msol
 * @return                                      Remnant mass in Msol
 */
double GiantBranch::CalculateFallbackBHMassMullerMandel(const double p_COCoreMass, const double p_HeCoreMass) {
    
    double remnantMass      = 0.0;
    std::size_t iterations  = 0;
    
    while (iterations++ < MULLERMANDEL_REMNANT_MASS_MAX_ITERATIONS &&
           (utils::Compare(remnantMass, OPTIONS->MaximumNeutronStarMass()) < 0 || utils::Compare(remnantMass, p_HeCoreMass) > 0)) {
        remnantMass = MULLERMANDEL_MUBH * p_COCoreMass + RAND->RandomGaussian(MULLERMANDEL_SIGMABH);
    }
    
    if (iterations >= MULLERMANDEL_REMNANT_MASS_MAX_ITERATIONS) // failure to find a solution implies a narrow range; just pick a midpoint in this case
        remnantMass = (OPTIONS->MaximumNeutronStarMass() + p_HeCoreMass) / 2.0;
    
    return remnantMass;
}


/*
 * Calculate remnant mass given Mass and COCoreMass per Muller et al. 2016 as presented in eq. B4 of Vigna-Gomez et al. 2018 (arXiv:1805.07974)
 *
 * double CalculateRemnantMassByMuller2016(const double p_Mass, const double p_COCoreMass)
 *
 * @param   [IN]    p_Mass                      Mass in Msol
 * @param   [IN]    p_COCoreMass                COCoreMass in Msol
 * @return                                      Remnant mass in Msol
 */
double GiantBranch::CalculateRemnantMassByMuller2016(const double p_Mass, const double p_COCoreMass) {
    double	remnantMass; 					                                                                        // Limit mass for a White Dwarf units Msun.

    if (utils::Compare(p_COCoreMass, 1.372) < 0) {
        // Not explicitly pointed out in Appendix B of Vigna-Gomez+2018 but assumed for continuity and simplicity
        // Muller+2016 didn't go as low as this in CO Core mass (see Figure A1 in that paper)
        remnantMass = 1.21;                         
    }
	else if (utils::Compare(p_COCoreMass, 1.49) < 0) { remnantMass = 1.21 - (0.4  * (p_COCoreMass - 1.372)); }
	else if (utils::Compare(p_COCoreMass, 1.65) < 0) { remnantMass = 1.16;                                   }
    else if (utils::Compare(p_COCoreMass, 2.4 ) < 0) { remnantMass = 1.32 + (0.3  * (p_COCoreMass - 1.65));  }
    else if (utils::Compare(p_COCoreMass, 3.2 ) < 0) { remnantMass = 1.42 + (0.7  * (p_COCoreMass - 2.4));   }
    else if (utils::Compare(p_COCoreMass, 3.6 ) < 0) { remnantMass = 1.32 + (0.25 * (p_COCoreMass - 3.2));   }
    else if (utils::Compare(p_COCoreMass, 4.05) < 0) { remnantMass = p_Mass * NEUTRINO_LOSS_FALLBACK_FACTOR; }      // Going to be a Black Hole
    else if (utils::Compare(p_COCoreMass, 4.6 ) < 0) { remnantMass = 1.5;                                    }
    else if (utils::Compare(p_COCoreMass, 5.7 ) < 0) { remnantMass = p_Mass * NEUTRINO_LOSS_FALLBACK_FACTOR; }      // Going to be a Black Hole
    else if (utils::Compare(p_COCoreMass, 6.0 ) < 0) { remnantMass = 1.64 - (0.2  * (p_COCoreMass - 5.7));   }
    else                                             { remnantMass = p_Mass * NEUTRINO_LOSS_FALLBACK_FACTOR; }      // Going to be a Black Hole

    return remnantMass;
}


/*
 * Calculate the baryonic mass of the remnant
 *
 * Fryer et al. 2012, eq 12
 *
 *
 * double CalculateBaryonicRemnantMass(const double p_ProtoMass, double p_FallbackMass)
 *
 * @param   [IN]    p_ProtoMass                 Mass of proto compact object in Msol
 * @param   [IN]    p_FallbackMass              Mass falling back onto proto compact object Mfb = fb*(MpreSN - Mproto)
 * @return                                      Baryonic mass of the remnant in Msol
 */
double GiantBranch::CalculateBaryonicRemnantMass(const double p_ProtoMass, double p_FallbackMass) {
    return p_ProtoMass + p_FallbackMass; //Mfb = fb*(MpreSN - Mproto);
}


/*
 * Calculate the gravitational mass of the remnant
 *
 * Fryer et al. 2012, eqs 13 & 14
 *
 *
 * double CalculateGravitationalRemnantMass(const double p_BaryonicRemnantMass)
 *
 * @param   [IN]    p_BaryonicRemnantMass       Baryonic remnant mass in Msol
 * @return                                      Gravitational mass of the remnant in Msol
 */
double GiantBranch::CalculateGravitationalRemnantMass(const double p_BaryonicRemnantMass) {

    ERROR  error = ERROR::NONE;

    double root;

    // decide whether to calculate GravitationalRemnantMass from Fryer+2012, Eq.13 for Neutron Star or Black Hole 
    // then calculate GravitationalRemnantMass 
    
    if (utils::Compare(p_BaryonicRemnantMass, m_BaryonicMassOfMaximumNeutronStarMass) < 0) {
        std::tie(error, root) = utils::SolveQuadratic(0.075, 1.0, -p_BaryonicRemnantMass);                 // Neutron Star
        if (error == ERROR::NO_REAL_ROOTS)
            THROW_ERROR(error);
    } 
    else {                                                                                                 // Black hole
        root = BH::CalculateNeutrinoMassLoss_Static(p_BaryonicRemnantMass);                                // Convert to gravitational mass due to neutrino mass loss
    }

    return root;
}


/*
 * Calculate the mass falling back onto the proto compact object
 *
 * Fryer et al. 2012, eq 11
 *
 *
 * double CalculateFallbackMass(const double p_PreSNMass, const double p_ProtoMass, const double p_Fallback)
 *
 * @param   [IN]    p_PreSNMass                 Pre supernova stellar mass in Msol
 * @param   [IN]    p_ProtoMass                 Pre supernova Carbon Oxygen (CO) core mass in Msol
 * @param   [IN]    p_Fallback                  Fraction of mass falling back onto proto object
 * @return                                      Mass falling back onto proto object
 */
double GiantBranch::CalculateFallbackMass(const double p_PreSNMass, const double p_ProtoMass, const double p_Fallback) {
    return p_Fallback * (p_PreSNMass - p_ProtoMass);
}


/*
 * Calculate the fraction of mass falling back onto the proto compact object using the rapid prescription
 *
 * Fryer et al. 2012, eq 15 & 16
 *
 *
 * double CalculateFallbackFractionRapid(const double p_PreSNMass, const double p_ProtoMass, const double p_COCoreMass)
 *
 * @param   [IN]    p_PreSNMass                 Pre supernova stellar mass in Msol
 * @param   [IN]    p_ProtoMass                 Fe/Ni proto compact object mass in Msol
 * @param   [IN]    p_COCoreMass                Pre supernova Carbon Oxygen (CO) core mass in Msol
 * @return                                      Fraction of mass falling back onto proto object
 */
double GiantBranch::CalculateFallbackFractionRapid(const double p_PreSNMass, const double p_ProtoMass, const double p_COCoreMass) {

    double fb;

    double spread = p_PreSNMass - p_ProtoMass;

    if (utils::Compare(p_COCoreMass, 2.5) < 0) {
        fb = 0.2 / spread;
    }
    else if (utils::Compare(p_COCoreMass, 6.0) < 0) {
        fb = ((0.286 * p_COCoreMass) - 0.514) / spread;
    }
    else if (utils::Compare(p_COCoreMass, 7.0) < 0) {
        fb = 1.0;
    }
    else if (utils::Compare(p_COCoreMass, 11.0) < 0) {
        double a1 = 0.25 - (1.275 / spread);
        double b1 = (-11.0 * a1) + 1.0;
        fb        = (a1 * p_COCoreMass) + b1;
    }
    else {
        fb = 1.0;
    }

    return std::max(0.0, std::min(1.0, fb));
}


/*
 * Calculate the mass of the proto core
 *
 * Fryer et al. 2012, eq 18
 *
 *
 * double CalculateProtoCoreMassDelayed(const double p_COCoreMass)
 *
 * @param   [IN]    p_COCoreMass                Pre supernova Carbon Oxygen (CO) core mass in Msol
 * @return                                      Mass of the Fe/Ni proto core in Msol
 */
double GiantBranch::CalculateProtoCoreMassDelayed(const double p_COCoreMass) {

    double protoMass;

         if (utils::Compare(p_COCoreMass,  3.5) < 0) { protoMass = 1.2; }
    else if (utils::Compare(p_COCoreMass,  6.0) < 0) { protoMass = 1.3; }
    else if (utils::Compare(p_COCoreMass, 11.0) < 0) { protoMass = 1.4; }
    else                                             { protoMass = 1.6; }

    return protoMass;
}


/*
 * Calculate the fraction of mass falling back onto the proto compact object using the delayed prescription
 *
 * Fryer et al. 2012, eq 19
 *
 *
 * double CalculateFallbackFractionDelayed(const double p_PreSNMass, const double p_ProtoMass, const double p_COCoreMass)
 *
 * @param   [IN]    p_PreSNMass                 Pre supernova stellar mass in Msol
 * @param   [IN]    p_ProtoMass                 Fe/Ni proto compact object mass in Msol
 * @param   [IN]    p_COCoreMass                Pre supernova Carbon Oxygen (CO) core mass in Msol
 * @return                                      Fraction of mass falling back onto proto object
 */
double GiantBranch::CalculateFallbackFractionDelayed(const double p_PreSNMass, const double p_ProtoMass, const double p_COCoreMass) {

    double fb;

    double spread = p_PreSNMass - p_ProtoMass;

    if (utils::Compare(p_COCoreMass, 2.5) < 0) {
        fb = 0.2 / spread;
    }
    else if (utils::Compare(p_COCoreMass, 3.5) < 0) {
        fb = ((0.5 * p_COCoreMass) - 1.05) / spread;
    }
    else if (utils::Compare(p_COCoreMass, 11.0) < 0) {
        double a2 = 0.133 - (0.093 / spread);
        double b2 = (-11.0 * a2) + 1.0;
        fb        = (a2 * p_COCoreMass) + b2;
    }
    else {
        fb = 1.0;
    }

    return std::max(0.0, std::min(1.0, fb));
}


/*
 * Calculate the remnant mass using one of the Fryer SN Engines
 *
 * Ref?  Fryer et al. 2012?
 *
 *
 * std::tuple<double, double> CalculateRemnantMassByFryer2012(const double p_Mass, const double p_COCoreMass)
 *
 * @param   [IN]    p_Mass                      Pre supernova mass in Msol
 * @param   [IN]    p_COCoreMass                Pre supernova Carbon Oxygen (CO) core mass in Msol
 * @return                                      Tuple containing Remnant mass in Msol and updated fraction of mass falling back onto compact object
 */
std::tuple<double, double> GiantBranch::CalculateRemnantMassByFryer2012(const double p_Mass, const double p_COCoreMass) {

    double mProto;
    double fallbackMass;
    double baryonicRemnantMass;

    double fallbackFraction         = 0.0;
    double gravitationalRemnantMass = 0.0;

    switch (OPTIONS->FryerSupernovaEngine()) {                                                  // which SN_ENGINE?

        case SN_ENGINE::DELAYED:

            mProto           = CalculateProtoCoreMassDelayed(p_COCoreMass);
            fallbackFraction = CalculateFallbackFractionDelayed(p_Mass, mProto, p_COCoreMass);
            fallbackMass     = CalculateFallbackMass(p_Mass, mProto, fallbackFraction);

            baryonicRemnantMass      = CalculateBaryonicRemnantMass(mProto, fallbackMass);
            gravitationalRemnantMass = CalculateGravitationalRemnantMass(baryonicRemnantMass);
            break;

        case SN_ENGINE::RAPID:

            mProto           = FRYER_PROTO_CORE_MASS_RAPID;
            fallbackFraction = CalculateFallbackFractionRapid(p_Mass, mProto, p_COCoreMass);
            fallbackMass     = CalculateFallbackMass(p_Mass, mProto, fallbackFraction);

            baryonicRemnantMass      = CalculateBaryonicRemnantMass(mProto, fallbackMass);
            gravitationalRemnantMass = CalculateGravitationalRemnantMass(baryonicRemnantMass);
            break;
    
        default:                                                                                // unknown prescription
            // the only way this can happen is if someone added an SN_ENGINE
            // and it isn't accounted for in this code.  We should not default here, with or without a warning.
            // We are here because the user chose a prescription this code doesn't account for, and that should
            // be flagged as an error and result in termination of the evolution of the star or binary.
            // The correct fix for this is to add code for the missing prescription or, if the missing
            // prescription is superfluous, remove it from the option.

            THROW_ERROR(ERROR::UNKNOWN_SN_ENGINE);                                              // throw error
    }

    return std::make_tuple(gravitationalRemnantMass, fallbackFraction);
}


/*
 * Calculate the remnant mass using the new Fryer prescription from 2022
 *
 * Fryer et al. 2022 eq. 5
 *
 *
 * std::tuple<double, double> CalculateRemnantMassByFryer2022(const double p_Mass, const double p_COCoreMass)
 *
 * @param   [IN]    p_Mass                      Pre supernova mass in Msol
 * @param   [IN]    p_COCoreMass                Pre supernova Carbon Oxygen (CO) core mass in Msol
 * @return                                      Tuple containing Remnant mass in Msol and updated fraction of mass falling back onto compact object
 */
std::tuple<double, double> GiantBranch::CalculateRemnantMassByFryer2022(const double p_Mass, const double p_COCoreMass) {


    double mProto;
    double fallbackMass;
    double baryonicRemnantMass;

    double fallbackFraction         = 0.0;
    double gravitationalRemnantMass = 0.0;

    double fmix         = OPTIONS->Fryer22fmix();
    double fact         = p_COCoreMass / fmix;
    baryonicRemnantMass = 1.2 + 0.05 * fmix + 0.01 * fact * fact + exp(fmix * (p_COCoreMass - OPTIONS->Fryer22Mcrit())) ; // equation 5. 
    baryonicRemnantMass = std::min(baryonicRemnantMass, p_Mass);                                // clamp baryonicRemnantMass to total mass

    // Now the proto mass, which is only used for the calculation of kicks, will
    // still be calculated using the DELAYED/RAPID prescriptions from Fryer 2012
    switch (OPTIONS->FryerSupernovaEngine()) {                                                  // which SN_ENGINE?

        case SN_ENGINE::DELAYED:  
            mProto                   = CalculateProtoCoreMassDelayed(p_COCoreMass);

            fallbackMass             = std::max(0.0, baryonicRemnantMass - mProto);             // fallbackMass larger than 0
            fallbackFraction         = fallbackMass / (p_Mass - mProto);
            fallbackFraction         = std::max(0.0, std::min(1.0, fallbackFraction));          // make sure the fb fraction lies between 0-1
            gravitationalRemnantMass = CalculateGravitationalRemnantMass(baryonicRemnantMass);
            break;

        case SN_ENGINE::RAPID:  
            mProto                   = FRYER_PROTO_CORE_MASS_RAPID;

            fallbackMass             = std::max(0.0, baryonicRemnantMass - mProto);             // fallbackMass larger than 0
            fallbackFraction         = fallbackMass / (p_Mass - mProto);
            fallbackFraction         = std::max(0.0, std::min(1.0, fallbackFraction));          // make sure the fb fraction lies between 0-1
            gravitationalRemnantMass = CalculateGravitationalRemnantMass(baryonicRemnantMass);
            break;
    
        default:                                                                                // unknown prescription
            // the only way this can happen is if someone added an SN_ENGINE
            // and it isn't accounted for in this code.  We should not default here, with or without a warning.
            // We are here because the user chose a prescription this code doesn't account for, and that should
            // be flagged as an error and result in termination of the evolution of the star or binary.
            // The correct fix for this is to add code for the missing prescription or, if the missing
            // prescription is superfluous, remove it from the option.

            THROW_ERROR(ERROR::UNKNOWN_SN_ENGINE);                                              // throw error
    }
                                   
    return std::make_tuple(gravitationalRemnantMass, fallbackFraction);
}



/*
 * Calculate fallback using linear interpolation
 *
 * Belczynski et al. 2002
 * Description?
 *
 *
 * double CalculateFallbackByBelczynski2002(const double p_COCoreMass)
 *
 * @param   [IN]    p_COCoreMass                Pre supernova Carbon Oxygen (CO) core mass in Msol
 * @return                                      Fraction of mass falling back onto compact object
 */
double GiantBranch::CalculateFallbackByBelczynski2002(const double p_COCoreMass) {
    return utils::Compare(p_COCoreMass, 5.0) <= 0 ? 0.0 : ((utils::Compare(p_COCoreMass, 7.6) < 0) ? (p_COCoreMass - 5.0) / 2.6 : 1.0);
}


/*
 * Calculate remnant mass
 *
 * Formula used by Hurley code not given in Hurley et al 2000 but in Belczynski et al. 2002
 *
 *
 * double CalculateRemnantMassByBelczynski2002(const double p_Mass, const double p_COCoreMass, const double p_FallbackFraction)
 *
 * @param   [IN]    p_COCoreMass                Pre supernova Carbon Oxygen (CO) core mass in Msol
 * @param   [IN]    p_FallbackFraction          Fraction of mass falling back onto compact object
 * @return                                      Remnant mass in Msol
 */
double GiantBranch::CalculateRemnantMassByBelczynski2002(const double p_Mass, const double p_COCoreMass, const double p_FallbackFraction) {
    double McFeNi = (utils::Compare(p_COCoreMass, 2.5) < 0) ? (0.161767 * p_COCoreMass) + 1.067055 : (0.314154 * p_COCoreMass) + 0.686088; // Iron core mass
    return McFeNi + (p_FallbackFraction * (p_Mass - McFeNi));
}


/*
 * Driver function for Core Collapse Supernovas
 *
 * This function determines which prescription is used for the core collapse SN (via program options)
 *
 * The function calls prescription functions that update the following parameters:
 *      Mass, stellarType, drawnKickMagnitude, kickMagnitude
 *
 * At the end of this function we set the following parameters which are (so far) independent of the
 * ccSN prescriptions (but do depend on the parameters above):
 *      Luminosity, Radius, Temperature, supernova events: current = SN, past = CCSN
 *
 *
 * STELLAR_TYPE ResolveCoreCollapseSN()
 *
 * @return                                      The stellar type to which the star should evolve
 */
STELLAR_TYPE GiantBranch::ResolveCoreCollapseSN() {

    STELLAR_TYPE stellarType = m_StellarType;
    double mass              = m_Mass;                                                                      // initial mass

    switch (OPTIONS->RemnantMassPrescription()) {                                                           // which prescription?

        case REMNANT_MASS_PRESCRIPTION::HURLEY2000:                                                         // Hurley 2000

            m_SupernovaDetails.fallbackFraction = 0.0;                                                      // not defined
            m_Mass                              = Remnants::CalculateRemnantMass_Static(m_COCoreMass);
            break;

        case REMNANT_MASS_PRESCRIPTION::BELCZYNSKI2002:                                                     // Belczynski 2002

            m_SupernovaDetails.fallbackFraction = CalculateFallbackByBelczynski2002(m_CoreMass);
            m_Mass                              = CalculateRemnantMassByBelczynski2002(m_Mass, m_CoreMass, m_SupernovaDetails.fallbackFraction);
            break;

        case REMNANT_MASS_PRESCRIPTION::FRYER2012:                                                          // Fryer 2012

            std::tie(m_Mass, m_SupernovaDetails.fallbackFraction) = CalculateRemnantMassByFryer2012(m_Mass, m_COCoreMass);
            break;


        case REMNANT_MASS_PRESCRIPTION::FRYER2022:                                                          // Fryer 2022

            std::tie(m_Mass, m_SupernovaDetails.fallbackFraction) = CalculateRemnantMassByFryer2022(m_Mass, m_COCoreMass);
            break;

        case REMNANT_MASS_PRESCRIPTION::MULLER2016:                                                         // Muller 2016

            m_SupernovaDetails.fallbackFraction = 0.0;                                                      // no subsequent kick adjustment by fallback fraction needed
            m_Mass                              = CalculateRemnantMassByMuller2016(m_Mass, m_COCoreMass);
            break;

        case REMNANT_MASS_PRESCRIPTION::MULLERMANDEL:                                                       // Mandel & Mueller, 2020

            m_SupernovaDetails.fallbackFraction = 0.0;                                                      // no subsequent kick adjustment by fallback fraction needed
            m_Mass                              = CalculateRemnantMassByMullerMandel(m_COCoreMass, m_HeCoreMass);
            break;

        case REMNANT_MASS_PRESCRIPTION::SCHNEIDER2020:                                                      // Schneider 2020

            m_Mass                              = CalculateRemnantMassBySchneider2020(m_COCoreMass);
            m_SupernovaDetails.fallbackFraction = utils::Compare(m_Mass, OPTIONS->MaximumNeutronStarMass() ) > 0 ? (m_Mass - NEUTRON_STAR_MASS) / (mass - NEUTRON_STAR_MASS) : 0.0;                                                                   // Fallback fraction of mass beyond proto-neutron-star for BH formation and kicks
            break;

        case REMNANT_MASS_PRESCRIPTION::SCHNEIDER2020ALT:                                                   // Schneider 2020, alternative

            m_Mass                              = CalculateRemnantMassBySchneider2020Alt(m_COCoreMass);
            m_SupernovaDetails.fallbackFraction = utils::Compare(m_Mass, OPTIONS->MaximumNeutronStarMass() ) > 0 ? (m_Mass - NEUTRON_STAR_MASS) / (mass - NEUTRON_STAR_MASS) : 0.0;                                                                   // Fallback fraction of mass beyond proto-neutron-star for BH formation and kicks
            break;
        
        case REMNANT_MASS_PRESCRIPTION::MALTSEV2024:                                                        // Maltsev+ 2024

            m_SupernovaDetails.fallbackFraction = 0.0;                                                      // no subsequent kick adjustment by fallback fraction needed; MULLERMANDEL kick prescription should be used
            m_Mass                              = CalculateRemnantMassByMaltsev2024(m_COCoreMass, m_HeCoreMass);
            break;
            
    
        default:                                                                                            // unknown prescription
            // the only way this can happen is if someone added a REMNANT_MASS_PRESCRIPTION
            // and it isn't accounted for in this code.  We should not default here, with or without a warning.
            // We are here because the user chose a prescription this code doesn't account for, and that should
            // be flagged as an error and result in termination of the evolution of the star or binary.
            // The correct fix for this is to add code for the missing prescription or, if the missing
            // prescription is superfluous, remove it from the option.

            THROW_ERROR(ERROR::UNKNOWN_REMNANT_MASS_PRESCRIPTION);                                          // throw error
    }
    
    // Set the stellar type to which the star should evolve (either use prescription or MAXIMUM_NS_MSS)
    if (OPTIONS->RemnantMassPrescription() == REMNANT_MASS_PRESCRIPTION::MULLER2016) {
        stellarType = CalculateRemnantTypeByMuller2016(m_COCoreMass);
    }
    else if (OPTIONS->RemnantMassPrescription() == REMNANT_MASS_PRESCRIPTION::MULLERMANDEL) {
        if (utils::Compare(m_Mass, OPTIONS->MaximumNeutronStarMass() ) > 0)
            stellarType = STELLAR_TYPE::BLACK_HOLE;
        else
            stellarType = STELLAR_TYPE::NEUTRON_STAR;
    }
    else if (OPTIONS->RemnantMassPrescription() == REMNANT_MASS_PRESCRIPTION::HURLEY2000) {
        stellarType = (utils::Compare(m_Mass, 1.8 ) > 0)
                        ? STELLAR_TYPE::BLACK_HOLE
                        : STELLAR_TYPE::NEUTRON_STAR;                                                       // Hurley+ 2000, Eq. (92)
    }
    else if (utils::Compare(m_Mass, OPTIONS->MaximumNeutronStarMass()) > 0) {
        std::tie(m_Luminosity, m_Radius, m_Temperature) = BH::CalculateCoreCollapseSNParams_Static(m_Mass);
        stellarType = STELLAR_TYPE::BLACK_HOLE;
    }
    else {
        std::tie(m_Luminosity, m_Radius, m_Temperature) = NS::CalculateCoreCollapseSNParams_Static(m_Mass);
        stellarType = STELLAR_TYPE::NEUTRON_STAR;
    }

    if (utils::Compare(mass,m_CoreMass) == 0 && utils::Compare(m_HeCoreMass, m_COCoreMass) == 0) {          // entire star is CO core, so this is a USSN
        SetSNCurrentEvent(SN_EVENT::USSN);                                                                  // flag ultra-stripped SN happening now
        SetSNPastEvent(SN_EVENT::USSN);                                                                     // ... and will be a past event
    }

    SetSNCurrentEvent(SN_EVENT::CCSN);                                                                      // flag core-collapse SN happening now
    SetSNPastEvent(SN_EVENT::CCSN);                                                                         // ... and will be a past event

    return stellarType;
}


/*
 * Resolve Electron Capture Supernova
 *
 * Calculate the mass of the remnant and set remnant type - always a Neutron Star
 * Updates attributes of star; sets SN flags
 *
 *
 * Short hand wavy story. The core ignites ONeMg and collapses through electron
 * capture (e.g., Nomoto 1984 for thorough discussion and a summary in Nomoto 1987).
 * The explosion is likely accompanied by a low natal kick
 *
 *
 * STELLAR_TYPE ResolveElectronCaptureSN()
 *
 * @return                                      Stellar type of remnant (always STELLAR_TYPE::NEUTRON_STAR)
 */
STELLAR_TYPE GiantBranch::ResolveElectronCaptureSN() {
    m_Mass       = MECS_REM;                                                                // defined in constants.h
    m_CoreMass   = m_Mass;
    m_HeCoreMass = m_Mass;
    m_COCoreMass = m_Mass;
    m_Mass0      = m_Mass;
        
    SetSNCurrentEvent(SN_EVENT::ECSN);                                                      // electron capture SN happening now
    SetSNPastEvent(SN_EVENT::ECSN);                                                         // ... and will be a past event

    return STELLAR_TYPE::NEUTRON_STAR;
}


/*
 * Resolve Pair-Instability Supernova
 *
 * Calculate the mass of the remnant and set remnant type according to mass
 * Updates attributes of star; sets SN events
 *
 *
 * Short handwavy story.  The core is hot and massive enough that there is a significant
 * amount of high-energetic gamma rays which can create an electron-positron pair.
 * When a significant amounts of pairs are created the photons taken away
 * reduce the radiative pressure. The core contracts becomes hotter creating
 * more pairs. This runaway process will end in a SN that explodes the entire core without
 * leaving a remnant.
 *
 *
 * STELLAR_TYPE ResolvePairInstabilitySN()
 *
 * @return                                      Stellar type of remnant
 */
STELLAR_TYPE GiantBranch::ResolvePairInstabilitySN() {

    m_Luminosity  = 0.0;
    m_Radius      = 0.0;
    m_Temperature = 0.0;
    m_Mass        = 0.0;
    m_CoreMass    = 0.0;
    m_COCoreMass  = 0.0;
    m_HeCoreMass  = 0.0;

    m_SupernovaDetails.drawnKickMagnitude = 0.0;
    m_SupernovaDetails.kickMagnitude      = 0.0;
    m_SupernovaDetails.fallbackFraction   = 0.0;

    SetSNCurrentEvent(SN_EVENT::PISN);                                                                  // pair instability SN happening now
    SetSNPastEvent(SN_EVENT::PISN);                                                                     // ... and will be a past event

    return STELLAR_TYPE::MASSLESS_REMNANT;
}


/*
 * Resolve Pulsational Pair-Instability Supernova
 *
 * Calculate the mass of the remnant and set remnant type according to mass
 * Updates attributes of star; sets SN events
 *
 *
 * STELLAR_TYPE ResolvePulsationalPairInstabilitySN()
 *
 * @return                                      Stellar type of remnant
 */
STELLAR_TYPE GiantBranch::ResolvePulsationalPairInstabilitySN() {

    STELLAR_TYPE stellarType = m_StellarType;

    double baryonicMass;
    switch (OPTIONS->PulsationalPairInstabilityPrescription()) {                                        // which prescription?

        case PPI_PRESCRIPTION::COMPAS:                                                                  // deprecated Feb 2025, to be removed
        case PPI_PRESCRIPTION::WOOSLEY:                                                                 // Woosley 2017 https://arxiv.org/abs/1608.08939
            baryonicMass = m_HeCoreMass;                                                                // strip off the hydrogen envelope if any was left
            m_Mass       = BH::CalculateNeutrinoMassLoss_Static(baryonicMass);                          // convert to gravitational mass due to neutrino mass loss
            break;

        case PPI_PRESCRIPTION::STARTRACK:                                                               // Belczynski et al. 2016 https://arxiv.org/abs/1607.03116
            baryonicMass = std::min(m_HeCoreMass, STARTRACK_PPISN_HE_CORE_MASS);                        // strip off the hydrogen envelope if any was left, limit to STARTRACK_PPISN_HE_CORE_MASS (default 45 Msun)
            m_Mass       = BH::CalculateNeutrinoMassLoss_Static(baryonicMass);                          // convert to gravitational mass due to neutrino mass loss
            break;

        case PPI_PRESCRIPTION::MARCHANT: {                                                              // Marchant et al. 2018 https://arxiv.org/abs/1810.13412

            // pow() is slow - use multiplication
            double HeCoreMass_2 = m_HeCoreMass * m_HeCoreMass;
            double HeCoreMass_3 = HeCoreMass_2 * m_HeCoreMass;
            double HeCoreMass_4 = HeCoreMass_2 * HeCoreMass_2;
            double HeCoreMass_5 = HeCoreMass_3 * HeCoreMass_2;
            double HeCoreMass_6 = HeCoreMass_3 * HeCoreMass_3;
            double HeCoreMass_7 = HeCoreMass_6 * m_HeCoreMass;

            double ratioOfRemnantToHeCoreMass = std::max(0.0, std::min(1.0, (-1.63057326E-08 * HeCoreMass_7) +
                                                                            ( 5.36316755E-06 * HeCoreMass_6) +
                                                                            (-7.52206933E-04 * HeCoreMass_5) +
                                                                            ( 5.83107626E-02 * HeCoreMass_4) +
                                                                            (-2.69801221E+00 * HeCoreMass_3) +
                                                                            ( 7.45060098E+01 * HeCoreMass_2) +
                                                                            (-1.13694590E+03 * m_HeCoreMass) +
                                                                              7.39643451E+03));

            baryonicMass = ratioOfRemnantToHeCoreMass * m_HeCoreMass;                                   // strip off the hydrogen envelope if any was left
            m_Mass       = BH::CalculateNeutrinoMassLoss_Static(baryonicMass);                          // convert to gravitational mass due to neutrino mass loss
            } break;

        case PPI_PRESCRIPTION::FARMER: {                                                                // Farmer et al. 2019 http://dx.doi.org/10.3847/1538-4357/ab518b
            double totalMassPrePPISN = m_Mass;                                                          // save the total stellar mass 
                                                                                                        // three cases:
            if (utils::Compare(m_COCoreMass, FARMER_PPISN_UPP_LIM_LIN_REGIME) < 0) {
                m_Mass = m_COCoreMass + 4.0;                                                            // a linear relation below CO core masses of 38 Msun
            }
            else if (utils::Compare(m_COCoreMass, FARMER_PPISN_UPP_LIM_QUAD_REGIME) < 0) {              // a quadratic relation in CO core mass for 38 =< CO_core < 60  JR: shouldn't have constants here (in the comment) - they may change
                const double a1 = -0.096;
                const double a2 = 8.564;
                const double a3 = -2.07;
                const double a4 = -152.97;

                m_Mass = a1 * PPOW(m_COCoreMass, 2.0) + a2 * m_COCoreMass + a3 * m_Log10Metallicity + a4;
            }
            else if (utils::Compare(m_COCoreMass, FARMER_PPISN_UPP_LIM_INSTABILLITY) < 0) {             // no remnant between 60 - 140 Msun  JR: shouldn't have constants here (in the comment) - they may change
                m_Mass = 0.0;
            }
            else {                                                                                      // BH mass becomes CO-core mass above the PISN gap
                m_Mass = m_COCoreMass;
            }

            m_Mass = std::min(totalMassPrePPISN, m_Mass);                                               // check if remnant mass is bigger than total mass    
            } break;
    
        case PPI_PRESCRIPTION::HENDRIKS: {    
            // Prescription from Hendriks et al. 2023 (https://arxiv.org/abs/2309.09339)
            // Based on Renzo et al. 2022 (https://iopscience.iop.org/article/10.3847/2515-5172/ac503e)
            // 
            // Suggest using --PPI-upper-limit 80.0 and --PISN-lower-limit 80.0

            double DeltaMPPICOShift = OPTIONS->PulsationalPairInstabilityCOCoreShiftHendriks();
            double DeltaMPPIExtraML = 0.0; 								// Make an option? Currently does nothing

            // Equation (6) of Hendricks et al. 2023			
            double PPIOnset        = m_COCoreMass - DeltaMPPICOShift - 34.8;
            double PPIOnsetSquared = PPIOnset * PPIOnset;
            double PPIOnsetCubed   = PPIOnsetSquared * PPIOnset;
            double firstTerm  = (0.0006 * m_Log10Metallicity + 0.0054) * PPIOnsetCubed;
            double secondTerm = 0.0013 * PPIOnsetSquared;
            double DeltaMPPI  = firstTerm - secondTerm + DeltaMPPIExtraML;
            
            DeltaMPPI = std::max(DeltaMPPI, 0.0);						// DeltaMPPI, the amount of the He core that's lost in pulsations, is non-negative
            m_Mass = std::max(m_HeCoreMass - DeltaMPPI, 0.0);			// Remnant mass should be non-negative		
            m_Mass = m_Mass > 10.0 ? m_Mass : 0.0;                      // If the predicted remnant mass is below 10 Msun, set it equal to 0 (assume a PISN)

        } break;

        default:                                                                                        // unknown prescription
            // the only way this can happen is if someone added a REMNANT_MASS_PRESCRIPTION
            // and it isn't accounted for in this code.  We should not default here, with or without a warning.
            // We are here because the user chose a prescription this code doesn't account for, and that should
            // be flagged as an error and result in termination of the evolution of the star or binary.
            // The correct fix for this is to add code for the missing prescription or, if the missing
            // prescription is superfluous, remove it from the option.

            THROW_ERROR(ERROR::UNKNOWN_PPI_PRESCRIPTION);                                               // throw error
    }

    if (utils::Compare(m_Mass, 0.0) <= 0) {                                                             // remnant mass <= 0?
        stellarType = ResolvePairInstabilitySN();                                                       // yes - PISN rather than PPISN
    }
    else {                                                                                              // no - PPISN
        SetSNCurrentEvent(SN_EVENT::PPISN);                                                             // pulsational pair instability SN happening now
        SetSNPastEvent(SN_EVENT::PPISN);                                                                // ... and will be a past event

        stellarType   = STELLAR_TYPE::BLACK_HOLE;                                                       // -> black hole
        
        m_Luminosity  = BH::CalculateLuminosityOnPhase_Static();                                        // black hole luminosity
        m_Radius      = BH::CalculateRadiusOnPhase_Static(m_Mass);                                      // Schwarzschild radius (not correct for rotating BH)
        m_Temperature = CalculateTemperatureOnPhase(m_Luminosity, m_Radius);

        m_SupernovaDetails.fallbackFraction = 1.0;                                                      // fraction of mass that falls back
    }

    return stellarType;
}


/*
 * The main supernova function
 *
 * This function determines the type of the supernova and calls the appropriate functions
 * to calculate attributes correctly, and to determine the type of remnant to which the
 * star should evolve.
 *
 *
 * STELLAR_TYPE ResolveSupernova()
 *
 * @return                                      Stellar type of remnant
 */
STELLAR_TYPE GiantBranch::ResolveSupernova() {

    STELLAR_TYPE stellarType = m_StellarType;

    if (IsSupernova()) {                                                                            // is going supernova
                                                                                                    // yes - resolve new supernova event
        // squirrel away some attributes before they get changed...
        m_SupernovaDetails.totalMassAtCOFormation  = m_Mass;
        m_SupernovaDetails.HeCoreMassAtCOFormation = m_HeCoreMass;
        m_SupernovaDetails.COCoreMassAtCOFormation = m_COCoreMass;
        m_SupernovaDetails.coreMassAtCOFormation   = m_CoreMass;

        double snMass = CalculateInitialSupernovaMass();                                            // calculate SN initial mass
        
        SetSNHydrogenContent();                                                                     // set H-rich or H-poor  JR: why don't we do this when we initialise the star at change of stellar type?

        if (                             OPTIONS->UsePulsationalPairInstability()              &&
            utils::Compare(m_HeCoreMass, OPTIONS->PulsationalPairInstabilityLowerLimit()) >= 0 &&
            utils::Compare(m_HeCoreMass, OPTIONS->PulsationalPairInstabilityUpperLimit()) <= 0) {   // Pulsational Pair Instability Supernova

            stellarType = ResolvePulsationalPairInstabilitySN();                                    // BH or MR
        }
        else if (                        OPTIONS->UsePairInstabilitySupernovae()    &&
            utils::Compare(m_HeCoreMass, OPTIONS->PairInstabilityLowerLimit()) >= 0 &&
            utils::Compare(m_HeCoreMass, OPTIONS->PairInstabilityUpperLimit()) <= 0) {              // Pair Instability Supernova

            stellarType = ResolvePairInstabilitySN();                                               // MR
        }
        else if (utils::Compare(snMass, MCBUR2) < 0 && (!m_MassTransferDonorHistory.empty() || OPTIONS->AllowNonStrippedECSN())) {
            stellarType = ResolveElectronCaptureSN();                                               // electron capture SN; requires progenitor to have been a MT donor unless non-stripped ECSN are allowed; forms NS
        }
        else {                                                                                      // Core Collapse Supernova
            stellarType = ResolveCoreCollapseSN();                                                  // BH or NS
        }
        
        // check if the SN actually happened
        if (utils::IsOneOf(stellarType, { STELLAR_TYPE::NEUTRON_STAR, STELLAR_TYPE::BLACK_HOLE, STELLAR_TYPE::MASSLESS_REMNANT })) {
                                                                                                    // SN happened
            if (utils::SNEventType(m_SupernovaDetails.events.current) != SN_EVENT::PISN && !utils::IsOneOf(stellarType, { STELLAR_TYPE::MASSLESS_REMNANT }))
                CalculateSNKickMagnitude(m_Mass, m_SupernovaDetails.totalMassAtCOFormation - m_Mass, stellarType);
        
            if (!utils::IsOneOf(stellarType, { STELLAR_TYPE::NEUTRON_STAR }))
                m_SupernovaDetails.rocketKickMagnitude = 0;                                         // only NSs can get rocket kicks

            // Stash SN details for later printing to the SSE Supernova log.
            // Only if SSE (BSE does its own SN printing), and only if not an ephemeral clone
            // Can't print it now because we may revert state (in Star::EvolveOneTimestep()).
            // Will be printed in Star::EvolveOneTimestep() after timestep is accepted (i.e. we don't revert state).
            // Need to record the stellar type to which the star will switch if we don't revert state.

            if (OPTIONS->EvolutionMode() == EVOLUTION_MODE::SSE && m_ObjectPersistence == OBJECT_PERSISTENCE::PERMANENT)
                StashSupernovaDetails(stellarType);
       }
    }

    return stellarType;
}
