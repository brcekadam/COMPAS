#ifndef __constants_h__
#define __constants_h__

#include <cmath>
#include <string>
#include <vector>
#include <map>
#include <unordered_map>
#include <fstream>
#include <csignal>
#include <limits>

// the defaults size of the boost list that handles variant types is 20 - so only 20 variant types are allowed
// we've exceeded that number - we're at 21 currently - so the size of the boost list needs to be increased
// we have to set the size of the list before we include the boost headers - otherwise boost redefines it
// boost only recognises values of 20, 30, 40, & 50: for now we set it to 30. 
#define BOOST_MPL_CFG_NO_PREPROCESSED_HEADERS
  
#if !defined(BOOST_MPL_LIMIT_LIST_SIZE)
    #if defined(BOOST_MPL_LIST_HPP_INCLUDED)
      # error "BOOST_MPL_LIMIT_LIST_SIZE must be set to accommodate the size of COMPAS_VARIABLE before the including the Boost headers"
    #endif

    #define BOOST_MPL_LIMIT_LIST_SIZE 30

#elif BOOST_MPL_LIMIT_LIST_SIZE < 30
    #if defined(BOOST_MPL_LIST_HPP_INCLUDED)
      # error "BOOST_MPL_LIMIT_LIST_SIZE must be set to accommodate the size of COMPAS_VARIABLE before the including the Boost headers"
    #else
      # error "BOOST_MPL_LIMIT_LIST_SIZE value is too low"
    #endif
    
#endif

#include <boost/variant.hpp>


/*
 * trick to allow SWITCH on literal strings
 * constexpr is (or can be) evaluated at compile-time, hence the ability to SWITCH using it.
 *
 * This function returns a hash value for a string - very small chance of collisions, but if a
 * collision happens the compiler will complain (because CASE values will be the same).
 *
 * Adapted from https://rextester.com/discussion/FGBUG88403/Switch-case-with-strings-in-C-11,
 * but variations available in many places with the right Google search terms...
 */

constexpr uint64_t StrHash(char const* p_Str, uint64_t p_Hash = 14695981039346656037ull) {
    return (*p_Str == 0) ? p_Hash : StrHash(p_Str + 1, (p_Hash * 1099511628211ull) ^ static_cast<uint64_t>(*p_Str));
}

constexpr uint64_t _(char const* p_Str) {
    return StrHash(p_Str);
}


// globals
// I don't really like this, but it's the only way (until we have a globals singleton)

typedef unsigned long int OBJECT_ID;
extern OBJECT_ID globalObjectId;    // used to uniquely identify objects - used primarily for error printing


// common type definitions
// easiest way of making them available globally is to put them here
typedef std::vector<std::string>                                        STR_VECTOR;
typedef std::vector<double>                                             DBL_VECTOR;
typedef std::vector<int>                                                INT_VECTOR;
typedef std::vector<bool>                                               BOOL_VECTOR;
typedef std::tuple <double, double>                                     DBL_DBL;
typedef std::tuple <double, double, double>                             DBL_DBL_DBL;
typedef std::tuple <double, double, double, double>                     DBL_DBL_DBL_DBL;
typedef std::tuple<std::string, std::string>                            STR_STR;
typedef std::tuple<std::string, std::string, std::string>               STR_STR_STR;
typedef std::tuple<std::string, std::string, std::string, std::string>  STR_STR_STR_STR;
typedef std::vector<std::tuple<DBL_VECTOR, DBL_VECTOR, DBL_VECTOR, DBL_VECTOR, DBL_VECTOR, DBL_VECTOR, DBL_VECTOR>> GE_QCRIT_RADII_QCRIT_VECTOR;
typedef std::tuple<DBL_VECTOR, GE_QCRIT_RADII_QCRIT_VECTOR> GE_QCRIT_TABLE; 
typedef std::vector<std::tuple<DBL_VECTOR, DBL_VECTOR>> GE_QCRIT_RADII_QCRIT_VECTOR_HE;
typedef std::tuple<DBL_VECTOR, GE_QCRIT_RADII_QCRIT_VECTOR_HE> GE_QCRIT_TABLE_HE; 

#include "typedefs.h"


// Constants in SI
// CPLB: Use CODATA values where applicable http://physics.nist.gov/cuu/Constants/index.html

// cmath (included file) provides the following pi-related constants:
//
// pi	        M_PI	    3.14159265358979323846
// pi/2         M_PI_2	    1.57079632679489661923
// pi/4         M_PI_4	    0.785398163397448309616
// 1/pi         M_1_PI	    0.318309886183790671538
// 2/pi         M_2_PI	    0.636619772367581343076
// 2/sqrt(pi)	M_2_SQRTPI	1.12837916709551257390
//
// I've added _2_PI and SQRT_M_2_PI below

#undef COMPARE_GLOBAL_TOLERANCE // define/undef this to compare floats with/without tolerance (see FLOAT_TOLERANCE_ABSOLUTE, FLOAT_TOLERANCE_RELATIVE and Compare() function)

constexpr double FLOAT_TOLERANCE_ABSOLUTE               = 0.0000005;                                                // absolute tolerance for floating-point comparisons if COMPARE_GLOBAL_TOLERANCE is defined
constexpr double FLOAT_TOLERANCE_RELATIVE               = 0.0000005;                                                // relative tolerance for floating-point comparisons if COMPARE_GLOBAL_TOLERANCE is defined

constexpr double ROOT_ABS_TOLERANCE                     = 1.0E-6;                                                   // absolute tolerance for root finder
constexpr double ROOT_REL_TOLERANCE                     = 1.0E-6;                                                   // relative tolerance for root finder

constexpr std::size_t MAX_STACK_TRACE_SIZE              = 64;                                                       // for debugging - overkill, but just in case

// initialisation constants
constexpr double DEFAULT_INITIAL_DOUBLE_VALUE           = 0.0;                                                      // default initial value for double variables
constexpr double DEFAULT_INITIAL_INTEGER_VALUE          = 0;                                                        // default initial value for int variables
constexpr double DEFAULT_INITIAL_ULONGINT_VALUE         = 0l;                                                       // default initial value for unsigned long int variables
constexpr double DEFAULT_INITIAL_BOOLEAN_VALUE          = false;                                                    // default initial value for bool variables


// conversion constants

// mass
constexpr double G_TO_KG                                = 1.0E-3;                                                   // convert grams to kg
constexpr double MSOL_TO_G                              = 1.98892E33;                                               // convert Solar Mass to g
constexpr double MSOL_TO_KG                             = MSOL_TO_G * G_TO_KG;                                      // convert Solar Mass to kg
constexpr double KG_TO_MSOL                             = 1.0 / MSOL_TO_KG;                                         // convert kg to Solar Mass

// length
constexpr double KM_TO_CM 					            = 1.0E5;									                // convert km to cm
constexpr double KM_TO_M                                = 1000.0;                                                   // convert km to m
constexpr double CM_TO_M                                = 1.0E-2;                                                   // convert cm to m

constexpr double RSOL_TO_KM                             = 6.957E5;                                                  // convert Solar Radius (RSOL) to km
constexpr double RSOL_TO_CM                             = 6.957E10;                                                 // convert Solar Radius (RSOL) to cm
constexpr double RSOL_TO_AU                             = 0.00465047;                                               // convert Solar Radius (RSOL) to AU

constexpr double AU_TO_CM                               = 14959787070000.0;                                         // convert Astronomical Units (AU) to cm
constexpr double AU_TO_RSOL				                = 1.0 / RSOL_TO_AU;                                         // convert Astronomical Units AU to Solar Radius RSOL
constexpr double AU_TO_KM                               = AU_TO_CM / 1.0E5;                                         // convert Astronomical Units AU to km

constexpr double KM_TO_RSOL					            = 1.0 / RSOL_TO_KM;						                    // convert km to Solar Radius (RSOL)
constexpr double KM_TO_AU                               = 1.0 / AU_TO_KM;                                           // convert km to Astronomical Units AU

// time
constexpr double DAYS_IN_QUAD                           = 1461.0;                                                   // number of days in any 4-year period
constexpr double DAYS_IN_YEAR                           = DAYS_IN_QUAD / 4.0;                                       // mean days per year, given DAYS_IN_QUAD
constexpr double SECONDS_IN_DAY                         = 24.0 * 60.0 * 60.0;                                       // number of seconds in 1 day
constexpr double SECONDS_IN_QUAD                        = DAYS_IN_QUAD * SECONDS_IN_DAY;                            // number of seconds in 1 quad
constexpr double SECONDS_IN_YEAR                        = SECONDS_IN_QUAD / 4.0;                                    // number of seconds in 1 year
constexpr double SECONDS_IN_MYR                         = SECONDS_IN_YEAR * 1.0E6;                                  // number of seconds in 1 Myr
constexpr double SECONDS_IN_MS                          = 1.0E-3;                                                   // number of seconds in 1 millisecond

constexpr double MYR_TO_YEAR                            = 1.0E6;                                                    // convert Myr to year
constexpr double YEAR_TO_MYR                            = 1.0E-6;                                                   // convert year to Myr

// energy
constexpr double JOULES_TO_ERG                          = 1.0E7;                                                    // convert Joules to Erg

// B field
constexpr double TESLA_TO_GAUSS                         = 1.0E4;					                                // convert Tesla to Gauss
constexpr double GAUSS_TO_TESLA                         = 1.0 / TESLA_TO_GAUSS;                                     // convert Gauss to Tesla

// systems
constexpr double CGS_SI                                 = G_TO_KG * CM_TO_M * CM_TO_M;                              // convert CGS to SI

// opacity
constexpr double OPACITY_CGS_TO_SI                      = 0.1;                                                      // cm^2 g^-1 to m^2 kg^-1

// constants

constexpr double _2_PI                                  = M_PI * 2.0;                                               // 2PI
constexpr double PI_2                                   = M_PI * M_PI;                                              // PI squared
constexpr double SQRT_M_2_PI                            = 0.79788456080286536;                                      // sqrt(2/PI)
constexpr double DEGREE                                 = M_PI / 180.0;                                             // 1 degree in radians

constexpr double GAMMA_E                                = 0.57721566490153286;                                      // Euler's Constant

constexpr double H0                                     = 67.8;                                                     // Hubble's Constant in km s^-1 Mpc^-1  (from Planck approx 67.80±0.77) CPLB: Use WMAP value
constexpr double H0SI                                   = H0 * 1000.0 / 3.0E22;                                     // Hubble's Constant in SI units, s^-1
constexpr double HUBBLE_TIME                            = 1 / H0SI;                                                 // Hubble time in s

constexpr double G                                      = 6.67E-11;                                                 // Gravitational constant in m^3 kg^-1 s^-2 (more accurately known as G M_sol)
constexpr double G_CGS                                  = 6.6743E-8;                                                // Gravitational constant in cm^3 g^-1 s^-2
constexpr double G_AU_Msol_yr                           = 4.0 * PI_2;                                               // Gravitational constant in AU^3 Msol^-1 yr^-2
constexpr double G_km_Msol_s                            = G * 1.0E-9 / KG_TO_MSOL;                                  // Gravitational constant in km^3 Msol^-1 s^-2
constexpr double G_SOLAR_YEAR                           = 3.14E7;                                                   // Gravitational constant in Lsol Rsol yr Msol^-2 for calculating photon tiring limit

constexpr double RSOL                                   = 6.957E8;                                                  // Solar Radius (in m)
constexpr double ZSOL                                   = 0.02;                                                     // Solar Metallicity used in scalings
constexpr double LOG10_ZSOL                             = -1.698970004336019;                                       // log10(ZSOL) - for performance
constexpr double ZSOL_ASPLUND                           = 0.0142;                                                   // Solar Metallicity (Asplund+ 2010) used in initial condition
constexpr double YSOL                                   = 0.2485;                                                   // Asplund+ 2009
constexpr double TSOL                                   = 5778.0;                                                   // Solar Temperature in kelvin
constexpr double LSOL                                   = 3.844E33;                                                 // Solar Luminosity in erg/s
constexpr double LSOLW                                  = 3.844E26;                                                 // Solar luminosity (in W)

constexpr double AU                                     = 149597870700.0;                                           // 1 AU (Astronomical Unit) in metres
constexpr double KM                                     = 1000.0;                                                   // 1 km (Kilometre) in metres
constexpr double C                                      = 299792458.0;                                              // Speed of light in m s^-1
constexpr double C_AU_yr                                = C / KM * KM_TO_AU * SECONDS_IN_YEAR;                      // Speed of light in AU yr^-1

constexpr double MU_0                                   = 4.0 * M_PI * 1.0E-7;                                      // Vacuum permeability in m kg s-2 A-2

constexpr double NEUTRINO_LOSS_FALLBACK_FACTOR          = 1.0;                                                      // Factor which accounts for mass loss in neutrino winds during a supernovae. Should be made a flag and added to pythonSubmit.py

constexpr double MC_L_C1                                = 9.20925E-5;                                               // Core Mass - Luminosity relation constant c1 (Hurley et al. 2000, eq 44)
constexpr double MC_L_C2                                = 5.402216;                                                 // Core Mass - Luminosity relation constant c2 (Hurley et al. 2000, eq 44)

constexpr double HE_RATE_CONSTANT                       = 7.66E-5;                                                  // Helium rate constant (Hurley et al. 2000, eq 68)
constexpr double HHE_RATE_CONSTANT                      = 1.27E-5;                                                  // Combined rate constant for both hydrogen and helium shell burning (Hurley et al. 2000, eq 71)

constexpr double EDDINGTON_PARAMETER_FACTOR             = HE_RATE_CONSTANT * 0.325;                                 // Eddington parameter factor

constexpr double BLACK_HOLE_LUMINOSITY                  = 1.0E-10;                                                  // Black Hole luminosity

constexpr double NEUTRON_STAR_MASS                      = 1.4;                                                      // Canonical NS mass in Msol
constexpr double NEUTRON_STAR_RADIUS                    = (1.0 / 7.0) * 1.0E-4;                                     // 10km in Rsol.  Hurley et al. 2000, just after eq 93

constexpr double MASSIVE_THRESHOLD                      = 8.0;                                                      // Mass (in solar masses) above which we consider stars to be "massive"
constexpr double HIGH_MASS_THRESHOLD                    = 12.0;                                                     // Mass (in solar masses) above which Hurley considers stars to be high mass stars
constexpr double VMS_MASS_THRESHOLD                     = 100.0;                                                    // Minimum mass for applying Very Massive (VMS) mass loss rates to be applied

constexpr double FRYER_PROTO_CORE_MASS_RAPID            = 1.0;                                                      // Proto neutron star core mass in the Fryer+ 2012 rapid prescription (Eq. 15)
constexpr double MCH                                    = 1.44;                                                     // Chandrasekhar mass
constexpr double MECS                                   = 1.38;                                                     // Mass of Neutron-Star (NS) formed in electron capture supernova (ECS). From Belczysnki+2008, before eq. 3.
constexpr double MECS_REM                               = 1.26;                                                     // Gravitational mass of Neutron-Star (NS) formed in electron capture supernova (ECS). From Belczysnki+2008, eq. 3
constexpr double MASS_LOSS_ETA                          = 0.5;                                                      // Mass loss efficiency -- can be set in the code as an option easily enough
constexpr double MCBUR1HURLEY					        = 1.6;							                            // Minimum core mass at base of the AGB to avoid fully degenerate CO core formation (Hurley value, Fryer+ and Belczynski+ use 1.83)
constexpr double MCBUR2					                = 2.25;							                            // Core mass at base of the AGB above which the CO core is completely non-degenerate

constexpr double NJ_MINIMUM_LUMINOSITY                  = 4.0E3;                                                    // Minimum luminosity in Lsun needed for Nieuwenhuijzen & de Jager wind mass loss
constexpr double VINK_MASS_LOSS_MINIMUM_TEMP            = 8.0E3;                                                    // Minimum temperature in K for Vink mass loss rates to be applied (12.5kK in Vink+Sander 2021)
constexpr double RSG_MAXIMUM_TEMP                       = 8.0E3;                                                    // Upper temperature in K for Red Supergiant (RSG) mass loss rates to be applied
constexpr double VINK_MASS_LOSS_BISTABILITY_TEMP        = 2.5E4;                                                    // Temperature in K for bistability jump in Vink mass loss (assumed to be 25000K following Belczysnki+2010)
constexpr double VINK_MASS_LOSS_MAXIMUM_TEMP            = 5.0E4;                                                    // Maximum temperature in K for Vink mass loss rates to be applied (show warning above this)
constexpr double LBV_LUMINOSITY_LIMIT_STARTRACK         = 6.0E5;                                                    // STARTRACK LBV luminosity limit
constexpr double LBV_LUMINOSITY_LIMIT_VANBEVEREN        = 3.0E5;                                                    // VANBEVEREN LBV luminosity limit

constexpr double CONVECTIVE_BOUNDARY_TEMPERATURE_BELCZYNSKI = 5.37E3;                                               // Threshold temperature for the star to develop a convective envelope, in Kelvin (10^3.73 K, from Belczynski+, 2008)

constexpr double MINIMUM_BLUE_LOOP_FRACTION             = 1.0E-10;                                                  // minimum fraction blue loop can be of He burning before we ignore it for HG radius calculation (see HG::CalculateRadiusOnPhase())

constexpr double TIMESTEP_QUANTUM                       = 1.0E-12;                                                  // Timestep quantum in Myr (=31.5576 seconds, given DAYS_IN_QUAD)
constexpr double ABSOLUTE_MINIMUM_TIMESTEP              = 3.0 * TIMESTEP_QUANTUM;                                   // In Myr (=94.6728 seconds, given TIMESTEP QUANTUM)
constexpr double NUCLEAR_MINIMUM_TIMESTEP               = 1.0E6 * TIMESTEP_QUANTUM;                                 // Minimum time step for nuclear evolution in My (= 1 year = 31557600 seconds, given TIMESTEP_QUANTUM)

constexpr unsigned int ABSOLUTE_MAXIMUM_TIMESTEPS       = 1000000;                                                  // Absolute maximum number of timesteps

constexpr int    MAX_BSE_INITIAL_CONDITIONS_ITERATIONS  = 100;                                                      // Maximum loop iterations looking for initial conditions for binary systems
constexpr int    MAX_TIMESTEP_RETRIES                   = 30;                                                       // Maximum retries to find a good timestep for stellar evolution

constexpr double MAXIMUM_MASS_LOSS_FRACTION             = 0.01;                                                     // Maximum allowable mass loss - 1.0% (of mass) expressed as a fraction
constexpr double MAXIMUM_RADIAL_CHANGE                  = 0.01;                                                     // Maximum allowable radial change - 1% (of radius) expressed as a fraction
constexpr double MINIMUM_MASS_SECONDARY                 = 4.0;                                                      // Minimum mass of secondary to evolve
constexpr double LAMBDA_NANJING_ZLIMIT                  = 0.0105;                                                   // Metallicity cutoff for Nanjing lambda calculations
constexpr double LAMBDA_NANJING_POPI_Z                  = 0.02;                                                     // Population I metallicity in Xu & Li (2010)
constexpr double LAMBDA_NANJING_POPII_Z                 = 0.001;                                                    // Population II metallicity in Xu & Li (2010)
constexpr double LAMBDA_NANJING_POPI_LOGZ               = -1.69897;                                                 // Population I log metallicity in Xu & Li (2010)
constexpr double LAMBDA_NANJING_POPII_LOGZ              = -3.0;                                                     // Population II log metallicity in Xu & Li (2010)
constexpr double LAMBDA_NANJING_MIN_MASS                = 1.0;                                                      // Minimum tabulated mass model in Xu & Li (2010)
constexpr double LAMBDA_NANJING_MAX_MASS                = 100.0;                                                    // Maximum tabulated mass model in Xu & Li (2010)

constexpr int    MAX_KEPLER_ITERATIONS                  = 1000;                                                     // Maximum number of iterations to solve Kepler's equation
constexpr double NEWTON_RAPHSON_EPSILON                 = 1.0E-5;                                                   // Accuracy for Newton-Raphson method

constexpr double EPSILON_PULSAR                         = 1.0;                                                      // Uncertainty due to mass loss

constexpr double MIN_HMXRB_STAR_TO_ROCHE_LOBE_RADIUS_RATIO = 0.8;                                                   // Minimum value of stellar radius | Roche Lobe radius for visible HMXRBs

constexpr double ADAPTIVE_RLOF_FRACTION_DONOR_GUESS     = 0.001;                                                    // Fraction of donor mass to use as guess in BaseBinaryStar::MassLossToFitInsideRocheLobe()
constexpr int    ADAPTIVE_RLOF_MAX_TRIES                = 30;                                                       // Maximum number of tries in BaseBinaryStar::MassLossToFitInsideRocheLobe()
constexpr int    ADAPTIVE_RLOF_MAX_ITERATIONS           = 50;                                                       // Maximum number of root finder iterations in BaseBinaryStar::MassLossToFitInsideRocheLobe()
constexpr double ADAPTIVE_RLOF_SEARCH_FACTOR_FRAC       = 1.0;                                                      // Search size factor (fractional part) in BaseBinaryStar::MassLossToFitInsideRocheLobe() (added to 1.0)

constexpr int    ADAPTIVE_MASS0_MAX_TRIES               = 30;                                                       // Maximum number of tries in HG::Mass0ToMatchDesiredCoreMass()
constexpr int    ADAPTIVE_MASS0_MAX_ITERATIONS          = 50;                                                       // Maximum number of iterations in HG::Mass0ToMatchDesiredCoreMass()
constexpr double ADAPTIVE_MASS0_SEARCH_FACTOR_FRAC      = 1.0;                                                      // Search size factor (fractional part) in HG::Mass0ToMatchDesiredCoreMass() (added to 1.0)

constexpr int    MULLERMANDEL_REMNANT_MASS_MAX_ITERATIONS = 1000;                                                   // Maximum number of iterations to find remnant mass in GiantBranch::CalculateRemnantMassByMullerMandel()

constexpr int    PULSAR_SPIN_ITERATIONS                 = 100;                                                      // Maximum number of iterations to find pulsar birth spin period in NS::CalculatePulsarBirthSpinPeriod()

constexpr int    SEMI_MAJOR_AXIS_SAMPLES                = 100;                                                      // Maximum number of samples when sampling period/semi-major axis in utils::SampleSemiMajorAxis()

constexpr int    TIDES_OMEGA_MAX_TRIES                  = 30;                                                       // Maximum number of tries in BaseBinaryStar::OmegaAfterCircularisation()
constexpr int    TIDES_OMEGA_MAX_ITERATIONS             = 50;                                                       // Maximum number of root finder iterations in BaseBinaryStar::OmegaAfterCircularisation()
constexpr double TIDES_OMEGA_SEARCH_FACTOR_FRAC         = 1.0;                                                      // Search size factor (fractional part) in BaseBinaryStar::OmegaAfterCircularisation() (added to 1.0)
constexpr double TIDES_MINIMUM_FRACTIONAL_EXTENT        = 1.0E-4;                                                   // Minimum fractional radius or mass of the stellar core or envelope, above which a given tidal dissipation mechanism is considered applicable
constexpr double TIDES_MAXIMUM_ORBITAL_CHANGE_FRAC      = 0.01;                                                     // Maximum allowed change in orbital and spin properties due to KAPIL2024 tides in a single timestep - 1% expressed as a fraction
constexpr double TIDES_MNIMUM_FRACTIONAL_NUCLEAR_TIME   = 0.001;                                                    // Minimum allowed timestep from tidal processes, as a fraction of the nuclear minimum time scale

constexpr double FARMER_PPISN_UPP_LIM_LIN_REGIME        = 38.0;                                                     // Maximum CO core mass to result in the linear remnant mass regime of the FARMER PPISN prescription
constexpr double FARMER_PPISN_UPP_LIM_QUAD_REGIME       = 60.0;                                                     // Maximum CO core mass to result in the quadratic remnant mass regime of the FARMER PPISN prescription
constexpr double FARMER_PPISN_UPP_LIM_INSTABILLITY      = 140.0;                                                    // Maximum CO core mass to result in PI (upper edge of PISN gap) from FARMER PPISN prescription
constexpr double STARTRACK_PPISN_HE_CORE_MASS           = 45.0;                                                     // Helium core mass remaining following PPISN as assumed in StarTrack (Belczynski et al. 2017 https://arxiv.org/abs/1607.03116)

constexpr double Q_CNO                                  = 9.9073E4;                                                 // Energy released per unit mass by hydrogen fusion via the CNO cycle in Lsol Myr Msol-1

// Initial mass of stars above which (including the limit) we allow convective core mass calculations from Shikauchi et al. (2024)
// Note that this value should always be > 0.7 Msol
constexpr double SHIKAUCHI_LOWER_MASS_LIMIT             = 15.0;

// logging constants

const LOGFILETYPE DEFAULT_LOGFILE_TYPE                  = LOGFILETYPE::HDF5;                                        // Default logfile type
const std::string DEFAULT_OUTPUT_CONTAINER_NAME         = "COMPAS_Output";                                          // Default name for output container (directory)
const std::string DEFAULT_HDF5_FILE_NAME                = "COMPAS_Output";                                          // Default name for HDF5 output file
const std::string DETAILED_OUTPUT_DIRECTORY_NAME        = "Detailed_Output";                                        // Name for detailed output directory within output container
const std::string RUN_DETAILS_FILE_NAME                 = "Run_Details";                                            // Name for run details output file within output container

constexpr int    HDF5_DEFAULT_CHUNK_SIZE                = 100000;                                                   // default HDF5 chunk size (number of dataset entries)
constexpr int    HDF5_DEFAULT_IO_BUFFER_SIZE            = 1;                                                        // number of HDF5 chunks to buffer for IO (per open dataset)
constexpr int    HDF5_MINIMUM_CHUNK_SIZE                = 1000;                                                     // minimum HDF5 chunk size (number of dataset entries)

// option constraints
// Use these constant to specify constraints that should be applied to program option values
// The values specified here should be checked in Options::OptionValues::CheckAndSetOptions()
// and in any relevant sampling functions

constexpr double MINIMUM_INITIAL_MASS                   = 0.00007;                                                  // Minimum initial mass (Msol) (~theoretical minimum? How low does COMPAS actually tolerate?)
constexpr double MAXIMUM_INITIAL_MASS                   = 150.0;                                                    // Maximum initial mass (Msol) (should actually be 100Msol?)

constexpr double MINIMUM_METALLICITY                    = 0.0001;                                                   // Minimum metallicity - Hurley equations known to fail for Z < 0.0001
constexpr double MAXIMUM_METALLICITY                    = 0.03;                                                     // Maximum metallicity (~> super-metal-rich?)


// IMF constants
constexpr double SALPETER_POWER                         = -2.35;
constexpr double SALPETER_MINIMUM                       = 0.5;
constexpr double SALPETER_MAXIMUM                       = 100.0;

// Kroupa IMF is a broken power law with three slopes
constexpr double KROUPA_POWER_1                         = -0.3;
constexpr double KROUPA_POWER_2                         = -1.3;
constexpr double KROUPA_POWER_3                         = -2.3;

// Declare some values here so we don't need to repeatedly calculate them in the code

// Often require the power law exponent plus one
constexpr double KROUPA_POWER_PLUS1_1                   = 0.7;
constexpr double KROUPA_POWER_PLUS1_2                   = -0.3;
constexpr double KROUPA_POWER_PLUS1_3                   = -1.3;

constexpr double ONE_OVER_KROUPA_POWER_1_PLUS1          = 1.0 / KROUPA_POWER_PLUS1_1;
constexpr double ONE_OVER_KROUPA_POWER_2_PLUS1          = 1.0 / KROUPA_POWER_PLUS1_2;
constexpr double ONE_OVER_KROUPA_POWER_3_PLUS1          = 1.0 / KROUPA_POWER_PLUS1_3;

// There are two breaks in the Kroupa power law -- they occur here (in solar masses)
constexpr double KROUPA_BREAK_1                         = 0.08;
constexpr double KROUPA_BREAK_2                         = 0.5;

// Some values that are really constants
constexpr double KROUPA_BREAK_1_PLUS1_1                 = 0.17067228025785934;                                      // pow(KROUPA_BREAK_1, KROUPA_POWER_PLUS1_1);
constexpr double KROUPA_BREAK_1_PLUS1_2                 = 2.13340350322324179;                                      // pow(KROUPA_BREAK_1, KROUPA_POWER_PLUS1_2);
constexpr double KROUPA_BREAK_1_POWER_1_2               = 0.08;                                                     // pow(KROUPA_BREAK_1, (KROUPA_POWER_1 - KROUPA_POWER_2));

constexpr double KROUPA_BREAK_2_PLUS1_2                 = 1.23114441334491628;                                      // pow(KROUPA_BREAK_2, KROUPA_POWER_PLUS1_2);
constexpr double KROUPA_BREAK_2_PLUS1_3                 = 2.46228882668983257;                                      // pow(KROUPA_BREAK_2, KROUPA_POWER_PLUS1_3);
constexpr double KROUPA_BREAK_2_POWER_2_3               = 0.5;                                                      // pow(KROUPA_BREAK_2, (KROUPA_POWER_2 - KROUPA_POWER_3));

constexpr double OPIKS_LAW_SEMIMAJOR_AXIS_DISTRIBUTION_POWER =  -1.0;

// Constants for the Muller and Mandel remnant mass and kick prescriptions
constexpr double MULLERMANDEL_M1                        = 2.0;	
constexpr double MULLERMANDEL_M2                        = 3.0; 
constexpr double MULLERMANDEL_M3                        = 7.0; 
constexpr double MULLERMANDEL_M4                        = 8.0; 
constexpr double MULLERMANDEL_MU1                       = 1.2;
constexpr double MULLERMANDEL_SIGMA1                    = 0.02;  
constexpr double MULLERMANDEL_MU2A                      = 1.4; 
constexpr double MULLERMANDEL_MU2B                      = 0.5;
constexpr double MULLERMANDEL_SIGMA2                    = 0.05;
constexpr double MULLERMANDEL_MU3A                      = 1.4;
constexpr double MULLERMANDEL_MU3B                      = 0.4;
constexpr double MULLERMANDEL_SIGMA3                    = 0.05;
constexpr double MULLERMANDEL_MUBH                    	= 0.8;
constexpr double MULLERMANDEL_SIGMABH                   = 0.5;
constexpr double MULLERMANDEL_MINNS                     = 1.13;
constexpr double MULLERMANDEL_KICKNS                    = 520.0;                                                    // As calibrated by Kapil+ 2023
constexpr double MULLERMANDEL_KICKBH                    = 200.0;
constexpr double MULLERMANDEL_SIGMAKICK                 = 0.3;

// Constants for the Maltsev+ 2024 SN remnant mass prescription
constexpr double MALTSEV2024_MMIN                       = 5.62;
constexpr double MALTSEV2024_MMAX                       = 16.18;
constexpr double MALTSEV2024_M1S                        = 6.6;
constexpr double MALTSEV2024_M1C                        = 6.6;
constexpr double MALTSEV2024_M1B                        = 7.7;
constexpr double MALTSEV2024_M1A                        = 7.4;
constexpr double MALTSEV2024_M2S                        = 7.2;
constexpr double MALTSEV2024_M2C                        = 7.1;
constexpr double MALTSEV2024_M2B                        = 8.3;
constexpr double MALTSEV2024_M2A                        = 8.4;
constexpr double MALTSEV2024_M3S                        = 12.9;
constexpr double MALTSEV2024_M3C                        = 13.2;
constexpr double MALTSEV2024_M3B                        = 15.2;
constexpr double MALTSEV2024_M3A                        = 15.4;
constexpr double MALTSEV2024_M1SZ01                     = 6.1;
constexpr double MALTSEV2024_M1CZ01                     = 6.3;
constexpr double MALTSEV2024_M1BZ01                     = 6.9;
constexpr double MALTSEV2024_M1AZ01                     = 7.0;
constexpr double MALTSEV2024_M2SZ01                     = 6.6;
constexpr double MALTSEV2024_M2CZ01                     = 7.1;
constexpr double MALTSEV2024_M2BZ01                     = 7.9;
constexpr double MALTSEV2024_M2AZ01                     = 7.4;
constexpr double MALTSEV2024_M3SZ01                     = 12.9;
constexpr double MALTSEV2024_M3CZ01                     = 12.3;
constexpr double MALTSEV2024_M3BZ01                     = 13.7;
constexpr double MALTSEV2024_M3AZ01                     = 13.7;

// Constants for WD evolution 

constexpr double COWD_LOG_MDOT_MIN_OFF_CENTER_IGNITION  = -5.688246139;                                             // Minimum log mass accretion rate for off center ignition in a CO WD. From Wang+ 2017. Log( 2.05 x 10^-6). 
constexpr double COWD_MASS_MIN_OFF_CENTER_IGNITION      = 1.33;                                                     // Minimum mass required for off center ignition, as shown in Wang, Podsiadlowski & Han (2017), sect 3.2.
constexpr double HEWD_HE_MDOT_CRIT                      = 2.0E-8;                                                   // Critical accretion rate for He WD accreting He-rich material. From Belczynski+ 2008, Mdot_crit2 in section 5.7.1.
constexpr double HEWD_MINIMUM_MASS_IGNITION             = 0.35;                                                     // Minimum mass for HeMS burning
constexpr double MASS_DOUBLE_DETONATION_CO              = 0.9;                                                      // Minimum mass for detonation which would yield something similar to SN Ia. Ruiter+ 2014.
constexpr double Q_HYDROGEN_BURNING                     = 6.4E18 * MSOL_TO_G / (SECONDS_IN_YEAR * LSOL);            // 6.4E18 is the energy yield of H burning in erg/g as given in Nomoto+ 2007 (2007ApJ...663.1269N)
constexpr double WD_HE_SHELL_MCRIT_DETONATION           = 0.05;                                                     // Minimum shell mass of He for detonation. Should be composed of helium (so, exclude burnt material), but not implemented yet. Ruiter+ 2014.
constexpr double WD_LOG_MT_LIMIT_PIERSANTI_RG_SS_0      = -6.84;
constexpr double WD_LOG_MT_LIMIT_PIERSANTI_RG_SS_1      = 1.349;
constexpr double WD_LOG_MT_LIMIT_PIERSANTI_SS_MF_0      = -8.115;
constexpr double WD_LOG_MT_LIMIT_PIERSANTI_SS_MF_1      = 2.29;
constexpr double WD_LOG_MT_LIMIT_PIERSANTI_SF_Dt_0      = -8.313;
constexpr double WD_LOG_MT_LIMIT_PIERSANTI_SF_Dt_1      = 1.018;
constexpr double WD_LOG_MT_LIMIT_NOMOTO_REDGIANT_0      = -8.33017155;
constexpr double WD_LOG_MT_LIMIT_NOMOTO_REDGIANT_1      = 2.88247131;
constexpr double WD_LOG_MT_LIMIT_NOMOTO_REDGIANT_2      = -0.98023471;
constexpr double WD_LOG_MT_LIMIT_NOMOTO_STABLE_0        = -9.21757267;
constexpr double WD_LOG_MT_LIMIT_NOMOTO_STABLE_1        = 3.57319872;
constexpr double WD_LOG_MT_LIMIT_NOMOTO_STABLE_2        = -1.2137735;


// coefficients for the calculation of initial angular frequency for Chemically Homogeneous Evolution
// Mandel from Butler 2018
const DBL_VECTOR CHE_Coefficients = { 5.7914E-04, -1.9196E-06, -4.0602E-07, 1.0150E-08, -9.1792E-11, 2.9051E-13 };

// WD effective baryon number lookup table
// unordered_map - key is stellar type
// Hurley et al. 2000, just after eq 90
const COMPASUnorderedMap<STELLAR_TYPE, double> WD_Baryon_Number = {
    {STELLAR_TYPE::HELIUM_WHITE_DWARF,         4.0},
    {STELLAR_TYPE::CARBON_OXYGEN_WHITE_DWARF, 15.0},
    {STELLAR_TYPE::OXYGEN_NEON_WHITE_DWARF,   17.0}
};

// symbolic names for term coefficients for Luminosity & Radius coefficients from Tout et al. 1996
enum class LR_TCoeff: int { a, b, c, d, e };
#define A LR_TCoeff::a
#define B LR_TCoeff::b
#define C LR_TCoeff::c
#define D LR_TCoeff::d
#define E LR_TCoeff::e

// symbolic names for luminosity coefficients (from Tout et al. 1996)
enum class L_Coeff: int { ALPHA, BETA, GAMMA, DELTA, EPSILON, ZETA, ETA };
#define ALPHA   L_Coeff::ALPHA
#define BETA    L_Coeff::BETA
#define GAMMA   L_Coeff::GAMMA
#define DELTA   L_Coeff::DELTA
#define EPSILON L_Coeff::EPSILON
#define ZETA    L_Coeff::ZETA
#define ETA     L_Coeff::ETA

// L (Luminosity) coefficients
// Table 1 in Tout et al 1996
// Key to map is L_Coeff.  Map element is unordered_map of term coefficient values.
const std::map<int, COMPASUnorderedMap<LR_TCoeff, double>> L_COEFF = {
    {static_cast<int>(ALPHA),   {{A, 0.39704170}, {B,  -0.32913574}, {C,  0.34776688}, {D,  0.37470851}, {E, 0.09011915}}},
    {static_cast<int>(BETA),    {{A, 8.52762600}, {B, -24.41225973}, {C, 56.43597107}, {D, 37.06152575}, {E, 5.45624060}}},
    {static_cast<int>(GAMMA),   {{A, 0.00025546}, {B,  -0.00123461}, {C, -0.00023246}, {D,  0.00045519}, {E, 0.00016176}}},
    {static_cast<int>(DELTA),   {{A, 5.43288900}, {B,  -8.62157806}, {C, 13.44202049}, {D, 14.51584135}, {E, 3.39793084}}},
    {static_cast<int>(EPSILON), {{A, 5.56357900}, {B, -10.32345224}, {C, 19.44322980}, {D, 18.97361347}, {E, 4.16903097}}},
    {static_cast<int>(ZETA),    {{A, 0.78866060}, {B,  -2.90870942}, {C,  6.54713531}, {D,  4.05606657}, {E, 0.53287322}}},
    {static_cast<int>(ETA),     {{A, 0.00586685}, {B,  -0.01704237}, {C,  0.03872348}, {D,  0.02570041}, {E, 0.00383376}}}
};

#undef ALPHA
#undef BETA
#undef GAMMA
#undef DELTA
#undef EPSILON
#undef ZETA
#undef ETA

// symbolic names for radius coefficients from Tout et al. 1996
enum class R_Coeff: int { THETA, IOTA, KAPPA, LAMBDA, MU, NU, XI, OMICRON, PI };
#define THETA   R_Coeff::THETA
#define IOTA    R_Coeff::IOTA
#define KAPPA   R_Coeff::KAPPA
#define LAMBDA  R_Coeff::LAMBDA
#define MU      R_Coeff::MU
#define NU      R_Coeff::NU
#define XI      R_Coeff::XI
#define OMICRON R_Coeff::OMICRON
#define Pi      R_Coeff::PI

// R (Radius) coefficients
// Table 2 in Tout et al. 1996
// Key to map is L_Coeff.  Map element is unordered_map of term coefficient values.
const std::map<int, COMPASUnorderedMap<LR_TCoeff, double>> R_COEFF = {
    {static_cast<int>(THETA),   {{A,  1.71535900}, {B,  0.62246212}, {C,  -0.92557761}, {D,  -1.16996966}, {E, -0.30631491}}},
    {static_cast<int>(IOTA),    {{A,  6.59778800}, {B, -0.42450044}, {C, -12.13339427}, {D, -10.73509484}, {E, -2.51487077}}},
    {static_cast<int>(KAPPA),   {{A, 10.08855000}, {B, -7.11727086}, {C, -31.67119479}, {D, -24.24848322}, {E, -5.33608972}}},
    {static_cast<int>(LAMBDA),  {{A,  1.01249500}, {B,  0.32699690}, {C,  -0.00923418}, {D,  -0.03876858}, {E, -0.00412750}}},
    {static_cast<int>(MU),      {{A,  0.07490166}, {B,  0.02410413}, {C,   0.07233664}, {D,   0.03040467}, {E,  0.00197741}}},
    {static_cast<int>(NU),      {{A,  0.01077422}, {B,  0.00000000}, {C,   0.00000000}, {D,   0.00000000}, {E,  0.00000000}}},
    {static_cast<int>(XI),      {{A,  3.08223400}, {B,  0.94472050}, {C,  -2.15200882}, {D,  -2.49219496}, {E, -0.63848738}}},
    {static_cast<int>(OMICRON), {{A, 17.84778000}, {B, -7.45345690}, {C, -48.96066856}, {D, -40.05386135}, {E, -9.09331816}}},
    {static_cast<int>(Pi),      {{A,  0.00022582}, {B, -0.00186899}, {C,   0.00388783}, {D,   0.00142402}, {E, -0.00007671}}}
};

#undef THETA
#undef IOTA
#undef KAPPA
#undef LAMBDA
#undef MU
#undef NU
#undef XI
#undef OMICRON
#undef Pi

#undef A
#undef B
#undef C
#undef D
#undef E

// symbolic names for term coefficients for A & B coefficients from Hurley et al. 2000
enum class AB_TCoeff: int { ALPHA, BETA, GAMMA, ETA, MU };
#define ALPHA AB_TCoeff::ALPHA
#define BETA  AB_TCoeff::BETA
#define GAMMA AB_TCoeff::GAMMA
#define ETA   AB_TCoeff::ETA
#define MU    AB_TCoeff::MU

// A coefficients
// Table in Appendix A of Hurley et al. 2000
// Key to map is n (A(n)).  Map element is unordered_map of term coefficient values.
const std::map<int, COMPASUnorderedMap<AB_TCoeff, double>> A_COEFF = {
    { 1, {{ALPHA,  1.593890E3 }, {BETA,  2.053038E3 }, {GAMMA,  1.231226E3 }, {ETA,  2.327785E2 }, {MU,  0.000000E0 }}},
    { 2, {{ALPHA,  2.706708E3 }, {BETA,  1.483131E3 }, {GAMMA,  5.772723E2 }, {ETA,  7.411230E1 }, {MU,  0.000000E0 }}},
    { 3, {{ALPHA,  1.466143E2 }, {BETA, -1.048442E2 }, {GAMMA, -6.795374E1 }, {ETA, -1.391127E1 }, {MU,  0.000000E0 }}},
    { 4, {{ALPHA,  4.141960E-2}, {BETA,  4.564888E-2}, {GAMMA,  2.958542E-2}, {ETA,  5.571483E-3}, {MU,  0.000000E0 }}},
    { 5, {{ALPHA,  3.426349E-1}, {BETA,  0.000000E0 }, {GAMMA,  0.000000E0 }, {ETA,  0.000000E0 }, {MU,  0.000000E0 }}},
    { 6, {{ALPHA,  1.949814E1 }, {BETA,  1.758178E0 }, {GAMMA, -6.008212E0 }, {ETA, -4.470533E0 }, {MU,  0.000000E0 }}},
    { 7, {{ALPHA,  4.903830E0 }, {BETA,  0.000000E0 }, {GAMMA,  0.000000E0 }, {ETA,  0.000000E0 }, {MU,  0.000000E0 }}},
    { 8, {{ALPHA,  5.212154E-2}, {BETA,  3.166411E-2}, {GAMMA, -2.750074E-3}, {ETA, -2.271549E-3}, {MU,  0.000000E0 }}},
    { 9, {{ALPHA,  1.312179E0 }, {BETA, -3.294936E-1}, {GAMMA,  9.231860E-2}, {ETA,  2.610989E-2}, {MU,  0.000000E0 }}},
    {10, {{ALPHA,  8.073972E-1}, {BETA,  0.000000E0 }, {GAMMA,  0.000000E0 }, {ETA,  0.000000E0 }, {MU,  0.000000E0 }}},

    {11, {{ALPHA,  1.031538E0 }, {BETA, -2.434480E-1}, {GAMMA,  7.732821E0 }, {ETA,  6.460705E0 }, {MU,  1.374484E0 }}},
    {12, {{ALPHA,  1.043715E0 }, {BETA, -1.577474E0 }, {GAMMA, -5.168234E0 }, {ETA, -5.596506E0 }, {MU, -1.299394E0 }}},
    {13, {{ALPHA,  7.859573E2 }, {BETA, -8.542048E0 }, {GAMMA, -2.642511E1 }, {ETA, -9.585707E0 }, {MU,  0.000000E0 }}},
    {14, {{ALPHA,  3.858911E3 }, {BETA,  2.459681E3 }, {GAMMA, -7.630093E1 }, {ETA, -3.486057E2 }, {MU, -4.861703E1 }}},
    {15, {{ALPHA,  2.888720E2 }, {BETA,  2.952979E2 }, {GAMMA,  1.850341E2 }, {ETA,  3.797254E1 }, {MU,  0.000000E0 }}},
    {16, {{ALPHA,  7.196580E0 }, {BETA,  5.613746E-1}, {GAMMA,  3.805871E-1}, {ETA,  8.398728E-2}, {MU,  0.000000E0 }}},
    {17, {{ALPHA,  0.000000E0 }, {BETA,  0.000000E0 }, {GAMMA,  0.000000E0 }, {ETA,  0.000000E0 }, {MU,  0.000000E0 }}},
    {18, {{ALPHA,  2.187715E-1}, {BETA, -2.154437E0 }, {GAMMA, -3.768678E0 }, {ETA, -1.975518E0 }, {MU, -3.021475E-1}}},
    {19, {{ALPHA,  1.466440E0 }, {BETA,  1.839725E0 }, {GAMMA,  6.442199E0 }, {ETA,  4.023635E0 }, {MU,  6.957529E-1}}},
    {20, {{ALPHA,  2.652091E1 }, {BETA,  8.178458E1 }, {GAMMA,  1.156058E2 }, {ETA,  7.633811E1 }, {MU,  1.950698E1 }}},

    {21, {{ALPHA,  1.472103E0 }, {BETA, -2.947609E0 }, {GAMMA, -3.312828E0 }, {ETA, -9.945065E-1}, {MU,  0.000000E0 }}},
    {22, {{ALPHA,  3.071048E0 }, {BETA, -5.679941E0 }, {GAMMA, -9.745523E0 }, {ETA, -3.594543E0 }, {MU,  0.000000E0 }}},
    {23, {{ALPHA,  2.617890E0 }, {BETA,  1.019135E0 }, {GAMMA, -3.292551E-2}, {ETA, -7.445123E-2}, {MU,  0.000000E0 }}},
    {24, {{ALPHA,  1.075567E-2}, {BETA,  1.773287E-2}, {GAMMA,  9.610479E-3}, {ETA,  1.732469E-3}, {MU,  0.000000E0 }}},
    {25, {{ALPHA,  1.476246E0 }, {BETA,  1.899331E0 }, {GAMMA,  1.195010E0 }, {ETA,  3.035051E-1}, {MU,  0.000000E0 }}},
    {26, {{ALPHA,  5.502535E0 }, {BETA, -6.601663E-2}, {GAMMA,  9.968707E-2}, {ETA,  3.599801E-2}, {MU,  0.000000E0 }}},
    {27, {{ALPHA,  9.511033E1 }, {BETA,  6.819618E1 }, {GAMMA, -1.045625E1 }, {ETA, -1.474939E1 }, {MU,  0.000000E0 }}},
    {28, {{ALPHA,  3.113458E1 }, {BETA,  1.012033E1 }, {GAMMA, -4.650511E0 }, {ETA, -2.463185E0 }, {MU,  0.000000E0 }}},
    {29, {{ALPHA,  1.413057E0 }, {BETA,  4.578814E-1}, {GAMMA, -6.850581E-2}, {ETA, -5.588658E-2}, {MU,  0.000000E0 }}},
    {30, {{ALPHA,  3.910862E1 }, {BETA,  5.196646E1 }, {GAMMA,  2.264970E1 }, {ETA,  2.873680E0 }, {MU,  0.000000E0 }}},

    {31, {{ALPHA,  4.597479E0 }, {BETA, -2.855179E-1}, {GAMMA,  2.709724E-1}, {ETA,  0.000000E0 }, {MU,  0.000000E0 }}},
    {32, {{ALPHA,  6.682518E0 }, {BETA,  2.827718E-1}, {GAMMA, -7.294429E-2}, {ETA,  0.000000E0 }, {MU,  0.000000E0 }}},
    {33, {{ALPHA,  0.000000E0 }, {BETA,  0.000000E0 }, {GAMMA,  0.000000E0 }, {ETA,  0.000000E0 }, {MU,  0.000000E0 }}},
    {34, {{ALPHA,  1.910302E-1}, {BETA,  1.158624E-1}, {GAMMA,  3.348990E-2}, {ETA,  2.599706E-3}, {MU,  0.000000E0 }}},
    {35, {{ALPHA,  3.931056E-1}, {BETA,  7.277637E-2}, {GAMMA, -1.366593E-1}, {ETA, -4.508946E-2}, {MU,  0.000000E0 }}},
    {36, {{ALPHA,  3.267776E-1}, {BETA,  1.204424E-1}, {GAMMA,  9.988332E-2}, {ETA,  2.455361E-2}, {MU,  0.000000E0 }}},
    {37, {{ALPHA,  5.990212E-1}, {BETA,  5.570264E-2}, {GAMMA,  6.207626E-2}, {ETA,  1.777283E-2}, {MU,  0.000000E0 }}},
    {38, {{ALPHA,  7.330122E-1}, {BETA,  5.192827E-1}, {GAMMA,  2.316416E-1}, {ETA,  8.346941E-3}, {MU,  0.000000E0 }}},
    {39, {{ALPHA,  1.172768E0 }, {BETA, -1.209262E-1}, {GAMMA, -1.193023E-1}, {ETA, -2.859837E-2}, {MU,  0.000000E0 }}},
    {40, {{ALPHA,  3.982622E-1}, {BETA, -2.296279E-1}, {GAMMA, -2.262539E-1}, {ETA, -5.219837E-2}, {MU,  0.000000E0 }}},

    {41, {{ALPHA,  3.571038E0 }, {BETA, -2.223635E-2}, {GAMMA, -2.611794E-2}, {ETA, -6.359648E-3}, {MU,  0.000000E0 }}},
    {42, {{ALPHA,  1.984800E0 }, {BETA,  1.138600E0 }, {GAMMA,  3.564000E-1}, {ETA,  0.000000E0 }, {MU,  0.000000E0 }}},
    {43, {{ALPHA,  6.300000E-2}, {BETA,  4.810000E-2}, {GAMMA,  9.840000E-3}, {ETA,  0.000000E0 }, {MU,  0.000000E0 }}},
    {44, {{ALPHA,  1.200000E0 }, {BETA,  2.450000E0 }, {GAMMA,  0.000000E0 }, {ETA,  0.000000E0 }, {MU,  0.000000E0 }}},
    {45, {{ALPHA,  2.321400E-1}, {BETA,  1.828075E-3}, {GAMMA, -2.232007E-2}, {ETA, -3.378734E-3}, {MU,  0.000000E0 }}},
    {46, {{ALPHA,  1.163659E-2}, {BETA,  3.427682E-3}, {GAMMA,  1.421393E-3}, {ETA, -3.710666E-3}, {MU,  0.000000E0 }}},
    {47, {{ALPHA,  1.048020E-2}, {BETA, -1.231921E-2}, {GAMMA, -1.686860E-2}, {ETA, -4.234354E-3}, {MU,  0.000000E0 }}},
    {48, {{ALPHA,  1.555590E0 }, {BETA, -3.223927E-1}, {GAMMA, -5.197429E-1}, {ETA, -1.066441E-1}, {MU,  0.000000E0 }}},
    {49, {{ALPHA,  9.770000E-2}, {BETA, -2.310000E-1}, {GAMMA, -7.530000E-2}, {ETA,  0.000000E0 }, {MU,  0.000000E0 }}},
    {50, {{ALPHA,  2.400000E-1}, {BETA,  1.800000E-1}, {GAMMA,  5.950000E-1}, {ETA,  0.000000E0 }, {MU,  0.000000E0 }}},

    {51, {{ALPHA,  3.300000E-1}, {BETA,  1.320000E-1}, {GAMMA,  2.180000E-1}, {ETA,  0.000000E0 }, {MU,  0.000000E0 }}},
    {52, {{ALPHA,  1.106400E0 }, {BETA,  4.150000E-1}, {GAMMA,  1.800000E-1}, {ETA,  0.000000E0 }, {MU,  0.000000E0 }}},
    {53, {{ALPHA,  1.190000E0 }, {BETA,  3.770000E-1}, {GAMMA,  1.760000E-1}, {ETA,  0.000000E0 }, {MU,  0.000000E0 }}},
    {54, {{ALPHA,  3.855707E-1}, {BETA, -6.104166E-1}, {GAMMA,  5.676742E0 }, {ETA,  1.060894E1 }, {MU,  5.284014E0 }}},
    {55, {{ALPHA,  3.579064E-1}, {BETA, -6.442936E-1}, {GAMMA,  5.494644E0 }, {ETA,  1.054952E1 }, {MU,  5.280991E0 }}},
    {56, {{ALPHA,  9.587587E-1}, {BETA,  8.777464E-1}, {GAMMA,  2.017321E-1}, {ETA,  0.000000E0 }, {MU,  0.000000E0 }}},
    {57, {{ALPHA,  1.513500E0 }, {BETA,  3.769000E-1}, {GAMMA,  0.000000E0 }, {ETA,  0.000000E0 }, {MU,  0.000000E0 }}},
    {58, {{ALPHA,  4.907546E-1}, {BETA, -1.683928E-1}, {GAMMA, -3.108742E-1}, {ETA, -7.202918E-2}, {MU,  0.000000E0 }}},
    {59, {{ALPHA,  4.537070E0 }, {BETA, -4.465455E0 }, {GAMMA, -1.612690E0 }, {ETA, -1.623246E0 }, {MU,  0.000000E0 }}},
    {60, {{ALPHA,  1.796220E0 }, {BETA,  2.814020E-1}, {GAMMA,  1.423325E0 }, {ETA,  3.421036E-1}, {MU,  0.000000E0 }}},

    {61, {{ALPHA,  2.256216E0 }, {BETA,  3.773400E-1}, {GAMMA,  1.537867E0 }, {ETA,  4.396373E-1}, {MU,  0.000000E0 }}},
    {62, {{ALPHA,  8.430000E-2}, {BETA, -4.750000E-2}, {GAMMA, -3.520000E-2}, {ETA,  0.000000E0 }, {MU,  0.000000E0 }}},
    {63, {{ALPHA,  7.360000E-2}, {BETA,  7.490000E-2}, {GAMMA,  4.426000E-2}, {ETA,  0.000000E0 }, {MU,  0.000000E0 }}},
    {64, {{ALPHA,  1.360000E-1}, {BETA,  3.520000E-2}, {GAMMA,  0.000000E0 }, {ETA,  0.000000E0 }, {MU,  0.000000E0 }}},
    {65, {{ALPHA,  1.564231E-3}, {BETA,  1.653042E-3}, {GAMMA, -4.439786E-3}, {ETA, -4.951011E-3}, {MU, -1.216530E-3}}},
    {66, {{ALPHA,  1.477000E0 }, {BETA,  2.960000E-1}, {GAMMA,  0.000000E0 }, {ETA,  0.000000E0 }, {MU,  0.000000E0 }}},
    {67, {{ALPHA,  5.210157E0 }, {BETA, -4.143695E0 }, {GAMMA, -2.120870E0 }, {ETA,  0.000000E0 }, {MU,  0.000000E0 }}},
    {68, {{ALPHA,  1.116000E0 }, {BETA,  1.660000E-1}, {GAMMA,  0.000000E0 }, {ETA,  0.000000E0 }, {MU,  0.000000E0 }}},
    {69, {{ALPHA,  1.071489E0 }, {BETA, -1.164852E-1}, {GAMMA, -8.623831E-2}, {ETA, -1.582349E-2}, {MU,  0.000000E0 }}},
    {70, {{ALPHA,  7.108492E-1}, {BETA,  7.935927E-1}, {GAMMA,  3.926983E-1}, {ETA,  3.622146E-2}, {MU,  0.000000E0 }}},

    {71, {{ALPHA,  3.478514E0 }, {BETA, -2.585474E-2}, {GAMMA, -1.512955E-2}, {ETA, -2.833691E-3}, {MU,  0.000000E0 }}},
    {72, {{ALPHA,  9.132108E-1}, {BETA, -1.653695E-1}, {GAMMA,  0.000000E0 }, {ETA,  3.636784E-2}, {MU,  0.000000E0 }}},
    {73, {{ALPHA,  3.969331E-3}, {BETA,  4.539076E-3}, {GAMMA,  1.720906E-3}, {ETA,  1.897857E-4}, {MU,  0.000000E0 }}},
    {74, {{ALPHA,  1.600000E0 }, {BETA,  7.640000E-1}, {GAMMA,  3.322000E-1}, {ETA,  0.000000E0 }, {MU,  0.000000E0 }}},
    {75, {{ALPHA,  8.109000E-1}, {BETA, -6.282000E-1}, {GAMMA,  0.000000E0 }, {ETA,  0.000000E0 }, {MU,  0.000000E0 }}},
    {76, {{ALPHA,  1.192334E-2}, {BETA,  1.083057E-2}, {GAMMA,  1.230969E0 }, {ETA,  1.551656E0 }, {MU,  0.000000E0 }}},
    {77, {{ALPHA, -1.668868E-1}, {BETA,  5.818123E-1}, {GAMMA, -1.105027E1 }, {ETA, -1.668070E1 }, {MU,  0.000000E0 }}},
    {78, {{ALPHA,  7.615495E-1}, {BETA,  1.068243E-1}, {GAMMA, -2.011333E-1}, {ETA, -9.371415E-2}, {MU,  0.000000E0 }}},
    {79, {{ALPHA,  9.409838E0 }, {BETA,  1.522928E0 }, {GAMMA,  0.000000E0 }, {ETA,  0.000000E0 }, {MU,  0.000000E0 }}},
    {80, {{ALPHA, -2.711000E-1}, {BETA, -5.756000E-1}, {GAMMA, -8.380000E-2}, {ETA,  0.000000E0 }, {MU,  0.000000E0 }}},

    {81, {{ALPHA,  2.493000E0 }, {BETA,  1.147500E0 }, {GAMMA,  0.000000E0 }, {ETA,  0.000000E0 }, {MU,  0.000000E0 }}}
};

// Critial mass ratios for a grid of masses and radii. These come from the team of Hongwei Ge, in a series of papers 
// titled "Adiabatic Mass Loss". The first two tables below are a subset of Tables A3 and A4 from Ge et al. 2024 
// (Paper V, arXiv:2408.16350), provided by private request to Hongwei Ge. This is the most up to date version of these 
// tables at the time of writing. Table A3 corresponds to low metallicity Z=0.001, Table A4 is for Z=0.02. Both of 
// these tables are represented below as two separate nested vectors, and we interpolate linearly in logM, logR, and logZ. 
//
// Both tables contain the Mass/Msol, log(R/Rsol), and critical mass ratios for mass transfer instability for a variety of model 
// variations, namely qCritST_full, qCritST_half, qCritST_nonc, qCritIC_full, qCritIC_half, qCritIC_nonc.
// The ST and IC suffixes distinguish their adiabatic and isentropic models, respectively (see Ge et al. 2020).
// The _full, _half, and _nonc suffixes distinguish fully conservative MT, half conservative MT, and fully non-conservative MT, 
// respectively, where non-conservative models assume isotropic re-emission AM loss. 
//
// The First entry in the tuple is the vector of unique mass values, second entry is the vector containing 7-tuples of vectors 
// for logR, and the six qCrits listed above (in order). Note that the radius may contract several times. These points have 
// been removed to facilitate the interpolation, so logR is monotonic.
//
// The third table below comes from Zhang et al. 2024 (Paper IV, arxiv:2406.13146). This table was computed exclusively for He
// stars, and does include the same variations as for the H-rich donors in the first two tables. There is only one critical mass
// ratio value per mass and radius, assuming fully conservative MT and a purely adiabatic (i.e non-istentropic) response to mass loss.
//
// In all cases, q is mAccretor/mDonor, which is inverted from the Ge et al. datatable. 
//
// Low Z = 0.001 table
const GE_QCRIT_TABLE QCRIT_GE_LOW_Z = {
    {1.0, 1.14, 1.3, 1.44, 1.6, 1.8, 2.0, 2.5, 3.2, 4.0, 5.0, 6.3, 8.0, 10.0, 13.0, 16.0, 20.0, 25.0, 32.0, 40.0, 50.0, 63.0, 80.0, 100.0},
    {
      {{-0.0391, -0.029, -0.0112, 0.0059, 0.0217, 0.0359, 0.0492, 0.0663, 0.1844, 0.2985, 0.4157, 0.5323, 0.6365, 0.7326, 0.8349, 0.9336, 1.0333, 1.1347, 1.2364, 1.3386, 1.4372, 1.5371, 1.6372, 1.7365, 1.834, 1.9394, 2.0102, 2.0483, 2.1567, 2.2588},  {0.425, 0.412, 0.392, 0.375, 0.360, 0.348, 0.338, 0.326, 0.269, 0.239, 0.242, 0.962, 1.085, 1.101, 1.092, 1.071, 1.042, 1.007, 0.970, 0.928, 0.881, 0.827, 0.767, 0.700, 0.625, 0.532, 0.463, 0.286, 0.221, 0.138},  {0.425, 0.414, 0.398, 0.384, 0.372, 0.361, 0.352, 0.342, 0.292, 0.265, 0.264, 0.817, 0.912, 0.926, 0.919, 0.903, 0.882, 0.856, 0.829, 0.797, 0.762, 0.720, 0.675, 0.623, 0.564, 0.491, 0.435, 0.286, 0.332, 0.332},  {0.455, 0.446, 0.433, 0.422, 0.412, 0.404, 0.396, 0.387, 0.345, 0.322, 0.315, 0.693, 0.768, 0.779, 0.775, 0.765, 0.750, 0.732, 0.713, 0.691, 0.666, 0.637, 0.603, 0.565, 0.521, 0.465, 0.423, 0.296, 0.340, 0.340},  {0.424, 0.411, 0.390, 0.373, 0.359, 0.347, 0.336, 0.324, 0.267, 0.237, 0.241, 0.913, 1.022, 1.029, 1.007, 0.974, 0.932, 0.883, 0.831, 0.773, 0.711, 0.643, 0.570, 0.492, 0.410, 0.318, 0.255, 0.149, 0.088, 0.022},  {0.424, 0.414, 0.397, 0.382, 0.370, 0.360, 0.351, 0.341, 0.291, 0.264, 0.263, 0.779, 0.867, 0.873, 0.858, 0.833, 0.801, 0.765, 0.726, 0.682, 0.635, 0.582, 0.525, 0.464, 0.398, 0.323, 0.270, 0.172, 0.117, 0.046},  {0.454, 0.446, 0.433, 0.421, 0.411, 0.403, 0.395, 0.386, 0.344, 0.321, 0.314, 0.668, 0.737, 0.743, 0.734, 0.717, 0.696, 0.671, 0.644, 0.613, 0.579, 0.541, 0.499, 0.453, 0.403, 0.344, 0.303, 0.209, 0.164, 0.099}},
      {{-0.0079, -0.0012, 0.0124, 0.0273, 0.0429, 0.0606, 0.0811, 0.1071, 0.1195, 0.1198, 0.2209, 0.3222, 0.4201, 0.5206, 0.6209, 0.7201, 0.8201, 0.9239, 1.0276, 1.1203, 1.2281, 1.3198, 1.4194, 1.5285, 1.627, 1.725, 1.8217, 1.9266, 2.0076, 2.0544, 2.1579, 2.2587, 2.3532},  {0.410, 0.401, 0.386, 0.372, 0.360, 0.348, 0.337, 0.326, 0.318, 0.309, 0.272, 0.246, 0.229, 0.250, 1.020, 1.094, 1.099, 1.082, 1.054, 1.022, 0.982, 0.950, 0.905, 0.847, 0.787, 0.721, 0.646, 0.553, 0.478, 0.363, 0.287, 0.201, 0.104},  {0.414, 0.407, 0.394, 0.382, 0.371, 0.361, 0.351, 0.342, 0.335, 0.327, 0.294, 0.272, 0.256, 0.269, 0.863, 0.921, 0.926, 0.913, 0.892, 0.869, 0.840, 0.816, 0.781, 0.738, 0.693, 0.642, 0.584, 0.511, 0.451, 0.352, 0.352, 0.352, 0.352},  {0.448, 0.441, 0.431, 0.421, 0.412, 0.403, 0.395, 0.387, 0.381, 0.374, 0.347, 0.328, 0.314, 0.317, 0.730, 0.776, 0.781, 0.773, 0.759, 0.743, 0.723, 0.707, 0.683, 0.652, 0.620, 0.583, 0.541, 0.486, 0.441, 0.353, 0.414, 0.414, 0.414},  {0.409, 0.400, 0.385, 0.372, 0.359, 0.347, 0.336, 0.325, 0.317, 0.309, 0.271, 0.245, 0.227, 0.248, 0.974, 1.036, 1.031, 1.002, 0.962, 0.918, 0.865, 0.820, 0.760, 0.688, 0.617, 0.541, 0.460, 0.367, 0.295, 0.201, 0.134, 0.058, 0.015},  {0.413, 0.406, 0.393, 0.381, 0.371, 0.360, 0.351, 0.341, 0.334, 0.327, 0.294, 0.271, 0.270, 0.268, 0.828, 0.879, 0.876, 0.855, 0.824, 0.792, 0.752, 0.718, 0.674, 0.619, 0.564, 0.505, 0.441, 0.366, 0.307, 0.221, 0.163, 0.090, 0.039},  {0.447, 0.441, 0.430, 0.420, 0.411, 0.403, 0.394, 0.386, 0.381, 0.374, 0.346, 0.327, 0.313, 0.316, 0.708, 0.748, 0.748, 0.734, 0.714, 0.692, 0.664, 0.641, 0.610, 0.571, 0.532, 0.488, 0.440, 0.383, 0.337, 0.256, 0.211, 0.150, 0.098}},
      {{0.0016, 0.0098, 0.0252, 0.0435, 0.0701, 0.0923, 0.1267, 0.1701, 0.1848, 0.1849, 0.2848, 0.3851, 0.4844, 0.5839, 0.6862, 0.7854, 0.882, 0.9848, 1.0819, 1.185, 1.2822, 1.3834, 1.4892, 1.5857, 1.682, 1.7881, 1.8817, 1.9827, 2.0905, 2.1867, 2.2913, 2.3847},  {0.417, 0.407, 0.392, 0.378, 0.361, 0.350, 0.336, 0.322, 0.315, 0.298, 0.270, 0.249, 0.233, 0.237, 1.054, 1.120, 1.122, 1.103, 1.075, 1.041, 1.005, 0.965, 0.915, 0.863, 0.804, 0.730, 0.655, 0.562, 0.412, 0.327, 0.227, 0.118},  {0.420, 0.411, 0.398, 0.386, 0.372, 0.362, 0.351, 0.338, 0.332, 0.317, 0.293, 0.275, 0.260, 0.261, 0.889, 0.942, 0.943, 0.929, 0.909, 0.883, 0.857, 0.828, 0.790, 0.751, 0.706, 0.650, 0.591, 0.519, 0.395, 0.327, 0.245, 0.152},  {0.452, 0.445, 0.434, 0.424, 0.412, 0.404, 0.394, 0.384, 0.379, 0.366, 0.345, 0.330, 0.318, 0.315, 0.751, 0.792, 0.795, 0.786, 0.773, 0.755, 0.737, 0.717, 0.691, 0.663, 0.631, 0.591, 0.548, 0.494, 0.391, 0.342, 0.280, 0.207},  {0.417, 0.407, 0.392, 0.377, 0.360, 0.349, 0.336, 0.336, 0.315, 0.298, 0.270, 0.248, 0.231, 0.235, 1.003, 1.058, 1.048, 1.018, 0.978, 0.931, 0.882, 0.828, 0.763, 0.698, 0.627, 0.543, 0.462, 0.370, 0.233, 0.159, 0.072, 0.016},  {0.419, 0.411, 0.398, 0.386, 0.371, 0.362, 0.350, 0.338, 0.332, 0.317, 0.293, 0.274, 0.259, 0.260, 0.852, 0.896, 0.890, 0.867, 0.838, 0.803, 0.766, 0.726, 0.677, 0.627, 0.573, 0.508, 0.443, 0.369, 0.250, 0.188, 0.107, 0.043},  {0.451, 0.445, 0.434, 0.424, 0.412, 0.404, 0.394, 0.383, 0.379, 0.365, 0.345, 0.329, 0.317, 0.314, 0.726, 0.762, 0.759, 0.745, 0.725, 0.701, 0.676, 0.649, 0.614, 0.579, 0.540, 0.492, 0.444, 0.388, 0.285, 0.238, 0.171, 0.109}},
      {{0.0094, 0.0178, 0.0431, 0.0596, 0.0885, 0.123, 0.1633, 0.2063, 0.2314, 0.2329, 0.3303, 0.4354, 0.5325, 0.6321, 0.7326, 0.8357, 0.9354, 1.0354, 1.1362, 1.2351, 1.336, 1.4396, 1.5329, 1.6374, 1.7418, 1.8345, 1.9551, 2.0503, 2.1428, 2.2423, 2.3321},  {0.423, 0.414, 0.393, 0.382, 0.366, 0.351, 0.336, 0.323, 0.313, 0.297, 0.271, 0.251, 0.236, 0.228, 1.042, 1.131, 1.134, 1.115, 1.086, 1.053, 1.013, 0.969, 0.924, 0.867, 0.802, 0.735, 0.632, 0.481, 0.405, 0.309, 0.206},  {0.425, 0.417, 0.399, 0.390, 0.376, 0.363, 0.351, 0.339, 0.331, 0.316, 0.294, 0.276, 0.263, 0.255, 0.879, 0.951, 0.953, 0.939, 0.918, 0.894, 0.864, 0.831, 0.797, 0.754, 0.705, 0.654, 0.574, 0.451, 0.391, 0.314, 0.228},  {0.455, 0.449, 0.434, 0.426, 0.415, 0.404, 0.394, 0.384, 0.377, 0.364, 0.346, 0.331, 0.320, 0.312, 0.744, 0.799, 0.803, 0.794, 0.779, 0.763, 0.743, 0.720, 0.697, 0.667, 0.632, 0.595, 0.537, 0.437, 0.393, 0.336, 0.271},  {0.423, 0.414, 0.393, 0.382, 0.366, 0.350, 0.336, 0.322, 0.313, 0.296, 0.270, 0.250, 0.235, 0.226, 0.988, 1.066, 1.057, 1.026, 0.985, 0.938, 0.886, 0.826, 0.768, 0.697, 0.618, 0.543, 0.435, 0.300, 0.228, 0.144, 0.061},  {0.424, 0.417, 0.399, 0.390, 0.376, 0.363, 0.350, 0.339, 0.330, 0.316, 0.293, 0.275, 0.262, 0.253, 0.840, 0.903, 0.897, 0.874, 0.844, 0.809, 0.770, 0.725, 0.682, 0.627, 0.567, 0.508, 0.422, 0.308, 0.249, 0.177, 0.098},  {0.455, 0.449, 0.434, 0.426, 0.415, 0.404, 0.394, 0.384, 0.377, 0.364, 0.345, 0.330, 0.319, 0.311, 0.718, 0.767, 0.765, 0.751, 0.730, 0.707, 0.680, 0.649, 0.618, 0.580, 0.537, 0.494, 0.430, 0.333, 0.289, 0.234, 0.169}},
      {{0.0188, 0.0303, 0.0557, 0.0841, 0.1175, 0.1574, 0.2007, 0.252, 0.2781, 0.338, 0.4376, 0.5394, 0.6339, 0.7372, 0.8377, 0.9379, 1.0397, 1.1356, 1.2395, 1.3355, 1.4396, 1.5353, 1.6336, 1.7332, 1.833, 1.867, 1.9053, 2.1012, 2.2144, 2.3178, 2.4107},  {0.430, 0.420, 0.402, 0.385, 0.368, 0.351, 0.337, 0.321, 0.310, 0.281, 0.258, 0.241, 0.228, 0.220, 1.032, 1.112, 1.125, 1.112, 1.085, 1.053, 1.002, 0.959, 0.909, 0.850, 0.782, 0.756, 0.619, 0.499, 0.388, 0.272, 0.145},  {0.430, 0.422, 0.406, 0.392, 0.377, 0.363, 0.351, 0.337, 0.328, 0.302, 0.282, 0.267, 0.255, 0.247, 0.873, 0.936, 0.947, 0.938, 0.917, 0.894, 0.856, 0.824, 0.787, 0.743, 0.692, 0.672, 0.560, 0.468, 0.379, 0.285, 0.179},  {0.459, 0.452, 0.439, 0.427, 0.415, 0.404, 0.394, 0.382, 0.375, 0.352, 0.336, 0.323, 0.313, 0.305, 0.740, 0.790, 0.799, 0.794, 0.781, 0.765, 0.739, 0.717, 0.691, 0.660, 0.623, 0.609, 0.519, 0.452, 0.387, 0.317, 0.238},  {0.430, 0.420, 0.402, 0.384, 0.368, 0.351, 0.336, 0.321, 0.310, 0.281, 0.258, 0.240, 0.226, 0.218, 0.972, 1.039, 1.038, 1.013, 0.972, 0.926, 0.862, 0.804, 0.740, 0.668, 0.589, 0.561, 0.442, 0.307, 0.207, 0.109, 0.019},  {0.430, 0.421, 0.406, 0.391, 0.377, 0.363, 0.351, 0.337, 0.328, 0.302, 0.282, 0.266, 0.254, 0.245, 0.828, 0.883, 0.884, 0.865, 0.835, 0.801, 0.753, 0.710, 0.661, 0.606, 0.546, 0.523, 0.424, 0.316, 0.233, 0.146, 0.049},  {0.459, 0.452, 0.439, 0.427, 0.415, 0.404, 0.393, 0.382, 0.374, 0.352, 0.335, 0.322, 0.312, 0.303, 0.710, 0.754, 0.757, 0.746, 0.726, 0.703, 0.669, 0.640, 0.606, 0.567, 0.523, 0.507, 0.423, 0.343, 0.281, 0.212, 0.125}},
      {{0.032, 0.0492, 0.0786, 0.1129, 0.1433, 0.1972, 0.2489, 0.3026, 0.3286, 0.3311, 0.428, 0.5312, 0.6328, 0.7301, 0.8284, 0.9305, 1.0295, 1.1287, 1.2304, 1.333, 1.4264, 1.4307, 1.5182, 1.5201, 1.6233, 1.7255, 1.8271, 1.9251, 2.034, 2.1287, 2.2229, 2.3246, 2.4363},  {0.438, 0.427, 0.408, 0.388, 0.374, 0.352, 0.335, 0.319, 0.310, 0.293, 0.271, 0.253, 0.239, 0.228, 0.224, 1.066, 1.127, 1.131, 1.116, 1.088, 1.055, 0.779, 0.787, 0.792, 0.779, 0.747, 0.705, 0.657, 0.595, 0.523, 0.422, 0.300, 0.135},  {0.437, 0.427, 0.411, 0.394, 0.382, 0.364, 0.349, 0.335, 0.327, 0.312, 0.293, 0.278, 0.265, 0.255, 0.250, 0.899, 0.948, 0.951, 0.941, 0.920, 0.896, 0.679, 0.687, 0.691, 0.682, 0.660, 0.628, 0.591, 0.543, 0.488, 0.408, 0.309, 0.172},  {0.464, 0.456, 0.443, 0.429, 0.419, 0.404, 0.392, 0.381, 0.374, 0.360, 0.344, 0.331, 0.321, 0.312, 0.306, 0.760, 0.798, 0.803, 0.796, 0.783, 0.766, 0.599, 0.607, 0.610, 0.606, 0.592, 0.570, 0.545, 0.511, 0.471, 0.412, 0.338, 0.236},  {0.438, 0.427, 0.408, 0.388, 0.374, 0.352, 0.335, 0.319, 0.310, 0.293, 0.271, 0.253, 0.238, 0.226, 0.221, 0.995, 1.042, 1.032, 1.002, 0.959, 0.912, 0.673, 0.666, 0.669, 0.640, 0.594, 0.539, 0.479, 0.405, 0.329, 0.238, 0.134, 0.017},  {0.436, 0.427, 0.411, 0.394, 0.382, 0.364, 0.349, 0.335, 0.327, 0.312, 0.293, 0.277, 0.264, 0.253, 0.247, 0.847, 0.885, 0.879, 0.858, 0.826, 0.791, 0.599, 0.596, 0.599, 0.578, 0.543, 0.501, 0.455, 0.397, 0.336, 0.261, 0.170, 0.047},  {0.464, 0.456, 0.443, 0.429, 0.419, 0.404, 0.392, 0.381, 0.374, 0.360, 0.344, 0.331, 0.320, 0.311, 0.305, 0.725, 0.756, 0.755, 0.741, 0.720, 0.696, 0.544, 0.545, 0.547, 0.534, 0.512, 0.482, 0.450, 0.407, 0.362, 0.305, 0.235, 0.127}},
      {{0.0484, 0.0698, 0.1022, 0.137, 0.1815, 0.2298, 0.2816, 0.3355, 0.3668, 0.4083, 0.5092, 0.6092, 0.7061, 0.8119, 0.9089, 0.9605, 1.0089, 1.1076, 1.206, 1.3058, 1.3497, 1.4207, 1.5092, 1.5443, 1.6049, 1.7053, 1.8076, 1.9069, 2.0042, 2.1102, 2.2111, 2.302, 2.4087},  {0.446, 0.433, 0.413, 0.394, 0.373, 0.355, 0.338, 0.322, 0.311, 0.285, 0.265, 0.249, 0.237, 0.226, 0.219, 0.754, 1.014, 1.101, 1.112, 1.104, 1.095, 0.990, 0.820, 0.820, 0.821, 0.797, 0.760, 0.715, 0.664, 0.597, 0.506, 0.395, 0.236},  {0.442, 0.432, 0.415, 0.399, 0.382, 0.366, 0.351, 0.337, 0.328, 0.305, 0.288, 0.274, 0.263, 0.252, 0.246, 0.652, 0.859, 0.928, 0.938, 0.932, 0.926, 0.845, 0.712, 0.713, 0.714, 0.698, 0.671, 0.637, 0.598, 0.546, 0.476, 0.387, 0.256},  {0.468, 0.460, 0.446, 0.433, 0.418, 0.405, 0.393, 0.382, 0.374, 0.354, 0.339, 0.328, 0.318, 0.309, 0.303, 0.567, 0.730, 0.784, 0.793, 0.791, 0.787, 0.727, 0.625, 0.627, 0.630, 0.620, 0.602, 0.578, 0.551, 0.515, 0.463, 0.397, 0.298},  {0.445, 0.433, 0.413, 0.394, 0.373, 0.355, 0.337, 0.322, 0.310, 0.285, 0.265, 0.249, 0.236, 0.223, 0.216, 0.696, 0.941, 1.011, 1.008, 0.984, 0.969, 0.865, 0.699, 0.693, 0.683, 0.643, 0.591, 0.532, 0.469, 0.390, 0.299, 0.201, 0.076},  {0.442, 0.432, 0.415, 0.399, 0.382, 0.366, 0.351, 0.337, 0.328, 0.305, 0.287, 0.273, 0.262, 0.251, 0.243, 0.609, 0.805, 0.862, 0.861, 0.844, 0.833, 0.753, 0.622, 0.618, 0.611, 0.582, 0.543, 0.498, 0.448, 0.386, 0.312, 0.230, 0.117},  {0.468, 0.460, 0.446, 0.433, 0.418, 0.405, 0.393, 0.382, 0.373, 0.353, 0.339, 0.327, 0.318, 0.308, 0.301, 0.536, 0.693, 0.740, 0.742, 0.732, 0.725, 0.665, 0.563, 0.562, 0.559, 0.540, 0.514, 0.482, 0.447, 0.402, 0.346, 0.284, 0.197}},
      {{0.0987, 0.1221, 0.161, 0.199, 0.2447, 0.3006, 0.3524, 0.4134, 0.4442, 0.4813, 0.579, 0.6841, 0.7724, 0.8804, 0.9697, 1.0836, 1.1571, 1.2582, 1.3564, 1.4338, 1.5189, 1.6085, 1.6216, 1.7217, 1.8201, 1.9193, 2.0202, 2.1313, 2.211, 2.3205, 2.3502, 2.4163, 2.527},  {0.460, 0.446, 0.422, 0.402, 0.381, 0.360, 0.343, 0.325, 0.314, 0.290, 0.270, 0.254, 0.242, 0.230, 0.222, 0.215, 0.634, 1.007, 1.046, 1.050, 0.838, 0.853, 0.858, 0.842, 0.809, 0.766, 0.715, 0.647, 0.580, 0.455, 0.417, 0.312, 0.097},  {0.454, 0.442, 0.422, 0.405, 0.388, 0.369, 0.355, 0.339, 0.330, 0.309, 0.292, 0.277, 0.267, 0.256, 0.249, 0.241, 0.556, 0.856, 0.887, 0.892, 0.726, 0.739, 0.743, 0.733, 0.709, 0.677, 0.638, 0.586, 0.535, 0.436, 0.406, 0.320, 0.132},  {0.477, 0.467, 0.451, 0.437, 0.422, 0.407, 0.395, 0.382, 0.375, 0.355, 0.342, 0.329, 0.321, 0.311, 0.305, 0.297, 0.491, 0.731, 0.757, 0.763, 0.635, 0.648, 0.652, 0.646, 0.631, 0.609, 0.583, 0.546, 0.509, 0.436, 0.413, 0.348, 0.185},  {0.460, 0.446, 0.422, 0.402, 0.381, 0.360, 0.342, 0.324, 0.314, 0.289, 0.270, 0.254, 0.242, 0.229, 0.220, 0.211, 0.575, 0.917, 0.939, 0.931, 0.726, 0.725, 0.727, 0.695, 0.647, 0.590, 0.525, 0.443, 0.373, 0.258, 0.222, 0.131, 0.131},  {0.454, 0.442, 0.422, 0.405, 0.387, 0.369, 0.355, 0.339, 0.330, 0.308, 0.292, 0.277, 0.267, 0.255, 0.247, 0.238, 0.511, 0.789, 0.808, 0.804, 0.643, 0.644, 0.646, 0.623, 0.587, 0.544, 0.494, 0.430, 0.374, 0.279, 0.249, 0.168, 0.168},  {0.477, 0.467, 0.451, 0.437, 0.422, 0.407, 0.395, 0.382, 0.374, 0.355, 0.341, 0.329, 0.321, 0.311, 0.303, 0.295, 0.458, 0.686, 0.704, 0.704, 0.578, 0.583, 0.585, 0.571, 0.548, 0.518, 0.483, 0.437, 0.395, 0.324, 0.300, 0.239, 0.239}},
      {{0.1633, 0.1852, 0.2184, 0.2609, 0.313, 0.3689, 0.4354, 0.4808, 0.5168, 0.5448, 0.6601, 0.7591, 0.8595, 0.9648, 1.0669, 1.1711, 1.2689, 1.3613, 1.4624, 1.5622, 1.6186, 1.6612, 1.7614, 1.7711, 1.8647, 1.9639, 2.0608, 2.1636, 2.2679, 2.3693, 2.46, 2.4818},  {0.476, 0.461, 0.440, 0.417, 0.392, 0.370, 0.348, 0.335, 0.322, 0.316, 0.281, 0.263, 0.249, 0.236, 0.226, 0.218, 0.211, 0.218, 0.910, 0.976, 0.984, 0.742, 0.792, 0.791, 0.789, 0.756, 0.709, 0.641, 0.553, 0.454, 0.298, 0.248},  {0.466, 0.454, 0.437, 0.417, 0.396, 0.378, 0.359, 0.347, 0.337, 0.330, 0.301, 0.285, 0.272, 0.261, 0.252, 0.244, 0.238, 0.243, 0.781, 0.835, 0.842, 0.651, 0.693, 0.693, 0.692, 0.668, 0.633, 0.582, 0.513, 0.436, 0.308, 0.267},  {0.486, 0.476, 0.462, 0.445, 0.428, 0.413, 0.397, 0.387, 0.378, 0.372, 0.348, 0.334, 0.324, 0.314, 0.306, 0.299, 0.293, 0.295, 0.676, 0.720, 0.728, 0.578, 0.615, 0.615, 0.618, 0.602, 0.579, 0.542, 0.493, 0.436, 0.339, 0.309},  {0.475, 0.461, 0.440, 0.416, 0.392, 0.370, 0.348, 0.334, 0.322, 0.316, 0.281, 0.263, 0.249, 0.236, 0.225, 0.216, 0.207, 0.211, 0.813, 0.859, 0.858, 0.636, 0.664, 0.662, 0.642, 0.594, 0.534, 0.458, 0.369, 0.274, 0.145, 0.108},  {0.466, 0.454, 0.436, 0.417, 0.396, 0.377, 0.359, 0.347, 0.336, 0.330, 0.300, 0.285, 0.272, 0.261, 0.251, 0.242, 0.234, 0.238, 0.710, 0.748, 0.748, 0.571, 0.597, 0.596, 0.582, 0.546, 0.500, 0.441, 0.371, 0.293, 0.181, 0.148},  {0.486, 0.476, 0.461, 0.445, 0.428, 0.413, 0.397, 0.387, 0.378, 0.372, 0.347, 0.334, 0.324, 0.314, 0.305, 0.297, 0.291, 0.291, 0.628, 0.663, 0.665, 0.523, 0.549, 0.549, 0.542, 0.518, 0.487, 0.444, 0.392, 0.335, 0.251, 0.227}},
      {{0.2241, 0.2542, 0.2888, 0.338, 0.3907, 0.4184, 0.4912, 0.5441, 0.5804, 0.6431, 0.739, 0.841, 0.9404, 1.0414, 1.1442, 1.2482, 1.3489, 1.4494, 1.5385, 1.5651, 1.6428, 1.7331, 1.7358, 1.8352, 1.925, 1.9373, 2.0361, 2.1339, 2.2353, 2.3371, 2.4347, 2.5323, 2.5624},  {0.490, 0.469, 0.447, 0.421, 0.396, 0.385, 0.360, 0.343, 0.331, 0.301, 0.283, 0.267, 0.255, 0.243, 0.233, 0.224, 0.217, 0.210, 0.211, 0.480, 0.853, 0.917, 0.131, 0.632, 0.726, 0.737, 0.727, 0.685, 0.621, 0.542, 0.455, 0.247, 0.176},  {0.477, 0.460, 0.442, 0.419, 0.399, 0.390, 0.368, 0.354, 0.343, 0.317, 0.301, 0.288, 0.277, 0.267, 0.257, 0.250, 0.243, 0.237, 0.238, 0.432, 0.737, 0.790, 0.162, 0.563, 0.641, 0.651, 0.645, 0.615, 0.566, 0.504, 0.437, 0.266, 0.205},  {0.494, 0.480, 0.465, 0.446, 0.429, 0.422, 0.403, 0.392, 0.383, 0.360, 0.347, 0.336, 0.327, 0.318, 0.310, 0.303, 0.297, 0.292, 0.291, 0.391, 0.644, 0.688, 0.221, 0.509, 0.577, 0.584, 0.584, 0.564, 0.530, 0.486, 0.438, 0.307, 0.260},  {0.489, 0.469, 0.447, 0.420, 0.396, 0.385, 0.359, 0.343, 0.330, 0.301, 0.283, 0.267, 0.254, 0.243, 0.233, 0.223, 0.214, 0.205, 0.202, 0.413, 0.747, 0.791, 0.123, 0.525, 0.591, 0.598, 0.570, 0.516, 0.444, 0.362, 0.276, 0.115, 0.014},  {0.477, 0.460, 0.442, 0.419, 0.399, 0.390, 0.368, 0.354, 0.343, 0.317, 0.301, 0.288, 0.277, 0.266, 0.257, 0.248, 0.241, 0.233, 0.230, 0.378, 0.658, 0.696, 0.155, 0.481, 0.539, 0.546, 0.526, 0.485, 0.430, 0.364, 0.295, 0.154, 0.045},  {0.494, 0.480, 0.465, 0.446, 0.429, 0.421, 0.403, 0.392, 0.383, 0.360, 0.347, 0.336, 0.326, 0.318, 0.310, 0.302, 0.295, 0.289, 0.286, 0.356, 0.591, 0.625, 0.216, 0.451, 0.506, 0.512, 0.502, 0.474, 0.435, 0.387, 0.336, 0.231, 0.144}},
      {{0.2862, 0.3086, 0.3534, 0.3958, 0.4409, 0.4959, 0.5541, 0.6068, 0.6442, 0.6483, 0.7405, 0.8402, 0.9436, 1.0401, 1.1427, 1.2402, 1.3473, 1.4372, 1.6098, 1.6389, 1.7395, 1.8416, 1.8981, 1.9484, 2.0229, 2.0483, 2.1425, 2.1914, 2.1935, 2.2423, 2.3402, 2.4403, 2.5407, 2.5649, 2.5650},  {0.505, 0.488, 0.459, 0.435, 0.413, 0.390, 0.369, 0.352, 0.338, 0.326, 0.302, 0.285, 0.270, 0.258, 0.246, 0.237, 0.227, 0.220, 0.206, 0.204, 0.207, 0.870, 0.912, 0.125, 0.671, 0.776, 0.666, 0.653, 0.660, 0.639, 0.581, 0.496, 0.357, 0.296, 0.296},  {0.489, 0.476, 0.451, 0.431, 0.412, 0.393, 0.376, 0.361, 0.349, 0.338, 0.317, 0.303, 0.290, 0.279, 0.269, 0.260, 0.252, 0.245, 0.234, 0.232, 0.234, 0.751, 0.785, 0.156, 0.594, 0.680, 0.597, 0.588, 0.593, 0.579, 0.534, 0.468, 0.356, 0.306, 0.306},  {0.503, 0.491, 0.471, 0.455, 0.439, 0.423, 0.408, 0.397, 0.386, 0.376, 0.359, 0.347, 0.336, 0.327, 0.319, 0.311, 0.304, 0.298, 0.289, 0.287, 0.289, 0.654, 0.683, 0.215, 0.534, 0.603, 0.548, 0.543, 0.547, 0.537, 0.506, 0.459, 0.374, 0.337, 0.337},  {0.504, 0.488, 0.458, 0.434, 0.412, 0.390, 0.369, 0.352, 0.338, 0.326, 0.302, 0.284, 0.270, 0.257, 0.246, 0.236, 0.226, 0.218, 0.201, 0.197, 0.196, 0.723, 0.747, 0.112, 0.518, 0.600, 0.517, 0.495, 0.499, 0.473, 0.405, 0.320, 0.205, 0.162, 0.162},  {0.489, 0.475, 0.451, 0.431, 0.412, 0.393, 0.375, 0.361, 0.349, 0.338, 0.317, 0.302, 0.289, 0.279, 0.268, 0.260, 0.251, 0.243, 0.229, 0.227, 0.225, 0.641, 0.663, 0.145, 0.477, 0.547, 0.484, 0.468, 0.471, 0.451, 0.398, 0.330, 0.234, 0.196, 0.196},  {0.502, 0.491, 0.471, 0.454, 0.439, 0.423, 0.408, 0.397, 0.386, 0.376, 0.359, 0.347, 0.336, 0.327, 0.319, 0.311, 0.303, 0.297, 0.285, 0.283, 0.282, 0.579, 0.600, 0.206, 0.451, 0.512, 0.469, 0.459, 0.461, 0.448, 0.411, 0.360, 0.289, 0.261, 0.261}},
      {{0.3507, 0.3756, 0.4147, 0.4564, 0.5074, 0.5621, 0.6203, 0.6728, 0.7091, 0.7744, 0.8732, 0.9762, 1.0737, 1.171, 1.2712, 1.3717, 1.4739, 1.5714, 1.6745, 1.7765, 1.8757, 1.9368, 1.9694, 2.0709, 2.1328, 2.2321, 2.3324, 2.4322, 2.5152},  {0.521, 0.502, 0.475, 0.450, 0.424, 0.401, 0.379, 0.362, 0.348, 0.319, 0.300, 0.284, 0.271, 0.259, 0.248, 0.238, 0.229, 0.220, 0.211, 0.202, 0.194, 0.459, 0.699, 0.812, 0.254, 0.623, 0.613, 0.553, 0.473},  {0.502, 0.486, 0.464, 0.443, 0.421, 0.402, 0.383, 0.369, 0.357, 0.332, 0.315, 0.301, 0.290, 0.280, 0.270, 0.262, 0.253, 0.246, 0.238, 0.231, 0.224, 0.416, 0.617, 0.710, 0.242, 0.563, 0.558, 0.512, 0.449},  {0.512, 0.499, 0.481, 0.463, 0.445, 0.429, 0.414, 0.402, 0.391, 0.370, 0.356, 0.345, 0.336, 0.327, 0.319, 0.312, 0.305, 0.298, 0.292, 0.286, 0.281, 0.382, 0.551, 0.629, 0.238, 0.520, 0.521, 0.490, 0.444},  {0.521, 0.502, 0.475, 0.450, 0.424, 0.400, 0.379, 0.362, 0.348, 0.319, 0.299, 0.284, 0.271, 0.259, 0.248, 0.238, 0.228, 0.219, 0.208, 0.196, 0.183, 0.369, 0.576, 0.654, 0.184, 0.477, 0.448, 0.379, 0.301},  {0.502, 0.486, 0.464, 0.443, 0.421, 0.401, 0.383, 0.369, 0.356, 0.332, 0.315, 0.301, 0.290, 0.280, 0.270, 0.261, 0.253, 0.244, 0.235, 0.225, 0.214, 0.344, 0.522, 0.591, 0.183, 0.450, 0.431, 0.377, 0.313},  {0.512, 0.499, 0.481, 0.463, 0.445, 0.429, 0.414, 0.402, 0.391, 0.370, 0.356, 0.345, 0.336, 0.327, 0.319, 0.312, 0.304, 0.297, 0.290, 0.282, 0.274, 0.332, 0.484, 0.548, 0.227, 0.441, 0.432, 0.394, 0.347}},
      {{0.4172, 0.4384, 0.4817, 0.5229, 0.5735, 0.6281, 0.6865, 0.7473, 0.7819, 0.7858, 0.8796, 0.9794, 1.0808, 1.1778, 1.28, 1.3839, 1.4769, 1.5757, 1.6749, 1.8857, 1.9817, 2.0831, 2.1795, 2.252, 2.2562, 2.3548, 2.4287, 2.4515, 2.5519, 2.5719, 2.5829},  {0.540, 0.523, 0.491, 0.465, 0.438, 0.413, 0.390, 0.370, 0.356, 0.345, 0.321, 0.302, 0.287, 0.274, 0.262, 0.250, 0.241, 0.231, 0.221, 0.200, 0.190, 0.182, 0.688, 0.734, 0.118, 0.560, 0.563, 0.562, 0.479, 0.456, 0.444},  {0.517, 0.503, 0.476, 0.455, 0.432, 0.411, 0.392, 0.375, 0.363, 0.353, 0.333, 0.317, 0.304, 0.293, 0.282, 0.272, 0.263, 0.255, 0.246, 0.229, 0.221, 0.213, 0.610, 0.650, 0.151, 0.512, 0.518, 0.518, 0.454, 0.436, 0.426},  {0.523, 0.511, 0.489, 0.472, 0.453, 0.436, 0.420, 0.406, 0.396, 0.386, 0.370, 0.357, 0.346, 0.337, 0.328, 0.320, 0.312, 0.305, 0.298, 0.284, 0.278, 0.273, 0.551, 0.586, 0.211, 0.480, 0.489, 0.491, 0.446, 0.432, 0.425},  {0.540, 0.522, 0.491, 0.465, 0.438, 0.413, 0.390, 0.370, 0.356, 0.345, 0.321, 0.302, 0.287, 0.274, 0.262, 0.250, 0.240, 0.230, 0.220, 0.195, 0.181, 0.166, 0.543, 0.568, 0.100, 0.414, 0.403, 0.397, 0.315, 0.295, 0.283},  {0.517, 0.503, 0.476, 0.455, 0.432, 0.411, 0.392, 0.375, 0.363, 0.353, 0.333, 0.317, 0.303, 0.292, 0.282, 0.272, 0.263, 0.254, 0.245, 0.225, 0.213, 0.200, 0.500, 0.524, 0.135, 0.399, 0.394, 0.390, 0.324, 0.308, 0.299},  {0.523, 0.511, 0.489, 0.472, 0.453, 0.436, 0.420, 0.406, 0.396, 0.386, 0.369, 0.356, 0.346, 0.337, 0.328, 0.319, 0.312, 0.305, 0.297, 0.281, 0.273, 0.262, 0.473, 0.499, 0.199, 0.400, 0.402, 0.401, 0.354, 0.341, 0.335}},
      {{0.4786, 0.5039, 0.5477, 0.583, 0.6468, 0.6954, 0.7625, 0.8173, 0.8565, 0.9311, 1.0241, 1.1271, 1.2236, 1.3247, 1.428, 1.5207, 1.6215, 1.7214, 1.8213, 1.927, 2.0055, 2.1314, 2.233, 2.3267, 2.4084, 2.4913, 2.5922, 2.63, 2.6607, 2.6765},  {0.561, 0.539, 0.505, 0.481, 0.445, 0.422, 0.396, 0.376, 0.360, 0.332, 0.313, 0.296, 0.282, 0.269, 0.257, 0.247, 0.236, 0.224, 0.212, 0.199, 0.189, 0.174, 0.163, 0.539, 0.708, 0.544, 0.456, 0.420, 0.386, 0.369},  {0.533, 0.515, 0.487, 0.468, 0.438, 0.419, 0.396, 0.380, 0.366, 0.342, 0.325, 0.311, 0.299, 0.288, 0.278, 0.268, 0.259, 0.249, 0.238, 0.227, 0.219, 0.206, 0.196, 0.487, 0.628, 0.497, 0.434, 0.406, 0.379, 0.365},  {0.535, 0.520, 0.497, 0.481, 0.457, 0.441, 0.422, 0.408, 0.397, 0.375, 0.362, 0.350, 0.341, 0.332, 0.323, 0.315, 0.307, 0.298, 0.290, 0.281, 0.274, 0.264, 0.257, 0.447, 0.568, 0.463, 0.429, 0.409, 0.390, 0.380},  {0.560, 0.538, 0.504, 0.481, 0.445, 0.422, 0.395, 0.376, 0.359, 0.332, 0.312, 0.295, 0.282, 0.269, 0.257, 0.246, 0.235, 0.223, 0.210, 0.196, 0.184, 0.165, 0.148, 0.388, 0.505, 0.358, 0.309, 0.278, 0.248, 0.232},  {0.533, 0.515, 0.487, 0.468, 0.438, 0.419, 0.396, 0.380, 0.365, 0.341, 0.325, 0.311, 0.299, 0.288, 0.277, 0.268, 0.258, 0.248, 0.237, 0.225, 0.215, 0.198, 0.183, 0.369, 0.473, 0.351, 0.318, 0.294, 0.269, 0.256},  {0.534, 0.520, 0.497, 0.481, 0.457, 0.441, 0.422, 0.408, 0.396, 0.375, 0.362, 0.350, 0.340, 0.331, 0.323, 0.315, 0.306, 0.298, 0.289, 0.279, 0.271, 0.258, 0.247, 0.360, 0.460, 0.358, 0.347, 0.329, 0.311, 0.302}},
      {{0.5496, 0.573, 0.6136, 0.6583, 0.7076, 0.7669, 0.8299, 0.8895, 0.9349, 0.981, 1.0742, 1.1774, 1.2751, 1.3772, 1.475, 1.5799, 1.6746, 1.7733, 1.7813, 1.8107, 1.9133, 2.0141, 2.1167, 2.225, 2.3157, 2.4002, 2.5056, 2.6112, 2.6608, 2.7104, 2.7218, 2.7496, 2.7827},  {0.589, 0.567, 0.532, 0.500, 0.469, 0.438, 0.411, 0.388, 0.368, 0.347, 0.326, 0.308, 0.293, 0.279, 0.264, 0.246, 0.228, 0.208, 0.206, 0.171, 0.158, 0.146, 0.135, 0.124, 0.115, 0.108, 0.099, 0.252, 0.369, 0.353, 0.343, 0.320, 0.282},  {0.556, 0.538, 0.509, 0.483, 0.457, 0.432, 0.409, 0.390, 0.372, 0.354, 0.336, 0.321, 0.308, 0.296, 0.284, 0.268, 0.252, 0.234, 0.231, 0.199, 0.187, 0.176, 0.166, 0.155, 0.148, 0.140, 0.132, 0.253, 0.361, 0.350, 0.342, 0.325, 0.294},  {0.552, 0.537, 0.514, 0.492, 0.471, 0.450, 0.431, 0.415, 0.400, 0.383, 0.369, 0.356, 0.346, 0.336, 0.326, 0.313, 0.299, 0.284, 0.282, 0.251, 0.240, 0.231, 0.222, 0.213, 0.206, 0.200, 0.193, 0.264, 0.367, 0.364, 0.359, 0.348, 0.326},  {0.589, 0.567, 0.532, 0.500, 0.469, 0.438, 0.411, 0.388, 0.368, 0.347, 0.326, 0.307, 0.292, 0.278, 0.264, 0.245, 0.227, 0.207, 0.204, 0.170, 0.157, 0.145, 0.133, 0.120, 0.108, 0.097, 0.082, 0.168, 0.250, 0.240, 0.231, 0.210, 0.180},  {0.556, 0.538, 0.509, 0.483, 0.457, 0.432, 0.409, 0.389, 0.372, 0.354, 0.336, 0.320, 0.308, 0.296, 0.283, 0.267, 0.251, 0.233, 0.230, 0.198, 0.186, 0.175, 0.164, 0.151, 0.141, 0.130, 0.116, 0.183, 0.265, 0.259, 0.253, 0.235, 0.210},  {0.552, 0.537, 0.514, 0.492, 0.471, 0.450, 0.431, 0.415, 0.400, 0.383, 0.368, 0.356, 0.346, 0.336, 0.326, 0.312, 0.299, 0.283, 0.281, 0.250, 0.240, 0.230, 0.221, 0.210, 0.201, 0.192, 0.180, 0.209, 0.297, 0.299, 0.294, 0.284, 0.265}},
      {{0.6043, 0.6267, 0.6689, 0.7159, 0.7652, 0.8228, 0.8934, 0.9629, 1.008, 1.0109, 1.1111, 1.213, 1.3147, 1.413, 1.5099, 1.6074, 1.7099, 1.8104, 1.8605, 1.9101, 2.011, 2.1119, 2.2151, 2.3229, 2.4121, 2.5158, 2.6128, 2.7109, 2.7607, 2.8101, 2.8489, 2.8561, 2.8681},  {0.618, 0.595, 0.556, 0.518, 0.484, 0.451, 0.418, 0.389, 0.368, 0.356, 0.332, 0.312, 0.293, 0.276, 0.257, 0.237, 0.214, 0.190, 0.175, 0.163, 0.148, 0.135, 0.123, 0.112, 0.103, 0.094, 0.086, 0.084, 0.226, 0.249, 0.209, 0.209, 0.201},  {0.579, 0.560, 0.528, 0.497, 0.469, 0.442, 0.413, 0.390, 0.371, 0.360, 0.340, 0.323, 0.307, 0.293, 0.276, 0.258, 0.238, 0.216, 0.202, 0.191, 0.177, 0.165, 0.154, 0.143, 0.135, 0.126, 0.118, 0.116, 0.238, 0.264, 0.233, 0.233, 0.226},  {0.569, 0.553, 0.527, 0.502, 0.479, 0.457, 0.433, 0.413, 0.397, 0.386, 0.369, 0.355, 0.342, 0.331, 0.317, 0.301, 0.284, 0.265, 0.252, 0.242, 0.230, 0.220, 0.210, 0.201, 0.193, 0.185, 0.178, 0.175, 0.263, 0.297, 0.279, 0.279, 0.275},  {0.618, 0.595, 0.556, 0.518, 0.484, 0.451, 0.418, 0.389, 0.367, 0.356, 0.332, 0.311, 0.293, 0.276, 0.256, 0.236, 0.213, 0.189, 0.173, 0.161, 0.146, 0.133, 0.121, 0.107, 0.095, 0.082, 0.069, 0.060, 0.144, 0.165, 0.140, 0.135, 0.129},  {0.579, 0.560, 0.528, 0.497, 0.469, 0.442, 0.413, 0.390, 0.371, 0.360, 0.340, 0.322, 0.307, 0.292, 0.275, 0.257, 0.237, 0.215, 0.201, 0.190, 0.176, 0.163, 0.152, 0.139, 0.128, 0.115, 0.102, 0.092, 0.169, 0.194, 0.174, 0.170, 0.165},  {0.569, 0.553, 0.527, 0.502, 0.479, 0.457, 0.433, 0.413, 0.397, 0.386, 0.369, 0.355, 0.342, 0.330, 0.316, 0.301, 0.283, 0.264, 0.251, 0.241, 0.229, 0.218, 0.208, 0.197, 0.188, 0.177, 0.165, 0.156, 0.209, 0.246, 0.235, 0.232, 0.229}},
      {{0.6619, 0.6829, 0.7264, 0.7703, 0.8246, 0.8804, 0.9533, 1.0364, 1.0886, 1.1214, 1.2227, 1.3212, 1.4203, 1.5191, 1.6227, 1.7289, 1.8236, 1.9235, 2.0233, 2.1219, 2.2232, 2.3241, 2.4243, 2.5217, 2.624, 2.7242, 2.8234, 2.873, 2.9294},  {0.658, 0.634, 0.587, 0.548, 0.506, 0.469, 0.430, 0.392, 0.366, 0.344, 0.319, 0.299, 0.280, 0.260, 0.237, 0.213, 0.192, 0.173, 0.155, 0.139, 0.124, 0.110, 0.099, 0.089, 0.080, 0.072, 0.068, 0.145, 0.169},  {0.610, 0.591, 0.553, 0.521, 0.486, 0.456, 0.422, 0.390, 0.368, 0.348, 0.328, 0.310, 0.294, 0.277, 0.257, 0.235, 0.217, 0.200, 0.183, 0.168, 0.154, 0.141, 0.130, 0.120, 0.111, 0.103, 0.099, 0.165, 0.198},  {0.593, 0.575, 0.545, 0.518, 0.490, 0.465, 0.436, 0.409, 0.391, 0.371, 0.355, 0.341, 0.328, 0.314, 0.297, 0.279, 0.263, 0.248, 0.234, 0.220, 0.207, 0.196, 0.187, 0.178, 0.170, 0.163, 0.157, 0.198, 0.247},  {0.658, 0.634, 0.587, 0.548, 0.506, 0.469, 0.430, 0.392, 0.366, 0.344, 0.319, 0.298, 0.279, 0.259, 0.236, 0.212, 0.191, 0.171, 0.153, 0.137, 0.121, 0.107, 0.094, 0.081, 0.068, 0.054, 0.044, 0.083, 0.107},  {0.610, 0.591, 0.553, 0.521, 0.486, 0.456, 0.422, 0.390, 0.368, 0.348, 0.327, 0.310, 0.293, 0.276, 0.256, 0.234, 0.216, 0.198, 0.181, 0.166, 0.151, 0.138, 0.125, 0.112, 0.099, 0.085, 0.073, 0.110, 0.143},  {0.593, 0.575, 0.545, 0.518, 0.490, 0.465, 0.436, 0.409, 0.391, 0.371, 0.355, 0.341, 0.328, 0.314, 0.296, 0.278, 0.262, 0.247, 0.232, 0.219, 0.206, 0.194, 0.183, 0.172, 0.161, 0.149, 0.137, 0.154, 0.205}},
      {{0.7182, 0.7405, 0.7812, 0.8362, 0.883, 0.9415, 1.0224, 1.1023, 1.1666, 1.203, 1.3002, 1.4053, 1.5046, 1.6048, 1.7044, 1.8057, 1.9078, 2.0084, 2.1051, 2.2023, 2.3, 2.4014, 2.5046, 2.6046, 2.7045, 2.8014, 2.9019, 2.9334, 2.966, 3.001, 3.0263, 3.0676, 3.0707},  {0.710, 0.680, 0.629, 0.572, 0.531, 0.487, 0.437, 0.396, 0.361, 0.338, 0.313, 0.289, 0.268, 0.248, 0.229, 0.208, 0.186, 0.166, 0.148, 0.130, 0.116, 0.103, 0.091, 0.081, 0.073, 0.065, 0.059, 0.061, 0.095, 0.110, 0.115, 0.092, 0.094},  {0.652, 0.627, 0.585, 0.539, 0.505, 0.468, 0.426, 0.392, 0.362, 0.341, 0.320, 0.299, 0.282, 0.265, 0.248, 0.229, 0.210, 0.192, 0.176, 0.159, 0.145, 0.133, 0.122, 0.112, 0.103, 0.095, 0.089, 0.090, 0.120, 0.141, 0.151, 0.134, 0.136},  {0.624, 0.604, 0.570, 0.532, 0.503, 0.473, 0.438, 0.409, 0.383, 0.361, 0.344, 0.327, 0.313, 0.300, 0.287, 0.271, 0.255, 0.239, 0.225, 0.211, 0.199, 0.189, 0.179, 0.171, 0.163, 0.156, 0.150, 0.149, 0.160, 0.190, 0.211, 0.206, 0.208},  {0.710, 0.680, 0.629, 0.572, 0.531, 0.487, 0.437, 0.396, 0.361, 0.338, 0.313, 0.289, 0.268, 0.248, 0.227, 0.206, 0.184, 0.164, 0.145, 0.127, 0.112, 0.099, 0.086, 0.073, 0.060, 0.047, 0.035, 0.036, 0.051, 0.059, 0.068, 0.057, 0.058},  {0.652, 0.627, 0.586, 0.539, 0.505, 0.468, 0.426, 0.392, 0.362, 0.341, 0.319, 0.299, 0.281, 0.264, 0.247, 0.228, 0.209, 0.190, 0.173, 0.156, 0.142, 0.130, 0.117, 0.104, 0.090, 0.077, 0.063, 0.063, 0.078, 0.092, 0.105, 0.097, 0.098},  {0.624, 0.604, 0.570, 0.531, 0.503, 0.473, 0.438, 0.409, 0.383, 0.361, 0.344, 0.327, 0.313, 0.300, 0.286, 0.270, 0.253, 0.238, 0.224, 0.209, 0.197, 0.186, 0.176, 0.165, 0.154, 0.143, 0.132, 0.131, 0.129, 0.148, 0.170, 0.172, 0.174}},
      {{0.7794, 0.8011, 0.8395, 0.8872, 0.9406, 1.0073, 1.0815, 1.1816, 1.2442, 1.2654, 1.3676, 1.4697, 1.5626, 1.6678, 1.7692, 1.8655, 1.9601, 2.0645, 2.1651, 2.2681, 2.3661, 2.465, 2.5652, 2.6639, 2.7652, 2.8685, 2.9681, 3.0195, 3.0374, 3.0663, 3.09, 3.1173},  {0.789, 0.749, 0.688, 0.625, 0.568, 0.508, 0.454, 0.396, 0.358, 0.329, 0.300, 0.274, 0.254, 0.232, 0.210, 0.190, 0.170, 0.149, 0.131, 0.115, 0.101, 0.089, 0.078, 0.069, 0.061, 0.054, 0.048, 0.046, 0.047, 0.062, 0.075, 0.095},  {0.714, 0.684, 0.633, 0.581, 0.534, 0.485, 0.440, 0.391, 0.358, 0.330, 0.305, 0.284, 0.266, 0.248, 0.229, 0.212, 0.194, 0.175, 0.158, 0.143, 0.130, 0.118, 0.108, 0.099, 0.090, 0.083, 0.076, 0.074, 0.074, 0.086, 0.105, 0.131},  {0.670, 0.650, 0.617, 0.581, 0.545, 0.508, 0.471, 0.426, 0.393, 0.347, 0.327, 0.310, 0.296, 0.282, 0.267, 0.253, 0.237, 0.222, 0.208, 0.197, 0.187, 0.178, 0.170, 0.167, 0.163, 0.150, 0.144, 0.142, 0.141, 0.139, 0.153, 0.190},  {0.789, 0.749, 0.688, 0.625, 0.568, 0.508, 0.454, 0.396, 0.358, 0.329, 0.299, 0.274, 0.253, 0.231, 0.209, 0.189, 0.168, 0.146, 0.128, 0.111, 0.097, 0.084, 0.073, 0.061, 0.049, 0.038, 0.027, 0.023, 0.023, 0.028, 0.036, 0.053},  {0.714, 0.684, 0.633, 0.581, 0.534, 0.485, 0.440, 0.391, 0.358, 0.330, 0.305, 0.284, 0.266, 0.247, 0.228, 0.210, 0.192, 0.172, 0.155, 0.139, 0.126, 0.114, 0.103, 0.091, 0.079, 0.066, 0.053, 0.048, 0.047, 0.051, 0.064, 0.088},  {0.670, 0.650, 0.617, 0.581, 0.546, 0.507, 0.467, 0.419, 0.387, 0.347, 0.327, 0.310, 0.296, 0.281, 0.267, 0.252, 0.236, 0.220, 0.206, 0.195, 0.185, 0.175, 0.167, 0.157, 0.148, 0.138, 0.127, 0.122, 0.121, 0.119, 0.115, 0.152}},
      {{0.8331, 0.8544, 0.8941, 0.9418, 0.9984, 1.0705, 1.149, 1.2455, 1.3263, 1.3405, 1.4392, 1.5377, 1.6439, 1.7362, 1.8397, 1.9419, 2.0404, 2.1364, 2.2413, 2.343, 2.4398, 2.5393, 2.6409, 2.7405, 2.84, 2.9392, 3.039, 3.1401, 3.1862, 3.2043, 3.2292},  {0.887, 0.847, 0.773, 0.697, 0.618, 0.535, 0.458, 0.394, 0.342, 0.301, 0.271, 0.246, 0.223, 0.204, 0.183, 0.162, 0.142, 0.125, 0.108, 0.094, 0.082, 0.071, 0.062, 0.055, 0.048, 0.042, 0.037, 0.036, 0.062, 0.084, 0.074},  {0.789, 0.759, 0.705, 0.649, 0.590, 0.526, 0.446, 0.399, 0.344, 0.304, 0.279, 0.258, 0.238, 0.222, 0.204, 0.185, 0.167, 0.151, 0.135, 0.122, 0.110, 0.100, 0.091, 0.083, 0.076, 0.070, 0.065, 0.060, 0.093, 0.123, 0.114},  {0.725, 0.706, 0.670, 0.632, 0.593, 0.549, 0.503, 0.453, 0.410, 0.330, 0.311, 0.295, 0.279, 0.267, 0.254, 0.239, 0.226, 0.215, 0.204, 0.194, 0.186, 0.178, 0.170, 0.162, 0.156, 0.150, 0.159, 0.139, 0.142, 0.189, 0.186},  {0.887, 0.847, 0.773, 0.697, 0.619, 0.531, 0.458, 0.394, 0.342, 0.301, 0.271, 0.246, 0.222, 0.203, 0.181, 0.160, 0.139, 0.122, 0.104, 0.089, 0.077, 0.066, 0.057, 0.047, 0.037, 0.028, 0.019, 0.014, 0.024, 0.048, 0.041},  {0.789, 0.759, 0.705, 0.649, 0.589, 0.520, 0.453, 0.390, 0.344, 0.304, 0.278, 0.257, 0.237, 0.220, 0.202, 0.183, 0.164, 0.148, 0.132, 0.118, 0.106, 0.096, 0.086, 0.076, 0.065, 0.055, 0.044, 0.034, 0.049, 0.084, 0.078},  {0.725, 0.706, 0.670, 0.632, 0.588, 0.536, 0.485, 0.436, 0.395, 0.330, 0.311, 0.295, 0.279, 0.267, 0.253, 0.238, 0.225, 0.213, 0.202, 0.191, 0.183, 0.174, 0.166, 0.156, 0.146, 0.136, 0.125, 0.113, 0.107, 0.154, 0.152}},
      {{0.887, 0.9129, 0.9459, 0.9935, 1.0491, 1.12, 1.2031, 1.3071, 1.3933, 1.4051, 1.5031, 1.5996, 1.703, 1.8052, 1.8988, 2.0097, 2.1111, 2.2099, 2.306, 2.4001, 2.5031, 2.6046, 2.7063, 2.8086, 2.9002, 3.0001, 3.1013, 3.2022, 3.2302, 3.2554, 3.2571, 3.2608},  {1.021, 0.967, 0.902, 0.816, 0.732, 0.633, 0.526, 0.418, 0.340, 0.287, 0.255, 0.230, 0.207, 0.187, 0.169, 0.149, 0.131, 0.114, 0.100, 0.087, 0.076, 0.065, 0.057, 0.050, 0.044, 0.039, 0.034, 0.031, 0.038, 0.050, 0.055, 0.075},  {0.891, 0.853, 0.804, 0.742, 0.678, 0.601, 0.517, 0.431, 0.362, 0.291, 0.264, 0.242, 0.223, 0.205, 0.190, 0.173, 0.156, 0.141, 0.127, 0.116, 0.105, 0.097, 0.089, 0.082, 0.076, 0.070, 0.065, 0.060, 0.059, 0.077, 0.084, 0.111},  {0.800, 0.775, 0.742, 0.702, 0.657, 0.600, 0.536, 0.472, 0.422, 0.329, 0.308, 0.290, 0.274, 0.260, 0.248, 0.234, 0.221, 0.210, 0.200, 0.192, 0.183, 0.175, 0.168, 0.161, 0.156, 0.151, 0.145, 0.140, 0.135, 0.129, 0.134, 0.171},  {1.021, 0.967, 0.902, 0.815, 0.716, 0.602, 0.499, 0.401, 0.339, 0.285, 0.254, 0.228, 0.205, 0.184, 0.166, 0.146, 0.127, 0.110, 0.096, 0.083, 0.070, 0.060, 0.051, 0.042, 0.033, 0.025, 0.017, 0.011, 0.013, 0.016, 0.023, 0.040},  {0.891, 0.853, 0.804, 0.739, 0.661, 0.573, 0.491, 0.412, 0.346, 0.290, 0.263, 0.241, 0.221, 0.203, 0.188, 0.170, 0.153, 0.137, 0.123, 0.112, 0.101, 0.092, 0.083, 0.074, 0.064, 0.054, 0.044, 0.034, 0.031, 0.037, 0.047, 0.072},  {0.800, 0.775, 0.742, 0.693, 0.636, 0.569, 0.508, 0.450, 0.404, 0.329, 0.307, 0.289, 0.273, 0.258, 0.247, 0.232, 0.219, 0.208, 0.197, 0.188, 0.179, 0.170, 0.163, 0.154, 0.145, 0.135, 0.124, 0.113, 0.108, 0.101, 0.107, 0.137}},
      {{0.9429, 0.9626, 1.0019, 1.0535, 1.1086, 1.1764, 1.2647, 1.3707, 1.4622, 1.4657, 1.5664, 1.66, 1.7589, 1.8626, 1.9613, 2.0658, 2.1671, 2.2651, 2.3637, 2.46, 2.5655, 2.6688, 2.76, 2.8675, 2.9798, 3.0742, 3.1577, 3.25, 3.2609, 3.293, 3.319, 3.3358},  {1.217, 1.172, 1.086, 0.989, 0.874, 0.741, 0.587, 0.454, 0.350, 0.282, 0.247, 0.221, 0.197, 0.177, 0.159, 0.142, 0.127, 0.113, 0.099, 0.087, 0.075, 0.065, 0.057, 0.050, 0.043, 0.038, 0.035, 0.032, 0.032, 0.036, 0.045, 0.047},  {1.039, 1.008, 0.944, 0.874, 0.786, 0.684, 0.563, 0.457, 0.375, 0.287, 0.256, 0.234, 0.214, 0.196, 0.180, 0.166, 0.152, 0.140, 0.128, 0.117, 0.107, 0.098, 0.091, 0.083, 0.076, 0.071, 0.066, 0.061, 0.061, 0.059, 0.071, 0.075},  {0.907, 0.886, 0.845, 0.797, 0.733, 0.657, 0.565, 0.487, 0.428, 0.344, 0.309, 0.290, 0.272, 0.257, 0.244, 0.232, 0.221, 0.211, 0.202, 0.192, 0.183, 0.186, 0.169, 0.162, 0.156, 0.151, 0.146, 0.142, 0.141, 0.139, 0.129, 0.130},  {1.217, 1.172, 1.076, 0.939, 0.803, 0.671, 0.546, 0.433, 0.336, 0.279, 0.244, 0.217, 0.194, 0.173, 0.155, 0.138, 0.122, 0.108, 0.094, 0.082, 0.069, 0.059, 0.051, 0.041, 0.031, 0.023, 0.017, 0.013, 0.012, 0.013, 0.015, 0.013},  {1.039, 1.008, 0.935, 0.831, 0.725, 0.623, 0.525, 0.436, 0.359, 0.285, 0.254, 0.231, 0.211, 0.193, 0.177, 0.162, 0.148, 0.136, 0.124, 0.113, 0.102, 0.092, 0.084, 0.074, 0.063, 0.053, 0.045, 0.037, 0.035, 0.034, 0.034, 0.034},  {0.907, 0.886, 0.832, 0.755, 0.677, 0.601, 0.528, 0.463, 0.409, 0.338, 0.308, 0.288, 0.271, 0.255, 0.242, 0.230, 0.219, 0.208, 0.198, 0.189, 0.179, 0.170, 0.163, 0.156, 0.143, 0.133, 0.132, 0.115, 0.113, 0.112, 0.102, 0.097}},
      {{1.0009, 1.0267, 1.0602, 1.1058, 1.1627, 1.2375, 1.3222, 1.4362, 1.5374, 1.5391, 1.6315, 1.7307, 1.8356, 1.9348, 2.0273, 2.1338, 2.2343, 2.3297, 2.4282, 2.5378, 2.6323, 2.7363, 2.8274, 2.9272, 3.0285, 3.1283, 3.2287, 3.33, 3.3507, 3.3707, 3.384, 3.3907},  {1.527, 1.462, 1.379, 1.272, 1.124, 0.923, 0.724, 0.514, 0.374, 0.260, 0.228, 0.201, 0.177, 0.158, 0.143, 0.127, 0.113, 0.101, 0.089, 0.077, 0.067, 0.058, 0.051, 0.045, 0.040, 0.035, 0.031, 0.028, 0.030, 0.033, 0.039, 0.045},  {1.271, 1.225, 1.165, 1.086, 0.976, 0.802, 0.616, 0.503, 0.392, 0.279, 0.250, 0.225, 0.203, 0.186, 0.172, 0.158, 0.145, 0.135, 0.124, 0.113, 0.104, 0.095, 0.088, 0.081, 0.075, 0.069, 0.064, 0.059, 0.058, 0.057, 0.063, 0.071},  {1.070, 1.041, 1.002, 0.949, 0.870, 0.760, 0.601, 0.519, 0.435, 0.360, 0.328, 0.302, 0.267, 0.251, 0.239, 0.226, 0.215, 0.206, 0.197, 0.187, 0.179, 0.171, 0.165, 0.159, 0.153, 0.148, 0.143, 0.138, 0.137, 0.136, 0.129, 0.127},  {1.524, 1.439, 1.299, 1.110, 0.912, 0.737, 0.599, 0.464, 0.353, 0.256, 0.224, 0.196, 0.173, 0.154, 0.138, 0.121, 0.108, 0.095, 0.083, 0.071, 0.061, 0.051, 0.044, 0.036, 0.028, 0.021, 0.015, 0.010, 0.010, 0.011, 0.012, 0.014},  {1.269, 1.205, 1.100, 0.957, 0.806, 0.672, 0.564, 0.458, 0.371, 0.276, 0.247, 0.222, 0.200, 0.183, 0.169, 0.154, 0.141, 0.130, 0.119, 0.107, 0.098, 0.088, 0.080, 0.071, 0.061, 0.051, 0.042, 0.033, 0.032, 0.032, 0.029, 0.034},  {1.067, 1.020, 0.945, 0.842, 0.730, 0.631, 0.551, 0.473, 0.412, 0.334, 0.319, 0.295, 0.265, 0.249, 0.236, 0.223, 0.212, 0.202, 0.192, 0.182, 0.174, 0.165, 0.158, 0.150, 0.140, 0.129, 0.120, 0.110, 0.109, 0.109, 0.101, 0.099}},
      {{1.054, 1.0742, 1.112, 1.1608, 1.2196, 1.2921, 1.3839, 1.5093, 1.6262, 1.6288, 1.7317, 1.8258, 1.9264, 2.025, 2.141, 2.2663, 2.361, 2.4589, 2.5693, 2.6659, 2.7595, 2.8628, 2.9668, 3.0644, 3.1682, 3.259, 3.3594, 3.4019, 3.4393},  {2.016, 1.968, 1.862, 1.718, 1.493, 1.206, 0.900, 0.601, 0.388, 0.247, 0.214, 0.187, 0.164, 0.145, 0.126, 0.108, 0.097, 0.086, 0.074, 0.066, 0.058, 0.051, 0.045, 0.040, 0.035, 0.031, 0.028, 0.026, 0.029},  {1.631, 1.597, 1.522, 1.387, 1.161, 1.037, 0.806, 0.525, 0.402, 0.277, 0.247, 0.223, 0.201, 0.183, 0.165, 0.148, 0.136, 0.125, 0.114, 0.105, 0.097, 0.089, 0.082, 0.076, 0.070, 0.066, 0.061, 0.059, 0.056},  {1.319, 1.297, 1.248, 1.156, 1.064, 0.914, 0.747, 0.574, 0.440, 0.331, 0.330, 0.283, 0.287, 0.248, 0.232, 0.217, 0.207, 0.198, 0.188, 0.180, 0.173, 0.166, 0.159, 0.154, 0.148, 0.143, 0.139, 0.137, 0.134},  {1.968, 1.848, 1.582, 1.266, 1.006, 0.801, 0.629, 0.468, 0.347, 0.241, 0.209, 0.182, 0.159, 0.140, 0.121, 0.103, 0.091, 0.080, 0.068, 0.059, 0.051, 0.043, 0.036, 0.028, 0.021, 0.016, 0.012, 0.009, 0.009},  {1.595, 1.506, 1.309, 1.072, 0.875, 0.718, 0.583, 0.453, 0.364, 0.272, 0.243, 0.218, 0.197, 0.179, 0.160, 0.143, 0.131, 0.119, 0.107, 0.098, 0.089, 0.080, 0.071, 0.061, 0.051, 0.043, 0.036, 0.032, 0.030},  {1.289, 1.228, 1.089, 0.920, 0.777, 0.660, 0.561, 0.461, 0.401, 0.345, 0.301, 0.296, 0.277, 0.262, 0.228, 0.233, 0.224, 0.192, 0.182, 0.173, 0.165, 0.157, 0.149, 0.139, 0.129, 0.120, 0.113, 0.108, 0.106}},
    }
};
//
// High Z = 0.02 table
const GE_QCRIT_TABLE QCRIT_GE_HIGH_Z = {
    {0.1, 0.13, 0.16, 0.2, 0.22, 0.25, 0.28, 0.32, 0.36, 0.4, 0.45, 0.5, 0.56, 0.63, 0.71, 0.8, 0.89, 1.0, 1.14, 1.3, 1.44, 1.6, 1.8, 2.0, 2.04, 2.5, 3.2, 4.0, 5.0, 6.3, 8.0, 10.0, 13.0, 16.0, 20.0, 25.0, 32.0, 40.0, 50.0, 63.0, 80.0, 100.0},
    {
      {{-0.8724, -0.8581, -0.8457, -0.8367, -0.8336},  {-334.448, -502.513, -502.513, 497.512, 332.226},  {-200.401, -250.627, -200.401, 990.099, 497.512},  {-90.992, -111.235, -100.100, -334.448, 100000.000},  {1.517, 1.531, 1.534, 1.536, 1.538},  {1.241, 1.252, 1.253, 1.255, 1.255},  {1.006, 1.013, 1.014, 1.015, 1.016}},
      {{-0.7771, -0.7662, -0.7581, -0.7545, -0.7546, -0.7295, -0.7239},  {-250.627, -1010.101, 100000.000, 100000.000, 100000.000, 0.941, 0.209},  {-143.062, -334.448, -502.513, -502.513, -334.448, -13.335, 199.601},  {-83.403, -143.062, -166.945, -143.062, -143.062, -12.822, 199.601},  {1.531, 1.534, 1.536, 1.538, 1.538, 0.828, 0.209},  {1.250, 1.253, 1.255, 1.256, 1.256, 0.708, 0.248},  {1.012, 1.014, 1.015, 1.017, 1.017, 0.605, 0.326}},
      {{-0.7068, -0.6996, -0.6958, -0.6953, -0.6949, -0.6642, -0.5622, -0.46, -0.3582, -0.256},  {332.226, -1010.101, 249.377, 100000.000, 497.512, 1.495, 0.600, 0.366, 0.204, 0.096},  {497.512, -334.448, 332.226, -334.448, -33.344, 249.377, 0.525, 0.330, 0.190, 0.147},  {-1010.101, -143.062, 990.099, -143.062, -26.323, 249.377, 0.459, 0.298, 0.264, 0.250},  {1.534, 1.538, 1.541, 1.543, 1.427, 1.186, 0.580, 0.354, 0.196, 0.096},  {1.252, 1.255, 1.258, 1.259, 1.171, 0.988, 0.509, 0.321, 0.183, 0.147},  {1.013, 1.016, 1.018, 1.019, 0.954, 0.819, 0.447, 0.291, 0.264, 0.250}},
      {{-0.6373, -0.6331, -0.6305, -0.6299, -0.6295, -0.6251, -0.6117, -0.5659, -0.4627, -0.3596, -0.256, -0.1531, -0.0464, 0.0491, 0.1517, 0.2617, 0.3627},  {-1010.101, 332.226, 100000.000, 497.512, 249.377, 100000.000, 249.377, 1.309, 0.823, 0.664, 0.575, 0.495, 0.415, 0.342, 0.258, 0.156, 0.042},  {-334.448, 497.512, -502.513, -83.403, -43.497, -35.727, -27.034, -16.396, -11.629, -11.629, 0.527, 0.527, 0.372, 0.527, 0.527, -11.629, 0.558},  {-143.062, 100000.000, -166.945, -52.659, -33.344, -27.785, -23.261, -15.154, -11.365, -11.365, 0.460, 0.460, 0.336, 0.460, 0.460, -11.365, 0.487},  {1.538, 1.543, 1.548, 1.550, 1.534, 1.484, 1.381, 1.151, 0.797, 0.649, 0.560, 0.479, 0.397, 0.323, 0.240, 0.142, 0.042},  {1.256, 1.259, 1.263, 1.266, 1.252, 1.215, 1.136, 0.960, 0.684, 0.565, 0.493, 0.427, 0.358, 0.296, 0.224, 0.136, 0.087},  {1.016, 1.019, 1.021, 1.024, 1.014, 0.987, 0.929, 0.798, 0.587, 0.494, 0.436, 0.382, 0.325, 0.272, 0.212, 0.204, 0.016}},
      {{-0.6087, -0.6047, -0.6021, -0.6015, -0.5993, -0.593, -0.5764, -0.5258, -0.4215, -0.3172, -0.2127, -0.1081, -0.0049, 0.0996, 0.2053, 0.3111, 0.4082, 0.5206, 0.6211},  {990.099, 497.512, 332.226, 100000.000, 249.377, 332.226, 990.099, 1.269, 0.870, 0.739, 0.674, 0.610, 0.543, 0.475, 0.403, 0.325, 0.247, 0.146, 0.030},  {-1010.101, 100000.000, -111.235, -47.642, 249.377, -32.268, -25.006, -15.876, 249.377, 0.637, 0.586, 0.534, 0.480, 0.424, 0.363, 0.297, 0.230, 0.139, 0.072},  {-200.401, -334.448, -62.539, -34.495, 249.377, -26.323, -21.281, -14.928, 249.377, 0.550, 0.510, 0.469, 0.426, 0.380, 0.330, 0.274, 0.215, 0.134, 0.022},  {1.541, 1.546, 1.550, 1.548, 1.522, 1.468, 1.364, 1.143, 0.844, 0.723, 0.659, 0.592, 0.523, 0.452, 0.378, 0.300, 0.224, 0.128, 0.030},  {1.258, 1.263, 1.266, 1.264, 1.244, 1.203, 1.124, 0.954, 0.720, 0.625, 0.574, 0.520, 0.464, 0.406, 0.344, 0.277, 0.211, 0.125, 0.071},  {1.018, 1.020, 1.024, 1.022, 1.008, 0.978, 0.920, 0.794, 0.615, 0.541, 0.502, 0.459, 0.414, 0.366, 0.315, 0.258, 0.200, 0.122, 0.019}},
      {{-0.5703, -0.5663, -0.5638, -0.5627, -0.5583, -0.5497, -0.5296, -0.4743, -0.3728, -0.2715, -0.1699, -0.0708, 0.0319, 0.1316, 0.2382, 0.3351, 0.4386, 0.546, 0.6432, 0.7393, 0.8399, 0.9457},  {332.226, 990.099, 332.226, 990.099, 332.226, 332.226, 71.378, 1.236, 0.931, 0.831, 0.789, 0.742, 0.685, 0.629, 0.566, 0.506, 0.439, 0.364, 0.292, 0.218, 0.136, 0.020},  {497.512, -1010.101, -50.025, -41.684, -35.727, -9.435, -9.435, 1.026, 0.789, 0.711, 0.678, 0.641, 0.596, 0.550, 0.500, 0.451, 0.395, 0.332, 0.271, 0.206, 0.132, 0.057},  {-1010.101, -200.401, -37.051, -32.268, -27.785, -9.435, -9.435, 0.846, 0.668, 0.609, 0.583, 0.555, 0.520, 0.484, 0.443, 0.404, 0.359, 0.306, 0.253, 0.196, 0.129, 0.021},  {1.546, 1.550, 1.550, 1.538, 1.506, 1.449, 1.344, 1.139, 0.904, 0.817, 0.775, 0.725, 0.665, 0.603, 0.537, 0.474, 0.404, 0.329, 0.258, 0.188, 0.113, 0.019},  {1.261, 1.266, 1.266, 1.256, 1.232, 1.189, 1.109, 0.951, 0.769, 0.700, 0.667, 0.628, 0.580, 0.531, 0.477, 0.426, 0.368, 0.304, 0.243, 0.181, 0.113, 0.055},  {1.020, 1.024, 1.024, 1.017, 0.999, 0.967, 0.908, 0.792, 0.653, 0.600, 0.575, 0.546, 0.509, 0.470, 0.427, 0.386, 0.338, 0.284, 0.232, 0.176, 0.113, 0.017}},
      {{-0.536, -0.5322, -0.5307, -0.5276, -0.5217, -0.5112, -0.4883, -0.4305, -0.3276, -0.2249, -0.1214, -0.0168, 0.0805, 0.1867, 0.2933, 0.3908, 0.4954, 0.593, 0.707, 0.8101, 0.9127, 1.0121, 1.1144, 1.2151},  {-1010.101, -1010.101, 990.099, 332.226, 100000.000, -1010.101, 1.634, 1.219, 0.963, 0.890, 0.864, 0.822, 0.774, 0.718, 0.660, 0.606, 0.544, 0.482, 0.407, 0.337, 0.263, 0.191, 0.115, 0.019},  {-334.448, -76.982, -47.642, -38.476, -33.344, -9.435, -9.435, 1.012, 0.815, 0.758, 0.737, 0.705, 0.667, 0.623, 0.577, 0.533, 0.483, 0.433, 0.371, 0.311, 0.247, 0.183, 0.114, 0.047},  {-143.062, -52.659, -34.495, -30.312, -27.034, -9.435, -9.435, 0.837, 0.688, 0.645, 0.630, 0.606, 0.577, 0.543, 0.507, 0.472, 0.432, 0.392, 0.340, 0.290, 0.235, 0.178, 0.114, 0.156},  {1.550, 1.550, 1.548, 1.527, 1.488, 1.431, 1.326, 1.139, 0.927, 0.876, 0.850, 0.806, 0.754, 0.693, 0.629, 0.569, 0.503, 0.439, 0.362, 0.292, 0.222, 0.155, 0.089, 0.013},  {1.264, 1.266, 1.264, 1.247, 1.218, 1.175, 1.095, 0.951, 0.786, 0.747, 0.727, 0.693, 0.652, 0.603, 0.553, 0.506, 0.452, 0.399, 0.335, 0.275, 0.213, 0.154, 0.092, 0.044},  {1.022, 1.024, 1.022, 1.010, 0.989, 0.958, 0.899, 0.792, 0.664, 0.637, 0.623, 0.597, 0.566, 0.529, 0.490, 0.452, 0.409, 0.367, 0.314, 0.263, 0.209, 0.155, 0.096, 0.016}},
      {{-0.4956, -0.4928, -0.4908, -0.4861, -0.4786, -0.4662, -0.4406, -0.3806, -0.2813, -0.1825, -0.085, 0.0174, 0.1168, 0.2151, 0.3155, 0.4085, 0.5101, 0.6091, 0.7101, 0.8133, 0.9053, 1.0115, 1.1044, 1.2088, 1.3081, 1.4077, 1.5028},  {100000.000, 497.512, 497.512, -502.513, 100000.000, 100000.000, 1.490, 1.209, 1.007, 0.944, 0.931, 0.903, 0.859, 0.811, 0.760, 0.712, 0.658, 0.607, 0.551, 0.488, 0.429, 0.359, 0.296, 0.224, 0.159, 0.092, 0.014},  {-9.347, -47.642, -41.684, -37.051, -32.268, -9.435, -9.435, 1.005, 0.849, 0.801, 0.791, 0.769, 0.735, 0.697, 0.658, 0.620, 0.576, 0.536, 0.491, 0.440, 0.391, 0.332, 0.278, 0.216, 0.157, 0.094, 0.039},  {-9.435, -35.727, -32.268, -29.420, -26.323, -9.435, -9.435, 0.831, 0.714, 0.678, 0.671, 0.655, 0.630, 0.601, 0.571, 0.542, 0.508, 0.477, 0.441, 0.400, 0.360, 0.311, 0.265, 0.211, 0.157, 0.098, 0.018},  {1.550, 1.550, 1.538, 1.511, 1.471, 1.412, 1.311, 1.143, 0.982, 0.929, 0.918, 0.889, 0.841, 0.787, 0.730, 0.676, 0.615, 0.559, 0.498, 0.432, 0.371, 0.301, 0.240, 0.175, 0.117, 0.062, 0.008},  {1.266, 1.266, 1.256, 1.235, 1.205, 1.160, 1.083, 0.954, 0.830, 0.789, 0.781, 0.758, 0.721, 0.679, 0.635, 0.592, 0.544, 0.499, 0.450, 0.396, 0.345, 0.286, 0.233, 0.175, 0.121, 0.068, 0.033},  {1.024, 1.024, 1.016, 1.002, 0.979, 0.947, 0.890, 0.794, 0.700, 0.669, 0.664, 0.648, 0.620, 0.589, 0.555, 0.523, 0.486, 0.450, 0.412, 0.368, 0.327, 0.276, 0.231, 0.179, 0.128, 0.076, 0.013}},
      {{-0.4601, -0.4588, -0.4549, -0.4491, -0.4404, -0.4267, -0.3993, -0.3391, -0.2412, -0.1446, -0.0441, 0.0538, 0.1524, 0.2466, 0.3437, 0.4443, 0.5423, 0.6364, 0.7361, 0.8389, 0.9307, 1.0371, 1.1311, 1.2251, 1.3315, 1.4219, 1.5193, 1.6162, 1.7182},  {-1010.101, 990.099, 100000.000, 497.512, -1010.101, 990.099, 1.447, 1.209, 1.039, 0.996, 0.992, 0.968, 0.926, 0.882, 0.836, 0.788, 0.737, 0.688, 0.634, 0.576, 0.521, 0.455, 0.394, 0.332, 0.262, 0.203, 0.143, 0.083, 0.011},  {-9.347, -9.347, -38.476, -35.727, -29.420, -9.347, -9.347, 1.005, 0.874, 0.841, 0.838, 0.820, 0.787, 0.754, 0.718, 0.680, 0.641, 0.602, 0.559, 0.513, 0.468, 0.414, 0.364, 0.312, 0.251, 0.199, 0.144, 0.088, 0.034},  {-9.347, -9.347, -30.312, -27.785, -24.396, -9.347, -9.347, 0.831, 0.734, 0.709, 0.707, 0.694, 0.670, 0.645, 0.618, 0.590, 0.560, 0.530, 0.497, 0.461, 0.426, 0.382, 0.341, 0.297, 0.245, 0.200, 0.149, 0.140, 0.134},  {1.553, 1.546, 1.524, 1.495, 1.456, 1.397, 1.300, 1.151, 1.015, 0.980, 0.979, 0.956, 0.909, 0.860, 0.808, 0.752, 0.694, 0.638, 0.578, 0.514, 0.456, 0.386, 0.325, 0.264, 0.198, 0.146, 0.094, 0.049, 0.004},  {1.266, 1.261, 1.247, 1.224, 1.193, 1.149, 1.075, 0.961, 0.855, 0.829, 0.828, 0.810, 0.775, 0.737, 0.696, 0.653, 0.608, 0.564, 0.517, 0.465, 0.417, 0.360, 0.309, 0.257, 0.199, 0.152, 0.103, 0.057, 0.023},  {1.025, 1.021, 1.010, 0.993, 0.971, 0.939, 0.885, 0.799, 0.720, 0.700, 0.700, 0.688, 0.662, 0.634, 0.604, 0.571, 0.538, 0.504, 0.467, 0.427, 0.389, 0.342, 0.300, 0.256, 0.205, 0.162, 0.115, 0.069, 0.010}},
      {{-0.4265, -0.4218, -0.4158, -0.4079, -0.3971, -0.3809, -0.3506, -0.2891, -0.1896, -0.0899, 0.0089, 0.1082, 0.2107, 0.3082, 0.4071, 0.5078, 0.6032, 0.7087, 0.8087, 0.9117, 1.0036, 1.1101, 1.2042, 1.3121, 1.4059, 1.5111, 1.5998, 1.8083, 1.9051},  {990.099, 249.377, 497.512, -502.513, 249.377, 1.490, 1.335, 1.148, 1.018, 1.004, 1.014, 0.996, 0.958, 0.918, 0.876, 0.831, 0.784, 0.731, 0.680, 0.624, 0.571, 0.507, 0.448, 0.378, 0.316, 0.247, 0.192, 0.069, 0.006},  {-9.347, -9.347, -9.347, 249.377, -9.347, -9.347, -17.244, 0.959, 0.858, 0.847, 0.855, 0.842, 0.813, 0.782, 0.750, 0.715, 0.679, 0.637, 0.597, 0.552, 0.510, 0.459, 0.410, 0.352, 0.301, 0.241, 0.193, 0.077, 0.030},  {-9.347, -9.347, -9.347, 249.377, -9.347, -9.347, -15.876, 0.797, 0.721, 0.714, 0.720, 0.712, 0.691, 0.668, 0.644, 0.618, 0.591, 0.559, 0.529, 0.495, 0.462, 0.421, 0.382, 0.335, 0.291, 0.240, 0.197, 0.135, 0.013},  {1.534, 1.508, 1.479, 1.445, 1.401, 1.340, 1.247, 1.109, 0.999, 0.990, 1.002, 0.982, 0.940, 0.894, 0.845, 0.792, 0.737, 0.676, 0.617, 0.554, 0.496, 0.428, 0.366, 0.296, 0.237, 0.174, 0.126, 0.031, 0.002},  {1.252, 1.235, 1.212, 1.185, 1.152, 1.106, 1.034, 0.928, 0.842, 0.837, 0.846, 0.832, 0.799, 0.764, 0.726, 0.685, 0.643, 0.595, 0.549, 0.499, 0.453, 0.397, 0.346, 0.287, 0.236, 0.180, 0.136, 0.041, 0.003},  {1.014, 1.001, 0.985, 0.965, 0.941, 0.907, 0.854, 0.775, 0.710, 0.706, 0.714, 0.705, 0.681, 0.656, 0.628, 0.598, 0.566, 0.530, 0.495, 0.457, 0.421, 0.421, 0.335, 0.286, 0.242, 0.193, 0.152, 0.056, 0.008}},
      {{-0.3833, -0.3768, -0.3686, -0.3584, -0.3449, -0.3253, -0.2916, -0.2299, -0.1287, -0.0276, 0.0708, 0.1718, 0.2754, 0.3766, 0.4772, 0.5775, 0.6809, 0.7754, 0.8745, 0.9768, 1.0813, 1.1741, 1.2811, 1.3751, 1.4821, 1.5875, 1.6899, 1.7877, 1.8895, 1.9868, 2.0867},  {1.869, 1.678, 1.582, 1.508, 1.427, 1.330, 1.206, 1.065, 0.979, 1.011, 1.034, 1.021, 0.990, 0.954, 0.915, 0.873, 0.824, 0.781, 0.733, 0.681, 0.623, 0.569, 0.503, 0.443, 0.373, 0.302, 0.236, 0.175, 0.114, 0.058, 0.005},  {2.045, 1.361, 1.289, 1.233, 1.171, 1.098, 1.003, 0.894, 0.828, 0.853, 0.871, 0.862, 0.838, 0.810, 0.781, 0.748, 0.711, 0.678, 0.641, 0.599, 0.554, 0.511, 0.458, 0.409, 0.351, 0.291, 0.235, 0.180, 0.125, 0.070, 0.026},  {1.555, 1.092, 1.039, 0.999, 0.954, 0.901, 0.830, 0.748, 0.699, 0.718, 0.733, 0.728, 0.711, 0.691, 0.668, 0.645, 0.617, 0.592, 0.565, 0.533, 0.499, 0.465, 0.424, 0.385, 0.337, 0.288, 0.240, 0.191, 0.141, 0.086, 0.126},  {1.479, 1.449, 1.412, 1.372, 1.321, 1.255, 1.160, 1.039, 0.965, 0.999, 1.021, 1.006, 0.968, 0.924, 0.878, 0.828, 0.770, 0.719, 0.663, 0.602, 0.538, 0.478, 0.408, 0.346, 0.276, 0.209, 0.150, 0.097, 0.051, 0.012, 0.001},  {1.211, 1.188, 1.161, 1.130, 1.092, 1.041, 0.968, 0.875, 0.817, 0.843, 0.861, 0.850, 0.822, 0.789, 0.753, 0.715, 0.671, 0.631, 0.588, 0.540, 0.489, 0.441, 0.384, 0.332, 0.273, 0.216, 0.163, 0.113, 0.067, 0.023, 0.012},  {0.984, 0.968, 0.948, 0.924, 0.896, 0.859, 0.804, 0.734, 0.691, 0.712, 0.727, 0.719, 0.700, 0.676, 0.650, 0.622, 0.589, 0.560, 0.528, 0.492, 0.452, 0.415, 0.370, 0.329, 0.279, 0.231, 0.183, 0.138, 0.091, 0.043, 0.009}},
      {{-0.3385, -0.3304, -0.3205, -0.3085, -0.2931, -0.2718, -0.2375, -0.1818, -0.1056, -0.0284, 0.0736, 0.1762, 0.2766, 0.3844, 0.4792, 0.5853, 0.6837, 0.7944, 0.8911, 0.9914, 1.0942, 1.1988, 1.3047, 1.3978, 1.504, 1.6092, 1.7122, 1.8121, 1.919, 2.0177, 2.1214, 2.2211},  {1.522, 1.468, 1.410, 1.350, 1.280, 1.199, 1.100, 0.997, 0.951, 0.989, 1.045, 1.053, 1.033, 1.001, 0.968, 0.928, 0.887, 0.838, 0.795, 0.747, 0.693, 0.636, 0.573, 0.514, 0.445, 0.373, 0.302, 0.238, 0.166, 0.109, 0.046, 0.004},  {1.242, 1.203, 1.159, 1.112, 1.059, 0.997, 0.921, 0.841, 0.806, 0.836, 0.879, 0.887, 0.872, 0.847, 0.822, 0.792, 0.760, 0.723, 0.690, 0.653, 0.611, 0.566, 0.516, 0.469, 0.413, 0.354, 0.294, 0.239, 0.177, 0.124, 0.061, 0.024},  {1.007, 0.978, 0.946, 0.911, 0.872, 0.826, 0.769, 0.708, 0.682, 0.706, 0.740, 0.746, 0.736, 0.718, 0.701, 0.679, 0.655, 0.628, 0.604, 0.576, 0.545, 0.511, 0.472, 0.436, 0.391, 0.344, 0.295, 0.248, 0.194, 0.145, 0.084, 0.121},  {1.399, 1.362, 1.323, 1.279, 1.225, 1.159, 1.072, 0.979, 0.939, 0.978, 1.032, 1.036, 1.012, 0.972, 0.931, 0.883, 0.832, 0.773, 0.722, 0.665, 0.603, 0.537, 0.468, 0.405, 0.333, 0.263, 0.197, 0.138, 0.080, 0.033, 0.008, 0.001},  {1.149, 1.124, 1.093, 1.058, 1.017, 0.966, 0.900, 0.828, 0.796, 0.827, 0.870, 0.874, 0.856, 0.826, 0.795, 0.758, 0.719, 0.675, 0.635, 0.591, 0.542, 0.490, 0.435, 0.384, 0.325, 0.265, 0.209, 0.155, 0.101, 0.051, 0.019, 0.011},  {0.940, 0.920, 0.898, 0.872, 0.842, 0.804, 0.754, 0.698, 0.675, 0.699, 0.733, 0.738, 0.725, 0.704, 0.682, 0.656, 0.628, 0.595, 0.566, 0.533, 0.497, 0.458, 0.415, 0.375, 0.327, 0.278, 0.230, 0.182, 0.132, 0.082, 0.043, 0.012}},
      {{-0.2837, -0.2743, -0.2632, -0.2502, -0.2344, -0.2139, -0.1855, -0.1451, -0.0687, 0.0091, 0.1072, 0.2088, 0.3064, 0.4031, 0.5016, 0.6035, 0.699, 0.7957, 0.9009, 0.9988, 1.0995, 1.1894, 1.2934, 1.3983, 1.4903, 1.595, 1.6855, 1.787, 1.8856, 1.9917, 2.0802},  {1.330, 1.284, 1.235, 1.183, 1.129, 1.067, 1.000, 0.941, 0.917, 0.979, 1.055, 1.075, 1.065, 1.042, 1.011, 0.977, 0.940, 0.897, 0.853, 0.810, 0.761, 0.715, 0.657, 0.594, 0.536, 0.465, 0.401, 0.329, 0.262, 0.187, 0.133},  {1.098, 1.062, 1.025, 0.986, 0.943, 0.895, 0.843, 0.797, 0.779, 0.828, 0.887, 0.904, 0.897, 0.879, 0.856, 0.831, 0.803, 0.769, 0.736, 0.703, 0.665, 0.629, 0.584, 0.535, 0.488, 0.432, 0.380, 0.319, 0.263, 0.198, 0.150},  {0.900, 0.874, 0.847, 0.818, 0.786, 0.750, 0.710, 0.674, 0.661, 0.700, 0.746, 0.760, 0.755, 0.742, 0.727, 0.708, 0.688, 0.664, 0.640, 0.615, 0.588, 0.561, 0.528, 0.490, 0.454, 0.410, 0.369, 0.321, 0.274, 0.219, 0.176},  {1.272, 1.235, 1.193, 1.149, 1.100, 1.044, 0.982, 0.927, 0.906, 0.967, 1.041, 1.058, 1.042, 1.011, 0.973, 0.929, 0.884, 0.832, 0.778, 0.726, 0.668, 0.613, 0.547, 0.477, 0.413, 0.340, 0.277, 0.209, 0.147, 0.085, 0.037},  {1.054, 1.026, 0.994, 0.960, 0.922, 0.878, 0.830, 0.787, 0.770, 0.819, 0.877, 0.891, 0.879, 0.857, 0.828, 0.796, 0.761, 0.721, 0.680, 0.639, 0.595, 0.552, 0.500, 0.444, 0.393, 0.333, 0.281, 0.223, 0.167, 0.109, 0.060},  {0.869, 0.847, 0.824, 0.799, 0.770, 0.737, 0.700, 0.667, 0.655, 0.693, 0.739, 0.751, 0.743, 0.728, 0.708, 0.684, 0.660, 0.631, 0.602, 0.572, 0.539, 0.508, 0.469, 0.426, 0.386, 0.338, 0.296, 0.247, 0.199, 0.147, 0.099}},
      {{-0.2237, -0.2149, -0.2052, -0.1944, -0.1821, -0.1667, -0.144, -0.1035, -0.027, 0.048, 0.1496, 0.2515, 0.3555, 0.456, 0.562, 0.6631, 0.7607, 0.8678, 0.9605, 1.0696, 1.1694, 1.2711, 1.3739, 1.4774, 1.5808, 1.6835, 1.7847, 1.8837, 1.9914, 2.0936, 2.1633, 2.2587, 2.3394, 2.3607},  {1.142, 1.107, 1.072, 1.036, 1.000, 0.962, 0.919, 0.873, 0.876, 0.961, 1.065, 1.096, 1.091, 1.070, 1.039, 1.007, 0.970, 0.924, 0.887, 0.840, 0.794, 0.741, 0.684, 0.622, 0.555, 0.483, 0.408, 0.334, 0.255, 0.180, 0.038, 0.011, 0.004, 0.002},  {0.953, 0.927, 0.899, 0.872, 0.843, 0.814, 0.780, 0.745, 0.746, 0.814, 0.895, 0.921, 0.917, 0.901, 0.879, 0.854, 0.826, 0.792, 0.763, 0.728, 0.692, 0.651, 0.608, 0.559, 0.506, 0.448, 0.387, 0.327, 0.260, 0.196, 0.053, 0.023, 0.011, 0.005},  {0.793, 0.773, 0.752, 0.732, 0.710, 0.687, 0.662, 0.635, 0.636, 0.689, 0.752, 0.772, 0.770, 0.760, 0.744, 0.727, 0.707, 0.682, 0.662, 0.636, 0.610, 0.581, 0.548, 0.512, 0.471, 0.427, 0.379, 0.331, 0.275, 0.224, 0.073, 0.044, 0.025, 0.013},  {1.114, 1.083, 1.050, 1.017, 0.983, 0.948, 0.907, 0.862, 0.863, 0.948, 1.049, 1.076, 1.064, 1.035, 0.996, 0.955, 0.908, 0.853, 0.806, 0.748, 0.691, 0.629, 0.562, 0.492, 0.419, 0.344, 0.271, 0.202, 0.130, 0.068, 0.007, 0.002, 0.001, 0.000},  {0.933, 0.908, 0.883, 0.858, 0.831, 0.803, 0.770, 0.736, 0.737, 0.804, 0.884, 0.906, 0.897, 0.876, 0.847, 0.816, 0.781, 0.738, 0.703, 0.659, 0.615, 0.567, 0.515, 0.459, 0.400, 0.339, 0.278, 0.219, 0.154, 0.095, 0.017, 0.019, 0.019, 0.019},  {0.778, 0.760, 0.741, 0.721, 0.701, 0.679, 0.654, 0.628, 0.629, 0.682, 0.744, 0.762, 0.757, 0.743, 0.723, 0.701, 0.676, 0.646, 0.621, 0.589, 0.557, 0.522, 0.483, 0.441, 0.396, 0.348, 0.299, 0.250, 0.195, 0.141, 0.039, 0.043, 0.043, 0.043}},
      {{-0.1759, -0.1695, -0.1619, -0.1529, -0.1415, -0.1263, -0.1025, -0.0658, 0.0099, 0.0862, 0.186, 0.2839, 0.3856, 0.4854, 0.5831, 0.6904, 0.7825, 0.891, 0.9806, 1.0863, 1.1837, 1.2832, 1.3841, 1.4857, 1.5876, 1.6889, 1.7891, 1.8876, 1.9836, 2.0873, 2.1719},  {1.001, 0.979, 0.955, 0.929, 0.901, 0.870, 0.833, 0.798, 0.812, 0.925, 1.066, 1.112, 1.116, 1.100, 1.076, 1.045, 1.014, 0.967, 0.934, 0.890, 0.847, 0.800, 0.747, 0.689, 0.625, 0.555, 0.480, 0.402, 0.328, 0.241, 0.096},  {0.844, 0.827, 0.808, 0.788, 0.766, 0.741, 0.712, 0.685, 0.696, 0.786, 0.897, 0.934, 0.936, 0.925, 0.907, 0.883, 0.861, 0.825, 0.800, 0.767, 0.735, 0.698, 0.657, 0.612, 0.562, 0.508, 0.448, 0.385, 0.325, 0.251, 0.115},  {0.710, 0.697, 0.683, 0.668, 0.651, 0.631, 0.609, 0.588, 0.597, 0.668, 0.753, 0.782, 0.785, 0.778, 0.766, 0.749, 0.733, 0.708, 0.690, 0.667, 0.643, 0.617, 0.588, 0.555, 0.517, 0.476, 0.430, 0.381, 0.334, 0.275, 0.142},  {0.985, 0.964, 0.942, 0.917, 0.889, 0.858, 0.821, 0.787, 0.799, 0.910, 1.048, 1.091, 1.087, 1.064, 1.032, 0.990, 0.951, 0.894, 0.851, 0.797, 0.743, 0.684, 0.620, 0.552, 0.480, 0.405, 0.329, 0.254, 0.185, 0.111, 0.019},  {0.832, 0.816, 0.797, 0.778, 0.756, 0.732, 0.703, 0.676, 0.686, 0.774, 0.883, 0.917, 0.915, 0.898, 0.874, 0.843, 0.814, 0.770, 0.738, 0.697, 0.656, 0.611, 0.562, 0.509, 0.452, 0.391, 0.329, 0.267, 0.206, 0.140, 0.038},  {0.701, 0.689, 0.675, 0.660, 0.643, 0.625, 0.602, 0.581, 0.590, 0.660, 0.745, 0.771, 0.771, 0.760, 0.743, 0.723, 0.702, 0.671, 0.649, 0.619, 0.590, 0.558, 0.522, 0.482, 0.439, 0.393, 0.344, 0.294, 0.244, 0.189, 0.073}},
      {{-0.138, -0.131, -0.1221, -0.1116, -0.0984, -0.0809, -0.0558, -0.0258, 0.0571, 0.1393, 0.2406, 0.3444, 0.4457, 0.5449, 0.644, 0.7523, 0.8506, 0.9542, 1.0547, 1.1482, 1.2566, 1.3549, 1.4544, 1.5544, 1.6544, 1.7662, 1.8643, 1.9605, 2.0653, 2.1653, 2.1926},  {0.890, 0.867, 0.842, 0.816, 0.787, 0.756, 0.724, 0.702, 0.732, 0.899, 1.081, 1.132, 1.135, 1.121, 1.098, 1.067, 1.033, 0.990, 0.953, 0.916, 0.868, 0.821, 0.769, 0.711, 0.647, 0.567, 0.489, 0.410, 0.323, 0.231, 0.144},  {0.757, 0.739, 0.719, 0.698, 0.676, 0.651, 0.625, 0.607, 0.631, 0.766, 0.908, 0.949, 0.952, 0.942, 0.924, 0.902, 0.876, 0.843, 0.816, 0.787, 0.751, 0.716, 0.676, 0.631, 0.581, 0.519, 0.457, 0.393, 0.322, 0.246, 0.163},  {0.643, 0.629, 0.614, 0.598, 0.580, 0.561, 0.540, 0.526, 0.546, 0.652, 0.762, 0.794, 0.797, 0.791, 0.779, 0.763, 0.745, 0.722, 0.703, 0.683, 0.658, 0.632, 0.603, 0.571, 0.534, 0.488, 0.440, 0.392, 0.335, 0.277, 0.193},  {0.878, 0.855, 0.831, 0.804, 0.776, 0.745, 0.712, 0.689, 0.716, 0.882, 1.060, 1.107, 1.103, 1.080, 1.048, 1.008, 0.963, 0.911, 0.863, 0.814, 0.755, 0.696, 0.633, 0.564, 0.492, 0.406, 0.327, 0.251, 0.170, 0.095, 0.042},  {0.748, 0.730, 0.711, 0.690, 0.667, 0.642, 0.616, 0.598, 0.620, 0.753, 0.894, 0.930, 0.928, 0.912, 0.887, 0.858, 0.824, 0.784, 0.748, 0.712, 0.667, 0.622, 0.573, 0.520, 0.463, 0.394, 0.329, 0.266, 0.195, 0.128, 0.067},  {0.637, 0.623, 0.608, 0.591, 0.574, 0.554, 0.534, 0.519, 0.538, 0.643, 0.752, 0.781, 0.781, 0.771, 0.755, 0.734, 0.711, 0.683, 0.658, 0.633, 0.600, 0.568, 0.533, 0.494, 0.451, 0.398, 0.348, 0.298, 0.240, 0.184, 0.112}},
      {{-0.0999, -0.091, -0.0802, -0.0675, -0.0515, -0.0326, -0.0084, 0.0141, 0.1013, 0.1875, 0.2858, 0.3838, 0.4842, 0.5821, 0.6848, 0.7751, 0.8806, 0.9732, 1.071, 1.1737, 1.2681, 1.3766, 1.4743, 1.5727, 1.6711, 1.7691, 1.8661, 1.9614, 2.0657, 2.1656, 2.2708},  {0.770, 0.740, 0.712, 0.683, 0.652, 0.623, 0.597, 0.586, 0.627, 0.864, 1.092, 1.148, 1.155, 1.142, 1.119, 1.094, 1.062, 1.024, 0.989, 0.950, 0.910, 0.860, 0.810, 0.755, 0.694, 0.626, 0.550, 0.468, 0.379, 0.279, 0.142},  {0.662, 0.637, 0.615, 0.592, 0.567, 0.543, 0.522, 0.513, 0.548, 0.738, 0.917, 0.961, 0.966, 0.957, 0.941, 0.923, 0.898, 0.869, 0.843, 0.814, 0.784, 0.746, 0.709, 0.667, 0.620, 0.567, 0.507, 0.442, 0.370, 0.288, 0.166},  {0.569, 0.550, 0.532, 0.513, 0.494, 0.475, 0.458, 0.450, 0.479, 0.631, 0.769, 0.803, 0.808, 0.803, 0.792, 0.779, 0.762, 0.742, 0.724, 0.704, 0.683, 0.657, 0.630, 0.600, 0.565, 0.526, 0.481, 0.432, 0.376, 0.312, 0.205},  {0.759, 0.728, 0.700, 0.671, 0.640, 0.610, 0.584, 0.573, 0.611, 0.845, 1.070, 1.120, 1.119, 1.098, 1.066, 1.033, 0.990, 0.942, 0.897, 0.846, 0.795, 0.732, 0.671, 0.604, 0.533, 0.457, 0.378, 0.298, 0.212, 0.129, 0.027},  {0.653, 0.629, 0.606, 0.582, 0.558, 0.534, 0.512, 0.503, 0.535, 0.724, 0.900, 0.940, 0.941, 0.925, 0.902, 0.877, 0.845, 0.809, 0.776, 0.737, 0.699, 0.651, 0.604, 0.553, 0.497, 0.437, 0.373, 0.307, 0.235, 0.162, 0.140},  {0.563, 0.543, 0.525, 0.507, 0.487, 0.468, 0.450, 0.443, 0.470, 0.621, 0.758, 0.789, 0.791, 0.782, 0.766, 0.749, 0.727, 0.702, 0.678, 0.652, 0.625, 0.591, 0.558, 0.521, 0.479, 0.434, 0.385, 0.335, 0.276, 0.217, 0.188}},
      {{-0.0471, -0.0352, -0.0217, -0.0063, 0.0105, 0.0286, 0.0467, 0.0627, 0.1611, 0.26, 0.3588, 0.4567, 0.5574, 0.6556, 0.7514, 0.8533, 0.9495, 1.0511, 1.1494, 1.2515, 1.3449, 1.4519, 1.5481, 1.6449, 1.7414, 1.8492, 1.9438, 2.0478, 2.1479, 2.2425},  {0.585, 0.541, 0.508, 0.477, 0.450, 0.428, 0.417, 0.416, 0.479, 0.898, 1.130, 1.171, 1.171, 1.156, 1.134, 1.106, 1.079, 1.038, 1.003, 0.964, 0.924, 0.875, 0.826, 0.770, 0.709, 0.630, 0.550, 0.456, 0.359, 0.256},  {0.512, 0.476, 0.449, 0.423, 0.406, 0.392, 0.381, 0.373, 0.426, 0.765, 0.947, 0.979, 0.980, 0.969, 0.952, 0.932, 0.912, 0.881, 0.855, 0.826, 0.796, 0.759, 0.722, 0.680, 0.633, 0.571, 0.509, 0.434, 0.355, 0.271},  {0.476, 0.464, 0.452, 0.441, 0.430, 0.419, 0.410, 0.403, 0.379, 0.652, 0.792, 0.818, 0.820, 0.812, 0.801, 0.787, 0.773, 0.752, 0.734, 0.714, 0.693, 0.668, 0.641, 0.611, 0.577, 0.532, 0.485, 0.429, 0.367, 0.304},  {0.573, 0.529, 0.497, 0.465, 0.437, 0.416, 0.404, 0.403, 0.463, 0.876, 1.104, 1.138, 1.131, 1.107, 1.076, 1.038, 1.000, 0.950, 0.904, 0.852, 0.801, 0.737, 0.676, 0.609, 0.537, 0.450, 0.368, 0.277, 0.187, 0.105},  {0.503, 0.467, 0.440, 0.420, 0.405, 0.392, 0.381, 0.373, 0.413, 0.748, 0.928, 0.955, 0.951, 0.933, 0.911, 0.883, 0.854, 0.816, 0.782, 0.743, 0.705, 0.657, 0.610, 0.558, 0.502, 0.433, 0.366, 0.292, 0.215, 0.142},  {0.476, 0.463, 0.452, 0.441, 0.430, 0.419, 0.410, 0.403, 0.371, 0.641, 0.779, 0.801, 0.799, 0.789, 0.773, 0.754, 0.735, 0.709, 0.685, 0.658, 0.632, 0.598, 0.564, 0.527, 0.486, 0.434, 0.383, 0.326, 0.263, 0.205}},
      {{0.0304, 0.0434, 0.0575, 0.0721, 0.0852, 0.097, 0.109, 0.1228, 0.1591, 0.1923, 0.2657, 0.3392, 0.4361, 0.5356, 0.6357, 0.731, 0.8312, 0.9342, 1.0339, 1.1293, 1.2283, 1.3304, 1.4231, 1.5288, 1.6233, 1.718, 1.824, 1.9174, 2.0207, 2.1209, 2.2162, 2.2694},  {0.411, 0.393, 0.380, 0.367, 0.358, 0.351, 0.345, 0.341, 0.335, 0.307, 0.363, 0.991, 1.171, 1.196, 1.190, 1.172, 1.148, 1.120, 1.082, 1.052, 1.016, 0.976, 0.936, 0.887, 0.838, 0.782, 0.711, 0.639, 0.545, 0.451, 0.343, 0.247},  {0.411, 0.397, 0.386, 0.375, 0.367, 0.361, 0.356, 0.352, 0.345, 0.321, 0.327, 0.838, 0.979, 0.999, 0.995, 0.981, 0.963, 0.942, 0.915, 0.892, 0.866, 0.836, 0.806, 0.769, 0.732, 0.690, 0.636, 0.579, 0.506, 0.431, 0.343, 0.262},  {0.441, 0.430, 0.422, 0.414, 0.407, 0.402, 0.397, 0.393, 0.385, 0.367, 0.354, 0.709, 0.817, 0.833, 0.831, 0.822, 0.810, 0.796, 0.777, 0.761, 0.743, 0.723, 0.702, 0.677, 0.651, 0.621, 0.582, 0.540, 0.486, 0.428, 0.360, 0.294},  {0.411, 0.393, 0.379, 0.367, 0.357, 0.350, 0.344, 0.340, 0.335, 0.306, 0.346, 0.965, 1.140, 1.157, 1.143, 1.117, 1.083, 1.045, 0.998, 0.956, 0.908, 0.856, 0.804, 0.740, 0.678, 0.611, 0.528, 0.449, 0.355, 0.263, 0.169, 0.103},  {0.411, 0.397, 0.385, 0.375, 0.367, 0.360, 0.356, 0.351, 0.344, 0.321, 0.316, 0.819, 0.956, 0.971, 0.961, 0.942, 0.917, 0.888, 0.853, 0.822, 0.786, 0.747, 0.708, 0.660, 0.613, 0.561, 0.496, 0.433, 0.358, 0.281, 0.201, 0.206},  {0.441, 0.430, 0.422, 0.413, 0.407, 0.401, 0.397, 0.393, 0.384, 0.366, 0.353, 0.696, 0.802, 0.814, 0.808, 0.796, 0.779, 0.759, 0.735, 0.714, 0.690, 0.663, 0.636, 0.602, 0.569, 0.532, 0.484, 0.437, 0.380, 0.320, 0.257, 0.288}},
      {{0.1105, 0.1305, 0.1549, 0.1838, 0.2138, 0.2392, 0.2629, 0.3128, 0.4025, 0.4475, 0.5442, 0.6436, 0.7432, 0.8343, 0.9409, 1.0351, 1.135, 1.2337, 1.3223, 1.424, 1.5156, 1.6196, 1.7123, 1.8163, 1.9081, 2.0099, 2.1091, 2.2035, 2.2922},  {0.355, 0.338, 0.325, 0.314, 0.308, 0.309, 0.317, 0.258, 1.039, 1.161, 1.217, 1.214, 1.200, 1.179, 1.151, 1.125, 1.089, 1.055, 1.022, 0.981, 0.941, 0.890, 0.840, 0.775, 0.709, 0.622, 0.527, 0.426, 0.281},  {0.369, 0.355, 0.343, 0.333, 0.327, 0.326, 0.330, 0.283, 0.876, 0.971, 1.014, 1.012, 1.003, 0.987, 0.967, 0.948, 0.921, 0.896, 0.871, 0.841, 0.811, 0.773, 0.735, 0.686, 0.635, 0.568, 0.493, 0.412, 0.292},  {0.413, 0.402, 0.392, 0.383, 0.377, 0.373, 0.373, 0.339, 0.739, 0.810, 0.844, 0.844, 0.838, 0.827, 0.814, 0.801, 0.782, 0.765, 0.748, 0.728, 0.708, 0.682, 0.655, 0.620, 0.584, 0.534, 0.478, 0.416, 0.321},  {0.354, 0.338, 0.324, 0.313, 0.307, 0.308, 0.316, 0.258, 1.010, 1.130, 1.178, 1.167, 1.145, 1.116, 1.078, 1.042, 0.995, 0.949, 0.905, 0.851, 0.798, 0.733, 0.670, 0.591, 0.515, 0.423, 0.328, 0.233, 0.126},  {0.368, 0.354, 0.342, 0.332, 0.326, 0.325, 0.329, 0.282, 0.855, 0.948, 0.985, 0.978, 0.962, 0.942, 0.912, 0.886, 0.852, 0.818, 0.785, 0.745, 0.705, 0.656, 0.607, 0.547, 0.488, 0.414, 0.337, 0.257, 0.161},  {0.413, 0.402, 0.392, 0.383, 0.376, 0.373, 0.373, 0.339, 0.724, 0.795, 0.825, 0.821, 0.811, 0.797, 0.778, 0.760, 0.736, 0.713, 0.691, 0.663, 0.636, 0.602, 0.568, 0.524, 0.481, 0.425, 0.366, 0.302, 0.223}},
      {{0.1554, 0.1773, 0.2091, 0.2456, 0.2816, 0.313, 0.3424, 0.3577, 0.3763, 0.4846, 0.5853, 0.6828, 0.7826, 0.8808, 0.9824, 1.0763, 1.1824, 1.2745, 1.3813, 1.4804, 1.5813, 1.6721, 1.7746, 1.8768, 1.978, 2.0771, 2.172, 2.2663, 2.3611},  {0.334, 0.321, 0.308, 0.297, 0.290, 0.288, 0.294, 0.296, 0.255, 1.024, 1.199, 1.221, 1.215, 1.196, 1.171, 1.132, 1.100, 1.068, 1.029, 0.989, 0.944, 0.899, 0.842, 0.776, 0.697, 0.607, 0.513, 0.352, 0.250},  {0.353, 0.341, 0.329, 0.319, 0.312, 0.309, 0.312, 0.312, 0.280, 0.864, 1.001, 1.018, 1.014, 1.001, 0.982, 0.953, 0.929, 0.907, 0.877, 0.847, 0.814, 0.781, 0.738, 0.687, 0.627, 0.557, 0.483, 0.351, 0.267},  {0.403, 0.393, 0.383, 0.374, 0.367, 0.362, 0.361, 0.359, 0.338, 0.730, 0.834, 0.848, 0.847, 0.838, 0.825, 0.805, 0.789, 0.773, 0.754, 0.734, 0.711, 0.688, 0.658, 0.623, 0.579, 0.527, 0.471, 0.367, 0.305},  {0.333, 0.320, 0.307, 0.296, 0.289, 0.287, 0.294, 0.295, 0.254, 0.991, 1.157, 1.172, 1.157, 1.129, 1.094, 1.045, 1.001, 0.960, 0.907, 0.853, 0.794, 0.736, 0.664, 0.584, 0.496, 0.402, 0.307, 0.182, 0.095},  {0.352, 0.340, 0.329, 0.318, 0.311, 0.309, 0.312, 0.311, 0.280, 0.840, 0.971, 0.982, 0.972, 0.951, 0.926, 0.890, 0.857, 0.826, 0.787, 0.747, 0.702, 0.658, 0.604, 0.543, 0.473, 0.398, 0.320, 0.212, 0.133},  {0.402, 0.393, 0.383, 0.373, 0.366, 0.362, 0.361, 0.359, 0.337, 0.714, 0.814, 0.824, 0.818, 0.805, 0.788, 0.763, 0.741, 0.720, 0.693, 0.666, 0.636, 0.605, 0.567, 0.522, 0.472, 0.414, 0.354, 0.266, 0.204}},
      {{0.1695, 0.1998, 0.2407, 0.2868, 0.3377, 0.3866, 0.4213, 0.4504, 0.5386, 0.54, 0.6345, 0.7287, 0.8225, 0.9176, 1.0088, 1.1055, 1.2078, 1.3049, 1.4075, 1.5142, 1.6236, 1.7235, 1.8239, 1.9349, 2.0334, 2.1394, 2.1921, 2.2893, 2.3807},  {0.338, 0.324, 0.309, 0.294, 0.280, 0.272, 0.273, 0.251, 0.301, 1.036, 1.203, 1.221, 1.217, 1.200, 1.179, 1.151, 1.114, 1.081, 1.045, 1.002, 0.953, 0.903, 0.847, 0.772, 0.692, 0.591, 0.490, 0.378, 0.269},  {0.356, 0.344, 0.330, 0.317, 0.304, 0.296, 0.295, 0.277, 0.306, 0.874, 1.004, 1.018, 1.015, 1.004, 0.988, 0.967, 0.940, 0.917, 0.889, 0.858, 0.822, 0.784, 0.742, 0.685, 0.624, 0.545, 0.463, 0.373, 0.284},  {0.406, 0.396, 0.384, 0.372, 0.361, 0.353, 0.349, 0.335, 0.341, 0.737, 0.837, 0.849, 0.847, 0.840, 0.830, 0.816, 0.797, 0.781, 0.763, 0.742, 0.717, 0.692, 0.662, 0.622, 0.578, 0.519, 0.453, 0.386, 0.319},  {0.338, 0.323, 0.308, 0.292, 0.278, 0.270, 0.272, 0.250, 0.299, 1.001, 1.160, 1.170, 1.156, 1.131, 1.101, 1.063, 1.014, 0.971, 0.921, 0.864, 0.799, 0.735, 0.663, 0.575, 0.486, 0.381, 0.295, 0.198, 0.107},  {0.356, 0.343, 0.329, 0.316, 0.303, 0.295, 0.294, 0.276, 0.305, 0.848, 0.972, 0.981, 0.972, 0.954, 0.931, 0.903, 0.867, 0.835, 0.797, 0.755, 0.707, 0.658, 0.604, 0.536, 0.466, 0.382, 0.308, 0.227, 0.146},  {0.406, 0.395, 0.383, 0.372, 0.360, 0.352, 0.349, 0.335, 0.340, 0.720, 0.816, 0.824, 0.818, 0.807, 0.792, 0.773, 0.748, 0.727, 0.702, 0.673, 0.640, 0.606, 0.568, 0.519, 0.467, 0.403, 0.341, 0.280, 0.217}},
      {{0.1864, 0.2229, 0.2683, 0.3215, 0.3781, 0.4395, 0.4987, 0.5316, 0.5618, 0.6374, 0.7315, 0.8336, 0.935, 1.0423, 1.1437, 1.2443, 1.263, 1.3573, 1.4604, 1.5595, 1.652, 1.7467, 1.8533, 1.9495, 2.045, 2.1238, 2.225, 2.3202, 2.4159},  {0.346, 0.331, 0.314, 0.297, 0.281, 0.267, 0.257, 0.250, 0.228, 0.276, 1.188, 1.203, 1.199, 1.185, 1.161, 1.134, 1.117, 1.087, 1.050, 1.012, 0.972, 0.926, 0.867, 0.804, 0.729, 0.612, 0.523, 0.401, 0.271},  {0.363, 0.350, 0.335, 0.319, 0.305, 0.292, 0.282, 0.275, 0.256, 0.288, 0.992, 1.005, 1.003, 0.992, 0.977, 0.956, 0.942, 0.921, 0.894, 0.866, 0.836, 0.802, 0.758, 0.710, 0.653, 0.560, 0.490, 0.392, 0.286},  {0.411, 0.400, 0.387, 0.374, 0.362, 0.350, 0.340, 0.333, 0.317, 0.329, 0.828, 0.840, 0.839, 0.833, 0.822, 0.809, 0.799, 0.785, 0.767, 0.748, 0.728, 0.705, 0.674, 0.641, 0.600, 0.527, 0.476, 0.402, 0.322},  {0.346, 0.331, 0.313, 0.296, 0.280, 0.266, 0.255, 0.248, 0.226, 0.274, 1.138, 1.145, 1.131, 1.106, 1.072, 1.033, 1.014, 0.973, 0.923, 0.870, 0.817, 0.756, 0.682, 0.606, 0.521, 0.410, 0.316, 0.211, 0.102},  {0.363, 0.349, 0.334, 0.319, 0.304, 0.291, 0.281, 0.274, 0.255, 0.287, 0.956, 0.962, 0.953, 0.935, 0.911, 0.881, 0.867, 0.837, 0.800, 0.761, 0.721, 0.676, 0.619, 0.561, 0.494, 0.403, 0.327, 0.239, 0.142},  {0.411, 0.400, 0.387, 0.374, 0.361, 0.349, 0.340, 0.333, 0.316, 0.328, 0.804, 0.811, 0.806, 0.795, 0.779, 0.759, 0.749, 0.729, 0.704, 0.678, 0.651, 0.620, 0.580, 0.538, 0.489, 0.417, 0.359, 0.291, 0.218}},
      {{0.2061, 0.2508, 0.2984, 0.3476, 0.407, 0.473, 0.5416, 0.5813, 0.6245, 0.7291, 0.8127, 0.9095, 1.0031, 1.0966, 1.1905, 1.2907, 1.3835, 1.4606, 1.5683, 1.6733, 1.7823, 1.8227, 1.9215, 2.0444, 2.1254, 2.2266, 2.3272, 2.4302},  {0.354, 0.335, 0.318, 0.301, 0.284, 0.269, 0.256, 0.247, 0.225, 0.260, 1.165, 1.183, 1.185, 1.176, 1.160, 1.135, 1.109, 1.072, 1.032, 0.988, 0.935, 0.835, 0.787, 0.718, 0.666, 0.579, 0.449, 0.296},  {0.369, 0.353, 0.338, 0.323, 0.308, 0.294, 0.281, 0.273, 0.253, 0.276, 0.976, 0.991, 0.992, 0.987, 0.975, 0.957, 0.937, 0.910, 0.880, 0.848, 0.810, 0.730, 0.695, 0.643, 0.603, 0.536, 0.432, 0.307},  {0.416, 0.403, 0.389, 0.377, 0.364, 0.351, 0.339, 0.332, 0.314, 0.321, 0.817, 0.829, 0.831, 0.829, 0.822, 0.810, 0.797, 0.778, 0.759, 0.737, 0.711, 0.650, 0.626, 0.590, 0.561, 0.512, 0.433, 0.338},  {0.353, 0.335, 0.317, 0.301, 0.284, 0.268, 0.254, 0.245, 0.224, 0.259, 1.112, 1.121, 1.112, 1.095, 1.068, 1.032, 0.992, 0.946, 0.892, 0.833, 0.765, 0.670, 0.608, 0.522, 0.458, 0.363, 0.244, 0.117},  {0.369, 0.353, 0.337, 0.323, 0.307, 0.293, 0.280, 0.272, 0.252, 0.274, 0.936, 0.945, 0.939, 0.927, 0.907, 0.881, 0.852, 0.818, 0.777, 0.733, 0.682, 0.607, 0.560, 0.493, 0.443, 0.367, 0.268, 0.153},  {0.416, 0.403, 0.389, 0.377, 0.363, 0.350, 0.339, 0.331, 0.313, 0.320, 0.791, 0.799, 0.797, 0.789, 0.777, 0.759, 0.741, 0.717, 0.690, 0.660, 0.625, 0.566, 0.533, 0.486, 0.449, 0.392, 0.315, 0.230}},
      {{0.2102, 0.2525, 0.3083, 0.3581, 0.418, 0.4792, 0.5492, 0.5886, 0.638, 0.747, 0.8297, 0.9266, 1.026, 1.1197, 1.2209, 1.3146, 1.4107, 1.4482, 1.47, 1.578, 1.6751, 1.7777, 1.8748, 1.9645, 2.0508, 2.1656, 2.2541, 2.3435, 2.4591},  {0.355, 0.338, 0.317, 0.300, 0.284, 0.269, 0.256, 0.247, 0.225, 0.260, 1.161, 1.181, 1.181, 1.172, 1.155, 1.131, 1.103, 0.977, 0.969, 0.937, 0.903, 0.861, 0.814, 0.763, 0.708, 0.620, 0.527, 0.402, 0.221},  {0.370, 0.355, 0.337, 0.322, 0.307, 0.294, 0.281, 0.273, 0.253, 0.275, 0.973, 0.988, 0.989, 0.984, 0.971, 0.954, 0.933, 0.835, 0.830, 0.807, 0.782, 0.751, 0.716, 0.678, 0.636, 0.568, 0.495, 0.395, 0.244},  {0.417, 0.404, 0.388, 0.376, 0.363, 0.351, 0.339, 0.332, 0.313, 0.320, 0.814, 0.827, 0.830, 0.827, 0.819, 0.808, 0.794, 0.721, 0.718, 0.703, 0.686, 0.666, 0.643, 0.616, 0.586, 0.537, 0.482, 0.406, 0.290},  {0.355, 0.337, 0.316, 0.300, 0.283, 0.268, 0.255, 0.246, 0.223, 0.258, 1.106, 1.116, 1.107, 1.089, 1.059, 1.025, 0.984, 0.867, 0.856, 0.811, 0.764, 0.708, 0.648, 0.586, 0.520, 0.420, 0.327, 0.220, 0.024},  {0.370, 0.355, 0.337, 0.322, 0.307, 0.293, 0.280, 0.272, 0.252, 0.274, 0.933, 0.942, 0.935, 0.923, 0.901, 0.876, 0.846, 0.754, 0.747, 0.713, 0.678, 0.637, 0.591, 0.544, 0.493, 0.414, 0.338, 0.247, 0.059},  {0.416, 0.404, 0.388, 0.376, 0.362, 0.351, 0.339, 0.331, 0.313, 0.319, 0.788, 0.796, 0.794, 0.786, 0.773, 0.756, 0.736, 0.667, 0.662, 0.639, 0.617, 0.588, 0.557, 0.524, 0.486, 0.428, 0.369, 0.300, 0.150}},
      {{0.2565, 0.3089, 0.3608, 0.4107, 0.4719, 0.5383, 0.611, 0.6533, 0.7003, 0.8159, 0.9302, 1.0234, 1.1269, 1.2307, 1.3342, 1.4388, 1.5429, 1.5919, 1.6698, 1.7978, 1.8888, 1.9812, 2.0912, 2.1877, 2.2821, 2.3811, 2.4758, 2.5726},  {0.369, 0.347, 0.327, 0.309, 0.291, 0.275, 0.260, 0.250, 0.229, 0.217, 0.251, 1.109, 1.131, 1.134, 1.122, 1.100, 1.070, 0.981, 0.956, 0.907, 0.867, 0.820, 0.757, 0.688, 0.592, 0.448, 0.283, 0.052},  {0.381, 0.362, 0.345, 0.330, 0.314, 0.299, 0.285, 0.276, 0.256, 0.244, 0.268, 0.934, 0.951, 0.954, 0.947, 0.931, 0.908, 0.841, 0.822, 0.787, 0.756, 0.721, 0.674, 0.622, 0.547, 0.432, 0.296, 0.067},  {0.425, 0.409, 0.394, 0.381, 0.367, 0.354, 0.342, 0.333, 0.316, 0.303, 0.313, 0.787, 0.802, 0.806, 0.803, 0.793, 0.778, 0.727, 0.716, 0.692, 0.672, 0.649, 0.616, 0.578, 0.524, 0.436, 0.328, 0.086},  {0.368, 0.347, 0.326, 0.309, 0.291, 0.275, 0.260, 0.250, 0.228, 0.216, 0.249, 1.046, 1.057, 1.048, 1.026, 0.990, 0.946, 0.858, 0.822, 0.756, 0.688, 0.641, 0.563, 0.481, 0.381, 0.251, 0.023, 0.481},  {0.381, 0.362, 0.345, 0.330, 0.313, 0.299, 0.285, 0.275, 0.256, 0.243, 0.266, 0.887, 0.897, 0.891, 0.876, 0.850, 0.818, 0.749, 0.723, 0.674, 0.623, 0.587, 0.527, 0.463, 0.383, 0.274, 0.060, 0.060},  {0.425, 0.409, 0.394, 0.381, 0.367, 0.354, 0.341, 0.333, 0.315, 0.302, 0.311, 0.756, 0.766, 0.765, 0.755, 0.739, 0.718, 0.666, 0.649, 0.616, 0.582, 0.556, 0.514, 0.468, 0.406, 0.321, 0.157, 0.157}},
      {{0.3185, 0.3701, 0.4211, 0.4749, 0.5368, 0.6028, 0.6745, 0.7197, 0.7518, 0.854, 0.9536, 1.0493, 1.1528, 1.2552, 1.3566, 1.4598, 1.565, 1.6659, 1.7359, 1.8392, 1.9346, 2.0384, 2.1341, 2.2172, 2.3138, 2.4199},  {0.386, 0.363, 0.341, 0.322, 0.302, 0.285, 0.269, 0.258, 0.239, 0.225, 0.218, 0.216, 0.252, 1.036, 1.075, 1.080, 1.075, 1.056, 0.961, 0.923, 0.884, 0.837, 0.787, 0.733, 0.645, 0.506},  {0.395, 0.376, 0.357, 0.340, 0.323, 0.307, 0.292, 0.282, 0.265, 0.251, 0.244, 0.241, 0.270, 0.878, 0.909, 0.914, 0.912, 0.898, 0.826, 0.798, 0.769, 0.735, 0.697, 0.656, 0.589, 0.479},  {0.436, 0.419, 0.403, 0.389, 0.374, 0.360, 0.347, 0.338, 0.322, 0.310, 0.301, 0.296, 0.312, 0.746, 0.772, 0.778, 0.779, 0.771, 0.718, 0.700, 0.681, 0.657, 0.632, 0.603, 0.556, 0.472},  {0.386, 0.362, 0.341, 0.321, 0.302, 0.285, 0.269, 0.258, 0.239, 0.224, 0.216, 0.213, 0.248, 0.959, 0.983, 0.975, 0.956, 0.922, 0.828, 0.776, 0.724, 0.660, 0.595, 0.529, 0.423, 0.295},  {0.395, 0.375, 0.357, 0.340, 0.322, 0.307, 0.292, 0.282, 0.265, 0.251, 0.243, 0.239, 0.266, 0.820, 0.842, 0.837, 0.824, 0.799, 0.727, 0.689, 0.649, 0.602, 0.552, 0.501, 0.418, 0.312},  {0.436, 0.419, 0.403, 0.388, 0.374, 0.360, 0.347, 0.338, 0.322, 0.309, 0.301, 0.295, 0.310, 0.708, 0.727, 0.726, 0.720, 0.705, 0.652, 0.626, 0.600, 0.567, 0.532, 0.496, 0.434, 0.351}},
      {{0.3767, 0.4256, 0.4783, 0.5342, 0.5945, 0.6617, 0.7348, 0.7745, 0.8044, 0.8996, 0.9923, 1.0906, 1.1855, 1.2812, 1.3658, 1.376, 1.4653, 1.5586, 1.6531, 1.7419, 1.8342, 1.8646, 1.9616, 2.0593, 2.1566, 2.2498, 2.3536, 2.4389, 2.541},  {0.403, 0.378, 0.355, 0.333, 0.313, 0.294, 0.276, 0.265, 0.248, 0.232, 0.223, 0.217, 0.216, 0.217, 0.262, 0.626, 0.937, 0.998, 1.007, 1.005, 0.990, 0.933, 0.901, 0.862, 0.814, 0.757, 0.669, 0.580, 0.424},  {0.409, 0.388, 0.368, 0.349, 0.332, 0.315, 0.299, 0.289, 0.273, 0.258, 0.249, 0.243, 0.240, 0.241, 0.279, 0.551, 0.802, 0.851, 0.859, 0.859, 0.850, 0.804, 0.781, 0.752, 0.717, 0.674, 0.607, 0.538, 0.412},  {0.446, 0.429, 0.412, 0.395, 0.380, 0.366, 0.351, 0.343, 0.328, 0.315, 0.306, 0.299, 0.294, 0.293, 0.321, 0.488, 0.691, 0.730, 0.739, 0.741, 0.737, 0.703, 0.687, 0.668, 0.645, 0.615, 0.568, 0.517, 0.418},  {0.403, 0.378, 0.355, 0.333, 0.313, 0.294, 0.276, 0.265, 0.247, 0.232, 0.222, 0.216, 0.214, 0.213, 0.257, 0.566, 0.852, 0.897, 0.892, 0.877, 0.848, 0.790, 0.743, 0.688, 0.625, 0.555, 0.459, 0.368, 0.368},  {0.409, 0.388, 0.368, 0.349, 0.331, 0.314, 0.299, 0.289, 0.272, 0.258, 0.248, 0.242, 0.239, 0.238, 0.275, 0.504, 0.738, 0.776, 0.775, 0.765, 0.744, 0.698, 0.664, 0.622, 0.575, 0.521, 0.446, 0.372, 0.372},  {0.446, 0.428, 0.411, 0.395, 0.380, 0.365, 0.351, 0.343, 0.328, 0.315, 0.306, 0.298, 0.293, 0.291, 0.318, 0.453, 0.648, 0.680, 0.682, 0.678, 0.666, 0.631, 0.608, 0.580, 0.547, 0.509, 0.455, 0.399, 0.399}},
      {{0.4348, 0.483, 0.5369, 0.5866, 0.6466, 0.7179, 0.7843, 0.8274, 0.8701, 0.9734, 1.0775, 1.1835, 1.2879, 1.396, 1.4978, 1.605, 1.7107, 1.814, 1.9138, 2.0197, 2.119, 2.2207, 2.3247, 2.429},  {0.421, 0.395, 0.369, 0.348, 0.326, 0.305, 0.287, 0.275, 0.254, 0.236, 0.225, 0.217, 0.214, 0.213, 0.212, 0.546, 0.879, 0.925, 0.929, 0.922, 0.828, 0.773, 0.700, 0.617},  {0.424, 0.402, 0.380, 0.362, 0.343, 0.324, 0.308, 0.297, 0.278, 0.262, 0.250, 0.243, 0.239, 0.238, 0.238, 0.487, 0.758, 0.796, 0.802, 0.797, 0.726, 0.685, 0.630, 0.566},  {0.457, 0.438, 0.420, 0.405, 0.388, 0.372, 0.359, 0.349, 0.332, 0.317, 0.307, 0.298, 0.293, 0.290, 0.289, 0.437, 0.660, 0.693, 0.700, 0.701, 0.648, 0.620, 0.582, 0.536},  {0.421, 0.395, 0.369, 0.348, 0.326, 0.304, 0.287, 0.275, 0.254, 0.236, 0.224, 0.217, 0.213, 0.210, 0.207, 0.479, 0.778, 0.806, 0.793, 0.766, 0.664, 0.597, 0.514, 0.423},  {0.424, 0.401, 0.379, 0.362, 0.343, 0.323, 0.308, 0.297, 0.278, 0.262, 0.250, 0.242, 0.238, 0.235, 0.233, 0.434, 0.683, 0.708, 0.700, 0.682, 0.602, 0.552, 0.488, 0.416},  {0.457, 0.438, 0.420, 0.404, 0.388, 0.372, 0.358, 0.349, 0.332, 0.317, 0.306, 0.298, 0.292, 0.288, 0.286, 0.397, 0.609, 0.633, 0.632, 0.622, 0.564, 0.529, 0.483, 0.430}},
      {{0.4941, 0.5384, 0.591, 0.647, 0.7067, 0.77, 0.8451, 0.8839, 0.9248, 1.0243, 1.1241, 1.2245, 1.3242, 1.4274, 1.5265, 1.6324, 1.7313, 1.8336, 1.9279, 2.0251, 2.1158, 2.2121, 2.2799, 2.3794, 2.4768, 2.5743, 2.5838, 2.5971, 2.6244, 2.6721},  {0.442, 0.415, 0.388, 0.362, 0.339, 0.318, 0.297, 0.284, 0.265, 0.245, 0.231, 0.221, 0.215, 0.211, 0.208, 0.205, 0.201, 0.524, 0.810, 0.845, 0.847, 0.830, 0.737, 0.682, 0.608, 0.511, 0.509, 0.496, 0.448, 0.319},  {0.441, 0.418, 0.395, 0.373, 0.353, 0.335, 0.316, 0.305, 0.287, 0.270, 0.257, 0.247, 0.240, 0.236, 0.234, 0.232, 0.229, 0.470, 0.706, 0.736, 0.740, 0.729, 0.656, 0.614, 0.558, 0.481, 0.481, 0.470, 0.431, 0.323},  {0.469, 0.450, 0.431, 0.413, 0.396, 0.381, 0.364, 0.355, 0.339, 0.324, 0.312, 0.303, 0.295, 0.290, 0.286, 0.284, 0.283, 0.426, 0.623, 0.650, 0.656, 0.652, 0.597, 0.568, 0.528, 0.471, 0.471, 0.463, 0.433, 0.345},  {0.442, 0.415, 0.388, 0.362, 0.339, 0.318, 0.296, 0.284, 0.264, 0.245, 0.231, 0.221, 0.214, 0.210, 0.206, 0.201, 0.193, 0.447, 0.694, 0.711, 0.697, 0.664, 0.570, 0.502, 0.420, 0.324, 0.324, 0.310, 0.265, 0.162},  {0.440, 0.418, 0.395, 0.373, 0.353, 0.335, 0.316, 0.305, 0.287, 0.269, 0.256, 0.247, 0.240, 0.235, 0.232, 0.228, 0.223, 0.409, 0.619, 0.635, 0.627, 0.604, 0.529, 0.477, 0.412, 0.335, 0.335, 0.323, 0.285, 0.194},  {0.469, 0.450, 0.431, 0.413, 0.396, 0.381, 0.364, 0.355, 0.339, 0.324, 0.312, 0.302, 0.295, 0.289, 0.285, 0.282, 0.279, 0.381, 0.563, 0.581, 0.580, 0.567, 0.509, 0.473, 0.426, 0.367, 0.367, 0.358, 0.327, 0.254}},
      {{0.5543, 0.6002, 0.6495, 0.7007, 0.7596, 0.8307, 0.9062, 0.945, 1.0021, 1.0933, 1.1798, 1.2815, 1.379, 1.4779, 1.57, 1.6655, 1.7667, 1.8709, 1.9713, 2.0648, 2.1493, 2.2414, 2.3271, 2.415, 2.4651, 2.5699, 2.6747},  {0.466, 0.436, 0.407, 0.382, 0.357, 0.331, 0.307, 0.294, 0.272, 0.253, 0.237, 0.224, 0.215, 0.207, 0.202, 0.198, 0.193, 0.187, 0.180, 0.508, 0.726, 0.745, 0.741, 0.725, 0.641, 0.555, 0.413},  {0.460, 0.434, 0.411, 0.389, 0.368, 0.346, 0.325, 0.313, 0.293, 0.276, 0.263, 0.250, 0.241, 0.233, 0.228, 0.224, 0.221, 0.216, 0.211, 0.459, 0.642, 0.659, 0.659, 0.649, 0.582, 0.515, 0.402},  {0.483, 0.462, 0.443, 0.425, 0.407, 0.388, 0.371, 0.361, 0.343, 0.328, 0.316, 0.305, 0.296, 0.288, 0.283, 0.278, 0.275, 0.272, 0.269, 0.421, 0.577, 0.595, 0.598, 0.595, 0.543, 0.495, 0.409},  {0.466, 0.435, 0.407, 0.382, 0.357, 0.331, 0.307, 0.293, 0.271, 0.252, 0.237, 0.224, 0.214, 0.207, 0.201, 0.196, 0.189, 0.180, 0.169, 0.417, 0.600, 0.600, 0.581, 0.546, 0.463, 0.371, 0.245},  {0.460, 0.434, 0.411, 0.389, 0.368, 0.346, 0.325, 0.313, 0.293, 0.276, 0.262, 0.250, 0.240, 0.233, 0.227, 0.223, 0.218, 0.211, 0.202, 0.388, 0.546, 0.550, 0.537, 0.512, 0.445, 0.372, 0.268},  {0.483, 0.462, 0.442, 0.425, 0.407, 0.388, 0.371, 0.361, 0.342, 0.328, 0.316, 0.305, 0.296, 0.288, 0.282, 0.277, 0.272, 0.267, 0.262, 0.368, 0.510, 0.519, 0.514, 0.500, 0.446, 0.393, 0.314}},
      {{0.6098, 0.6536, 0.7051, 0.7616, 0.821, 0.8847, 0.961, 1.0081, 1.0592, 1.1435, 1.248, 1.3425, 1.4449, 1.5472, 1.6379, 1.7372, 1.8431, 1.9394, 2.0354, 2.1389, 2.2337, 2.3265, 2.4095, 2.4984, 2.5874, 2.6803, 2.7853},  {0.491, 0.459, 0.427, 0.397, 0.370, 0.345, 0.319, 0.301, 0.281, 0.263, 0.243, 0.230, 0.218, 0.208, 0.202, 0.196, 0.190, 0.184, 0.177, 0.170, 0.482, 0.680, 0.674, 0.643, 0.602, 0.452, 0.293},  {0.480, 0.454, 0.427, 0.402, 0.379, 0.357, 0.335, 0.319, 0.301, 0.285, 0.268, 0.255, 0.244, 0.235, 0.228, 0.222, 0.217, 0.213, 0.208, 0.203, 0.440, 0.607, 0.605, 0.582, 0.552, 0.433, 0.305},  {0.498, 0.476, 0.454, 0.434, 0.415, 0.397, 0.378, 0.365, 0.348, 0.335, 0.320, 0.309, 0.299, 0.290, 0.283, 0.277, 0.272, 0.268, 0.265, 0.262, 0.408, 0.552, 0.555, 0.541, 0.523, 0.432, 0.336},  {0.490, 0.459, 0.427, 0.397, 0.370, 0.344, 0.318, 0.300, 0.280, 0.262, 0.242, 0.228, 0.216, 0.207, 0.200, 0.194, 0.187, 0.179, 0.169, 0.157, 0.381, 0.539, 0.522, 0.481, 0.421, 0.287, 0.156},  {0.480, 0.453, 0.427, 0.401, 0.378, 0.357, 0.334, 0.319, 0.301, 0.285, 0.267, 0.254, 0.243, 0.234, 0.227, 0.221, 0.215, 0.209, 0.201, 0.192, 0.360, 0.499, 0.488, 0.458, 0.413, 0.303, 0.192},  {0.498, 0.476, 0.454, 0.433, 0.414, 0.397, 0.378, 0.365, 0.348, 0.334, 0.320, 0.308, 0.298, 0.289, 0.282, 0.276, 0.270, 0.265, 0.260, 0.254, 0.348, 0.478, 0.474, 0.455, 0.425, 0.339, 0.257}},
      {{0.6744, 0.7197, 0.7696, 0.8232, 0.886, 0.9574, 1.0331, 1.0896, 1.1286, 1.2266, 1.3272, 1.4232, 1.5288, 1.625, 1.7269, 1.8239, 1.9293, 2.0288, 2.1165, 2.2179, 2.3212, 2.4224, 2.5305, 2.6371, 2.7467, 2.7916, 2.8864},  {0.525, 0.488, 0.454, 0.422, 0.389, 0.358, 0.329, 0.306, 0.289, 0.266, 0.246, 0.230, 0.216, 0.206, 0.196, 0.189, 0.181, 0.175, 0.169, 0.162, 0.155, 0.408, 0.584, 0.525, 0.434, 0.368, 0.201},  {0.507, 0.477, 0.448, 0.421, 0.394, 0.368, 0.343, 0.324, 0.308, 0.288, 0.270, 0.256, 0.243, 0.233, 0.224, 0.217, 0.210, 0.204, 0.200, 0.195, 0.190, 0.381, 0.533, 0.489, 0.419, 0.365, 0.227},  {0.518, 0.494, 0.470, 0.448, 0.426, 0.404, 0.384, 0.367, 0.352, 0.335, 0.321, 0.309, 0.297, 0.288, 0.279, 0.272, 0.265, 0.260, 0.256, 0.253, 0.250, 0.363, 0.498, 0.470, 0.423, 0.381, 0.280},  {0.525, 0.488, 0.453, 0.421, 0.388, 0.357, 0.328, 0.304, 0.287, 0.263, 0.243, 0.228, 0.214, 0.203, 0.194, 0.186, 0.179, 0.171, 0.163, 0.152, 0.140, 0.309, 0.436, 0.384, 0.286, 0.220, 0.097},  {0.507, 0.477, 0.448, 0.421, 0.394, 0.367, 0.342, 0.322, 0.306, 0.286, 0.268, 0.255, 0.242, 0.231, 0.222, 0.215, 0.208, 0.201, 0.195, 0.186, 0.176, 0.302, 0.418, 0.379, 0.302, 0.246, 0.138},  {0.518, 0.493, 0.470, 0.448, 0.425, 0.403, 0.383, 0.366, 0.351, 0.334, 0.319, 0.308, 0.296, 0.287, 0.278, 0.271, 0.264, 0.257, 0.252, 0.246, 0.240, 0.303, 0.417, 0.392, 0.339, 0.296, 0.218}},
      {{0.7249, 0.7702, 0.8212, 0.8771, 0.9439, 1.0217, 1.1114, 1.1683, 1.2003, 1.3054, 1.4086, 1.5074, 1.6067, 1.713, 1.8162, 1.9163, 2.0197, 2.1185, 2.2206, 2.3245, 2.4243, 2.5264, 2.6285, 2.7105, 2.794, 2.8816, 2.9489, 3.0202, 3.026},  {0.557, 0.517, 0.477, 0.440, 0.401, 0.364, 0.327, 0.301, 0.284, 0.258, 0.237, 0.219, 0.203, 0.189, 0.176, 0.166, 0.156, 0.148, 0.141, 0.133, 0.126, 0.119, 0.268, 0.377, 0.342, 0.240, 0.144, 0.071, 0.073},  {0.533, 0.499, 0.467, 0.436, 0.404, 0.372, 0.341, 0.319, 0.303, 0.281, 0.262, 0.246, 0.232, 0.218, 0.207, 0.197, 0.187, 0.179, 0.172, 0.166, 0.160, 0.154, 0.265, 0.365, 0.340, 0.259, 0.179, 0.115, 0.118},  {0.537, 0.510, 0.483, 0.457, 0.431, 0.405, 0.378, 0.360, 0.345, 0.327, 0.312, 0.299, 0.286, 0.274, 0.264, 0.255, 0.246, 0.238, 0.231, 0.224, 0.219, 0.215, 0.269, 0.367, 0.354, 0.299, 0.243, 0.201, 0.203},  {0.556, 0.515, 0.476, 0.438, 0.399, 0.361, 0.324, 0.298, 0.280, 0.254, 0.232, 0.215, 0.199, 0.185, 0.173, 0.162, 0.153, 0.144, 0.134, 0.124, 0.112, 0.100, 0.187, 0.261, 0.236, 0.135, 0.071, 0.030, 0.031},  {0.532, 0.499, 0.466, 0.434, 0.402, 0.370, 0.338, 0.316, 0.299, 0.277, 0.258, 0.242, 0.228, 0.215, 0.203, 0.193, 0.184, 0.175, 0.167, 0.157, 0.147, 0.136, 0.197, 0.272, 0.255, 0.171, 0.112, 0.068, 0.069},  {0.537, 0.509, 0.482, 0.457, 0.430, 0.404, 0.377, 0.358, 0.342, 0.324, 0.309, 0.296, 0.284, 0.272, 0.262, 0.252, 0.243, 0.235, 0.227, 0.218, 0.209, 0.201, 0.217, 0.298, 0.293, 0.234, 0.190, 0.156, 0.158}},
      {{0.779, 0.8239, 0.8777, 0.936, 1.0066, 1.094, 1.1909, 1.2612, 1.2973, 1.3991, 1.5013, 1.6029, 1.7035, 1.8074, 1.9068, 2.0087, 2.1065, 2.2163, 2.3061, 2.4227, 2.5158, 2.6188, 2.7163, 2.8192, 2.9088, 2.9939, 3.0828},  {0.599, 0.553, 0.505, 0.461, 0.415, 0.367, 0.323, 0.291, 0.269, 0.243, 0.221, 0.201, 0.184, 0.168, 0.155, 0.142, 0.131, 0.119, 0.110, 0.101, 0.094, 0.088, 0.083, 0.143, 0.213, 0.170, 0.075},  {0.567, 0.528, 0.489, 0.452, 0.414, 0.374, 0.336, 0.308, 0.288, 0.266, 0.246, 0.229, 0.214, 0.199, 0.186, 0.174, 0.163, 0.152, 0.143, 0.133, 0.126, 0.120, 0.116, 0.155, 0.229, 0.199, 0.117},  {0.563, 0.531, 0.498, 0.468, 0.437, 0.405, 0.374, 0.350, 0.329, 0.311, 0.296, 0.282, 0.269, 0.257, 0.246, 0.236, 0.226, 0.215, 0.207, 0.198, 0.191, 0.184, 0.178, 0.187, 0.259, 0.251, 0.194},  {0.598, 0.551, 0.503, 0.457, 0.410, 0.362, 0.317, 0.283, 0.261, 0.235, 0.213, 0.194, 0.178, 0.162, 0.149, 0.136, 0.125, 0.113, 0.104, 0.092, 0.082, 0.071, 0.062, 0.087, 0.134, 0.112, 0.025},  {0.565, 0.527, 0.487, 0.450, 0.410, 0.369, 0.331, 0.302, 0.281, 0.259, 0.240, 0.223, 0.208, 0.194, 0.182, 0.169, 0.158, 0.147, 0.137, 0.125, 0.115, 0.104, 0.095, 0.106, 0.161, 0.148, 0.059},  {0.562, 0.529, 0.497, 0.466, 0.435, 0.402, 0.370, 0.346, 0.324, 0.306, 0.291, 0.278, 0.266, 0.254, 0.243, 0.232, 0.223, 0.212, 0.203, 0.192, 0.183, 0.173, 0.164, 0.164, 0.207, 0.212, 0.139}},
      {{0.8329, 0.8774, 0.9302, 0.9949, 1.0673, 1.1574, 1.2675, 1.3554, 1.3727, 1.4708, 1.5744, 1.6662, 1.7714, 1.8734, 1.9704, 2.0724, 2.17, 2.2604, 2.3655, 2.4746, 2.5585, 2.6581, 2.7628, 2.859, 2.96, 3.05, 3.1447, 3.237},  {0.653, 0.596, 0.541, 0.482, 0.428, 0.372, 0.317, 0.275, 0.255, 0.229, 0.205, 0.187, 0.168, 0.152, 0.138, 0.125, 0.113, 0.103, 0.092, 0.082, 0.076, 0.069, 0.063, 0.059, 0.064, 0.108, 0.068, 0.041},  {0.609, 0.561, 0.516, 0.469, 0.424, 0.377, 0.331, 0.294, 0.273, 0.251, 0.231, 0.215, 0.198, 0.184, 0.171, 0.158, 0.146, 0.136, 0.126, 0.116, 0.109, 0.102, 0.095, 0.089, 0.094, 0.137, 0.107, 0.081},  {0.595, 0.562, 0.527, 0.491, 0.457, 0.420, 0.382, 0.348, 0.314, 0.296, 0.280, 0.267, 0.254, 0.243, 0.232, 0.222, 0.213, 0.205, 0.197, 0.189, 0.183, 0.176, 0.169, 0.164, 0.158, 0.183, 0.179, 0.165},  {0.650, 0.592, 0.535, 0.475, 0.419, 0.361, 0.305, 0.262, 0.240, 0.215, 0.193, 0.176, 0.159, 0.143, 0.130, 0.117, 0.105, 0.095, 0.084, 0.073, 0.065, 0.056, 0.045, 0.038, 0.036, 0.061, 0.021, 0.011},  {0.606, 0.558, 0.512, 0.463, 0.417, 0.369, 0.321, 0.284, 0.261, 0.240, 0.221, 0.206, 0.190, 0.176, 0.164, 0.151, 0.140, 0.130, 0.119, 0.108, 0.099, 0.088, 0.077, 0.067, 0.065, 0.092, 0.051, 0.037},  {0.593, 0.559, 0.525, 0.488, 0.451, 0.411, 0.371, 0.339, 0.305, 0.288, 0.273, 0.261, 0.248, 0.237, 0.227, 0.218, 0.209, 0.201, 0.192, 0.183, 0.175, 0.166, 0.156, 0.147, 0.139, 0.144, 0.123, 0.117}},
      {{0.8929, 0.9375, 0.9895, 1.0552, 1.1378, 1.2309, 1.3559, 1.4615, 1.4757, 1.5758, 1.6801, 1.7799, 1.8742, 1.9722, 2.0736, 2.1782, 2.2748, 2.3743, 2.478, 2.5789, 2.6792, 2.7749, 2.8788, 2.9741, 3.073, 3.1565, 3.2419, 3.3267},  {0.728, 0.660, 0.590, 0.516, 0.443, 0.378, 0.310, 0.256, 0.230, 0.204, 0.182, 0.163, 0.147, 0.133, 0.119, 0.106, 0.095, 0.085, 0.076, 0.068, 0.061, 0.055, 0.049, 0.045, 0.049, 0.067, 0.053, 0.034},  {0.668, 0.617, 0.564, 0.507, 0.448, 0.391, 0.327, 0.279, 0.251, 0.229, 0.209, 0.192, 0.178, 0.165, 0.152, 0.140, 0.130, 0.120, 0.111, 0.102, 0.095, 0.088, 0.082, 0.076, 0.075, 0.097, 0.092, 0.072},  {0.639, 0.605, 0.569, 0.531, 0.488, 0.444, 0.395, 0.356, 0.296, 0.279, 0.264, 0.250, 0.239, 0.229, 0.219, 0.210, 0.201, 0.194, 0.186, 0.179, 0.173, 0.167, 0.161, 0.156, 0.151, 0.147, 0.164, 0.153},  {0.722, 0.651, 0.581, 0.503, 0.427, 0.359, 0.290, 0.236, 0.208, 0.184, 0.164, 0.148, 0.134, 0.121, 0.108, 0.096, 0.086, 0.076, 0.066, 0.058, 0.049, 0.040, 0.031, 0.025, 0.022, 0.033, 0.014, 0.010},  {0.664, 0.611, 0.557, 0.496, 0.431, 0.372, 0.311, 0.263, 0.233, 0.212, 0.194, 0.179, 0.166, 0.154, 0.143, 0.131, 0.121, 0.111, 0.102, 0.092, 0.082, 0.073, 0.062, 0.053, 0.046, 0.060, 0.041, 0.035},  {0.636, 0.601, 0.561, 0.514, 0.465, 0.421, 0.373, 0.336, 0.286, 0.269, 0.254, 0.242, 0.232, 0.222, 0.212, 0.203, 0.195, 0.187, 0.179, 0.171, 0.163, 0.154, 0.144, 0.136, 0.127, 0.126, 0.110, 0.112}},
      {{0.9471, 0.9955, 1.0498, 1.1179, 1.2027, 1.3139, 1.3822, 1.4588, 1.5268, 1.5976, 1.6047, 1.703, 1.7972, 1.8952, 1.9864, 2.0918, 2.1888, 2.2771, 2.377, 2.4785, 2.5704, 2.678, 2.7651, 2.8719, 2.9565, 3.0585, 3.1547, 3.2663, 3.3783},  {0.818, 0.740, 0.663, 0.580, 0.491, 0.394, 0.346, 0.300, 0.264, 0.224, 0.189, 0.167, 0.149, 0.133, 0.120, 0.106, 0.095, 0.086, 0.076, 0.068, 0.061, 0.053, 0.049, 0.043, 0.039, 0.035, 0.033, 0.045, 0.032},  {0.738, 0.680, 0.622, 0.559, 0.488, 0.410, 0.370, 0.332, 0.303, 0.270, 0.216, 0.197, 0.181, 0.167, 0.155, 0.142, 0.132, 0.123, 0.114, 0.105, 0.098, 0.090, 0.085, 0.079, 0.074, 0.069, 0.064, 0.072, 0.070},  {0.691, 0.652, 0.613, 0.566, 0.512, 0.451, 0.421, 0.393, 0.370, 0.347, 0.284, 0.267, 0.253, 0.240, 0.230, 0.218, 0.209, 0.201, 0.193, 0.185, 0.179, 0.173, 0.166, 0.160, 0.156, 0.151, 0.146, 0.150, 0.150},  {0.806, 0.724, 0.635, 0.537, 0.444, 0.354, 0.309, 0.267, 0.234, 0.198, 0.162, 0.143, 0.128, 0.114, 0.103, 0.092, 0.082, 0.073, 0.064, 0.056, 0.049, 0.042, 0.035, 0.028, 0.022, 0.017, 0.012, 0.016, 0.009},  {0.729, 0.666, 0.597, 0.519, 0.445, 0.371, 0.335, 0.300, 0.273, 0.242, 0.196, 0.179, 0.164, 0.152, 0.141, 0.130, 0.120, 0.111, 0.102, 0.093, 0.086, 0.078, 0.070, 0.060, 0.053, 0.045, 0.037, 0.037, 0.032},  {0.684, 0.637, 0.583, 0.524, 0.467, 0.412, 0.384, 0.358, 0.338, 0.317, 0.269, 0.253, 0.240, 0.228, 0.218, 0.208, 0.199, 0.191, 0.183, 0.175, 0.168, 0.160, 0.152, 0.143, 0.135, 0.126, 0.117, 0.113, 0.107}},
      {{1.0029, 1.0499, 1.1095, 1.1814, 1.2743, 1.3952, 1.4725, 1.5544, 1.6361, 1.7219, 1.7321, 1.8285, 1.9297, 2.0342, 2.1316, 2.231, 2.3319, 2.4341, 2.5361, 2.6397, 2.7391, 2.8383, 2.9398, 3.0485, 3.1463, 3.2481, 3.3292, 3.4049},  {0.934, 0.851, 0.758, 0.656, 0.541, 0.423, 0.365, 0.313, 0.271, 0.225, 0.162, 0.143, 0.126, 0.111, 0.099, 0.088, 0.078, 0.069, 0.061, 0.054, 0.048, 0.043, 0.039, 0.034, 0.031, 0.029, 0.038, 0.039},  {0.826, 0.766, 0.697, 0.617, 0.526, 0.431, 0.384, 0.342, 0.307, 0.270, 0.201, 0.183, 0.167, 0.152, 0.140, 0.129, 0.119, 0.110, 0.101, 0.093, 0.087, 0.081, 0.075, 0.069, 0.065, 0.061, 0.061, 0.072},  {0.756, 0.716, 0.667, 0.606, 0.537, 0.465, 0.429, 0.397, 0.371, 0.345, 0.273, 0.256, 0.242, 0.228, 0.225, 0.215, 0.198, 0.199, 0.182, 0.186, 0.168, 0.163, 0.157, 0.152, 0.147, 0.142, 0.154, 0.134},  {0.907, 0.799, 0.672, 0.552, 0.443, 0.343, 0.294, 0.253, 0.218, 0.178, 0.132, 0.116, 0.102, 0.091, 0.081, 0.072, 0.063, 0.055, 0.047, 0.041, 0.035, 0.029, 0.023, 0.017, 0.013, 0.010, 0.012, 0.010},  {0.804, 0.720, 0.621, 0.527, 0.439, 0.358, 0.318, 0.284, 0.255, 0.223, 0.175, 0.160, 0.146, 0.133, 0.123, 0.113, 0.104, 0.095, 0.086, 0.078, 0.071, 0.063, 0.054, 0.046, 0.040, 0.033, 0.031, 0.032},  {0.734, 0.670, 0.595, 0.522, 0.455, 0.393, 0.362, 0.337, 0.315, 0.295, 0.267, 0.251, 0.224, 0.213, 0.203, 0.194, 0.185, 0.176, 0.168, 0.161, 0.154, 0.146, 0.137, 0.127, 0.120, 0.112, 0.110, 0.093}},
      {{1.0629, 1.1151, 1.1753, 1.2546, 1.3564, 1.4811, 1.5716, 1.6672, 1.7614, 1.8518, 1.9506, 2.0475, 2.1469, 2.2375, 2.3402, 2.4447, 2.5372, 2.6403, 2.733, 2.8354, 2.938, 3.0352, 3.1225, 3.2232, 3.3239, 3.3965, 3.4704},  {1.091, 0.987, 0.864, 0.714, 0.560, 0.429, 0.359, 0.300, 0.254, 0.215, 0.128, 0.113, 0.100, 0.089, 0.078, 0.069, 0.062, 0.055, 0.049, 0.043, 0.039, 0.035, 0.031, 0.028, 0.026, 0.033, 0.037},  {0.947, 0.870, 0.775, 0.658, 0.537, 0.433, 0.376, 0.332, 0.294, 0.258, 0.172, 0.156, 0.143, 0.132, 0.121, 0.111, 0.103, 0.095, 0.088, 0.082, 0.076, 0.071, 0.066, 0.062, 0.059, 0.056, 0.072},  {0.844, 0.789, 0.718, 0.631, 0.539, 0.459, 0.414, 0.382, 0.353, 0.329, 0.247, 0.233, 0.220, 0.210, 0.200, 0.191, 0.183, 0.176, 0.170, 0.163, 0.158, 0.153, 0.148, 0.143, 0.139, 0.135, 0.137},  {0.999, 0.831, 0.672, 0.532, 0.412, 0.315, 0.264, 0.224, 0.191, 0.158, 0.098, 0.086, 0.076, 0.068, 0.059, 0.051, 0.044, 0.038, 0.033, 0.028, 0.023, 0.018, 0.014, 0.011, 0.009, 0.010, 0.010},  {0.871, 0.742, 0.617, 0.506, 0.409, 0.330, 0.288, 0.255, 0.227, 0.201, 0.142, 0.130, 0.119, 0.110, 0.101, 0.091, 0.083, 0.076, 0.069, 0.062, 0.055, 0.048, 0.042, 0.036, 0.031, 0.030, 0.033},  {0.776, 0.679, 0.584, 0.499, 0.423, 0.362, 0.330, 0.305, 0.285, 0.269, 0.246, 0.235, 0.224, 0.216, 0.207, 0.173, 0.165, 0.158, 0.152, 0.145, 0.137, 0.129, 0.122, 0.115, 0.108, 0.106, 0.095}},
      {{1.1279, 1.1796, 1.2481, 1.335, 1.4428, 1.5145, 1.5912, 1.6953, 1.7934, 1.9019, 2.0188, 2.1012, 2.203, 2.3081, 2.4017, 2.5082, 2.6022, 2.6946, 2.7988, 2.8955, 2.9905, 3.0905, 3.1926, 3.2958, 3.3927, 3.4605, 3.5279},  {1.307, 1.147, 0.941, 0.733, 0.552, 0.472, 0.405, 0.323, 0.273, 0.233, 0.189, 0.115, 0.100, 0.087, 0.077, 0.068, 0.060, 0.054, 0.047, 0.042, 0.038, 0.034, 0.030, 0.027, 0.025, 0.028, 0.035},  {1.109, 0.987, 0.831, 0.669, 0.527, 0.463, 0.409, 0.341, 0.299, 0.269, 0.232, 0.159, 0.144, 0.131, 0.120, 0.110, 0.102, 0.094, 0.087, 0.081, 0.075, 0.070, 0.065, 0.060, 0.057, 0.055, 0.066},  {0.956, 0.868, 0.752, 0.632, 0.523, 0.474, 0.432, 0.388, 0.344, 0.325, 0.300, 0.236, 0.222, 0.210, 0.200, 0.190, 0.183, 0.176, 0.169, 0.163, 0.157, 0.152, 0.147, 0.141, 0.136, 0.133, 0.124},  {1.045, 0.826, 0.624, 0.471, 0.355, 0.302, 0.259, 0.216, 0.185, 0.156, 0.124, 0.079, 0.069, 0.059, 0.052, 0.044, 0.038, 0.033, 0.028, 0.024, 0.020, 0.016, 0.013, 0.010, 0.008, 0.009, 0.010},  {0.903, 0.734, 0.574, 0.452, 0.358, 0.314, 0.278, 0.242, 0.216, 0.192, 0.166, 0.128, 0.115, 0.104, 0.095, 0.085, 0.077, 0.070, 0.064, 0.058, 0.052, 0.045, 0.039, 0.034, 0.030, 0.029, 0.030},  {0.794, 0.667, 0.544, 0.449, 0.374, 0.340, 0.311, 0.284, 0.264, 0.247, 0.232, 0.231, 0.221, 0.211, 0.203, 0.194, 0.187, 0.180, 0.174, 0.168, 0.161, 0.153, 0.146, 0.139, 0.106, 0.104, 0.127}},
      {{1.1921, 1.2529, 1.3279, 1.4202, 1.5448, 1.6258, 1.7109, 1.8276, 1.9456, 2.0581, 2.164, 2.2702, 2.3409, 2.4455, 2.5366, 2.6464, 2.7418, 2.8473, 2.94, 3.0334, 3.1361, 3.2417, 3.3381, 3.4403, 3.513, 3.5774},  {1.520, 1.232, 0.942, 0.697, 0.481, 0.423, 0.342, 0.289, 0.239, 0.199, 0.171, 0.143, 0.094, 0.081, 0.072, 0.063, 0.056, 0.049, 0.043, 0.039, 0.035, 0.031, 0.027, 0.024, 0.023, 0.033},  {1.264, 1.048, 0.828, 0.637, 0.464, 0.419, 0.351, 0.302, 0.267, 0.236, 0.210, 0.187, 0.138, 0.124, 0.115, 0.105, 0.097, 0.089, 0.082, 0.077, 0.071, 0.066, 0.061, 0.057, 0.054, 0.062},  {1.062, 0.907, 0.745, 0.601, 0.464, 0.432, 0.375, 0.349, 0.316, 0.290, 0.270, 0.255, 0.215, 0.203, 0.195, 0.185, 0.177, 0.169, 0.163, 0.158, 0.153, 0.147, 0.141, 0.135, 0.132, 0.125},  {1.009, 0.720, 0.517, 0.378, 0.273, 0.230, 0.197, 0.164, 0.137, 0.116, 0.099, 0.082, 0.054, 0.046, 0.041, 0.034, 0.030, 0.026, 0.023, 0.019, 0.016, 0.013, 0.010, 0.008, 0.007, 0.009},  {0.873, 0.647, 0.484, 0.372, 0.284, 0.248, 0.220, 0.192, 0.169, 0.150, 0.136, 0.122, 0.102, 0.092, 0.084, 0.076, 0.070, 0.064, 0.059, 0.054, 0.048, 0.042, 0.036, 0.030, 0.028, 0.028},  {0.766, 0.595, 0.467, 0.377, 0.306, 0.276, 0.255, 0.233, 0.216, 0.202, 0.192, 0.186, 0.201, 0.192, 0.186, 0.178, 0.172, 0.166, 0.161, 0.156, 0.149, 0.143, 0.137, 0.132, 0.130, 0.127}},
    }
};
//
// He star table
const GE_QCRIT_TABLE_HE QCRIT_GE_HE_STAR = {
    {0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0, 2.1, 2.2, 2.3, 2.4, 2.5, 2.6, 2.8, 3.0, 3.2, 3.5, 3.8, 4.2, 4.6, 5.0, 5.4, 5.8, 6.3, 6.9, 7.5, 8.0, 8.6, 9.3, 10.0},
    {
      {{-2.6886, -2.6486, -2.609, -2.5716, -2.5363, -2.5009, -2.4674, -2.4347, -2.4031, -2.3756, -2.3446},  {0.516, 0.506, 0.495, 0.485, 0.474, 0.463, 0.451, 0.438, 0.424, 0.347, 0.329}},
      {{-2.4785, -2.4489, -2.418, -2.3886, -2.3599, -2.3328, -2.3052, -2.2797, -2.2539, -2.2288, -2.2039},  {0.490, 0.483, 0.475, 0.466, 0.458, 0.449, 0.440, 0.431, 0.421, 0.409, 0.391}},
      {{-2.3271, -2.3001, -2.2704, -2.2431, -2.2163, -2.1905, -2.165, -2.1402, -2.1164, -2.0925, -2.0689},  {0.478, 0.471, 0.463, 0.455, 0.447, 0.439, 0.430, 0.421, 0.412, 0.401, 0.383}},
      {{-2.207, -2.1808, -2.1527, -2.126, -2.1008, -2.0757, -2.0509, -2.0277, -2.0044, -1.9818, -1.9589},  {0.472, 0.465, 0.457, 0.449, 0.441, 0.433, 0.424, 0.416, 0.407, 0.397, 0.379}},
      {{-2.1068, -2.084, -2.0552, -2.0321, -2.0058, -1.9804, -1.9563, -1.9345, -1.911, -1.8891, -1.8664},  {0.468, 0.462, 0.454, 0.447, 0.439, 0.430, 0.422, 0.414, 0.404, 0.394, 0.377}},
      {{-2.0215, -1.9984, -1.9739, -1.9483, -1.9222, -1.899, -1.875, -1.8517, -1.8301, -1.808, -1.7866, -1.701, -1.5877},  {0.467, 0.461, 0.454, 0.446, 0.437, 0.429, 0.421, 0.412, 0.403, 0.393, 0.376, 0.194, 0.160}},
      {{-1.9419, -1.9231, -1.899, -1.8755, -1.8507, -1.8262, -1.8033, -1.7829, -1.7589, -1.7382, -1.7162, -1.6156, -1.4951, -1.3803, -1.2697},  {0.466, 0.460, 0.453, 0.445, 0.437, 0.429, 0.421, 0.413, 0.403, 0.394, 0.376, 0.203, 0.179, 0.154, 0.129}},
      {{-1.8794, -1.8607, -1.8361, -1.8116, -1.7862, -1.7609, -1.7404, -1.7173, -1.6956, -1.6747, -1.654, -1.6536, -1.5397, -1.4296, -1.3075, -1.2006, -1.0765, -0.9684, -0.8665},  {0.468, 0.462, 0.454, 0.446, 0.438, 0.429, 0.422, 0.413, 0.404, 0.394, 0.377, 0.229, 0.208, 0.188, 0.166, 0.148, 0.128, 0.111, 0.095}},
      {{-1.819, -1.8011, -1.776, -1.7509, -1.728, -1.7035, -1.6813, -1.6605, -1.6392, -1.6183, -1.5977, -1.5808, -1.4799, -1.3623, -1.2452, -1.1258, -1.02, -0.8915, -0.7769, -0.6609, -0.5465, -0.458, -0.3784, -0.2617, -0.095},  {0.469, 0.463, 0.455, 0.447, 0.439, 0.431, 0.422, 0.414, 0.406, 0.396, 0.378, 0.230, 0.212, 0.192, 0.173, 0.155, 0.140, 0.123, 0.109, 0.096, 0.085, 0.076, 0.069, 0.058, 0.043}},
      {{-1.7639, -1.7436, -1.7218, -1.6993, -1.6758, -1.6504, -1.6292, -1.608, -1.587, -1.5673, -1.5465, -1.5149, -1.4125, -1.3089, -1.1871, -1.0796, -0.9488, -0.851, -0.728, -0.6204, -0.509, -0.4177, -0.2778, -0.1588, -0.0564, 0.0297, 0.089, 0.1673, 0.272, 0.3828, 0.4952, 0.6091, 0.7244, 0.8414, 0.9604, 1.0808, 1.2, 1.3145, 1.4514, 1.5367, 1.6549},  {0.470, 0.463, 0.456, 0.449, 0.441, 0.432, 0.424, 0.416, 0.407, 0.397, 0.380, 0.231, 0.213, 0.196, 0.178, 0.162, 0.145, 0.133, 0.119, 0.108, 0.098, 0.090, 0.080, 0.071, 0.065, 0.060, 0.057, 0.053, 0.048, 0.043, 0.039, 0.035, 0.032, 0.028, 0.025, 0.023, 0.020, 0.018, 0.015, 0.014, 0.012}},
      {{-1.7133, -1.6933, -1.6713, -1.6484, -1.6275, -1.6027, -1.5813, -1.5622, -1.5399, -1.5209, -1.4996, -1.4322, -1.3636, -1.2419, -1.1355, -1.0253, -0.9099, -0.7886, -0.6823, -0.5707, -0.4536, -0.3315, -0.2305, -0.1008, 0.0057, 0.1148, 0.2263, 0.3405, 0.4575, 0.5246, 0.63, 0.7058, 0.8049, 0.9377, 1.0881, 1.2368, 1.3784, 1.4811, 1.6236, 1.7508, 1.8321, 1.9048, 2.0157, 2.0924, 2.1893, 2.3081, 2.429, 2.5533, 2.6861, 2.8422, 3.0746, 3.1325, 3.2661, 3.3366, 3.4452, 3.5505, 3.615, 3.8234, 3.9139, 4.0265, 4.1342},  {0.472, 0.465, 0.458, 0.450, 0.443, 0.434, 0.426, 0.418, 0.409, 0.400, 0.382, 0.228, 0.216, 0.197, 0.181, 0.166, 0.151, 0.137, 0.126, 0.115, 0.105, 0.095, 0.088, 0.079, 0.073, 0.067, 0.062, 0.057, 0.052, 0.050, 0.046, 0.044, 0.041, 0.037, 0.033, 0.030, 0.027, 0.026, 0.023, 0.021, 0.020, 0.019, 0.018, 0.017, 0.016, 0.015, 0.014, 0.013, 0.012, 0.011, 0.010, 0.009, 0.009, 0.008, 0.008, 0.007, 0.007, 0.006, 0.006, 0.006, 0.005}},
      {{-1.6673, -1.645, -1.6196, -1.6, -1.5789, -1.5537, -1.5375, -1.5143, -1.4959, -1.4761, -1.456, -1.3885, -1.3027, -1.1979, -1.0725, -0.9601, -0.8419, -0.7382, -0.6296, -0.5155, -0.3956, -0.2699, -0.1656, -0.0582, 0.0521, 0.1653, 0.2812, 0.4, 0.5227, 0.6179, 0.7501, 0.8539, 0.9604, 1.101, 1.2028, 1.3326, 1.4339, 1.4818, 1.6173, 1.7137, 1.8375, 1.9618, 2.088, 2.2114, 2.3606, 2.4767, 2.5992, 2.6492, 2.7927, 2.8364, 2.9098, 2.992, 3.1166, 3.2766, 3.4377, 3.6175, 3.8265, 3.975, 4.0799, 4.1604, 4.2135, 4.3199, 4.4107, 4.9056},  {0.474, 0.466, 0.458, 0.451, 0.444, 0.434, 0.428, 0.419, 0.411, 0.401, 0.383, 0.231, 0.216, 0.200, 0.182, 0.167, 0.153, 0.141, 0.130, 0.119, 0.109, 0.100, 0.092, 0.086, 0.079, 0.073, 0.068, 0.063, 0.058, 0.054, 0.050, 0.047, 0.043, 0.040, 0.037, 0.035, 0.032, 0.032, 0.029, 0.028, 0.026, 0.024, 0.022, 0.021, 0.019, 0.018, 0.017, 0.016, 0.015, 0.015, 0.014, 0.013, 0.012, 0.011, 0.010, 0.010, 0.009, 0.008, 0.008, 0.007, 0.007, 0.007, 0.006, 0.007}},
      {{-1.5823, -1.5648, -1.5423, -1.5217, -1.4998, -1.4795, -1.4576, -1.4381, -1.4182, -1.3985, -1.3782, -1.3768, -1.3452, -1.3126, -1.2094, -1.0846, -0.9738, -0.8581, -0.7572, -0.6297, -0.5179, -0.4011, -0.3037, -0.1769, -0.0713, 0.0379, 0.1504, 0.2661, 0.3848, 0.507, 0.6017, 0.7329, 0.836, 0.9436, 1.0531, 1.1957, 1.298, 1.4282, 1.5221, 1.644, 1.7636, 1.882, 1.9997, 2.1167, 2.2337, 2.3514, 2.4719, 2.5668, 2.6705, 2.7513, 2.8585, 2.9339, 3.0304, 3.215, 3.3561, 3.4834, 3.6162, 3.7229, 3.8339, 3.9434, 4.0745, 4.2056, 4.3276, 4.4406, 4.5533, 4.6742, 4.7647, 5.1304, 5.2533},  {0.478, 0.471, 0.464, 0.456, 0.449, 0.441, 0.432, 0.424, 0.415, 0.405, 0.387, 0.247, 0.241, 0.235, 0.218, 0.199, 0.184, 0.169, 0.157, 0.144, 0.133, 0.123, 0.115, 0.106, 0.099, 0.092, 0.086, 0.080, 0.075, 0.069, 0.066, 0.061, 0.057, 0.054, 0.051, 0.047, 0.044, 0.041, 0.039, 0.036, 0.034, 0.032, 0.030, 0.028, 0.026, 0.025, 0.023, 0.022, 0.021, 0.020, 0.019, 0.018, 0.017, 0.016, 0.015, 0.014, 0.013, 0.012, 0.012, 0.011, 0.010, 0.010, 0.009, 0.011, 0.019, 0.025, 0.026, 0.011, 0.006}},
      {{-1.5077, -1.3828, -1.3647, -1.3483, -1.3294, -1.3104, -1.2435, -1.1432, -1.0233, -0.9167, -0.8047, -0.6859, -0.5595, -0.448, -0.331, -0.2334, -0.1065, -0.001, 0.1076, 0.2193, 0.3342, 0.4521, 0.5732, 0.6671, 0.7971, 0.8995, 1.0065, 1.1527, 1.2599, 1.363, 1.4941, 1.6195, 1.711, 1.831, 1.95, 2.0683, 2.1859, 2.3036, 2.4221, 2.5127, 2.6404, 2.7471, 2.8331, 2.8869, 3.0404, 3.1239, 3.2777, 3.4238, 3.5606, 3.6741, 3.786, 3.9029, 4.0197, 4.1549, 4.2745, 4.3945, 4.5145, 4.6318, 4.7498, 4.8564, 5.3118},  {0.483, 0.435, 0.427, 0.419, 0.409, 0.393, 0.240, 0.224, 0.206, 0.191, 0.176, 0.163, 0.150, 0.139, 0.129, 0.121, 0.112, 0.105, 0.098, 0.092, 0.086, 0.081, 0.075, 0.072, 0.067, 0.063, 0.059, 0.055, 0.052, 0.049, 0.046, 0.043, 0.041, 0.038, 0.036, 0.034, 0.032, 0.030, 0.028, 0.027, 0.025, 0.024, 0.023, 0.022, 0.021, 0.020, 0.018, 0.017, 0.016, 0.015, 0.014, 0.013, 0.013, 0.012, 0.012, 0.019, 0.036, 0.049, 0.049, 0.041, 0.007}},
      {{-1.4406, -1.4271, -1.4044, -1.3836, -1.3645, -1.3439, -1.3251, -1.3058, -1.2869, -1.2697, -1.2495, -1.1782, -1.0626, -0.9442, -0.8386, -0.7273, -0.609, -0.5046, -0.3945, -0.2789, -0.158, -0.0318, 0.0728, 0.1804, 0.291, 0.4047, 0.5213, 0.6411, 0.7652, 0.8618, 0.9969, 1.1036, 1.2126, 1.3198, 1.4568, 1.5549, 1.6805, 1.8025, 1.9225, 2.012, 2.1306, 2.2486, 2.3666, 2.4854, 2.6073, 2.7039, 2.8102, 2.8946, 3.0142, 3.0984, 3.2527, 3.3334, 3.4877, 3.637, 3.7531, 3.8651, 3.9814, 4.0967, 4.2278, 4.3558, 4.4754, 4.5844, 4.7017, 4.8192, 4.937, 5.3836, 5.4979},  {0.487, 0.480, 0.472, 0.465, 0.457, 0.449, 0.441, 0.433, 0.424, 0.415, 0.397, 0.243, 0.224, 0.207, 0.192, 0.178, 0.165, 0.154, 0.144, 0.134, 0.124, 0.115, 0.109, 0.102, 0.096, 0.090, 0.085, 0.079, 0.074, 0.070, 0.066, 0.062, 0.059, 0.056, 0.052, 0.049, 0.046, 0.043, 0.041, 0.039, 0.037, 0.035, 0.033, 0.031, 0.029, 0.028, 0.026, 0.025, 0.024, 0.023, 0.021, 0.020, 0.019, 0.017, 0.016, 0.015, 0.015, 0.014, 0.013, 0.014, 0.028, 0.053, 0.068, 0.065, 0.052, 0.010, 0.005}},
      {{-1.3799, -1.3656, -1.3467, -1.3258, -1.3066, -1.2861, -1.2674, -1.2505, -1.2317, -1.2131, -1.1943, -1.1067, -0.991, -0.8555, -0.7491, -0.6369, -0.538, -0.4125, -0.3244, -0.2095, -0.0402, 0.0608, 0.1913, 0.2717, 0.3814, 0.4941, 0.6096, 0.7283, 0.8512, 0.9469, 1.0463, 1.186, 1.2938, 1.3999, 1.5358, 1.6332, 1.7582, 1.8797, 1.9993, 2.1183, 2.2074, 2.3254, 2.4435, 2.5624, 2.6848, 2.782, 2.8895, 2.9766, 3.0331, 3.1913, 3.3343, 3.4193, 3.5645, 3.7097, 3.8261, 3.9394, 4.0555, 4.1939, 4.3138, 4.4318, 4.5515, 4.6707, 4.7884, 4.8982, 5.0166, 5.4675, 5.5836},  {0.492, 0.484, 0.477, 0.469, 0.462, 0.454, 0.446, 0.438, 0.429, 0.419, 0.400, 0.243, 0.225, 0.205, 0.191, 0.177, 0.166, 0.154, 0.146, 0.136, 0.123, 0.116, 0.108, 0.103, 0.097, 0.092, 0.086, 0.081, 0.076, 0.072, 0.069, 0.064, 0.061, 0.058, 0.054, 0.051, 0.048, 0.045, 0.043, 0.040, 0.039, 0.036, 0.034, 0.032, 0.030, 0.029, 0.028, 0.026, 0.026, 0.024, 0.022, 0.021, 0.020, 0.018, 0.017, 0.016, 0.016, 0.015, 0.014, 0.015, 0.037, 0.070, 0.082, 0.074, 0.057, 0.010, 0.005}},
      {{-1.3241, -1.3024, -1.2843, -1.2674, -1.2486, -1.2313, -1.2129, -1.1965, -1.1801, -1.162, -1.1444, -1.0679, -0.9539, -0.8393, -0.7203, -0.6123, -0.4973, -0.3748, -0.2668, -0.1535, -0.0353, 0.0876, 0.1893, 0.2937, 0.428, 0.5386, 0.6519, 0.768, 0.8877, 0.9806, 1.1096, 1.2114, 1.3168, 1.4573, 1.5591, 1.689, 1.783, 1.9048, 2.0243, 2.1429, 2.261, 2.3786, 2.4964, 2.615, 2.7061, 2.8333, 2.9387, 3.0212, 3.1347, 3.2126, 3.3678, 3.4453, 3.6035, 3.72, 3.8582, 3.9739, 4.0934, 4.21, 4.3464, 4.4702, 4.5872, 4.7094, 4.823, 4.937, 5.0521, 5.5163, 5.6259},  {0.497, 0.485, 0.478, 0.472, 0.464, 0.457, 0.449, 0.441, 0.433, 0.422, 0.406, 0.249, 0.230, 0.213, 0.197, 0.183, 0.171, 0.158, 0.148, 0.139, 0.129, 0.121, 0.114, 0.108, 0.101, 0.095, 0.090, 0.084, 0.079, 0.076, 0.071, 0.068, 0.064, 0.060, 0.057, 0.053, 0.051, 0.048, 0.045, 0.043, 0.040, 0.038, 0.036, 0.034, 0.033, 0.031, 0.029, 0.028, 0.026, 0.025, 0.024, 0.023, 0.021, 0.020, 0.019, 0.018, 0.017, 0.016, 0.015, 0.015, 0.033, 0.075, 0.091, 0.084, 0.066, 0.012, 0.006}},
      {{-1.2748, -1.2617, -1.243, -1.2222, -1.2062, -1.1857, -1.1696, -1.1506, -1.1342, -1.1161, -1.0983, -1.0982, -1.0048, -0.8907, -0.7758, -0.6568, -0.5491, -0.4349, -0.3133, -0.2063, -0.0944, 0.0223, 0.1434, 0.2435, 0.3723, 0.4781, 0.5866, 0.6978, 0.8115, 0.9283, 1.0492, 1.1434, 1.2752, 1.3784, 1.5165, 1.6172, 1.7462, 1.8394, 1.9606, 2.0795, 2.1975, 2.3151, 2.4324, 2.5497, 2.6679, 2.7579, 2.884, 2.9872, 3.0673, 3.1718, 3.2434, 3.4066, 3.4862, 3.6381, 3.7908, 3.9187, 4.0436, 4.1386, 4.2619, 4.3981, 4.5202, 4.6456, 4.7649, 4.8767, 4.9929, 5.1107, 5.5712, 5.6898},  {0.500, 0.493, 0.486, 0.477, 0.471, 0.462, 0.455, 0.446, 0.438, 0.427, 0.410, 0.409, 0.248, 0.230, 0.213, 0.197, 0.184, 0.171, 0.159, 0.149, 0.140, 0.131, 0.123, 0.116, 0.109, 0.103, 0.097, 0.092, 0.087, 0.082, 0.077, 0.074, 0.069, 0.066, 0.061, 0.059, 0.055, 0.053, 0.050, 0.047, 0.044, 0.042, 0.040, 0.037, 0.035, 0.034, 0.032, 0.030, 0.029, 0.028, 0.027, 0.025, 0.024, 0.022, 0.021, 0.019, 0.018, 0.017, 0.016, 0.015, 0.015, 0.023, 0.072, 0.095, 0.089, 0.069, 0.013, 0.007}},
      {{-1.2282, -1.2134, -1.1914, -1.1777, -1.1586, -1.1412, -1.1252, -1.1085, -1.0901, -1.0739, -1.0556, -0.9596, -0.8452, -0.7317, -0.6154, -0.5103, -0.3985, -0.2796, -0.1536, -0.0433, 0.0713, 0.1902, 0.2883, 0.4144, 0.5179, 0.6239, 0.7598, 0.8712, 0.9855, 1.1045, 1.1977, 1.3272, 1.4282, 1.53, 1.6629, 1.7907, 1.8831, 2.0032, 2.1209, 2.2378, 2.3544, 2.4709, 2.5875, 2.7049, 2.8246, 2.9183, 3.019, 3.1355, 3.243, 3.3166, 3.4013, 3.5653, 3.7081, 3.848, 3.9804, 4.0804, 4.182, 4.312, 4.4301, 4.5686, 4.6919, 4.8085, 4.9281, 5.042, 5.1558, 5.62, 5.7341},  {0.505, 0.496, 0.488, 0.482, 0.474, 0.466, 0.459, 0.451, 0.442, 0.432, 0.411, 0.250, 0.231, 0.215, 0.199, 0.186, 0.174, 0.162, 0.151, 0.141, 0.133, 0.124, 0.118, 0.111, 0.105, 0.099, 0.093, 0.088, 0.083, 0.078, 0.075, 0.070, 0.067, 0.064, 0.060, 0.056, 0.054, 0.051, 0.048, 0.046, 0.043, 0.041, 0.039, 0.037, 0.034, 0.033, 0.031, 0.030, 0.028, 0.027, 0.026, 0.024, 0.023, 0.021, 0.020, 0.019, 0.018, 0.017, 0.016, 0.015, 0.016, 0.041, 0.091, 0.091, 0.072, 0.015, 0.008}},
      {{-1.1861, -1.173, -1.1544, -1.1372, -1.1181, -1.1006, -1.0847, -1.0681, -1.0499, -1.0339, -1.0164, -0.9102, -0.7957, -0.6808, -0.563, -0.4575, -0.3462, -0.2285, -0.1253, 0.0041, 0.1164, 0.2327, 0.3284, 0.4514, 0.5777, 0.6812, 0.787, 0.9222, 1.0334, 1.1479, 1.2367, 1.3599, 1.4894, 1.5884, 1.7183, 1.8128, 1.9346, 2.0527, 2.1691, 2.2848, 2.4004, 2.5158, 2.6314, 2.7478, 2.8662, 2.9905, 3.0916, 3.1678, 3.2648, 3.329, 3.4874, 3.6441, 3.7595, 3.9044, 4.0054, 4.112, 4.2238, 4.37, 4.468, 4.6143, 4.7423, 4.854, 4.969, 5.0823, 5.2022},  {0.509, 0.502, 0.494, 0.487, 0.479, 0.471, 0.464, 0.456, 0.446, 0.436, 0.417, 0.251, 0.232, 0.215, 0.200, 0.187, 0.175, 0.163, 0.154, 0.143, 0.134, 0.126, 0.120, 0.113, 0.106, 0.101, 0.096, 0.089, 0.085, 0.080, 0.077, 0.073, 0.068, 0.065, 0.061, 0.059, 0.056, 0.053, 0.050, 0.047, 0.045, 0.042, 0.040, 0.038, 0.036, 0.034, 0.032, 0.031, 0.030, 0.029, 0.027, 0.025, 0.024, 0.022, 0.021, 0.020, 0.019, 0.018, 0.017, 0.016, 0.017, 0.054, 0.087, 0.090, 0.075}},
      {{-1.1391, -1.1296, -1.1048, -1.0857, -1.0706, -1.0566, -1.0413, -1.0254, -1.0101, -0.9963, -0.98, -0.9799, -0.8779, -0.763, -0.6477, -0.5301, -0.4254, -0.3157, -0.1999, -0.0783, 0.0275, 0.1594, 0.2732, 0.3669, 0.4871, 0.6103, 0.7111, 0.8399, 0.9452, 1.0527, 1.1631, 1.2771, 1.3957, 1.5201, 1.6159, 1.7429, 1.8667, 1.9862, 2.1027, 2.2176, 2.3319, 2.4461, 2.5598, 2.674, 2.7892, 2.9065, 2.9976, 3.1282, 3.2022, 3.2923, 3.4239, 3.501, 3.6543, 3.7816, 3.9384, 4.0383, 4.1733, 4.2771, 4.4096, 4.5277, 4.6302, 4.7659, 4.8901, 5.0054, 5.1315},  {0.515, 0.505, 0.494, 0.486, 0.479, 0.473, 0.465, 0.457, 0.449, 0.440, 0.422, 0.421, 0.254, 0.235, 0.217, 0.202, 0.189, 0.177, 0.165, 0.154, 0.145, 0.136, 0.128, 0.122, 0.114, 0.108, 0.102, 0.096, 0.091, 0.087, 0.082, 0.078, 0.074, 0.070, 0.067, 0.063, 0.060, 0.056, 0.053, 0.051, 0.048, 0.046, 0.043, 0.041, 0.039, 0.037, 0.035, 0.033, 0.032, 0.031, 0.029, 0.028, 0.026, 0.024, 0.023, 0.022, 0.020, 0.019, 0.018, 0.017, 0.017, 0.016, 0.033, 0.074, 0.087}},
      {{-1.1002, -1.0899, -1.0697, -1.0486, -1.0336, -1.0197, -1.0046, -0.9891, -0.9758, -0.9607, -0.9457, -0.8334, -0.7182, -0.6025, -0.5019, -0.3805, -0.2714, -0.1571, -0.0375, 0.0661, 0.1949, 0.306, 0.4201, 0.5372, 0.6571, 0.7549, 0.8797, 0.9815, 1.1113, 1.2181, 1.3279, 1.4418, 1.5608, 1.6834, 1.7754, 1.8951, 2.0121, 2.1257, 2.2662, 2.3777, 2.4892, 2.6013, 2.7139, 2.8273, 2.9426, 3.0623, 3.159, 3.2684, 3.363, 3.4257, 3.5679, 3.7112, 3.8697, 3.9732, 4.0717, 4.2033, 4.3044, 4.4347, 4.5521, 4.675, 4.8174},  {0.520, 0.508, 0.499, 0.490, 0.483, 0.477, 0.469, 0.461, 0.453, 0.443, 0.426, 0.254, 0.235, 0.218, 0.204, 0.189, 0.177, 0.166, 0.155, 0.147, 0.137, 0.129, 0.122, 0.115, 0.108, 0.103, 0.097, 0.093, 0.087, 0.083, 0.079, 0.075, 0.071, 0.067, 0.065, 0.061, 0.058, 0.055, 0.052, 0.049, 0.047, 0.044, 0.042, 0.040, 0.038, 0.036, 0.034, 0.033, 0.031, 0.030, 0.028, 0.027, 0.025, 0.024, 0.022, 0.021, 0.020, 0.019, 0.018, 0.017, 0.018}},
      {{-1.0724, -1.0605, -1.0318, -1.0252, -1.0095, -0.9922, -0.9766, -0.9604, -0.9466, -0.9295, -0.9136, -0.8008, -0.6847, -0.5848, -0.45, -0.3462, -0.2388, -0.1267, -0.0295, 0.1119, 0.2165, 0.3458, 0.4566, 0.57, 0.6859, 0.8041, 0.9002, 1.0224, 1.1472, 1.2491, 1.3537, 1.489, 1.6021, 1.7185, 1.835, 1.9492, 2.0614, 2.1708, 2.2797, 2.3877, 2.5242, 2.6332, 2.7432, 2.8547, 2.9681, 3.0857, 3.2122, 3.2839, 3.3691, 3.4888, 3.6226, 3.6904, 3.8745, 4.0199, 4.1171, 4.2455, 4.3909},  {0.523, 0.515, 0.502, 0.499, 0.492, 0.484, 0.476, 0.468, 0.460, 0.448, 0.430, 0.256, 0.236, 0.221, 0.203, 0.190, 0.179, 0.168, 0.159, 0.147, 0.139, 0.130, 0.123, 0.116, 0.110, 0.104, 0.099, 0.094, 0.089, 0.085, 0.081, 0.076, 0.072, 0.069, 0.065, 0.062, 0.059, 0.056, 0.054, 0.051, 0.048, 0.046, 0.044, 0.042, 0.039, 0.037, 0.035, 0.034, 0.033, 0.031, 0.029, 0.028, 0.026, 0.024, 0.023, 0.022, 0.021}},
      {{-1.0288, -1.0144, -0.9971, -0.9811, -0.9635, -0.9477, -0.9312, -0.9154, -0.8998, -0.883, -0.7793, -0.6636, -0.5472, -0.4293, -0.3082, -0.2006, -0.0892, 0.0261, 0.1453, 0.2473, 0.373, 0.4802, 0.6117, 0.7235, 0.8372, 0.9527, 1.0698, 1.1649, 1.2858, 1.4099, 1.512, 1.6175, 1.7542, 1.8643, 1.9736, 2.081, 2.2143, 2.3212, 2.4297, 2.5418, 2.7041},  {0.529, 0.513, 0.506, 0.498, 0.490, 0.482, 0.473, 0.464, 0.453, 0.433, 0.259, 0.239, 0.221, 0.205, 0.191, 0.179, 0.168, 0.158, 0.148, 0.140, 0.131, 0.125, 0.117, 0.111, 0.105, 0.100, 0.094, 0.090, 0.086, 0.081, 0.078, 0.074, 0.070, 0.067, 0.064, 0.061, 0.058, 0.055, 0.053, 0.050, 0.047}},
      {{-0.9957, -0.9864, -0.9735, -0.9549, -0.9359, -0.9224, -0.9076, -0.8947, -0.8818, -0.8675, -0.854, -0.8535, -0.7384, -0.6225, -0.5224, -0.4041, -0.2833, -0.1768, -0.0483, 0.0654, 0.1625, 0.2819, 0.404, 0.5077, 0.6343, 0.7415, 0.85, 0.982, 1.0932, 1.2056, 1.3196, 1.4359, 1.5556, 1.655, 1.7851, 1.8949},  {0.534, 0.519, 0.513, 0.504, 0.495, 0.488, 0.481, 0.473, 0.465, 0.455, 0.437, 0.281, 0.258, 0.238, 0.223, 0.207, 0.192, 0.180, 0.168, 0.158, 0.150, 0.141, 0.132, 0.126, 0.118, 0.112, 0.107, 0.101, 0.096, 0.091, 0.087, 0.083, 0.079, 0.075, 0.072, 0.068}},
      {{-0.9642, -0.942, -0.9276, -0.9153, -0.9011, -0.8907, -0.8763, -0.8638, -0.8515, -0.8393, -0.8264, -0.8244, -0.7094, -0.5935, -0.4764, -0.3748, -0.2541, -0.1482, -0.0213, 0.0905, 0.2048, 0.3214, 0.4401, 0.5403, 0.6621, 0.7852, 0.8886, 1.0135, 1.1183, 1.245, 1.3521, 1.4622, 1.606},  {0.539, 0.516, 0.509, 0.504, 0.496, 0.491, 0.483, 0.476, 0.468, 0.458, 0.441, 0.282, 0.259, 0.239, 0.221, 0.207, 0.193, 0.181, 0.168, 0.158, 0.149, 0.140, 0.132, 0.126, 0.119, 0.112, 0.107, 0.102, 0.097, 0.092, 0.088, 0.084, 0.079}},
      {{-0.9367, -0.9367, -0.918, -0.9015, -0.8859, -0.8682, -0.8475, -0.8371, -0.8246, -0.8125, -0.7996, -0.6873, -0.5711, -0.4537, -0.3348, -0.2314, -0.1084, -0.0007, 0.1091, 0.2396, 0.3532, 0.4682, 0.5842, 0.6814, 0.7985, 0.9157, 1.033, 1.1505, 1.2696, 1.4009},  {0.532, 0.532, 0.522, 0.514, 0.507, 0.498, 0.486, 0.480, 0.471, 0.462, 0.444, 0.261, 0.240, 0.222, 0.206, 0.194, 0.180, 0.170, 0.160, 0.149, 0.141, 0.133, 0.126, 0.120, 0.114, 0.108, 0.103, 0.098, 0.093, 0.088}},
      {{-0.9055, -0.8809, -0.87, -0.8576, -0.8437, -0.8335, -0.8218, -0.8096, -0.7976, -0.7861, -0.7743, -0.6605, -0.544, -0.4261, -0.3069, -0.2035, -0.081, 0.0256, 0.1518, 0.2612, 0.2612, 0.4829, 0.6129, 0.7242, 0.8349, 0.945, 1.073, 1.1978},  {0.547, 0.522, 0.516, 0.510, 0.502, 0.497, 0.490, 0.482, 0.474, 0.464, 0.447, 0.261, 0.240, 0.222, 0.206, 0.194, 0.180, 0.170, 0.159, 0.150, 0.150, 0.134, 0.126, 0.120, 0.114, 0.109, 0.103, 0.098}},
      {{-0.8807, -0.8592, -0.8406, -0.8279, -0.8176, -0.8079, -0.7964, -0.7844, -0.7727, -0.7627, -0.7504, -0.7497, -0.6351, -0.5189, -0.401, -0.2817, -0.1786, -0.0569, 0.0485, 0.1726, 0.2798, 0.4053, 0.5127, 0.6373, 0.7428, 0.8642, 0.9842, 1.1092},  {0.539, 0.527, 0.518, 0.511, 0.506, 0.500, 0.493, 0.485, 0.477, 0.469, 0.450, 0.285, 0.261, 0.241, 0.223, 0.206, 0.194, 0.181, 0.170, 0.159, 0.151, 0.142, 0.134, 0.127, 0.121, 0.115, 0.109, 0.104}},
      {{-0.848, -0.844, -0.825, -0.8074, -0.7892, -0.7752, -0.7625, -0.7468, -0.7337, -0.7197, -0.7048, -0.6004, -0.4836, -0.3656, -0.2466, -0.1438, -0.0236, 0.0967, 0.1996, 0.319, 0.437, 0.5527, 0.6651, 0.789, 0.9164},  {0.556, 0.552, 0.542, 0.533, 0.524, 0.516, 0.508, 0.499, 0.490, 0.478, 0.457, 0.265, 0.244, 0.225, 0.209, 0.196, 0.183, 0.171, 0.162, 0.152, 0.144, 0.136, 0.129, 0.123, 0.116}},
      {{-0.8019, -0.7969, -0.7866, -0.765, -0.7496, -0.7311, -0.7185, -0.7029, -0.6901, -0.6766, -0.6624, -0.6488, -0.5331, -0.4157, -0.3141, -0.1949, -0.0757, 0.0431, 0.1609, 0.2769, 0.3902, 0.4998, 0.6192, 0.7319},  {0.564, 0.559, 0.553, 0.542, 0.534, 0.523, 0.515, 0.505, 0.496, 0.485, 0.463, 0.283, 0.260, 0.239, 0.223, 0.207, 0.192, 0.180, 0.168, 0.158, 0.150, 0.142, 0.135, 0.128}},
      {{-0.7613, -0.7548, -0.7435, -0.7317, -0.7108, -0.6931, -0.6769, -0.662, -0.6515, -0.6384, -0.6245, -0.529, -0.4122, -0.2943, -0.1762, -0.0586, 0.0578, 0.1721, 0.2832, 0.4049, 0.5192, 0.6373},  {0.573, 0.568, 0.561, 0.555, 0.543, 0.532, 0.522, 0.512, 0.504, 0.492, 0.470, 0.267, 0.246, 0.227, 0.210, 0.195, 0.183, 0.171, 0.162, 0.152, 0.144, 0.137}},
      {{-0.6991, -0.6867, -0.6743, -0.6643, -0.6447, -0.6301, -0.6185, -0.6064, -0.595, -0.5832, -0.5711, -0.513, -0.4467, -0.3799, -0.3128, -0.2456, -0.1786, -0.1119, -0.0292, 0.036, 0.1003, 0.1788, 0.2398, 0.3134, 0.3834, 0.4494, 0.5226, 0.5896},  {0.582, 0.574, 0.567, 0.561, 0.549, 0.539, 0.531, 0.522, 0.513, 0.501, 0.479, 0.275, 0.262, 0.249, 0.238, 0.227, 0.217, 0.208, 0.198, 0.190, 0.183, 0.175, 0.170, 0.163, 0.158, 0.153, 0.148, 0.144}},
      {{-0.6517, -0.6416, -0.6064, -0.5957, -0.585, -0.5707, -0.5587, -0.5474, -0.5354, -0.5227, -0.4789, -0.4122, -0.3451, -0.2781, -0.2113, -0.1283, -0.0626, 0.0022, 0.0659, 0.1433, 0.2178, 0.2885, 0.3548, 0.4275, 0.4924},  {0.598, 0.590, 0.569, 0.561, 0.554, 0.544, 0.535, 0.525, 0.513, 0.488, 0.278, 0.264, 0.251, 0.239, 0.229, 0.216, 0.207, 0.199, 0.192, 0.183, 0.176, 0.169, 0.164, 0.158, 0.153}},
      {{-0.5886, -0.5764, -0.5603, -0.5442, -0.5284, -0.5163, -0.5069, -0.4869, -0.4757, -0.464, -0.4309, -0.3643, -0.2975, -0.2308, -0.1478, -0.0821, -0.0175, 0.0613, 0.1221, 0.1945, 0.2622, 0.336, 0.4102, 0.4787},  {0.615, 0.606, 0.595, 0.584, 0.573, 0.563, 0.556, 0.538, 0.525, 0.502, 0.280, 0.266, 0.252, 0.240, 0.227, 0.217, 0.208, 0.198, 0.191, 0.183, 0.176, 0.169, 0.163, 0.159}},
      {{-0.5367, -0.5315, -0.5175, -0.5042, -0.4918, -0.4778, -0.4633, -0.451, -0.4377, -0.4246, -0.4109, -0.3919, -0.3084, -0.2414, -0.1745, -0.1081, -0.0428, 0.0366, 0.0976, 0.1696, 0.2484, 0.3176, 0.3923, 0.4691},  {0.638, 0.631, 0.621, 0.612, 0.603, 0.592, 0.580, 0.569, 0.556, 0.541, 0.520, 0.282, 0.264, 0.251, 0.239, 0.228, 0.218, 0.206, 0.198, 0.190, 0.181, 0.175, 0.169, 0.164}},
      {{-0.4853, -0.4799, -0.4691, -0.4554, -0.4401, -0.4266, -0.4147, -0.4008, -0.3881, -0.3755, -0.3622, -0.3405, -0.2732, -0.2056, -0.1382, -0.0717, 0.0093, 0.0713, 0.1444, 0.2112, 0.2921, 0.3643, 0.4407},  {0.652, 0.650, 0.640, 0.630, 0.618, 0.606, 0.596, 0.582, 0.571, 0.559, 0.539, 0.281, 0.266, 0.252, 0.239, 0.228, 0.215, 0.206, 0.197, 0.189, 0.181, 0.175, 0.169}},
      {{-0.4429, -0.4377, -0.4252, -0.4095, -0.3944, -0.3811, -0.3695, -0.3558, -0.3433, -0.3299, -0.3171, -0.3099, -0.2419, -0.1734, -0.1053, -0.0385, 0.0263, 0.1028, 0.1732, 0.2476, 0.3174, 0.3917},  {0.675, 0.670, 0.660, 0.647, 0.634, 0.622, 0.612, 0.601, 0.590, 0.578, 0.558, 0.281, 0.266, 0.251, 0.238, 0.227, 0.216, 0.205, 0.196, 0.187, 0.180, 0.173}},
      {{-0.3454, -0.3328, -0.3263, -0.3165, -0.3051, -0.2959, -0.2853, -0.2751, -0.2739, -0.204, -0.1339, -0.0642, 0.0037, 0.0847, 0.1454, 0.2271, 0.2951, 0.3703, 0.447},  {0.702, 0.636, 0.630, 0.622, 0.613, 0.605, 0.595, 0.577, 0.283, 0.266, 0.251, 0.238, 0.225, 0.212, 0.203, 0.193, 0.185, 0.177, 0.169}},
      {{-0.3374, -0.3294, -0.3203, -0.3082, -0.2947, -0.2829, -0.2727, -0.2605, -0.2494, -0.2374, -0.226, -0.2184, -0.1636, -0.0905, -0.0184, 0.0518, 0.1179, 0.1947, 0.2624, 0.3381, 0.4155},  {0.707, 0.698, 0.690, 0.680, 0.669, 0.659, 0.651, 0.641, 0.631, 0.619, 0.600, 0.278, 0.265, 0.249, 0.234, 0.222, 0.211, 0.199, 0.190, 0.181, 0.173}},
      {{-0.2814, -0.2759, -0.2681, -0.2523, -0.2377, -0.2255, -0.2137, -0.2041, -0.1924, -0.1821, -0.1709, -0.1308, -0.0735, 0.0029, 0.0769, 0.1464, 0.2112, 0.2845, 0.3553, 0.4296, 0.5077},  {0.736, 0.729, 0.723, 0.709, 0.698, 0.687, 0.677, 0.668, 0.657, 0.647, 0.627, 0.265, 0.252, 0.236, 0.221, 0.209, 0.199, 0.189, 0.180, 0.171, 0.161}},
      {{-0.2439, -0.2434, -0.2267, -0.2112, -0.1951, -0.1832, -0.1708, -0.157, -0.1448, -0.1329, -0.1198, -0.0686, 0.0082, 0.0858, 0.1407, 0.2272, 0.292, 0.3613, 0.4354, 0.5159, 0.5721},  {0.781, 0.776, 0.760, 0.746, 0.733, 0.722, 0.712, 0.700, 0.689, 0.677, 0.654, 0.256, 0.239, 0.224, 0.214, 0.199, 0.190, 0.180, 0.171, 0.162, 0.158}},
      {{-0.2039, -0.1953, -0.1899, -0.1771, -0.1556, -0.1434, -0.1293, -0.1175, -0.104, -0.0915, -0.0789, -0.0129, 0.0599, 0.1337, 0.2068, 0.2765, 0.3422, 0.4126, 0.4884, 0.5674},  {0.801, 0.793, 0.789, 0.776, 0.757, 0.746, 0.734, 0.724, 0.711, 0.698, 0.675, 0.252, 0.232, 0.218, 0.205, 0.194, 0.184, 0.175, 0.165, 0.157}},
      {{-0.1677, -0.1592, -0.1421, -0.1298, -0.1159, -0.1013, -0.0868, -0.0731, -0.0593, -0.0454, -0.0312, 0.017, 0.0972, 0.1599, 0.2406, 0.2979, 0.3702, 0.4504, 0.5203, 0.5983, 0.6735},  {0.841, 0.831, 0.814, 0.803, 0.791, 0.778, 0.765, 0.752, 0.740, 0.725, 0.699, 0.262, 0.241, 0.226, 0.209, 0.197, 0.184, 0.172, 0.164, 0.155, 0.147}},
      {{-0.1202, -0.1111, -0.0941, -0.0788, -0.0626, -0.0486, -0.0346, -0.0195, -0.0049, 0.0095, 0.0238, 0.0382, 0.0953, 0.1756, 0.2393, 0.3211, 0.381, 0.4533, 0.5217, 0.5951, 0.674, 0.7481, 0.8196},  {0.878, 0.867, 0.850, 0.835, 0.821, 0.808, 0.796, 0.782, 0.768, 0.752, 0.725, 0.290, 0.271, 0.248, 0.233, 0.215, 0.203, 0.189, 0.177, 0.165, 0.151, 0.142, 0.135}},
      {{-0.0689, -0.0609, -0.0474, -0.0336, -0.0198, -0.0, 0.0233, 0.0385, 0.052, 0.0657, 0.0802, 0.1029, 0.1627, 0.232, 0.3074, 0.3674, 0.4451, 0.5197, 0.5885, 0.6534, 0.7327, 0.8048, 0.8844},  {0.911, 0.899, 0.885, 0.872, 0.858, 0.840, 0.819, 0.804, 0.791, 0.775, 0.743, 0.319, 0.292, 0.267, 0.245, 0.230, 0.213, 0.198, 0.185, 0.174, 0.160, 0.147, 0.133}},
    }
};

// B coefficients
// Table in Appendix A of Hurley et al. 2000
// Key to map is n (B(n)).  Map element is unordered_map of term coefficient values.
const std::map<int, COMPASUnorderedMap<AB_TCoeff, double>> B_COEFF = {
    { 1, {{ALPHA,  3.970000E-1}, {BETA,  2.882600E-1}, {GAMMA,  5.293000E-1}, {ETA,  0.000000E0 }, {MU,  0.000000E0 }}},
    { 2, {{ALPHA,  0.000000E0 }, {BETA,  0.000000E0 }, {GAMMA,  0.000000E0 }, {ETA,  0.000000E0 }, {MU,  0.000000E0 }}},
    { 3, {{ALPHA,  0.000000E0 }, {BETA,  0.000000E0 }, {GAMMA,  0.000000E0 }, {ETA,  0.000000E0 }, {MU,  0.000000E0 }}},
    { 4, {{ALPHA,  9.960283E-1}, {BETA,  8.164393E-1}, {GAMMA,  2.383830E0 }, {ETA,  2.223436E0 }, {MU,  8.638115E-1}}},
    { 5, {{ALPHA,  2.561062E-1}, {BETA,  7.072646E-2}, {GAMMA, -5.444596E-2}, {ETA, -5.798167E-2}, {MU, -1.349129E-2}}},
    { 6, {{ALPHA,  1.157338E0 }, {BETA,  1.467883E0 }, {GAMMA,  4.299661E0 }, {ETA,  3.130500E0 }, {MU,  6.992080E-1}}},
    { 7, {{ALPHA,  4.022765E-1}, {BETA,  3.050010E-1}, {GAMMA,  9.962137E-1}, {ETA,  7.914079E-1}, {MU,  1.728098E-1}}},
    { 8, {{ALPHA,  0.000000E0 }, {BETA,  0.000000E0 }, {GAMMA,  0.000000E0 }, {ETA,  0.000000E0 }, {MU,  0.000000E0 }}},
    { 9, {{ALPHA,  2.751631E3 }, {BETA,  3.557098E2 }, {GAMMA,  0.000000E0 }, {ETA,  0.000000E0 }, {MU,  0.000000E0 }}},
    {10, {{ALPHA, -3.820831E-2}, {BETA,  5.872664E-2}, {GAMMA,  0.000000E0 }, {ETA,  0.000000E0 }, {MU,  0.000000E0 }}},

    {11, {{ALPHA,  1.071738E2 }, {BETA, -8.970339E1 }, {GAMMA, -3.949739E1 }, {ETA,  0.000000E0 }, {MU,  0.000000E0 }}},
    {12, {{ALPHA,  7.348793E2 }, {BETA, -1.531020E2 }, {GAMMA, -3.793700E1 }, {ETA,  0.000000E0 }, {MU,  0.000000E0 }}},
    {13, {{ALPHA,  9.219293E0 }, {BETA, -2.005865E0 }, {GAMMA, -5.561309E-1}, {ETA,  0.000000E0 }, {MU,  0.000000E0 }}},
    {14, {{ALPHA,  2.917412E0 }, {BETA,  1.575290E0 }, {GAMMA,  5.751814E-1}, {ETA,  0.000000E0 }, {MU,  0.000000E0 }}},
    {15, {{ALPHA,  3.629118E0 }, {BETA, -9.112722E-1}, {GAMMA,  1.042291E0 }, {ETA,  0.000000E0 }, {MU,  0.000000E0 }}},
    {16, {{ALPHA,  4.916389E0 }, {BETA,  2.862149E0 }, {GAMMA,  7.844850E-1}, {ETA,  0.000000E0 }, {MU,  0.000000E0 }}},
    {17, {{ALPHA,  0.000000E0 }, {BETA,  0.000000E0 }, {GAMMA,  0.000000E0 }, {ETA,  0.000000E0 }, {MU,  0.000000E0 }}},
    {18, {{ALPHA,  5.496045E1 }, {BETA, -1.289968E1 }, {GAMMA,  6.385758E0 }, {ETA,  0.000000E0 }, {MU,  0.000000E0 }}},
    {19, {{ALPHA,  1.832694E0 }, {BETA, -5.766608E-2}, {GAMMA,  5.696128E-2}, {ETA,  0.000000E0 }, {MU,  0.000000E0 }}},
    {20, {{ALPHA,  1.211104E2 }, {BETA,  0.000000E0 }, {GAMMA,  0.000000E0 }, {ETA,  0.000000E0 }, {MU,  0.000000E0 }}},

    {21, {{ALPHA,  2.214088E2 }, {BETA,  2.187113E2 }, {GAMMA,  1.170177E1 }, {ETA, -2.635340E1 }, {MU,  0.000000E0 }}},
    {22, {{ALPHA,  2.063983E0 }, {BETA,  7.363827E-1}, {GAMMA,  2.654323E-1}, {ETA, -6.140719E-2}, {MU,  0.000000E0 }}},
    {23, {{ALPHA,  2.003160E0 }, {BETA,  9.388871E-1}, {GAMMA,  9.656450E-1}, {ETA,  2.362266E-1}, {MU,  0.000000E0 }}},
    {24, {{ALPHA,  1.609901E1 }, {BETA,  7.391573E0 }, {GAMMA,  2.277010E1 }, {ETA,  8.334227E0 }, {MU,  0.000000E0 }}},
    {25, {{ALPHA,  1.747500E-1}, {BETA,  6.271202E-2}, {GAMMA, -2.324229E-2}, {ETA, -1.844559E-2}, {MU,  0.000000E0 }}},
    {26, {{ALPHA,  0.000000E0 }, {BETA,  0.000000E0 }, {GAMMA,  0.000000E0 }, {ETA,  0.000000E0 }, {MU,  0.000000E0 }}},
    {27, {{ALPHA,  2.752869E0 }, {BETA,  2.729201E-2}, {GAMMA,  4.996927E-1}, {ETA,  2.496551E-1}, {MU,  0.000000E0 }}},
    {28, {{ALPHA,  3.518506E0 }, {BETA,  1.112440E0 }, {GAMMA, -4.556216E-1}, {ETA, -2.179426E-1}, {MU,  0.000000E0 }}},
    {29, {{ALPHA,  1.626062E2 }, {BETA, -1.168838E1 }, {GAMMA, -5.498343E0 }, {ETA,  0.000000E0 }, {MU,  0.000000E0 }}},
    {30, {{ALPHA,  3.336833E-1}, {BETA, -1.458043E-1}, {GAMMA, -2.011751E-2}, {ETA,  0.000000E0 }, {MU,  0.000000E0 }}},

    {31, {{ALPHA,  7.425137E1 }, {BETA,  1.790236E1 }, {GAMMA,  3.033910E1 }, {ETA,  1.018259E1 }, {MU,  0.000000E0 }}},
    {32, {{ALPHA,  9.268325E2 }, {BETA, -9.739859E1 }, {GAMMA, -7.702152E1 }, {ETA, -3.158268E1 }, {MU,  0.000000E0 }}},
    {33, {{ALPHA,  2.474401E0 }, {BETA,  3.892972E-1}, {GAMMA,  0.000000E0 }, {ETA,  0.000000E0 }, {MU,  0.000000E0 }}},
    {34, {{ALPHA,  1.127018E1 }, {BETA,  1.622158E0 }, {GAMMA, -1.443664E0 }, {ETA, -9.474699E-1}, {MU,  0.000000E0 }}},
    {35, {{ALPHA,  0.000000E0 }, {BETA,  0.000000E0 }, {GAMMA,  0.000000E0 }, {ETA,  0.000000E0 }, {MU,  0.000000E0 }}},
    {36, {{ALPHA,  1.445216E-1}, {BETA, -6.180219E-2}, {GAMMA,  3.093878E-2}, {ETA,  1.567090E-2}, {MU,  0.000000E0 }}},
    {37, {{ALPHA,  1.304129E0 }, {BETA,  1.395919E-1}, {GAMMA,  4.142455E-3}, {ETA, -9.732503E-3}, {MU,  0.000000E0 }}},
    {38, {{ALPHA,  5.114149E-1}, {BETA, -1.160850E-2}, {GAMMA,  0.000000E0 }, {ETA,  0.000000E0 }, {MU,  0.000000E0 }}},
    {39, {{ALPHA,  1.314955E2 }, {BETA,  2.009258E1 }, {GAMMA, -5.143082E-1}, {ETA, -1.379140E0 }, {MU,  0.000000E0 }}},
    {40, {{ALPHA,  1.823973E1 }, {BETA, -3.074559E0 }, {GAMMA, -4.307878E0 }, {ETA,  0.000000E0 }, {MU,  0.000000E0 }}},

    {41, {{ALPHA,  2.327037E0 }, {BETA,  2.403445E0 }, {GAMMA,  1.208407E0 }, {ETA,  2.087263E-1}, {MU,  0.000000E0 }}},
    {42, {{ALPHA,  1.997378E0 }, {BETA, -8.126205E-1}, {GAMMA,  0.000000E0 }, {ETA,  0.000000E0 }, {MU,  0.000000E0 }}},
    {43, {{ALPHA,  1.079113E-1}, {BETA,  1.762409E-2}, {GAMMA,  1.096601E-2}, {ETA,  3.058818E-3}, {MU,  0.000000E0 }}},
    {44, {{ALPHA,  2.327409E0 }, {BETA,  6.901582E-1}, {GAMMA, -2.158431E-1}, {ETA, -1.084117E-1}, {MU,  0.000000E0 }}},
    {45, {{ALPHA,  0.000000E0 }, {BETA,  0.000000E0 }, {GAMMA,  0.000000E0 }, {ETA,  0.000000E0 }, {MU,  0.000000E0 }}},
    {46, {{ALPHA,  2.214315E0 }, {BETA, -1.975747E0 }, {GAMMA,  0.000000E0 }, {ETA,  0.000000E0 }, {MU,  0.000000E0 }}},
    {47, {{ALPHA,  0.000000E0 }, {BETA,  0.000000E0 }, {GAMMA,  0.000000E0 }, {ETA,  0.000000E0 }, {MU,  0.000000E0 }}},
    {48, {{ALPHA,  5.072525E0 }, {BETA,  1.146189E1 }, {GAMMA,  6.961724E0 }, {ETA,  1.316965E0 }, {MU,  0.000000E0 }}},
    {49, {{ALPHA,  5.139740E0 }, {BETA,  0.000000E0 }, {GAMMA,  0.000000E0 }, {ETA,  0.000000E0 }, {MU,  0.000000E0 }}},
    {50, {{ALPHA,  0.000000E0 }, {BETA,  0.000000E0 }, {GAMMA,  0.000000E0 }, {ETA,  0.000000E0 }, {MU,  0.000000E0 }}},

    {51, {{ALPHA,  1.125124E0 }, {BETA,  1.306486E0 }, {GAMMA,  3.622359E0 }, {ETA,  2.601976E0 }, {MU,  3.031270E-1}}},
    {52, {{ALPHA,  3.349489E-1}, {BETA,  4.531269E-3}, {GAMMA,  1.131793E-1}, {ETA,  2.300156E-1}, {MU,  7.632745E-2}}},
    {53, {{ALPHA,  1.467794E0 }, {BETA,  2.798142E0 }, {GAMMA,  9.455580E0 }, {ETA,  8.963904E0 }, {MU,  3.339719E0 }}},
    {54, {{ALPHA,  4.658512E-1}, {BETA,  2.597451E-1}, {GAMMA,  9.048179E-1}, {ETA,  7.394505E-1}, {MU,  1.607092E-1}}},
    {55, {{ALPHA,  1.042200E0 }, {BETA,  1.315600E-1}, {GAMMA,  4.500000E-2}, {ETA,  0.000000E0 }, {MU,  0.000000E0 }}},
    {56, {{ALPHA,  1.110866E0 }, {BETA,  9.623856E-1}, {GAMMA,  2.735487E0 }, {ETA,  2.445602E0 }, {MU,  8.826352E-1}}},
    {57, {{ALPHA, -1.584333E-1}, {BETA, -1.728865E-1}, {GAMMA, -4.461431E-1}, {ETA, -3.925259E-1}, {MU, -1.276203E-1}}}
};

#undef ALPHA
#undef BETA
#undef GAMMA
#undef ETA
#undef MU

// C coefficients
// Key to map is n (C(n)).  Map element is unordered_map of term coefficient values.
const std::unordered_map<int, double> C_COEFF = {{1, -8.672073E-2}, {2, 9.301992E0}, {3, 4.637345E0}};

// CDF from Table 7 in Dufton et al 2013 https://arxiv.org/abs/1212.2424
// There is an assumption in the code that this function is monitonically increasing - it is now
// and should remain so if the map is modified.
const std::map<double, double> BStarRotationalVelocityCDFTable = {
    {000.0, 0.000}, {020.0, 0.046}, {040.0, 0.094}, {060.0, 0.144}, {080.0, 0.192}, {100.0, 0.239},
    {120.0, 0.253}, {140.0, 0.270}, {160.0, 0.288}, {180.0, 0.322}, {200.0, 0.377}, {220.0, 0.435},
    {240.0, 0.492}, {260.0, 0.548}, {280.0, 0.609}, {300.0, 0.674}, {320.0, 0.739}, {340.0, 0.796},
    {360.0, 0.841}, {380.0, 0.879}, {400.0, 0.912}, {420.0, 0.938}, {440.0, 0.956}, {460.0, 0.971},
    {480.0, 0.983}, {500.0, 0.990}, {520.0, 0.993}, {540.0, 0.995}, {560.0, 0.996}, {580.0, 0.997}
};

// These neutron star (NS) equations-of-state (EOS) are taken from the review Ozel & Freire 2016,
// Masses, Radii, and Equation of State of Neutron Stars,
// Annual Reviews of Astronomy and Astrophysics,
// https://arxiv.org/abs/1603.02698, downloaded from
// their website http://xtreme.as.arizona.edu/NeutronStars/
//
// for now we choose one example EOS ARP3 from
// Akmal et al 1998 https://arxiv.org/abs/nucl-th/9804027
const std::map<double, double> ARP3MassRadiusRelation = {
    {0.184 , 16.518}, {0.188 , 16.292}, {0.192 , 16.067}, {0.195 , 15.857}, {0.199 , 15.658}, {0.203 , 15.46 }, {0.207 , 15.277}, {0.212, 15.102}, {0.216, 14.933},
    {0.221 , 14.774}, {0.225 , 14.619}, {0.23  , 14.473}, {0.235 , 14.334}, {0.24  , 14.199}, {0.245 , 14.073}, {0.251 , 13.951}, {0.256, 13.834}, {0.262, 13.725},
    {0.268 , 13.618}, {0.273 , 13.52 }, {0.28  , 13.423}, {0.286 , 13.332}, {0.292 , 13.245}, {0.299 , 13.162}, {0.306 , 13.084}, {0.313, 13.009}, {0.32 , 12.94 },
    {0.327 , 12.871}, {0.335 , 12.806}, {0.342 , 12.747}, {0.35  , 12.691}, {0.358 , 12.638}, {0.366 , 12.586}, {0.374 , 12.538}, {0.383, 12.493}, {0.391, 12.451},
    {0.4   , 12.409}, {0.409 , 12.371}, {0.418 , 12.336}, {0.427 , 12.302}, {0.438 , 12.269}, {0.448 , 12.239}, {0.458 , 12.211}, {0.468, 12.184}, {0.479, 12.16 },
    {0.49  , 12.136}, {0.501 , 12.116}, {0.512 , 12.096}, {0.524 , 12.078}, {0.535 , 12.061}, {0.547 , 12.046}, {0.559 , 12.031}, {0.572, 12.018}, {0.585, 12.007},
    {0.598 , 11.997}, {0.611 , 11.987}, {0.625 , 11.979}, {0.638 , 11.972}, {0.652 , 11.966}, {0.666 , 11.96 }, {0.681 , 11.955}, {0.695, 11.952}, {0.71 , 11.949},
    {0.725 , 11.947}, {0.74  , 11.946}, {0.756 , 11.945}, {0.772 , 11.945}, {0.788 , 11.945}, {0.804 , 11.946}, {0.82  , 11.947}, {0.837, 11.949}, {0.854, 11.952},
    {0.871 , 11.955}, {0.888 , 11.957}, {0.906 , 11.961}, {0.923 , 11.964}, {0.941 , 11.968}, {0.959 , 11.972}, {0.977 , 11.977}, {0.995, 11.981}, {1.014, 11.985},
    {1.032 , 11.99 }, {1.05  , 11.994}, {1.069 , 11.999}, {1.088 , 12.004}, {1.107 , 12.009}, {1.126 , 12.013}, {1.145 , 12.018}, {1.164, 12.022}, {1.184, 12.027},
    {1.203 , 12.031}, {1.222 , 12.035}, {1.242 , 12.039}, {1.261 , 12.043}, {1.281 , 12.047}, {1.3   , 12.05 }, {1.32  , 12.053}, {1.339, 12.056}, {1.358, 12.058},
    {1.378 , 12.061}, {1.397 , 12.063}, {1.416 , 12.064}, {1.436 , 12.066}, {1.455 , 12.067}, {1.474 , 12.068}, {1.493 , 12.068}, {1.512, 12.068}, {1.531, 12.068},
    {1.549 , 12.067}, {1.568 , 12.066}, {1.586 , 12.065}, {1.604 , 12.063}, {1.623 , 12.06 }, {1.64  , 12.058}, {1.658 , 12.055}, {1.676, 12.052}, {1.693, 12.048},
    {1.71  , 12.044}, {1.727 , 12.039}, {1.744 , 12.034}, {1.761 , 12.029}, {1.777 , 12.024}, {1.793 , 12.017}, {1.809 , 12.011}, {1.825, 12.004}, {1.84 , 11.997},
    {1.856 , 11.989}, {1.871 , 11.981}, {1.886 , 11.973}, {1.9   , 11.965}, {1.915 , 11.956}, {1.929 , 11.946}, {1.943 , 11.937}, {1.956, 11.927}, {1.969, 11.916},
    {1.982 , 11.906}, {1.995 , 11.895}, {2.008 , 11.884}, {2.02  , 11.827}, {2.032 , 11.86 }, {2.044 , 11.848}, {2.056 , 11.836}, {2.067, 11.823}, {2.078, 11.81 },
    {2.089 , 11.797}, {2.099 , 11.784}, {2.109 , 11.77 }, {2.119 , 11.756}, {2.129 , 11.742}, {2.139 , 11.727}, {2.148 , 11.713}, {2.157, 11.698}, {2.166, 11.683},
    {2.174 , 11.668}, {2.182 , 11.652}, {2.19  , 11.637}, {2.198 , 11.621}, {2.206 , 11.605}, {2.213 , 11.589}, {2.221 , 11.573}, {2.227, 11.556}, {2.234, 11.54 },
    {2.241 , 11.523}, {2.247 , 11.506}, {2.253 , 11.49 }, {2.259 , 11.473}, {2.264 , 11.456}, {2.27  , 11.438}, {2.275 , 11.421}, {2.28 , 11.404}, {2.285, 11.386},
    {2.29  , 11.369}, {2.294 , 11.351}, {2.299 , 11.333}, {2.303 , 11.316}, {2.307 , 11.299}, {2.31  , 11.281}, {2.314 , 11.263}, {2.317, 11.245}, {2.321, 11.227},
    {2.324 , 11.209}, {2.327 , 11.191}, {2.33  , 11.173}, {2.332 , 11.155}, {2.335 , 11.136}, {2.337 , 11.119}, {2.339 , 11.101}, {2.342, 11.083}, {2.344, 11.065},
    {2.345 , 11.046}, {2.347 , 11.028}, {2.349 , 11.01 }, {2.35  , 10.992}, {2.352 , 10.974}, {2.353 , 10.956}, {2.354 , 10.938}, {2.355, 10.92 }, {2.356, 10.902},
    {2.3571, 10.885}, {2.3572, 10.866}, {2.3581, 10.849}, {2.3582, 10.831}, {2.3591, 10.813}, {2.3592, 10.795}, {2.3593, 10.777}, {2.361, 10.76 }, {2.362, 10.742}
};

// mass / Msun of stellar models computed by Xu & Li (2010), with additional unpublished 50 and 100 Msun models
const DBL_VECTOR NANJING_MASSES = { 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 12.0, 14.0, 16.0, 20.0, 50.0, 100.0 };

// mass / Msun bin edges of Xu & Li (2010) lambda prescription as implemented in StarTrack. These are the midpoints between the masses in NANJING_MASSES
const DBL_VECTOR NANJING_MASSES_MIDPOINTS = { 1.5, 2.5, 3.5, 4.5, 5.5, 6.5, 7.5, 8.5, 9.5, 11.0, 13.0, 15.0, 18.0, 35.0, 75.0 };

// coefficients for calculating binding and recombination energy as described in Loveridge et al. 2011
// electronic tables, program and further information in: http://astro.ru.nl/~sluys/index.php?title=BE

// struct LoveridgeCoefficients
// m, r, alpha(m,r) used in Loveridge et al. 2011, eq 5
struct LoveridgeCoefficients {
    int    m;
    int    r;
    double alpha_mr;
};

// symbolic names for metallicities described in Loveridge et al., 2011
// These are used as indices into the loveridgeCoefficients multi-dimensional vector (described below)
// This is a bit of a hack until I figure out how to (elegantly) iterate over enum classes...
enum class LOVERIDGE_METALLICITY: int { Z000010, Z000100, Z001000, Z001500, Z002000, Z003000, COUNT };
const std::vector<std::tuple<LOVERIDGE_METALLICITY, double>> LOVERIDGE_METALLICITY_VALUE = {
    { LOVERIDGE_METALLICITY::Z000010, 0.00010 },
    { LOVERIDGE_METALLICITY::Z000100, 0.00100 },
    { LOVERIDGE_METALLICITY::Z001000, 0.01000 },
    { LOVERIDGE_METALLICITY::Z001500, 0.01500 },
    { LOVERIDGE_METALLICITY::Z002000, 0.02000 },
    { LOVERIDGE_METALLICITY::Z003000, 0.03000 }
};

// values for the division between High Mass and Low Mass, indexed by metallicity (LOVERIDGE_METALLICITY)
// Values are given in Msol
// From Loveridge et al. 2011, table 1
const DBL_VECTOR LOVERIDGE_LM_HM_CUTOFFS = { 11.7000, 11.7000, 10.2000, 11.7000, 11.7000, 13.4000 };

// coefficients for the division between Low Mass RGB 1 and Low Mass RGB 2, indexed by metallicity (LOVERIDGE_METALLICITY)
// From Loveridge et al. 2011, table 2 & eq 4
const std::vector<DBL_VECTOR> LOVERIDGE_LM1_LM2_CUTOFFS = {
    { 6.047230E-01,  1.397240E+00, -8.249630E-02,  1.141430E+00,  0.000000E+00 },       // Metallicity Z000010 (0.00010)
    { 4.566000E-01,  1.186600E+00,  2.390880E+00, -3.054040E+00,  1.404900E+00 },       // Metallicity Z000100 (0.00100)
    { 2.821740E-01,  1.149380E+00,  1.884450E+00, -1.082300E+00,  0.000000E+00 },       // Metallicity Z001000 (0.01000)
    { 2.518180E-01,  1.210490E+00,  1.634900E+00, -8.369100E-01,  0.000000E+00 },       // Metallicity Z001500 (0.01500)
    { 2.406370E-01,  1.089220E+00,  1.953180E+00, -1.032180E+00,  0.000000E+00 },       // Metallicity Z002000 (0.02000)
    { 2.348880E-01,  8.972940E-01,  2.519950E+00, -1.414110E+00,  0.000000E+00 }        // Metallicity Z003000 (0.03000)
};

// multi-dimensional vector of struct LoveridgeCoefficients, indexed by
// metallicity (LOVERIDGE_METALLICITY) and evolutionary stage (LOVERIDGE_GROUP)
// This vector records the coefficients (m, r, alpha(m,r) used in Loveridge et al. 2001, eq 5
// Using a vector indexed by metallicity and evolutionary stage because it's faster than a map, and
// the lambdas could be calculated at every timestep.
const std::vector<std::vector<std::vector<LoveridgeCoefficients>>> LOVERIDGE_COEFFICIENTS = {
    {                                                                                   // Metallicity Z000010 (0.00010)
        {                                                                               // LMR1 (Z000010)
            { 0,     0,     1.49884369566236408389E+01},
            { 0,     1,     3.55674019216888570583E+00},
            { 0,     2,    -1.50579325323499482181E+01},
            { 0,     3,     2.74507278637946647848E+01},
            { 0,     4,    -2.40420132204742671433E+01},
            { 0,     5,     6.08559902401751795509E+00},

            { 1,     0,     4.77689517615753889146E+00},
            { 1,     1,    -3.52448257631879471319E+01},
            { 1,     2,     1.26181509166749165729E+02},
            { 1,     3,    -2.05760448139415075275E+02},
            { 1,     4,     1.66826460262414656199E+02},
            { 1,     5,    -4.13681688626000720888E+01},

            { 2,     0,    -4.03038048915965774199E+00},
            { 2,     1,     5.94984473818373302834E+01},
            { 2,     2,    -3.15891070807037635859E+02},
            { 2,     3,     5.63465810280387586317E+02},
            { 2,     4,    -4.56401018552895436642E+02},
            { 2,     5,     1.12357256780075772440E+02},

            { 3,     0,    -2.22258127619659440199E+00},
            { 3,     1,     1.21818087660567726971E+02},
            { 3,     2,    -6.06793401690043339158E+01},
            { 3,     3,    -2.73750098046708558286E+02},
            { 3,     4,     3.80968083978251002009E+02},
            { 3,     5,    -1.07210865421446229107E+02},

            { 4,     0,    -4.95448396548353997559E+01},
            { 4,     1,    -1.05676079281290000722E+02},
            { 4,     2,     3.74254532751612941865E+02},
            { 4,     3,    -2.84755814237885886087E+02},
            { 4,     4,     5.32060692031168436245E+00},
            { 4,     5,     1.94841031059088862776E+01},

            { 5,     0,     5.85417149247924371025E+01},
            { 5,     1,    -6.46713025038344397899E+01},
            { 5,     2,    -8.04514035300949643670E+01},
            { 5,     3,     1.54013846765123219029E+02},
            { 5,     4,    -6.62783052076742649206E+01},
            { 5,     5,     9.83910595056972248074E+00}
        },
        {                                                                               // LMR2 (Z000010)
            { 0,     0,     2.10206064943832124925E+01},
            { 0,     1,    -2.39940628010456791230E+01},
            { 0,     2,     3.67437434259622861532E+01},
            { 0,     3,    -2.87504026348741277275E+01},
            { 0,     4,     1.10696952815601967757E+01},
            { 0,     5,    -1.67354101724841819454E+00},

            { 1,     0,     6.24726695402092602194E+01},
            { 1,     1,    -2.25859701401090774198E+02},
            { 1,     2,     3.25693445380178616233E+02},
            { 1,     3,    -2.28906270354160255920E+02},
            { 1,     4,     7.82835291167177160787E+01},
            { 1,     5,    -1.04409269263635096081E+01},

            { 2,     0,     1.68774936141528343114E+02},
            { 2,     1,    -4.70922534725343496120E+02},
            { 2,     2,     5.20150477052292671942E+02},
            { 2,     3,    -2.82942436111233064366E+02},
            { 2,     4,     7.54607477257930696624E+01},
            { 2,     5,    -7.80062541052705249456E+00},

            { 3,     0,     1.26323501968766254322E+03},
            { 3,     1,    -5.43724065618109580100E+03},
            { 3,     2,     9.47031538171058127773E+03},
            { 3,     3,    -8.20344328990647773026E+03},
            { 3,     4,     3.48253888526251739677E+03},
            { 3,     5,    -5.75361752664876235031E+02},

            { 4,     0,     1.45320316532362594444E+04},
            { 4,     1,    -6.10692503818239565589E+04},
            { 4,     2,     9.45752483181984280236E+04},
            { 4,     3,    -6.92033750093292765087E+04},
            { 4,     4,     2.43234260768021413242E+04},
            { 4,     5,    -3.32540856427475091550E+03},

            { 5,     0,    -7.83727239733487567719E+03},
            { 5,     1,     6.87101874631883547409E+04},
            { 5,     2,    -1.42788737041162559763E+05},
            { 5,     3,     1.25369407779255416244E+05},
            { 5,     4,    -5.05985607497797464021E+04},
            { 5,     5,     7.77505329663658358186E+03}
        },
        {                                                                               // LMA (Z000010)
            { 0,     0,     1.77846204423370872973E+04},
            { 0,     1,    -1.11008631122171675088E+05},
            { 0,     2,     3.07385212080689030699E+05},
            { 0,     3,    -4.97253519789625774138E+05},
            { 0,     4,     5.20899845929651521146E+05},
            { 0,     5,    -3.69562230436008889228E+05},
            { 0,     6,     1.79995285036839370150E+05},
            { 0,     7,    -5.94766776453754428076E+04},
            { 0,     8,     1.27704226161695205519E+04},
            { 0,     9,    -1.60998313917297696207E+03},
            { 0,    10,     9.05540938508377593053E+01},

            { 1,     0,     1.27448576992469941615E+05},
            { 1,     1,    -8.26126519162579439580E+05},
            { 1,     2,     2.40616127883097669110E+06},
            { 1,     3,    -4.13147857106406055391E+06},
            { 1,     4,     4.61530595326700154692E+06},
            { 1,     5,    -3.49471896526116924360E+06},
            { 1,     6,     1.81238210394326946698E+06},
            { 1,     7,    -6.34629941238884348422E+05},
            { 1,     8,     1.43452093324876565021E+05},
            { 1,     9,    -1.88911385051759207272E+04},
            { 1,    10,     1.10039680760221062883E+03},

            { 2,     0,    -1.11280780545824207366E+06},
            { 2,     1,     6.52773804973363596946E+06},
            { 2,     2,    -1.68385778483916968107E+07},
            { 2,     3,     2.51369743430132977664E+07},
            { 2,     4,    -2.40291083278050050139E+07},
            { 2,     5,     1.53511958359827548265E+07},
            { 2,     6,    -6.62599811194045469165E+06},
            { 2,     7,     1.90266653405042551458E+06},
            { 2,     8,    -3.46290151645659178030E+05},
            { 2,     9,     3.57968178594517958118E+04},
            { 2,    10,    -1.57403299302352661471E+03},

            { 3,     0,     3.61468220664994791150E+06},
            { 3,     1,    -1.42838660574260130525E+07},
            { 3,     2,     9.89719916261141002178E+06},
            { 3,     3,     4.17630757517836764455E+07},
            { 3,     4,    -1.13186791305614486337E+08},
            { 3,     5,     1.34723130784819722176E+08},
            { 3,     6,    -9.39780759352868050337E+07},
            { 3,     7,     4.08028988015334084630E+07},
            { 3,     8,    -1.08834827876409199089E+07},
            { 3,     9,     1.63703844878768618219E+06},
            { 3,    10,    -1.06502153903636033647E+05},

            { 4,     0,    -2.18383920460389368236E+07},
            { 4,     1,     1.01377685264262214303E+08},
            { 4,     2,    -1.24736986550756111741E+08},
            { 4,     3,    -1.38097782211961090565E+08},
            { 4,     4,     5.82118970734395384789E+08},
            { 4,     5,    -7.72188668410225749016E+08},
            { 4,     6,     5.69788365736976385117E+08},
            { 4,     7,    -2.56651440166880398989E+08},
            { 4,     8,     7.03184175257203429937E+07},
            { 4,     9,    -1.07993168413460906595E+07},
            { 4,    10,     7.14464107997456681915E+05},

            { 5,     0,     6.22013266083969771862E+07},
            { 5,     1,    -3.70892432569035887718E+08},
            { 5,     2,     7.32722076455112814903E+08},
            { 5,     3,    -3.41863748162672758102E+08},
            { 5,     4,    -8.72008743590860724449E+08},
            { 5,     5,     1.70439955004952502251E+09},
            { 5,     6,    -1.44412099096650505066E+09},
            { 5,     7,     7.01467443390604257584E+08},
            { 5,     8,    -2.01846242185972064734E+08},
            { 5,     9,     3.21032091475058645010E+07},
            { 5,    10,    -2.18098983308966364712E+06},

            { 6,     0,     8.32790288301304075867E+06},
            { 6,     1,     1.67320728489836782217E+08},
            { 6,     2,    -4.90329729719172537327E+08},
            { 6,     3,    -3.57015713222805708647E+07},
            { 6,     4,     1.73606751974490427971E+09},
            { 6,     5,    -2.95511773960293579102E+09},
            { 6,     6,     2.49140620948586368561E+09},
            { 6,     7,    -1.22675662513774180412E+09},
            { 6,     8,     3.58779851358991682529E+08},
            { 6,     9,    -5.79478609330464825034E+07},
            { 6,    10,     3.99136670739177288488E+06},

            { 7,     0,    -3.09375949266583919525E+08},
            { 7,     1,     1.30927392519327545166E+09},
            { 7,     2,    -2.71972201258040809631E+09},
            { 7,     3,     4.17056345501154565811E+09},
            { 7,     4,    -5.31121138472141742706E+09},
            { 7,     5,     5.18211883778930091858E+09},
            { 7,     6,    -3.52386318383868503571E+09},
            { 7,     7,     1.57889703104470300674E+09},
            { 7,     8,    -4.41837538483527064323E+08},
            { 7,     9,     6.98762237560149580240E+07},
            { 7,    10,    -4.76660235680679418147E+06},

            { 8,     0,     5.32124163954892635345E+08},
            { 8,     1,    -2.58332422589960527420E+09},
            { 8,     2,     5.78740993511894130707E+09},
            { 8,     3,    -8.18639627587050056458E+09},
            { 8,     4,     8.33603336255734443665E+09},
            { 8,     5,    -6.37392318348361968994E+09},
            { 8,     6,     3.59443787565530967712E+09},
            { 8,     7,    -1.42294472536891078949E+09},
            { 8,     8,     3.68482395798513412476E+08},
            { 8,     9,    -5.55465651649229973555E+07},
            { 8,    10,     3.67803603909488115460E+06},

            { 9,     0,    -3.67132831968976199627E+08},
            { 9,     1,     1.87202460446496915817E+09},
            { 9,     2,    -4.32871663092040729523E+09},
            { 9,     3,     6.09032624990053653717E+09},
            { 9,     4,    -5.87735082626143074036E+09},
            { 9,     5,     4.10221144979333591461E+09},
            { 9,     6,    -2.08821944738185024261E+09},
            { 9,     7,     7.54238944109313964844E+08},
            { 9,     8,    -1.81710771632690608501E+08},
            { 9,     9,     2.59757891273681670427E+07},
            { 9,    10,    -1.65624886172299087048E+06},

            {10,     0,     9.30481611941920518875E+07},
            {10,     1,    -4.87096365846687316895E+08},
            {10,     2,     1.14762701275433015823E+09},
            {10,     3,    -1.62114284406450057030E+09},
            {10,     4,     1.53847292382848095894E+09},
            {10,     5,    -1.03373922282777369022E+09},
            {10,     6,     4.99196402871803343296E+08},
            {10,     7,    -1.70266350400250434875E+08},
            {10,     8,     3.88883041814338341355E+07},
            {10,     9,    -5.31391571130599640310E+06},
            {10,    10,     3.26864794701321807224E+05}
        },
        {                                                                               // HM (Z000010)
            { 0,     0,     3.04240080558681453113E+05},
            { 0,     1,    -7.24150360611511170864E+07},
            { 0,     2,     2.70646920133822739124E+08},
            { 0,     3,    -4.33003387311595380306E+08},
            { 0,     4,     4.03091375320128858089E+08},
            { 0,     5,    -2.49456192293578892946E+08},
            { 0,     6,     9.70010504787636250257E+07},
            { 0,     7,    -8.23838906384931178764E+05},
            { 0,     8,    -3.28722313181499056518E+07},
            { 0,     9,     2.78880451327459998429E+07},
            { 0,    10,    -1.40338722991762422025E+07},
            { 0,    11,     5.06007529781199898571E+06},
            { 0,    12,    -1.37089679529935028404E+06},
            { 0,    13,     2.68858224094756355044E+05},
            { 0,    14,    -3.36690685426766431192E+04},
            { 0,    15,     1.98948270157709225714E+03},

            { 1,     0,     5.37679127182313203812E+07},
            { 1,     1,     1.46388718178314149380E+08},
            { 1,     2,    -9.35395212179676651955E+08},
            { 1,     3,     1.55915548629785561562E+09},
            { 1,     4,    -1.34035874177675724030E+09},
            { 1,     5,     7.60583187397939205170E+08},
            { 1,     6,    -3.48267020550069510937E+08},
            { 1,     7,     1.27565538716441661119E+08},
            { 1,     8,    -2.52623368581141643226E+07},
            { 1,     9,    -2.38211218436619313434E+06},
            { 1,    10,     2.83049418278661277145E+06},
            { 1,    11,    -1.12351439162609330378E+06},
            { 1,    12,     5.12258502489926875569E+05},
            { 1,    13,    -1.87547968379564146744E+05},
            { 1,    14,     3.78142640863094275119E+04},
            { 1,    15,    -3.14819542193592815238E+03},

            { 2,     0,    -2.72150296501208603382E+08},
            { 2,     1,     2.54191144958757966757E+08},
            { 2,     2,     1.09221455545737147331E+09},
            { 2,     3,    -2.27038845949890232086E+09},
            { 2,     4,     1.74805773803573489189E+09},
            { 2,     5,    -6.91637455687782287598E+08},
            { 2,     6,     2.04606336638418078423E+08},
            { 2,     7,    -7.53842029037321358919E+07},
            { 2,     8,     1.98339755878359973431E+07},
            { 2,     9,    -2.51659290353266568854E+06},
            { 2,    10,     2.78171596346601564437E+06},
            { 2,    11,    -2.17292116712976712734E+06},
            { 2,    12,     6.62271515154829132371E+05},
            { 2,    13,    -6.85233401404126780108E+04},
            { 2,    14,    -7.24533940448010434920E+03},
            { 2,    15,     1.77091708853070326768E+03},

            { 3,     0,     5.84549742654761433601E+08},
            { 3,     1,    -1.22755266303496289253E+09},
            { 3,     2,    -1.99841301217404045165E+07},
            { 3,     3,     1.62942775427007746696E+09},
            { 3,     4,    -1.29122969498875904083E+09},
            { 3,     5,     2.81517572141482591629E+08},
            { 3,     6,     3.29046875321928299963E+07},
            { 3,     7,    -2.18449690962994331494E+06},
            { 3,     8,    -2.63601564748177304864E+06},
            { 3,     9,    -4.48942869859682396054E+06},
            { 3,    10,     2.44109624343642685562E+06},
            { 3,    11,    -4.88764745659762818832E+05},
            { 3,    12,     5.56227389265404126490E+04},
            { 3,    13,    -2.70979608753431966761E+03},
            { 3,    14,     1.47652876079566203771E+03},
            { 3,    15,    -5.79993594759073516798E+02},

            { 4,     0,    -6.80110049336398720741E+08},
            { 4,     1,     1.85585044735255384445E+09},
            { 4,     2,    -1.28215713107603287697E+09},
            { 4,     3,    -3.01044179228078424931E+08},
            { 4,     4,     5.51112592108273148537E+08},
            { 4,     5,    -7.78949714611403048038E+07},
            { 4,     6,    -3.93153402725159674883E+07},
            { 4,     7,    -1.19108830926259383559E+07},
            { 4,     8,     1.43100595134035050869E+07},
            { 4,     9,    -4.91528834089644718915E+06},
            { 4,    10,     1.68842217521909135394E+06},
            { 4,    11,    -1.75640131265274452744E+05},
            { 4,    12,    -1.40626545874571602326E+05},
            { 4,    13,     3.77159267968944986933E+04},
            { 4,    14,    -1.39635312769678466793E+03},
            { 4,    15,    -1.88574192421971709166E+01},

            { 5,     0,     4.31829985381816327572E+08},
            { 5,     1,    -1.41004882914948272705E+09},
            { 5,     2,     1.38554457511777806282E+09},
            { 5,     3,    -3.82333992940614879131E+08},
            { 5,     4,    -7.42225503488553911448E+07},
            { 5,     5,     1.56423537414295095950E+07},
            { 5,     6,     1.26731311307457205839E+06},
            { 5,     7,     8.38170131250077579170E+06},
            { 5,     8,     1.69866128198041976430E+06},
            { 5,     9,    -2.63807965752872545272E+06},
            { 5,    10,    -1.30410926100970245898E+05},
            { 5,    11,     2.79376974434808362275E+05},
            { 5,    12,    -2.47335903090722194975E+03},
            { 5,    13,    -8.09120832153944047604E+03},
            { 5,    14,    -1.18437059561102614680E+03},
            { 5,    15,     2.22406017739857361448E+02},

            { 6,     0,    -1.05728919394894704223E+08},
            { 6,     1,     4.78614212386758863926E+08},
            { 6,     2,    -5.47275303584568023682E+08},
            { 6,     3,     1.81002869935932725668E+08},
            { 6,     4,     2.37412980302353836596E+07},
            { 6,     5,    -1.49761239025783911347E+07},
            { 6,     6,     6.94548275970081239939E+06},
            { 6,     7,    -7.08765067577983532101E+06},
            { 6,     8,     6.15334952530000591651E+05},
            { 6,     9,     1.04307962681467283983E+06},
            { 6,    10,    -1.66822670094482309651E+05},
            { 6,    11,    -3.14983610603470224305E+04},
            { 6,    12,    -8.67230359049541766581E+03},
            { 6,    13,     3.35310141926492042330E+03},
            { 6,    14,     5.20457690967627286227E+02},
            { 6,    15,    -9.46886765048879937012E+01},

            { 7,     0,    -4.00641674232298433781E+07},
            { 7,     1,     3.53151559413958713412E+07},
            { 7,     2,    -3.32651288367051668465E+07},
            { 7,     3,     8.78685073832835257053E+07},
            { 7,     4,    -8.83196869217660874128E+07},
            { 7,     5,     3.53058029766001477838E+07},
            { 7,     6,    -2.53460763685103179887E+06},
            { 7,     7,    -3.27478091701123584062E+06},
            { 7,     8,     1.51544590128077217378E+06},
            { 7,     9,    -1.97898065117523074150E+05},
            { 7,    10,    -4.86252464679210970644E+04},
            { 7,    11,     3.00813760869401539821E+03},
            { 7,    12,     9.57301967925047574681E+03},
            { 7,    13,    -2.60470610160503747466E+03},
            { 7,    14,     1.71621252115364512747E+02},
            { 7,    15,    -4.79603274832587800347E+00},

            { 8,     0,     3.64528375028086900711E+07},
            { 8,     1,    -7.85125069899140596390E+07},
            { 8,     2,     7.79785535495323836803E+07},
            { 8,     3,    -5.67020461767953187227E+07},
            { 8,     4,     3.75835219815801903605E+07},
            { 8,     5,    -2.08943791429358497262E+07},
            { 8,     6,     7.66625923235765285790E+06},
            { 8,     7,    -1.06592743443907494657E+06},
            { 8,     8,    -3.14079768184977932833E+05},
            { 8,     9,     1.45893387703247601166E+05},
            { 8,    10,    -2.00583693587448651670E+04},
            { 8,    11,    -1.07287338880777724626E+03},
            { 8,    12,     1.86369211872262280849E+03},
            { 8,    13,    -5.41388893973812855620E+02},
            { 8,    14,     3.82187077215814881015E+01},
            { 8,    15,     4.76195893580739504358E+00},

            { 9,     0,    -9.92938956532538309693E+06},
            { 9,     1,     2.04695142809341102839E+07},
            { 9,     2,    -1.02214338619423378259E+07},
            { 9,     3,    -9.93348938946167938411E+06},
            { 9,     4,     1.64549253709395211190E+07},
            { 9,     5,    -9.86318830775812268257E+06},
            { 9,     6,     3.13746389613276068121E+06},
            { 9,     7,    -6.67069913326343754306E+05},
            { 9,     8,     2.16856145024514931720E+05},
            { 9,     9,    -1.08820304041782830609E+05},
            { 9,    10,     3.55738569172507050098E+04},
            { 9,    11,    -3.12188240990893382332E+03},
            { 9,    12,    -1.97305440808249090878E+03},
            { 9,    13,     7.01633632749421281005E+02},
            { 9,    14,    -8.03244389802174652004E+01},
            { 9,    15,     1.98516265048327333886E+00},

            {10,     0,     1.05106274467385723256E+06},
            {10,     1,    -1.90074146260506426916E+06},
            {10,     2,    -5.14866259700076945592E+05},
            {10,     3,     4.43767094488164782524E+06},
            {10,     4,    -5.41318389280360378325E+06},
            {10,     5,     3.27827314338852465153E+06},
            {10,     6,    -1.00315844098674086854E+06},
            {10,     7,     3.13072588745724533510E+04},
            {10,     8,     1.04236830211061693262E+05},
            {10,     9,    -4.78932228588938814937E+04},
            {10,    10,     1.31090834075820221187E+04},
            {10,    11,    -3.39479313977342053477E+03},
            {10,    12,     9.29207146820721618496E+02},
            {10,    13,    -1.86795988806825192796E+02},
            {10,    14,     1.98510067854014131683E+01},
            {10,    15,    -7.40917959154837491020E-01}
        },
        {                                                                               // RECOM (Z000010)
            { 0,     0,     1.43249838796439163957E+01},
            { 0,     1,    -1.34095062880606157307E+01},
            { 0,     2,     6.08550621146791144156E+01},
            { 0,     3,    -1.46577101984033276949E+02},
            { 0,     4,     2.10671492276514470632E+02},
            { 0,     5,    -1.92066044872448998149E+02},
            { 0,     6,     1.12528193269927839992E+02},
            { 0,     7,    -4.10571661215909387010E+01},
            { 0,     8,     8.47360078103331737509E+00},
            { 0,     9,    -7.53719632151268914555E-01},

            { 1,     0,    -4.69105628993968171159E+00},
            { 1,     1,     8.90187555533843095645E+01},
            { 1,     2,    -5.17082466752443224323E+02},
            { 1,     3,     1.49287044908038524227E+03},
            { 1,     4,    -2.42034351477910740869E+03},
            { 1,     5,     2.35977318225110639105E+03},
            { 1,     6,    -1.41747263752418666627E+03},
            { 1,     7,     5.14675357902906625895E+02},
            { 1,     8,    -1.03788126354235686222E+02},
            { 1,     9,     8.93875687034091015448E+00},

            { 2,     0,    -1.98086570192114912459E+01},
            { 2,     1,     9.02725480979547114657E+01},
            { 2,     2,     3.72962652775568585639E+02},
            { 2,     3,    -2.86443864935814417549E+03},
            { 2,     4,     6.54791371890192021965E+03},
            { 2,     5,    -7.57542293321982197085E+03},
            { 2,     6,     4.98863219830037542124E+03},
            { 2,     7,    -1.90003291421173662457E+03},
            { 2,     8,     3.91370534933327633098E+02},
            { 2,     9,    -3.38888734262211954729E+01},

            { 3,     0,    -7.61727485013823990556E+00},
            { 3,     1,    -8.99676221452425437519E+01},
            { 3,     2,     2.39429289955121817002E+02},
            { 3,     3,     2.05771307220595781473E+03},
            { 3,     4,    -7.66672108820838730026E+03},
            { 3,     5,     1.07610646452268938447E+04},
            { 3,     6,    -7.83052958018395565887E+03},
            { 3,     7,     3.14853537506878637942E+03},
            { 3,     8,    -6.66805677756379054699E+02},
            { 3,     9,     5.84483701602106719974E+01},

            { 4,     0,     6.84689749428262302899E+02},
            { 4,     1,    -6.41337423137322639377E+03},
            { 4,     2,     2.14353806903193471953E+04},
            { 4,     3,    -3.96806322513512059231E+04},
            { 4,     4,     4.68443039997480664169E+04},
            { 4,     5,    -3.66080543837553850608E+04},
            { 4,     6,     1.87025281715490273200E+04},
            { 4,     7,    -5.95491842750383239036E+03},
            { 4,     8,     1.06698781883290007499E+03},
            { 4,     9,    -8.22032529086119865269E+01},
            { 5,     0,    -1.08938704804548956417E+03},
            { 5,     1,     1.48146628563555968867E+04},
            { 5,     2,    -5.75561987341030326206E+04},
            { 5,     3,     1.11350161916194265359E+05},
            { 5,     4,    -1.26655054015904039261E+05},
            { 5,     5,     9.02577656921312009217E+04},
            { 5,     6,    -4.08346371174871092080E+04},
            { 5,     7,     1.13850781915498100716E+04},
            { 5,     8,    -1.78117510051617318823E+03},
            { 5,     9,     1.19593664869470444501E+02},

            { 6,     0,     7.71158520059956913428E+01},
            { 6,     1,    -1.14063205364834775537E+04},
            { 6,     2,     5.61224719137911379221E+04},
            { 6,     3,    -1.19398022502862586407E+05},
            { 6,     4,     1.41535388629036053317E+05},
            { 6,     5,    -1.01977152037983221817E+05},
            { 6,     6,     4.57325986781135798083E+04},
            { 6,     7,    -1.24671193954050559114E+04},
            { 6,     8,     1.88748088780042576218E+03},
            { 6,     9,    -1.21427200114438761602E+02},

            { 7,     0,     7.51859271237372013275E+02},
            { 7,     1,     2.60919535160936538887E+03},
            { 7,     2,    -2.43703099252169959072E+04},
            { 7,     3,     6.04402419973716023378E+04},
            { 7,     4,    -7.66051639442676823819E+04},
            { 7,     5,     5.70820587259594394709E+04},
            { 7,     6,    -2.60263704770619915507E+04},
            { 7,     7,     7.14161111890670417779E+03},
            { 7,     8,    -1.08138895268239821235E+03},
            { 7,     9,     6.92421978398739810245E+01},

            { 8,     0,    -4.88596016673188898949E+02},
            { 8,     1,     6.06324096291598948483E+02},
            { 8,     2,     4.18032368029227109218E+03},
            { 8,     3,    -1.41908964780536971375E+04},
            { 8,     4,     1.99120584653974874527E+04},
            { 8,     5,    -1.55743938039613858564E+04},
            { 8,     6,     7.29395282901822793065E+03},
            { 8,     7,    -2.03371603812545799883E+03},
            { 8,     8,     3.11077689089091336427E+02},
            { 8,     9,    -2.00523533590040692332E+01},

            { 9,     0,     9.23231057243035166948E+01},
            { 9,     1,    -2.49558777666313233112E+02},
            { 9,     2,    -9.91094586704329287841E+01},
            { 9,     3,     1.18411592281756065859E+03},
            { 9,     4,    -1.97229180010892332575E+03},
            { 9,     5,     1.65110538528935694558E+03},
            { 9,     6,    -8.01140239767901903178E+02},
            { 9,     7,     2.28187127507413094918E+02},
            { 9,     8,    -3.54063664898269081505E+01},
            { 9,     9,     2.30675599774240458473E+00}
        }
    },
    {                                                                                   // Metallicity Z000100 (0.00100)
        {                                                                               // LMR1 (Z000100)
            { 0,     0,     1.50346099257731538046E+01},
            { 0,     1,     2.82604496789430559289E+00},
            { 0,     2,    -9.81648621516211328242E+00},
            { 0,     3,     9.56667166719136474740E+00},
            { 0,     4,    -3.13477224724504077713E+00},

            { 1,     0,     2.07944042796325678779E+00},
            { 1,     1,    -5.06461279977239176020E+00},
            { 1,     2,     2.03394884719189121824E+01},
            { 1,     3,    -1.93735154009608017134E+01},
            { 1,     4,     5.75610929539014914980E+00},

            { 2,     0,    -1.05127268925442436398E+01},
            { 2,     1,     3.68098490829458810936E+01},
            { 2,     2,    -5.91474113982236531228E+01},
            { 2,     3,     3.16463783616280913691E+01},
            { 2,     4,    -3.65289167448958851381E+00},

            { 3,     0,     2.22020213350119419715E+01},
            { 3,     1,    -7.77940541900198496705E+01},
            { 3,     2,     1.05860178274410074550E+02},
            { 3,     3,    -4.94197227165773327329E+01},
            { 3,     4,     3.82693132506987110375E+00},

            { 4,     0,    -1.12988524726900649853E+01},
            { 4,     1,     4.13367947020142807446E+01},
            { 4,     2,    -5.60688828217752686101E+01},
            { 4,     3,     2.71194601089813360772E+01},
            { 4,     4,    -2.72130850772848109642E+00}
        },
        {                                                                               // LMR2 (Z000100)
            { 0,     0,     1.88453218130641353412E+01},
            { 0,     1,    -1.66303719796976672285E+01},
            { 0,     2,     2.63222496911415824172E+01},
            { 0,     3,    -2.13131516790098238801E+01},
            { 0,     4,     8.39818275012596870965E+00},
            { 0,     5,    -1.28795907731821102082E+00},

            { 1,     0,     4.31329741689262391446E+01},
            { 1,     1,    -1.56619568143149109574E+02},
            { 1,     2,     2.26032050176760577642E+02},
            { 1,     3,    -1.59897885791484213769E+02},
            { 1,     4,     5.55606561823480475937E+01},
            { 1,     5,    -7.60112702919777394328E+00},

            { 2,     0,     2.46987072824178248709E+01},
            { 2,     1,     4.20629378936399973554E+01},
            { 2,     2,    -2.65767178117826063044E+02},
            { 2,     3,     3.25730032378355303990E+02},
            { 2,     4,    -1.58371785651665220485E+02},
            { 2,     5,     2.75938397603546903269E+01},

            { 3,     0,    -9.91381643467949629667E+02},
            { 3,     1,     3.38104557879559251887E+03},
            { 3,     2,    -4.11533190831385491038E+03},
            { 3,     3,     2.33454985382554878015E+03},
            { 3,     4,    -6.14157763172000272789E+02},
            { 3,     5,     5.77076386143955843977E+01},

            { 4,     0,     1.43797791534433781635E+03},
            { 4,     1,    -5.28916096082312014914E+03},
            { 4,     2,     6.91057399431488920527E+03},
            { 4,     3,    -4.26364756616708109505E+03},
            { 4,     4,     1.25783448603752435702E+03},
            { 4,     5,    -1.41377167956629875789E+02},

            { 5,     0,    -5.01761058636202108119E+02},
            { 5,     1,     1.99253292955415395227E+03},
            { 5,     2,    -2.72417925456980492527E+03},
            { 5,     3,     1.74907541356862247994E+03},
            { 5,     4,    -5.38448511479686771963E+02},
            { 5,     5,     6.36637946720314573668E+01}
        },
        {                                                                               // LMA (Z000100)
            { 0,     0,     7.69284893177311052568E+05},
            { 0,     1,    -4.60964236746346391737E+06},
            { 0,     2,     9.24477105650694482028E+06},
            { 0,     3,     2.96809810085771558806E+06},
            { 0,     4,    -5.47904367331531867385E+07},
            { 0,     5,     1.38449097585996419191E+08},
            { 0,     6,    -2.05590189063783049583E+08},
            { 0,     7,     2.11734651675154030323E+08},
            { 0,     8,    -1.60435957671923041344E+08},
            { 0,     9,     9.17609013885079324245E+07},
            { 0,    10,    -3.99674217866430878639E+07},
            { 0,    11,     1.32154540634083170444E+07},
            { 0,    12,    -3.26794649031273694709E+06},
            { 0,    13,     5.86210662802961654961E+05},
            { 0,    14,    -7.21158724758123280481E+04},
            { 0,    15,     5.44588834510656215571E+03},
            { 0,    16,    -1.90445008267364357835E+02},

            { 2,     0,     2.21188523097329214215E+07},
            { 2,     1,    -1.37731912987798571587E+08},
            { 2,     2,     3.57687151988543748856E+08},
            { 2,     3,    -4.63590327267392992973E+08},
            { 2,     4,     2.01845788407205879688E+08},
            { 2,     5,     2.87663499872765600681E+08},
            { 2,     6,    -5.59995404294795393944E+08},
            { 2,     7,     4.26777807943069577217E+08},
            { 2,     8,    -1.24264635254063099623E+08},
            { 2,     9,    -6.53086944807844161987E+07},
            { 2,    10,     9.24923524751316756010E+07},
            { 2,    11,    -5.24022445798896402121E+07},
            { 2,    12,     1.83838309755647853017E+07},
            { 2,    13,    -4.26215335990773420781E+06},
            { 2,    14,     6.40775684917965554632E+05},
            { 2,    15,    -5.69637563759163822397E+04},
            { 2,    16,     2.28325195741889183410E+03},

            { 4,     0,    -5.86271602527818381786E+07},
            { 4,     1,     4.78369991991747081280E+08},
            { 4,     2,    -1.74361930468868136406E+09},
            { 4,     3,     3.76621819705243778229E+09},
            { 4,     4,    -5.37751677961362648010E+09},
            { 4,     5,     5.35077163876746082306E+09},
            { 4,     6,    -3.80998322478477764130E+09},
            { 4,     7,     1.97099712467356157303E+09},
            { 4,     8,    -7.60832823835394382477E+08},
            { 4,     9,     2.42272232218582451344E+08},
            { 4,    10,    -8.05926269997112601995E+07},
            { 4,    11,     3.11124469145792499185E+07},
            { 4,    12,    -1.09693358577220626175E+07},
            { 4,    13,     2.84038850674086064100E+06},
            { 4,    14,    -4.84199551332909730263E+05},
            { 4,    15,     4.84970404654026060598E+04},
            { 4,    16,    -2.17011384837293599048E+03},

            { 6,     0,     9.83702780789473533630E+08},
            { 6,     1,    -6.23115263299602222443E+09},
            { 6,     2,     1.74950104067897109985E+10},
            { 6,     3,    -2.85244593670344696045E+10},
            { 6,     4,     2.95578184637604980469E+10},
            { 6,     5,    -1.97934299705980529785E+10},
            { 6,     6,     8.07257493996787166595E+09},
            { 6,     7,    -1.50554167974234604836E+09},
            { 6,     8,    -1.18252537929885953665E+08},
            { 6,     9,    -3.70554775613283040002E+06},
            { 6,    10,     1.23448486340704277158E+08},
            { 6,    11,    -7.41264279707325249910E+07},
            { 6,    12,     2.09071950705411061645E+07},
            { 6,    13,    -2.96454796610304014757E+06},
            { 6,    14,     1.18142086018433197751E+05},
            { 6,    15,     1.92799522383320254448E+04},
            { 6,    16,    -1.90271167151171880505E+03},

            { 8,     0,     6.69561843877016976476E+07},
            { 8,     1,     1.56256926101374864578E+09},
            { 8,     2,    -8.75745768682143020630E+09},
            { 8,     3,     1.94779956266501197815E+10},
            { 8,     4,    -2.33126852312486190796E+10},
            { 8,     5,     1.55144418171162052155E+10},
            { 8,     6,    -4.34538685387376117706E+09},
            { 8,     7,    -1.24875014824374175072E+09},
            { 8,     8,     1.36300258540026307106E+09},
            { 8,     9,    -2.47726977709512680769E+08},
            { 8,    10,    -1.37758820791601359844E+08},
            { 8,    11,     7.42609580452102422714E+07},
            { 8,    12,    -5.21190464590375591069E+06},
            { 8,    13,    -5.56763955730174109340E+06},
            { 8,    14,     1.98916589357259985991E+06},
            { 8,    15,    -2.82567601932666904759E+05},
            { 8,    16,     1.54458629248594224919E+04},

            {10,     0,    -2.91743516387244606018E+09},
            {10,     1,     1.43474200087861461639E+10},
            {10,     2,    -2.97121547997276992798E+10},
            {10,     3,     3.29244879152800636292E+10},
            {10,     4,    -1.97049122193543739319E+10},
            {10,     5,     4.44795698088213920593E+09},
            {10,     6,     1.68826391791642093658E+09},
            {10,     7,    -1.34669317233685207367E+09},
            {10,     8,     3.03101724866939246655E+08},
            {10,     9,    -4.99869293617860823870E+07},
            {10,    10,     1.42133397069574892521E+07},
            {10,    11,     2.69422055173848085105E+07},
            {10,    12,    -3.11643383659875914454E+07},
            {10,    13,     1.40837340470498204231E+07},
            {10,    14,    -3.33681038689361186698E+06},
            {10,    15,     4.13678001847658306360E+05},
            {10,    16,    -2.12956128078225046920E+04},

            {12,     0,     3.83066159035058796406E+08},
            {12,     1,    -3.57522953798943567276E+09},
            {12,     2,     1.00626903613678798676E+10},
            {12,     3,    -1.31025349462869167328E+10},
            {12,     4,     7.85384311503630828857E+09},
            {12,     5,    -2.41011892096044540405E+08},
            {12,     6,    -2.56408374643986368179E+09},
            {12,     7,     1.39171229309356403351E+09},
            {12,     8,    -1.44009585297834306955E+08},
            {12,     9,    -8.26904684745683073997E+07},
            {12,    10,    -2.38144140223466372117E+06},
            {12,    11,     1.79516107555324509740E+07},
            {12,    12,    -1.88674114853504090570E+06},
            {12,    13,    -2.42738823937244340777E+06},
            {12,    14,     9.99039133402326726355E+05},
            {12,    15,    -1.56469409353959286818E+05},
            {12,    16,     9.15292969974886182172E+03},

            {14,     0,     3.09020185410522365570E+09},
            {14,     1,    -1.36212178749835147858E+10},
            {14,     2,     2.66881548670487785339E+10},
            {14,     3,    -3.03946970897722930908E+10},
            {14,     4,     2.17784943226391410828E+10},
            {14,     5,    -9.40680629512597274780E+09},
            {14,     6,     1.39847537299111199379E+09},
            {14,     7,     1.10993825469740843773E+09},
            {14,     8,    -8.73110848697505235672E+08},
            {14,     9,     2.50292736143364071846E+08},
            {14,    10,     1.06462010817974358797E+07},
            {14,    11,    -3.22772786734968498349E+07},
            {14,    12,     1.01494283884754441679E+07},
            {14,    13,    -1.02655365919480670709E+06},
            {14,    14,    -1.37982127262780006276E+05},
            {14,    15,     4.30180039417772932211E+04},
            {14,    16,    -3.02473938756877396372E+03},

            {16,     0,    -1.56852265694742107391E+09},
            {16,     1,     7.15821216061697673798E+09},
            {16,     2,    -1.43193921312494201660E+10},
            {16,     3,     1.61786574433440303802E+10},
            {16,     4,    -1.08238279920512008667E+10},
            {16,     5,     3.66439506765549850464E+09},
            {16,     6,     3.04446730119516968727E+08},
            {16,     7,    -9.95706614031146168709E+08},
            {16,     8,     5.34070064296932458878E+08},
            {16,     9,    -1.68596081582052469254E+08},
            {16,    10,     4.50922229884207323194E+07},
            {16,    11,    -1.64262273564462922513E+07},
            {16,    12,     6.76394867549139633775E+06},
            {16,    13,    -2.05885415643938211724E+06},
            {16,    14,     3.96860804045895056333E+05},
            {16,    15,    -4.34586026756422943436E+04},
            {16,    16,     2.06736711765395784823E+03}
        },
        {                                                                               // HM (Z000100)
            { 0,     0,     2.66920628896922789863E+04},
            { 0,     1,     5.69598144516321923584E+06},
            { 0,     2,    -3.36106590539333373308E+07},
            { 0,     3,     9.25393487313766926527E+07},
            { 0,     4,    -1.64881303894646197557E+08},
            { 0,     5,     2.12786252918928653002E+08},
            { 0,     6,    -2.02235432385967493057E+08},
            { 0,     7,     1.36366495759085386992E+08},
            { 0,     8,    -5.89415461827230826020E+07},
            { 0,     9,     1.05569686857345979661E+07},
            { 0,    10,     4.26705694591859634966E+06},
            { 0,    11,    -3.31228471217951877043E+06},
            { 0,    12,     6.44460078807694022544E+05},
            { 0,    13,     1.79727468496846966445E+05},
            { 0,    14,    -1.14602002866395123419E+05},
            { 0,    15,     8.73345725996870532981E+03},
            { 0,    16,     1.03006575184149296547E+04},
            { 0,    17,    -4.45924381157494190120E+03},
            { 0,    18,     8.75712682690710266797E+02},
            { 0,    19,    -8.97432677867381869419E+01},
            { 0,    20,     3.90904107825454749658E+00},

            { 2,     0,    -4.66175954969746712595E+06},
            { 2,     1,     1.44553042693169433624E+07},
            { 2,     2,    -1.38532479348339065909E+07},
            { 2,     3,     1.61449323716102633625E+07},
            { 2,     4,    -5.04766978015377074480E+07},
            { 2,     5,     8.52236769195135086775E+07},
            { 2,     6,    -7.52658732323246151209E+07},
            { 2,     7,     3.77911912456169873476E+07},
            { 2,     8,    -1.05537492595233544707E+07},
            { 2,     9,     1.23416606663387943991E+06},
            { 2,    10,     1.07657456006634638470E+04},
            { 2,    11,     1.43573359232634276850E+05},
            { 2,    12,    -1.42859383483003708534E+05},
            { 2,    13,     5.63552937440058449283E+04},
            { 2,    14,    -1.09635084253284767328E+04},
            { 2,    15,     1.32853574780430949431E+03},
            { 2,    16,    -4.14076353943859373885E+02},
            { 2,    17,     1.47253804236653735416E+02},
            { 2,    18,    -2.47527050904476979554E+01},
            { 2,    19,     1.91210492011224597597E+00},
            { 2,    20,    -7.05892167374512213840E-02},

            { 4,     0,     1.02998189728514533490E+07},
            { 4,     1,    -3.48644858077521771193E+07},
            { 4,     2,     2.46925453657950051129E+07},
            { 4,     3,     2.65148321726851537824E+07},
            { 4,     4,    -3.81589329059887528419E+07},
            { 4,     5,     6.93028380426791380160E+05},
            { 4,     6,     1.99752902745822705328E+07},
            { 4,     7,    -1.14559117107986155897E+07},
            { 4,     8,     2.16656442928818566725E+06},
            { 4,     9,    -2.71810134198548854329E+05},
            { 4,    10,     2.68526744210913369898E+05},
            { 4,    11,    -9.77884134982192044845E+04},
            { 4,    12,    -1.53021212986694445135E+03},
            { 4,    13,     1.19037018096100418916E+04},
            { 4,    14,    -5.85794015733223022835E+03},
            { 4,    15,     1.56346892706033872855E+03},
            { 4,    16,    -2.75046129788926521087E+02},
            { 4,    17,     1.01298210263538095433E+02},
            { 4,    18,    -3.87720363616409713359E+01},
            { 4,    19,     6.66209009612612579332E+00},
            { 4,    20,    -3.93432972239732559050E-01},

            { 6,     0,    -1.23277652244743425399E+07},
            { 6,     1,     5.31219698912802562118E+07},
            { 6,     2,    -6.70881175458954200149E+07},
            { 6,     3,     1.38345048209509141743E+07},
            { 6,     4,     2.72927585462916977704E+07},
            { 6,     5,    -1.06984214016320873052E+07},
            { 6,     6,    -1.07550707413996290416E+07},
            { 6,     7,     9.12851204341692477465E+06},
            { 6,     8,    -2.45932699648965429515E+06},
            { 6,     9,     3.79086684429724293295E+05},
            { 6,    10,    -1.15384063428619745537E+05},
            { 6,    11,     9.37208910898514295695E+03},
            { 6,    12,     5.50001460585494169209E+03},
            { 6,    13,     8.99582878278114890236E+02},
            { 6,    14,     2.17976963642488072992E+02},
            { 6,    15,    -2.90083472695319642298E+02},
            { 6,    16,    -1.06912609256554702597E+02},
            { 6,    17,     7.44274626977634596869E+01},
            { 6,    18,    -1.32865712570078340349E+01},
            { 6,    19,     1.15176452822245067864E+00},
            { 6,    20,    -7.52331325998904648644E-02},

            { 8,     0,     6.81678287605022825301E+06},
            { 8,     1,    -3.98179300438016653061E+07},
            { 8,     2,     6.65612044507603123784E+07},
            { 8,     3,    -3.76181603687587901950E+07},
            { 8,     4,    -6.56919684995194338262E+06},
            { 8,     5,     1.43847081179449260235E+07},
            { 8,     6,    -4.06374374651811597869E+06},
            { 8,     7,     2.91022107704004331026E+05},
            { 8,     8,    -5.34165458237739163451E+05},
            { 8,     9,     2.58318759598450968042E+05},
            { 8,    10,    -2.04026808616285597964E+03},
            { 8,    11,    -1.91223168630281179503E+03},
            { 8,    12,    -3.69962055497544270111E+03},
            { 8,    13,    -9.23920259961619422029E+02},
            { 8,    14,     4.11256827202935312471E+02},
            { 8,    15,     1.26735718583368523582E+02},
            { 8,    16,    -4.68765072065054724249E+01},
            { 8,    17,     6.13281037162783970729E+00},
            { 8,    18,    -8.50877614045823293942E-01},
            { 8,    19,    -1.89283136262390905280E-01},
            { 8,    20,     5.75261353234685565705E-02},

            {10,     0,    -3.55154449666948523372E+05},
            {10,     1,     1.30929190946150850505E+07},
            {10,     2,    -3.15488824688208103180E+07},
            {10,     3,     2.70812807114166915417E+07},
            {10,     4,    -6.85649736623421218246E+06},
            {10,     5,    -2.40983373484842292964E+06},
            {10,     6,     1.08043478886758117005E+06},
            {10,     7,     1.11299945145412624697E+05},
            {10,     8,     4.58449855702385902987E+04},
            {10,     9,    -5.80789547527073082165E+04},
            {10,    10,    -4.28841803810808505659E+03},
            {10,    11,     4.99567543157751151739E+03},
            {10,    12,    -5.59329152251473260549E+02},
            {10,    13,     3.66269193209130548894E+02},
            {10,    14,    -4.12453739939884940213E+01},
            {10,    15,    -1.03224939821626904290E+01},
            {10,    16,    -1.26744942581002923987E+01},
            {10,    17,     2.50242050167653262704E+00},
            {10,    18,     4.89212292141251114952E-01},
            {10,    19,     7.41773856600715748508E-04},
            {10,    20,    -2.14966363670156049293E-02},

            {12,     0,    -1.71770314663222013041E+06},
            {12,     1,     1.02925354974316677544E+06},
            {12,     2,     4.75587077003400865942E+06},
            {12,     3,    -6.94997276824742369354E+06},
            {12,     4,     2.96480201554762385786E+06},
            {12,     5,     1.05305805605570232728E+05},
            {12,     6,    -2.38596544252316496568E+05},
            {12,     7,    -5.79872394392845744733E+04},
            {12,     8,     3.01928438393556607480E+04},
            {12,     9,     1.89217764968022629546E+03},
            {12,    10,    -4.88373167221754613365E+02},
            {12,    11,     7.01114489071818979937E+02},
            {12,    12,    -3.17674857603029806796E+02},
            {12,    13,    -9.41372505874503957557E+01},
            {12,    14,     1.83240842852757204184E+01},
            {12,    15,     1.08114635042327797976E+01},
            {12,    16,    -1.76741740095267596544E+00},
            {12,    17,     6.94852470861354931664E-01},
            {12,    18,    -2.89336670178484411942E-01},
            {12,    19,    -3.66873283517264768896E-03},
            {12,    20,     7.64308669905610499340E-03},

            {14,     0,     9.98775663821288966574E+05},
            {14,     1,    -2.18451837998499209061E+06},
            {14,     2,     1.40499697486835205927E+06},
            {14,     3,     4.46811799451349797891E+04},
            {14,     4,    -2.46225148374797863653E+05},
            {14,     5,    -7.05764615857048920589E+04},
            {14,     6,     9.16282262372330296785E+04},
            {14,     7,    -1.90152534930542642542E+04},
            {14,     8,     4.91854753979303222877E+03},
            {14,     9,    -2.98215458946964145071E+03},
            {14,    10,    -1.24575514575863394384E+02},
            {14,    11,     3.07045703026585499629E+02},
            {14,    12,     2.56707666795900628642E+01},
            {14,    13,    -2.55503462019516938142E+01},
            {14,    14,    -7.83210099348742638803E-01},
            {14,    15,     2.75956632611307561831E+00},
            {14,    16,    -6.89935042869725956294E-01},
            {14,    17,    -9.58434272104521506330E-02},
            {14,    18,     3.10995999077311284509E-02},
            {14,    19,     1.11788217615811057148E-02},
            {14,    20,    -2.50872008751619425190E-03},

            {16,     0,    -2.56910131149193679448E+05},
            {16,     1,     6.59413942707619862631E+05},
            {16,     2,    -5.99113926101273740642E+05},
            {16,     3,     1.79299729173374042148E+05},
            {16,     4,     3.57433247453356452752E+04},
            {16,     5,    -1.71963241067330600345E+04},
            {16,     6,    -8.07591195715742651373E+03},
            {16,     7,     1.64482614549100276236E+03},
            {16,     8,     1.58385752216435230366E+03},
            {16,     9,    -3.77870980619565330016E+02},
            {16,    10,    -1.99848538469132250839E+01},
            {16,    11,    -2.99325126666875540593E+01},
            {16,    12,     1.39774853358364321565E+01},
            {16,    13,     2.38884594919343573594E+00},
            {16,    14,    -1.77569226381652423008E+00},
            {16,    15,     2.51906789777950501641E-01},
            {16,    16,    -9.35456748727026554668E-02},
            {16,    17,     4.71639212726435025358E-02},
            {16,    18,    -2.76687708834074250132E-05},
            {16,    19,    -3.76375892191212375881E-03},
            {16,    20,     5.43491168280245645454E-04},

            {18,     0,     3.26573367579413752537E+04},
            {18,     1,    -8.28040252641195256729E+04},
            {18,     2,     6.32083466612035263097E+04},
            {18,     3,     1.10908655478984019283E+04},
            {18,     4,    -4.24696541649424107163E+04},
            {18,     5,     2.13154570304911503626E+04},
            {18,     6,    -4.01758333934846973534E+02},
            {18,     7,    -3.02668804705343745809E+03},
            {18,     8,     9.00317782135503534846E+02},
            {18,     9,    -1.71590190988977440156E+01},
            {18,    10,    -7.58773236769787295941E+00},
            {18,    11,    -1.47039559521070568593E+01},
            {18,    12,     5.86203446077681977755E+00},
            {18,    13,    -4.35902509094680423729E-01},
            {18,    14,    -3.96411294687243842549E-01},
            {18,    15,     1.90872936553090594147E-01},
            {18,    16,    -3.63349815421271857274E-02},
            {18,    17,     4.18942750214783450613E-03},
            {18,    18,    -1.95696242731476692869E-03},
            {18,    19,     6.55720172656096495292E-04},
            {18,    20,    -6.89531819695781068207E-05},

            {20,     0,    -1.79757648979984833204E+03},
            {20,     1,     4.21453422482446967479E+03},
            {20,     2,    -1.45713195934336681603E+03},
            {20,     3,    -5.27886800091239820176E+03},
            {20,     4,     8.07499233626505065331E+03},
            {20,     5,    -5.33458308189604667859E+03},
            {20,     6,     1.84109344300976272280E+03},
            {20,     7,    -2.89109084486973813455E+02},
            {20,     8,     5.54488971681655407053E+00},
            {20,     9,    -2.55037857514060917197E+00},
            {20,    10,     3.28120362746293148248E+00},
            {20,    11,    -1.43670757857304920435E+00},
            {20,    12,     8.25401835544457118665E-01},
            {20,    13,    -3.05993525251011089239E-01},
            {20,    14,     2.31463397238021693914E-02},
            {20,    15,     2.54799453700065604844E-02},
            {20,    16,    -1.16997904598278249649E-02},
            {20,    17,     2.34703472519215026668E-03},
            {20,    18,    -1.64174625118059316908E-04},
            {20,    19,    -1.87233906354765485458E-05},
            {20,    20,     3.01646763959664224869E-06}
        },
        {                                                                               // RECOM (Z000100)
            { 0,     0,     1.29041924400426033515E+01},
            { 0,     1,     1.46247169347020422592E+00},
            { 0,     2,    -2.56099379621667111451E+00},
            { 0,     3,     1.53157682928077298889E+00},
            { 0,     4,    -3.12651080854979501744E-01},

            { 1,     0,     1.13224784437094316836E+00},
            { 1,     1,    -1.21990615381324185584E+00},
            { 1,     2,     4.03103598240615390580E+00},
            { 1,     3,    -2.98697934479165105870E+00},
            { 1,     4,     7.11864900599255223668E-01},

            { 2,     0,    -1.84989977300083419109E+00},
            { 2,     1,     2.98072414316297784609E+00},
            { 2,     2,    -4.44057678245443998577E+00},
            { 2,     3,     2.62017135972970738322E+00},
            { 2,     4,    -6.33529954639474812694E-01},

            { 3,     0,     3.17264901652564645929E+00},
            { 3,     1,    -4.83562740400180235412E+00},
            { 3,     2,     4.06335321968123874825E+00},
            { 3,     3,    -1.38974858743467533095E+00},
            { 3,     4,     2.62122104574736836113E-01},

            { 4,     0,    -2.02663383998561608124E+00},
            { 4,     1,     3.18617136072649298484E+00},
            { 4,     2,    -2.29903125672929053991E+00},
            { 4,     3,     5.58176645345621169625E-01},
            { 4,     4,    -5.91094511981828871217E-02},

            { 5,     0,     3.82436170169052624956E-01},
            { 5,     1,    -6.15679347898177353748E-01},
            { 5,     2,     4.45802375081275292779E-01},
            { 5,     3,    -1.00232559556148334567E-01},
            { 5,     4,     6.38538257034518188376E-03}
        }
    },
    {                                                                                   // Metallicity Z001000 (0.01000)
        {                                                                               // LMR1 (Z001000)
            { 0,     0,     1.61732010479488863552E+01},
            { 0,     1,    -6.07511035663581200339E+00},
            { 0,     2,     8.58287622972541264232E+00},
            { 0,     3,    -3.10585034200754250833E+00},
            { 0,     4,    -7.06748010502642909358E-01},

            { 1,     0,    -1.90469269583824840630E+00},
            { 1,     1,     4.00007719345960026658E+01},
            { 1,     2,    -8.46518605722360462096E+01},
            { 1,     3,     5.91797386241528400319E+01},
            { 1,     4,    -1.13107588613914398223E+01},

            { 2,     0,    -8.70230071618597023075E+00},
            { 2,     1,    -2.89910577605233434895E+01},
            { 2,     2,     1.24597859970597440338E+02},
            { 2,     3,    -1.16207323210657051504E+02},
            { 2,     4,     3.01232264093183879083E+01},

            { 3,     0,     1.54664218136763143008E+01},
            { 3,     1,    -1.42381111546532466150E+01},
            { 3,     2,    -5.53996704254122320776E+01},
            { 3,     3,     7.85240667339935924929E+01},
            { 3,     4,    -2.52702301461299612129E+01},

            { 4,     0,    -4.23582571635512561414E+00},
            { 4,     1,     9.21228272981120710483E+00},
            { 4,     2,     6.44366319286792332832E+00},
            { 4,     3,    -1.82412432120631002874E+01},
            { 4,     4,     7.14935319397919499806E+00}
        },
        {                                                                               // LMR2 (Z001000)
            { 0,     0,     1.59314592243044277353E+01},
            { 0,     1,    -4.60401610759330548461E+00},
            { 0,     2,     6.24185962640975233739E+00},
            { 0,     3,    -4.97550125746145077699E+00},
            { 0,     4,     1.93048999642347540728E+00},
            { 0,     5,    -2.89477547622687292339E-01},

            { 1,     0,     1.01061304060731380616E+01},
            { 1,     1,    -3.46086484084549255158E+01},
            { 1,     2,     5.27667822465917879526E+01},
            { 1,     3,    -3.86837936682783691822E+01},
            { 1,     4,     1.36020700266615026663E+01},
            { 1,     5,    -1.85012309699758370485E+00},

            { 2,     0,     2.00677620782971253277E+01},
            { 2,     1,    -2.26769827650059667690E+01},
            { 2,     2,    -1.15679690058618760418E+01},
            { 2,     3,    -5.80771237470029455530E-01},
            { 2,     4,     1.37546460492662188102E+01},
            { 2,     5,    -4.53156102738701349608E+00},

            { 3,     0,    -5.09057297132129690453E+02},
            { 3,     1,     1.15813297915374073455E+03},
            { 3,     2,    -8.67528490581213873156E+02},
            { 3,     3,     3.72884696292590604116E+02},
            { 3,     4,    -1.38863311991984033966E+02},
            { 3,     5,     2.83694073141493667833E+01},

            { 4,     0,     2.61623056004032514466E+03},
            { 4,     1,    -6.70644555846752700745E+03},
            { 4,     2,     6.13584072176503832452E+03},
            { 4,     3,    -2.70631835697935139251E+03},
            { 4,     4,     6.36762028562752561811E+02},
            { 4,     5,    -7.14561574114333950547E+01},

            { 5,     0,    -2.63489815083014082120E+03},
            { 5,     1,     7.38504520674821560533E+03},
            { 5,     2,    -7.50607780122314488835E+03},
            { 5,     3,     3.63058087806650701168E+03},
            { 5,     4,    -8.69070310941866864596E+02},
            { 5,     5,     8.58816125753537704668E+01}
        },
        {                                                                               // LMA (Z001000)
            { 0,     0,    -2.09773228588053607382E+05},
            { 0,     1,     1.43389594778738403693E+06},
            { 0,     2,    -7.75133130475684534758E+06},
            { 0,     3,     3.63601069582714810967E+07},
            { 0,     4,    -1.18813444845486640930E+08},
            { 0,     5,     2.60283723461131542921E+08},
            { 0,     6,    -3.91364485037142515182E+08},
            { 0,     7,     4.10038222587301969528E+08},
            { 0,     8,    -2.94145223515916526318E+08},
            { 0,     9,     1.30304661949904114008E+08},
            { 0,    10,    -1.72744424556804075837E+07},
            { 0,    11,    -2.00820756134503111243E+07},
            { 0,    12,     1.57965446997256800532E+07},
            { 0,    13,    -5.16442762537559401244E+06},
            { 0,    14,     2.88797206094442051835E+05},
            { 0,    15,     4.85840419889588956721E+05},
            { 0,    16,    -2.37449358671291265637E+05},
            { 0,    17,     5.89009159964685095474E+04},
            { 0,    18,    -8.73970189435798602062E+03},
            { 0,    19,     7.40523972585050728412E+02},
            { 0,    20,    -2.77978725470881755655E+01},

            { 2,     0,     1.24556972785087689757E+08},
            { 2,     1,    -1.05485290802952396870E+09},
            { 2,     2,     4.02473985730900144577E+09},
            { 2,     3,    -9.11583314025122642517E+09},
            { 2,     4,     1.35512409280278778076E+10},
            { 2,     5,    -1.37244791677059879303E+10},
            { 2,     6,     9.40152573164973068237E+09},
            { 2,     7,    -4.02953830579527378082E+09},
            { 2,     8,     7.38075188234021663666E+08},
            { 2,     9,     2.04097261853698253632E+08},
            { 2,    10,    -1.31296867969154432416E+08},
            { 2,    11,    -3.72527087734079034999E+06},
            { 2,    12,     2.13105948168305270374E+07},
            { 2,    13,    -5.90374779784860368818E+06},
            { 2,    14,    -1.58750756547301629325E+05},
            { 2,    15,     1.73802705629212403437E+05},
            { 2,    16,     1.34389908452468196629E+05},
            { 2,    17,    -8.46886743101308384212E+04},
            { 2,    18,     2.01200486051142797805E+04},
            { 2,    19,    -2.33582333582660248794E+03},
            { 2,    20,     1.10848048996528689258E+02},

            { 4,     0,     5.05357352773523867130E+08},
            { 4,     1,    -3.31354820566867399216E+09},
            { 4,     2,     9.55477016233197402954E+09},
            { 4,     3,    -1.56750020973316497803E+10},
            { 4,     4,     1.55897307926384067535E+10},
            { 4,     5,    -8.83530815627405929565E+09},
            { 4,     6,     1.63784678158393263817E+09},
            { 4,     7,     1.26310745307742381096E+09},
            { 4,     8,    -9.06822750568983674049E+08},
            { 4,     9,     1.67932512482709646225E+08},
            { 4,    10,    -1.23911428774325661361E+07},
            { 4,    11,     5.71523280062625110149E+07},
            { 4,    12,    -4.26682323901077881455E+07},
            { 4,    13,     8.90929416129896789789E+06},
            { 4,    14,     1.67166402933667809702E+06},
            { 4,    15,    -4.69427729580903833266E+05},
            { 4,    16,    -5.00381548112256336026E+05},
            { 4,    17,     3.01400804692653589882E+05},
            { 4,    18,    -7.32120482670434430474E+04},
            { 4,    19,     8.80328726870056379994E+03},
            { 4,    20,    -4.34033934541470671320E+02},

            { 6,     0,    -1.07852724296717691422E+09},
            { 6,     1,     6.84582537129478645325E+09},
            { 6,     2,    -1.98636935069379844666E+10},
            { 6,     3,     3.49891808203243942261E+10},
            { 6,     4,    -4.16472301939004364014E+10},
            { 6,     5,     3.48742082858454284668E+10},
            { 6,     6,    -2.03515912933958663940E+10},
            { 6,     7,     7.63230475433247947693E+09},
            { 6,     8,    -1.32354384835169196129E+09},
            { 6,     9,    -1.69552074669079929590E+08},
            { 6,    10,     6.97440331351650208235E+07},
            { 6,    11,     4.60175066012795269489E+07},
            { 6,    12,    -2.68479745975134037435E+07},
            { 6,    13,     2.40770592574207345024E+06},
            { 6,    14,     2.26369602661154372618E+06},
            { 6,    15,    -1.54705866766018536873E+06},
            { 6,    16,     7.80276827808928559534E+05},
            { 6,    17,    -2.90525562987123674247E+05},
            { 6,    18,     6.73388625471522245789E+04},
            { 6,    19,    -8.55506609125114846393E+03},
            { 6,    20,     4.57366957264935876992E+02},

            { 8,     0,     6.01887171994113349915E+09},
            { 8,     1,    -3.59631408617360916138E+10},
            { 8,     2,     9.21373711482781219482E+10},
            { 8,     3,    -1.30859993208645553589E+11},
            { 8,     4,     1.09869343601253295898E+11},
            { 8,     5,    -5.19391019925100326538E+10},
            { 8,     6,     9.93795400740454673767E+09},
            { 8,     7,     1.17822075295209312439E+09},
            { 8,     8,     2.08547094828250437975E+08},
            { 8,     9,    -7.77891584686658740044E+08},
            { 8,    10,     9.79231656343366652727E+07},
            { 8,    11,     1.46029283089236587286E+08},
            { 8,    12,    -5.63426838521486073732E+07},
            { 8,    13,    -4.34994308839864679612E+04},
            { 8,    14,     2.67452041592313116416E+06},
            { 8,    15,     2.48601794890810997458E+05},
            { 8,    16,    -2.88661329611176857725E+05},
            { 8,    17,     3.44704079574566349038E+04},
            { 8,    18,     4.71065540904905537900E+03},
            { 8,    19,    -1.22791766059891961049E+03},
            { 8,    20,     6.52659303871330394031E+01},

            {10,     0,     4.09839446260637617111E+09},
            {10,     1,    -1.19243749172475032806E+10},
            {10,     2,     6.44260733488419055939E+08},
            {10,     3,     4.01784270278003311157E+10},
            {10,     4,    -7.02370927844391479492E+10},
            {10,     5,     5.53928774156074829102E+10},
            {10,     6,    -1.99839479763696517944E+10},
            {10,     7,     1.21178171327214643359E+08},
            {10,     8,     2.12195873823345303535E+09},
            {10,     9,    -3.04125284478327631950E+08},
            {10,    10,    -1.36583732859478831291E+08},
            {10,    11,     1.00033640716232210398E+07},
            {10,    12,     2.11985158940840736032E+07},
            {10,    13,    -1.41918619870038777590E+07},
            {10,    14,     1.06037309161468278617E+07},
            {10,    15,    -4.77449659335303772241E+06},
            {10,    16,     7.63050677152805146761E+05},
            {10,    17,     1.56993307696408010088E+05},
            {10,    18,    -8.52750815847208868945E+04},
            {10,    19,     1.33620566100046889915E+04},
            {10,    20,    -7.54749310106248799457E+02},

            {12,     0,    -2.77498112127752447128E+09},
            {12,     1,     1.99816442713947868347E+10},
            {12,     2,    -4.77683150712246398926E+10},
            {12,     3,     5.04117814886737747192E+10},
            {12,     4,    -1.67946036119296798706E+10},
            {12,     5,    -1.38693190414022674561E+10},
            {12,     6,     1.56277095143661766052E+10},
            {12,     7,    -4.57458119846564769745E+09},
            {12,     8,    -8.68468246746188163757E+08},
            {12,     9,     7.75119369125685334206E+08},
            {12,    10,    -1.16491136086922392249E+08},
            {12,    11,     1.47421104143672753125E+07},
            {12,    12,    -1.52535673836318179965E+07},
            {12,    13,    -1.70087987937960936688E+06},
            {12,    14,     4.67476437233580090106E+06},
            {12,    15,    -1.17463899375570146367E+06},
            {12,    16,    -1.10623359788225745433E+05},
            {12,    17,     8.96177599209738109494E+04},
            {12,    18,    -1.26272792247867328115E+04},
            {12,    19,     2.57183382652074556063E+02},
            {12,    20,     4.53946488958070730746E+01},

            {14,     0,    -1.19404109194643745422E+10},
            {14,     1,     4.75026388992112045288E+10},
            {14,     2,    -8.02570569421333770752E+10},
            {14,     3,     7.41607800924110260010E+10},
            {14,     4,    -4.05307739822724685669E+10},
            {14,     5,     1.41796423462739467621E+10},
            {14,     6,    -4.56807621318898487091E+09},
            {14,     7,     1.94941044600867557526E+09},
            {14,     8,    -5.25762401084450244904E+08},
            {14,     9,    -4.20587376462206020951E+07},
            {14,    10,     1.47395545866844896227E+07},
            {14,    11,     5.00008831927335187793E+07},
            {14,    12,    -3.55091173088149204850E+07},
            {14,    13,     1.31608539435279052705E+07},
            {14,    14,    -4.59113670743147470057E+06},
            {14,    15,     1.81956762125828582793E+06},
            {14,    16,    -5.94762161201574141160E+05},
            {14,    17,     1.30033906684856090578E+05},
            {14,    18,    -1.75126293747036870627E+04},
            {14,    19,     1.30486252871918713936E+03},
            {14,    20,    -4.01670777875890010478E+01}
        },
        {                                                                               // HM (Z001000)
            { 0,     0,     2.17405636868515895912E+05},
            { 0,     1,     2.06852476160948835313E+07},
            { 0,     2,    -6.14611410296131018549E+06},
            { 0,     3,    -3.71233464909802228212E+07},
            { 0,     4,     2.59646896356047354639E+07},
            { 0,     5,     2.70401430864754803479E+07},
            { 0,     6,    -5.03161212555342614651E+07},
            { 0,     7,     3.58598481040194258094E+07},
            { 0,     8,    -1.78006177931661754847E+07},
            { 0,     9,     9.30532460730251483619E+06},
            { 0,    10,    -5.15848465171259269118E+06},
            { 0,    11,     2.19027753214642358944E+06},
            { 0,    12,    -5.81644839415979455225E+05},
            { 0,    13,     7.99240321203274652362E+04},
            { 0,    14,    -2.51621643905049268142E+03},
            { 0,    15,    -3.81475220523696918917E+02},

            { 1,     0,    -3.06765020921669416130E+07},
            { 1,     1,    -7.19287141298301517963E+07},
            { 1,     2,    -4.77435954108010306954E+07},
            { 1,     3,     4.06090998470959842205E+08},
            { 1,     4,    -4.94869498500309586525E+08},
            { 1,     5,     2.69898549761012136936E+08},
            { 1,     6,    -7.97774313283056169748E+07},
            { 1,     7,     4.26528677017149850726E+07},
            { 1,     8,    -4.68106602932079061866E+07},
            { 1,     9,     3.02614077510188445449E+07},
            { 1,    10,    -1.02789476280537806451E+07},
            { 1,    11,     1.48880611471109278500E+06},
            { 1,    12,     6.37803643199263970018E+04},
            { 1,    13,    -1.54060450772956992296E+04},
            { 1,    14,    -1.20700325670922356949E+04},
            { 1,    15,     2.24656015193673965769E+03},

            { 2,     0,     1.75066912655571192503E+08},
            { 2,     1,     6.56024659875666722655E+07},
            { 2,     2,     1.23625106866149380803E+08},
            { 2,     3,    -8.66427535415050864220E+08},
            { 2,     4,     1.04647074302398228645E+09},
            { 2,     5,    -5.43304198956838369370E+08},
            { 2,     6,     6.62721930982138589025E+07},
            { 2,     7,     6.96526135511116981506E+07},
            { 2,     8,    -3.78209700669894739985E+07},
            { 2,     9,     5.84715094197653792799E+06},
            { 2,    10,     1.18073224759165989235E+06},
            { 2,    11,    -3.39494064604910614435E+05},
            { 2,    12,    -1.30009674988549915724E+05},
            { 2,    13,     2.04151331991039405693E+04},
            { 2,    14,     1.63426397272586054896E+04},
            { 2,    15,    -3.46623073630780800158E+03},

            { 3,     0,    -4.40378749002870142460E+08},
            { 3,     1,     1.90811958049112737179E+08},
            { 3,     2,    -3.46994414631571888924E+08},
            { 3,     3,     1.03839183062096667290E+09},
            { 3,     4,    -1.09191434318746423721E+09},
            { 3,     5,     6.03376176521021485329E+08},
            { 3,     6,    -1.86538522271839350462E+08},
            { 3,     7,     3.01234692302216663957E+07},
            { 3,     8,    -2.60876296034886920825E+06},
            { 3,     9,    -5.91722311695919837803E+05},
            { 3,    10,    -1.04200386924736900255E+05},
            { 3,    11,     4.46675275850050500594E+05},
            { 3,    12,    -9.51964972779181844089E+04},
            { 3,    13,     3.58457144797106002443E+03},
            { 3,    14,    -1.11371310452127327153E+04},
            { 3,    15,     2.85845040236972772618E+03},

            { 4,     0,     5.72344695093964815140E+08},
            { 4,     1,    -5.13427680762670874596E+08},
            { 4,     2,     7.89372994735144257545E+08},
            { 4,     3,    -1.00183048339205336571E+09},
            { 4,     4,     4.85858157195736885071E+08},
            { 4,     5,    -3.68810635419888198376E+07},
            { 4,     6,    -6.77671612100235223770E+07},
            { 4,     7,     3.51244092417839094996E+07},
            { 4,     8,    -4.04172301070910505950E+06},
            { 4,     9,     4.42372389482865459286E+05},
            { 4,    10,    -1.11637046404088055715E+06},
            { 4,    11,     1.32348000868789851665E+05},
            { 4,    12,     5.77673805733840563335E+04},
            { 4,    13,     1.39889883355121310160E+03},
            { 4,    14,     1.84537734558708234545E+03},
            { 4,    15,    -1.29546033045933745598E+03},

            { 5,     0,    -3.11964641996447086334E+08},
            { 5,     1,     3.20795370867759823799E+08},
            { 5,     2,    -9.73902886296924829483E+08},
            { 5,     3,     1.11109912862654781342E+09},
            { 5,     4,    -3.31644242757628619671E+08},
            { 5,     5,    -7.20401027258160263300E+07},
            { 5,     6,     5.49449315561505854130E+07},
            { 5,     7,    -1.50211363970660679042E+07},
            { 5,     8,     4.37729396379113662988E+05},
            { 5,     9,     1.78953864822479896247E+06},
            { 5,    10,    -9.64483419783796125557E+04},
            { 5,    11,    -1.03372025982109145843E+05},
            { 5,    12,     2.27654966719349886262E+03},
            { 5,    13,    -9.21556404868084791815E+03},
            { 5,    14,     2.29694125905770533791E+03},
            { 5,    15,     2.69172643454379056038E+02},

            { 6,     0,    -1.47363671236134380102E+08},
            { 6,     1,     3.98069493547490179539E+08},
            { 6,     2,     3.80600399916491806507E+08},
            { 6,     3,    -7.77899828440285444260E+08},
            { 6,     4,     2.13971953739619493484E+08},
            { 6,     5,     7.50766140657877773046E+07},
            { 6,     6,    -1.24408092643524613231E+07},
            { 6,     7,    -1.19486997694680951536E+07},
            { 6,     8,     1.72403185042079538107E+06},
            { 6,     9,     5.82624710496194384177E+04},
            { 6,    10,     6.62807449280618020566E+04},
            { 6,    11,    -5.73858576309664931614E+03},
            { 6,    12,     7.87238279678164053621E+03},
            { 6,    13,     5.58027583166682688898E+02},
            { 6,    14,    -7.80191280094344847384E+02},
            { 6,    15,    -4.41524171396966949033E+01},

            { 7,     0,     3.70057909980202019215E+08},
            { 7,     1,    -8.43582286640240669250E+08},
            { 7,     2,     3.06515664902774870396E+08},
            { 7,     3,     2.79680007567675352097E+08},
            { 7,     4,    -1.31316622141129702330E+08},
            { 7,     5,    -4.78340494412834197283E+07},
            { 7,     6,     2.05878316230771690607E+07},
            { 7,     7,     4.84660386411372851580E+06},
            { 7,     8,    -1.28697308684668410569E+06},
            { 7,     9,    -4.26216149967036908492E+05},
            { 7,    10,     1.35944377558333886554E+05},
            { 7,    11,    -2.78481544974478456425E+04},
            { 7,    12,     6.38116677538639805789E+02},
            { 7,    13,     1.86324844948431268676E+03},
            { 7,    14,    -3.54990721623104548144E+02},
            { 7,    15,     4.62881801006983621960E+01},

            { 8,     0,    -2.73364609220364511013E+08},
            { 8,     1,     6.70469859252789497375E+08},
            { 8,     2,    -4.83480623472343206406E+08},
            { 8,     3,     4.03385611767915189266E+07},
            { 8,     4,     7.35016723277738541365E+07},
            { 8,     5,    -1.62650907308186870068E+07},
            { 8,     6,    -2.72791540988075407222E+06},
            { 8,     7,     3.87529004919559811242E+04},
            { 8,     8,     2.21605940650244883727E+05},
            { 8,     9,    -4.99119267431606931495E+04},
            { 8,    10,     3.88533899694829306100E+04},
            { 8,    11,    -4.12454335026239914441E+03},
            { 8,    12,    -5.45806706724830746680E+02},
            { 8,    13,    -7.47011879871399173680E+02},
            { 8,    14,     2.58488648387669854856E+02},
            { 8,    15,    -2.59121452073608971034E+01},

            { 9,     0,     1.00383462413229718804E+08},
            { 9,     1,    -2.74271434515541374683E+08},
            { 9,     2,     2.71244504862839102745E+08},
            { 9,     3,    -1.11115631573666274548E+08},
            { 9,     4,     4.91877530747146159410E+06},
            { 9,     5,     1.25292149804861973971E+07},
            { 9,     6,    -5.59249804659485630691E+06},
            { 9,     7,     1.28330904890557797626E+06},
            { 9,     8,    -1.04317198392582082306E+05},
            { 9,     9,    -9.62038930517230619444E+03},
            { 9,    10,    -3.97809489634518195089E+03},
            { 9,    11,    -5.57959747341926686204E+02},
            { 9,    12,     8.89096577259679747840E+02},
            { 9,    13,    -5.60425945702165773099E+00},
            { 9,    14,    -4.78179067842277234490E+01},
            { 9,    15,     5.74915010119235514452E+00},

            {10,     0,    -1.55431711502586677670E+07},
            {10,     1,     4.79521485383154302835E+07},
            {10,     2,    -6.02665872499209195375E+07},
            {10,     3,     4.05979536634513586760E+07},
            {10,     4,    -1.60619130310957916081E+07},
            {10,     5,     3.70530953097213804722E+06},
            {10,     6,    -3.78917655865714070387E+05},
            {10,     7,    -1.84983833769938682963E+04},
            {10,     8,    -6.31914341540330042335E+03},
            {10,     9,     7.89696836505760438740E+03},
            {10,    10,    -2.17729950811668868482E+03},
            {10,    11,     7.29602008514560338881E+02},
            {10,    12,    -2.38483209776185589135E+02},
            {10,    13,     2.80674740177570534172E+01},
            {10,    14,     1.83297935389278832119E+00},
            {10,    15,    -4.35700622586662367208E-01}
        },
        {                                                                               // RECOM (Z001000)
            { 0,     0,     1.29690611505892015032E+01},
            { 0,     1,     1.15133146967729094179E+00},
            { 0,     2,    -2.03151799780861574973E+00},
            { 0,     3,     1.18911797645565364689E+00},
            { 0,     4,    -2.35730192355203121979E-01},

            { 1,     0,     6.77126934573034411358E-01},
            { 1,     1,     8.12545635631524709730E-01},
            { 1,     2,     1.02571355665047869721E+00},
            { 1,     3,    -1.20203984384456630252E+00},
            { 1,     4,     3.31810706642633568286E-01},

            { 2,     0,    -2.98789230039933073613E+00},
            { 2,     1,     4.40181013672451193486E+00},
            { 2,     2,    -3.18890331653727177041E+00},
            { 2,     3,     8.83079598994844006121E-01},
            { 2,     4,    -1.51558396747922291548E-01},

            { 3,     0,     6.22545012651770335310E+00},
            { 3,     1,    -1.26234063625864720848E+01},
            { 3,     2,     8.89677911090983108977E+00},
            { 3,     3,    -2.01593702143950714856E+00},
            { 3,     4,     1.24221816733169190816E-01},

            { 4,     0,    -3.48677279439661846894E+00},
            { 4,     1,     8.37556424986163960966E+00},
            { 4,     2,    -6.53850783620592235224E+00},
            { 4,     3,     1.69898655557246569536E+00},
            { 4,     4,    -1.21208024376684805890E-01},

            { 5,     0,     4.04481505949380593101E-01},
            { 5,     1,    -1.37770027286069152161E+00},
            { 5,     2,     1.26623104358849225548E+00},
            { 5,     3,    -3.70847513213560797674E-01},
            { 5,     4,     3.04915286309964846112E-02}
        }
    },
    {                                                                                   // Metallicity Z001500 (0.01500)
        {                                                                               // LMR1 (Z001500)
            { 0,     0,     1.60726718444739766767E+01},
            { 0,     1,    -6.24715135137397137299E+00},
            { 0,     2,     1.07672344037789748938E+01},
            { 0,     3,    -7.24052976085269506257E+00},
            { 0,     4,     1.46349058730338077439E+00},

            { 1,     0,    -3.63506353750313060402E+00},
            { 1,     1,     5.36864271459196729097E+01},
            { 1,     2,    -1.16942397088811162575E+02},
            { 1,     3,     9.34035699421548457622E+01},
            { 1,     4,    -2.50874003951161199666E+01},

            { 2,     0,    -7.91754787628987521941E+00},
            { 2,     1,    -5.94653927755674089894E+01},
            { 2,     2,     2.05203056353080341978E+02},
            { 2,     3,    -1.96558176123197796414E+02},
            { 2,     4,     5.97010203155139578257E+01},

            { 3,     0,     2.27278594676717915490E+01},
            { 3,     1,    -2.40820672531560520113E+00},
            { 3,     2,    -1.20196505857006002316E+02},
            { 3,     3,     1.49881801286317909216E+02},
            { 3,     4,    -5.13533879653931180087E+01},

            { 4,     0,    -1.08895531675941441563E+01},
            { 4,     1,     1.55183170137950732226E+01},
            { 4,     2,     1.95899488178189216114E+01},
            { 4,     3,    -3.88388622741310385322E+01},
            { 4,     4,     1.51805521241851977265E+01}
        },
        {                                                                               // LMR2 (Z001500)
            { 0,     0,     1.57661900292269407942E+01},
            { 0,     1,    -3.86739001681062477545E+00},
            { 0,     2,     4.90951789020641982120E+00},
            { 0,     3,    -3.83335877824608228792E+00},
            { 0,     4,     1.46458915608047957058E+00},
            { 0,     5,    -2.16702239649122857523E-01},

            { 1,     0,     8.37711942202392734202E+00},
            { 1,     1,    -2.70768859464224540545E+01},
            { 1,     2,     3.93888011339611523454E+01},
            { 1,     3,    -2.70047525691479570753E+01},
            { 1,     4,     8.70098442072441713435E+00},
            { 1,     5,    -1.06476013712186490245E+00},

            { 2,     0,     3.06826409038962175657E+01},
            { 2,     1,    -1.09405785435920606119E+02},
            { 2,     2,     1.81370175730107348500E+02},
            { 2,     3,    -1.81704582752789093547E+02},
            { 2,     4,     9.02892160295599097708E+01},
            { 2,     5,    -1.65129464015789579889E+01},

            { 3,     0,    -3.62814204845711003600E+02},
            { 3,     1,     1.08689918143087788849E+03},
            { 3,     2,    -1.31054018838510614842E+03},
            { 3,     3,     9.90790478329023699189E+02},
            { 3,     4,    -4.33419474916350452531E+02},
            { 3,     5,     7.62679478243772734913E+01},

            { 4,     0,     1.40418643511127538659E+03},
            { 4,     1,    -3.88233656649431759433E+03},
            { 4,     2,     3.84683355224100250780E+03},
            { 4,     3,    -2.03701250178350619535E+03},
            { 4,     4,     6.48142655530930028362E+02},
            { 4,     5,    -9.67319650941358304408E+01},

            { 5,     0,    -1.22407686126082307965E+03},
            { 5,     1,     3.76978222294538090864E+03},
            { 5,     2,    -3.95869928580390433126E+03},
            { 5,     3,     1.99118186826012993151E+03},
            { 5,     4,    -5.21844465216713615519E+02},
            { 5,     5,     6.04450098185751301116E+01}
        },
        {                                                                               // LMA (Z001500)
            { 0,     0,    -5.15111813953318560380E+04},
            { 0,     1,    -1.38581324034049781039E+06},
            { 0,     2,     1.39308014554255567491E+07},
            { 0,     3,    -5.71603108284442424774E+07},
            { 0,     4,     1.37981165851579815149E+08},
            { 0,     5,    -2.21159494412106066942E+08},
            { 0,     6,     2.47918615635263592005E+08},
            { 0,     7,    -1.97267821930349111557E+08},
            { 0,     8,     1.08701847878496184945E+08},
            { 0,     9,    -3.66657088282169476151E+07},
            { 0,    10,     2.47653496909504802898E+06},
            { 0,    11,     5.12285355178847070783E+06},
            { 0,    12,    -3.49757309776654699817E+06},
            { 0,    13,     1.39230419173041009344E+06},
            { 0,    14,    -4.29369586936748528387E+05},
            { 0,    15,     1.17250418664394310326E+05},
            { 0,    16,    -2.87188898192620908958E+04},
            { 0,    17,     5.69885720592272627982E+03},
            { 0,    18,    -8.01726305735812047715E+02},
            { 0,    19,     6.87272631810855187950E+01},
            { 0,    20,    -2.67041856670744914837E+00},

            { 2,     0,    -1.45339113482951931655E+07},
            { 2,     1,     1.12667823825035437942E+08},
            { 2,     2,    -4.17006587079701125622E+08},
            { 2,     3,     9.85341220899901628494E+08},
            { 2,     4,    -1.66272967677597117424E+09},
            { 2,     5,     2.10006772395531940460E+09},
            { 2,     6,    -2.00687866295555543900E+09},
            { 2,     7,     1.43029664485314631462E+09},
            { 2,     8,    -7.29760192023199796677E+08},
            { 2,     9,     2.43147726715557068586E+08},
            { 2,    10,    -4.05304287569378465414E+07},
            { 2,    11,    -4.01306759306349034887E+05},
            { 2,    12,    -7.32600274554128758609E+05},
            { 2,    13,     1.17855480548669607379E+06},
            { 2,    14,     2.50996201334879733622E+05},
            { 2,    15,    -6.29675976741036050953E+05},
            { 2,    16,     3.31786518519565521274E+05},
            { 2,    17,    -9.37220926397035946138E+04},
            { 2,    18,     1.56585215689283104439E+04},
            { 2,    19,    -1.46922968079042675527E+03},
            { 2,    20,     6.01559001927688683509E+01},

            { 4,     0,     2.33350108446424514055E+08},
            { 4,     1,    -1.60289968489695382118E+09},
            { 4,     2,     4.81990833858025169373E+09},
            { 4,     3,    -8.16360325217743110657E+09},
            { 4,     4,     8.14584602763201999664E+09},
            { 4,     5,    -4.07331865929028415680E+09},
            { 4,     6,    -5.06252506179879188538E+08},
            { 4,     7,     2.29763469151120567322E+09},
            { 4,     8,    -1.62128710050702285767E+09},
            { 4,     9,     5.32701252654437601566E+08},
            { 4,    10,    -4.48733723450238630176E+07},
            { 4,    11,    -2.11723976412126012146E+07},
            { 4,    12,     2.32885320360103854910E+06},
            { 4,    13,     2.03951553722685994580E+06},
            { 4,    14,    -3.13630035119212174322E+05},
            { 4,    15,    -7.05223482959375251085E+04},
            { 4,    16,    -5.37405402453973874799E+04},
            { 4,    17,     4.94362966298164537875E+04},
            { 4,    18,    -1.40518968918354457855E+04},
            { 4,    19,     1.83732556677333059270E+03},
            { 4,    20,    -9.50453907494923839749E+01},

            { 6,     0,    -1.76169219592734664679E+08},
            { 6,     1,     1.42698092226201272011E+09},
            { 6,     2,    -5.28067056308736801147E+09},
            { 6,     3,     1.15336708877598838806E+10},
            { 6,     4,    -1.61229222116397609711E+10},
            { 6,     5,     1.45899208559092102051E+10},
            { 6,     6,    -8.01445678044835853577E+09},
            { 6,     7,     1.90927120293207097054E+09},
            { 6,     8,     4.88043666148680567741E+08},
            { 6,     9,    -3.87662269468736469746E+08},
            { 6,    10,    -3.33959551647167652845E+07},
            { 6,    11,     9.24862097769977152348E+07},
            { 6,    12,    -2.25346824999888129532E+07},
            { 6,    13,    -5.86599007108959276229E+06},
            { 6,    14,     3.70074268507561367005E+06},
            { 6,    15,    -2.56417061140370758949E+05},
            { 6,    16,    -2.29390710704398312373E+05},
            { 6,    17,     6.94915698487437912263E+04},
            { 6,    18,    -7.11327988466724036698E+03},
            { 6,    19,     3.57413795845535204876E+01},
            { 6,    20,     2.80444597204764924925E+01},

            { 8,     0,     1.24245267843930959702E+09},
            { 8,     1,    -6.44322780962794589996E+09},
            { 8,     2,     1.38199616798107013702E+10},
            { 8,     3,    -1.50973626273141441345E+10},
            { 8,     4,     7.46748678501209640503E+09},
            { 8,     5,     8.19630542308683276176E+08},
            { 8,     6,    -3.01758840394945383072E+09},
            { 8,     7,     1.45032130975556683540E+09},
            { 8,     8,    -3.13138358723539471626E+08},
            { 8,     9,     1.84361824896876424551E+08},
            { 8,    10,    -1.69216581016204208136E+08},
            { 8,    11,     5.96086349759788066149E+07},
            { 8,    12,     3.41292303451553406194E+06},
            { 8,    13,    -8.38948524171856231987E+06},
            { 8,    14,     2.24497033832136029378E+06},
            { 8,    15,    -2.11810850483575311955E+05},
            { 8,    16,     5.84158368683356893598E+04},
            { 8,    17,    -4.11205502318392173038E+04},
            { 8,    18,     1.17356696779575449909E+04},
            { 8,    19,    -1.52952212651557533718E+03},
            { 8,    20,     7.85072147480235713601E+01},

            {10,     0,    -1.81932512063875246048E+09},
            {10,     1,     1.17609472824396915436E+10},
            {10,     2,    -3.26657704396900253296E+10},
            {10,     3,     4.99244642700828399658E+10},
            {10,     4,    -4.43242869808500061035E+10},
            {10,     5,     2.08982330227992744446E+10},
            {10,     6,    -2.02708169583341908455E+09},
            {10,     7,    -2.88544961080271482468E+09},
            {10,     8,     1.12173715222722339630E+09},
            {10,     9,     1.71997591831600487232E+08},
            {10,    10,    -1.91453234115329504013E+08},
            {10,    11,     3.25976767552792169154E+07},
            {10,    12,    -2.81013549671907734592E+05},
            {10,    13,     1.62781335264503862709E+06},
            {10,    14,    -1.02482288018334569642E+05},
            {10,    15,    -4.78252977657809446100E+05},
            {10,    16,     1.08856413307630777126E+05},
            {10,    17,     3.74121964994652880705E+04},
            {10,    18,    -2.00568317658710147953E+04},
            {10,    19,     3.30817698426615561402E+03},
            {10,    20,    -1.97037754182892768995E+02},

            {12,     0,     3.35585658860964119434E+08},
            {12,     1,     6.04561877929950594902E+08},
            {12,     2,    -5.13752404817218589783E+09},
            {12,     3,     1.06028901045264606476E+10},
            {12,     4,    -1.25273394138480796814E+10},
            {12,     5,     1.06453110355505924225E+10},
            {12,     6,    -6.95095817367120552063E+09},
            {12,     7,     3.16043713721437406540E+09},
            {12,     8,    -7.21320971392456412315E+08},
            {12,     9,    -5.26415434883427247405E+07},
            {12,    10,     3.88067712231941595674E+07},
            {12,    11,     2.20965666929983980954E+07},
            {12,    12,    -1.00345439357733502984E+07},
            {12,    13,    -2.83703418567746179178E+06},
            {12,    14,     2.26345678725755680352E+06},
            {12,    15,    -2.51391728097161831101E+05},
            {12,    16,    -1.16427576675173724652E+05},
            {12,    17,     3.66684894778185043833E+04},
            {12,    18,    -2.41298286705503824123E+03},
            {12,    19,    -3.40666737232470040908E+02},
            {12,    20,     4.25259705893395718590E+01},

            {14,     0,    -1.88270337231409168243E+09},
            {14,     1,     6.39532526253386974335E+09},
            {14,     2,    -7.82324355639520549774E+09},
            {14,     3,     2.03781762196393775940E+09},
            {14,     4,     5.26969578325209617615E+09},
            {14,     5,    -7.26948054175921058655E+09},
            {14,     6,     4.65417136459286022186E+09},
            {14,     7,    -1.63391954770838069916E+09},
            {14,     8,     1.87346108643729388714E+08},
            {14,     9,     8.89076657232203483582E+07},
            {14,    10,    -2.42219275225403495133E+07},
            {14,    11,    -1.37909139873882699758E+07},
            {14,    12,     6.82268681301416642964E+06},
            {14,    13,     7.02051292148833046667E+05},
            {14,    14,    -6.92968335637202719226E+05},
            {14,    15,    -3.38601851769478351343E+05},
            {14,    16,     3.62546876010872249026E+05},
            {14,    17,    -1.30131463530392415123E+05},
            {14,    18,     2.47381054441249470983E+04},
            {14,    19,    -2.51595863192902470473E+03},
            {14,    20,     1.08731772056314042629E+02}
        },
        {                                                                               // HM (Z001500)
            { 0,     0,    -6.40222930769867496565E+05},
            { 0,     1,     1.10169653380964145064E+08},
            { 0,     2,    -1.83288798818755358458E+08},
            { 0,     3,     2.34685038301938652992E+08},
            { 0,     4,    -2.90868321715029954910E+08},
            { 0,     5,     2.58803902657175153494E+08},
            { 0,     6,    -1.63288883819993436337E+08},
            { 0,     7,     8.19842075092330276966E+07},
            { 0,     8,    -3.74597286079473644495E+07},
            { 0,     9,     1.73063130301320850849E+07},
            { 0,    10,    -7.45417360680353827775E+06},
            { 0,    11,     2.34362116585996653885E+06},
            { 0,    12,    -3.93361441283391148318E+05},
            { 0,    13,    -1.10781031619090754248E+02},
            { 0,    14,     1.09339428886468467681E+04},
            { 0,    15,    -1.22318915545462755290E+03},

            { 1,     0,    -2.80689302632132731378E+06},
            { 1,     1,    -8.98291909237491011620E+08},
            { 1,     2,     1.18592188686839365959E+09},
            { 1,     3,    -7.65219587451042652130E+08},
            { 1,     4,     5.70295193113373160362E+08},
            { 1,     5,    -3.92352877095889389515E+08},
            { 1,     6,     1.43709948555291593075E+08},
            { 1,     7,     3.04550315389612223953E+06},
            { 1,     8,    -3.80446096340993121266E+07},
            { 1,     9,     2.31604841598522365093E+07},
            { 1,    10,    -5.65692487392537854612E+06},
            { 1,    11,    -3.33505242560792976292E+03},
            { 1,    12,     1.18223884409625112312E+05},
            { 1,    13,     8.29620922556901932694E+04},
            { 1,    14,    -3.61397359909771475941E+04},
            { 1,    15,     3.94064567463254252289E+03},

            { 2,     0,     2.01327877490426152945E+08},
            { 2,     1,     3.10118016671073770523E+09},
            { 2,     2,    -4.22208827383246850967E+09},
            { 2,     3,     1.99760734988433146477E+09},
            { 2,     4,    -7.30736768670006632805E+08},
            { 2,     5,     4.63059271373849034309E+08},
            { 2,     6,    -2.61651098515476316214E+08},
            { 2,     7,     1.29333847253416240215E+08},
            { 2,     8,    -5.09144396244669929147E+07},
            { 2,     9,     9.08995892725536786020E+06},
            { 2,    10,    -4.72198901409551501274E+05},
            { 2,    11,     8.76993455861743539572E+05},
            { 2,    12,    -4.06936925562300020829E+05},
            { 2,    13,    -6.55970599629438947886E+03},
            { 2,    14,     3.07809196263461344643E+04},
            { 2,    15,    -4.44326151571276750474E+03},

            { 3,     0,    -1.20734569145530676842E+09},
            { 3,     1,    -5.34256358175416564941E+09},
            { 3,     2,     8.35253308946118927002E+09},
            { 3,     3,    -3.55132392736585330963E+09},
            { 3,     4,     3.64000813302139878273E+08},
            { 3,     5,     6.22636854320012405515E+07},
            { 3,     6,    -7.19534735248802900314E+07},
            { 3,     7,     3.05025745163594335318E+07},
            { 3,     8,     7.68708163315938971937E+06},
            { 3,     9,    -3.07800794374816305935E+06},
            { 3,    10,    -2.14286733640326838940E+06},
            { 3,    11,     9.30921378395954845473E+05},
            { 3,    12,    -5.63446497619864894659E+04},
            { 3,    13,     1.63174566022897852235E+03},
            { 3,    14,    -1.16542300609449248441E+04},
            { 3,    15,     2.36843518043800349915E+03},

            { 4,     0,     3.21429576530229854584E+09},
            { 4,     1,     4.39998454027949714661E+09},
            { 4,     2,    -9.63131293606875228882E+09},
            { 4,     3,     4.22524392135507345200E+09},
            { 4,     4,    -1.44018430614179670811E+08},
            { 4,     5,    -1.74848804728251099586E+08},
            { 4,     6,     3.49971447999607846141E+07},
            { 4,     7,    -2.03914002418685518205E+07},
            { 4,     8,    -1.16793722923461743630E+06},
            { 4,     9,     4.07216598363398248330E+06},
            { 4,    10,    -6.63453501301569864154E+05},
            { 4,    11,    -2.07226327457722654799E+05},
            { 4,    12,     7.70487153054169175448E+04},
            { 4,    13,     3.74511172163199034912E+03},
            { 4,    14,    -2.35924085496946099738E+03},
            { 4,    15,    -2.00074150360817185401E+02},

            { 5,     0,    -4.86051768030780982971E+09},
            { 5,     1,    -2.52283411110093832016E+08},
            { 5,     2,     6.15744921080406188965E+09},
            { 5,     3,    -2.94612427147580671310E+09},
            { 5,     4,    -1.64832248763824313879E+08},
            { 5,     5,     2.56579223709281235933E+08},
            { 5,     6,    -1.66799195149700418115E+07},
            { 5,     7,     1.29958687433429018711E+05},
            { 5,     8,    -1.39862347051455453038E+06},
            { 5,     9,     6.31508472462175530382E+05},
            { 5,    10,     1.11229199938330348232E+04},
            { 5,    11,    -1.09408809112410061061E+05},
            { 5,    12,    -3.31703747407727632890E+03},
            { 5,    13,    -7.42683722320846641196E+02},
            { 5,    14,     3.67269071514491452035E+03},
            { 5,    15,    -4.51573620652630495442E+02},

            { 6,     0,     4.54118409271736431122E+09},
            { 6,     1,    -2.90767159976795244217E+09},
            { 6,     2,    -1.58025838498863267899E+09},
            { 6,     3,     1.16398375723541212082E+09},
            { 6,     4,     1.09676718899973109365E+08},
            { 6,     5,    -1.17450609424565017223E+08},
            { 6,     6,     3.66726234801336098462E+06},
            { 6,     7,     1.73320626005022414029E+06},
            { 6,     8,    -2.36895981676343455911E+06},
            { 6,     9,     6.50776589120257762261E+05},
            { 6,    10,     9.80770387879816262284E+04},
            { 6,    11,     6.08450624306344252545E+03},
            { 6,    12,     1.44103835104128138482E+04},
            { 6,    13,    -8.31901750284624540654E+03},
            { 6,    14,    -3.29760727015796192063E+02},
            { 6,    15,     2.22378172851269994226E+02},

            { 7,     0,    -2.66951759357458686829E+09},
            { 7,     1,     2.78244247697885799408E+09},
            { 7,     2,    -5.23326483537968039513E+08},
            { 7,     3,    -1.32406331732717037201E+08},
            { 7,     4,    -7.73866632845929916948E+06},
            { 7,     5,    -1.31289335217684842646E+07},
            { 7,     6,     8.68883711662754230201E+06},
            { 7,     7,     6.03206328956170473248E+06},
            { 7,     8,    -1.68968455437217419967E+06},
            { 7,     9,    -1.77911507059704919811E+05},
            { 7,    10,    -2.13170445228394455626E+04},
            { 7,    11,    -1.98704865429771484742E+03},
            { 7,    12,     1.19148875731268594791E+03},
            { 7,    13,     2.79295327351346804790E+03},
            { 7,    14,    -2.78996535146859685028E+02},
            { 7,    15,    -4.49988682514610260910E+01},

            { 8,     0,     9.50171459722543239594E+08},
            { 8,     1,    -1.19109509072384023666E+09},
            { 8,     2,     4.45397440506503641605E+08},
            { 8,     3,    -4.46003071982323154807E+07},
            { 8,     4,    -3.80750898013454861939E+06},
            { 8,     5,     1.37090059550021272153E+07},
            { 8,     6,    -5.08256188045835960656E+06},
            { 8,     7,    -1.65700011583761894144E+06},
            { 8,     8,     4.63444510281368449796E+05},
            { 8,     9,     6.44070571055476248148E+04},
            { 8,    10,     5.85348885756574745756E+04},
            { 8,    11,    -2.07048623592222938896E+04},
            { 8,    12,    -1.63741310696323444063E+03},
            { 8,    13,     6.19764619828629179210E+02},
            { 8,    14,    -5.07010754652976132206E+01},
            { 8,    15,     1.04668911802231079378E+01},

            { 9,     0,    -1.81045244850725114346E+08},
            { 9,     1,     2.23702214978076577187E+08},
            { 9,     2,    -5.39576641499614790082E+07},
            { 9,     3,    -3.48971727291814684868E+07},
            { 9,     4,     2.05198001371338181198E+07},
            { 9,     5,    -1.21495342617135704495E+06},
            { 9,     6,    -3.16323564470732118934E+06},
            { 9,     7,     1.45230267412367230281E+06},
            { 9,     8,     4.31289414334815301117E+04},
            { 9,     9,    -1.48747640901944541838E+05},
            { 9,    10,     8.63897175520341988886E+03},
            { 9,    11,     7.41770223687096313370E+03},
            { 9,    12,    -3.19056811887575008768E+02},
            { 9,    13,    -3.21467405033000488856E+02},
            { 9,    14,     5.29898436774840320140E+01},
            { 9,    15,    -3.45842797011066638291E+00},

            {10,     0,     1.29428707331854440272E+07},
            {10,     1,    -8.40693462242027930915E+06},
            {10,     2,    -1.66950000229333285242E+07},
            {10,     3,     2.57591525679362528026E+07},
            {10,     4,    -1.54869040150702800602E+07},
            {10,     5,     4.65623680094749014825E+06},
            {10,     6,    -3.66297825326654070523E+05},
            {10,     7,    -1.74355945057241042377E+05},
            {10,     8,     1.17860134753315287526E+04},
            {10,     9,     1.88480662367682598415E+04},
            {10,    10,    -1.91529993314948478655E+03},
            {10,    11,    -1.26547072451270787496E+03},
            {10,    12,     2.41500958409452209708E+02},
            {10,    13,     1.39658146862240180042E+01},
            {10,    14,    -6.16658750062571581196E+00},
            {10,    15,     4.30482965248756621612E-01}
        },
        {                                                                               // RECOM (Z001500)
            { 0,     0,     1.30114084816216948326E+01},
            { 0,     1,     9.82927123720802486950E-01},
            { 0,     2,    -1.78831205264574788494E+00},
            { 0,     3,     1.04428283556912560037E+00},
            { 0,     4,    -2.05873948734389328186E-01},

            { 1,     0,    -3.88724524733189344405E-02},
            { 1,     1,     3.79408475484787910403E+00},
            { 1,     2,    -2.91118101352734681697E+00},
            { 1,     3,     7.77250458270365651714E-01},
            { 1,     4,    -5.93043121201380050989E-03},

            { 2,     0,    -1.19760400814100687050E+00},
            { 2,     1,    -3.70516957950696035340E+00},
            { 2,     2,     8.13081389223339456862E+00},
            { 2,     3,    -4.84039188270623732535E+00},
            { 2,     4,     8.01677869943414500575E-01},

            { 3,     0,     5.25028472760498932104E+00},
            { 3,     1,    -6.49412032630474111983E+00},
            { 3,     2,    -1.37460399330354388070E+00},
            { 3,     3,     3.58589212553425529251E+00},
            { 3,     4,    -8.38094449459952040016E-01},

            { 4,     0,    -3.35142128249576209953E+00},
            { 4,     1,     6.88602644973469768530E+00},
            { 4,     2,    -2.99992638566578184722E+00},
            { 4,     3,    -4.88499106441473263107E-01},
            { 4,     4,     2.76990054008695862908E-01},

            { 5,     0,     3.32597771543929221494E-01},
            { 5,     1,    -1.23233281888843992924E+00},
            { 5,     2,     8.33374957715134034864E-01},
            { 5,     3,    -6.82770912855755612858E-02},
            { 5,     4,    -2.81426477843974842674E-02}
        }
    },
    {                                                                                   // Metallicity Z002000 (0.02000)
        {                                                                               // LMR1 (Z002000)
            { 0,     0,     1.50305864307172658556E+01},
            { 0,     1,     4.96110465703095582235E-01},
            { 0,     2,    -9.16608282566077403608E-01},
            { 0,     3,     2.48033316991016244968E-01},

            { 1,     0,     1.77238406649739155263E+00},
            { 1,     1,    -6.52669943496246296455E-01},
            { 1,     2,     6.23616025599144307989E-01},
            { 1,     3,    -1.77703681033998805994E-01}
        },
        {                                                                               // LMR2 (Z002000)
            { 0,     0,     1.56617552673034907684E+01},
            { 0,     1,    -3.38271409511212528543E+00},
            { 0,     2,     3.99712000854973803499E+00},
            { 0,     3,    -3.02506056136319223526E+00},
            { 0,     4,     1.12651646333780397491E+00},
            { 0,     5,    -1.62914377442334834534E-01},

            { 1,     0,     7.36910246144982838956E+00},
            { 1,     1,    -2.30880577686855161801E+01},
            { 1,     2,     3.32294193400783512971E+01},
            { 1,     3,    -2.22327077943441828722E+01},
            { 1,     4,     6.85986578867425578210E+00},
            { 1,     5,    -7.85596865167521363205E-01},

            { 2,     0,     3.04288311387925105578E+01},
            { 2,     1,    -1.38697534163042433875E+02},
            { 2,     2,     2.72764619820308155340E+02},
            { 2,     3,    -2.82250602598139892052E+02},
            { 2,     4,     1.37077309352478550863E+02},
            { 2,     5,    -2.43363801583426564434E+01},

            { 3,     0,    -2.40391883810966504598E+02},
            { 3,     1,     7.20123937632566253342E+02},
            { 3,     2,    -9.56047390568349555906E+02},
            { 3,     3,     8.80126339460319968566E+02},
            { 3,     4,    -4.43060386994394036719E+02},
            { 3,     5,     8.34534841118213392974E+01},

            { 4,     0,     1.19040215556355087756E+03},
            { 4,     1,    -3.06320152050599153881E+03},
            { 4,     2,     2.78453473219604893529E+03},
            { 4,     3,    -1.43582241909097660937E+03},
            { 4,     4,     4.99720580375658244066E+02},
            { 4,     5,    -8.44179485778146130315E+01},

            { 5,     0,    -1.27279708339606577283E+03},
            { 5,     1,     3.68563309227265335721E+03},
            { 5,     2,    -3.66200350857788362191E+03},
            { 5,     3,     1.74828529177925724980E+03},
            { 5,     4,    -4.43446115665101842751E+02},
            { 5,     5,     5.16038082522078980219E+01}
        },
        {                                                                               // LMA (Z002000)
            { 0,     0,    -1.51049254914130983707E+03},
            { 0,     1,     5.90804592923560267081E+03},
            { 0,     2,    -9.78411552163013584504E+03},
            { 0,     3,     9.03281646038073631644E+03},
            { 0,     4,    -5.04733341717106668511E+03},
            { 0,     5,     1.70840439977420237483E+03},
            { 0,     6,    -3.16251664756395882705E+02},
            { 0,     7,     1.74202707579704139107E+01},
            { 0,     8,     3.67713293499369164863E+00},
            { 0,     9,    -5.09055045574688058707E-01},

            { 1,     0,     1.18642296961242536781E+05},
            { 1,     1,    -4.22882345998080738354E+05},
            { 1,     2,     6.02728568485423107632E+05},
            { 1,     3,    -4.09618748060796875507E+05},
            { 1,     4,     8.90979215331482701004E+04},
            { 1,     5,     5.95414095173196401447E+04},
            { 1,     6,    -5.18590927506359803374E+04},
            { 1,     7,     1.72331732723822569824E+04},
            { 1,     8,    -2.81277261133411593619E+03},
            { 1,     9,     1.86487606226860776815E+02},

            { 2,     0,    -1.40854420295521779917E+06},
            { 2,     1,     3.16920022432578727603E+06},
            { 2,     2,     4.61300052596119523514E+05},
            { 2,     3,    -8.61105917229603044689E+06},
            { 2,     4,     1.23556398056155946106E+07},
            { 2,     5,    -8.94036121833491511643E+06},
            { 2,     6,     3.81112897666227677837E+06},
            { 2,     7,    -9.71158692551792715676E+05},
            { 2,     8,     1.37401910653425060445E+05},
            { 2,     9,    -8.32876128728674302693E+03},

            { 3,     0,     1.72458181532576158643E+07},
            { 3,     1,    -4.52772539129179418087E+07},
            { 3,     2,     2.01148809934781342745E+07},
            { 3,     3,     6.09227038679622411728E+07},
            { 3,     4,    -1.07679655401528432965E+08},
            { 3,     5,     8.29801717072715312243E+07},
            { 3,     6,    -3.64493890278426855803E+07},
            { 3,     7,     9.44902866139182634652E+06},
            { 3,     8,    -1.35153547941456316039E+06},
            { 3,     9,     8.25342644314703647979E+04},

            { 4,     0,    -1.36692876826326847076E+08},
            { 4,     1,     4.51025336598177015781E+08},
            { 4,     2,    -5.31011502233843207359E+08},
            { 4,     3,     1.64806847028889596462E+08},
            { 4,     4,     2.15493288928440541029E+08},
            { 4,     5,    -2.68977365806154668331E+08},
            { 4,     6,     1.38776598747677832842E+08},
            { 4,     7,    -3.89913542043956518173E+07},
            { 4,     8,     5.84939274544744472951E+06},
            { 4,     9,    -3.68464314099946233910E+05},

            { 5,     0,     5.88007632627574801445E+08},
            { 5,     1,    -2.19546849066643857956E+09},
            { 5,     2,     3.29578755621284484863E+09},
            { 5,     3,    -2.44640117334454774857E+09},
            { 5,     4,     7.63562415504432678223E+08},
            { 5,     5,     1.37457055430822163820E+08},
            { 5,     6,    -2.06529556624937295914E+08},
            { 5,     7,     7.50013034409245699644E+07},
            { 5,     8,    -1.26735707259846013039E+07},
            { 5,     9,     8.54653687692407402210E+05},

            { 6,     0,    -1.38327254361235642433E+09},
            { 6,     1,     5.60614067713224506378E+09},
            { 6,     2,    -9.43342679709010314941E+09},
            { 6,     3,     8.54102479910387897491E+09},
            { 6,     4,    -4.42476607837071800232E+09},
            { 6,     5,     1.22015302616584563255E+09},
            { 6,     6,    -9.17636537403069436550E+07},
            { 6,     7,    -4.19402333332933112979E+07},
            { 6,     8,     1.20023414749641995877E+07},
            { 6,     9,    -9.83676207032694481313E+05},

            { 7,     0,     1.70433842457238912582E+09},
            { 7,     1,    -7.46564615829306602478E+09},
            { 7,     2,     1.36213117744959812164E+10},
            { 7,     3,    -1.36449069424619369507E+10},
            { 7,     4,     8.23382455328104782104E+09},
            { 7,     5,    -3.05178841976632213593E+09},
            { 7,     6,     6.64362036529520750046E+08},
            { 7,     7,    -7.09507497840798497200E+07},
            { 7,     8,     7.89097510117984027602E+05},
            { 7,     9,     3.47019468185705947690E+05},

            { 8,     0,    -9.03133463246678471565E+08},
            { 8,     1,     4.47621684993862724304E+09},
            { 8,     2,    -8.96349174100218391418E+09},
            { 8,     3,     9.77159178772765159607E+09},
            { 8,     4,    -6.46222307786082077026E+09},
            { 8,     5,     2.69448710911580371857E+09},
            { 8,     6,    -7.04153944276203036308E+08},
            { 8,     7,     1.09216364710150659084E+08},
            { 8,     8,    -8.71251345965249091387E+06},
            { 8,     9,     2.35094345904298679670E+05},

            { 9,     0,     7.47826861721657812595E+07},
            { 9,     1,    -6.76538905537033677101E+08},
            { 9,     2,     1.73090915363686656952E+09},
            { 9,     3,    -2.18473390758292007446E+09},
            { 9,     4,     1.61084727759012627602E+09},
            { 9,     5,    -7.39637661531803250313E+08},
            { 9,     6,     2.13872151698962450027E+08},
            { 9,     7,    -3.76329373731744587421E+07},
            { 9,     8,     3.63438629714747751132E+06},
            { 9,     9,    -1.44132561512217595009E+05}
        },
        {                                                                               // HM (Z002000)
            { 0,    -4,     2.53729226092276048660E+09},
            { 0,    -3,    -9.90007540851090812683E+09},
            { 0,    -2,     1.52593834113249797821E+10},
            { 0,    -1,    -1.01549528384376564026E+10},
            { 0,     0,     2.99829578010063916445E+07},
            { 0,     1,     5.78893490124735832214E+09},
            { 0,     2,    -6.14609898490778446198E+09},
            { 0,     3,     4.04278241615104389191E+09},
            { 0,     4,    -1.80654362590145325661E+09},
            { 0,     5,     4.91334348953023314476E+08},
            { 0,     6,    -5.44061721887778788805E+07},
            { 0,     7,    -6.47383359300961066037E+06},
            { 0,     8,     1.76558650155847985297E+06},
            { 0,     9,     1.00096546886151510989E+05},
            { 0,    10,    -3.39961914804177431506E+04},

            { 1,    -4,    -1.43765508369811935425E+10},
            { 1,    -3,     5.53386381611402206421E+10},
            { 1,    -2,    -8.64972441405090179443E+10},
            { 1,    -1,     6.18263781388091888428E+10},
            { 1,     0,    -1.27031002818969459534E+10},
            { 1,     1,    -1.22741674009294738770E+10},
            { 1,     2,     1.32545543793032627106E+10},
            { 1,     3,    -7.72103742705006504059E+09},
            { 1,     4,     2.99797494889310359955E+09},
            { 1,     5,    -4.51316994491591095924E+08},
            { 1,     6,    -2.12359218253416419029E+08},
            { 1,     7,     1.28106564253871873021E+08},
            { 1,     8,    -2.52927841362378038466E+07},
            { 1,     9,     1.52832625190069107339E+06},
            { 1,    10,     5.33991439763792805024E+04},

            { 2,    -4,     3.53888765741276016235E+10},
            { 2,    -3,    -1.30949004957460281372E+11},
            { 2,    -2,     2.03740123625964263916E+11},
            { 2,    -1,    -1.45485633356229827881E+11},
            { 2,     0,     3.38169397075351562500E+10},
            { 2,     1,     1.75117599864639778137E+10},
            { 2,     2,    -1.68526525631832427979E+10},
            { 2,     3,     7.96671769498775196075E+09},
            { 2,     4,    -3.11421864218404531479E+09},
            { 2,     5,     6.72407024087293863297E+08},
            { 2,     6,     1.69479333713295638561E+08},
            { 2,     7,    -1.48711963495587915182E+08},
            { 2,     8,     3.25634279867792949080E+07},
            { 2,     9,    -1.75731226257729995996E+06},
            { 2,    10,    -1.35079742770175478654E+05},

            { 3,    -4,    -5.02259172277233886719E+10},
            { 3,    -3,     1.70340081645104858398E+11},
            { 3,    -2,    -2.59178892024369842529E+11},
            { 3,    -1,     1.80456668127881774902E+11},
            { 3,     0,    -4.10109717020726623535E+10},
            { 3,     1,    -1.59603847360727577209E+10},
            { 3,     2,     1.16582282107490043640E+10},
            { 3,     3,    -2.69128104888696050644E+09},
            { 3,     4,     7.08683177906285405159E+08},
            { 3,     5,    -4.34911479088057935238E+08},
            { 3,     6,     1.01044743218226939440E+08},
            { 3,     7,     1.50179933915116582066E+07},
            { 3,     8,    -4.86460542015785723925E+06},
            { 3,     9,    -1.18190563420071476139E+06},
            { 3,    10,     2.78623554417385661509E+05},

            { 4,    -4,     4.64705908425810546875E+10},
            { 4,    -3,    -1.30875109087104644775E+11},
            { 4,    -2,     1.84384613263450225830E+11},
            { 4,    -1,    -1.14315528459969711304E+11},
            { 4,     0,     1.54155640604285240173E+10},
            { 4,     1,     1.32641401056199913025E+10},
            { 4,     2,    -5.51530197341053485870E+09},
            { 4,     3,    -1.70800161550766736269E+08},
            { 4,     4,     7.78830644995285868645E+08},
            { 4,     5,    -2.06842292132320821285E+08},
            { 4,     6,    -2.58832894926396012306E+07},
            { 4,     7,     2.22092620407946482301E+07},
            { 4,     8,    -6.95661957240117993206E+06},
            { 4,     9,     1.92690144883910729550E+06},
            { 4,    10,    -2.33277229786988784326E+05},

            { 5,    -4,    -3.03660798334242973328E+10},
            { 5,    -3,     5.92434377201234664917E+10},
            { 5,    -2,    -6.54677950676404571533E+10},
            { 5,    -1,     1.98006473040911102295E+10},
            { 5,     0,     1.80741102193409042358E+10},
            { 5,     1,    -1.24110750642199649811E+10},
            { 5,     2,     2.06936088729867076874E+09},
            { 5,     3,    -5.55200935775747060776E+08},
            { 5,     4,     3.31793995932478666306E+08},
            { 5,     5,    -3.77282162409072965384E+07},
            { 5,     6,     3.77271593664872169029E+05},
            { 5,     7,    -2.39515931142189726233E+06},
            { 5,     8,     4.78976734628147096373E+05},
            { 5,     9,    -2.24492564827198133571E+05},
            { 5,    10,     4.79792411042873136466E+04},

            { 6,    -4,     1.39020977034095058441E+10},
            { 6,    -3,    -1.24214230056445674896E+10},
            { 6,    -2,     3.40159755029327201843E+09},
            { 6,    -1,     1.62527040079761524200E+10},
            { 6,     0,    -2.40572617307258148193E+10},
            { 6,     1,     1.01576765639793376923E+10},
            { 6,     2,    -2.79329356560366511345E+08},
            { 6,     3,    -3.70225353787005186081E+08},
            { 6,     4,    -1.05583829363182082772E+08},
            { 6,     5,     4.22791698590517193079E+07},
            { 6,     6,    -2.37510881991530815139E+06},
            { 6,     7,    -2.48551971484777750447E+06},
            { 6,     8,     2.30980219498265953735E+06},
            { 6,     9,    -5.97702223240277264267E+05},
            { 6,    10,     4.54666582402248095605E+04},

            { 7,    -4,    -2.83252506215792465210E+09},
            { 7,    -3,    -5.87107593619644260406E+09},
            { 7,    -2,     9.35640963237767410278E+09},
            { 7,    -1,    -9.10900080282122421265E+09},
            { 7,     0,     9.58587735224136924744E+09},
            { 7,     1,    -5.49030619645441913605E+09},
            { 7,     2,     1.36275859178348875046E+09},
            { 7,     3,    -8.57005601597409099340E+07},
            { 7,     4,    -7.86904883546570241451E+07},
            { 7,     5,     4.99692720418420732021E+07},
            { 7,     6,    -7.09569033777875266969E+06},
            { 7,     7,    -8.55338083130266284570E+05},
            { 7,     8,    -5.61264059584917849861E+05},
            { 7,     9,     3.08073506627650989685E+05},
            { 7,    10,    -3.28288523737060459098E+04},

            { 8,    -4,    -1.46525129680925726891E+09},
            { 8,    -3,     9.31328032341473007202E+09},
            { 8,    -2,    -1.05576564874556827545E+10},
            { 8,    -1,     4.15304574768818807602E+09},
            { 8,     0,    -5.00505353118861734867E+08},
            { 8,     1,     2.38979281271329969168E+08},
            { 8,     2,    -3.43075090044980287552E+08},
            { 8,     3,     1.53399431073772042990E+08},
            { 8,     4,     6.15683430398282781243E+06},
            { 8,     5,    -2.10229696774827726185E+07},
            { 8,     6,     6.68991617223582346924E+05},
            { 8,     7,     2.39900342987180408090E+06},
            { 8,     8,    -4.43480899599777883850E+05},
            { 8,     9,    -2.42983373542564295349E+04},
            { 8,    10,     7.76394675149368777056E+03},

            { 9,    -4,     1.21311579654507946968E+09},
            { 9,    -3,    -5.24070703990287113190E+09},
            { 9,    -2,     7.06700234139572620392E+09},
            { 9,    -1,    -4.37710070889741706848E+09},
            { 9,     0,     1.25442565132006978989E+09},
            { 9,     1,    -4.98846050957186818123E+07},
            { 9,     2,    -6.82854203030566722155E+07},
            { 9,     3,     3.25620574902522787452E+07},
            { 9,     4,    -2.20218804932029694319E+07},
            { 9,     5,     9.02162431089685671031E+06},
            { 9,     6,    -2.40492887758824275807E+05},
            { 9,     7,    -8.68585806973894243129E+05},
            { 9,     8,     2.25402435254328796873E+05},
            { 9,     9,    -1.55503066541077005240E+04},
            { 9,    10,    -3.94672900362825259890E+02},

            {10,    -4,    -2.52924874342063546181E+08},
            {10,    -3,     1.08433769070899581909E+09},
            {10,    -2,    -1.75274841405356764793E+09},
            {10,    -1,     1.53857146715462470055E+09},
            {10,     0,    -8.55304471302122354507E+08},
            {10,     1,     3.31495823209429860115E+08},
            {10,     2,    -9.56754638305779546499E+07},
            {10,     3,     1.92680795062012895942E+07},
            {10,     4,    -1.01626479920090991072E+06},
            {10,     5,    -6.69255818740246933885E+05},
            {10,     6,     1.56196811545844484499E+04},
            {10,     7,     9.97457630392971332185E+04},
            {10,     8,    -2.92981386761484645831E+04},
            {10,     9,     2.94914409981776043423E+03},
            {10,    10,    -6.54138557526380992613E+01}
        },
        {                                                                               // RECOM (Z002000)
            { 0,     0,     1.34142729153601880654E+01},
            { 0,     1,    -8.33977084929894418863E-01},
            { 0,     2,     6.31803234438807592710E-01},
            { 0,     3,    -1.70217751115103232973E-01},

            { 1,     0,     1.57297096125805579980E+00},
            { 1,     1,     3.04879789211891349954E-01},
            { 1,     2,    -6.10228805816398045536E-01},
            { 1,     3,     2.31090540934666771600E-01},

            { 2,     0,    -1.41108730414326188907E+00},
            { 2,     1,     1.23004856793249794933E+00},
            { 2,     2,    -2.83350980192163703908E-01},
            { 2,     3,    -4.69140139332027000796E-02},

            { 3,     0,     5.62788156340448986192E-01},
            { 3,     1,    -6.89182559641456582433E-01},
            { 3,     2,     2.69502539954297459790E-01},
            { 3,     3,    -2.14046891159626988255E-02}
        }
    },
    {                                                                                   // Metallicity Z003000 (0.03000)
        {                                                                               // LMR1 (Z003000)
            { 0,     0,     1.54866491842634630416E+01},
            { 0,     1,    -2.35824863480211543987E+00},
            { 0,     2,     1.34711604769709447638E+00},
            { 0,     3,     2.79306201502264705994E+00},
            { 0,     4,    -2.56152597553130023655E+00},

            { 1,     0,     7.89806670922946318925E-01},
            { 1,     1,     1.42719958374201905116E+01},
            { 1,     2,    -2.47713064947769225910E+01},
            { 1,     3,     8.04278071310951681028E+00},
            { 1,     4,     3.66647228234742073028E+00},

            { 2,     0,     1.66175550737136745738E+01},
            { 2,     1,    -7.35725916220352331720E+01},
            { 2,     2,     1.08031105163848295092E+02},
            { 2,     3,    -5.83176668396733290933E+01},
            { 2,     4,     6.99282806748704732769E+00},

            { 3,     0,    -4.02973895611078916090E+01},
            { 3,     1,     1.40473419123728177738E+02},
            { 3,     2,    -1.80465312924060071964E+02},
            { 3,     3,     9.78861437573207098239E+01},
            { 3,     4,    -1.77262651464499150222E+01},

            { 4,     0,     2.12553651132979481986E+01},
            { 4,     1,    -7.07553640430039791909E+01},
            { 4,     2,     8.71360130418007798880E+01},
            { 4,     3,    -4.65547494341429839437E+01},
            { 4,     4,     9.00395775876992843223E+00}
        },
        {                                                                               // LMR2 (Z003000)
            { 0,     0,     1.56146519862202897144E+01},
            { 0,     1,    -3.25295555117365031705E+00},
            { 0,     2,     3.80218538552061113833E+00},
            { 0,     3,    -2.85992446872646643996E+00},
            { 0,     4,     1.05223371503530960247E+00},
            { 0,     5,    -1.49446662816108466476E-01},

            { 1,     0,     7.77113365437422576321E+00},
            { 1,     1,    -2.69242148162415411150E+01},
            { 1,     2,     4.21377981796297760297E+01},
            { 1,     3,    -3.05791586287080576767E+01},
            { 1,     4,     1.03058811476673746199E+01},
            { 1,     5,    -1.30770786874792355192E+00},

            { 2,     0,     4.37708329067558921111E+00},
            { 2,     1,     4.05698553767332636966E+00},
            { 2,     2,    -4.39975579696197893753E+00},
            { 2,     3,    -4.02226205989940837071E+01},
            { 2,     4,     3.97982102495747582793E+01},
            { 2,     5,    -9.68255076300243366916E+00},

            { 3,     0,    -1.14171738557163791938E+02},
            { 3,     1,     1.44654983520949940612E+02},
            { 3,     2,     1.25658014849041876460E+02},
            { 3,     3,    -7.61233422914715873731E+01},
            { 3,     4,    -4.97277512863255921616E+01},
            { 3,     5,     2.27566957849081354937E+01},

            { 4,     0,     8.52284180915415390700E+02},
            { 4,     1,    -2.11793768130346688849E+03},
            { 4,     2,     1.43778303702111816165E+03},
            { 4,     3,    -3.53485239324663382376E+02},
            { 4,     4,     5.96099116904753358881E+01},
            { 4,     5,    -1.52841095719764492600E+01},

            { 5,     0,    -7.27744971361332773085E+02},
            { 5,     1,     2.37447873535868984618E+03},
            { 5,     2,    -2.29824328800873217915E+03},
            { 5,     3,     9.63017826362198434254E+02},
            { 5,     4,    -1.94941240897338843752E+02},
            { 5,     5,     1.82013102751265911650E+01}
        },
        {                                                                               // LMA (Z003000)
            { 0,    -5,    -2.57403227906218804419E+07},
            { 0,    -4,     1.91919278581367194653E+08},
            { 0,    -3,    -6.42722772365061283112E+08},
            { 0,    -2,     1.26557774018711972237E+09},
            { 0,    -1,    -1.59807275300796556473E+09},
            { 0,     0,     1.29431773149255251884E+09},
            { 0,     1,    -5.73752193441836237907E+08},
            { 0,     2,    -3.44428973121863901615E+07},
            { 0,     3,     2.58171971437299102545E+08},
            { 0,     4,    -2.08805258131750792265E+08},
            { 0,     5,     9.92256309368441253901E+07},
            { 0,     6,    -3.16719607180158011615E+07},
            { 0,     7,     6.90780604186993371695E+06},
            { 0,     8,    -9.94125170976186404005E+05},
            { 0,     9,     8.54546806862898811232E+04},
            { 0,    10,    -3.33231922729193092891E+03},

            { 2,    -5,     1.34762385503116726875E+08},
            { 2,    -4,    -1.74399176370349168777E+09},
            { 2,    -3,     8.78154209435013771057E+09},
            { 2,    -2,    -2.48007449760757026672E+10},
            { 2,    -1,     4.53496317537334213257E+10},
            { 2,     0,    -5.77108539579762878418E+10},
            { 2,     1,     5.31677144188315505981E+10},
            { 2,     2,    -3.61960228471166305542E+10},
            { 2,     3,     1.83377540802361297607E+10},
            { 2,     4,    -6.87766287682024002075E+09},
            { 2,     5,     1.87236466119424223900E+09},
            { 2,     6,    -3.54625421447190940380E+08},
            { 2,     7,     4.26550618172946497798E+07},
            { 2,     8,    -2.50770450767832808197E+06},
            { 2,     9,    -2.72280137764262872224E+04},
            { 2,    10,     8.76822961747765839391E+03},

            { 4,    -5,     3.64113230917176294327E+09},
            { 4,    -4,    -2.19993132639651412964E+10},
            { 4,    -3,     5.51823742774377899170E+10},
            { 4,    -2,    -6.74001478202328414917E+10},
            { 4,    -1,     1.97503194124485588074E+10},
            { 4,     0,     6.55027963168653793335E+10},
            { 4,     1,    -1.20047754730772918701E+11},
            { 4,     2,     1.10973078984599426270E+11},
            { 4,     3,    -6.64446966999905242920E+10},
            { 4,     4,     2.71468286742218437195E+10},
            { 4,     5,    -7.48815090425585556030E+09},
            { 4,     6,     1.28782471484118795395E+09},
            { 4,     7,    -9.90502311621036529541E+07},
            { 4,     8,    -7.30136450831011310220E+06},
            { 4,     9,     2.22677327622929215431E+06},
            { 4,    10,    -1.42107178068599314429E+05},

            { 6,    -5,    -1.69094733263322086334E+10},
            { 6,    -4,     1.18528670819833084106E+11},
            { 6,    -3,    -3.71969875084951599121E+11},
            { 6,    -2,     6.88453453818612182617E+11},
            { 6,    -1,    -8.33162353980749511719E+11},
            { 6,     0,     6.92006751423001220703E+11},
            { 6,     1,    -4.06654644455712707520E+11},
            { 6,     2,     1.78505518636014801025E+11},
            { 6,     3,    -6.90610161810206146240E+10},
            { 6,     4,     3.03147898201140594482E+10},
            { 6,     5,    -1.43313685446296138763E+10},
            { 6,     6,     5.54512398337787055969E+09},
            { 6,     7,    -1.50349138969701838493E+09},
            { 6,     8,     2.63695637427561759949E+08},
            { 6,     9,    -2.68951201610311269760E+07},
            { 6,    10,     1.21421668571477057412E+06},

            { 8,    -5,     1.62520989450171852112E+10},
            { 8,    -4,    -1.22380441325664489746E+11},
            { 8,    -3,     4.03718140234380981445E+11},
            { 8,    -2,    -7.62315076723765625000E+11},
            { 8,    -1,     8.94726561199124023438E+11},
            { 8,     0,    -6.46755386562959716797E+11},
            { 8,     1,     2.37624479000744293213E+11},
            { 8,     2,     2.53444076942211761475E+10},
            { 8,     3,    -7.69330323678240509033E+10},
            { 8,     4,     3.81364515969743728638E+10},
            { 8,     5,    -6.76366483839152717590E+09},
            { 8,     6,    -1.57825522737156033516E+09},
            { 8,     7,     1.18057807852275967598E+09},
            { 8,     8,    -2.92124731642549633980E+08},
            { 8,     9,     3.56642162956294342875E+07},
            { 8,    10,    -1.79988696399625157937E+06},

            {10,    -5,    -5.29522123969188880920E+09},
            {10,    -4,     4.37924632110345230103E+10},
            {10,    -3,    -1.66864842737637054443E+11},
            {10,    -2,     3.77193811645369873047E+11},
            {10,    -1,    -5.45755946513809143066E+11},
            {10,     0,     5.12213537621273925781E+11},
            {10,     1,    -2.95784493345500854492E+11},
            {10,     2,     7.77097973981596527100E+10},
            {10,     3,     2.15638580927770462036E+10},
            {10,     4,    -2.78122413775628700256E+10},
            {10,     5,     1.08470967326972980499E+10},
            {10,     6,    -1.76167227877269959450E+09},
            {10,     7,    -1.19039692308903947473E+08},
            {10,     8,     1.03663300992385908961E+08},
            {10,     9,    -1.72756292958132103086E+07},
            {10,    10,     1.02701810218507179525E+06},

            {12,    -5,     1.80343253078517799377E+10},
            {12,    -4,    -1.26691173309104354858E+11},
            {12,    -3,     3.56351695791545471191E+11},
            {12,    -2,    -5.53758351468571411133E+11},
            {12,    -1,     5.36138039982130004883E+11},
            {12,     0,    -3.39866845315202819824E+11},
            {12,     1,     1.44163375900454193115E+11},
            {12,     2,    -4.49762794375830764771E+10},
            {12,     3,     1.71580814489644603729E+10},
            {12,     4,    -1.12449749792033348083E+10},
            {12,     5,     6.78490694496251773834E+09},
            {12,     6,    -2.84463375643244552612E+09},
            {12,     7,     7.96841892772628188133E+08},
            {12,     8,    -1.43851289638218492270E+08},
            {12,     9,     1.52044649291066508740E+07},
            {12,    10,    -7.16406403363849385642E+05},

            {14,    -5,     1.16520079605483184814E+11},
            {14,    -4,    -5.39465539126750244141E+11},
            {14,    -3,     1.09764790470851672363E+12},
            {14,    -2,    -1.26552677553650268555E+12},
            {14,    -1,     8.70273333877262084961E+11},
            {14,     0,    -3.15431150532989440918E+11},
            {14,     1,     9.53251283628360748291E+08},
            {14,     2,     5.57105743806749877930E+10},
            {14,     3,    -2.28806497057574310303E+10},
            {14,     4,     3.86158310274266898632E+08},
            {14,     5,     3.09330777606074094772E+09},
            {14,     6,    -1.37702519450405025482E+09},
            {14,     7,     3.09963461585896790028E+08},
            {14,     8,    -3.97062345451950877905E+07},
            {14,     9,     2.66685145412934198976E+06},
            {14,    10,    -6.71761975842174288118E+04}
        },
        {                                                                               // HM (Z003000)
            { 0,     0,     3.12154248604355193675E+06},
            { 0,     1,     1.40693554095653563738E+08},
            { 0,     2,    -3.71595276305998325348E+08},
            { 0,     3,     1.76670755437760323286E+08},
            { 0,     4,     2.02951692746982872486E+08},
            { 0,     5,    -2.92336468998508751392E+08},
            { 0,     6,     1.68110673419797122478E+08},
            { 0,     7,    -5.15604220053785145283E+07},
            { 0,     8,     4.97172830721258837730E+06},
            { 0,     9,     1.57284540400977246463E+06},
            { 0,    10,    -1.75988885345202215831E+04},
            { 0,    11,    -2.97943821613070613239E+05},
            { 0,    12,     8.61773562932996719610E+04},
            { 0,    13,    -9.04684915798255497066E+03},
            { 0,    14,     5.27682281201942714688E+02},
            { 0,    15,    -4.83809114868389968933E+01},

            { 1,     0,    -7.35394878249169588089E+07},
            { 1,     1,    -5.76316277886640667915E+08},
            { 1,     2,     2.41011072754778337479E+09},
            { 1,     3,    -2.36860622483675336838E+09},
            { 1,     4,     7.09952326698126316071E+08},
            { 1,     5,     2.10905643436775207520E+08},
            { 1,     6,    -2.58844355418170630932E+08},
            { 1,     7,     1.12699839262904345989E+08},
            { 1,     8,    -2.48447791184037849307E+07},
            { 1,     9,     3.05833750104070967063E+06},
            { 1,    10,    -2.27715413729118229821E+06},
            { 1,    11,     1.28163629808912402950E+06},
            { 1,    12,    -2.45961475246419315226E+05},
            { 1,    13,     7.71686543784326659079E+03},
            { 1,    14,     9.24986993266994090845E+02},
            { 1,    15,     1.19680958374041750858E+02},

            { 2,     0,     3.01712491925559639931E+08},
            { 2,     1,     1.35997955868716746569E+08},
            { 2,     2,    -4.33322825757982063293E+09},
            { 2,     3,     5.50388031191586112976E+09},
            { 2,     4,    -2.64093562924157762527E+09},
            { 2,     5,     6.12809232748553872108E+08},
            { 2,     6,    -9.04877006798645257950E+07},
            { 2,     7,     9.55706229663263447583E+06},
            { 2,     8,    -1.24463357289487235248E+07},
            { 2,     9,     1.14754493813732806593E+07},
            { 2,    10,    -2.62512806448984425515E+06},
            { 2,    11,    -1.38746075201159896096E+05},
            { 2,    12,    -3.06434607562738347042E+04},
            { 2,    13,     5.51760989160078097484E+04},
            { 2,    14,    -7.38420266268990508252E+03},
            { 2,    15,    -7.88663649547120257921E+01},

            { 3,     0,    -2.52314621471064984798E+08},
            { 3,     1,     1.59695673223666906357E+09},
            { 3,     2,     2.56902829277027750015E+09},
            { 3,     3,    -4.80529198174804210663E+09},
            { 3,     4,     2.24889814265354108810E+09},
            { 3,     5,    -3.74476785226586341858E+08},
            { 3,     6,     5.43274903410231973976E+06},
            { 3,     7,     4.29081287132954373956E+07},
            { 3,     8,    -2.75819523820120133460E+07},
            { 3,     9,     5.40615162239633128047E+06},
            { 3,    10,    -1.34637905180512159131E+06},
            { 3,    11,     8.50683737928409944288E+05},
            { 3,    12,    -1.24430090013735956745E+05},
            { 3,    13,    -2.52241309959800564684E+04},
            { 3,    14,     3.73821107091213934837E+03},
            { 3,    15,     3.32248890148269254041E+02},

            { 4,     0,    -5.60517200779632091522E+08},
            { 4,     1,    -1.89984697471257662773E+09},
            { 4,     2,     6.79333540273263901472E+07},
            { 4,     3,     1.90352540671050238609E+09},
            { 4,     4,    -7.50747156377030849457E+08},
            { 4,     5,    -4.72567503656732067466E+07},
            { 4,     6,     3.00467742178344316781E+07},
            { 4,     7,    -1.93436903588892333210E+06},
            { 4,     8,     2.95066336620633816347E+06},
            { 4,     9,     1.70540967107301065698E+06},
            { 4,    10,    -1.16727812393627688289E+06},
            { 4,    11,     4.34517006904071604367E+04},
            { 4,    12,     2.35401513264981185785E+04},
            { 4,    13,     6.11857844834931529476E+03},
            { 4,    14,     1.03171863933184098983E+02},
            { 4,    15,    -3.75092728367472147966E+02},

            { 5,     0,     1.24819667430232143402E+09},
            { 5,     1,     2.19967497561418384314E+08},
            { 5,     2,    -7.89014919695381075144E+07},
            { 5,     3,    -5.97878301337307453156E+08},
            { 5,     4,     1.59211666115665256977E+08},
            { 5,     5,     5.16454561360483840108E+07},
            { 5,     6,    -4.18833435321484040469E+06},
            { 5,     7,    -3.39808681376820290461E+06},
            { 5,     8,    -7.22078069451165269129E+05},
            { 5,     9,     1.99512798423410247779E+05},
            { 5,    10,     9.07387635003017785493E+04},
            { 5,    11,    -1.52441052433743388974E+04},
            { 5,    12,    -8.25416524310134082043E+03},
            { 5,    13,     2.18147893849961155865E+03},
            { 5,    14,    -6.64708537198572457783E+02},
            { 5,    15,     1.44256443204359840138E+02},

            { 6,     0,    -8.01485488430039882660E+08},
            { 6,     1,     5.38297769089608192444E+08},
            { 6,     2,    -6.06655185404267787933E+08},
            { 6,     3,     4.95703300265223979950E+08},
            { 6,     4,    -1.09494887044486209750E+08},
            { 6,     5,     1.03584531923860888928E+07},
            { 6,     6,    -5.11998278185982257128E+06},
            { 6,     7,    -3.64656091820526169613E+06},
            { 6,     8,     1.75311863341458095238E+06},
            { 6,     9,    -4.11223165191869309638E+05},
            { 6,    10,     1.88387970340288331499E+05},
            { 6,    11,    -3.47362162359701178502E+04},
            { 6,    12,     8.57684754506128047069E+03},
            { 6,    13,    -2.74638773621986092621E+03},
            { 6,    14,     9.61819122641777823901E+01},
            { 6,    15,     2.69118539619979557642E+01},

            { 7,     0,    -7.32038980764418244362E+07},
            { 7,     1,     2.00004937755068242550E+08},
            { 7,     2,     1.87085510840455114841E+08},
            { 7,     3,    -2.40693202493527889252E+08},
            { 7,     4,     2.91464472802970111370E+07},
            { 7,     5,     1.80457588328361362219E+07},
            { 7,     6,    -5.22623080352768115699E+06},
            { 7,     7,     6.93950939136536209844E+05},
            { 7,     8,     1.29864879686568258330E+06},
            { 7,     9,    -6.70193084813015302643E+05},
            { 7,    10,     9.20749454564202897018E+04},
            { 7,    11,    -2.00941969415327257593E+04},
            { 7,    12,     4.96906965956355270464E+03},
            { 7,    13,    -3.33430820278479188801E+02},
            { 7,    14,     2.60022419015300613410E+02},
            { 7,    15,    -5.26054399949432109906E+01},

            { 8,     0,     3.49899361381577372551E+08},
            { 8,     1,    -6.16363983678707957268E+08},
            { 8,     2,     2.87744997564267575741E+08},
            { 8,     3,    -5.80734918315911106765E+06},
            { 8,     4,    -7.49858021421681251377E+06},
            { 8,     5,    -1.02456923512990288436E+07},
            { 8,     6,     5.04230622996888868511E+06},
            { 8,     7,     8.30524694711205520434E+04},
            { 8,     8,    -8.61025903787584858947E+05},
            { 8,     9,     1.34900060086739453254E+05},
            { 8,    10,     8.92692328603777714306E+04},
            { 8,    11,    -2.89001283128855830000E+04},
            { 8,    12,     2.12396475300008751219E+03},
            { 8,    13,     2.95241329416678070174E+02},
            { 8,    14,    -1.65338632433199251182E+02},
            { 8,    15,     2.41591961058125228590E+01},

            { 9,     0,    -1.71260510420660376549E+08},
            { 9,     1,     3.28761719424411475658E+08},
            { 9,     2,    -2.11206060614314973354E+08},
            { 9,     3,     3.46205389136431142688E+07},
            { 9,     4,     2.04971128339806981385E+07},
            { 9,     5,    -1.15401356420353483409E+07},
            { 9,     6,     2.01258777025233558379E+06},
            { 9,     7,    -1.79581004342673142673E+05},
            { 9,     8,     2.07729967156824888662E+05},
            { 9,     9,    -9.74314817908821714809E+04},
            { 9,    10,     2.22068544668985268800E+04},
            { 9,    11,    -8.59795432872551464243E+03},
            { 9,    12,     3.31306345412375048909E+03},
            { 9,    13,    -7.04637490543811736643E+02},
            { 9,    14,     9.43370988450807317349E+01},
            { 9,    15,    -7.12454991962896322377E+00},

            {10,     0,     2.84390209118411280215E+07},
            {10,     1,    -5.95417063166125416756E+07},
            {10,     2,     4.45506571666083186865E+07},
            {10,     3,    -9.78308669950580969453E+06},
            {10,     4,    -5.62624247386193368584E+06},
            {10,     5,     4.15237079235649015754E+06},
            {10,     6,    -5.87503508675953838974E+05},
            {10,     7,    -4.26206334621473739389E+05},
            {10,     8,     2.55988526607063336996E+05},
            {10,     9,    -6.22110169340419161017E+04},
            {10,    10,     3.20104591341956165707E+03},
            {10,    11,     3.22804174274844626780E+03},
            {10,    12,    -1.24817289214048719259E+03},
            {10,    13,     2.26605355126834183466E+02},
            {10,    14,    -2.29672306377457609017E+01},
            {10,    15,     1.15281981127687060962E+00}
        },
        {                                                                               // RECOM (Z003000)
            { 0,     0,     1.34060093558913830947E+01},
            { 0,     1,    -7.67393959247841483950E-01},
            { 0,     2,     5.67829538231681030247E-01},
            { 0,     3,    -1.53301717168738610431E-01},

            { 1,     0,     1.68198984963847353313E+00},
            { 1,     1,     2.75310597472883875070E-02},
            { 1,     2,    -3.85904737173993372945E-01},
            { 1,     3,     1.74248736565002115828E-01},

            { 2,     0,    -1.42654332247593096383E+00},
            { 2,     1,     1.33518023284770448456E+00},
            { 2,     2,    -3.89425375461592448989E-01},
            { 2,     3,    -1.47490049820737292169E-02},

            { 3,     0,     5.53979839492381387345E-01},
            { 3,     1,    -6.97392793919068609831E-01},
            { 3,     2,     2.81780811937450081928E-01},
            { 3,     3,    -2.63299277947934735888E-02}
        }
    }
};


// Coefficients for determining Main Sequence core mass
// from Shikauchi et al. (2024), https://arxiv.org/abs/2409.00460
// Table 2
const std::vector<DBL_VECTOR> SHIKAUCHI_ALPHA_COEFFICIENTS = {
    {0.45, -0.0557105,  -0.86589929},       // 0.1*Z_Sun
    {0.45, -0.06968022, -0.73688164},       // 1/3*Z_Sun
    {0.45, -0.05878711, -0.84646162}        // Solar metallicity Z_Sun
};
// Table 3
const std::vector<DBL_VECTOR> SHIKAUCHI_FMIX_COEFFICIENTS = {
    {0.86914766, -0.60815098, 37.20654856},     // 0.1*Z_Sun
    {0.86269445, -0.62623353, 35.74630996},     // 1/3*Z_Sun
    {0.86605495, -0.64960375, 35.57019104}      // Solar metallicity Z_Sun
};
// Table 4
const std::vector<DBL_VECTOR> SHIKAUCHI_L_COEFFICIENTS = {
    {3.2555795,  1.84666823, -0.79986388, -0.75728099, -0.38831172, 0.08223542, 0.49543834, 0.31314176, -0.36705796, 1.72200581},   // 0.1*Z_Sun
    {3.35622529, 1.96904931, -0.88894808, -0.81112488, -0.47925922, 0.09056925, 0.53094768, 0.33971972, -0.35581284, 1.65390003},   // 1/3*Z_Sun
    {3.27883249, 1.79370338, -0.71413866, -0.77019351, -0.3898752,  0.07499563, 0.5920458,  0.33846556, -0.49649838, 1.71263853}    // Solar metallicity Z_Sun
};

#endif // __constants_h__
