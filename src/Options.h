/******************************************************************************************/
/*                                                                                        */
/*                See Options.cpp for instructions for adding a new option                */
/*                                                                                        */
/******************************************************************************************/


#ifndef __Options_H__
#define __Options_H__

#define OPTIONS Options::Instance()

#include <iostream>
#include <string>
#include <sstream>
#include <typeinfo>
#include <typeindex>
#include <iterator>

#include "constants.h"

#include <boost/algorithm/string.hpp>   // Boost string manipulation
#include <boost/program_options.hpp>    // Boost command line options tools
#include <boost/filesystem.hpp>         // Boost filesystem tools for handling paths etc.

#include <boost/any.hpp>

#include "typedefs.h"
#include "ErrorCatalog.h"
#include "profiling.h"
#include "utils.h"
#include "Rand.h"
#include "changelog.h"

namespace po = boost::program_options;


const std::string NOT_PROVIDED = std::to_string(255);


// OPT_VALUE macro
//
// Getter functions return the value of the class member variable - the class
// member variable is set to a value depending upon the value of the corresponding
// option entered by the user.
// 
// Since users now specify grid line values using options, getter functions need to
// know which option value to return - the one specified on the commandline (if in
// fact the option was specified on the commandline), or the one specified on the
// grid line (if in fact the option was specified on the grid line).
//
// The general idea is to use the value specified by the user on the grid line (if
// the use actually specified the option on the grid line) in preference to the
// value specified by the the use on the commandline (if the use actually specified
// the option on the commandline).  That's what the OPT_VALUE macro defined below
// does - if the grid line exists (i.e. if a grid file is being used), the macro will
// check whether the user specified the option on the grid line, and if they did return
// that value, and if the option was not specified on the grid line, it will return
// either the default value for the option (if the 'fallback' parameter is 'false'),
// or the commandline value (if the 'fallback' parameter is 'true').  Note that the
// commandline value for an option will be the value specified by the user on the 
// commandline if in fact the option was specified on the commandline, and it will be
// the default value for the option if the option was not specified on the commandline.
//
// To reiterate: by using the OPT_VALUE macro, the value of the option returned will
// be (in order of priority):
//
//    1. the value specified on the grid line IFF the user specified the option on the 
//       grid line
//
//    2. if 'fallback' is 'true':
//           the value specified on the commandline IFF the user did not specify the 
//           option on the grid line but did specify the option on the commandline
//       else if 'fallback' is 'false':
//           the default value for the option
//
//    3. the default value for the option
//
// For most options we will use
//
//     OPT_VALUE("option-name", m_OptionVar, true)
//
// so that we fallback to the commandline value for an option if the user did not
// specify the option on the grid line, but there may be options for which we may not
// want to fall back to the commandline value, even if the user did not specify a value
// on the grid line - we may instead want to return the default value for the option 
// rather than return the value specified on the commandline (if any).  For those 
// options, the getter should use
//
//     OPT_VALUE("option-name", m_OptionVar, false)
//
// OPT_VALUE relies on the option name ("option-name") being exactly as specified in
// Options::AddOptions() - if it isn't, it will fail
//
// Finally, there are some options that we may always want to return the value specified
// on the commandline.  For those options, the getter should return
//
//     m_CmdLine.optionValues.m_ClassMemberVariable
//
// In that case the value entered by the user on the commandline will be returned IFF
// the user specified the option on the commandline, otherwise the default value for
// the option will be returned.
//
// Note that the `optName` argument to this macros should be the option name without
// the "-" or "--" prefix.


#define OPT_VALUE(optName, optValue, fallback)  ((m_GridLine.optionValues.m_Populated && \
                                                 (!m_GridLine.optionValues.m_VM[optName].defaulted() || !fallback)) \
                                                    ? m_GridLine.optionValues.optValue \
                                                    : m_CmdLine.optionValues.optValue)


// OPT_DEFAULTED macro
//
// Use to determine, for a given option, whether an option value was specified or it defaulted
// to the COMPAS default.  Note that this is different from the question of whether the option
// value is equal to the COMPAS default - this macro indicates how the option value was set:
// was it provided by the user, or was it set to the default value because a value was not
// provided by the user.
//
// This is a reasonable proxy for the Option::OptionSpecified() function, but only if the
// `optName` argument is actually a valid option name.
//
// The boost defaulted() function will return:
//
//     (a) TRUE  if `optName` is a valid option name and a value was not specified by the user
//         (so Option::OptionSpecified() would return FALSE)
//     (b) FALSE if `optName` is a valid option name and a value was specified by the user
//         (so Option::OptionSpecified() would return TRUE)
//     (c) FALSE if `optName` is not found in the stored list of valid option names
//         (i.e not a valid option name)
//
// In the (a) and (b) cases the boost defaulted() function (and so this macro) is a valid proxy for 
// Option::OptionSpecified(), but in the (c) case, while the result is technically correct (i.e the 
// default value was not set for the option), it is not a valid proxy for Option::OptionSpecified()
// (in this case, Option::OptionSpecified() would return FALSE)
//
// When `optName` is known to be a valid option name, this macro is a valid proxy for, and is ~7 times
// faster than, Option::OptionSpecified() (because it doesn't need to search the list of valid option
// names for `optName`): just take the NOT of the OPT_DEFAULTED macro.
//
// Note that the `optName` argument to this macros should be the option name without the "-" or "--" prefix.


#define OPT_DEFAULTED(optName) (m_GridLine.optionValues.m_Populated \
                                ? m_GridLine.optionValues.m_VM[optName].defaulted() \
                                : m_CmdLine.optionValues.m_VM[optName].defaulted())


/*
 * Options Singleton
 *
 * Holds program options
 *
 * Singletons are sometimes frowned-upon, but doing it this way means
 * the program options don't need to be passed around to all and sundry.
 * I think convenience and clarity sometimes trump dogma.
 */

class Options {

private:


    // The following vectors are used to specify deprecated option strings, option values,
    // and their replacements (if applicable).
    //
    // The vectors below need to be updated whenever we deprecate an option or an option value,
    // and the option (or value) eventually removed when the deprecation notice period is over.
    // 
    //
    // "deprecatedOptionStrings" vector
    // --------------------------------
    //
    // Each tuple in the "deprecatedOptionStrings" vector records an option that has been deprecated,
    // but is still available for users to specify.  The tuple entries are:
    // 
    //     - the option string for the deprecated option (just the option string - no leading "--")
    //     - the option string for any replacement for the deprecated option (just the option string,
    //       no leading "--").  If there is no replacement (i.e. the deprecated option will be removed
    //       and no replacement option implemented), the replacement option string should be the empty
    //       string ("")
    //     - a boolean flag to indicate if the deprecation notice for the option has been shown - should
    //       be "false" in the vector, and will be set true if and when the deprecation notice for that
    //       option is shown the first time in a COMPAS run (a deprecation notice for a deprecated option
    //       is only shown once per COMPAS run).
    // 
    //
    // "deprecatedOptionValues" vector
    // -------------------------------
    //
    // Sometimes we may want to deprecate an option value (e.g. one of the possible mass loss prescriptions).
    // We may want to do this to rename an option value, or we might want to remove it completely (without
    // replacement).
    // 
    // Each tuple in the "deprecatedOptionValues" vector records an option value that has been deprecated,
    // but is still available for users to specify.  The tuple entries are:
    // 
    //     - the option string for which a value is to be deprecated (just the option string - no leading "--")
    //     - the value string for the value to be deprecated (e.g. for the value QCRIT_PRESCRIPTION::CLAEYS for
    //       the option "critical-mass-ratio-prescription", specify "CLAEYS" in the vector)
    //     - the value string for any replacement value for the deprecated value (e.g. if the value
    //       QCRIT_PRESCRIPTION::CLAEYS for the option "critical-mass-ratio-prescription", is to be replaced
    //       with QCRIT_PRESCRIPTION::CLAEYS123, specify "CLAEYS123" in the vector).  If there is no replacement
    //       (i.e. the deprecated value will be removed and no replacement value implemented), the replacement
    //       value string should be the empty string ("")
    //     - a boolean flag to indicate if the deprecation notice for the option value has been shown - should
    //       be "false" in the vector, and will be set true if and when the deprecation notice for that option
    //       value is shown the first time in a COMPAS run (a deprecation notice for a deprecated option value
    //       is only shown once per COMPAS run).

    std::vector<std::tuple<std::string, std::string, bool>> deprecatedOptionStrings = {
        { "black-hole-kicks",                            "black-hole-kicks-mode",                           false },
        { "chemically-homogeneous-evolution",            "chemically-homogeneous-evolution-mode",           false },
        { "kick-direction",                              "kick-direction-distribution",                     false },
        { "luminous-blue-variable-prescription",         "LBV-mass-loss-prescription",                      false },
        { "mass-transfer",                               "use-mass-transfer",                               false },
        { "mass-transfer-thermal-limit-accretor",        "mass-transfer-thermal-limit-accretor-multiplier", false },
        { "OB-mass-loss",                                "OB-mass-loss-prescription",                       false },
        { "retain-core-mass-during-caseA-mass-transfer", "",                                                false },
        { "RSG-mass-loss",                               "RSG-mass-loss-prescription",                      false },
        { "VMS-mass-loss",                               "VMS-mass-loss-prescription",                      false },
        { "WR-mass-loss",                                "WR-mass-loss-prescription",                       false }
    };

    std::vector<std::tuple<std::string, std::string, std::string, bool>> deprecatedOptionValues = {
        { "critical-mass-ratio-prescription",    "GE20", "GE", false },
        { "critical-mass-ratio-prescription",    "GE20_IC", "GE_IC", false },
        { "LBV-mass-loss-prescription",          "NONE", "ZERO", false },
        { "luminous-blue-variable-prescription", "NONE", "ZERO", false },
        { "pulsational-pair-instability-prescription", "COMPAS", "WOOSLEY", false},
        { "OB-mass-loss",                        "NONE", "ZERO", false },
        { "OB-mass-loss-prescription",           "NONE", "ZERO", false },
        { "RSG-mass-loss",                       "NONE", "ZERO", false },
        { "RSG-mass-loss-prescription",          "NONE", "ZERO", false },
        { "VMS-mass-loss",                       "NONE", "ZERO", false },
        { "VMS-mass-loss-prescription",          "NONE", "ZERO", false },
        { "WR-mass-loss",                        "NONE", "ZERO", false },
        { "WR-mass-loss-prescription",           "NONE", "ZERO", false }
    };

    // the following vector is used to replace deprecated options in the logfile-definitions file
    std::vector<std::tuple<std::string, std::string, bool>> deprecatedOptionProperties = {
        { "black_hole_kicks", "black_hole_kicks_mode",      false },
        { "lbv_prescription", "LBV-mass-loss-prescription", false }
    };


    // The following vectors are used to constrain which options can be specified
    // when:
    //
    // m_ShorthandAllowed records option strings that may be specified using the
    // shorthand notation described in Options::PreprocessOptionValues()
    //
    // m_GridLineExcluded records option strings that may not be specified on a grid line
    //
    // m_SSEOnly records option strings that apply to SSE only
    // m_BSEOnly records option strings that apply to BSE only
    //
    // m_RangeExcluded records option strings for which a range may not be specified
    // m_SetExcluded records option strings for which a set may not be specified
    //
    // Each of these is described in more detail below


    // m_ShorthandAllowed records option strings that may be specified using the
    // shorthand notation described in Options::ExpandShorthandOptionValues() and
    // in Log.h.
    //
    // Furthermore, the vector records whether option values for such options can be
    // defaulted (i.e. some or all values need not be specified), and if so, what
    // string should be substituted for the unspecified values (at this stage there is
    // no type associated with the options - what is being manipulated is the command-line
    // or grid-line string that will be parsed to determine options specified and their
    // values).
    //
    // This vector is checked immediately prior to the command line or grid line being parsed
    // for ranges and sets - i.e. before the command line or grid line is passed to boost
    // for final parsing.  If any option strings not in this vector are specified using
    // shorthand notation, boost will parse them as usual and likely (though not necessarily)
    // complain (boost will only complain if the option/value pair is malformed or unknown,
    // which would almost certainly be the case - but it isn't guaranteed to be). 

    typedef std::tuple<std::string, bool, std::string> SHORTHAND_ENTRY;         // option name, default allowed (i.e. can be omitted), default string
    std::vector<SHORTHAND_ENTRY> m_ShorthandAllowed = {

        // trying to keep entries alphabetical so easier to find specific entries

        // option name          default allowed     default string
        { "debug-classes",      false,              "" },                       // don't allow defaults - we don't know how many classes to specify

        { "log-classes",        false,              "" },                       // don't allow defaults - we don't know how many classes to specify

        { "notes",              true,               "" },                       // allow defaults - number of notes is 0..#notes-hdrs
        { "notes-hdrs",         false,              "" }                        // don't allow defaults - we don't know how many headers to specify
    };


    // m_GridLineExcluded records option strings that may not be specified on a grid line
    //
    // This vector is checked when the grid line is parsed - if any option strings are 
    // specified but excluded from the grid line, a warning will be issued and the option
    // ignored (processing will continue). The reasons we might exclude options from the 
    // grid file should be obvious upon reading the list of excluded options.  I can't 
    // think of a good reason to exclude options from the commandline, so I haven't 
    // implemented that functionality (though it wouldn't be too difficult to add it).
    //
    // I could probably have done this using a different set of options in Boost for
    // the commandline and gridfile, but in the end I decided this way was actually
    // easier, cleaner, and gives us a bit more control.

    std::vector<std::string> m_GridLineExcluded = {

        // trying to keep entries alphabetical so easier to find specific entries

        "add-options-to-sysparms",

        "create-yaml-file",

        "debug-level",
        "debug-classes",
        "debug-to-file",
        "detailed-output",

        "enable-warnings",
        "errors-to-file",

        "fp-error-mode",

        "grid",
        "grid-start-line",
        "grid-num-lines",

        "hdf5-buffer-size",
        "hdf5-chunk-size",
        "help", "h",
        "hmxr-binaries",

        "log-level", 
        "log-classes",

        "logfile-common-envelopes",
        "logfile-common-envelopes-record-types",
        "logfile-definitions",
        "logfile-detailed-output",
        "logfile-detailed-output-record-types",
        "logfile-double-compact-objects",
        "logfile-double-compact-objects-record-types",
        "logfile-name-prefix",
        "logfile-pulsar-evolution",
        "logfile-pulsar-evolution-record-types",
        "logfile-rlof-parameters",
        "logfile-rlof-parameters-record-types",
        "logfile-supernovae",
        "logfile-supernovae-record-types",
        "logfile-switch-log",
        "logfile-system-parameters",
        "logfile-system-parameters-record-types",
        "logfile-type",

        "mode",

        "notes-hdrs",
        "number-of-systems",

        "output-container", "c",
        "outputPath", "o",

        "population-data-printing",
        "print-bool-as-string",

        "quiet", 

        "rlof-printing",

        "store-input-files",
        "switch-log",

        "timestep-multiplier",

        "version", "v",

        "yaml-template"
    };

    
    // m_SSEOnly records option strings that apply to SSE only
    // m_BSEOnly records option strings that apply to BSE only
    //
    // These vectors are checked when the commandline or grid line is parsed for
    // ranges and sets.  Ranges and sets are played out, and stars/binaries evolved
    // based on the grid of options defined by any ranges and sets specified by the
    // user.
    //
    // A problem arises when a user is (say) evolving single stars (in SSE mode)
    // and (perhaps inadvertently) specifies a range or set for an option that
    // applies only to BSE.  Before ranges and sets were implemented, this was not
    // a problem - the SSE code just ignores any BSE-only options (and vice-versa).
    // But with ranges and sets implemented, we need to know whether we should play
    // out the range or set.  If a range for BSE-only only option is specified the
    // SSE code will happily ignore the option, but unless we know not to, we will
    // still play out the range of values (only for them to be ignored) - so we will
    // evolve as many stars as there are values in the range, and they will all be
    // the same (because the SSE code will ignore the BSE-only option each time 
    // through the loop while the range is playing out).
    //
    // To get around this we specify in the following vectors the names of any
    // options that are SSE only or BSE only - that way we can choose not to play
    // them out as required.
    //
    // It's not the end of the world if we forget to put some entries in these
    // vectors - the worst thing that will happen is that duplicate stars/binaries
    // will be evolved as ranges/sets of options that will be ignored are played out.
    // In that case the user has a simple remedy: don't specify ranges for BSE options
    // when evolving in SSE mode (and vice-versa).  All we're doing here with these
    // vectors is helping the user avoid duplicating stars/binaries if they specify
    // inconsistent options.

    std::vector<std::string> m_SSEOnly = {

        // trying to keep enties alphabetical so easier to find specific entries

        "initial-mass",

        "kick-magnitude",
        "kick-magnitude-random",

        "rotational-frequency"
    };

    std::vector<std::string> m_BSEOnly = {

        // trying to keep entries alphabetical so easier to find specific entries

        "allow-rlof-at-birth",
        "allow-touching-at-birth",
        "angular-momentum-conservation-during-circularisation", 

        "case-BB-stability-prescription",
        "circularise-binary-during-mass-transfer",
        "common-envelope-allow-main-sequence-survive",
        "common-envelope-alpha", 
        "common-envelope-alpha-thermal",
        "common-envelope-formalism",
        "common-envelope-lambda",
        "common-envelope-lambda-multiplier",
        "common-envelope-lambda-prescription",
        "common-envelope-mass-accretion-constant",
        "common-envelope-mass-accretion-max",
        "common-envelope-mass-accretion-min",
        "common-envelope-mass-accretion-prescription",
        "common-envelope-recombination-energy-density",
        "common-envelope-slope-kruckow",

        "eccentricity", "e",
        "eccentricity-distribution",
        "eccentricity-max",
        "eccentricity-min",
        "emit-gravitational-radiation",
        "enhance-CHE-lifetimes-luminosities",
        "evolve-double-white-dwarfs",
        "evolve-pulsars",
        "evolve-unbound-systems",

        "include-WD-binaries-as-DCO",
        "initial-mass-1",
        "initial-mass-2",

        "kick-magnitude-1",
        "kick-magnitude-2",
        "kick-magnitude-random-1",
        "kick-magnitude-random-2",
        "kick-mean-anomaly-1",
        "kick-mean-anomaly-2",
        "kick-phi-1",
        "kick-phi-2",
        "kick-theta-1",
        "kick-theta-2",

        "logfile-common-envelopes",
        "logfile-common-envelopes-record-types",
        "logfile-double-compact-objects",
        "logfile-double-compact-objects-record-types",
        "logfile-pulsar-evolution",
        "logfile-pulsar-evolution-record-types",
        "logfile-rlof-parameters",
        "logfile-rlof-parameters-record-types",

        "mass-ratio", "q",
        "mass-ratio-max",
        "mass-ratio-min",
        "mass-ratio-distribution",
        "mass-transfer",
        "mass-transfer-fa",
        "mass-transfer-jloss",
        "mass-transfer-jloss-macleod-linear-fraction-degen",
        "mass-transfer-jloss-macleod-linear-fraction-non-degen",
        "mass-transfer-accretion-efficiency-prescription",
        "mass-transfer-angular-momentum-loss-prescription",
        "mass-transfer-rejuvenation-prescription",
        "mass-transfer-thermal-limit-accretor-multiplier",
        "mass-transfer-thermal-limit-C",
        "maximum-mass-donor-nandez-ivanova",
        "minimum-secondary-mass",

        "orbital-period",
        "orbital-period-distribution",
        "orbital-period-max",
        "orbital-period-min",

        "rlof-printing",
        "rotational-frequency-1",
        "rotational-frequency-2",

        "rocket-kick-magnitude-1",
        "rocket-kick-magnitude-2",
        "rocket-kick-phi-1", 
        "rocket-kick-phi-2", 
        "rocket-kick-theta-1",
        "rocket-kick-theta-2",

        "scale-CHE-mass-loss-with-surface-helium-abundance",
        "semi-major-axis", "a",
        "semi-major-axis-distribution",
        "semi-major-axis-max",
        "semi-major-axis-min",
    };

    
    // m_RangeExcluded records option strings that cannot be ranges
    // m_SetExcluded records option strings that cannot be sets
    //
    // These vectors are checked when the commandline or grid line is parsed for
    // ranges and sets.  Ranges can only be specified for numerical options - other
    // data types are not ordered, so ranges don't make sense (what would the
    // increment be...).  Sets can be specified for options of all data types,
    // but sets (and ranges) don't make sense for some options (things like "help",
    // "quiet", logfile names etc....)
  

    std::vector<std::string> m_RangeExcluded = {

        // trying to keep entries alphabetical so easier to find specific entries

        "add-options-to-sysparms",

        "allow-non-stripped-ECSN"
        "allow-rlof-at-birth",
        "allow-touching-at-birth",
        "angular-momentum-conservation-during-circularisation",

        "black-hole-kicks-mode",

        "case-BB-stability-prescription",
        "check-photon-tiring-limit",
        "chemically-homogeneous-evolution-mode",
        "circularise-binary-during-mass-transfer",
        "common-envelope-allow-main-sequence-survive",
        "common-envelope-formalism",
        "common-envelope-lambda-prescription",
        "common-envelope-mass-accretion-prescription",
        "create-yaml-file",
        "critical-mass-ratio-prescription",

        "debug-classes",
        "debug-level",
        "debug-to-file",
        "detailed-output",

        "eccentricity-distribution",
        "emit-gravitational-radiation",
        "enable-warnings",
        "enable-rotationally-enhanced-mass-loss",
        "enhance-CHE-lifetimes-luminosities",
        "envelope-state-prescription",
        "errors-to-file",
        "evolve-double-white-dwarfs"
        "evolve-pulsars",
        "evolve-unbound-systems",

        "fp-error-mode",
        "fryer-supernova-engine",

        "grid",
        "grid-start-line",
        "grid-num-lines",

        "hdf5-buffer-size",
        "hdf5-chunk-size",
        "help", "h",
        "hmxr-binaries",

        "include-WD-binaries-as-DCO",
        "initial-mass-function", "i",

        "kick-direction-distribution",
        "kick-magnitude-distribution", 

        "LBV-mass-loss-prescription",

        "log-level", 
        "log-classes",

        "logfile-common-envelopes",
        "logfile-common-envelopes-record-types",
        "logfile-definitions",
        "logfile-detailed-output",
        "logfile-double-compact-objects",
        "logfile-double-compact-objects-record-types",
        "logfile-name-prefix",
        "logfile-pulsar-evolution",
        "logfile-pulsar-evolution-record-types",
        "logfile-rlof-parameters",
        "logfile-rlof-parameters-record-types",
        "logfile-supernovae",
        "logfile-supernovae-record-types",
        "logfile-switch-log",
        "logfile-system-parameters",
        "logfile-system-parameters-record-types",
        "logfile-type",

        "main-sequence-core-mass-prescription",
        "mass-change-fraction",
        "mass-loss-prescription",
        "mass-ratio-distribution",
        "mass-transfer-accretion-efficiency-prescription",
        "mass-transfer-angular-momentum-loss-prescription",
        "mass-transfer-rejuvenation-prescription",
        "mass-transfer-thermal-limit-accretor-multiplier",
        "metallicity-distribution",
        "mode",

        "natal-kick-for-PPISN",
        "notes",
        "notes-hdrs",
        "neutrino-mass-loss-BH-formation",
        "neutron-star-equation-of-state",

        "OB-mass-loss-prescription",
        "orbital-period-distribution",
        "output-container", "c",
        "outputPath", "o",

        "pair-instability-supernovae",
        "population-data-printing",
        "print-bool-as-string",
        "pulsar-birth-magnetic-field-distribution",
        "pulsar-birth-spin-period-distribution",
        "pulsational-pair-instability",
        "pulsational-pair-instability-prescription",

        "quiet", 

        "RSG-mass-loss-prescription",
        "radial-change-fraction",
        "random-seed",
        "remnant-mass-prescription",
        "revised-energy-formalism-nandez-ivanova",
        "rlof-printing",
        "rotational-velocity-distribution",

        "scale-CHE-mass-loss-with-surface-helium-abundance",
        "semi-major-axis-distribution",
        "stellar-zeta-prescription",
        "store-input-files",
        "switch-log",

        "tides-prescription",

        "timesteps-filename",

        "use-mass-loss",
        "use-mass-transfer",

        "VMW-mass-loss-prescription",
        "version", "v",

        "WR-mass-loss-prescription",

        "yaml-template"
    };
    
    std::vector<std::string> m_SetExcluded = {

        // trying to keep entries alphabetical so easier to find specific entries

        "add-options-to-sysparms",

        "create-yaml-file",

        "debug-classes",
        "debug-level",
        "debug-to-file",
        "detailed-output",

        "enable-warnings",
        "errors-to-file",

        "fp-error-mode",

        "grid",
        "grid-start-line",
        "grid-num-lines",

        "hdf5-buffer-size",
        "hdf5-chunk-size",
        "help", "h",
        "hmxr-binaries",

        "log-classes",
        "log-level", 

        "logfile-common-envelopes",
        "logfile-common-envelopes-record-types",
        "logfile-definitions",
        "logfile-detailed-output",
        "logfile-detailed-output-record-types",
        "logfile-double-compact-objects",
        "logfile-double-compact-objects-record-types",
        "logfile-name-prefix",
        "logfile-pulsar-evolution",
        "logfile-pulsar-evolution-record-types",
        "logfile-rlof-parameters",
        "logfile-rlof-parameters-record-types",
        "logfile-supernovae",
        "logfile-supernovae-record-types",
        "logfile-switch-log",
        "logfile-system-parameters",
        "logfile-system-parameters-record-types",
        "logfile-type",

        "mass-change-fraction",
        "mode",

        "natal-kick-for-PPISN",
        "notes",
        "notes-hdrs",

        "output-container", "c",
        "outputPath", "o",

        "population-data-printing",
        "print-bool-as-string",

        "quiet",

        "radial-change-fraction",
        "random-seed",
        "rlof-printing",

        "store-input-files",
        "switch-log",

        "timesteps-filename",

        "version", "v",

        "yaml-template"
    };


public:
    
    // The OptionsValues class holds the values for the options.  This allows the Options class
    // to hold values for both the commandline options (the options specified by the user on the
    // commandline) and the grid file options (the options specified by the user on a grid file
    // record - on a per object (str/binary) basis).  When using grid files, the grid file options
    // take precedence over the commandline options: the grid file option values are set, then
    // options that were not specified in the grid file record are set from the command line
    // options.  Options that were not specified in the grid file record, and that were not
    // specified on the commandline, are set to the COMPAS default values.

    class OptionValues {

        friend class Options;                                                                                                   // So the Options class can access members directly

        private:

            template <typename T> 
            struct ENUM_OPT { 
                std::string typeString; 
                T           type; 
            }; 
            
            // member variables - alphabetically in groups (sort of...)

            bool                                                m_AllowNonStrippedECSN;                                         // Indicates whether single stars should undergo ECSNe if they were not stripped by a companion
            bool                                                m_AllowRLOFAtBirth;                                             // Indicates whether binaries that have one or both stars in RLOF at birth are allowed to evolve
            bool                                                m_AllowTouchingAtBirth;                                         // Indicates whether binaries that are touching at birth are allowed to evolve

            bool                                                m_DebugToFile;                                                  // Flag used to determine whether debug statements should also be written to a log file
            bool                                                m_ErrorsToFile;                                                 // Flag used to determine whether error statements should also be written to a log file

            bool                                                m_EnableWarnings;                                               // Flag used to determine if warnings (via SHOW_WARN macros) should be displayed
            ENUM_OPT<FP_ERROR_MODE>                             m_FPErrorMode;                                                  // Specifies the mode for floating-point error handling

            std::vector<std::string>                            m_Notes;                                                        // Notes contents - for user-defined annotations
            std::vector<std::string>                            m_NotesHdrs;                                                    // Notes header strings - for user-defined annotations

            bool                                                m_EvolveDoubleWhiteDwarfs;                                      // Whether to evolve double white dwarfs or not
            bool                                                m_EvolveMainSequenceMergers;                                    // Option to evolve binaries in which two stars merged on the main sequence
            bool                                                m_EvolvePulsars;                                                // Whether to evolve pulsars or not
	        bool                                                m_EvolveUnboundSystems;                                         // Option to chose if unbound systems are evolved until death or the evolution stops after the system is unbound during a SN.

            bool                                                m_EmitGravitationalRadiation;                                   // Option to emit gravitational radiation for each timestep of binary evolution

            bool                                                m_HMXRBinaries;                                                 // Flag if we want to store HMXRBs in RLOF output file
            bool                                                m_WDBinariesAsDCO;                                              // Flag if we want to store WD binariess in DCO output file

            bool                                                m_NatalKickForPPISN;                                            // Flag if PPISN remnant should receive a non-zero natal kick

            bool                                                m_DetailedOutput;                                               // Print detailed output details to file (default = false)
            bool                                                m_PopulationDataPrinting;                                       // Print certain data for small populations, but not for larger one
            bool                                                m_PrintBoolAsString;                                            // Flag used to indicate that boolean properties should be printed as "TRUE" or "FALSE" (default is 1 or 0)
            bool                                                m_Quiet;                                                        // Suppress some output
            bool                                                m_RlofPrinting;                                                 // RLOF printing

            bool                                                m_ShortHelp;                                                    // Flag to indicate whether user wants short help ('-h', just option names) or long help ('--help', plus descriptions)

            bool                                                m_StoreInputFiles;                                              // Store input files in output container (default = true)

            bool                                                m_SwitchLog;                                                    // Print switch log details to file (default = false)


            // Miscellaneous evolution variables

            ENUM_OPT<EVOLUTION_MODE>                            m_EvolutionMode;                                                // Mode of evolution: SSE or BSE

            int                                                 m_ObjectsToEvolve;                                              // Number of stars (SSE) or binaries (BSE) to evolve
            bool                                                m_FixedRandomSeed;                                              // Whether to use a fixed random seed given by options.randomSeed (set to true if --random-seed is passed on command line)
            unsigned long int                                   m_RandomSeed;                                                   // Random seed to use
    
            double                                              m_MaxEvolutionTime;                                             // Maximum time to evolve a binary by
            unsigned long int                                   m_MaxNumberOfTimestepIterations;                                // Maximum number of timesteps to evolve binary for before giving up
            double                                              m_TimestepMultiplier;                                           // Multiplier for time step size (<1 -- shorter timesteps, >1 -- longer timesteps)
   
            double m_MassChangeFraction;                                                                                        // Approximate goal for fractional radial change per timestep
            double m_RadialChangeFraction;                                                                                      // Approximate goal for fractional radial change per timestep

            std::streamsize                                     m_GridStartLine;                                                // The grid file line to start processing (0-based)
            std::streamsize                                     m_GridLinesToProcess;                                           // The number of grid file lines to process (starting at m_GridStartLine)

            std::string                                         m_TimestepsFileName;                                            // The name of the timesteps file

            // Initial distribution variables

            double                                              m_InitialMass;                                                  // Initial mass of single star (SSE)
            double                                              m_InitialMass1;                                                 // Initial mass of primary (BSE)
            double                                              m_InitialMass2;                                                 // Initial mass of secondary (BSE)

            ENUM_OPT<INITIAL_MASS_FUNCTION>                     m_InitialMassFunction;                                          // Which initial mass function
            double                                              m_InitialMassFunctionMin;                                       // Minimum mass to generate in Msol
            double                                              m_InitialMassFunctionMax;                                       // Maximum mass to generate in Msol
            double                                              m_InitialMassFunctionPower;                                     // single IMF power law set manually

            // Mass ratio
            double                                              m_MassRatio;                                                    // Mass ratio for BSE
            ENUM_OPT<MASS_RATIO_DISTRIBUTION>                   m_MassRatioDistribution;                                        // Which mass ratio distribution
            double                                              m_MassRatioDistributionMin;                                     // Minimum initial mass ratio when using a distribution
            double                                              m_MassRatioDistributionMax;                                     // Maximum initial mass ratio when using a distribution

            double                                              m_MinimumMassSecondary;                                         // Minimum mass of secondary to draw (in Msol)

            // Semi major axis
            double                                              m_SemiMajorAxis;                                                // Semi-major axis
            ENUM_OPT<SEMI_MAJOR_AXIS_DISTRIBUTION>              m_SemiMajorAxisDistribution;                                    // Which semi-major axis distribution
            double                                              m_SemiMajorAxisDistributionMin;                                 // Minimum a in AU
            double                                              m_SemiMajorAxisDistributionMax;                                 // Maximum a in AU

            // Orbital period
            double                                              m_OrbitalPeriod;                                                // Orbital period in days
            ENUM_OPT<ORBITAL_PERIOD_DISTRIBUTION>               m_OrbitalPeriodDistribution;                                    // Which orbital period distribution
            double                                              m_OrbitalPeriodDistributionMin;                                 // Minimum initial period in days
            double                                              m_OrbitalPeriodDistributionMax;                                 // Maximum initial period in days

            // Wind mass loss
            bool                                                m_EnableRotationallyEnhancedMassLoss;                           // Whether to enable rotationally enhanced mass loss
            double                                              m_CoolWindMassLossMultiplier;                                   // Multiplication factor to reduce cool wind mass loss rate at each timestep
            double                                              m_OverallWindMassLossMultiplier;                                // Multiplication factor to reduce the overall wind mass loss rate at each timestep
            double                                              m_ScaleTerminalWindVelocityWithMetallicityPower;                // Power with which to scale terminal wind velocity with metallicity (v_inf ~ Z^x)

            // Eccentricity
            double                                              m_Eccentricity;                                                 // Eccentricity
            ENUM_OPT<ECCENTRICITY_DISTRIBUTION>                 m_EccentricityDistribution;                                     // Which eccentricity distribution
            double                                              m_EccentricityDistributionMin;                                  // Minimum initial eccentricity when using a distribution
            double                                              m_EccentricityDistributionMax;                                  // Maximum initial eccentricity when using a distribution

            // Kick options
            ENUM_OPT<KICK_MAGNITUDE_DISTRIBUTION>               m_KickMagnitudeDistribution;                                    // Which kick magnitude distribution
            double                                              m_KickMagnitudeDistributionSigmaCCSN_NS;                        // Kick magnitude sigma in km s^-1 for neutron stars (default = "250" )
            double                                              m_KickMagnitudeDistributionSigmaCCSN_BH;                        // Kick magnitude sigma in km s^-1 for black holes (default = "250" )
            double                                              m_KickMagnitudeDistributionMaximum;                             // Maximum kick magnitude to draw. If negative, no maximum
	        double                                              m_KickMagnitudeDistributionSigmaForECSN;			            // Kick magnitude sigma for ECSN in km s^-1 (default = "0" )
	        double                                              m_KickMagnitudeDistributionSigmaForUSSN;			            // Kick magnitude sigma for USSN in km s^-1 (default = "20" )
	        double                                              m_KickScalingFactor;								            // Arbitrary factor for scaling kicks

            // Kick direction options
            ENUM_OPT<KICK_DIRECTION_DISTRIBUTION>               m_KickDirectionDistribution;                                    // Kick direction distribution
            double                                              m_KickDirectionPower;                                           // Exponent

            // User-specified supernova parameter values
            double                                              m_KickMagnitude;                                                // Supernova kick magnitude - SSE
            double                                              m_KickMagnitude1;                                               // Supernova kick magnitude - BSE primary star
            double                                              m_KickMagnitude2;                                               // Supernova kick magnitude - BSE secondary star
            double                                              m_KickMagnitudeRandom;                                          // Random number U(0,1) for choosing the supernova kick magnitude - SSE
            double                                              m_KickMagnitudeRandom1;                                         // Random number U(0,1) for choosing the supernova kick magnitude - BSE primary star
            double                                              m_KickMagnitudeRandom2;                                         // Random number U(0,1) for choosing the supernova kick magnitude - BSE secondary star
            double                                              m_KickMeanAnomaly1;                                             // Mean anomaly at instantaneous time of the SN - uniform in [0, 2pi] - BSE primary star
            double                                              m_KickMeanAnomaly2;                                             // Mean anomaly at instantaneous time of the SN - uniform in [0, 2pi] - BSE secondary star
            double                                              m_KickPhi1;                                                     // Angle between 'x' and 'y', both in the orbital plane of supernovae vector (rad) - BSE primary star
            double                                              m_KickPhi2;                                                     // Angle between 'x' and 'y', both in the orbital plane of supernovae vector (rad) - BSE secondary star
            double                                              m_KickTheta1;                                                   // Angle between the orbital plane and the 'z' axis of supernovae vector (rad) - BSE primary star
            double                                              m_KickTheta2;                                                   // Angle between the orbital plane and the 'z' axis of supernovae vector (rad) - BSE secondar star

            double                                              m_MullerMandelKickBH;                                           // Multiplier for BH kicks per Mandel and Mueller, 2020
            double                                              m_MullerMandelKickNS;                                           // Multiplier for NS kicks per Mandel and Mueller, 2020
            double                                              m_MullerMandelSigmaKick;                                        // Scatter for kicks per Mandel and Mueller, 2020

            // Black hole kicks
            ENUM_OPT<BLACK_HOLE_KICKS_MODE>                     m_BlackHoleKicksMode;                                           // Which black hole kicks mode

            // Rocket kicks
            double                                              m_RocketKickMagnitude1;                                         // Rocket kick magnitude primary - only for neutron stars
            double                                              m_RocketKickMagnitude2;                                         // Rocket kick magnitude secondary - only for neutron stars
            double                                              m_RocketKickPhi1;                                               // Rocket kick phi angle primary
            double                                              m_RocketKickPhi2;                                               // Rocket kick phi angle secondary
            double                                              m_RocketKickTheta1;                                             // Rocket kick theta angle primary
            double                                              m_RocketKickTheta2;                                             // Rocket kick theta angle secondary
                                                                                                                                
            // CHE - Chemically Homogeneous Evolution
            ENUM_OPT<CHE_MODE>                                  m_CheMode;                                                      // Which Chemically Homogeneous Evolution mode
            bool                                                m_EnhanceCHELifetimesLuminosities;                              // Whether to enhance the lifetimes and luminosities of CHE stars relative to SSE MS stars
            bool                                                m_ScaleCHEMassLossWithSurfaceHeliumAbundance;                   // Whether to transition between OB and WR mass loss rates for CHE stars on the MS

            // Supernova remnant mass
            ENUM_OPT<REMNANT_MASS_PRESCRIPTION>                 m_RemnantMassPrescription;                                      // Which remnant mass prescription

            ENUM_OPT<SN_ENGINE>                                 m_FryerSupernovaEngine;                                         // Which Fryer et al. supernova engine

            ENUM_OPT<NEUTRINO_MASS_LOSS_PRESCRIPTION>           m_NeutrinoMassLossAssumptionBH;                                 // Which neutrino mass loss assumption for BH formation
            double                                              m_NeutrinoMassLossValueBH;                                      // Value (corresponding to assumption) for neutrino mass loss for BH formation


            double                                              m_Fryer22fmix;                                                  // Parameter describing the mixing growth time when using Fryer 2022 remnant mass presc. 
            double                                              m_Fryer22Mcrit;                                                 // Critical mass for black hole formation when using Fryer 2022 remnant mass presc. 

            // Fixed uk options
            bool                                                m_UseFixedUK;                                                   // Whether to fix uk to a certain value (default is to NOT fix uk)
            double                                              m_FixedUK;                                                      // Dimensionless value to fix the kick magnitude to

            // Pair instability and pulsational pair instability mass loss
            bool                                                m_UsePairInstabilitySupernovae;                                 // Whether to use pair instability supernovae (PISN)
            double                                              m_PairInstabilityLowerLimit;                                    // Minimum core mass leading to PISN
            double                                              m_PairInstabilityUpperLimit;                                    // Maximum core mass leading to PISN

            bool                                                m_UsePulsationalPairInstability;                                // Whether to use pulsational pair instability (PPI)
            double                                              m_PulsationalPairInstabilityLowerLimit;                         // Maximum core mass leading to PPI
            double                                              m_PulsationalPairInstabilityUpperLimit;                         // Minimum core mass leading to PPI

            double                                              m_PulsationalPairInstabilityCOCoreShiftHendriks;                // Shift in CO Core mass for PPI from Hendriks+23
            
            ENUM_OPT<PPI_PRESCRIPTION>                          m_PulsationalPairInstabilityPrescription;                       // Which PPI prescription

	        double                                              m_MaximumNeutronStarMass;						                // Maximum mass of a neutron star allowed, set to default in StarTrack

            // Setup default output directory and desired output directory
            std::string                                         m_OutputPathString;                                             // String to hold the output directory
            std::string                                         m_OutputContainerName;                                          // Name of output container (directory)

            // Mass loss options
            bool                                                m_UseMassLoss;                                                  // Whether to activate mass loss (default = True)
            bool                                                m_CheckPhotonTiringLimit;                                       // Whether to check the photon tiring limit for wind mass loss

            // Can also have options for modifying strength of winds etc here

            ENUM_OPT<MASS_LOSS_PRESCRIPTION>                    m_MassLossPrescription;                                         // Which mass loss prescription

            ENUM_OPT<LBV_MASS_LOSS_PRESCRIPTION>                m_LBVMassLossPrescription;                                      // Which LBV mass loss prescription to use
            double                                              m_LuminousBlueVariableFactor;                                   // Multiplicitive factor for luminous blue variable (LBV) mass loss rates when using Belczynski’s prescription
            double                                              m_WolfRayetFactor;                                              // Multiplicitive factor for Wolf-Rayet (WR) wind mass loss rates

            ENUM_OPT<OB_MASS_LOSS_PRESCRIPTION>                 m_OBMassLossPrescription;                                       // Which OB mass loss prescrioption
            ENUM_OPT<VMS_MASS_LOSS_PRESCRIPTION>                m_VMSMassLossPrescription;                                      // Which VMS mass loss prescription for M > 100 Msol        
            ENUM_OPT<RSG_MASS_LOSS_PRESCRIPTION>                m_RSGMassLossPrescription;                                      // Which RSG mass loss prescription to use for RSG       
            ENUM_OPT<WR_MASS_LOSS_PRESCRIPTION>                 m_WRMassLossPrescription;                                       // Which WR mass loss prescription to use for WR       


            // Mass transfer options
            bool                                                m_UseMassTransfer;                                              // Whether to use mass transfer (default = true)
	        bool                                                m_CirculariseBinaryDuringMassTransfer;						    // Whether to circularise binary when it starts (default = true)
	        bool                                                m_AngularMomentumConservationDuringCircularisation;			    // Whether to conserve angular momentum while circularising or circularise to periastron (default = false)
            double                                              m_ConvectiveEnvelopeTemperatureThreshold;                       // The boundary between convective and radiative envelopes for HG and Giant stars
        
            bool                                                m_ExpelConvectiveEnvelopeAboveLuminosityThreshold;              // Whether to expel the convective envelope in a pulsation when log_10(L/M) reaches the threshold defined by m_LuminosityToMassThreshold
            double                                              m_LuminosityToMassThreshold;                                    // Threshold value of log_10(L/M) above which the convective envelope is expelled in a pulsation
        
            bool                                                m_RetainCoreMassDuringCaseAMassTransfer;                        // Whether to retain the approximate core mass of a case A donor as a minimum core at end of MS or HeMS (default = false)

            ENUM_OPT<CORE_MASS_PRESCRIPTION>                    m_MainSequenceCoreMassPrescription;                             // Which MS core prescription
        
            ENUM_OPT<CASE_BB_STABILITY_PRESCRIPTION>            m_CaseBBStabilityPrescription;									// Which prescription for the stability of case BB/BC mass transfer


            ENUM_OPT<MT_ACCRETION_EFFICIENCY_PRESCRIPTION>      m_MassTransferAccretionEfficiencyPrescription;                  // Which accretion efficiency prescription

            double                                              m_MassTransferFractionAccreted;                                 // In mass transfer, amount of mass transferred that is accreted. 1 for conservative, 0 for fully-non conservative.
            double                                              m_MassTransferCParameter;                                       // Detailed model parameter used in mass transfer
            double                                              m_EddingtonAccretionFactor;                                     // Multiplication factor for eddington accretion for NS & BH
                                                                                                                                // i.e. >1 is super-eddington
                                                                                                                                //       0. is no accretion

	        ENUM_OPT<MT_THERMALLY_LIMITED_VARIATION>            m_MassTransferThermallyLimitedVariation;                        // Choose how to deal with mass transfer if it is set as thermally limited.

            double                                              m_MassTransferJloss;                                            // Specific angular momentum of the material leaving the system (not accreted)
            double                                              m_MassTransferJlossMacLeodLinearFractionDegen;                  // Linear interpolation fraction for jloss for degenerate accretors, between accretor and L2 position 
            double                                              m_MassTransferJlossMacLeodLinearFractionNonDegen;               // Linear interpolation fraction for jloss for non-degenerate accretors, between accretor and L2 position 
            ENUM_OPT<MT_ANGULAR_MOMENTUM_LOSS_PRESCRIPTION>     m_MassTransferAngularMomentumLossPrescription;                  // Which mass transfer angular momentum loss prescription

            // Mass transfer rejuvenation prescription
            ENUM_OPT<MT_REJUVENATION_PRESCRIPTION>              m_MassTransferRejuvenationPrescription;                         // Which mass transfer rejuvenation prescription

            // Mass transfer critical mass ratios
            ENUM_OPT<QCRIT_PRESCRIPTION>                        m_QCritPrescription;                                            // The critical mass ratio prescription, if any
            double                                              m_MassTransferCriticalMassRatioMSLowMassNonDegenerateAccretor;  // Critical mass ratio for MT from a MS low mass star
            double                                              m_MassTransferCriticalMassRatioMSLowMassDegenerateAccretor;     // Critical mass ratio for MT from a MS low mass star on to a degenerate accretor

            double                                              m_MassTransferCriticalMassRatioMSHighMassNonDegenerateAccretor; // Critical mass ratio for MT from a MS high mass star
            double                                              m_MassTransferCriticalMassRatioMSHighMassDegenerateAccretor;    // Critical mass ratio for MT from a MS high mass star on to a degenerate accretor

            double                                              m_MassTransferCriticalMassRatioGiantNonDegenerateAccretor;      // Critical mass ratio for MT from a giant
            double                                              m_MassTransferCriticalMassRatioGiantDegenerateAccretor;         // Critical mass ratio for MT from a giant on to a degenerate accretor

            double                                              m_MassTransferCriticalMassRatioHGNonDegenerateAccretor;         // Critical mass ratio for MT from a HG star
            double                                              m_MassTransferCriticalMassRatioHGDegenerateAccretor;            // Critical mass ratio for MT from a HG star on to a degenerate accretor

            double                                              m_MassTransferCriticalMassRatioHeliumMSNonDegenerateAccretor;   // Critical mass ratio for MT from a Helium MS star
            double                                              m_MassTransferCriticalMassRatioHeliumMSDegenerateAccretor;      // Critical mass ratio for MT from a Helium MS star on to a degenerate accretor

            double                                              m_MassTransferCriticalMassRatioHeliumHGNonDegenerateAccretor;   // Critical mass ratio for MT from a Helium HG star
            double                                              m_MassTransferCriticalMassRatioHeliumHGDegenerateAccretor;      // Critical mass ratio for MT from a Helium HG star on to a degenerate accretor

            double                                              m_MassTransferCriticalMassRatioHeliumGiantNonDegenerateAccretor;// Critical mass ratio for MT from a helium giant
            double                                              m_MassTransferCriticalMassRatioHeliumGiantDegenerateAccretor;   // Critical mass ratio for MT from a helium giant on to a degenerate accretor

            double                                              m_MassTransferCriticalMassRatioWhiteDwarfNonDegenerateAccretor; // Critical mass ratio for MT from a white dwarf
            double                                              m_MassTransferCriticalMassRatioWhiteDwarfDegenerateAccretor;    // Critical mass ratio for MT from a white dwarf on to a degenerate accretor

            // Common Envelope options
            double                                              m_CommonEnvelopeAlpha;                                          // Common envelope efficiency alpha parameter (default = X)
            double                                              m_CommonEnvelopeLambda;                                         // Common envelope Lambda parameter (default = X)
	        double                                              m_CommonEnvelopeSlopeKruckow;									// Common envelope power factor for Kruckow fit normalized according to Kruckow+2016, Fig. 1
            double                                              m_CommonEnvelopeAlphaThermal;                                   // lambda = alpha_th*lambda_b + (1-alpha_th)*lambda_g
            double                                              m_CommonEnvelopeLambdaMultiplier;                               // Multiply common envelope lambda by some constant
            bool                                                m_CommonEnvelopeLambdaNanjingEnhanced;                          // Use Nanjing lambda's with enhanced extrapolation in stellar radius
            bool                                                m_CommonEnvelopeLambdaNanjingInterpolateInMass;                 // Use Nanjing lambda's with mass interpolation (only used when using enhanced Nanjing lambda's)
            bool                                                m_CommonEnvelopeLambdaNanjingInterpolateInMetallicity;          // Use Nanjing lambda's with metallicity interpolation (only used when using enhanced Nanjing lambda's)
            bool                                                m_CommonEnvelopeLambdaNanjingUseRejuvenatedMass;                // Whether or not to use mass after rejuvenation (m_Mass0) instead of true birth mass when calculating Nanjing lambda's
            bool                                                m_AllowMainSequenceStarToSurviveCommonEnvelope;                 // Whether or not to allow a main sequence star to survive a common envelope event
            bool                                                m_AllowRadiativeEnvelopeStarToSurviveCommonEnvelope;            // Whether or not to allow a radiative-envelope star to survive a common envelope event
            bool                                                m_AllowImmediateRLOFpostCEToSurviveCommonEnvelope;              // Whether or not to allow Roche Lobe Overflow immediately after a CE to survive a common envelope event
    
            // Prescription for envelope state (radiative or convective)
            ENUM_OPT<ENVELOPE_STATE_PRESCRIPTION>               m_EnvelopeStatePrescription;

            // Accretion during common envelope
            ENUM_OPT<CE_ACCRETION_PRESCRIPTION>                 m_CommonEnvelopeMassAccretionPrescription;
            double                                              m_CommonEnvelopeMassAccretionMin;
            double                                              m_CommonEnvelopeMassAccretionMax;
            double                                              m_CommonEnvelopeMassAccretionConstant;

            // Common envelope formalism
            ENUM_OPT<CE_FORMALISM>                              m_CommonEnvelopeFormalism;                                      // Formalism for CE evolution
        
	        // Common envelope lambda prescription
	        ENUM_OPT<CE_LAMBDA_PRESCRIPTION>                    m_CommonEnvelopeLambdaPrescription;							    // Prescription to use for CE lambda

	        // Common envelope Nandez and Ivanova energy formalism
	        bool                                                m_RevisedEnergyFormalismNandezIvanova;			                // Use the revised energy formalism from Nandez & Ivanova 2016 (default = false)
	        double                                              m_MaximumMassDonorNandezIvanova;								// Maximum mass allowed to use the revised energy formalism in Msol (default = 2.0)
	        double                                              m_CommonEnvelopeRecombinationEnergyDensity;					    // Factor using to calculate the binding energy depending on the mass of the envelope. (default = 1.5x10^13 erg/g)


            // Tides
            ENUM_OPT<TIDES_PRESCRIPTION>                        m_TidesPrescription;                                             // Which tides prescription (default = NONE)


            // Zetas
            ENUM_OPT<ZETA_PRESCRIPTION>                         m_StellarZetaPrescription;                                 	    // Prescription to use for calculating stellar zetas (default = SOBERMAN)

	        double                                              m_ZetaAdiabaticArbitrary;
	        double                                              m_ZetaMainSequence;
            double                                              m_ZetaRadiativeEnvelopeGiant;


            // Metallicity options
            double                                              m_Metallicity;                                                  // Metallicity
            ENUM_OPT<METALLICITY_DISTRIBUTION>                  m_MetallicityDistribution;                                      // Which metallicity distribution
            double                                              m_MetallicityDistributionMin;                                   // Minimum initial metallicity when using a distribution
            double                                              m_MetallicityDistributionMax;                                   // Maximum initial metallicity when using a distribution

            double                                              m_mCBUR1;                                                       // Minimum core mass at base of the AGB to avoid fully degenerate CO core formation


            // Neutron star equation of state
            ENUM_OPT<NS_EOS>                                    m_NeutronStarEquationOfState;                                   // NS EOS


            // Pulsar birth magnetic field distribution string
            ENUM_OPT<PULSAR_BIRTH_MAGNETIC_FIELD_DISTRIBUTION>  m_PulsarBirthMagneticFieldDistribution;                         // Birth magnetic field distribution for pulsars
            double                                              m_PulsarBirthMagneticFieldDistributionMin;                      // Minimum birth magnetic field (log10 B/G)
            double                                              m_PulsarBirthMagneticFieldDistributionMax;                      // Maximum birth magnetic field (log10 B/G)

            // Pulsar birth spin period distribution string
            ENUM_OPT<PULSAR_BIRTH_SPIN_PERIOD_DISTRIBUTION>     m_PulsarBirthSpinPeriodDistribution;                            // Birth spin period distribution for pulsars
            double                                              m_PulsarBirthSpinPeriodDistributionMin;                         // Minimum birth spin period (ms)
            double                                              m_PulsarBirthSpinPeriodDistributionMax;                         // Maximum birth spin period (ms)

            double                                              m_PulsarMagneticFieldDecayTimescale;                            // Timescale on which magnetic field decays (Myr)
            double                                              m_PulsarMagneticFieldDecayMassscale;                            // Mass scale on which magnetic field decays during accretion (solar masses)
            double                                              m_PulsarLog10MinimumMagneticField;                              // log10 of the minimum pulsar magnetic field in Gauss


            // Rotational Velocity distribution options
            ENUM_OPT<ROTATIONAL_VELOCITY_DISTRIBUTION>          m_RotationalVelocityDistribution;                               // Rotational velocity distribution
            double                                              m_RotationalFrequency;                                          // Rotational frequency for single star (SSE)
            double                                              m_RotationalFrequency1;                                         // Rotational frequency for primary (BSE)
            double                                              m_RotationalFrequency2;                                         // Rotational frequency for secondary (BSE)

	        // Grids

            std::string                                         m_GridFilename;                                                 // Grid filename


            // Debug and logging options

            int                                                 m_DebugLevel;                                                   // Debug level - used to determine which debug statements are actually written
            std::vector<std::string>                            m_DebugClasses;                                                 // Debug classes - used to determine which debug statements are actually written

            int                                                 m_LogLevel;                                                     // Logging level - used to determine which logging statements are actually written
            std::vector<std::string>                            m_LogClasses;                                                   // Logging classes - used to determine which logging statements are actually written


            // Logfiles
            std::string                                         m_LogfileDefinitionsFilename;                                   // Filename for the logfile record definitions
            std::string                                         m_LogfileNamePrefix;                                            // Prefix for log file names
            ENUM_OPT<LOGFILETYPE>                               m_LogfileType;                                                  // File type log files

            std::string                                         m_LogfileSystemParameters;                                      // output file name: system parameters
            std::string                                         m_LogfileDetailedOutput;                                        // output file name: detailed output
            std::string                                         m_LogfileDoubleCompactObjects;                                  // output file name: double compact objects
            std::string                                         m_LogfileSupernovae;                                            // output file name: supernovae
            std::string                                         m_LogfileCommonEnvelopes;                                       // output file name: common envelopes
            std::string                                         m_LogfileRLOFParameters;                                        // output file name: Roche Lobe overflow
            std::string                                         m_LogfilePulsarEvolution;                                       // output file name: pulsar evolution
            std::string                                         m_LogfileSwitchLog;                                             // output file name: switch log

            int                                                 m_LogfileSystemParametersRecordTypes;                           // enabled record types: system parameters
            int                                                 m_LogfileDetailedOutputRecordTypes;                             // enabled record types: detailed output
            int                                                 m_LogfileDoubleCompactObjectsRecordTypes;                       // enabled record types: double compact objects
            int                                                 m_LogfileSupernovaeRecordTypes;                                 // enabled record types: supernovae
            int                                                 m_LogfileCommonEnvelopesRecordTypes;                            // enabled record types: common envelopes
            int                                                 m_LogfileRLOFParametersRecordTypes;                             // enabled record types: Roche Lobe overflow
            int                                                 m_LogfilePulsarEvolutionRecordTypes;                            // enabled record types: pulsar evolution

            ENUM_OPT<ADD_OPTIONS_TO_SYSPARMS>                   m_AddOptionsToSysParms;                                         // Whether/when to add program option columns to BSE/SSE sysparms file

            int                                                 m_HDF5BufferSize;                                               // HDF5 file IO buffer size (number of chunks)
            int                                                 m_HDF5ChunkSize;                                                // HDF5 file chunk size (number of dataset entries)


            // YAML file

            std::string                                         m_YAMLfilename;                                                 // filename of YAML file to be created
            std::string                                         m_YAMLtemplate;                                                 // filename of user-supplied YAML template file


            // the boost variables map
            // this holds information on the options as specified by the user

            po::variables_map m_VM;

            bool m_Populated;                                                                                                   // flag to indicate whether we're using a grid line


            // member functions

            std::string CheckAndSetOptions();

            void        Initialise();

            template<class T>
            void ModifyVariableMap(std::map<std::string, po::variable_value>& vm, const std::string& opt, const T& val) { 
                vm[opt].value() = boost::any(val);
            }

            std::string SetCalculatedOptionDefaults(const BOOST_MAP p_BoostMap);

        public:

    };  // class OptionValues


    // complex option values are values for options that the user has supplied as ranges or sets
    //
    // complex option values are described by a tuple containing:
    //
    //     optionName       (std::string)               the name of the option
    //     complexValue     (RangeOrSetDescriptorT)     the complex option value - a RANGE or a SET
    //
    // ranges and sets are described by the RangeOrSetDescriptorT struct (see below)
    // the struct elements are described as:
    //
    //     type         (INT)                           type indicates whether the entry refers to a RANGE (type 0) or SET (type 1)
    //     dataType     (TYPENAME)                      the data type of the option to which the RangeOrSetDescriptorT pertaines
    //     parameters   (std::vector<std::string>)      a vector of strings that hold the parameters as they were supplied by the user
    //                                                  for a RANGE there must be exactly 3 parameters: start, count, increment
    //                                                  a SET must have at least one parameter (element); there is no maximum number of elements
    //     rangeParms   (std::vector<RangeParameterT>)  numerical values for range parameters (see RangeParameter struct)
    //     currPos      (INT)                           the current iterator position (the code iterates over the range or set)

    enum class COMPLEX_TYPE: int {NONE, RANGE, SET};

    typedef union RangeParameter {
        // these are the only types we need at the moment
        // we can add more if we ever need them
        long double   ldVal;    // LONG DOUBLE
        double        dVal;     // FLOAT/DOUBLE
        unsigned long ulVal;    // UNSIGNED LONG (INT)
        long          lVal;     // LONG (INT)
        int           iVal;     // INT
    } RangeParameterT; 

    typedef struct RangeOrSetDescriptor {
        COMPLEX_TYPE                 type;                                              // RANGE or SET
        TYPENAME                     dataType;                                          // the option datatype
        std::vector<std::string>     parameters;                                        // the range or set parameters
        std::vector<RangeParameterT> rangeParms;                                        // range parameters numerical values
        int                          currPos;                                           // current position of iterator - count for RANGE, pos for SET                                             
    } RangeOrSetDescriptorT;

    typedef std::vector<std::tuple<std::string, RangeOrSetDescriptorT>> COMPLEX_OPTION_VALUES;

    typedef std::tuple<TYPENAME, bool, std::string, std::string> ATTR;                  // <dataType, defaulted, typeStr, valueStr>

    typedef STR_STR_STR_STR OPTIONSTR;                                                  // option strings for specified options: <asEntered, asEnteredDownshifted, longName, shortName>

    // we have two structs:
    //    one for the commandline (program-level) options, and 
    //    one for the grid file line (evolving object-level) options
    //
    // each struct contains:
    //
    //    an OPTIONS_ORIGIN variable to indicate whether this struct is for command-line or grid file options (so the struct can be queried)
    //    an OptionValues object - holds the values of the options 
    //    a  Boost options_descriptions object
    //    a  COMPLEX_OPTION_VALUES object - holds the complex option values (ranges, sets)
    //    a  struct containing the option strings of the specified options

    typedef struct OptionsDescriptor {
        OPTIONS_ORIGIN          optionsOrigin;
        OptionValues            optionValues;
        po::options_description optionDescriptions;
        COMPLEX_OPTION_VALUES   complexOptionValues;
        std::vector<OPTIONSTR>  optionsSpecified;
    } OptionsDescriptorT;


// class Options

private:

    Options() {};
    Options(Options const&) = delete;
    Options& operator = (Options const&) = delete;

    static Options* m_Instance;


    // member variables

    GridfileT                   m_Gridfile = {"", ERROR::EMPTY_FILENAME};

    OptionsDescriptorT          m_CmdLine = {OPTIONS_ORIGIN::CMDLINE, {}, {}, {}, {}};
    OptionsDescriptorT          m_GridLine = {OPTIONS_ORIGIN::GRIDFILE, {}, {}, {}, {}};

    std::vector<OptionDetailsT> m_CmdLineOptionsDetails;                                                    // for Run_Details and YAML files
    std::map<std::string, std::string> m_optionDefaults;                                                    // for Run_Details and YAML files

    // member functions

    bool                        AddOptions(OptionValues *p_Options, po::options_description *p_OptionsDescription);
    std::vector<std::string>    AllowedOptionValues(const std::string p_OptionString);
    std::string                 AllowedOptionValuesFormatted(const std::string p_OptionString);
    int                         AdvanceOptionVariation(OptionsDescriptorT &p_OptionsDescriptor);

    void                        BuildDefaultsMap(po::options_description *p_OptionsDescription);

    std::tuple<std::string, int, std::vector<std::string>> ExpandShorthandOptionValues(int p_ArgCount, char *p_ArgStrings[]);

    bool                        IsSupportedNumericDataType(TYPENAME p_TypeName);

    ATTR                        OptionAttributes(const po::variables_map p_VM, const po::variables_map::const_iterator p_IT);

    PROGRAM_STATUS              ParseCommandLineOptions(int argc, char * argv[]);
    std::string                 ParseOptionValues(int p_ArgCount, char *p_ArgStrings[], OptionsDescriptorT &p_OptionsDescriptor);

    std::vector<OptionDetailsT> OptionDetails(const OptionsDescriptorT &p_Options);

public:

    static Options*             Instance();

    int                         AdvanceCmdLineOptionValues()  { return AdvanceOptionVariation(m_CmdLine); }
    int                         AdvanceGridLineOptionValues() { return AdvanceOptionVariation(m_GridLine); }

    int                         ApplyNextGridLine();

    std::string                 CheckDeprecatedOptionProperty(const std::string p_OptionProperty);
    std::string                 CheckDeprecatedOptionString(const std::string p_OptionString);
    std::string                 CheckDeprecatedOptionValue(const std::string p_OptionString, const std::string p_OptionValue);
    void                        CloseGridFile() { m_Gridfile.handle.close(); m_Gridfile.filename = ""; m_Gridfile.error = ERROR::EMPTY_FILENAME; }

    bool                        Initialise(int p_OptionCount, char *p_OptionStrings[]);
    bool                        InitialiseEvolvingObject(const std::string p_OptionsString);

    ERROR                       OpenGridFile(const std::string p_GridFilename);
    bool                        OptionDefaulted(const std::string p_OptionString) const                                 { return OPT_DEFAULTED(p_OptionString); }
    bool                        OptionSpecified(const std::string p_OptionString);

    COMPAS_VARIABLE             OptionValue(const T_ANY_PROPERTY p_Property) const;

    void                        PrintOptionHelp(const bool p_Verbose);

    ERROR                       RewindGridFile();

    ERROR                       SeekToGridFileLine(const unsigned int p_Line);

    std::string                 SetRandomSeed(const unsigned long int p_RandomSeed, const OPTIONS_ORIGIN p_OptionsSet);

    // getters

    ADD_OPTIONS_TO_SYSPARMS                     AddOptionsToSysParms() const                                            { return m_CmdLine.optionValues.m_AddOptionsToSysParms.type; }

    bool                                        AllowNonStrippedECSN() const                                            { return OPT_VALUE("allow-non-stripped-ECSN", m_AllowNonStrippedECSN, true); }
    bool                                        AllowMainSequenceStarToSurviveCommonEnvelope() const                    { return OPT_VALUE("common-envelope-allow-main-sequence-survive", m_AllowMainSequenceStarToSurviveCommonEnvelope, true); }
    bool                                        AllowRadiativeEnvelopeStarToSurviveCommonEnvelope() const               { return OPT_VALUE("common-envelope-allow-radiative-envelope-survive", m_AllowRadiativeEnvelopeStarToSurviveCommonEnvelope, true); }
    bool                                        AllowImmediateRLOFpostCEToSurviveCommonEnvelope() const                 { return OPT_VALUE("common-envelope-allow-immediate-rlof-post-ce-survive", m_AllowImmediateRLOFpostCEToSurviveCommonEnvelope, true); }
    bool                                        AllowRLOFAtBirth() const                                                { return OPT_VALUE("allow-rlof-at-birth", m_AllowRLOFAtBirth, true); }
    bool                                        AllowTouchingAtBirth() const                                            { return OPT_VALUE("allow-touching-at-birth", m_AllowTouchingAtBirth, true); }
    bool                                        AngularMomentumConservationDuringCircularisation() const                { return OPT_VALUE("angular-momentum-conservation-during-circularisation", m_AngularMomentumConservationDuringCircularisation, true); }


    BLACK_HOLE_KICKS_MODE                       BlackHoleKicksMode() const                                              { return OPT_VALUE("black-hole-kicks-mode", m_BlackHoleKicksMode.type, true); }
    
    CASE_BB_STABILITY_PRESCRIPTION              CaseBBStabilityPrescription() const                                     { return OPT_VALUE("case-BB-stability-prescription", m_CaseBBStabilityPrescription.type, true); }
    
    bool                                        CheckPhotonTiringLimit() const                                          { return OPT_VALUE("check-photon-tiring-limit", m_CheckPhotonTiringLimit, true); }

    CHE_MODE                                    CHEMode() const                                                         { return OPT_VALUE("chemically-homogeneous-evolution-mode", m_CheMode.type, true); }

    bool                                        CirculariseBinaryDuringMassTransfer() const                             { return OPT_VALUE("circularise-binary-during-mass-transfer", m_CirculariseBinaryDuringMassTransfer, true); }

    std::vector<OptionDetailsT>                 CmdLineOptionsDetails() const                                           { return m_CmdLineOptionsDetails; }

    bool                                        CommandLineGrid() const                                                 { return m_CmdLine.complexOptionValues.size() != 0; }
    
    double                                      CommonEnvelopeAlpha() const                                             { return OPT_VALUE("common-envelope-alpha", m_CommonEnvelopeAlpha, true); }
    double                                      CommonEnvelopeAlphaThermal() const                                      { return OPT_VALUE("common-envelope-alpha-thermal", m_CommonEnvelopeAlphaThermal, true); }
    CE_FORMALISM                                CommonEnvelopeFormalism() const                                         { return OPT_VALUE("common-envelope-formalism", m_CommonEnvelopeFormalism.type, true); }
    double                                      CommonEnvelopeLambda() const                                            { return OPT_VALUE("common-envelope-lambda", m_CommonEnvelopeLambda, true); }
    double                                      CommonEnvelopeLambdaMultiplier() const                                  { return OPT_VALUE("common-envelope-lambda-multiplier", m_CommonEnvelopeLambdaMultiplier, true); }
    bool                                        CommonEnvelopeLambdaNanjingEnhanced() const                             { return OPT_VALUE("common-envelope-lambda-nanjing-enhanced", m_CommonEnvelopeLambdaNanjingEnhanced, true); }
    bool                                        CommonEnvelopeLambdaNanjingInterpolateInMass() const                    { return OPT_VALUE("common-envelope-lambda-nanjing-interpolate-in-mass", m_CommonEnvelopeLambdaNanjingInterpolateInMass, true); }
    bool                                        CommonEnvelopeLambdaNanjingInterpolateInMetallicity() const             { return OPT_VALUE("common-envelope-lambda-nanjing-interpolate-in-metallicity", m_CommonEnvelopeLambdaNanjingInterpolateInMetallicity, true); }
    bool                                        CommonEnvelopeLambdaNanjingUseRejuvenatedMass() const                   { return OPT_VALUE("common-envelope-lambda-nanjing-use-rejuvenated-mass", m_CommonEnvelopeLambdaNanjingUseRejuvenatedMass, true); }
    CE_LAMBDA_PRESCRIPTION                      CommonEnvelopeLambdaPrescription() const                                { return OPT_VALUE("common-envelope-lambda-prescription", m_CommonEnvelopeLambdaPrescription.type, true); }
    double                                      CommonEnvelopeMassAccretionConstant() const                             { return OPT_VALUE("common-envelope-mass-accretion-constant", m_CommonEnvelopeMassAccretionConstant, true); }
    double                                      CommonEnvelopeMassAccretionMax() const                                  { return OPT_VALUE("common-envelope-mass-accretion-max", m_CommonEnvelopeMassAccretionMax, true); }
    double                                      CommonEnvelopeMassAccretionMin() const                                  { return OPT_VALUE("common-envelope-mass-accretion-min", m_CommonEnvelopeMassAccretionMin, true); }
    CE_ACCRETION_PRESCRIPTION                   CommonEnvelopeMassAccretionPrescription() const                         { return OPT_VALUE("common-envelope-mass-accretion-prescription", m_CommonEnvelopeMassAccretionPrescription.type, true); }
    double                                      CommonEnvelopeRecombinationEnergyDensity() const                        { return OPT_VALUE("common-envelope-recombination-energy-density", m_CommonEnvelopeRecombinationEnergyDensity, true); }
    double                                      CommonEnvelopeSlopeKruckow() const                                      { return OPT_VALUE("common-envelope-slope-kruckow", m_CommonEnvelopeSlopeKruckow, true); }

    double                                      ConvectiveEnvelopeTemperatureThreshold() const                          { return OPT_VALUE("convective-envelope-temperature-threshold", m_ConvectiveEnvelopeTemperatureThreshold, true); }

    double                                      CoolWindMassLossMultiplier() const                                      { return OPT_VALUE("cool-wind-mass-loss-multiplier", m_CoolWindMassLossMultiplier, true); }

    std::vector<std::string>                    DebugClasses() const                                                    { return m_CmdLine.optionValues.m_DebugClasses; }
    int                                         DebugLevel() const                                                      { return m_CmdLine.optionValues.m_DebugLevel; }
    bool                                        DebugToFile() const                                                     { return m_CmdLine.optionValues.m_DebugToFile; }
    bool                                        DetailedOutput() const                                                  { return m_CmdLine.optionValues.m_DetailedOutput; }

    bool                                        EnableRotationallyEnhancedMassLoss() const                              { return OPT_VALUE("enable-rotationally-enhanced-mass-loss", m_EnableRotationallyEnhancedMassLoss, true); }
    bool                                        EnhanceCHELifetimesLuminosities() const                                 { return OPT_VALUE("enhance-CHE-lifetimes-luminosities", m_EnhanceCHELifetimesLuminosities, false); }
    bool                                        EnableWarnings() const                                                  { return m_CmdLine.optionValues.m_EnableWarnings; }
    bool                                        ErrorsToFile() const                                                    { return m_CmdLine.optionValues.m_ErrorsToFile; }
    
    double                                      Eccentricity() const                                                    { return OPT_VALUE("eccentricity", m_Eccentricity, true); }
    ECCENTRICITY_DISTRIBUTION                   EccentricityDistribution() const                                        { return OPT_VALUE("eccentricity-distribution", m_EccentricityDistribution.type, true); }
    double                                      EccentricityDistributionMax() const                                     { return OPT_VALUE("eccentricity-distribution-max", m_EccentricityDistributionMax, true); }
    double                                      EccentricityDistributionMin() const                                     { return OPT_VALUE("eccentricity-distribution-min", m_EccentricityDistributionMin, true); }
    double                                      EddingtonAccretionFactor() const                                        { return OPT_VALUE("eddington-accretion-factor", m_EddingtonAccretionFactor, true); }
    ENVELOPE_STATE_PRESCRIPTION                 EnvelopeStatePrescription() const                                       { return OPT_VALUE("envelope-state-prescription", m_EnvelopeStatePrescription.type, true); }
    EVOLUTION_MODE                              EvolutionMode() const                                                   { return m_CmdLine.optionValues.m_EvolutionMode.type; }
    bool                                        EvolveDoubleWhiteDwarfs() const                                         { return OPT_VALUE("evolve-double-white-dwarfs", m_EvolveDoubleWhiteDwarfs, true); }
    bool                                        EvolveMainSequenceMergers() const                                       { return OPT_VALUE("evolve-main-sequence-mergers", m_EvolveMainSequenceMergers, true); }
    bool                                        EvolvePulsars() const                                                   { return OPT_VALUE("evolve-pulsars", m_EvolvePulsars, true); }
    bool                                        EvolveUnboundSystems() const                                            { return OPT_VALUE("evolve-unbound-systems", m_EvolveUnboundSystems, true); }
    bool                                        EmitGravitationalRadiation() const                                      { return OPT_VALUE("emit-gravitational-radiation", m_EmitGravitationalRadiation, true); }
    bool                                        ExpelConvectiveEnvelopeAboveLuminosityThreshold() const                 { return OPT_VALUE("expel-convective-envelope-above-luminosity-threshold", m_ExpelConvectiveEnvelopeAboveLuminosityThreshold, true); }

    bool                                        FixedRandomSeedCmdLine() const                                          { return m_CmdLine.optionValues.m_FixedRandomSeed; }
    bool                                        FixedRandomSeedGridLine() const                                         { return m_GridLine.optionValues.m_FixedRandomSeed; }
    double                                      FixedUK() const                                                         { return m_GridLine.optionValues.m_UseFixedUK || m_CmdLine.optionValues.m_FixedUK; }
    FP_ERROR_MODE                               FPErrorMode() const                                                     { return m_CmdLine.optionValues.m_FPErrorMode.type; }
    SN_ENGINE                                   FryerSupernovaEngine() const                                            { return OPT_VALUE("fryer-supernova-engine", m_FryerSupernovaEngine.type, true); }
    double                                      Fryer22fmix() const                                                     { return OPT_VALUE("fryer-22-fmix", m_Fryer22fmix, true); }
    double                                      Fryer22Mcrit() const                                                    { return OPT_VALUE("fryer-22-mcrit", m_Fryer22Mcrit, true); }

    std::string                                 GridFilename() const                                                    { return m_CmdLine.optionValues.m_GridFilename; }
    std::streamsize                             GridStartLine() const                                                   { return m_CmdLine.optionValues.m_GridStartLine; }
    std::streamsize                             GridLinesToProcess() const                                              { return m_CmdLine.optionValues.m_GridLinesToProcess; }

    size_t                                      HDF5ChunkSize() const                                                   { return m_CmdLine.optionValues.m_HDF5ChunkSize; }
    size_t                                      HDF5BufferSize() const                                                  { return m_CmdLine.optionValues.m_HDF5BufferSize; }
    bool                                        HMXRBinaries() const                                                    { return OPT_VALUE("hmxr-binaries", m_HMXRBinaries, true); }

    bool                                        IncludeWDBinariesAsDCO() const                                          { return OPT_VALUE("include-WD-binaries-as-DCO", m_WDBinariesAsDCO, true); }
    double                                      InitialMass() const                                                     { return OPT_VALUE("initial-mass", m_InitialMass, true); }
    double                                      InitialMass1() const                                                    { return OPT_VALUE("initial-mass-1", m_InitialMass1, true); }
    double                                      InitialMass2() const                                                    { return OPT_VALUE("initial-mass-2", m_InitialMass2, true); }

    INITIAL_MASS_FUNCTION                       InitialMassFunction() const                                             { return OPT_VALUE("initial-mass-function", m_InitialMassFunction.type, true); }
    double                                      InitialMassFunctionMax() const                                          { return OPT_VALUE("initial-mass-max", m_InitialMassFunctionMax, true); }
    double                                      InitialMassFunctionMin() const                                          { return OPT_VALUE("initial-mass-min", m_InitialMassFunctionMin, true); }
    double                                      InitialMassFunctionPower() const                                        { return OPT_VALUE("initial-mass-power", m_InitialMassFunctionPower, true); }

    KICK_DIRECTION_DISTRIBUTION                 KickDirectionDistribution() const                                       { return OPT_VALUE("kick-direction-distribution", m_KickDirectionDistribution.type, true); }
    double                                      KickDirectionPower() const                                              { return OPT_VALUE("kick-direction-power", m_KickDirectionPower, true); }
    double                                      KickScalingFactor() const                                               { return OPT_VALUE("kick-scaling-factor", m_KickScalingFactor, true); }
    KICK_MAGNITUDE_DISTRIBUTION                 KickMagnitudeDistribution() const                                       { return OPT_VALUE("kick-magnitude-distribution", m_KickMagnitudeDistribution.type, true); }

    double                                      KickMagnitudeDistributionMaximum() const                                { return OPT_VALUE("kick-magnitude-max", m_KickMagnitudeDistributionMaximum, true); }

    double                                      KickMagnitudeDistributionSigmaCCSN_BH() const                           { return OPT_VALUE("kick-magnitude-sigma-CCSN-BH", m_KickMagnitudeDistributionSigmaCCSN_BH, true); }
    double                                      KickMagnitudeDistributionSigmaCCSN_NS() const                           { return OPT_VALUE("kick-magnitude-sigma-CCSN-NS", m_KickMagnitudeDistributionSigmaCCSN_NS, true); }
    double                                      KickMagnitudeDistributionSigmaForECSN() const                           { return OPT_VALUE("kick-magnitude-sigma-ECSN", m_KickMagnitudeDistributionSigmaForECSN, true); }
    double                                      KickMagnitudeDistributionSigmaForUSSN() const                           { return OPT_VALUE("kick-magnitude-sigma-USSN", m_KickMagnitudeDistributionSigmaForUSSN, true); }

    double                                      KickMagnitude() const                                                   { return OPT_VALUE("kick-magnitude", m_KickMagnitude, true); }
    double                                      KickMagnitude1() const                                                  { return OPT_VALUE("kick-magnitude-1", m_KickMagnitude1, true); }
    double                                      KickMagnitude2() const                                                  { return OPT_VALUE("kick-magnitude-2", m_KickMagnitude2, true); }

    double                                      KickMagnitudeRandom() const                                             { return OPT_VALUE("kick-magnitude-random", m_KickMagnitudeRandom, true); }
    double                                      KickMagnitudeRandom1() const                                            { return OPT_VALUE("kick-magnitude-random-1", m_KickMagnitudeRandom1, true); }
    double                                      KickMagnitudeRandom2() const                                            { return OPT_VALUE("kick-magnitude-random-2", m_KickMagnitudeRandom2, true); }

    std::vector<std::string>                    LogClasses() const                                                      { return m_CmdLine.optionValues.m_LogClasses; }
    std::string                                 LogfileCommonEnvelopes() const                                          { return m_CmdLine.optionValues.m_LogfileCommonEnvelopes; }
    int                                         LogfileCommonEnvelopesRecordTypes() const                               { return m_CmdLine.optionValues.m_LogfileCommonEnvelopesRecordTypes; }
    std::string                                 LogfileDefinitionsFilename() const                                      { return m_CmdLine.optionValues.m_LogfileDefinitionsFilename; }
    std::string                                 LogfileDetailedOutput() const                                           { return m_CmdLine.optionValues.m_Populated && !m_CmdLine.optionValues.m_VM["logfile-detailed-output"].defaulted()
                                                                                                                                    ? m_CmdLine.optionValues.m_LogfileDetailedOutput
                                                                                                                                    : (m_CmdLine.optionValues.m_EvolutionMode.type == EVOLUTION_MODE::SSE
                                                                                                                                        ? std::get<0>(LOGFILE_DESCRIPTOR.at(LOGFILE::SSE_DETAILED_OUTPUT))
                                                                                                                                        : std::get<0>(LOGFILE_DESCRIPTOR.at(LOGFILE::BSE_DETAILED_OUTPUT))
                                                                                                                                      );
                                                                                                                        }
    int                                         LogfileDetailedOutputRecordTypes() const                                { return m_CmdLine.optionValues.m_LogfileDetailedOutputRecordTypes; }
    std::string                                 LogfileDoubleCompactObjects() const                                     { return m_CmdLine.optionValues.m_LogfileDoubleCompactObjects; }
    int                                         LogfileDoubleCompactObjectsRecordTypes() const                          { return m_CmdLine.optionValues.m_LogfileDoubleCompactObjectsRecordTypes; }
    std::string                                 LogfileNamePrefix() const                                               { return m_CmdLine.optionValues.m_LogfileNamePrefix; }
    std::string                                 LogfilePulsarEvolution() const                                          { return m_CmdLine.optionValues.m_Populated && !m_CmdLine.optionValues.m_VM["logfile-pulsar-evolution"].defaulted()
                                                                                                                                    ? m_CmdLine.optionValues.m_LogfilePulsarEvolution
                                                                                                                                    : (m_CmdLine.optionValues.m_EvolutionMode.type == EVOLUTION_MODE::SSE
                                                                                                                                        ? std::get<0>(LOGFILE_DESCRIPTOR.at(LOGFILE::SSE_PULSAR_EVOLUTION))
                                                                                                                                        : std::get<0>(LOGFILE_DESCRIPTOR.at(LOGFILE::BSE_PULSAR_EVOLUTION))
                                                                                                                                      );
                                                                                                                        }
    int                                         LogfilePulsarEvolutionRecordTypes() const                               { return m_CmdLine.optionValues.m_LogfilePulsarEvolutionRecordTypes; }
    std::string                                 LogfileRLOFParameters() const                                           { return m_CmdLine.optionValues.m_LogfileRLOFParameters; }
    int                                         LogfileRLOFParametersRecordTypes() const                                { return m_CmdLine.optionValues.m_LogfileRLOFParametersRecordTypes; }
    std::string                                 LogfileSupernovae() const                                               { return m_CmdLine.optionValues.m_Populated && !m_CmdLine.optionValues.m_VM["logfile-supernovae"].defaulted()
                                                                                                                                    ? m_CmdLine.optionValues.m_LogfileSupernovae
                                                                                                                                    : (m_CmdLine.optionValues.m_EvolutionMode.type == EVOLUTION_MODE::SSE
                                                                                                                                        ? std::get<0>(LOGFILE_DESCRIPTOR.at(LOGFILE::SSE_SUPERNOVAE))
                                                                                                                                        : std::get<0>(LOGFILE_DESCRIPTOR.at(LOGFILE::BSE_SUPERNOVAE))
                                                                                                                                      );
                                                                                                                        }
    int                                         LogfileSupernovaeRecordTypes() const                                    { return m_CmdLine.optionValues.m_LogfileSupernovaeRecordTypes; }
    std::string                                 LogfileSwitchLog() const                                                { return m_CmdLine.optionValues.m_Populated && !m_CmdLine.optionValues.m_VM["logfile-switch-log"].defaulted()
                                                                                                                                    ? m_CmdLine.optionValues.m_LogfileSwitchLog
                                                                                                                                    : (m_CmdLine.optionValues.m_EvolutionMode.type == EVOLUTION_MODE::SSE
                                                                                                                                        ? std::get<0>(LOGFILE_DESCRIPTOR.at(LOGFILE::SSE_SWITCH_LOG))
                                                                                                                                        : std::get<0>(LOGFILE_DESCRIPTOR.at(LOGFILE::BSE_SWITCH_LOG))
                                                                                                                                      );
                                                                                                                        }
    std::string                                 LogfileSystemParameters() const                                         { return m_CmdLine.optionValues.m_Populated && !m_CmdLine.optionValues.m_VM["logfile-system-parameters"].defaulted()
                                                                                                                                    ? m_CmdLine.optionValues.m_LogfileSystemParameters
                                                                                                                                    : (m_CmdLine.optionValues.m_EvolutionMode.type == EVOLUTION_MODE::SSE
                                                                                                                                        ? std::get<0>(LOGFILE_DESCRIPTOR.at(LOGFILE::SSE_SYSTEM_PARAMETERS))
                                                                                                                                        : std::get<0>(LOGFILE_DESCRIPTOR.at(LOGFILE::BSE_SYSTEM_PARAMETERS))
                                                                                                                                      );
                                                                                                                        }
    int                                         LogfileSystemParametersRecordTypes() const                              { return m_CmdLine.optionValues.m_LogfileSystemParametersRecordTypes; }
    LOGFILETYPE                                 LogfileType() const                                                     { return m_CmdLine.optionValues.m_LogfileType.type; }
    std::string                                 LogfileTypeString() const                                               { return m_CmdLine.optionValues.m_LogfileType.typeString; }
    int                                         LogLevel() const                                                        { return m_CmdLine.optionValues.m_LogLevel; }
    
    double                                      LuminosityToMassThreshold() const                                       { return OPT_VALUE("luminosity-to-mass-threshold", m_LuminosityToMassThreshold, true); }

    double                                      LuminousBlueVariableFactor() const                                      { return OPT_VALUE("luminous-blue-variable-multiplier", m_LuminousBlueVariableFactor, true); }
    LBV_MASS_LOSS_PRESCRIPTION                  LBVMassLossPrescription() const                                         { return OPT_VALUE("LBV-mass-loss-prescription", m_LBVMassLossPrescription.type, true); }
    
    CORE_MASS_PRESCRIPTION                      MainSequenceCoreMassPrescription() const                                { return OPT_VALUE("main-sequence-core-mass-prescription", m_MainSequenceCoreMassPrescription.type, true); }
    
    double                                      MassChangeFraction() const                                              { return m_CmdLine.optionValues.m_MassChangeFraction; }
    
    MASS_LOSS_PRESCRIPTION                      MassLossPrescription() const                                            { return OPT_VALUE("mass-loss-prescription", m_MassLossPrescription.type, true); }

    double                                      MassRatio() const                                                       { return OPT_VALUE("mass-ratio", m_MassRatio, true); }
    MASS_RATIO_DISTRIBUTION                     MassRatioDistribution() const                                           { return OPT_VALUE("mass-ratio-distribution", m_MassRatioDistribution.type, true); }
    double                                      MassRatioDistributionMax() const                                        { return OPT_VALUE("mass-ratio-max", m_MassRatioDistributionMax, true); }
    double                                      MassRatioDistributionMin() const                                        { return OPT_VALUE("mass-ratio-min", m_MassRatioDistributionMin, true); }

    MT_ACCRETION_EFFICIENCY_PRESCRIPTION        MassTransferAccretionEfficiencyPrescription() const                     { return OPT_VALUE("mass-transfer-accretion-efficiency-prescription", m_MassTransferAccretionEfficiencyPrescription.type, true); }
    MT_ANGULAR_MOMENTUM_LOSS_PRESCRIPTION       MassTransferAngularMomentumLossPrescription() const                     { return OPT_VALUE("mass-transfer-angular-momentum-loss-prescription", m_MassTransferAngularMomentumLossPrescription.type, true); }
    double                                      MassTransferCParameter() const                                          { return OPT_VALUE("mass-transfer-thermal-limit-C", m_MassTransferCParameter, true); }

    double                                      MassTransferCriticalMassRatioMSLowMassDegenerateAccretor() const        { return OPT_VALUE("critical-mass-ratio-ms-low-mass-degenerate-accretor", m_MassTransferCriticalMassRatioMSLowMassDegenerateAccretor, true); }
    double                                      MassTransferCriticalMassRatioMSLowMassNonDegenerateAccretor() const     { return OPT_VALUE("critical-mass-ratio-ms-low-mass-non-degenerate-accretor", m_MassTransferCriticalMassRatioMSLowMassNonDegenerateAccretor, true); }
    double                                      MassTransferCriticalMassRatioMSHighMassDegenerateAccretor() const       { return OPT_VALUE("critical-mass-ratio-ms-high-mass-degenerate-accretor", m_MassTransferCriticalMassRatioMSHighMassDegenerateAccretor, true); }
    double                                      MassTransferCriticalMassRatioMSHighMassNonDegenerateAccretor() const    { return OPT_VALUE("critical-mass-ratio-ms-high-mass-non-degenerate-accretor", m_MassTransferCriticalMassRatioMSHighMassNonDegenerateAccretor, true); }
    double                                      MassTransferCriticalMassRatioGiantDegenerateAccretor() const            { return OPT_VALUE("critical-mass-ratio-giant-degenerate-accretor", m_MassTransferCriticalMassRatioGiantDegenerateAccretor, true); }
    double                                      MassTransferCriticalMassRatioGiantNonDegenerateAccretor() const         { return OPT_VALUE("critical-mass-ratio-giant-non-degenerate-accretor", m_MassTransferCriticalMassRatioGiantNonDegenerateAccretor, true); }
    double                                      MassTransferCriticalMassRatioHGDegenerateAccretor() const               { return OPT_VALUE("critical-mass-ratio-hg-degenerate-accretor", m_MassTransferCriticalMassRatioHGDegenerateAccretor, true); }
    double                                      MassTransferCriticalMassRatioHGNonDegenerateAccretor() const            { return OPT_VALUE("critical-mass-ratio-hg-non-degenerate-accretor", m_MassTransferCriticalMassRatioHGNonDegenerateAccretor, true); }
    double                                      MassTransferCriticalMassRatioHeliumGiantDegenerateAccretor() const      { return OPT_VALUE("critical-mass-ratio-helium-giant-degenerate-accretor", m_MassTransferCriticalMassRatioHeliumGiantDegenerateAccretor, true); }
    double                                      MassTransferCriticalMassRatioHeliumGiantNonDegenerateAccretor() const   { return OPT_VALUE("critical-mass-ratio-helium-giant-non-degenerate-accretor", m_MassTransferCriticalMassRatioHeliumGiantNonDegenerateAccretor, true); }
    double                                      MassTransferCriticalMassRatioHeliumHGDegenerateAccretor() const         { return OPT_VALUE("critical-mass-ratio-helium-hg-degenerate-accretor", m_MassTransferCriticalMassRatioHeliumHGDegenerateAccretor, true); }
    double                                      MassTransferCriticalMassRatioHeliumHGNonDegenerateAccretor() const      { return OPT_VALUE("critical-mass-ratio-helium-hg-non-degenerate-accretor", m_MassTransferCriticalMassRatioHeliumHGNonDegenerateAccretor, true); }
    double                                      MassTransferCriticalMassRatioHeliumMSDegenerateAccretor() const         { return OPT_VALUE("critical-mass-ratio-helium-ms-degenerate-accretor", m_MassTransferCriticalMassRatioHeliumMSDegenerateAccretor, true); }
    double                                      MassTransferCriticalMassRatioHeliumMSNonDegenerateAccretor() const      { return OPT_VALUE("critical-mass-ratio-helium-ms-non-degenerate-accretor", m_MassTransferCriticalMassRatioHeliumMSNonDegenerateAccretor, true); }
    double                                      MassTransferCriticalMassRatioWhiteDwarfDegenerateAccretor() const       { return OPT_VALUE("critical-mass-ratio-white-dwarf-degenerate-accretor", m_MassTransferCriticalMassRatioWhiteDwarfDegenerateAccretor, true); }
    double                                      MassTransferCriticalMassRatioWhiteDwarfNonDegenerateAccretor() const    { return OPT_VALUE("critical-mass-ratio-white-dwarf-non-degenerate-accretor", m_MassTransferCriticalMassRatioWhiteDwarfNonDegenerateAccretor, true); }

    double                                      MassTransferFractionAccreted() const                                    { return OPT_VALUE("mass-transfer-fa", m_MassTransferFractionAccreted, true); }
    double                                      MassTransferJloss() const                                               { return OPT_VALUE("mass-transfer-jloss", m_MassTransferJloss, true); }
    double                                      MassTransferJlossMacLeodLinearFractionDegen() const                     { return OPT_VALUE("mass-transfer-jloss-macleod-linear-fraction-degen", m_MassTransferJlossMacLeodLinearFractionDegen, true); }
    double                                      MassTransferJlossMacLeodLinearFractionNonDegen() const                  { return OPT_VALUE("mass-transfer-jloss-macleod-linear-fraction-non-degen", m_MassTransferJlossMacLeodLinearFractionNonDegen, true); }
    MT_REJUVENATION_PRESCRIPTION                MassTransferRejuvenationPrescription() const                            { return OPT_VALUE("mass-transfer-rejuvenation-prescription", m_MassTransferRejuvenationPrescription.type, true); }
    MT_THERMALLY_LIMITED_VARIATION              MassTransferThermallyLimitedVariation() const                           { return OPT_VALUE("mass-transfer-thermal-limit-accretor-multiplier", m_MassTransferThermallyLimitedVariation.type, true); }
    double                                      MaxEvolutionTime() const                                                { return OPT_VALUE("maximum-evolution-time", m_MaxEvolutionTime, true); }
    double                                      MaximumNeutronStarMass() const                                          { return OPT_VALUE("maximum-neutron-star-mass", m_MaximumNeutronStarMass, true); }
    unsigned long int                           MaxNumberOfTimestepIterations() const                                   { return OPT_VALUE("maximum-number-timestep-iterations", m_MaxNumberOfTimestepIterations, true); }
    double                                      MaximumDonorMass() const                                                { return OPT_VALUE("maximum-mass-donor-nandez-ivanova", m_MaximumMassDonorNandezIvanova, true); }
    double                                      MCBUR1() const                                                          { return OPT_VALUE("mcbur1", m_mCBUR1, true); }

    double                                      Metallicity() const                                                     { return OPT_VALUE("metallicity", m_Metallicity, true); }
    METALLICITY_DISTRIBUTION                    MetallicityDistribution() const                                         { return OPT_VALUE("metallicity-distribution", m_MetallicityDistribution.type, true); }
    double                                      MetallicityDistributionMax() const                                      { return OPT_VALUE("metallicity-distribution-max", m_MetallicityDistributionMax, true); }
    double                                      MetallicityDistributionMin() const                                      { return OPT_VALUE("metallicity-distribution-min", m_MetallicityDistributionMin, true); }

    double                                      MinimumMassSecondary() const                                            { return OPT_VALUE("minimum-secondary-mass", m_MinimumMassSecondary, true); }

    double                                      MullerMandelKickMultiplierBH() const                                    { return OPT_VALUE("muller-mandel-kick-multiplier-BH", m_MullerMandelKickBH, true); }
    double                                      MullerMandelKickMultiplierNS() const                                    { return OPT_VALUE("muller-mandel-kick-multiplier-NS", m_MullerMandelKickNS, true); }
    double                                      MullerMandelSigmaKick() const                                           { return OPT_VALUE("muller-mandel-sigma-kick", m_MullerMandelSigmaKick, true); }

    bool                                        NatalKickForPPISN() const                                               { return OPT_VALUE("natal-kick-for-PPISN", m_NatalKickForPPISN, false); }
    NEUTRINO_MASS_LOSS_PRESCRIPTION             NeutrinoMassLossAssumptionBH() const                                    { return OPT_VALUE("neutrino-mass-loss-BH-formation", m_NeutrinoMassLossAssumptionBH.type, true); }
    double                                      NeutrinoMassLossValueBH() const                                         { return OPT_VALUE("neutrino-mass-loss-BH-formation-value", m_NeutrinoMassLossValueBH, true); }

    NS_EOS                                      NeutronStarEquationOfState() const                                      { return OPT_VALUE("neutron-star-equation-of-state", m_NeutronStarEquationOfState.type, true); }

    std::string                                 Notes(const size_t p_Idx) const                                         { return OPT_VALUE("notes", m_Notes[p_Idx], true); }
    std::vector<std::string>                    Notes() const                                                           { return OPT_VALUE("notes", m_Notes, true); }
    std::string                                 NotesHdrs(const size_t p_Idx) const                                     { return m_CmdLine.optionValues.m_NotesHdrs[p_Idx]; }
    std::vector<std::string>                    NotesHdrs() const                                                       { return m_CmdLine.optionValues.m_NotesHdrs; }
 
    size_t                                      nObjectsToEvolve() const                                                { return m_CmdLine.optionValues.m_ObjectsToEvolve; }
    OB_MASS_LOSS_PRESCRIPTION                   OBMassLossPrescription() const                                          { return OPT_VALUE("OB-mass-loss-prescription", m_OBMassLossPrescription.type, true); }
    bool                                        OptimisticCHE() const                                                   { return CHEMode() == CHE_MODE::OPTIMISTIC; }

    double                                      OrbitalPeriod() const                                                   { return OPT_VALUE("orbital-period", m_OrbitalPeriod, true); }
    ORBITAL_PERIOD_DISTRIBUTION                 OrbitalPeriodDistribution() const                                       { return OPT_VALUE("orbital-period-distribution", m_OrbitalPeriodDistribution.type, true); }
    double                                      OrbitalPeriodDistributionMax() const                                    { return OPT_VALUE("orbital-period-max", m_OrbitalPeriodDistributionMax, true); }
    double                                      OrbitalPeriodDistributionMin() const                                    { return OPT_VALUE("orbital-period-min", m_OrbitalPeriodDistributionMin, true); }

    std::string                                 OutputContainerName() const                                             { return m_CmdLine.optionValues.m_OutputContainerName; }
    std::string                                 OutputPathString() const                                                { return m_CmdLine.optionValues.m_OutputPathString; }

    double                                      OverallWindMassLossMultiplier() const                                   { return OPT_VALUE("overall-wind-mass-loss-multiplier", m_OverallWindMassLossMultiplier, true); }

    double                                      PairInstabilityLowerLimit() const                                       { return OPT_VALUE("PISN-lower-limit", m_PairInstabilityLowerLimit, true); }
    double                                      PairInstabilityUpperLimit() const                                       { return OPT_VALUE("PISN-upper-limit", m_PairInstabilityUpperLimit, true); }

    bool                                        PopulationDataPrinting() const                                          { return m_CmdLine.optionValues.m_PopulationDataPrinting; }
    bool                                        PrintBoolAsString() const                                               { return m_CmdLine.optionValues.m_PrintBoolAsString; }

    PULSAR_BIRTH_MAGNETIC_FIELD_DISTRIBUTION    PulsarBirthMagneticFieldDistribution() const                            { return OPT_VALUE("pulsar-birth-magnetic-field-distribution", m_PulsarBirthMagneticFieldDistribution.type, true); }
    double                                      PulsarBirthMagneticFieldDistributionMax() const                         { return OPT_VALUE("pulsar-birth-magnetic-field-distribution-max", m_PulsarBirthMagneticFieldDistributionMax, true); }
    double                                      PulsarBirthMagneticFieldDistributionMin() const                         { return OPT_VALUE("pulsar-birth-magnetic-field-distribution-min", m_PulsarBirthMagneticFieldDistributionMin, true); }

    PULSAR_BIRTH_SPIN_PERIOD_DISTRIBUTION       PulsarBirthSpinPeriodDistribution() const                               { return OPT_VALUE("pulsar-birth-spin-period-distribution", m_PulsarBirthSpinPeriodDistribution.type, true); }
    double                                      PulsarBirthSpinPeriodDistributionMax() const                            { return OPT_VALUE("pulsar-birth-spin-period-distribution-max", m_PulsarBirthSpinPeriodDistributionMax, true); }
    double                                      PulsarBirthSpinPeriodDistributionMin() const                            { return OPT_VALUE("pulsar-birth-spin-period-distribution-min", m_PulsarBirthSpinPeriodDistributionMin, true); }

    double                                      PulsarLog10MinimumMagneticField() const                                 { return OPT_VALUE("pulsar-minimum-magnetic-field", m_PulsarLog10MinimumMagneticField, true); }

    double                                      PulsarMagneticFieldDecayMassscale() const                               { return OPT_VALUE("pulsar-magnetic-field-decay-massscale", m_PulsarMagneticFieldDecayMassscale, true); }
    double                                      PulsarMagneticFieldDecayTimescale() const                               { return OPT_VALUE("pulsar-magnetic-field-decay-timescale", m_PulsarMagneticFieldDecayTimescale, true); }

    PPI_PRESCRIPTION                            PulsationalPairInstabilityPrescription() const                          { return OPT_VALUE("pulsational-pair-instability-prescription", m_PulsationalPairInstabilityPrescription.type, true); }
    double                                      PulsationalPairInstabilityLowerLimit() const                            { return OPT_VALUE("PPI-lower-limit", m_PulsationalPairInstabilityLowerLimit, true); }
    double                                      PulsationalPairInstabilityUpperLimit() const                            { return OPT_VALUE("PPI-upper-limit", m_PulsationalPairInstabilityUpperLimit, true); }
    double                                      PulsationalPairInstabilityCOCoreShiftHendriks() const                   { return OPT_VALUE("PPI-CO-Core-Shift-Hendriks", m_PulsationalPairInstabilityCOCoreShiftHendriks, true); }
    
    QCRIT_PRESCRIPTION                          QCritPrescription() const                                               { return OPT_VALUE("critical-mass-ratio-prescription", m_QCritPrescription.type, true); }

    bool                                        Quiet() const                                                           { return m_CmdLine.optionValues.m_Quiet; }

    double                                      RadialChangeFraction() const                                            { return m_CmdLine.optionValues.m_RadialChangeFraction; }
    
    unsigned long int                           RandomSeed() const                                                      { return OPT_VALUE("random-seed", m_RandomSeed, true); }
    unsigned long int                           RandomSeedCmdLine() const                                               { return m_CmdLine.optionValues.m_RandomSeed; }
    unsigned long int                           RandomSeedGridLine() const                                              { return m_GridLine.optionValues.m_RandomSeed; }

    REMNANT_MASS_PRESCRIPTION                   RemnantMassPrescription() const                                         { return OPT_VALUE("remnant-mass-prescription", m_RemnantMassPrescription.type, true); }
    
    bool                                        RequestedHelp() const                                                   { return m_CmdLine.optionValues.m_VM["help"].as<bool>(); }
    bool                                        RequestedVersion() const                                                { return m_CmdLine.optionValues.m_VM["version"].as<bool>(); }
    
    bool                                        RetainCoreMassDuringCaseAMassTransfer() const                           { return m_CmdLine.optionValues.m_RetainCoreMassDuringCaseAMassTransfer; }
    
    bool                                        RLOFPrinting() const                                                    { return m_CmdLine.optionValues.m_RlofPrinting; }

    double                                      RocketKickMagnitude1() const                                            { return OPT_VALUE("rocket-kick-magnitude-1", m_RocketKickMagnitude1, true); }
    double                                      RocketKickMagnitude2() const                                            { return OPT_VALUE("rocket-kick-magnitude-2", m_RocketKickMagnitude2, true); }
    double                                      RocketKick_Phi1() const                                                 { return OPT_VALUE("rocket-kick-phi-1", m_RocketKickPhi1, true); }
    double                                      RocketKick_Phi2() const                                                 { return OPT_VALUE("rocket-kick-phi-2", m_RocketKickPhi2, true); }
    double                                      RocketKick_Theta1() const                                               { return OPT_VALUE("rocket-kick-theta-1", m_RocketKickTheta1, true); }
    double                                      RocketKick_Theta2() const                                               { return OPT_VALUE("rocket-kick-theta-2", m_RocketKickTheta2, true); }

    ROTATIONAL_VELOCITY_DISTRIBUTION            RotationalVelocityDistribution() const                                  { return OPT_VALUE("rotational-velocity-distribution", m_RotationalVelocityDistribution.type, true); }
    double                                      RotationalFrequency() const                                             { return OPT_VALUE("rotational-frequency", m_RotationalFrequency, true); }
    double                                      RotationalFrequency1() const                                            { return OPT_VALUE("rotational-frequency-1", m_RotationalFrequency1, true); }
    double                                      RotationalFrequency2() const                                            { return OPT_VALUE("rotational-frequency-2", m_RotationalFrequency2, true); }
    RSG_MASS_LOSS_PRESCRIPTION                  RSGMassLossPrescription() const                                         { return OPT_VALUE("RSG-mass-loss-prescription", m_RSGMassLossPrescription.type, true); }

    bool                                        ScaleCHEMassLossWithSurfaceHeliumAbundance() const                      { return OPT_VALUE("scale-CHE-mass-loss-with-surface-helium-abundance", m_ScaleCHEMassLossWithSurfaceHeliumAbundance, false); }
    double                                      ScaleTerminalWindVelocityWithMetallicityPower() const                   { return OPT_VALUE("scale-terminal-wind-velocity-with-metallicity-power", m_ScaleTerminalWindVelocityWithMetallicityPower, true);}
    double                                      SemiMajorAxis() const                                                   { return OPT_VALUE("semi-major-axis", m_SemiMajorAxis, true); }
    SEMI_MAJOR_AXIS_DISTRIBUTION                SemiMajorAxisDistribution() const                                       { return OPT_VALUE("semi-major-axis-distribution", m_SemiMajorAxisDistribution.type, true); }
    double                                      SemiMajorAxisDistributionMax() const                                    { return OPT_VALUE("semi-major-axis-max", m_SemiMajorAxisDistributionMax, true); }
    double                                      SemiMajorAxisDistributionMin() const                                    { return OPT_VALUE("semi-major-axis-min", m_SemiMajorAxisDistributionMin, true); }

    void                                        ShowHelp()                                                              { PrintOptionHelp(!m_CmdLine.optionValues.m_ShortHelp); }

    double                                      SN_MeanAnomaly1() const                                                 { return OPT_VALUE("kick-mean-anomaly-1", m_KickMeanAnomaly1, true); }
    double                                      SN_MeanAnomaly2() const                                                 { return OPT_VALUE("kick-mean-anomaly-2", m_KickMeanAnomaly2, true); }
    double                                      SN_Phi1() const                                                         { return OPT_VALUE("kick-phi-1", m_KickPhi1, true); }
    double                                      SN_Phi2() const                                                         { return OPT_VALUE("kick-phi-2", m_KickPhi2, true); }
    double                                      SN_Theta1() const                                                       { return OPT_VALUE("kick-theta-1", m_KickTheta1, true); }
    double                                      SN_Theta2() const                                                       { return OPT_VALUE("kick-theta-2", m_KickTheta2, true); }

    bool                                        StoreInputFiles() const                                                 { return m_CmdLine.optionValues.m_StoreInputFiles; }
    bool                                        SwitchLog() const                                                       { return m_CmdLine.optionValues.m_SwitchLog; }

    ZETA_PRESCRIPTION                           StellarZetaPrescription() const                                         { return OPT_VALUE("stellar-zeta-prescription", m_StellarZetaPrescription.type, true); }

    TIDES_PRESCRIPTION                          TidesPrescription() const                                               { return OPT_VALUE("tides-prescription", m_TidesPrescription.type, true); }

    std::string                                 TimestepsFileName() const                                               { return OPT_VALUE("timesteps-filename", m_TimestepsFileName, true); }
    double                                      TimestepMultiplier() const                                              { return m_CmdLine.optionValues.m_TimestepMultiplier; }

    bool                                        UseFixedUK() const                                                      { return (m_GridLine.optionValues.m_UseFixedUK || m_CmdLine.optionValues.m_UseFixedUK); }
    bool                                        UseMassLoss() const                                                     { return OPT_VALUE("use-mass-loss", m_UseMassLoss, true); }
    bool                                        UseMassTransfer() const                                                 { return OPT_VALUE("use-mass-transfer", m_UseMassTransfer, true); }
    bool                                        UsePairInstabilitySupernovae() const                                    { return OPT_VALUE("pair-instability-supernovae", m_UsePairInstabilitySupernovae, true); }
    bool                                        UsePulsationalPairInstability() const                                   { return OPT_VALUE("pulsational-pair-instability", m_UsePulsationalPairInstability, true); }

    VMS_MASS_LOSS_PRESCRIPTION                  VMSMassLossPrescription() const                                         { return OPT_VALUE("VMS-mass-loss-prescription", m_VMSMassLossPrescription.type, true); }
    double                                      WolfRayetFactor() const                                                 { return OPT_VALUE("wolf-rayet-multiplier", m_WolfRayetFactor, true); }
    WR_MASS_LOSS_PRESCRIPTION                   WRMassLossPrescription() const                                          { return OPT_VALUE("WR-mass-loss-prescription", m_WRMassLossPrescription.type, true); }
    std::string                                 YAMLfilename() const                                                    { return m_CmdLine.optionValues.m_YAMLfilename; }
    std::string                                 YAMLtemplate() const                                                    { return m_CmdLine.optionValues.m_YAMLtemplate; }

    double                                      ZetaRadiativeEnvelopeGiant() const                                      { return OPT_VALUE("zeta-radiative-envelope-giant", m_ZetaRadiativeEnvelopeGiant, true); }
    double                                      ZetaMainSequence() const                                                { return OPT_VALUE("zeta-main-sequence", m_ZetaMainSequence, true); }
    double                                      ZetaAdiabaticArbitrary() const                                          { return OPT_VALUE("zeta-adiabatic-arbitrary", m_ZetaAdiabaticArbitrary, true); }

};

#endif // __Options_H__
