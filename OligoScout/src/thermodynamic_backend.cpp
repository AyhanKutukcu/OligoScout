/*
 * THERMODYNAMIC_BACKEND_V1
 */

#include "primerpair/thermodynamic_backend.hpp"

#include <cerrno>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <filesystem>
#include <iomanip>
#include <sstream>
#include <stdexcept>
#include <string>
#include <string_view>

#include <sys/wait.h>
#include <unistd.h>

namespace primerpair {

namespace {

bool executable_file(
    const std::filesystem::path& path
) {
    return
        std::filesystem::is_regular_file(path) &&
        ::access(
            path.c_str(),
            X_OK
        ) == 0;
}


std::string shell_quote(
    const std::string_view value
) {
    std::string output;

    output.reserve(
        value.size() + 2
    );

    output.push_back('\'');

    for (const char c : value) {

        if (c == '\'') {
            output += "'\"'\"'";
        } else {
            output.push_back(c);
        }
    }

    output.push_back('\'');

    return output;
}


void validate_sequence(
    const std::string_view sequence
) {
    if (sequence.empty()) {
        throw std::invalid_argument(
            "Thermodynamic sequence cannot be empty."
        );
    }

    for (const char c : sequence) {

        switch (c) {
            case 'A':
            case 'C':
            case 'G':
            case 'T':
            case 'a':
            case 'c':
            case 'g':
            case 't':
                break;

            default:
                throw std::invalid_argument(
                    "Thermodynamic sequence must contain "
                    "only canonical DNA bases A/C/G/T."
                );
        }
    }
}


void validate_conditions(
    const ThermodynamicConditions& c
) {
    if (
        !std::isfinite(c.monovalent_mM) ||
        !std::isfinite(c.divalent_mM) ||
        !std::isfinite(c.dntp_mM) ||
        !std::isfinite(c.dna_nM) ||
        !std::isfinite(c.temperature_celsius)
    ) {
        throw std::invalid_argument(
            "Thermodynamic conditions must be finite."
        );
    }

    if (
        c.monovalent_mM < 0.0 ||
        c.divalent_mM < 0.0 ||
        c.dntp_mM < 0.0
    ) {
        throw std::invalid_argument(
            "Ion and dNTP concentrations cannot "
            "be negative."
        );
    }

    if (c.dna_nM <= 0.0) {
        throw std::invalid_argument(
            "DNA concentration must be > 0."
        );
    }

    if (
        c.temperature_celsius <= -273.15
    ) {
        throw std::invalid_argument(
            "Temperature must be above "
            "absolute zero."
        );
    }

    if (c.max_loop > 30) {
        throw std::invalid_argument(
            "ntthal max_loop must be <= 30."
        );
    }

    if (
        c.tm_method < 0 ||
        c.tm_method > 2
    ) {
        throw std::invalid_argument(
            "Invalid oligotm Tm method."
        );
    }

    if (
        c.salt_correction < 0 ||
        c.salt_correction > 2
    ) {
        throw std::invalid_argument(
            "Invalid oligotm salt correction."
        );
    }
}


std::string alignment_name(
    const ThermodynamicAlignment alignment
) {
    switch (alignment) {
        case ThermodynamicAlignment::Any:
            return "ANY";

        case ThermodynamicAlignment::End1:
            return "END1";

        case ThermodynamicAlignment::End2:
            return "END2";
    }

    throw std::logic_error(
        "Unknown thermodynamic alignment type."
    );
}


double parse_numeric_output(
    const std::string& raw
) {
    std::size_t consumed = 0;

    double value = 0.0;

    try {
        value = std::stod(
            raw,
            &consumed
        );
    } catch (
        const std::exception&
    ) {
        throw std::runtime_error(
            "Thermodynamic backend returned "
            "non-numeric output: " +
            raw
        );
    }

    for (
        std::size_t i = consumed;
        i < raw.size();
        ++i
    ) {
        const char c = raw[i];

        if (
            c != ' ' &&
            c != '\t' &&
            c != '\r' &&
            c != '\n'
        ) {
            throw std::runtime_error(
                "Unexpected trailing thermodynamic "
                "backend output: " +
                raw
            );
        }
    }

    if (!std::isfinite(value)) {
        throw std::runtime_error(
            "Thermodynamic backend returned "
            "non-finite value."
        );
    }

    return value;
}


double run_numeric_command(
    const std::string& command
) {
    const std::string wrapped =
        "LC_ALL=C " +
        command +
        " 2>&1";

    FILE* pipe =
        ::popen(
            wrapped.c_str(),
            "r"
        );

    if (pipe == nullptr) {
        throw std::runtime_error(
            "Failed to start thermodynamic backend."
        );
    }

    std::string output;

    char buffer[4096];

    while (
        std::fgets(
            buffer,
            static_cast<int>(
                sizeof(buffer)
            ),
            pipe
        ) != nullptr
    ) {
        output += buffer;
    }

    const int status =
        ::pclose(pipe);

    if (
        status == -1 ||
        !WIFEXITED(status) ||
        WEXITSTATUS(status) != 0
    ) {
        throw std::runtime_error(
            "Thermodynamic backend command failed. "
            "Output: " +
            output
        );
    }

    return
        parse_numeric_output(
            output
        );
}


std::string config_path_with_separator(
    const std::filesystem::path& path
) {
    std::string value =
        path.string();

    if (
        !value.empty() &&
        value.back() != '/'
    ) {
        value.push_back('/');
    }

    return value;
}


std::string numeric(
    const double value
) {
    std::ostringstream out;

    out
        << std::setprecision(17)
        << value;

    return out.str();
}

}  // namespace


Primer3ThermodynamicBackend::
Primer3ThermodynamicBackend(
    std::filesystem::path ntthal_path,
    std::filesystem::path oligotm_path,
    std::filesystem::path
        thermodynamic_parameters_path
)
    : ntthal_path_(
          std::move(ntthal_path)
      ),
      oligotm_path_(
          std::move(oligotm_path)
      ),
      thermodynamic_parameters_path_(
          std::move(
              thermodynamic_parameters_path
          )
      )
{
    if (!executable_file(ntthal_path_)) {
        throw std::invalid_argument(
            "ntthal executable is unavailable: " +
            ntthal_path_.string()
        );
    }

    if (!executable_file(oligotm_path_)) {
        throw std::invalid_argument(
            "oligotm executable is unavailable: " +
            oligotm_path_.string()
        );
    }

    if (
        !std::filesystem::is_directory(
            thermodynamic_parameters_path_
        )
    ) {
        throw std::invalid_argument(
            "Primer3 thermodynamic parameter "
            "directory is unavailable: " +
            thermodynamic_parameters_path_
                .string()
        );
    }
}


bool Primer3ThermodynamicBackend::
available() const noexcept
{
    return
        executable_file(ntthal_path_) &&
        executable_file(oligotm_path_) &&
        std::filesystem::is_directory(
            thermodynamic_parameters_path_
        );
}


double Primer3ThermodynamicBackend::
oligo_tm(
    const std::string_view sequence,
    const ThermodynamicConditions& conditions
) const
{
    validate_sequence(sequence);
    validate_conditions(conditions);

    /*
     * Current Primer3 oligotm CLI supports
     * oligos from 2 through 36 bases.
     */
    if (
        sequence.size() < 2 ||
        sequence.size() > 36
    ) {
        throw std::invalid_argument(
            "oligotm sequence length must "
            "be in the range 2..36."
        );
    }

    std::ostringstream command;

    command
        << shell_quote(
               oligotm_path_.string()
           )
        << " -mv "
        << numeric(
               conditions.monovalent_mM
           )
        << " -dv "
        << numeric(
               conditions.divalent_mM
           )
        << " -n "
        << numeric(
               conditions.dntp_mM
           )
        << " -d "
        << numeric(
               conditions.dna_nM
           )
        << " -tp "
        << conditions.tm_method
        << " -sc "
        << conditions.salt_correction
        << ' '
        << shell_quote(sequence);

    return
        run_numeric_command(
            command.str()
        );
}


double Primer3ThermodynamicBackend::
duplex_tm(
    const std::string_view sequence1,
    const std::string_view sequence2,
    const ThermodynamicAlignment alignment,
    const ThermodynamicConditions& conditions
) const
{
    validate_sequence(sequence1);
    validate_sequence(sequence2);
    validate_conditions(conditions);

    std::ostringstream command;

    command
        << shell_quote(
               ntthal_path_.string()
           )
        << " -r"
        << " -path "
        << shell_quote(
               config_path_with_separator(
                   thermodynamic_parameters_path_
               )
           )
        << " -mv "
        << numeric(
               conditions.monovalent_mM
           )
        << " -dv "
        << numeric(
               conditions.divalent_mM
           )
        << " -n "
        << numeric(
               conditions.dntp_mM
           )
        << " -d "
        << numeric(
               conditions.dna_nM
           )
        << " -t "
        << numeric(
               conditions.temperature_celsius
           )
        << " -maxloop "
        << conditions.max_loop
        << " -a "
        << alignment_name(
               alignment
           )
        << " -s1 "
        << shell_quote(sequence1)
        << " -s2 "
        << shell_quote(sequence2);

    return
        run_numeric_command(
            command.str()
        );
}


double Primer3ThermodynamicBackend::
hairpin_tm(
    const std::string_view sequence,
    const ThermodynamicConditions& conditions
) const
{
    validate_sequence(sequence);
    validate_conditions(conditions);

    std::ostringstream command;

    command
        << shell_quote(
               ntthal_path_.string()
           )
        << " -r"
        << " -path "
        << shell_quote(
               config_path_with_separator(
                   thermodynamic_parameters_path_
               )
           )
        << " -mv "
        << numeric(
               conditions.monovalent_mM
           )
        << " -dv "
        << numeric(
               conditions.divalent_mM
           )
        << " -n "
        << numeric(
               conditions.dntp_mM
           )
        << " -d "
        << numeric(
               conditions.dna_nM
           )
        << " -t "
        << numeric(
               conditions.temperature_celsius
           )
        << " -maxloop "
        << conditions.max_loop
        << " -a HAIRPIN"
        << " -s1 "
        << shell_quote(sequence);

    return
        run_numeric_command(
            command.str()
        );
}


PrimerThermodynamicProfile
Primer3ThermodynamicBackend::
profile_primer(
    const std::string_view sequence,
    const ThermodynamicConditions& conditions
) const
{
    PrimerThermodynamicProfile output;

    output.oligo_tm_celsius =
        oligo_tm(
            sequence,
            conditions
        );

    output.hairpin_tm_celsius =
        hairpin_tm(
            sequence,
            conditions
        );

    output.homodimer_any_tm_celsius =
        duplex_tm(
            sequence,
            sequence,
            ThermodynamicAlignment::Any,
            conditions
        );

    output.homodimer_end1_tm_celsius =
        duplex_tm(
            sequence,
            sequence,
            ThermodynamicAlignment::End1,
            conditions
        );

    output.homodimer_end2_tm_celsius =
        duplex_tm(
            sequence,
            sequence,
            ThermodynamicAlignment::End2,
            conditions
        );

    return output;
}


PrimerPairThermodynamicProfile
Primer3ThermodynamicBackend::
profile_pair(
    const std::string_view left_primer,
    const std::string_view right_primer,
    const ThermodynamicConditions& conditions
) const
{
    PrimerPairThermodynamicProfile output;

    output.left =
        profile_primer(
            left_primer,
            conditions
        );

    output.right =
        profile_primer(
            right_primer,
            conditions
        );

    output.heterodimer_any_tm_celsius =
        duplex_tm(
            left_primer,
            right_primer,
            ThermodynamicAlignment::Any,
            conditions
        );

    output.heterodimer_end1_tm_celsius =
        duplex_tm(
            left_primer,
            right_primer,
            ThermodynamicAlignment::End1,
            conditions
        );

    output.heterodimer_end2_tm_celsius =
        duplex_tm(
            left_primer,
            right_primer,
            ThermodynamicAlignment::End2,
            conditions
        );

    return output;
}

}  // namespace primerpair
