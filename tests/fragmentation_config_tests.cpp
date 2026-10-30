#include "config_io.hpp"

#include <iostream>

namespace {

bool check(bool condition, const char* message)
{
    if (!condition) std::cerr << "FAILED: " << message << '\n';
    return condition;
}

void write_file(const fs::path& path, const std::string& text)
{
    std::ofstream output(path);
    output << text;
    if (!output) throw std::runtime_error("cannot write test input: " + path.string());
}

} // namespace

int main(int argc, char** argv)
{
    if (argc != 9) return 1;
    try {
        const fs::path directory(argv[1]);
        fs::create_directories(directory);
        write_file(directory / "rocks.txt", "unused.pts\n");
        write_file(directory / "starts.txt", "0 0 0\n");
        const fs::path scenario = directory / "scenario.xml";
        const std::string base = R"xml(
<rockfall3d/>
<terrain dem="unused.tif"><interaction type="coulomb_friction" mu="0.5"/></terrain>
<rocks list="rocks.txt" density="2500" random_orientations="1"/>
<starts list="starts.txt"/>
<simulation dt="0.002" output_dir="output"/>
<output><trajectory format="binary" metric="time"/></output>
<stop max_time="1" max_steps="1000" settle_min_steps="100" settle_linear_speed="0.1" settle_angular_speed="0.1"/>
<integrator type="schweizer"/>
<solver type="sor-prox" alpha="1" atol="1e-8" rtol="1e-6" maxiter="1000"/>
)xml";
        auto load = [&](const std::string& model) {
            write_file(scenario, base + model);
            return load_batch_simulation_config(scenario);
        };
        const auto described = load(argv[2]);
        bool ok = check(described.fragmentation.has_value(), "advertised XML enables fragmentation");
        if (!described.fragmentation) return 1;
        const auto& config = *described.fragmentation;
        ok &= check(config.thresholdSpecificEnergy == parse_number<Real>(argv[3], "threshold"), "threshold XML mapping");
        ok &= check(config.halfSpecificEnergy == parse_number<Real>(argv[4], "saturation"), "saturation XML mapping");
        ok &= check(config.exponent == parse_number<Real>(argv[5], "exponent"), "exponent XML mapping");
        ok &= check(config.maxRemovedFraction == parse_number<Real>(argv[6], "fraction"), "fraction XML mapping");
        ok &= check(config.shapeOutputDirectory == fs::path("fragment_shapes"), "default shape directory");
        for (const auto field : { OutputFieldId::Mass, OutputFieldId::ShapeIndex }) {
            ok &= check(std::find(described.output.fields.begin(), described.output.fields.end(), field)
                != described.output.fields.end(), "automatic fragmentation output fields");
        }

        const auto defaults = load("<fragmentation threshold_specific_energy=\"0\" half_specific_energy=\"1\"/>");
        ok &= check(defaults.fragmentation->exponent == parse_number<Real>(argv[7], "default exponent"), "advertised exponent default");
        ok &= check(defaults.fragmentation->maxRemovedFraction == parse_number<Real>(argv[8], "default fraction"), "advertised fraction default");
        const auto disabled = load("");
        ok &= check(!disabled.fragmentation, "omitted element disables fragmentation");
        ok &= check(disabled.output.fields.size() == 1, "disabled model does not add output fields");
        const auto boundary = load("<fragmentation threshold_specific_energy=\"0\" half_specific_energy=\"1e-12\" exponent=\"1e-12\" max_removed_fraction=\"0.5\"/>");
        ok &= check(boundary.fragmentation.has_value(), "valid inclusive/exclusive boundaries accepted");

        for (const std::string invalid : {
            "<fragmentation half_specific_energy=\"100\"/>",
            "<fragmentation threshold_specific_energy=\"0\"/>",
            "<fragmentation threshold_specific_energy=\"-1\" half_specific_energy=\"100\"/>",
            "<fragmentation threshold_specific_energy=\"1\" half_specific_energy=\"1\"/>",
            "<fragmentation threshold_specific_energy=\"2\" half_specific_energy=\"1\"/>",
            "<fragmentation threshold_specific_energy=\"0\" half_specific_energy=\"1\" exponent=\"0\"/>",
            "<fragmentation threshold_specific_energy=\"0\" half_specific_energy=\"1\" max_removed_fraction=\"0\"/>",
            "<fragmentation threshold_specific_energy=\"0\" half_specific_energy=\"1\" max_removed_fraction=\"0.500001\"/>",
            "<fragmentation threshold_specific_energy=\"nan\" half_specific_energy=\"100\"/>",
            "<fragmentation threshold_specific_energy=\"0\" half_specific_energy=\"inf\"/>",
            "<fragmentation threshold_specific_energy=\"0\" half_specific_energy=\"1\" exponent=\"inf\"/>",
            "<fragmentation threshold_specific_energy=\"0\" half_specific_energy=\"1\" max_removed_fraction=\"nan\"/>"
        }) {
            bool rejected = false;
            try { load(invalid); }
            catch (const std::runtime_error&) { rejected = true; }
            ok &= check(rejected, "invalid fragmentation parameters rejected by real parser");
        }
        return ok ? 0 : 1;
    }
    catch (const std::exception& error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
