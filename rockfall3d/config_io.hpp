#pragma once
#include "math_alias.hpp"
#include "terrain.hpp"
#include "time_stepper.hpp"
#include "trajectory_format.hpp"

#include <algorithm>
#include <cctype>
#include <filesystem>
#include <fstream>
#include <optional>
#include <random>
#include <regex>
#include <sstream>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <vector>

// XML file parser

namespace fs = std::filesystem;

// solver parameters
struct SolverConfig
{
    ContactSolverType type{ ContactSolverType::JorProx };
    Real alpha{ Real(0.5) };
    Real atol{ Real(1e-8) };
    Real rtol{ Real(1e-6) };
    std::uint32_t maxiter{ 10000 };
};

struct StartPoint
{
    Real x{ Real(0) };
    Real y{ Real(0) };
    Real z{ Real(0) };
    std::optional<Real> dz;

    bool isTerrainRelative() const noexcept
    {
        return dz.has_value();
    }
};

struct BatchSimulationConfig
{
    // Resolved file-system inputs and parsed simulation settings.
    fs::path configPath;
    fs::path demPath;
    fs::path rockListPath;
    fs::path startPointsPath;
    fs::path outputDir;
    std::optional<fs::path> roughnessRasterPath;
    std::optional<fs::path> terrainExportPath;
    std::optional<Real> terrainExportResolution;
    bool writeTimingReport{ false };
    bool rasterOutputEnabled{ false };
    Real rasterResolution{ Real(2.0) };

    std::vector<fs::path> rockPaths;
    std::vector<StartPoint> startPoints;
    StepFunc roughnessStepFunc{ StepFunc::Smoothstep };

    std::uint32_t randomOrientations{ 1 };
    Real rockDensity{ Real(2500.0) };
    Real mu{ Real(0.5) };
    Real dt{ Real(0.002) };
    std::optional<Real> dumpStep;
    std::optional<std::uint64_t> randomSeed;
    std::optional<std::uint64_t> orientationSeed;
    std::optional<std::uint64_t> roughnessSeed;
    std::uint32_t parallelThreads{ 1 };
    OutputSettings output;

    SolverConfig solver;
    TimeStepperStopCriteria stop;
};

// removes whitespace at the beginning and end of a string and returns cleaned copy
inline std::string trim_copy(const std::string& value)
{
    const auto begin = std::find_if_not(value.begin(), value.end(), [](unsigned char c) { return std::isspace(c) != 0; });
    const auto end = std::find_if_not(value.rbegin(), value.rend(), [](unsigned char c) { return std::isspace(c) != 0; }).base();
    return begin < end ? std::string(begin, end) : std::string();
}

// converts all characters in a string to lowercase (used for case-insensitive parsing)
inline std::string to_lower_copy(std::string value) 
{
    std::transform(value.begin(), value.end(), value.begin(), [](unsigned char c) {
        return static_cast<char>(std::tolower(c));
        });
    return value;
}

// detects Windows server paths
inline bool looks_like_unc_path(const std::string& raw) noexcept
{
    return raw.rfind("//", 0) == 0 || raw.rfind("\\\\", 0) == 0;
}

inline fs::path parse_config_path(const std::string& raw)
{
    std::string value = trim_copy(raw);
#ifdef _WIN32
    if (looks_like_unc_path(value)) {
        std::replace(value.begin(), value.end(), '/', '\\');
    }
#endif
    return fs::path(value);
}

// Resolve config paths relative to the XML location while still accepting
// absolute paths and UNC/network paths.
inline fs::path resolve_relative(const fs::path& baseDir, const std::string& raw)
{
    fs::path path = parse_config_path(raw);
    if (path.empty()) {
        return path;
    }
    if (path.is_absolute() || looks_like_unc_path(trim_copy(raw))) {
        return path.lexically_normal();
    }
    return (baseDir / path).lexically_normal();
}

inline std::string read_text_file(const fs::path& path)
{
    std::ifstream in(path);
    if (!in) {
        throw std::runtime_error("cannot open config file: " + path.string());
    }
    std::ostringstream buffer;
    buffer << in.rdbuf();
    return buffer.str();
}

// Find self-closing XML tags i.e. <simulation dt="0.002" output_dir="output"/>
// Returns all attributes as a map
// The XML format is intentionally small, so a lightweight regex-based parser is
// sufficient here and keeps dependencies low for the command-line tool.
inline std::optional<std::unordered_map<std::string, std::string>> find_tag_attributes(
    const std::string& xml,
    const std::string& tag)
{
    const std::regex tagExpr("<\\s*" + tag + R"(\b([^>]*)/\s*>)", std::regex::icase);
    std::smatch tagMatch;
    if (!std::regex_search(xml, tagMatch, tagExpr)) {
        return std::nullopt;
    }

    std::unordered_map<std::string, std::string> attrs;
    const std::regex attrExpr("([A-Za-z0-9_:\\-]+)\\s*=\\s*\"([^\"]*)\"");
    const std::string attrText = tagMatch[1].str();
    for (std::sregex_iterator it(attrText.begin(), attrText.end(), attrExpr), end; it != end; ++it) {
        attrs[to_lower_copy((*it)[1].str())] = (*it)[2].str();
    }
    return attrs;
}

// fetches required attribute from a parsed XML tag
// if attribute exists, it returns attribute value
// if attribute does not exist, it throws error
inline const std::string& require_attr(const std::unordered_map<std::string, std::string>& attrs,
    const std::string& name,
    const std::string& tag)
{
    const auto it = attrs.find(to_lower_copy(name));
    if (it == attrs.end()) {
        throw std::runtime_error("missing attribute '" + name + "' on <" + tag + "/> tag");
    }
    return it->second;
}

// fetches optional XML attribute
// if attribute exists, return value
// if it does not exist, return std::nullopt
inline std::optional<std::string> optional_attr(const std::unordered_map<std::string, std::string>& attrs,
    const std::string& name)
{
    const auto it = attrs.find(to_lower_copy(name));
    if (it == attrs.end()) {
        return std::nullopt;
    }
    return it->second;
}

// converts a string into a number of type T and validates that the whole string is a valid number
template <typename T>
inline T parse_number(const std::string& text, const std::string& fieldName)
{
    std::istringstream in(trim_copy(text));
    T value{};
    in >> value;
    if (!in || !in.eof()) {
        throw std::runtime_error("invalid numeric value for " + fieldName + ": " + text);
    }
    return value;
}

// converts a string into a boolean and validates boolean
inline bool parse_bool(const std::string& text, const std::string& fieldName)
{
    const std::string value = to_lower_copy(trim_copy(text));
    if (value == "1" || value == "true" || value == "yes" || value == "on") {
        return true;
    }
    if (value == "0" || value == "false" || value == "no" || value == "off") {
        return false;
    }
    throw std::runtime_error("invalid boolean value for " + fieldName + ": " + text);
}

// parse step function keyword
inline StepFunc parse_step_func(const std::string& text)
{
    const std::string value = to_lower_copy(trim_copy(text));
    if (value == "linear") return StepFunc::Linear;
    if (value == "smoothstep") return StepFunc::Smoothstep;
    if (value == "smootherstep") return StepFunc::Smootherstep;
    throw std::runtime_error("unknown roughness step function: " + text);
}

// parse solver type keyword
inline ContactSolverType parse_solver_type(const std::string& text)
{
    const std::string value = to_lower_copy(trim_copy(text));
    if (value == "jor" || value == "jor-prox" || value == "jor_prox" || value == "jorprox") return ContactSolverType::JorProx;
    if (value == "sor" || value == "sor-prox" || value == "sor_prox" || value == "sorprox") return ContactSolverType::SorProx;
    throw std::runtime_error("unknown solver type: " + text);
}

// parse output format (txt vs. binary)
inline OutputFormat parse_output_format(const std::string& text)
{
    const std::string value = to_lower_copy(trim_copy(text));
    if (value.empty() || value == "text" || value == "legacy" || value == "legacy-text" || value == "legacy_text") {
        return OutputFormat::LegacyText;
    }
    if (value == "binary" || value == "binary-sharded" || value == "binary_sharded" || value == "binary-shard" || value == "binary_shard") {
        return OutputFormat::BinarySharded;
    }
    throw std::runtime_error("unknown output format: " + text);
}

// parses which fields should be dumped to binary trajectory files
inline void append_output_field_alias(std::vector<OutputFieldId>& fields, const std::string& token)
{
    const std::string value = to_lower_copy(trim_copy(token));
    if (value.empty()) {
        return;
    }

    auto append_unique = [&](OutputFieldId field) {
        if (std::find(fields.begin(), fields.end(), field) == fields.end()) {
            fields.push_back(field);
        }
    };

    if (value == "time" || value == "t") { append_unique(OutputFieldId::Time); return; }
    if (value == "x") { append_unique(OutputFieldId::X); return; }
    if (value == "y") { append_unique(OutputFieldId::Y); return; }
    if (value == "z") { append_unique(OutputFieldId::Z); return; }
    if (value == "qw") { append_unique(OutputFieldId::Qw); return; }
    if (value == "qx") { append_unique(OutputFieldId::Qx); return; }
    if (value == "qy") { append_unique(OutputFieldId::Qy); return; }
    if (value == "qz") { append_unique(OutputFieldId::Qz); return; }
    if (value == "e" || value == "q" || value == "quat" || value == "quaternion") {
        append_unique(OutputFieldId::Qw);
        append_unique(OutputFieldId::Qx);
        append_unique(OutputFieldId::Qy);
        append_unique(OutputFieldId::Qz);
        return;
    }
    if (value == "vx") { append_unique(OutputFieldId::Vx); return; }
    if (value == "vy") { append_unique(OutputFieldId::Vy); return; }
    if (value == "vz") { append_unique(OutputFieldId::Vz); return; }
    if (value == "vtot") { append_unique(OutputFieldId::Vtot); return; }
    if (value == "wx") { append_unique(OutputFieldId::Wx); return; }
    if (value == "wy") { append_unique(OutputFieldId::Wy); return; }
    if (value == "wz") { append_unique(OutputFieldId::Wz); return; }
    if (value == "wtot") { append_unique(OutputFieldId::Wtot); return; }
    if (value == "ekin_trans" || value == "ekintrans") { append_unique(OutputFieldId::EkinTrans); return; }
    if (value == "ekin_rot" || value == "ekinrot") { append_unique(OutputFieldId::EkinRot); return; }
    if (value == "ekin_tot" || value == "ekin" || value == "ekintot") { append_unique(OutputFieldId::EkinTot); return; }
    if (value == "jumph" || value == "jump_h") { append_unique(OutputFieldId::JumpH); return; }

    throw std::runtime_error("unknown output variable: " + token);
}

inline std::vector<OutputFieldId> parse_output_fields(const std::string& text)
{
    std::vector<OutputFieldId> fields;
    std::string token;
    std::stringstream ss(text);
    while (std::getline(ss, token, ',')) {
        append_output_field_alias(fields, token);
    }
    if (fields.empty()) {
        throw std::runtime_error("expected at least one output variable");
    }
    return fields;
}

inline std::vector<std::string> read_nonempty_lines(const fs::path& path)
{
    std::ifstream in(path);
    if (!in) {
        throw std::runtime_error("cannot open list file: " + path.string());
    }

    std::vector<std::string> lines;
    std::string line;
    while (std::getline(in, line)) {
        const std::size_t commentPos = line.find('#');
        if (commentPos != std::string::npos) line.erase(commentPos);
        line = trim_copy(line);
        if (!line.empty()) lines.push_back(line);
    }
    return lines;
}

// reads a list file containing paths (e.g. rocks.txt)
inline std::vector<fs::path> load_path_list_file(const fs::path& path)
{
    const auto lines = read_nonempty_lines(path);
    std::vector<fs::path> result;
    result.reserve(lines.size());
    for (const auto& line : lines) result.push_back(resolve_relative(path.parent_path(), line));
    return result;
}

// Reads start-point files in one of two forms:
// - X Y Z       : use absolute world-space release coordinates directly
// - X Y Z dz    : ignore the given Z and instead place the rock dz above the
//                 terrain (including optional roughness) at X/Y
inline std::vector<StartPoint> load_start_points_file(const fs::path& path)
{
    const auto lines = read_nonempty_lines(path);
    std::vector<StartPoint> result;
    result.reserve(lines.size());
    for (const auto& line : lines) {
        std::istringstream in(line);
        std::vector<Real> values;
        Real value{};
        while (in >> value) {
            values.push_back(value);
        }
        if (!in.eof() || (values.size() != 3 && values.size() != 4)) {
            throw std::runtime_error("invalid start point line in " + path.string() + ": " + line);
        }

        StartPoint point;
        point.x = values[0];
        point.y = values[1];
        point.z = values[2];
        if (values.size() == 4) {
            point.dz = values[3];
        }
        result.push_back(point);
    }
    return result;
}

inline Eigen::Matrix<Real, 4, 1> make_random_quaternion_wxyz(std::mt19937_64& gen)
{
    std::uniform_real_distribution<Real> dist(Real(0), Real(1));
    const Real u1 = dist(gen);
    const Real u2 = dist(gen);
    const Real u3 = dist(gen);
    const Real twoPi = Real(2) * Real(3.1415926535897932384626433832795);

    const Real s1 = std::sqrt(Real(1) - u1);
    const Real s2 = std::sqrt(u1);
    const Real theta1 = twoPi * u2;
    const Real theta2 = twoPi * u3;

    const Real qx = s1 * std::sin(theta1);
    const Real qy = s1 * std::cos(theta1);
    const Real qz = s2 * std::sin(theta2);
    const Real qw = s2 * std::cos(theta2);

    return (Eigen::Matrix<Real, 4, 1>() << qw, qx, qy, qz).finished();
}

inline BatchSimulationConfig load_batch_simulation_config(const fs::path& configPath)
{
    BatchSimulationConfig cfg;
    cfg.configPath = configPath;
    const fs::path baseDir = configPath.parent_path();
    const std::string xml = read_text_file(configPath);

    const auto terrainTag = find_tag_attributes(xml, "terrain");
    const auto rocksTag = find_tag_attributes(xml, "rocks");
    const auto startsTag = find_tag_attributes(xml, "starts");
    const auto simulationTag = find_tag_attributes(xml, "simulation");
    const auto stopTag = find_tag_attributes(xml, "stop");
    const auto solverTag = find_tag_attributes(xml, "solver");
    const auto roughnessTag = find_tag_attributes(xml, "roughness");
    const auto outputTag = find_tag_attributes(xml, "output");
    const auto rasterTag = find_tag_attributes(xml, "raster");
    const auto debugTag = find_tag_attributes(xml, "debug");

    if (!terrainTag || !rocksTag || !startsTag || !simulationTag || !stopTag || !solverTag) {
        throw std::runtime_error("config XML must contain self-closing <terrain/>, <rocks/>, <starts/>, <simulation/>, <stop/>, and <solver/> tags");
    }

    cfg.demPath = resolve_relative(baseDir, require_attr(*terrainTag, "dem", "terrain"));
    cfg.rockListPath = resolve_relative(baseDir, require_attr(*rocksTag, "list", "rocks"));
    cfg.rockDensity = parse_number<Real>(require_attr(*rocksTag, "density", "rocks"), "rock density");
    cfg.randomOrientations = parse_number<std::uint32_t>(require_attr(*rocksTag, "random_orientations", "rocks"), "random_orientations");
    if (cfg.randomOrientations == 0) throw std::runtime_error("random_orientations must be at least 1");

    cfg.startPointsPath = resolve_relative(baseDir, require_attr(*startsTag, "list", "starts"));

    cfg.mu = parse_number<Real>(require_attr(*simulationTag, "mu", "simulation"), "mu");
    cfg.dt = parse_number<Real>(require_attr(*simulationTag, "dt", "simulation"), "dt");
    if (const auto dumpAttr = optional_attr(*simulationTag, "dump_step")) {
        cfg.dumpStep = parse_number<Real>(*dumpAttr, "dump_step");
        if (*cfg.dumpStep <= Real(0)) {
            throw std::runtime_error("dump_step must be positive");
        }
    }
    cfg.outputDir = resolve_relative(baseDir, require_attr(*simulationTag, "output_dir", "simulation"));
    if (const auto workersAttr = optional_attr(*simulationTag, "parallel_threads")) {
        cfg.parallelThreads = parse_number<std::uint32_t>(*workersAttr, "parallel_threads");
    }
    else if (const auto workersAttr = optional_attr(*simulationTag, "parallel_processes")) {
        cfg.parallelThreads = parse_number<std::uint32_t>(*workersAttr, "parallel_processes");
    }
    if (const auto seedAttr = optional_attr(*simulationTag, "random_seed")) {
        cfg.randomSeed = parse_number<std::uint64_t>(*seedAttr, "random_seed");
    }
    if (const auto seedAttr = optional_attr(*simulationTag, "orientation_seed")) {
        cfg.orientationSeed = parse_number<std::uint64_t>(*seedAttr, "orientation_seed");
    }
    if (const auto seedAttr = optional_attr(*simulationTag, "roughness_seed")) {
        cfg.roughnessSeed = parse_number<std::uint64_t>(*seedAttr, "roughness_seed");
    }
    if (!cfg.orientationSeed && cfg.randomSeed) cfg.orientationSeed = cfg.randomSeed;
    if (!cfg.roughnessSeed && cfg.randomSeed) cfg.roughnessSeed = cfg.randomSeed;
    if (cfg.parallelThreads == 0) cfg.parallelThreads = 1;

    cfg.stop.maxTime = parse_number<Real>(require_attr(*stopTag, "max_time", "stop"), "max_time");
    cfg.stop.maxSteps = parse_number<std::uint32_t>(require_attr(*stopTag, "max_steps", "stop"), "max_steps");
    cfg.stop.settleMinSteps = parse_number<std::uint32_t>(require_attr(*stopTag, "settle_min_steps", "stop"), "settle_min_steps");
    cfg.stop.settleLinearSpeed = parse_number<Real>(require_attr(*stopTag, "settle_linear_speed", "stop"), "settle_linear_speed");
    cfg.stop.settleAngularSpeed = parse_number<Real>(require_attr(*stopTag, "settle_angular_speed", "stop"), "settle_angular_speed");

    const std::string solverTypeText = optional_attr(*solverTag, "type")
        .value_or(optional_attr(*solverTag, "method")
            .value_or(std::string()));
    if (solverTypeText.empty()) {
        throw std::runtime_error("missing attribute 'type' or 'method' on <solver/> tag");
    }
    cfg.solver.type = parse_solver_type(solverTypeText);
    cfg.solver.alpha = parse_number<Real>(require_attr(*solverTag, "alpha", "solver"), "solver alpha");
    cfg.solver.atol = parse_number<Real>(require_attr(*solverTag, "atol", "solver"), "solver atol");
    cfg.solver.rtol = parse_number<Real>(require_attr(*solverTag, "rtol", "solver"), "solver rtol");
    cfg.solver.maxiter = parse_number<std::uint32_t>(require_attr(*solverTag, "maxiter", "solver"), "solver maxiter");

    if (roughnessTag) {
        cfg.roughnessRasterPath = resolve_relative(
            baseDir,
            require_attr(*roughnessTag, "raster", "roughness"));
        if (const auto stepAttr = optional_attr(*roughnessTag, "step_func")) {
            cfg.roughnessStepFunc = parse_step_func(*stepAttr);
        }
    }

    if (outputTag) {
        cfg.output.enabled = true;
        cfg.output.format = parse_output_format(require_attr(*outputTag, "format", "output"));
        cfg.output.fields = parse_output_fields(require_attr(*outputTag, "variables", "output"));
    }
    else if (!rasterTag) {
        cfg.output.enabled = true;
    }

    if (rasterTag) {
        cfg.rasterOutputEnabled = true;
    }

    if (debugTag) {
        if (const auto pathAttr = optional_attr(*debugTag, "terrain_export")) {
            cfg.terrainExportPath = resolve_relative(baseDir, *pathAttr);
        }
        if (const auto resAttr = optional_attr(*debugTag, "terrain_export_resolution")) {
            cfg.terrainExportResolution = parse_number<Real>(*resAttr, "terrain_export_resolution");
        }
        if (cfg.terrainExportPath.has_value() != cfg.terrainExportResolution.has_value()) {
            throw std::runtime_error("<debug/> must specify both terrain_export and terrain_export_resolution when terrain export is requested");
        }
        if (const auto reportAttr = optional_attr(*debugTag, "timing_report")) {
            cfg.writeTimingReport = parse_bool(*reportAttr, "timing_report");
        }
    }

    cfg.rockPaths = load_path_list_file(cfg.rockListPath);
    cfg.startPoints = load_start_points_file(cfg.startPointsPath);
    if (cfg.rockPaths.empty()) throw std::runtime_error("rock list file contains no rock paths: " + cfg.rockListPath.string());
    if (cfg.startPoints.empty()) throw std::runtime_error("start points file contains no start points: " + cfg.startPointsPath.string());

    if (roughnessTag) {
        if (!cfg.roughnessRasterPath) {
            throw std::runtime_error("<roughness/> requires raster=\"...\"");
        }
    }

    return cfg;
}
