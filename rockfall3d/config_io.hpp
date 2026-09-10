#pragma once
#include "deposition_output.hpp"
#include "fragmentation.hpp"
#include "math_alias.hpp"
#include "raster_output.hpp"
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

// Resolved file-system inputs and parsed simulation settings.
struct BatchSimulationConfig
{
    fs::path configPath;
    fs::path demPath;
    fs::path rockListPath;
    fs::path startPointsPath;
    fs::path outputDir;
    std::optional<fs::path> roughnessRasterPath;
    std::optional<fs::path> terrainExportPath;
    std::optional<Real> terrainExportResolution;
    bool writeTimingReport{ false };
    std::optional<RasterOutputConfig> rasterOutput;
    std::optional<DepositionPointOutputConfig> depositionPointOutput;
    std::optional<FragmentationConfig> fragmentation;

    std::vector<fs::path> rockPaths;
    std::vector<StartPoint> startPoints;
    StepFunc roughnessStepFunc{ StepFunc::Smoothstep };

    std::uint32_t randomOrientations{ 1 };
    Real rockDensity{ Real(2500.0) };
    Real epsilonN{ Real(0) };
    Real epsilonT{ Real(0) };
    Real dt{ Real(0.002) };
    std::optional<Real> dumpStep;
    std::optional<std::uint64_t> randomSeed;
    std::optional<std::uint64_t> orientationSeed;
    std::optional<std::uint64_t> roughnessSeed;
    std::uint32_t parallelThreads{ 1 };
    OutputSettings output;

    SolverConfig solver;
    RotationalIntegratorSettings integrator;
    TerrainInteractionSettings terrainInteractions;
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

struct XmlTagBlock
{
    std::unordered_map<std::string, std::string> attrs;
    std::string innerXml;
};

inline std::optional<XmlTagBlock> find_tag_block(const std::string& xml,
    const std::string& tag)
{
    const std::regex tagExpr("<\\s*" + tag + R"(\b([^>]*)>([\s\S]*?)<\s*/\s*)" + tag + R"(\s*>)",
        std::regex::icase);
    std::smatch tagMatch;
    if (!std::regex_search(xml, tagMatch, tagExpr)) {
        return std::nullopt;
    }

    XmlTagBlock block;
    block.innerXml = tagMatch[2].str();

    const std::regex attrExpr("([A-Za-z0-9_:\\-]+)\\s*=\\s*\"([^\"]*)\"");
    const std::string attrText = tagMatch[1].str();
    for (std::sregex_iterator it(attrText.begin(), attrText.end(), attrExpr), end; it != end; ++it) {
        block.attrs[to_lower_copy((*it)[1].str())] = (*it)[2].str();
    }
    return block;
}

inline std::vector<std::unordered_map<std::string, std::string>> find_all_tag_attributes(
    const std::string& xml,
    const std::string& tag)
{
    const std::regex tagExpr("<\\s*" + tag + R"(\b([^>]*)/\s*>)", std::regex::icase);
    const std::regex attrExpr("([A-Za-z0-9_:\\-]+)\\s*=\\s*\"([^\"]*)\"");
    std::vector<std::unordered_map<std::string, std::string>> tags;

    for (std::sregex_iterator tagIt(xml.begin(), xml.end(), tagExpr), end; tagIt != end; ++tagIt) {
        std::unordered_map<std::string, std::string> attrs;
        const std::string attrText = (*tagIt)[1].str();
        for (std::sregex_iterator attrIt(attrText.begin(), attrText.end(), attrExpr), attrEnd; attrIt != attrEnd; ++attrIt) {
            attrs[to_lower_copy((*attrIt)[1].str())] = (*attrIt)[2].str();
        }
        tags.push_back(std::move(attrs));
    }
    return tags;
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

inline RotationalIntegratorType parse_integrator_type(const std::string& text)
{
    const std::string value = to_lower_copy(trim_copy(text));
    if (value == "schweizer") return RotationalIntegratorType::Schweizer;
    if (value == "capobianco") return RotationalIntegratorType::Capobianco;
    throw std::runtime_error(
        "unsupported integrator type '" + text + "' (expected 'schweizer' or 'capobianco')");
}

// parse output format (txt vs. binary)
inline OutputFormat parse_output_format(const std::string& text)
{
    const std::string value = to_lower_copy(trim_copy(text));
    if (value == "text") {
        return OutputFormat::Text;
    }
    if (value == "binary") {
        return OutputFormat::Binary;
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
    if (value == "mass") { append_unique(OutputFieldId::Mass); return; }
    if (value == "shape_index" || value == "shapeindex") { append_unique(OutputFieldId::ShapeIndex); return; }

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

inline std::vector<OutputFieldId> default_output_fields()
{
    return std::vector<OutputFieldId>{
        OutputFieldId::Time,
        OutputFieldId::X,
        OutputFieldId::Y,
        OutputFieldId::Z,
        OutputFieldId::Qw,
        OutputFieldId::Qx,
        OutputFieldId::Qy,
        OutputFieldId::Qz,
        OutputFieldId::Vx,
        OutputFieldId::Vy,
        OutputFieldId::Vz,
        OutputFieldId::Wx,
        OutputFieldId::Wy,
        OutputFieldId::Wz
    };
}

inline RasterMetric parse_raster_metric(const std::string& text)
{
    const std::string value = to_lower_copy(trim_copy(text));
    if (value == "x") return RasterMetric::X;
    if (value == "y") return RasterMetric::Y;
    if (value == "z") return RasterMetric::Z;
    if (value == "vx") return RasterMetric::Vx;
    if (value == "vy") return RasterMetric::Vy;
    if (value == "vz") return RasterMetric::Vz;
    if (value == "vtot") return RasterMetric::Vtot;
    if (value == "wx") return RasterMetric::Wx;
    if (value == "wy") return RasterMetric::Wy;
    if (value == "wz") return RasterMetric::Wz;
    if (value == "wtot") return RasterMetric::Wtot;
    if (value == "ekin_trans" || value == "ekintrans") return RasterMetric::EkinTrans;
    if (value == "ekin_rot" || value == "ekinrot") return RasterMetric::EkinRot;
    if (value == "ekin_tot" || value == "ekin" || value == "ekintot") return RasterMetric::EkinTot;
    if (value == "jumph" || value == "jump_h") return RasterMetric::JumpH;
    throw std::runtime_error("unknown raster metric: " + text);
}

inline RasterStatistic parse_raster_statistic(const std::string& text)
{
    const std::string value = to_lower_copy(trim_copy(text));
    if (value == "count") return RasterStatistic::Count;
    if (value == "min") return RasterStatistic::Min;
    if (value == "max") return RasterStatistic::Max;
    if (value == "sum") return RasterStatistic::Sum;
    if (value == "mean") return RasterStatistic::Mean;
    if (value == "median") return RasterStatistic::Median;
    if (value == "percentile") return RasterStatistic::Percentile;
    throw std::runtime_error("unknown raster statistic: " + text);
}

inline RasterLayerSpec parse_raster_layer_config(
    const std::unordered_map<std::string, std::string>& rasterTag)
{
    RasterLayerSpec layer;
    const std::string metricText = require_attr(rasterTag, "metric", "raster");
    const std::string metricKey = to_lower_copy(trim_copy(metricText));

    if (metricKey == "count") {
        layer.statistic = RasterStatistic::Count;
        layer.metric = RasterMetric::None;

        if (const auto statisticAttr = optional_attr(rasterTag, "statistic")) {
            const RasterStatistic parsed = parse_raster_statistic(*statisticAttr);
            if (parsed != RasterStatistic::Count) {
                throw std::runtime_error("<raster metric=\"count\"/> must not define a statistic other than count");
            }
        }
    }
    else {
        layer.metric = parse_raster_metric(metricText);
        layer.statistic = parse_raster_statistic(require_attr(rasterTag, "statistic", "raster"));

        if (layer.statistic == RasterStatistic::Percentile) {
            const std::string pText = optional_attr(rasterTag, "p")
                .value_or(optional_attr(rasterTag, "percentile")
                    .value_or(std::string()));
            if (pText.empty()) {
                throw std::runtime_error("<raster statistic=\"percentile\"/> requires p=\"...\"");
            }
            layer.percentile = parse_number<Real>(pText, "raster percentile");
            if (!(*layer.percentile > Real(0) && *layer.percentile <= Real(100))) {
                throw std::runtime_error("raster percentile must be in (0, 100]");
            }
        }
    }

    if (const auto nameAttr = optional_attr(rasterTag, "name")) {
        layer.outputStem = trim_copy(*nameAttr);
    }

    return layer;
}

inline void parse_output_block(const XmlTagBlock& outputBlock,
    BatchSimulationConfig& cfg)
{
    const auto trajectoryTags = find_all_tag_attributes(outputBlock.innerXml, "trajectory");
    const auto rasterTags = find_all_tag_attributes(outputBlock.innerXml, "raster");
    const auto shapeTags = find_all_tag_attributes(outputBlock.innerXml, "shape");

    if (trajectoryTags.size() > 1) {
        throw std::runtime_error("<output> may contain at most one <trajectory/> tag");
    }
    if (shapeTags.size() > 1) {
        throw std::runtime_error("<output> may contain at most one <shape/> tag");
    }
    if (trajectoryTags.empty() && rasterTags.empty() && shapeTags.empty()) {
        throw std::runtime_error("<output> must contain at least one <trajectory/>, <raster/>, or <shape/> tag");
    }

    if (!trajectoryTags.empty()) {
        const auto& trajectoryTag = trajectoryTags.front();
        cfg.output.enabled = true;
        cfg.output.format = parse_output_format(require_attr(trajectoryTag, "format", "trajectory"));
        cfg.output.fields = parse_output_fields(require_attr(trajectoryTag, "metric", "trajectory"));
    }

    if (!rasterTags.empty()) {
        RasterOutputConfig rasterCfg;
        std::optional<Real> resolution;
        rasterCfg.layers.reserve(rasterTags.size());
        for (const auto& rasterTag : rasterTags) {
            const Real layerResolution = parse_number<Real>(
                require_attr(rasterTag, "resolution", "raster"),
                "raster resolution");
            if (!resolution.has_value()) {
                resolution = layerResolution;
                rasterCfg.resolution = layerResolution;
            }
            else if (std::abs(layerResolution - *resolution) > Real(1e-12)) {
                throw std::runtime_error("all <raster/> tags inside <output> must use the same resolution");
            }

            rasterCfg.layers.push_back(parse_raster_layer_config(rasterTag));
        }
        cfg.rasterOutput = std::move(rasterCfg);
    }

    if (!shapeTags.empty()) {
        const auto& shapeTag = shapeTags.front();
        const std::string metric = to_lower_copy(trim_copy(require_attr(shapeTag, "metric", "shape")));
        if (metric != "deposition") {
            throw std::runtime_error("<shape/> currently only supports metric=\"deposition\"");
        }

        DepositionPointOutputConfig depositionCfg;
        if (const auto stemAttr = optional_attr(shapeTag, "name")) {
            const std::string stem = trim_copy(*stemAttr);
            if (!stem.empty()) {
                depositionCfg.outputStem = stem;
            }
        }
        cfg.depositionPointOutput = std::move(depositionCfg);
    }
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
// X Y Z       : use absolute world-space release coordinates directly
// X Y Z dz    : ignore the given Z and instead place the rock dz above the
//               terrain (including optional roughness and initial orientation) at X/Y
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

    const auto terrainSelfClosingTag = find_tag_attributes(xml, "terrain");
    const auto terrainBlock = find_tag_block(xml, "terrain");
    std::optional<std::unordered_map<std::string, std::string>> terrainTag;
    if (terrainSelfClosingTag && terrainBlock) {
        throw std::runtime_error("config XML must contain only one <terrain> definition");
    }
    if (terrainSelfClosingTag) terrainTag = *terrainSelfClosingTag;
    if (terrainBlock) terrainTag = terrainBlock->attrs;
    const auto rocksTag = find_tag_attributes(xml, "rocks");
    const auto startsTag = find_tag_attributes(xml, "starts");
    const auto simulationTag = find_tag_attributes(xml, "simulation");
    const auto stopTag = find_tag_attributes(xml, "stop");
    const auto solverTag = find_tag_attributes(xml, "solver");
    const auto integratorTag = find_tag_attributes(xml, "integrator");
    const auto roughnessTag = find_tag_attributes(xml, "roughness");
    const auto fragmentationTag = find_tag_attributes(xml, "fragmentation");
    const auto outputTag = find_tag_block(xml, "output");
    const auto debugTag = find_tag_attributes(xml, "debug");

    if (!terrainTag || !rocksTag || !startsTag || !simulationTag || !stopTag || !solverTag || !integratorTag) {
        throw std::runtime_error("config XML must contain <terrain>, <rocks/>, <starts/>, <simulation/>, <stop/>, <integrator/>, and <solver/> tags");
    }

    cfg.demPath = resolve_relative(baseDir, require_attr(*terrainTag, "dem", "terrain"));
    bool foundCoulombFriction = false;
    if (terrainBlock) {
        bool foundRollingFriction = false;
        for (const auto& interactionTag : find_all_tag_attributes(
            terrainBlock->innerXml, "interaction")) {
            const std::string type = to_lower_copy(trim_copy(
                require_attr(interactionTag, "type", "interaction")));
            if (type == "rolling_friction" || type == "rolling-friction"
                || type == "rolling friction") {
                if (foundRollingFriction) {
                    throw std::runtime_error("<terrain> may contain only one rolling_friction interaction");
                }
                cfg.terrainInteractions.rollingFrictionCoefficient = parse_number<Real>(
                    require_attr(interactionTag, "mu_r", "rolling_friction interaction"),
                    "rolling friction mu_r");
                if (!(cfg.terrainInteractions.rollingFrictionCoefficient >= Real(0))) {
                    throw std::runtime_error("rolling friction mu_r must be non-negative");
                }
                foundRollingFriction = true;
            }
            else if (type == "coulomb_friction" || type == "coulomb-friction"
                || type == "coulomb friction") {
                if (foundCoulombFriction) {
                    throw std::runtime_error("<terrain> may contain only one coulomb_friction interaction");
                }

                const auto muAttr = optional_attr(interactionTag, "mu");
                const auto muMinAttr = optional_attr(interactionTag, "mumin");
                const auto muMaxAttr = optional_attr(interactionTag, "mumax");
                const auto kappaAttr = optional_attr(interactionTag, "kappa");
                const auto betaAttr = optional_attr(interactionTag, "beta");
                const bool hasAnySlippageAttribute = muMinAttr || muMaxAttr
                    || kappaAttr || betaAttr;

                auto& friction = cfg.terrainInteractions.coulombFriction;
                if (muAttr && !hasAnySlippageAttribute) {
                    friction.model = CoulombFrictionModel::Constant;
                    friction.mu = parse_number<Real>(*muAttr, "coulomb friction mu");
                    if (!std::isfinite(friction.mu) || friction.mu < Real(0)) {
                        throw std::runtime_error("coulomb friction mu must be finite and non-negative");
                    }
                }
                else if (!muAttr && muMinAttr && muMaxAttr && kappaAttr && betaAttr) {
                    friction.model = CoulombFrictionModel::SlippageDependent;
                    friction.muMin = parse_number<Real>(*muMinAttr, "coulomb friction mumin");
                    friction.muMax = parse_number<Real>(*muMaxAttr, "coulomb friction mumax");
                    friction.kappa = parse_number<Real>(*kappaAttr, "coulomb friction kappa");
                    friction.beta = parse_number<Real>(*betaAttr, "coulomb friction beta");
                    if (!std::isfinite(friction.muMin) || friction.muMin < Real(0)) {
                        throw std::runtime_error("coulomb friction mumin must be finite and non-negative");
                    }
                    if (!std::isfinite(friction.muMax) || friction.muMax < friction.muMin) {
                        throw std::runtime_error("coulomb friction mumax must be finite and at least mumin");
                    }
                    if (!std::isfinite(friction.kappa) || friction.kappa < Real(0)) {
                        throw std::runtime_error("coulomb friction kappa must be finite and non-negative");
                    }
                    if (!std::isfinite(friction.beta) || friction.beta < Real(0)) {
                        throw std::runtime_error("coulomb friction beta must be finite and non-negative");
                    }
                }
                else {
                    throw std::runtime_error(
                        "coulomb_friction requires either only mu, or all of mumin, mumax, kappa, and beta");
                }
                foundCoulombFriction = true;
            }
            else {
                throw std::runtime_error("unknown terrain interaction type: " + type);
            }
        }
    }
    cfg.rockListPath = resolve_relative(baseDir, require_attr(*rocksTag, "list", "rocks"));
    cfg.rockDensity = parse_number<Real>(require_attr(*rocksTag, "density", "rocks"), "rock density");
    cfg.randomOrientations = parse_number<std::uint32_t>(require_attr(*rocksTag, "random_orientations", "rocks"), "random_orientations");
    if (cfg.randomOrientations == 0) throw std::runtime_error("random_orientations must be at least 1");

    cfg.startPointsPath = resolve_relative(baseDir, require_attr(*startsTag, "list", "starts"));

    if (const auto legacyMuAttr = optional_attr(*simulationTag, "mu")) {
        if (foundCoulombFriction) {
            throw std::runtime_error(
                "define Coulomb friction either in <simulation mu=\"...\"> or as a terrain interaction, not both");
        }
        auto& friction = cfg.terrainInteractions.coulombFriction;
        friction.model = CoulombFrictionModel::Constant;
        friction.mu = parse_number<Real>(*legacyMuAttr, "mu");
        if (!std::isfinite(friction.mu) || friction.mu < Real(0)) {
            throw std::runtime_error("mu must be finite and non-negative");
        }
    }
    else if (!foundCoulombFriction) {
        throw std::runtime_error(
            "Coulomb friction must be defined by <simulation mu=\"...\"> or a coulomb_friction terrain interaction");
    }
    if (const auto enAttr = optional_attr(*simulationTag, "en")) {
        cfg.epsilonN = parse_number<Real>(*enAttr, "en");
        if (!(cfg.epsilonN >= Real(0) && cfg.epsilonN <= Real(1))) {
            throw std::runtime_error("en must be in [0, 1]");
        }
    }
    if (const auto etAttr = optional_attr(*simulationTag, "et")) {
        cfg.epsilonT = parse_number<Real>(*etAttr, "et");
        if (!(cfg.epsilonT >= Real(0) && cfg.epsilonT <= Real(1))) {
            throw std::runtime_error("et must be in [0, 1]");
        }
    }
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

    cfg.integrator.type = parse_integrator_type(require_attr(*integratorTag, "type", "integrator"));
    if (const auto value = optional_attr(*integratorTag, "atol")) {
        cfg.integrator.atol = parse_number<Real>(*value, "integrator atol");
    }
    if (const auto value = optional_attr(*integratorTag, "rtol")) {
        cfg.integrator.rtol = parse_number<Real>(*value, "integrator rtol");
    }
    if (const auto value = optional_attr(*integratorTag, "maxiter")) {
        cfg.integrator.maxiter = parse_number<std::uint32_t>(*value, "integrator maxiter");
    }
    if (cfg.integrator.atol < Real(0) || cfg.integrator.rtol < Real(0)) {
        throw std::runtime_error("integrator tolerances must be non-negative");
    }
    if (cfg.integrator.maxiter == 0) {
        throw std::runtime_error("integrator maxiter must be at least 1");
    }

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
        parse_output_block(*outputTag, cfg);
    }
    else {
        cfg.output.enabled = true;
        cfg.output.format = OutputFormat::Text;
        cfg.output.fields = default_output_fields();
    }

    if (fragmentationTag) {
        FragmentationConfig fragmentation;
        fragmentation.thresholdSpecificEnergy = parse_number<Real>(
            require_attr(*fragmentationTag, "threshold_specific_energy", "fragmentation"),
            "fragmentation threshold_specific_energy");
        fragmentation.halfSpecificEnergy = parse_number<Real>(
            require_attr(*fragmentationTag, "half_specific_energy", "fragmentation"),
            "fragmentation half_specific_energy");
        if (const auto exponent = optional_attr(*fragmentationTag, "exponent")) {
            fragmentation.exponent = parse_number<Real>(*exponent, "fragmentation exponent");
        }
        if (const auto maxFraction = optional_attr(*fragmentationTag, "max_removed_fraction")) {
            fragmentation.maxRemovedFraction = parse_number<Real>(
                *maxFraction, "fragmentation max_removed_fraction");
        }
        if (const auto shapeDirectory = optional_attr(*fragmentationTag, "shape_output_dir")) {
            fragmentation.shapeOutputDirectory = parse_config_path(trim_copy(*shapeDirectory));
        }

        if (!(fragmentation.thresholdSpecificEnergy >= Real(0))) {
            throw std::runtime_error("fragmentation threshold_specific_energy must be non-negative");
        }
        if (!(fragmentation.halfSpecificEnergy > fragmentation.thresholdSpecificEnergy)) {
            throw std::runtime_error("fragmentation half_specific_energy must be greater than threshold_specific_energy");
        }
        if (!(fragmentation.exponent > Real(0))) {
            throw std::runtime_error("fragmentation exponent must be positive");
        }
        if (!(fragmentation.maxRemovedFraction > Real(0)
            && fragmentation.maxRemovedFraction <= Real(0.5))) {
            throw std::runtime_error("fragmentation max_removed_fraction must be in (0, 0.5]");
        }
        if (fragmentation.shapeOutputDirectory.empty()
            || fragmentation.shapeOutputDirectory.is_absolute()
            || looks_like_unc_path(fragmentation.shapeOutputDirectory.string())) {
            throw std::runtime_error("fragmentation shape_output_dir must be a non-empty relative path");
        }
        cfg.fragmentation = std::move(fragmentation);
    }

    if (cfg.fragmentation && cfg.output.enabled) {
        append_output_field_alias(cfg.output.fields, "mass");
        append_output_field_alias(cfg.output.fields, "shape_index");
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
