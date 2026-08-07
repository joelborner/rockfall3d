#pragma once

#include <algorithm>
#include <cctype>
#include <cstdint>
#include <stdexcept>
#include <string>
#include <vector>

// trajectory output format
enum class OutputFormat : std::uint8_t
{
    Text = 0,
    Binary = 1
};

// trajectory output field IDs
enum class OutputFieldId : std::uint16_t
{
    Time = 1,
    X = 2,
    Y = 3,
    Z = 4,
    Qw = 5,
    Qx = 6,
    Qy = 7,
    Qz = 8,
    Vx = 9,
    Vy = 10,
    Vz = 11,
    Vtot = 12,
    Wx = 13,
    Wy = 14,
    Wz = 15,
    Wtot = 16,
    EkinTrans = 17,
    EkinRot = 18,
    EkinTot = 19,
    JumpH = 20,
    ScarDepth = 21,
    ScarLayerDepth = 22,
    DragFx = 23,
    DragFy = 24,
    DragFz = 25,
    DragFtot = 26,
    DragTx = 27,
    DragTy = 28,
    DragTz = 29,
    DragTtot = 30,
    ProjectedArea = 31,
    DragPowerTrans = 32,
    DragPowerRot = 33,
    DragPowerTotal = 34,
    DragWorkTrans = 35,
    DragWorkRot = 36,
    DragWorkTotal = 37,
    ContactSteps = 38,
    SolverIterations = 39
};

enum class OutputScalarType : std::uint8_t
{
    Float32 = 1,
    Float64 = 2
};

struct OutputFieldDescriptor
{
    OutputFieldId id{};
    OutputScalarType scalarType{ OutputScalarType::Float32 };
    std::uint32_t offset{};
};

struct OutputSettings
{
    bool enabled{ false };
    OutputFormat format{ OutputFormat::Text };
    std::vector<OutputFieldId> fields;
};

// turns a string to lowercase in-place
inline std::string output_to_lower(std::string value)
{
    std::transform(value.begin(), value.end(), value.begin(), [](unsigned char c) {
        return static_cast<char>(std::tolower(c));
        });
    return value;
}

// returns the string for an output field
inline std::string output_field_name(OutputFieldId field)
{
    switch (field) {
    case OutputFieldId::Time: return "time";
    case OutputFieldId::X: return "x";
    case OutputFieldId::Y: return "y";
    case OutputFieldId::Z: return "z";
    case OutputFieldId::Qw: return "qw";
    case OutputFieldId::Qx: return "qx";
    case OutputFieldId::Qy: return "qy";
    case OutputFieldId::Qz: return "qz";
    case OutputFieldId::Vx: return "vx";
    case OutputFieldId::Vy: return "vy";
    case OutputFieldId::Vz: return "vz";
    case OutputFieldId::Vtot: return "vtot";
    case OutputFieldId::Wx: return "wx";
    case OutputFieldId::Wy: return "wy";
    case OutputFieldId::Wz: return "wz";
    case OutputFieldId::Wtot: return "wtot";
    case OutputFieldId::EkinTrans: return "ekin_trans";
    case OutputFieldId::EkinRot: return "ekin_rot";
    case OutputFieldId::EkinTot: return "ekin_tot";
    case OutputFieldId::JumpH: return "jumph";
    case OutputFieldId::ScarDepth: return "scar_depth";
    case OutputFieldId::ScarLayerDepth: return "scar_layer_depth";
    case OutputFieldId::DragFx: return "drag_fx";
    case OutputFieldId::DragFy: return "drag_fy";
    case OutputFieldId::DragFz: return "drag_fz";
    case OutputFieldId::DragFtot: return "drag_ftot";
    case OutputFieldId::DragTx: return "drag_tx";
    case OutputFieldId::DragTy: return "drag_ty";
    case OutputFieldId::DragTz: return "drag_tz";
    case OutputFieldId::DragTtot: return "drag_ttot";
    case OutputFieldId::ProjectedArea: return "projected_area";
    case OutputFieldId::DragPowerTrans: return "drag_power_trans";
    case OutputFieldId::DragPowerRot: return "drag_power_rot";
    case OutputFieldId::DragPowerTotal: return "drag_power_total";
    case OutputFieldId::DragWorkTrans: return "drag_work_trans";
    case OutputFieldId::DragWorkRot: return "drag_work_rot";
    case OutputFieldId::DragWorkTotal: return "drag_work_total";
    case OutputFieldId::ContactSteps: return "contact_steps";
    case OutputFieldId::SolverIterations: return "solver_iterations";
    }
    return "unknown";
}

inline OutputScalarType output_field_scalar_type(OutputFieldId field)
{
    switch (field) {
    case OutputFieldId::X:
    case OutputFieldId::Y:
    case OutputFieldId::ContactSteps:
    case OutputFieldId::SolverIterations:
        return OutputScalarType::Float64;
    default:
        return OutputScalarType::Float32;
    }
}

inline std::uint32_t output_scalar_size(OutputScalarType type)
{
    switch (type) {
    case OutputScalarType::Float32: return 4u;
    case OutputScalarType::Float64: return 8u;
    }
    return 4u;
}

// builds a memory layout for the specified output fields, returning a list of field descriptors and the total record size in bytes
inline std::vector<OutputFieldDescriptor> build_output_layout(
    const std::vector<OutputFieldId>& fields,
    std::uint32_t& recordSize)
{
    std::vector<OutputFieldDescriptor> layout;
    layout.reserve(fields.size());
    std::uint32_t offset = 0;
    for (const auto field : fields) {
        const OutputScalarType scalarType = output_field_scalar_type(field);
        layout.push_back(OutputFieldDescriptor{ field, scalarType, offset });
        offset += output_scalar_size(scalarType);
    }
    recordSize = offset;
    return layout;
}

// returns true if the specified output fields include any that require the minimum gap to be tracked during simulation (currently only JumpH)
inline bool output_fields_require_min_gap(const std::vector<OutputFieldId>& fields)
{
    return std::find(fields.begin(), fields.end(), OutputFieldId::JumpH) != fields.end();
}
