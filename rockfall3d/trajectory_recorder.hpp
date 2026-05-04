#pragma once

#include "math_alias.hpp"
#include "rock.hpp"
#include "trajectory_format.hpp"

#include <cmath>
#include <cstdint>
#include <cstdio>
#include <filesystem>
#include <fstream>
#include <limits>
#include <optional>
#include <stdexcept>
#include <string>
#include <vector>

namespace fs = std::filesystem;

// Identifiers stored in the binary trajectory header.
struct TrajectoryMetadata
{
    std::size_t runIndex{};
    std::size_t rockIndex{};
    std::size_t startIndex{};
    std::uint32_t orientationIndex{};
    std::string caseName;
};

// Minimal state required by the recorder layer. Derived output variables
// are computed from this payload when a sample is actually written.
struct RecordedStepData
{
    Real t{ Real(0) };// time
    Eigen::Matrix<Real, 7, 1> q = (Eigen::Matrix<Real, 7, 1>() << 0, 0, 0, 1, 0, 0, 0).finished();// generalised coordinates
    Vec6 u = Vec6::Zero();// generalised velocities
    std::optional<Real> minGap;// jump height
};

struct RecordedTrajectoryLocation
{
    std::string storageFile;
    std::uint64_t byteOffset{ 0 };
    std::uint64_t sampleCount{ 0 };
};

class TrajectoryRecorder
{
public:
    virtual ~TrajectoryRecorder() = default;

    // Begin one trajectory stream. initialStep is the t=0 sample.
    virtual void beginTrajectory(const TrajectoryMetadata& meta,
        const Rock& rock,
        const RecordedStepData& initialStep) = 0;

    // Called for every accepted step; the recorder itself decides whether the
    // step falls onto a dump interval and therefore needs to be written.
    virtual void onAcceptedStep(const RecordedStepData& step) = 0;

    // Finalize the trajectory and return where it was stored.
    virtual RecordedTrajectoryLocation endTrajectory(const RecordedStepData& finalStep,
        double simSeconds) = 0;

    // Reset recorder state if a trajectory aborts before normal finalization.
    virtual void abortActiveTrajectory() noexcept = 0;
};

// Compute the value of a recorded field from the raw step data. This is used
// to populate the output fields during trajectory recording.
inline double recorded_field_value(OutputFieldId field,
    const Rock& rock,
    const RecordedStepData& step)
{
    switch (field) {
    case OutputFieldId::Time: return static_cast<double>(step.t);
    case OutputFieldId::X: return static_cast<double>(step.q(0));
    case OutputFieldId::Y: return static_cast<double>(step.q(1));
    case OutputFieldId::Z: return static_cast<double>(step.q(2));
    case OutputFieldId::Qw: return static_cast<double>(step.q(3));
    case OutputFieldId::Qx: return static_cast<double>(step.q(4));
    case OutputFieldId::Qy: return static_cast<double>(step.q(5));
    case OutputFieldId::Qz: return static_cast<double>(step.q(6));
    case OutputFieldId::Vx: return static_cast<double>(step.u(0));
    case OutputFieldId::Vy: return static_cast<double>(step.u(1));
    case OutputFieldId::Vz: return static_cast<double>(step.u(2));
    case OutputFieldId::Vtot: return static_cast<double>(step.u.head<3>().norm());
    case OutputFieldId::Wx: return static_cast<double>(step.u(3));
    case OutputFieldId::Wy: return static_cast<double>(step.u(4));
    case OutputFieldId::Wz: return static_cast<double>(step.u(5));
    case OutputFieldId::Wtot: return static_cast<double>(step.u.tail<3>().norm());
    case OutputFieldId::EkinTrans:
        return static_cast<double>(Real(0.5) * rock.mass() * step.u.head<3>().squaredNorm());
    case OutputFieldId::EkinRot:
        return static_cast<double>(Real(0.5) * step.u.tail<3>().dot(rock.inertiaK() * step.u.tail<3>()));
    case OutputFieldId::EkinTot:
        return static_cast<double>(
            Real(0.5) * rock.mass() * step.u.head<3>().squaredNorm()
            + Real(0.5) * step.u.tail<3>().dot(rock.inertiaK() * step.u.tail<3>()));
    case OutputFieldId::JumpH:
        if (step.minGap.has_value()) {
            return static_cast<double>(*step.minGap);
        }
        return std::numeric_limits<double>::quiet_NaN();
    }
    return std::numeric_limits<double>::quiet_NaN();
}

// Trajectory recorder implementation that writes binary files in a custom format. 
// The format is designed for efficient storage and fast sequential reading, but is not intended to be human-readable.
class BinaryRecorder final : public TrajectoryRecorder
{
public:
    BinaryRecorder(const fs::path& outputDir,
        std::uint32_t shardIndex,
        Real dumpStep,
        const OutputSettings& settings)
        : dumpStep_(dumpStep),
        fields_(settings.fields)
    {
        if (fields_.empty()) {
            throw std::runtime_error("binary output requires at least one selected trajectory metric");
        }
        if (!(dumpStep_ > Real(0))) {
            throw std::runtime_error("binary output requires a positive dump step");
        }

        layout_ = build_output_layout(fields_, recordSize_);

        fs::create_directories(outputDir);
        shardFileName_ = makeShardFileName(shardIndex);
        shardPath_ = outputDir / shardFileName_;
        out_.open(shardPath_, std::ios::binary | std::ios::trunc);
        if (!out_) {
            throw std::runtime_error("cannot open binary output file: " + shardPath_.string());
        }
        writeFileHeader();
    }

    void beginTrajectory(const TrajectoryMetadata& meta,
        const Rock& rock,
        const RecordedStepData& initialStep) override
    {
        if (active_) {
            throw std::runtime_error("binary recorder beginTrajectory called while another trajectory is still active");
        }

        currentRock_ = &rock;
        currentHeaderOffset_ = tellp_u64();
        currentSampleCount_ = 0;
        currentSampleCountPos_ = std::streampos(-1);
        currentSimSecondsPos_ = std::streampos(-1);
        lastWrittenTime_ = std::numeric_limits<Real>::quiet_NaN();
        nextDumpTime_ = initialStep.t + dumpStep_;
        active_ = true;

        writeTrajectoryHeader(meta);
        writeSample(initialStep);
    }

    void onAcceptedStep(const RecordedStepData& step) override
    {
        requireActive();

        const Real eps = Real(1e-12);
        if (step.t + eps < nextDumpTime_) {
            return;
        }

        writeSample(step);
        nextDumpTime_ = step.t + dumpStep_;
    }

	// Finalize the trajectory by writing any pending sample, updating the sample count and simulation seconds in the header, 
    // and returning the location of the recorded trajectory in the output file.
    RecordedTrajectoryLocation endTrajectory(const RecordedStepData& finalStep,
        double simSeconds) override
    {
        requireActive();

        const Real eps = Real(1e-12);
        if (!std::isfinite(lastWrittenTime_) || std::abs(finalStep.t - lastWrittenTime_) > eps) {
            writeSample(finalStep);
        }

        const std::streampos endPos = out_.tellp();
        out_.seekp(currentSampleCountPos_);
        writePod(currentSampleCount_);
        out_.seekp(currentSimSecondsPos_);
        writePod(simSeconds);
        out_.seekp(endPos);

        active_ = false;
        currentRock_ = nullptr;

        return RecordedTrajectoryLocation{
            shardFileName_,
            currentHeaderOffset_,
            currentSampleCount_
        };
    }

    const std::string& shardFileName() const noexcept
    {
        return shardFileName_;
    }

    void abortActiveTrajectory() noexcept override
    {
        active_ = false;
        currentRock_ = nullptr;
        currentSampleCount_ = 0;
        currentSampleCountPos_ = std::streampos(-1);
        currentSimSecondsPos_ = std::streampos(-1);
        nextDumpTime_ = Real(0);
        lastWrittenTime_ = std::numeric_limits<Real>::quiet_NaN();
    }

private:
    static std::string makeShardFileName(std::uint32_t shardIndex)
    {
        char buffer[32];
        std::snprintf(buffer, sizeof(buffer), "trajectories_%03u.rfbin", shardIndex + 1u);
        return std::string(buffer);
    }

    template <typename T>
    void writePod(const T& value)
    {
        out_.write(reinterpret_cast<const char*>(&value), sizeof(T));
        if (!out_) {
            throw std::runtime_error("failed while writing binary output file: " + shardPath_.string());
        }
    }

    std::uint64_t tellp_u64()
    {
        const auto pos = out_.tellp();
        if (pos < 0) {
            throw std::runtime_error("failed to query binary output position for: " + shardPath_.string());
        }
        return static_cast<std::uint64_t>(pos);
    }

	// Write the file header, which includes a magic number, version, field layout, and dump step. 
    // This header is written once per output file and describes how the trajectory records are structured.
    void writeFileHeader()
    {
        const char magic[8] = { 'R', 'F', '3', 'D', 'B', 'I', 'N', '1' };
        out_.write(magic, sizeof(magic));

        const std::uint32_t version = 1u;
        const std::uint32_t fieldCount = static_cast<std::uint32_t>(layout_.size());
        const std::uint32_t reserved = 0u;
        const double dumpStep = static_cast<double>(dumpStep_);

        writePod(version);
        writePod(fieldCount);
        writePod(recordSize_);
        writePod(reserved);
        writePod(dumpStep);

        for (const auto& field : layout_) {
            const std::uint16_t id = static_cast<std::uint16_t>(field.id);
            const std::uint8_t scalarType = static_cast<std::uint8_t>(field.scalarType);
            const std::uint8_t fieldReserved = 0u;
            writePod(id);
            writePod(scalarType);
            writePod(fieldReserved);
            writePod(field.offset);
        }
    }

	// Write the trajectory header, which includes identifiers for the trajectory and placeholders for the sample count and simulation seconds.
    void writeTrajectoryHeader(const TrajectoryMetadata& meta)
    {
        const std::uint32_t marker = 0x314A5254u; // "TRJ1"
        const std::uint32_t headerBytes = 44u;
        const std::uint64_t runIndex = static_cast<std::uint64_t>(meta.runIndex + 1u);
        const std::uint32_t rockIndex = static_cast<std::uint32_t>(meta.rockIndex + 1u);
        const std::uint32_t startIndex = static_cast<std::uint32_t>(meta.startIndex + 1u);
        const std::uint32_t orientationIndex = meta.orientationIndex + 1u;
        const std::uint64_t sampleCountPlaceholder = 0u;
        const double simSecondsPlaceholder = 0.0;

        writePod(marker);
        writePod(headerBytes);
        writePod(runIndex);
        writePod(rockIndex);
        writePod(startIndex);
        writePod(orientationIndex);
        currentSampleCountPos_ = out_.tellp();
        writePod(sampleCountPlaceholder);
        currentSimSecondsPos_ = out_.tellp();
        writePod(simSecondsPlaceholder);
    }

    static float toFloat32(Real value)
    {
        return static_cast<float>(value);
    }

	// Write a single sample of trajectory data according to the selected output fields and their layout.
    void writeSample(const RecordedStepData& step)
    {
        requireActive();
        if (currentRock_ == nullptr) {
            throw std::runtime_error("binary recorder has no active rock");
        }

        for (const auto& field : layout_) {
            const double value = recorded_field_value(field.id, *currentRock_, step);
            if (field.scalarType == OutputScalarType::Float64) {
                writePod(value);
            }
            else {
                const float value32 = static_cast<float>(value);
                writePod(value32);
            }
        }

        ++currentSampleCount_;
        lastWrittenTime_ = step.t;
    }

    void requireActive() const
    {
        if (!active_) {
            throw std::runtime_error("binary recorder received step data without an active trajectory");
        }
    }

private:
    fs::path shardPath_;
    std::string shardFileName_;
    std::ofstream out_;
    Real dumpStep_{ Real(0) };
    std::vector<OutputFieldId> fields_;
    std::vector<OutputFieldDescriptor> layout_;
    std::uint32_t recordSize_{ 0 };
    bool active_{ false };
    const Rock* currentRock_{ nullptr };
    std::uint64_t currentHeaderOffset_{ 0 };
    std::uint64_t currentSampleCount_{ 0 };
    std::streampos currentSampleCountPos_{};
    std::streampos currentSimSecondsPos_{};
    Real nextDumpTime_{ Real(0) };
    Real lastWrittenTime_{ std::numeric_limits<Real>::quiet_NaN() };
};

// Trajectory recorder implementation that writes human-readable text files.
class TextRecorder final : public TrajectoryRecorder
{
public:
    TextRecorder(const fs::path& outputDir,
        std::uint32_t shardIndex,
        Real dumpStep,
        const OutputSettings& settings)
        : dumpStep_(dumpStep),
        fields_(settings.fields)
    {
        if (fields_.empty()) {
            throw std::runtime_error("text output requires at least one selected output variable");
        }
        if (!(dumpStep_ > Real(0))) {
            throw std::runtime_error("text output requires a positive dump step");
        }

        fs::create_directories(outputDir);
        shardFileName_ = makeShardFileName(shardIndex);
        shardPath_ = outputDir / shardFileName_;
        out_.open(shardPath_);
        if (!out_) {
            throw std::runtime_error("cannot open text output file: " + shardPath_.string());
        }
        out_.setf(std::ios::fixed, std::ios::floatfield);
        out_.precision(9);
        writeFileHeader();
    }

	// Begin a new trajectory, resetting the internal state and writing the initial step.
    void beginTrajectory(const TrajectoryMetadata& meta,
        const Rock& rock,
        const RecordedStepData& initialStep) override
    {
        if (active_) {
            throw std::runtime_error("text recorder beginTrajectory called while another trajectory is still active");
        }

        currentRock_ = &rock;
        currentHeaderOffset_ = tellp_u64();
        currentSampleCount_ = 0;
        lastWrittenTime_ = std::numeric_limits<Real>::quiet_NaN();
        nextDumpTime_ = initialStep.t + dumpStep_;
        active_ = true;

        writeTrajectoryHeader(meta);
        writeSample(initialStep);
    }

	// Called for every accepted step. the recorder itself decides whether the step falls onto a dump interval and therefore needs to be written.
    void onAcceptedStep(const RecordedStepData& step) override
    {
        requireActive();
        const Real eps = Real(1e-12);
        if (step.t + eps < nextDumpTime_) {
            return;
        }
        writeSample(step);
        nextDumpTime_ = step.t + dumpStep_;
    }

	// Finalize the trajectory by writing any pending sample and returning the location of the recorded trajectory in the output file.
    RecordedTrajectoryLocation endTrajectory(const RecordedStepData& finalStep,
        double) override
    {
        requireActive();
        const Real eps = Real(1e-12);
        if (!std::isfinite(lastWrittenTime_) || std::abs(finalStep.t - lastWrittenTime_) > eps) {
            writeSample(finalStep);
        }

        active_ = false;
        currentRock_ = nullptr;
        out_.flush();

        return RecordedTrajectoryLocation{
            shardFileName_,
            currentHeaderOffset_,
            currentSampleCount_
        };
    }

    void abortActiveTrajectory() noexcept override
    {
        active_ = false;
        currentRock_ = nullptr;
        currentSampleCount_ = 0;
        nextDumpTime_ = Real(0);
        lastWrittenTime_ = std::numeric_limits<Real>::quiet_NaN();
    }

private:
    static std::string makeShardFileName(std::uint32_t shardIndex)
    {
        char buffer[32];
        std::snprintf(buffer, sizeof(buffer), "trajectories_%03u.txt", shardIndex + 1u);
        return std::string(buffer);
    }

    std::uint64_t tellp_u64()
    {
        const auto pos = out_.tellp();
        if (pos < 0) {
            throw std::runtime_error("failed to query text output position for: " + shardPath_.string());
        }
        return static_cast<std::uint64_t>(pos);
    }

    void writeFileHeader()
    {
        out_ << "# rockfall3d text trajectory shard\n";
        out_ << "# fields";
        for (const auto field : fields_) {
            out_ << ' ' << output_field_name(field);
        }
        out_ << '\n';
    }

    void writeTrajectoryHeader(const TrajectoryMetadata& meta)
    {
        out_ << "# trajectory"
            << " run_index=" << (meta.runIndex + 1u)
            << " rock_index=" << (meta.rockIndex + 1u)
            << " start_index=" << (meta.startIndex + 1u)
            << " orientation_index=" << (meta.orientationIndex + 1u)
            << " case_name=\"" << meta.caseName << "\"\n";
        if (!out_) {
            throw std::runtime_error("failed while writing text output file: " + shardPath_.string());
        }
    }

    void writeSample(const RecordedStepData& step)
    {
        requireActive();
        if (currentRock_ == nullptr) {
            throw std::runtime_error("text recorder has no active rock");
        }

        for (std::size_t i = 0; i < fields_.size(); ++i) {
            if (i > 0) {
                out_ << ' ';
            }
            out_ << recorded_field_value(fields_[i], *currentRock_, step);
        }
        out_ << '\n';
        if (!out_) {
            throw std::runtime_error("failed while writing text output file: " + shardPath_.string());
        }

        ++currentSampleCount_;
        lastWrittenTime_ = step.t;
    }

    void requireActive() const
    {
        if (!active_) {
            throw std::runtime_error("text recorder received step data without an active trajectory");
        }
    }

private:
    fs::path shardPath_;
    std::string shardFileName_;
    std::ofstream out_;
    Real dumpStep_{ Real(0) };
    std::vector<OutputFieldId> fields_;
    bool active_{ false };
    const Rock* currentRock_{ nullptr };
    std::uint64_t currentHeaderOffset_{ 0 };
    std::uint64_t currentSampleCount_{ 0 };
    Real nextDumpTime_{ Real(0) };
    Real lastWrittenTime_{ std::numeric_limits<Real>::quiet_NaN() };
};
