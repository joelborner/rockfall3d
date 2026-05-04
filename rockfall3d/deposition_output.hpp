#pragma once

#include "math_alias.hpp"

#include <filesystem>
#include <optional>
#include <stdexcept>
#include <string>
#include <vector>

#ifdef WITH_GDAL
#include "gdal.h"
#include "ogr_api.h"
#include "ogr_srs_api.h"
#endif

namespace fs = std::filesystem;

struct DepositionPointOutputConfig
{
    std::string outputStem{ "deposition_points" };
};

struct DepositionPointRecord
{
    bool valid{ false };
    std::size_t runIndex{ 0 };
    std::size_t rockIndex{ 0 };
    std::string rockPath;
    std::size_t startIndex{ 0 };
    std::uint32_t orientationIndex{ 0 };
    Real startX{ Real(0) };
    Real startY{ Real(0) };
    Real startZ{ Real(0) };
    Real x{ Real(0) };
    Real y{ Real(0) };
    Real z{ Real(0) };
    std::string status;
    std::string stopCriterion;
};

inline void add_shapefile_field(OGRLayerH layer,
    const char* name,
    OGRFieldType type,
    int width = 0,
    int precision = 0)
{
#ifndef WITH_GDAL
    (void)layer;
    (void)name;
    (void)type;
    (void)width;
    (void)precision;
    throw std::runtime_error("deposition point output requires a build with GDAL");
#else
    OGRFieldDefnH field = OGR_Fld_Create(name, type);
    if (!field) {
        throw std::runtime_error(std::string("failed to create shapefile field definition: ") + name);
    }
    if (width > 0) {
        OGR_Fld_SetWidth(field, width);
    }
    if (precision > 0) {
        OGR_Fld_SetPrecision(field, precision);
    }
    if (OGR_L_CreateField(layer, field, TRUE) != OGRERR_NONE) {
        OGR_Fld_Destroy(field);
        throw std::runtime_error(std::string("failed to create shapefile field: ") + name);
    }
    OGR_Fld_Destroy(field);
#endif
}

// Writes deposition point records to a shapefile in the specified output directory, using the given projection for spatial reference.
// The shapefile will contain one point feature per valid record, with attributes containing trajectory information.
inline void write_deposition_points_shapefile(const fs::path& outputDir,
    const DepositionPointOutputConfig& config,
    const std::string& projectionWkt,
    const std::vector<DepositionPointRecord>& records)
{
#ifndef WITH_GDAL
    (void)outputDir;
    (void)config;
    (void)projectionWkt;
    (void)records;
    throw std::runtime_error("deposition point output requires a build with GDAL");
#else
    GDALAllRegister();

    const fs::path shpPath = outputDir / (config.outputStem + ".shp");
    GDALDriverH driver = GDALGetDriverByName("ESRI Shapefile");
    if (!driver) {
        throw std::runtime_error("ESRI Shapefile driver is not available");
    }

    if (fs::exists(shpPath)) {
        GDALDeleteDataset(driver, shpPath.string().c_str());
    }

    GDALDatasetH ds = GDALCreate(driver, shpPath.string().c_str(), 0, 0, 0, GDT_Unknown, nullptr);
    if (!ds) {
        throw std::runtime_error("failed to create deposition point shapefile: " + shpPath.string());
    }

    OGRSpatialReferenceH srs = nullptr;
    if (!projectionWkt.empty()) {
        srs = OSRNewSpatialReference(nullptr);
        char* wkt = const_cast<char*>(projectionWkt.c_str());
        if (OSRImportFromWkt(srs, &wkt) != OGRERR_NONE) {
            OSRDestroySpatialReference(srs);
            GDALClose(ds);
            throw std::runtime_error("failed to import DEM projection for deposition point shapefile");
        }
    }

    OGRLayerH layer = GDALDatasetCreateLayer(
        ds,
        config.outputStem.c_str(),
        srs,
        wkbPoint25D,
        nullptr);

    if (srs) {
        OSRDestroySpatialReference(srs);
    }
    if (!layer) {
        GDALClose(ds);
        throw std::runtime_error("failed to create deposition point layer: " + shpPath.string());
    }

    add_shapefile_field(layer, "run_idx", OFTInteger64, 18);
    add_shapefile_field(layer, "rock_idx", OFTInteger64, 18);
    add_shapefile_field(layer, "rock_path", OFTString, 254);
    add_shapefile_field(layer, "start_idx", OFTInteger64, 18);
    add_shapefile_field(layer, "orient_idx", OFTInteger64, 18);
    add_shapefile_field(layer, "start_x", OFTReal, 24, 12);
    add_shapefile_field(layer, "start_y", OFTReal, 24, 12);
    add_shapefile_field(layer, "start_z", OFTReal, 24, 12);
    add_shapefile_field(layer, "status", OFTString, 16);
    add_shapefile_field(layer, "stop_crit", OFTString, 32);

    OGRFeatureDefnH layerDefn = OGR_L_GetLayerDefn(layer);
    if (!layerDefn) {
        GDALClose(ds);
        throw std::runtime_error("failed to query deposition point layer definition");
    }

    const int runIdxField = OGR_FD_GetFieldIndex(layerDefn, "run_idx");
    const int rockIdxField = OGR_FD_GetFieldIndex(layerDefn, "rock_idx");
    const int rockPathField = OGR_FD_GetFieldIndex(layerDefn, "rock_path");
    const int startIdxField = OGR_FD_GetFieldIndex(layerDefn, "start_idx");
    const int orientIdxField = OGR_FD_GetFieldIndex(layerDefn, "orient_idx");
    const int startXField = OGR_FD_GetFieldIndex(layerDefn, "start_x");
    const int startYField = OGR_FD_GetFieldIndex(layerDefn, "start_y");
    const int startZField = OGR_FD_GetFieldIndex(layerDefn, "start_z");
    const int statusField = OGR_FD_GetFieldIndex(layerDefn, "status");
    const int stopCritField = OGR_FD_GetFieldIndex(layerDefn, "stop_crit");

    for (const auto& record : records) {
        if (!record.valid) {
            continue;
        }

        OGRFeatureH feature = OGR_F_Create(layerDefn);
        if (!feature) {
            GDALClose(ds);
            throw std::runtime_error("failed to allocate deposition point feature");
        }

        OGR_F_SetFieldInteger64(feature, runIdxField, static_cast<GIntBig>(record.runIndex + 1u));
        OGR_F_SetFieldInteger64(feature, rockIdxField, static_cast<GIntBig>(record.rockIndex + 1u));
        OGR_F_SetFieldString(feature, rockPathField, record.rockPath.c_str());
        OGR_F_SetFieldInteger64(feature, startIdxField, static_cast<GIntBig>(record.startIndex + 1u));
        OGR_F_SetFieldInteger64(feature, orientIdxField, static_cast<GIntBig>(record.orientationIndex + 1u));
        OGR_F_SetFieldDouble(feature, startXField, static_cast<double>(record.startX));
        OGR_F_SetFieldDouble(feature, startYField, static_cast<double>(record.startY));
        OGR_F_SetFieldDouble(feature, startZField, static_cast<double>(record.startZ));
        OGR_F_SetFieldString(feature, statusField, record.status.c_str());
        OGR_F_SetFieldString(feature, stopCritField, record.stopCriterion.c_str());

        OGRGeometryH point = OGR_G_CreateGeometry(wkbPoint25D);
        if (!point) {
            OGR_F_Destroy(feature);
            GDALClose(ds);
            throw std::runtime_error("failed to allocate deposition point geometry");
        }
        OGR_G_SetPoint(point, 0,
            static_cast<double>(record.x),
            static_cast<double>(record.y),
            static_cast<double>(record.z));
        OGR_F_SetGeometryDirectly(feature, point);

        if (OGR_L_CreateFeature(layer, feature) != OGRERR_NONE) {
            OGR_F_Destroy(feature);
            GDALClose(ds);
            throw std::runtime_error("failed to write deposition point feature");
        }
        OGR_F_Destroy(feature);
    }

    GDALClose(ds);
#endif
}
