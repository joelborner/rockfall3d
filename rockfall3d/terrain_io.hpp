#pragma once
#include "terrain.hpp"
#include <fstream>
#include <sstream>
#include <vector>
#include <cctype>
#include <limits>
#include <stdexcept>

// Load terrain in ASCII grid format
inline Terrain load_asc(const std::string& path)
{
    std::ifstream in(path);
    if (!in) throw std::runtime_error("cannot open .asc file: " + path);

    std::size_t ncols = 0, nrows = 0;
    Real xll = 0, yll = 0, cell = 0, nodata = -9999.0;

    for (int i = 0; i < 6; ++i)
    {
        std::string key;
        in >> key;
        if (key == "ncols")         in >> ncols;
        else if (key == "nrows")    in >> nrows;
        else if (key == "xllcorner" || key == "xllcenter") in >> xll;
        else if (key == "yllcorner" || key == "yllcenter") in >> yll;
        else if (key == "cellsize") in >> cell;
        else if (key == "NODATA_value") in >> nodata;
    }
    if (ncols == 0 || nrows == 0)
        throw std::runtime_error("bad ASC header in " + path);

    Terrain dem;
    dem.xll = xll;
    dem.yll = yll;
    dem.cell = cell;
    dem.invCell = Real(1) / cell;
    dem.nx = static_cast<unsigned>(ncols);
    dem.ny = static_cast<unsigned>(nrows);
    dem.xmax = dem.xll + Real(dem.nx - 1) * dem.cell;
    dem.ymax = dem.yll + Real(dem.ny - 1) * dem.cell;
    dem.z.resize(nrows, ncols);

    for (std::size_t r = 0; r < nrows; ++r)
        for (std::size_t c = 0; c < ncols; ++c)
        {
            Real v;
            in >> v;
            dem.z(r, c) = (v == nodata)
                ? std::numeric_limits<Real>::quiet_NaN()
                : v;
        }

    dem.gradx = (dem.z.block(0, 1, nrows, ncols - 1) -
        dem.z.block(0, 0, nrows, ncols - 1)) / cell;
    dem.grady = (dem.z.block(1, 0, nrows - 1, ncols) -
        dem.z.block(0, 0, nrows - 1, ncols)) / cell;
    return dem;
}

#ifdef WITH_GDAL
#include "gdal.h"
#include <iostream>

// load terrain in GEOTIFF format with GDAL
inline Terrain load_tif(const std::string& path)
{
    GDALAllRegister();
    GDALDatasetH ds = GDALOpen(path.c_str(), GA_ReadOnly);
    if (!ds) throw std::runtime_error("GDALOpen failed for terrain file: " + path);

    const int nx = GDALGetRasterXSize(ds);
    const int ny = GDALGetRasterYSize(ds);

    double gt[6]{};
    if (GDALGetGeoTransform(ds, gt) != CE_None)
        throw std::runtime_error("no geotransform in terrain file: " + path);

    const bool east_rev = gt[1] < 0;   // pixel width sign
    const bool north_up = gt[5] < 0;   // pixel height sign

    // allocate final matrix: row-major, south-west origin
    MatX z(ny, nx);                    // row = south->north

    GDALRasterBandH band = GDALGetRasterBand(ds, 1);

    std::vector<double> scan(nx);      // one scan-line buffer

    // nodata handling
    int hasND = 0;
    const double ndv = GDALGetRasterNoDataValue(band, &hasND);
    const double NaN = std::numeric_limits<double>::quiet_NaN();

    for (int jSrc = 0; jSrc < ny; ++jSrc)
    {
        // read source row jSrc (0 = top of file)
        if (GDALRasterIO(band, GF_Read, 0, jSrc, nx, 1,
            scan.data(), nx, 1,
            GDT_Float64, 0, 0) != CE_None)
            throw std::runtime_error("GDALRasterIO row failed for terrain file: " + path);

        // destination row (south-west origin)
        int jDst = north_up ? (ny - 1 - jSrc) : jSrc;

        // copy/flip columns
        for (int iSrc = 0; iSrc < nx; ++iSrc)
        {
            int iDst = east_rev ? (nx - 1 - iSrc) : iSrc;
            double v = scan[iSrc];
            if (hasND && v == ndv) v = NaN;
            z(jDst, iDst) = static_cast<Real>(v);
        }
    }

    // build terrain
    Terrain dem;
    dem.cell = static_cast<Real>(std::abs(gt[1]));
    dem.invCell = Real(1) / dem.cell;
    dem.nx = nx; dem.ny = ny; dem.z = std::move(z);
    const char* projection = GDALGetProjectionRef(ds);
    if (projection != nullptr) {
        dem.projectionWkt = projection;
    }

    dem.xll = east_rev ? static_cast<Real>(gt[0] + (nx - 1) * gt[1])
        : static_cast<Real>(gt[0]);
    dem.yll = north_up ? static_cast<Real>(gt[3] + ny * gt[5])
        : static_cast<Real>(gt[3]);
    dem.xmax = dem.xll + Real(dem.nx - 1) * dem.cell;
    dem.ymax = dem.yll + Real(dem.ny - 1) * dem.cell;

    dem.gradx = (dem.z.block(0, 1, ny, nx - 1) -
        dem.z.block(0, 0, ny, nx - 1)) / dem.cell;
    dem.grady = (dem.z.block(1, 0, ny - 1, nx) -
        dem.z.block(0, 0, ny - 1, nx)) / dem.cell;

    GDALClose(ds);
    return dem;
}
#endif

// dispatcher to load terrain
inline Terrain load_terrain(const std::string& src)
{
    std::size_t dot = src.find_last_of('.');
    if (dot == std::string::npos)
        throw std::runtime_error("file has no extension");

    std::string ext = src.substr(dot + 1);
    for (char& ch : ext) ch = static_cast<char>(std::tolower(ch));

    if (ext == "asc")
        return load_asc(src);

#ifdef WITH_GDAL
    if (ext == "tif" || ext == "tiff")
        return load_tif(src);
#endif

    throw std::runtime_error("unsupported terrain file type: " + ext);
}

