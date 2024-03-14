#pragma once

#include <geometry.hh>

Geometry::Point2D param(size_t n, const Geometry::Point2D &);
bool loadCache(std::string filename);
bool saveCache(std::string filename);
