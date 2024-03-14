#pragma once

enum class BlendType { Linear, G1, G2, Bump, ERBS };
double blendFunction(double x, BlendType type);
