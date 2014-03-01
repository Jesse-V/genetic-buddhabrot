
#ifndef MAIN
#define MAIN

#include <algorithm>
#include <vector>
#include <string>

typedef std::vector<std::vector<unsigned long>> Matrix2D;
typedef std::vector<std::string> DNA;

std::string mate(std::string a, std::string b);
std::string mutate(const std::string& str);
std::mt19937 getMersenneTwister();

Matrix2D readRawMatrix(std::string filename);
Matrix2D compressMatrix(const Matrix2D& matrix);
std::pair<float, float> getMetadata(const Matrix2D& matrix);
Matrix2D getGradient(const Matrix2D& matrix, float sum, float maxDensity);
void writeGradient(const Matrix2D& gradientMatrix);
float getIntensity(float value, float lowerRange, float upperRange, bool maxFull);

#endif
