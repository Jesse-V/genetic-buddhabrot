
#ifndef MAIN
#define MAIN

#include <algorithm>
#include <vector>
#include <string>

typedef std::vector<std::vector<std::size_t>> Matrix2D;
typedef std::vector<std::string> DNA;
typedef std::pair<std::size_t, std::string> ScoredDNA;

void cycle(std::vector<std::string>& population, const Matrix2D& image);
long score(const std::string& str, const Matrix2D& image);
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
