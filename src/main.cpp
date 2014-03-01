
#include "main.hpp"
#include <complex>
#include <fstream>
#include <future>
#include <chrono>
#include <iostream>


const std::size_t IMAGE_SIZE = 4096;


int main(int argc, char** argv)
{
    Matrix2D matrix = compressMatrix(readRawMatrix("anti_histogram.txt"));
    auto metadata = getMetadata(matrix); //don't worry, it's just metadata...
    auto gradientMatrix = getGradient(matrix, metadata.first, metadata.second);
    writeGradient(gradientMatrix);

    return EXIT_SUCCESS;
    //return EXIT_FAILURE;
}



Matrix2D readRawMatrix(std::string filename)
{
    Matrix2D rawMatrix;

    std::cout << "Reading... ";
    std::cout.flush();

    std::ifstream fin;
    fin.open(filename, std::ofstream::out);

    for (std::size_t j = 0; j < IMAGE_SIZE; j++)
    {
        std::vector<unsigned long> row;
        for (std::size_t k = 0; k < IMAGE_SIZE; k++)
        {
            unsigned long value;
            fin >> value;
            row.push_back(value);
        }
        rawMatrix.push_back(row);
    }
    fin.close();

    std::cout << "done." << std::endl;
    std::cout.flush();

    return rawMatrix;
}



Matrix2D compressMatrix(const Matrix2D& raw)
{
    Matrix2D matrix;

    std::cout << "Compressing... ";
    std::cout.flush();

    for (std::size_t j = 0; j < IMAGE_SIZE; j += 8)
    {
        std::vector<unsigned long> row;
        for (std::size_t k = 0; k < IMAGE_SIZE; k += 8)
        {
            unsigned long sum = 0;
            for (std::size_t x = 0; x < 8; x++)
                for (std::size_t y = 0; y < 8; y++)
                    sum += raw[j + x][k + y];
            row.push_back(sum);
        }
        matrix.push_back(row);
    }

    std::cout << "done." << std::endl;
    std::cout.flush();

    return matrix;
}



std::pair<float, float> getMetadata(const Matrix2D& matrix)
{
    unsigned long temp = 0;
    for (const auto &row : matrix)
        for (const auto &cell : row)
            temp += cell;
    float sum = (float)temp;

    float maxDensity = 0;
    for (const auto &row : matrix)
        for (const auto &cell : row)
            if (cell / sum > maxDensity)
                maxDensity = cell / sum;

    return std::make_pair(sum, maxDensity);
}



Matrix2D getGradient(const Matrix2D& matrix, float sum, float maxDensity)
{
    Matrix2D gradientMatrix;

    for (const auto &row : matrix)
    {
        std::vector<unsigned long> array;
        for (const auto &cell : row)
        {
            float gradiant = getIntensity(cell / sum, 0, maxDensity * 0.008f, true);
            array.push_back((unsigned long)(gradiant * 255));
        }
        gradientMatrix.push_back(array);
    }

    return gradientMatrix;
}



void writeGradient(const Matrix2D& gradientMatrix)
{
    std::cout << "Writing... ";
    std::cout.flush();

    std::ofstream fout;
    fout.open(std::string("image.ppm"), std::ofstream::out);

    fout << "P2 512 512 255" << std::endl;
    for (const auto &row : gradientMatrix)
    {
        for (const auto &cell : row)
            fout << cell << " ";
        fout << std::endl;
    }
    fout.close();

    std::cout << "done." << std::endl;
    std::cout.flush();
}



float getIntensity(float value, float lowerRange, float upperRange, bool maxFull)
{
    float intensity;
    if (value < lowerRange)
        intensity = 0;
    else if (value > upperRange)
        intensity = maxFull ? 1 : 0;
    else
        intensity = (value - lowerRange) / (upperRange - lowerRange);

    return intensity;

}
