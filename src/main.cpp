
#include "main.hpp"
#include <complex>
#include <fstream>
#include <future>
#include <chrono>
#include <iostream>


const std::size_t IMAGE_SIZE = 4096;
const std::size_t POPULATION_SIZE = 30;
const std::size_t DNA_LENGTH = 32;
const std::size_t MUTATE_RATE = 3;



int main(int argc, char** argv)
{
    Matrix2D matrix = compressMatrix(readRawMatrix("anti_histogram.txt"));
    auto metadata = getMetadata(matrix); //don't worry, it's just metadata...
    auto gradientMatrix = getGradient(matrix, metadata.first, metadata.second);
    writeGradient(gradientMatrix, "image.ppm");

    auto mersenneTwister = getMersenneTwister();
    std::uniform_int_distribution<int> randomInt(1, 4);

    std::vector<std::string> population;
    for (int j = 0; j < POPULATION_SIZE; j++)
    {
        std::string str(DNA_LENGTH, '+');
        for (std::size_t k = 0; k < DNA_LENGTH; k++)
        {
            int type = randomInt(mersenneTwister);
            switch (type)
            {
                case 1 :
                    str[k] = '+';
                    break;
                case 2 :
                    str[k] = '-';
                    break;
                case 3 :
                    str[k] = '/';
                    break;
                case 4 :
                    str[k] = '*';
                    break;
            }
        }

        population.push_back(str);
    }

    cycle(population, gradientMatrix);

    return EXIT_SUCCESS;
    //return EXIT_FAILURE;
}



void cycle(std::vector<std::string>& population, const Matrix2D& image)
{
    std::vector<ScoredDNA> scores;
    for (const auto &individual : population)
    {
        auto value = score(individual, image);
        scores.push_back(std::make_pair(value, individual));
    }

    std::sort(scores.begin(), scores.end(), [](ScoredDNA a, ScoredDNA b) {
        return a.first < b.first;
    });

    //for (const auto &score : scores)
    //    std::cout << score.first << "," << score.second << std::endl;

    std::vector<std::string> newPopulation;
    for (int j = 0 ; j < POPULATION_SIZE / 2; j++)
        newPopulation.push_back(scores[j].second);

    for (int j = 0; j < POPULATION_SIZE / 2; j++)
        newPopulation.push_back(mutate(mate(scores[j].second, scores[j + 1].second)));

    //for (const auto &k : newPopulation)
    //    std::cout << k << std::endl;

    std::cout << newPopulation.size() << std::endl;

    std::cout << scores[0].first << ", " << scores[0].second << std::endl;
    std::ostringstream oss("");
    oss << "bestResult-" << scores[0].first << ".ppm";
    visualizeIndividual(scores[0].second, oss.str());
}



long score(const std::string& str, const Matrix2D& image)
{
    std::size_t error = 0;

    for (std::size_t x = 0; x < 512; x++)
    {
        for (std::size_t y = 0; y < 512; y++)
        {
            auto a = calculate(x, y, str);
            error += a > image[x][y] ? a - image[x][y] : image[x][y] - a;
        }
    }

    return error;
}



long calculate(std::size_t x, std::size_t y, const std::string& str)
{
    long a = x, b = y;

    bool state = true;
    for (std::size_t j = 0; j < str.size(); j++)
    {
        if (str[j] == '+')
            a = a + b;
        if (str[j] == '-')
            a = a < b ? 0 : a - b;
        if (str[j] == '*')
            a = a * b < 0 ? 2147483648 : a * b;
        if (str[j] == '/')
            a = b == 0 ? 2147483648 : a / b;

        b = state ? x : y;
        state = !state;
    }

    return a;
}



std::string mate(std::string a, std::string b)
{
    std::string child = a;

    for (int j = 0; j < DNA_LENGTH; j++)
    {
             if (a[j] == '+' && b[j] == '+')
            child[j] = '+';
        else if (a[j] == '+' && b[j] == '-')
            child[j] = '-';
        else if (a[j] == '+' && b[j] == '*')
            child[j] = '*';
        else if (a[j] == '+' && b[j] == '/')
            child[j] = '/';
        else if (a[j] == '-' && b[j] == '+')
            child[j] = '-';
        else if (a[j] == '-' && b[j] == '-')
            child[j] = '+';
        else if (a[j] == '-' && b[j] == '*')
            child[j] = '/';
        else if (a[j] == '-' && b[j] == '/')
            child[j] = '*';
        else if (a[j] == '*' && b[j] == '+')
            child[j] = '*';
        else if (a[j] == '*' && b[j] == '-')
            child[j] = '/';
        else if (a[j] == '*' && b[j] == '*')
            child[j] = '+';
        else if (a[j] == '*' && b[j] == '/')
            child[j] = '-';
        else if (a[j] == '/' && b[j] == '+')
            child[j] = '/';
        else if (a[j] == '/' && b[j] == '-')
            child[j] = '*';
        else if (a[j] == '/' && b[j] == '*')
            child[j] = '-';
        else if (a[j] == '/' && b[j] == '/')
            child[j] = '+';
    }

    return child;
}



std::string mutate(const std::string& str)
{
    auto mersenneTwister = getMersenneTwister();
    std::uniform_int_distribution<int> randomInt(1, 4);

    std::string mutated = str;
    for (std::size_t j = 0; j < MUTATE_RATE; j++)
    {
        int type = randomInt(mersenneTwister);
        switch (type)
        {
            case 1 :
                mutated[j] = '+';
                break;
            case 2 :
                mutated[j] = '-';
                break;
            case 3 :
                mutated[j] = '/';
                break;
            case 4 :
                mutated[j] = '*';
                break;
        }
    }

    return mutated;
}



std::mt19937 getMersenneTwister()
{
    std::array<int, std::mt19937::state_size> seed_data;
    std::random_device r;
    std::generate_n(seed_data.data(), seed_data.size(), std::ref(r));
    std::seed_seq seq(std::begin(seed_data), std::end(seed_data));
    std::mt19937 mersenneTwister(seq);
    return mersenneTwister;
}



void visualizeIndividual(const std::string& str, const std::string& filename)
{
    Matrix2D matrix;
    for (std::size_t j = 0; j < 512; j++)
    {
        std::vector<unsigned long> row;
        for (std::size_t k = 0; k < 512; k++)
            row.push_back(calculate(j, k, str));
        matrix.push_back(row);
    }

    float max = 0;
    for (const auto &row : matrix)
        for (const auto &cell : row)
            if (cell > max)
                max = cell;

    std::cout << "Writing result... ";
    std::cout.flush();

    std::ofstream fout;
    fout.open(filename, std::ofstream::out);

    fout << "P2 512 512 255" << std::endl;
    for (std::size_t j = 0; j < 512; j++)
    {
        for (std::size_t k = 0; k < 512; k++)
            fout << (int)(matrix[j][k] * 255 / max) << " ";
        fout << std::endl;
    }

    fout.close();

    std::cout << "done." << std::endl;
    std::cout.flush();
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
    std::cout << sum << ", " << maxDensity << std::endl;

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



void writeGradient(const Matrix2D& gradientMatrix, const std::string& filename)
{
    std::cout << "Writing... ";
    std::cout.flush();

    std::ofstream fout;
    fout.open(filename, std::ofstream::out);

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
