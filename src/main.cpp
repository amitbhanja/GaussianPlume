/*
 * uFVM-CPP:
 *      Parallel unstructured 3D fluid flow CFD package 2020
 *
 * Developed by:
 *      Fadl Moukalled
 *      Mhamad Mahdi Alloush
 *      Adam Fares
 *
 * File:
 *      main.cpp
*/

// Built-in libraries

// Code libraries
#include "realm.h"
#include <unordered_map>
#include <string_view>
#include <string>
#include <fstream>
#include <optional>
#include <algorithm>
#include <cmath>
#include <queue>
#include <future>
#include <iostream>

static const double molar_volume_ideal_gas = 24.45; // litres per mole

static double molar_mass = 0;
static double flow_rate = 0;
static ufvm::Point stack_outlet(0,0,0);

std::string trim(const std::string& str)
{
    size_t first = str.find_first_not_of(" \t");
    size_t last = str.find_last_not_of(" \t");
    return (first == std::string::npos || last == std::string::npos) ? "" : str.substr(first, last - first + 1);
}

std::optional<double> to_double(std::string_view str)
{
    try
    {
        size_t idx;
        double val = std::stod(std::string(str), &idx);
        if(idx == str.length())
            return val;
    }
    catch (const std::invalid_argument& ex)
    {
        std::cerr << "Invalid Argument Error " << ex.what() << std::endl;
    }
    catch (const std::out_of_range& ex)
    {
        std::cerr << "Out of Range Error " << ex.what() << std::endl;
    }

    return std::nullopt;
}

bool setStackOutlet(const std::string& str)
{
    //std::cout << "Stack Outlet str " << str << std::endl;
    size_t idX = 1;
    int currPoint = 0;
    while(idX < str.length())
    {
        size_t commaPos = str.find(',', idX);
        //std::cout << "Comma Pos " << commaPos << " idx " << idX << std::endl;
        if(commaPos == std::string::npos)
            commaPos = str.length() - 1;
        if(auto valueDouble = to_double(str.substr(idX, commaPos - idX)))
        {
            if(currPoint == 0)
                stack_outlet.x() = *valueDouble;
            else if (currPoint == 1)
                stack_outlet.y() = *valueDouble;
            else
                stack_outlet.z() = *valueDouble;
        }
        else
        {
            return false;
        }
        idX = commaPos + 1;
        ++currPoint;
    }
    return true;
}

bool readInputParams()
{
    const std::string inputParamsPath = "../Input/inputParams";
    std::ifstream inputFile(inputParamsPath);

    if(!inputFile.is_open())
    {
        std::cerr << "Unable to open file: " << inputParamsPath << std::endl;
        return false;
    }

    std::string line;
    while(std::getline(inputFile, line))
    {
        // Find the comment position
        size_t commentPos = line.find("//");
        if(commentPos != std::string::npos)
        {
            line = line.substr(0, commentPos);
        }

        if(line.empty())
            continue;

        size_t separatorPos = line.find('=');
        if(separatorPos != std::string::npos)
        {
            std::string key = trim(line.substr(0, separatorPos));
            std::string value = trim(line.substr(separatorPos + 1));
            if(key == "stack_outlet")
            {
                if(!setStackOutlet(value))
                    return false;
            }
            else if(key == "molar_mass")
            {
                if(auto valueDouble = to_double(value))
                {
                    molar_mass = *valueDouble;
                }
                else
                {
                    return false;
                }
            }
            else if(key == "flow_rate")
            {
                if(auto valueDouble = to_double(value))
                {
                    flow_rate = *valueDouble;
                    flow_rate /= 3600;
                }
                else
                {
                    return false;
                }
            }
        }
    }
    std::cout << "Molar Mass " << molar_mass << std::endl;
    std::cout << "Stack Outlet " << stack_outlet << std::endl;
    return true;
}

std::vector<std::pair<std::string, std::string>> process_chunk(const std::string& filePath, std::streampos startPos, size_t linesToRead)
{
    std::ifstream ifs(filePath);
    if(!ifs.is_open())
    {
        std::cerr << "Unable to open file for reading " << std::endl;
        throw std::runtime_error("Unable to open file for reading ");
    }

    ifs.seekg(startPos);
    std::vector<std::string> output;

    std::string line;
    ufvm::realm realm;
    realm.setStackOutlet(stack_outlet);

    std::cout << "LinestoRead " << linesToRead << std::endl;

    for(size_t i = 0; i < linesToRead && std::getline(ifs, line); ++i)
    {
        size_t strIdx = 0;
        int currCnt = 0;

        std::cout << "Reading " << line << std::endl;
        while(strIdx < line.length())
        {
            size_t commaPos = 0;
            std::string curr = "";
            if(currCnt == 0)
            {
                ufvm::Point point;
                commaPos = line.find(',', strIdx);
                curr = line.substr(strIdx, commaPos - strIdx);
                strIdx = commaPos + 1;
                point.x() = std::stod(curr);
                commaPos = line.find(',', strIdx);
                curr = line.substr(strIdx, commaPos - strIdx);
                strIdx = commaPos + 1;
                point.y() = std::stod(curr);
                commaPos = line.find(',', strIdx);
                curr = line.substr(strIdx, commaPos - strIdx);
                point.z() = std::stod(curr);
                strIdx = commaPos + 1;
                realm.addreceiver_point(point);
                //std::cout << "Point " << point << std::endl;
            }
            else if (currCnt == 1)
            {
                commaPos = line.find(',', strIdx);
                curr = line.substr(strIdx, commaPos - strIdx);
                double ppm = std::stod(curr);
                realm.addrecevier_measurement(((double)(ppm * 0.72)/(double)1000));
                strIdx = commaPos + 1;
                //std::cout << "Receiver Measurment " << ppm << std::endl;
            }
            else if (currCnt == 2)
            {
                commaPos = line.find(',', strIdx);
                curr = line.substr(strIdx, commaPos - strIdx);
                strIdx = commaPos + 1;
                double windVelocity = std::stod(curr);
                commaPos = line.find(',', strIdx);
                if(commaPos == std::string::npos)
                    commaPos = line.length();
                curr = line.substr(strIdx, commaPos - strIdx);
                double windAngle = std::stod(curr);
                strIdx = commaPos + 1;
                ufvm::Velocity v(windVelocity, windAngle);
                realm.addwind_velocity(v);
                //std::cout << "Velocity point " << v << std::endl;
            }
            ++currCnt;
        }
    }

    realm.setFlowRate(flow_rate);
    realm.run();
    return realm.getConcentrations();
}

std::vector<std::streampos> calculate_offsets(const std::string& filePath, size_t linesPerChunk, size_t headerlines = 1)
{
    std::ifstream ifs(filePath);
    if(!ifs.is_open())
    {
        throw std::runtime_error("Unable to open file");
    }

    std::vector<std::streampos> offsets;

    std::string header;
    for(size_t i = 0; i < headerlines; ++i)
    {
        std::getline(ifs, header);
    }

    offsets.push_back(ifs.tellg());

    size_t lineCount = 0;
    std::string line;

    while(std::getline(ifs, line))
    {
        if(++lineCount % linesPerChunk == 0)
        {
            std::cout << "Line Count " << lineCount << std::endl;
            offsets.push_back(ifs.tellg());
        }
    }

    ifs.close();
    return offsets;
}

int main(int argc, char *argv[])
{
    std::cout << "Read Input Params " << std::endl;
    if(!readInputParams())
    {
        std::cerr << "Unable to move foreward due to invalid input params" << std::endl;
        return -1;
    }

    std::cout << "Read Input Params End " << std::endl;
    const size_t linesPerChunk = 100;
    const size_t maxFutures = 5;
    const size_t headerLines = 1;
    const std::string filePath = "../Input/data.csv";

    std::cout << "Calculate offsets " << std::endl;
    std::vector<std::streampos> offsets = calculate_offsets(filePath, linesPerChunk, headerLines);

    size_t totalChunks = offsets.size();
    //process_chunk(filePath, offsets[0], linesPerChunk);


    std::queue<std::future<std::vector<std::pair<std::string, std::string> >>> futureQueue;
    std::vector<std::pair<std::string, std::string> > allLines;

    for(size_t i = 0; i < totalChunks; ++i)
    {
        std::cout << "Offset i " << offsets[i] << std::endl;
        auto lines = process_chunk(filePath, offsets[i], linesPerChunk);
        allLines.insert(allLines.end(), lines.begin(), lines.end());
//        futureQueue.push(std::async(std::launch::async, process_chunk, filePath, offsets[i], linesPerChunk));

//        if(futureQueue.size() >= maxFutures)
//        {
//            auto lines = futureQueue.front().get();
//            futureQueue.pop();
//            allLines.insert(allLines.end(), lines.begin(), lines.end());
//        }
    }
/*
    while(!futureQueue.empty()) {
        auto lines = futureQueue.front().get();
        futureQueue.pop();
        allLines.insert(allLines.end(), lines.begin(), lines.end());
    }
*/
    
    const std::string outfilePath = "../Output/forward_model_data.csv";
    std::ofstream ofs(outfilePath);

    if(ofs.is_open())
    {
        ofs << "Latitude,Longitude,Height,Stability,Sensor Measurement,Forward-Model-Concentration\n";

        for(const auto& pair : allLines)
        {
            ofs << pair.first << "," << pair.second << "\n";
        }

        ofs.close();
        std::cout << "Data Written to csv file " << std::endl;
    }
    else
    {
        std::cerr << "Unable to open file to write " << std::endl;
        return -1;
    }

    return 0;
}
