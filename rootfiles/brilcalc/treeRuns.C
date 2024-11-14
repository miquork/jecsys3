#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <unordered_map>
#include <stdexcept>
#include <ctime>
#include <TFile.h>
#include <TTree.h>

// Function to convert date and time strings to an integer in YYMMDDHHMM format
long long convertToTimeInt(const std::string &dateStr, const std::string &timeStr) {
    struct tm tm = {};
    strptime((dateStr + " " + timeStr).c_str(), "%m/%d/%y %H:%M:%S", &tm);

    // Format as YYMMDDHHMM
    char buffer[11];
    strftime(buffer, sizeof(buffer), "%y%m%d%H%M", &tm);
    return std::stoll(buffer);
}

// Function to convert date and time to minutes since 2022-01-01 00:00
long long getMinutesSince2022(const std::string &dateStr, const std::string &timeStr) {
    struct tm tm = {};
    strptime((dateStr + " " + timeStr).c_str(), "%m/%d/%y %H:%M:%S", &tm);

    // Set the start of 2022 in time format
    struct tm tm_2022 = {};
    strptime("01/01/22 00:00:00", "%m/%d/%y %H:%M:%S", &tm_2022);
    
    // Compute difference in seconds and convert to minutes
    time_t t_run = mktime(&tm);
    time_t t_start = mktime(&tm_2022);
    return (t_run - t_start) / 60;
}

// Function to load data from a file and fill vectors
void loadData(const std::string &filename, std::vector<int> &years, std::vector<int> &fills, 
              std::vector<int> &runs, std::vector<float> &lumis, std::vector<long long> &times,
              std::vector<long long> &minutes_since_2022, std::unordered_map<int, float> &fill_lumi_map) {
    std::ifstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Failed to open file: " << filename << std::endl;
        return;
    }

    std::string line;
    // Skip four header rows
    for (int i = 0; i < 4; ++i) std::getline(file, line);

    while (std::getline(file, line)) {
        // Use sscanf to check for valid 'run:fill' format at the beginning
        int run, fill;
        char dateStr[20], timeStr[20];
        float lumi;
        
        // Parse only lines matching the data format (run:fill, date, time, and recorded lumi)
        if (sscanf(line.c_str(), " | %d:%d | %10s %8s | %*d | %*d | %*f | %f |", &run, &fill, dateStr, timeStr, &lumi) != 5) {
            std::cerr << "Skipping line due to formatting: " << line << std::endl;
            continue;
        }

        // Convert date and time
        long long timeVal = convertToTimeInt(dateStr, timeStr);
        long long minutesSince2022 = getMinutesSince2022(dateStr, timeStr);

        // Determine year based on run
        int year = (run < 363380) ? 2022 : (run < 376370) ? 2023 : 2024;

        // Store parsed values
        years.push_back(year);
        fills.push_back(fill);
        runs.push_back(run);
        lumis.push_back(lumi);
        times.push_back(timeVal);
        minutes_since_2022.push_back(minutesSince2022);

        // Accumulate luminosity per fill
        fill_lumi_map[fill] += lumi;
    }
    file.close();
}

void treeRuns() {
    // Vectors to store data
    std::vector<int> years, fills, runs;
    std::vector<float> lumis;
    std::vector<long long> times, minutes_since_2022;
    std::unordered_map<int, float> fill_lumi_map;  // Cumulative luminosity per fill

    // Load data from each file (assuming files are in the same directory)
    loadData("2022.txt", years, fills, runs, lumis, times, minutes_since_2022, fill_lumi_map);
    loadData("2023.txt", years, fills, runs, lumis, times, minutes_since_2022, fill_lumi_map);
    loadData("2024.txt", years, fills, runs, lumis, times, minutes_since_2022, fill_lumi_map);

    // Create ROOT file and TTree
    TFile *file = new TFile("luminosity_data.root", "RECREATE");
    TTree *tree = new TTree("LumiTree", "Lumi data per run and fill");

    // Variables to hold data for each entry
    int year, fill, run;
    float lumi, fill_lumi, cumulative_lumi = 0.0;
    long long time, mins_since_2022;

    // Set branches
    tree->Branch("year", &year, "year/I");
    tree->Branch("fill", &fill, "fill/I");
    tree->Branch("run", &run, "run/I");
    tree->Branch("lumi", &lumi, "lumi/F");
    tree->Branch("time", &time, "time/L");
    tree->Branch("fill_lumi", &fill_lumi, "fill_lumi/F");
    tree->Branch("cumulative_lumi", &cumulative_lumi, "cumulative_lumi/F");
    tree->Branch("minutes_since_2022", &mins_since_2022, "minutes_since_2022/L");

    // Fill TTree with data
    std::unordered_map<int, bool> fill_recorded;  // Track if a fill's lumi has been recorded
    for (size_t i = 0; i < years.size(); ++i) {
        year = years[i];
        fill = fills[i];
        run = runs[i];
        lumi = lumis[i];
        time = times[i];
        mins_since_2022 = minutes_since_2022[i];

        // Only set `fill_lumi` once per fill to avoid duplication
        if (!fill_recorded[fill]) {
            fill_lumi = fill_lumi_map[fill];
            fill_recorded[fill] = true;
        } else {
            fill_lumi = 0.0;  // Set to zero for duplicate entries of the same fill
        }

        // Update cumulative luminosity
        cumulative_lumi += lumi;

        // Fill the tree with the updated values
        tree->Fill();
    }

    // Write TTree to file and close
    tree->Write();
    file->Close();

    std::cout << "Data has been written to luminosity_data.root as a TTree." << std::endl;
}
