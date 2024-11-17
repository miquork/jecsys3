#include <iostream>
#include <vector>
#include <unordered_set>
#include <algorithm>
#include <queue>
#include <TFile.h>
#include <TTree.h>
#include <TH1D.h>
#include <TCanvas.h>
#include <cmath>
#include <set>
#include <fstream>

struct RunRange {
    int start_run;
    int end_run;
    float luminosity;
    long long start_time;
    long long end_time;
    int n_runs;
    int n_ls;

    RunRange(int start, int end, float lumi, long long start_t, long long end_t, int nruns, int nls) 
      : start_run(start), end_run(end), luminosity(lumi), start_time(start_t), end_time(end_t), n_runs(nruns), n_ls(nls) {}
};

struct Pair {
    int index1, index2;
    float combined_lumi;
    float min_lumi;
    int n_runs;
    int n_ls;

Pair(int i1, int i2, float lumi, float min, int nruns, int nls) : index1(i1), index2(i2), combined_lumi(lumi), min_lumi(min), n_runs(nruns), n_ls(nls) {}

    bool operator<(const Pair &other) const {
      return combined_lumi > other.combined_lumi; // Min-heap priority queue (by combined luminosity)
      //return min_lumi > other.min_lumi; // Min-heap priority queue (by minimum individual luminosity)
    }
};

// Function to check if two run ranges can be merged
bool canMerge(const RunRange &left, const RunRange &right, const std::unordered_set<int> &break_points, float max_lumi, long long max_time_diff) {
    // Check that the left and right runs do not cross a break point
    // Also ensure no overlap and check maximum luminosity and time difference criteria
    return (left.end_run < right.start_run) &&
           (break_points.find(left.end_run) == break_points.end()) &&
           (left.luminosity + right.luminosity <= max_lumi) &&
           (std::abs(right.start_time - left.end_time) <= max_time_diff);
}

void addEraBoundaryBreakPoints(std::vector<RunRange> &runs, const std::vector<std::pair<int, int>> &era_boundaries, std::unordered_set<int> &break_points, std::unordered_map<int, std::string> &break_point_sources) {
    for (const auto &boundary : era_boundaries) {
        int era_start = boundary.first;
        int era_end = boundary.second;

        // Find where era boundaries occur between existing runs
        for (size_t i = 0; i < runs.size() - 1; ++i) {
            if (runs[i].end_run < era_start && runs[i + 1].start_run > era_start) {
                break_points.insert(runs[i].end_run);
                break_point_sources[runs[i].end_run] = "Added by addEraBoundaryBreakPoints (Start Boundary)";
            }
            if (runs[i].end_run < era_end && runs[i + 1].start_run > era_end) {
                break_points.insert(runs[i].end_run);
                break_point_sources[runs[i].end_run] = "Added by addEraBoundaryBreakPoints (End Boundary)";
            }
        }
    }
}

std::vector<RunRange> clusterRuns(std::vector<RunRange> &runs, const std::unordered_set<int> &break_points, float max_lumi, long long max_time_diff) {
    // Priority queue to keep track of all possible merges by smallest combined luminosity (or smallest individual one)
    std::priority_queue<Pair> merge_queue;

    // Generate initial list of valid pairs (only consecutive pairs)
    for (size_t i = 0; i < runs.size() - 1; ++i) {
        if (canMerge(runs[i], runs[i + 1], break_points, max_lumi, max_time_diff)) {
	  merge_queue.emplace(i, i + 1, runs[i].luminosity + runs[i + 1].luminosity, min(runs[i].luminosity, runs[i +1].luminosity), runs[i].n_runs + runs[i+1].n_runs, runs[i].n_ls + runs[i+1].n_ls);
        }
    }

    // Iteratively merge clusters with the smallest combined luminosity
    while (!merge_queue.empty()) {
        // Get the pair with the smallest combined luminosity
        Pair to_merge = merge_queue.top();
        merge_queue.pop();

        int index1 = to_merge.index1;
        int index2 = to_merge.index2;

        // Merge the two clusters
        RunRange &left = runs[index1];
        RunRange &right = runs[index2];

        // Create the merged RunRange
        RunRange merged_range(left.start_run, right.end_run, left.luminosity + right.luminosity, left.start_time, right.end_time, left.n_runs + right.n_runs, left.n_ls + right.n_ls);

        // Replace the first cluster with the merged one
        runs[index1] = merged_range;

        // Remove the second cluster from the runs vector
        runs.erase(runs.begin() + index2);

        // Clear and regenerate the merge_queue with updated pairs (only consecutive pairs)
        std::priority_queue<Pair> new_merge_queue;
        for (size_t i = 0; i < runs.size() - 1; ++i) {
            if (canMerge(runs[i], runs[i + 1], break_points, max_lumi, max_time_diff)) {
	      new_merge_queue.emplace(i, i + 1, runs[i].luminosity + runs[i + 1].luminosity, min(runs[i].luminosity, runs[i+1].luminosity), runs[i].n_runs + runs[i+1].n_runs, runs[i].n_ls + runs[i+1].n_ls);
            }
        }
        merge_queue = std::move(new_merge_queue);
    }

    return runs;
}

void clusterRuns() {
    // Open the ROOT file and get the TTree
    TFile *file = new TFile("luminosity_data.root", "READ");
    TTree *tree = (TTree*)file->Get("LumiTree");

    // Variables to read from the tree
    int run, nls;
    float lumi;
    long long time;
    tree->SetBranchAddress("run", &run);
    tree->SetBranchAddress("lumi", &lumi);
    tree->SetBranchAddress("nls", &nls);
    //tree->SetBranchAddress("time", &time); // this is YYMMDDHHMM format
    tree->SetBranchAddress("minutes_since_2022", &time);

    // Load runs and luminosities into a vector of RunRange structures
    std::vector<RunRange> runs;
    for (int i = 0; i < tree->GetEntries(); i++) {
        tree->GetEntry(i);
        runs.emplace_back(run, run, lumi, time, time, 1, nls); // Each run starts as its own range
    }

    // Define era boundaries as pairs of (start_run, end_run)
    std::vector<std::pair<int, int>> era_boundaries = {
        // 2022
        {355100, 355793}, {355794, 357486}, {357487, 357733}, {357734, 358219}, {358220, 359021}, {359022, 360331},
        {360332, 362180}, {362181, 362349}, {362350, 362760},
        // 2023
        {366442, 367079}, {367080, 367515}, {367516, 367620}, {367621, 367763}, {367765, 369802}, {369803, 370602}, {370603, 370790},
        // 2024
        //{378971, 379411}, {379412, 380252}, {380253, 380947}, {380948, 381383}, {381384, 381943}, {381944, 383779}, // v1
	{378971, 378970}, {378971, 379411}, {379412, 380252}, {380253, 380947}, {380948, 381383}, {381384, 381943}, {381944, 383779},
        {383780, 385813}, {385814, 386408}, {386409, 386951} // v2
    };

    // Additional HCAL breaks
    // https://twiki.cern.ch/twiki/bin/viewauth/CMS/HcalRespCorrsTagsRun3
    // https://twiki.cern.ch/twiki/bin/view/CMS/HcalRespCorrsTags2011
    //eras.push_back(pair<int,string>(383195,"24_v2.0")); // HLT
    //eras.push_back(pair<int,string>(383219,"24_v2.1")); // HLT
    //eras.push_back(pair<int,string>(386401,"24_v3.0")); // HLT, eraH
    era_boundaries.push_back(pair<int,int>(382287,382287)); // 382262, 382298
    era_boundaries.push_back(pair<int,int>(383219,383219)); // 383175, 383247
    era_boundaries.push_back(pair<int,int>(386401,386401)); // 386319, 386478
    // => https://cms-talk.web.cern.ch/t/fast-track-validation-hlt-prompt-hcal-respcorrs-condition-update-from-hb-time-adjustment/25302/5
    era_boundaries.push_back(pair<int,int>(368775,368775)); // 368765, 368822
    
    // Additional ECAL intercalibration updates
    // https://cms-talk.web.cern.ch/c/ppd/alca/108
    // => https://cms-talk.web.cern.ch/t/gt-online-hlt-express-prompt-update-of-ecal-pedestals-conditions-run-384719-w34/46354/10
    era_boundaries.push_back(pair<int,int>(384719,384710)); // mid-24G!  384644, 384933
    // => https://cms-talk.web.cern.ch/t/l1-pre-announcement-of-ecal-intercalibration-update-at-l1-run-386025/57306
    era_boundaries.push_back(pair<int,int>(386025-1,386025-1)); // shift by 1 to get this right
    //era_boundaries.push_back(pair<int,int>(387576-1,387576-1));
    // => https://cms-talk.web.cern.ch/t/gt-online-hlt-express-prompt-update-of-ecal-pedestals-conditions-run-368782-w24/25407
    //era_boundaries.push_back(pair<int,int>(368919,368919)); // nah, matches 24D break

    // => https://cms-talk.web.cern.ch/t/full-track-validation-hlt-prompt-ecalintercalibconstants-conditions-from-runs-378981-379616/39588/4
    era_boundaries.push_back(pair<int,int>(379956,379956)); // mid-24C 379866, 379984
    
    // Define break points for the era boundaries and map their sources
    std::unordered_set<int> break_points;
    std::unordered_set<int> era_break_points;
    std::unordered_map<int, std::string> break_point_sources;
    for (const auto &boundary : era_boundaries) {
        break_points.insert(boundary.second);
        era_break_points.insert(boundary.second);
        break_point_sources[boundary.second] = "Original Era Boundary";
    }
    addEraBoundaryBreakPoints(runs, era_boundaries, break_points, break_point_sources);

    // Add break points for time gaps greater than 7 days (10080 minutes)
    long long max_time_gap = 10080;  // Maximum time gap in minutes (7 days)
    for (size_t i = 0; i < runs.size() - 1; ++i) {
        if (std::abs(runs[i + 1].start_time - runs[i].end_time) > max_time_gap) {
            break_points.insert(runs[i].end_run);
            break_point_sources[runs[i].end_run] = "Added due to 7-day gap requirement";
        }
    }

    // Define maximum luminosity threshold for merging
    float max_lumi = 1.5;  // Maximum luminosity per cluster (in fb^-1)

    // Perform initial clustering
    auto merged_runs = clusterRuns(runs, break_points, max_lumi, max_time_gap);

    // Perform second pass clustering, reusing the same function but with a reduced set of break points (only era boundaries)
    //merged_runs = clusterRuns(merged_runs, era_break_points, max_lumi, std::numeric_limits<long long>::max());

    // Print out the results
    std::ofstream outfile("run_clustering_output.txt");
    std::cout << "Merged Run Ranges ("<<merged_runs.size()<<" runs):" << std::endl;
    outfile << "Merged Run Ranges ("<<merged_runs.size()<<" runs):" << std::endl;
    for (const auto &range : merged_runs) {
        std::cout << "Runs " << range.start_run << " to " << range.end_run 
                  << ", Luminosity: " << range.luminosity << " fb^-1"
                  << ", Time Span: " << range.start_time << " to " << range.end_time << " (minutes)" << std::endl;
        outfile << "Runs " << range.start_run << " to " << range.end_run 
                << ", Luminosity: " << range.luminosity << " fb^-1"
                << ", Time Span: " << range.start_time << " to " << range.end_time << " (minutes)" << std::endl;
        if (break_points.find(range.end_run) != break_points.end()) {
            std::cout << "Break Point after Run: " << range.end_run << " (" << break_point_sources[range.end_run] << ")" << std::endl;
            outfile << "Break Point after Run: " << range.end_run << " (" << break_point_sources[range.end_run] << ")" << std::endl;
        }
    }
    outfile.close();

    // Prepare the arrays for TH1D histogram
    const int nrun = merged_runs.size();
    //double vrun[nrun + 1];
    //double vlum[nrun];
    vector<double> vrun;
    vector<double> vlum;

    for (int i = 0; i < nrun; ++i) {
      /*
        vrun[i] = merged_runs[i].start_run;
        vlum[i] = merged_runs[i].luminosity;
        if (break_points.find(merged_runs[i].end_run) != break_points.end()) {
            vrun[i + 1] = merged_runs[i].end_run + 1;
            vlum[i] = 0.0; // Insert a bin with zero luminosity to represent the era break
        }
      */
      vrun.push_back(merged_runs[i].start_run);
      vlum.push_back(merged_runs[i].luminosity);
      if (break_points.find(merged_runs[i].end_run) != break_points.end()) {
	vrun.push_back(merged_runs[i].end_run + 1);
	vlum.push_back(0.0); // Insert a bin with zero luminosity to represent the era break
      }
    }
    //vrun[nrun] = merged_runs.back().end_run + 1;
    vrun.push_back(merged_runs.back().end_run + 1);

    // Print arrays for verification
    int nruns = vrun.size()-1;
    std::ofstream arrays_outfile("run_histogram_data.txt");
    arrays_outfile << "double vrun["<<nruns+1<<"] = {";
    for (int i = 0; i <= nruns; ++i) {
        arrays_outfile << vrun[i] << (i < nruns ? ", " : "};\n");
    }

    arrays_outfile << "double vlum["<<nruns<<"] = {";
    for (int i = 0; i < nruns; ++i) {
        arrays_outfile << vlum[i] << (i < nruns - 1 ? ", " : "};\n");
    }
    arrays_outfile.close();

    // Print fibs and nibs
    std::ofstream fibs("fibs.txt");
    fibs << "[run1, run2] | name = [year][era]-nib[N]-fib[M] | lum(/fb) |  start and end time  | runs  | number of LS" << std::endl;
    int fib = 0;
    int nib = 1;
    int current_year = -1;
    std::string current_era;
    //std::unordered_set<int> nib_breaks = {370602, 380252}; // Example list of nib break points
    std::unordered_set<int> nib_breaks = {368822, 382298, 383247, 384933}; // 2023-2024 list

    // Define era boundaries as tuples of (start_run, end_run, "[year][era]")
    std::vector<std::tuple<int, int, std::string>> era_boundaries2 = {
        // 2022
      //{355100, 355793, "2022A"}, {355794, 357486, "2022B"}, /*{357487, 357733, "2022Dv1"}, {357734, 358219, "2022Dv2"}, {358220, 359021, "2022Dv3"},*/ // v1
      {355100, 355793, "2022B"}, {355794, 357486, "2022C"}, /*{357487, 357733, "2022Dv1"}, {357734, 358219, "2022Dv2"}, {358220, 359021, "2022Dv3"},*/ // v2
      {357487, 359021, "2022D"},
      {359022, 360331, "2022E"}, {360332, 362180, "2022F"}, {362181, 362349, "2022HI"}, {362350, 362760, "2022G"},
      // 2023
      {366442, 367079, "2023B"}, /*{367080, 367515, "2023Cv1"}, {367516, 367620, "2023Cv2"}, {367621, 367763, "2023Cv3"},*/
      {367080, 367763, "2023C"},
      {367765, 369802, "2023Cv4"}, /*{369803, 370602, "2023Dv1"}, {370603, 370790, "2023Dv2"},*/
      {369803, 370790, "2023D"},
      // 2024
      //{378971, 379411, "2024A"}, {379412, 380252, "2024B"}, {380253, 380947, "2024C"}, {380948, 381383, "2024D"}, {381384, 381943, "2024E"}, {381944, 383779, "2024F"}, // v1
      {378971, 379411, "2024B"}, {379412, 380252, "2024C"}, {380253, 380947, "2024D"}, {380948, 381383, "2024Ev1"}, {381384, 381943, "2024Ev2"}, {381944, 383779, "2024F"}, // v2
      {383780, 385813, "2024G"}, {385814, 386408, "2024H"}, {386409, 386951, "2024I"}
    };
    
    for (const auto &range : merged_runs) {
        // Determine the current era and year from era_boundaries
        for (const auto &boundary : era_boundaries2) {
            if (range.start_run >= std::get<0>(boundary) && range.end_run <= std::get<1>(boundary)) {
                if (current_year != stoi(std::get<2>(boundary).substr(0, 4)) || current_era != std::get<2>(boundary)) {
                    // Reset fib counter if era changes
                    current_year = stoi(std::get<2>(boundary).substr(0, 4));
                    current_era = std::get<2>(boundary);
                    fib = 0;
                    nib = 1;
                }
                break;
            }
        }

        // Increment nib if a nib break point is crossed
        if (nib_breaks.find(range.start_run) != nib_breaks.end()) {
            nib++;
        }

        // Increment fib number for non-zero luminosity runs
        if (range.luminosity > 0) {
            fib++;
        }

        // Print fib information
        fibs << "[" << range.start_run << ", " << range.end_run << "]  |  "
             << setw(7) << right << current_era << "-nib" << nib << "-fib" << setw(2) << left << fib << "  |  "
             << setw(9) << right << setprecision(4) << range.luminosity << " fb^-1  |  "
	     << setw(7) << range.start_time << " to " << setw(7) << range.end_time << "  |  "
	     << setw(3) << range.n_runs << "  |  " << range.n_ls  << std::endl;
    }
    fibs.close();


    // Create the histogram and draw it
    TCanvas *c1 = new TCanvas("c1", "Run Range Luminosity", 800, 600);
    TH1D *h = new TH1D("h", ";Run range;Run Range Lum (/fb)", nruns, &vrun[0]);
    for (int i = 0; i != nruns; ++i) {
        h->SetBinContent(i + 1, vlum[i]);
	if (vlum[i]!=0) h->SetBinError(i + 1, 0.);
	//if (vlum[i]!=0) h->SetBinError(i + 1, 0.015);
    }
    h->Draw("HE");
    c1->SaveAs("run_luminosity_vs_run.pdf");

    TCanvas *c2 = new TCanvas("c2", "Run Range Luminosity", 800, 600);
    TH1D *h2 = new TH1D("h2", ";Run Range Lum (/fb);Number of Clusters", 40, 0, 2);
    for (int i = 0; i != nruns; ++i) {
      if (vlum[i]!=0)
        h2->Fill(vlum[i]);
    }
    h2->Draw("HIST");
    c2->SaveAs("run_luminosity.pdf");

    // Fill up also simple lumi per run histogram
    int nrun1 = h->GetBinLowEdge(1);
    int nrun2 = h->GetBinLowEdge(h->GetNbinsX()+1);
    TH1D *hr = new TH1D("hr",";Run;Lum (/fb)",nrun2-nrun1,nrun1,nrun2);
    for (int i = 0; i < tree->GetEntries(); i++) {
      tree->GetEntry(i);
      int j = hr->FindBin(run);
      hr->SetBinContent(j, lumi);
    }

    // Do the same for cumulative luminosities
    TH1D *hc = (TH1D*)h->Clone("hc"); hc->Reset();
    double lumsum(0);
    for (int i = 1; i != hc->GetNbinsX()+1; ++i) {
      lumsum += h->GetBinContent(i);
      hc->SetBinContent(i, lumsum);
    }
    TH1D *hrc = (TH1D*)hr->Clone("hrc"); hrc->Reset();
    double lumsumr(0);
    for (int i = 1; i != hrc->GetNbinsX()+1; ++i) {
      lumsumr += hr->GetBinContent(i);
      hrc->SetBinContent(i, lumsumr);
    }

    // Then binning in terms of cumulative luminosity
    vector<double> v(1,0);
    for (int i = 1; i != hc->GetNbinsX()+1; ++i) {
      double cumlum = hc->GetBinContent(i);
      if (v[v.size()-1]!=cumlum) v.push_back(cumlum);
    }
    TH1D *hcb = new TH1D("hcb",";Cum. Lum. (/fb);Cum. Lum. (/fb);",v.size()-1,&v[0]);
    for (int i = 1; i != hc->GetNbinsX()+1; ++i) {
      double cumlum = hc->GetBinContent(i);
      int j = hcb->FindBin(cumlum)-1;
      hcb->SetBinContent(j, cumlum);
    }
    
    TCanvas *c3 = new TCanvas("c3", "Run Range Cumulative Luminosity", 800, 600);
    TH1D *h3 = new TH1D("h3", ";Run Range;Run Range Cum. Lum (/fb)", nruns, &vrun[0]);
    hrc->Draw("HIST");
    hc->Draw("HIST SAME");
    c3->SaveAs("run_cumulative_luminosity.pdf");

    TFile *fout = new TFile("clusterRuns.root","RECREATE");
    hr->Write("hlum",TObject::kOverwrite);
    h->Write("hlum2",TObject::kOverwrite);
    hrc->Write("hcumlum",TObject::kOverwrite);
    hc->Write("hcumlum2",TObject::kOverwrite);
    hcb->Write("hcumlumbins2",TObject::kOverwrite);
    fout->Write();
    fout->Close();
    
    // Clean up
    file->Close();
    delete file;
}
