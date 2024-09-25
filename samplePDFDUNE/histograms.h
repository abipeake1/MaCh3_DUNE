#ifndef HISTOGRAMS_H
#define HISTOGRAMS_H

#include <iostream>
#include <vector>
#include <algorithm>  // For std::unique
#include <memory>     // For std::unique_ptr
#include "TH1D.h"
#include "TH2D.h"
#include "TH3D.h"     // ROOT histograms

// Function to generate linearly spaced values similar to numpy's linspace
template <typename T>
std::vector<T> linspace(T a, T b, size_t N) {
    T h = (b - a) / static_cast<T>(N);
    std::vector<T> xs(N);
    typename std::vector<T>::iterator x;
    T val;
    for (x = xs.begin(), val = a; x != xs.end(); ++x, val += h)
        *x = val;
    return xs;
}

// Function to extend a vector with linearly spaced values, ensuring no duplicates
void ExtendLinspace(std::vector<double>& in, double low, double high, size_t nbins) {
    for (auto x : linspace(low, high, nbins)) {
        in.push_back(x);
    }
    // Remove duplicate values in the vector
    auto x = std::unique(in.begin(), in.end());
    in.erase(x, in.end());
}

// Class to hold histograms and manage their binning
class Multidim_HistogramManager {
public:
    // Declare unique pointers to 2D and 3D histograms
    std::unique_ptr<TH3D> Abis3DHistogram;
    std::unique_ptr<TH3D> OA3DHistogram;
    std::unique_ptr<TH2D> OA2DHistogram;

    // Constructor that initializes the histograms
    Multidim_HistogramManager() {
        // Define bin ranges using linspace and ExtendLinspace
        std::vector<double> ELep_bins;
        ExtendLinspace(ELep_bins, 0, 20, 20);
        
        std::vector<double> theta_bins;
        ExtendLinspace(theta_bins, 0, 3, 10);
        ExtendLinspace(theta_bins, 3, 9, 3);
        ExtendLinspace(theta_bins, 9, 100, 2);
        
        std::vector<double> ENuReco_bins;
        ExtendLinspace(ENuReco_bins, 0, 4, 8);
        ExtendLinspace(ENuReco_bins, 4, 10, 6);
        
        std::vector<double> offaxis_position;
        ExtendLinspace(offaxis_position, -35, 0, 10);
        
        // Print bin sizes
        std::cout << "ELep_bins.size() = " << ELep_bins.size() << std::endl;
        std::cout << "theta_bins.size() = " << theta_bins.size() << std::endl;
        std::cout << "ENuReco_bins.size() = " << ENuReco_bins.size() << std::endl;

        // Initialize the 3D histogram for Abis
        Abis3DHistogram = std::make_unique<TH3D>("Abis3DHistogram", "", 
                                                 ELep_bins.size() - 1, ELep_bins.data(),
                                                 theta_bins.size() - 1, theta_bins.data(),
                                                 ENuReco_bins.size() - 1, ENuReco_bins.data());

        std::cout << "Number of bins in 3D Abis histogram = " << Abis3DHistogram->GetNcells() << std::endl;

        // Initialize the 3D histogram for Off-Axis
        OA3DHistogram = std::make_unique<TH3D>("OA3DHistogram", "", 
                                               ELep_bins.size() - 1, ELep_bins.data(),
                                               ENuReco_bins.size() - 1, ENuReco_bins.data(),
                                               offaxis_position.size() - 1, offaxis_position.data());
        
        std::cout << "Number of bins in  OA3DHistogram histogram = " <<  OA3DHistogram->GetNcells() << std::endl;

        // Initialize the 2D histogram for Off-Axis and ENuReco
        OA2DHistogram = std::make_unique<TH2D>("OA2DHistogram", "", 
                                               offaxis_position.size() - 1, offaxis_position.data(),
                                               ENuReco_bins.size() - 1, ENuReco_bins.data());

        std::cout << "Number of bins in  OA2DHistogram histogram = " <<  OA2DHistogram->GetNcells() << std::endl;
    }
};

#endif // HISTOGRAMS_H
