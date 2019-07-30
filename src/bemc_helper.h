#ifndef JETREADER_READER_BEMC_HELPER_H
#define JETREADER_READER_BEMC_HELPER_H

#include <vector>

namespace jetreader {

class BemcHelper {
public:
  BemcHelper();

  // provides access to a tower's absolute location with respect to the center
  // of the detector. Pseudorapidity (eta) is not corrected for the vertex
  // position along the beam line. The uncorrected pseudorapidity is useful for
  // QA, but should not be used for physics.
  double towerEta(unsigned tow_id);
  double towerPhi(unsigned tow_id);

  // calculates the pseudorapidity of the track with respect to the vertex -
  // this should be used for physics analyses.
  double vertexCorrectedEta(unsigned tow_id, double vz);

  // get module, eta position in module and phi position in module
  void hardwareLocation(unsigned soft_id, unsigned &module, unsigned &eta,
                        unsigned &phi);

private:
  // detector layout: 120 modules - 60 in phi x 2 in eta
  // each module is subivided into 40 towers - 2 in phi x 20 in eta
  // gives a total count of 4800 towers - 120 in phi x 40 in eta

  unsigned modules_;            // total number of module
  unsigned tow_per_module_phi_; // number of towers in module in phi direction
  unsigned tow_per_module_eta_; // number of towers in module in eta direction
  unsigned tow_per_module_;     // total number of towers per module
  unsigned towers_;             // total number of towers in detector

  // detector length constants
  double barrel_radius_;     // radius of barrel
  double half_width_module_; // half-width of module for the given radius

  // tower boundaries and centers
  std::vector<double> eta_bounds_; // edges for each tower (one side of TPC)
  std::vector<double>
      eta_center_; // the center of each tower in eta (one side of TPC)
  std::vector<double> phi_offset_; // phi offset from zero for first module
  std::vector<double> phi_module_step_width_; // width between modules in phi
  std::vector<double> phi_center_y_; // phi center for towers inside a module (in cm)
  std::vector<double>
      phi_center_; // phi center for towers inside a module (in degrees)

  // internal and external boundaries for the modules in eta
  double eta_min_;
  double eta_max_;
};

} // namespace jetreader

#endif // JETREADER_READER_BEMC_HELPER_H