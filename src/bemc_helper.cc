#include "bemc_helper.h"

#include <exception>
#include <cmath>

namespace jetreader {

BemcHelper::BemcHelper()
    : modules_(120), tow_per_module_phi_(2), tow_per_module_eta_(20),
      eta_min_(0.0035), eta_max_(0.984), barrel_radius_(225.405),
      half_width_module_(11.174) {

  tow_per_module_ = tow_per_module_eta_ * tow_per_module_phi_;
  towers_ = tow_per_module_ * modules_;

  eta_bounds_.resize(tow_per_module_eta_ + 1);
  for (int i = 0; i <= tow_per_module_eta_; ++i)
    eta_bounds_[i] = 0.05 * i;
  eta_bounds_[0] = eta_min_;
  eta_bounds_[tow_per_module_eta_] = eta_max_;

  eta_center_.resize(tow_per_module_eta_);
  for (int i = 0; i < tow_per_module_eta_; ++i)
    eta_center_[i] = (eta_bounds_[i] + eta_bounds_[i + 1]) / 2.0;

  const double pi = M_PI;
  phi_offset_ = {72.0 / 180.0 * pi, 108.0 / 180.0 * pi};

  phi_module_step_width_ = {-pi * 2.0 / (modules_ / 2.0),
                            pi * 2.0 / (modules_ / 2.0)};

  phi_center_y_ = {-1.0 * half_width_module_ / 2.0, half_width_module_ / 2.0};

  phi_center_ = {atan2(phi_center_y_[0], barrel_radius_),
                 atan2(phi_center_y_[1], barrel_radius_)};
}

double BemcHelper::towerEta(unsigned tow_id) {
  unsigned module_idx, module_eta, module_phi;
  hardwareLocation(tow_id, module_idx, module_eta, module_phi);

  double eta = eta_center_[module_eta];
  if (module_idx >= modules_ / 2.0)
    eta *= -1.0;
  return eta;
}

double BemcHelper::towerPhi(unsigned tow_id) {
  unsigned module_idx, module_eta, module_phi;
  hardwareLocation(tow_id, module_idx, module_eta, module_phi);

  int detector_side, module_on_side;
  double phi = phi_center_[module_phi];
  if (module_idx < modules_ / 2) {
    phi *= -1.0;
    detector_side = 0;
    module_on_side = module_idx;
  } else {
    detector_side = 1;
    module_on_side = module_idx - modules_ / 2;
  }
  phi += phi_offset_[detector_side];
  phi += phi_module_step_width_[detector_side] * module_on_side;

  double pi = M_PI;
  while (phi < -1.0 * pi)
    phi += 2.0 * pi;
  while (phi >= pi)
    phi -= 2.0 * pi;

  return phi;
}

double BemcHelper::vertexCorrectedEta(unsigned tow_id, double vz) {
  double tower_eta = towerEta(tow_id);
  double tower_theta = 2.0 * atan(exp(-tower_eta));
  double z = 0.0;
  if (tower_eta != 0.0)
    z = barrel_radius_ / tan(tower_theta);
  double z_diff = z - vz;
  double theta_corr = atan2(barrel_radius_, z_diff);
  double eta_corr = -log(tan(theta_corr / 2.0));
  return eta_corr;
}

void BemcHelper::hardwareLocation(unsigned soft_id, unsigned &module,
                                  unsigned &eta, unsigned &phi) {

  if (soft_id <= 0 || soft_id > towers_)
    std::exception("tower index out of bounds: ", soft_id,
                    " requested, but tower index range is [1, ", towers_, "]");

  int tower_idx = soft_id - 1;
  module = tower_idx / tow_per_module_;
  int local_tower_id = tower_idx - tow_per_module_ * module;
  phi = local_tower_id / tow_per_module_eta_;
  eta = local_tower_id % tow_per_module_eta_;
}

} // namespace jetreader
