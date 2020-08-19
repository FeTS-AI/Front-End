/*=========================================================================

  Program:   ALFABIS fast medical image registration programs
  Language:  C++
  Website:   github.com/pyushkevich/greedy
  Copyright (c) Paul Yushkevich, University of Pennsylvania. All rights reserved.

  This program is part of ALFABIS: Adaptive Large-Scale Framework for
  Automatic Biomedical Image Segmentation.

  ALFABIS development is funded by the NIH grant R01 EB017255.

  ALFABIS is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  ALFABIS is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with ALFABIS.  If not, see <http://www.gnu.org/licenses/>.

=========================================================================*/
#include "GreedyParameters.h"
#include "CommandLineHelper.h"


void
GreedyParameters
::SetToDefaults(GreedyParameters &param)
{
  param.dim = 2;
  param.mode = GreedyParameters::GREEDY;
  param.flag_dump_moving = false;
  param.flag_debug_deriv = false;
  param.flag_debug_aff_obj = false;
  param.dump_frequency = 1;
  param.epsilon_per_level = 1.0;
  param.sigma_pre.sigma = sqrt(3.0);
  param.sigma_pre.physical_units = false;
  param.sigma_post.sigma = sqrt(0.5);
  param.sigma_post.physical_units = false;
  param.threads = 0;
  param.metric = GreedyParameters::SSD;
  param.time_step_mode = GreedyParameters::SCALE;
  param.deriv_epsilon = 1e-4;
  param.flag_powell = false;
  param.warp_exponent = 6;
  param.warp_precision = 0.1;
  param.ncc_noise_factor = 0.001;
  param.affine_init_mode = VOX_IDENTITY;
  param.affine_dof = GreedyParameters::DOF_AFFINE;
  param.affine_jitter = 0.5;
  param.flag_float_math = false;
  param.flag_stationary_velocity_mode = false;
  param.flag_stationary_velocity_mode_use_lie_bracket = false;
  param.background = 0.0;
  param.current_weight = 1.0;

  param.iter_per_level.push_back(100);
  param.iter_per_level.push_back(100);

  // Moments of inertia parameters
  param.moments_flip_determinant = 0;
  param.flag_moments_id_covariance = false;
  param.moments_order = 1;
}

bool GreedyParameters::ParseCommandLine(const std::string &cmd, CommandLineHelper &cl)
{
  if(cmd == "-d")
    {
    this->dim = cl.read_integer();
    }
  else if(cmd == "-float")
    {
    this->flag_float_math = true;
    }
  else if(cmd == "-n")
    {
    this->iter_per_level = cl.read_int_vector();
    }
  else if(cmd == "-w")
    {
    this->current_weight = cl.read_double();
    }
  else if(cmd == "-e")
    {
    this->epsilon_per_level = cl.read_double_vector();
    }
  else if(cmd == "-m")
    {
    std::string metric_name = cl.read_string();
    if(metric_name == "NCC" || metric_name == "ncc")
      {
      this->metric = GreedyParameters::NCC;
      this->metric_radius = cl.read_int_vector();
      }
    else if(metric_name == "MI" || metric_name == "mi")
      {
      this->metric = GreedyParameters::MI;
      }
    else if(metric_name == "NMI" || metric_name == "nmi")
      {
      this->metric = GreedyParameters::NMI;
      }
    else if(metric_name == "MAHAL" || metric_name == "mahal")
      {
      this->metric = GreedyParameters::MAHALANOBIS;
      }
    }
  else if(cmd == "-tscale")
    {
    std::string mode = cl.read_string();
    if(mode == "SCALE" || mode == "scale")
      this->time_step_mode = GreedyParameters::SCALE;
    else if(mode == "SCALEDOWN" || mode == "scaledown")
      this->time_step_mode = GreedyParameters::SCALEDOWN;
    }
  else if(cmd == "-noise")
    {
    this->ncc_noise_factor = cl.read_double();
    }
  else if(cmd == "-s")
    {
    this->sigma_pre.sigma = cl.read_scalar_with_units(this->sigma_pre.physical_units);
    this->sigma_post.sigma = cl.read_scalar_with_units(this->sigma_post.physical_units);
    }
  else if(cmd == "-i")
    {
    ImagePairSpec ip;
    ip.weight = this->current_weight;
    ip.fixed = cl.read_existing_filename();
    ip.moving = cl.read_existing_filename();
    this->inputs.push_back(ip);
    }
  else if(cmd == "-id")
    {
    this->initial_warp = cl.read_existing_filename();
    }
  else if(cmd == "-ia")
    {
    this->affine_init_mode = RAS_FILENAME;
    this->affine_init_transform = cl.read_transform_spec();
    }
  else if(cmd == "-ia-identity" || cmd == "-iaid" || cmd == "-ia-id")
    {
    this->affine_init_mode = RAS_IDENTITY;
    }
  else if(cmd == "-ia-image-centers" || cmd == "-iaic" || cmd == "-ia-ic")
    {
    this->affine_init_mode = IMG_CENTERS;
    }
  else if(cmd == "-dof")
    {
    int dof = cl.read_integer();
    if(dof == 6)
      this->affine_dof = GreedyParameters::DOF_RIGID;
    else if(dof == 12)
        this->affine_dof = GreedyParameters::DOF_AFFINE;
    else throw GreedyException("DOF parameter only accepts 6 and 12 as values");
    }
  else if(cmd == "-jitter")
    {
    this->affine_jitter = cl.read_double();
    }
  else if(cmd == "-search")
    {
    this->rigid_search.iterations = cl.read_integer();

    std::string angle_cmd = cl.read_string();

    if(angle_cmd == "any" || angle_cmd == "ANY")
      {
      this->rigid_search.mode = ANY_ROTATION;
      this->rigid_search.sigma_angle = 0.0;
      }
    else if(angle_cmd == "flip" || angle_cmd == "FLIP")
      {
      this->rigid_search.mode = ANY_ROTATION_AND_FLIP;
      this->rigid_search.sigma_angle = 0.0;
      }
    else
      {
      this->rigid_search.mode = RANDOM_NORMAL_ROTATION;
      this->rigid_search.sigma_angle = atof(angle_cmd.c_str());
      }

    this->rigid_search.sigma_xyz = cl.read_double();
    }
  else if(cmd == "-it")
    {
    int nFiles = cl.command_arg_count();
    for(int i = 0; i < nFiles; i++)
      this->moving_pre_transforms.push_back(cl.read_transform_spec());
    }
  else if(cmd == "-gm")
    {
    this->gradient_mask = cl.read_existing_filename();
    }
  else if(cmd == "-gm-trim")
    {
    this->gradient_mask_trim_radius = cl.read_int_vector();
    }
  else if(cmd == "-mm")
    {
    this->moving_mask = cl.read_existing_filename();
    }
  else if(cmd == "-o")
    {
    this->output = cl.read_output_filename();
    }
  else if(cmd == "-dump-moving")
    {
    this->flag_dump_moving = true;
    }
  else if(cmd == "-powell")
    {
    this->flag_powell = true;
    }
  else if(cmd == "-dump-frequency" || cmd == "-dump-freq")
    {
    this->dump_frequency = cl.read_integer();
    }
  else if(cmd == "-debug-deriv")
    {
    this->flag_debug_deriv = true;
    }
  else if(cmd == "-debug-deriv-eps")
    {
    this->deriv_epsilon = cl.read_double();
    }
  else if(cmd == "-debug-aff-obj")
    {
    this->flag_debug_aff_obj = true;
    }
  else if(cmd == "-threads")
    {
    this->threads = cl.read_integer();
    }
  else if(cmd == "-a")
    {
    this->mode = GreedyParameters::AFFINE;
    }
  else if(cmd == "-moments")
    {
    this->mode = GreedyParameters::MOMENTS;

    // For backward compatibility allow no parameter, which defaults to order 1
    this->moments_order = cl.command_arg_count() > 0 ? cl.read_integer() : 2;
    if(this->moments_order != 1 && this->moments_order != 2)
      throw GreedyException("Parameter to -moments must be 1 or 2");
    }
  else if(cmd == "-brute")
    {
    this->mode = GreedyParameters::BRUTE;
    this->brute_search_radius = cl.read_int_vector();
    }
  else if(cmd == "-r")
    {
    this->mode = GreedyParameters::RESLICE;
    int nFiles = cl.command_arg_count();
    for(int i = 0; i < nFiles; i++)
      this->reslice_param.transforms.push_back(cl.read_transform_spec());
    }
  else if(cmd == "-iw")
    {
    this->mode = GreedyParameters::INVERT_WARP;
    this->invwarp_param.in_warp = cl.read_existing_filename();
    this->invwarp_param.out_warp = cl.read_output_filename();
    }
  else if(cmd == "-jac")
    {
    this->mode = GreedyParameters::JACOBIAN_WARP;
    this->jacobian_param.in_warp = cl.read_existing_filename();
    this->jacobian_param.out_det_jac = cl.read_output_filename();
    }
  else if(cmd == "-root")
    {
    this->mode = GreedyParameters::ROOT_WARP;
    this->warproot_param.in_warp = cl.read_existing_filename();
    this->warproot_param.out_warp = cl.read_output_filename();
    }

  else if(cmd == "-rm")
    {
    ResliceSpec rp;
    rp.interp = this->current_interp;
    rp.moving = cl.read_existing_filename();
    rp.output = cl.read_output_filename();
    this->reslice_param.images.push_back(rp);
    }
  else if(cmd == "-rs")
    {
    ResliceMeshSpec rp;
    rp.fixed = cl.read_existing_filename();
    rp.output = cl.read_output_filename();
    this->reslice_param.meshes.push_back(rp);
    }
  else if(cmd == "-rf")
    {
    this->reslice_param.ref_image = cl.read_existing_filename();
    }
  else if(cmd == "-rc")
    {
    this->reslice_param.out_composed_warp = cl.read_output_filename();
    }
  else if(cmd == "-rj")
    {
    this->reslice_param.out_jacobian_image = cl.read_output_filename();
    }
  else if(cmd == "-oinv")
    {
    this->inverse_warp = cl.read_output_filename();
    }
  else if(cmd == "-oroot")
    {
    this->root_warp = cl.read_output_filename();
    }
  else if(cmd == "-exp")
    {
    this->warp_exponent = cl.read_integer();
    }
  else if(cmd == "-sv")
    {
    this->flag_stationary_velocity_mode = true;
    this->flag_stationary_velocity_mode_use_lie_bracket = false;
    }
  else if(cmd == "-svlb")
    {
    this->flag_stationary_velocity_mode = true;
    this->flag_stationary_velocity_mode_use_lie_bracket = true;
    }
  else if(cmd == "-ri")
    {
    std::string mode = cl.read_string();
    if(mode == "nn" || mode == "NN" || mode == "0")
      {
      this->current_interp.mode = InterpSpec::NEAREST;
      }
    else if(mode == "linear" || mode == "LINEAR" || mode == "1")
      {
      this->current_interp.mode = InterpSpec::LINEAR;
      }
    else if(mode == "label" || mode == "LABEL")
      {
      this->current_interp.mode = InterpSpec::LABELWISE;
      this->current_interp.sigma.sigma = cl.read_scalar_with_units(
                                           this->current_interp.sigma.physical_units);
      }
    else
      {
      std::cerr << "Unknown interpolation mode" << std::endl;
      }
    }
  else if(cmd == "-rb")
    {
    this->current_interp.outside_value = cl.read_double();
    }
  else if(cmd == "-wp")
    {
    this->warp_precision = cl.read_double();
    }
  else if(cmd == "-det")
    {
    int det_value = cl.read_integer();
    if(det_value != -1 && det_value != 1)
      throw GreedyException("Incorrect -det parameter value %f", det_value);
    this->moments_flip_determinant = det_value;
    }
  else if(cmd == "-cov-id")
    {
    this->flag_moments_id_covariance = true;
    }
  else
    {
    return false;
    }

  return true;
}
