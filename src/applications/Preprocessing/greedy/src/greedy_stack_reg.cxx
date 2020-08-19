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
#include "CommandLineHelper.h"
#include "ShortestPath.h"

#include <iostream>
#include <sstream>
#include <cstdio>
#include <vector>
#include <string>
#include <algorithm>
#include <numeric>
#include <cerrno>

#include "itkMatrixOffsetTransformBase.h"
#include "itkImageAlgorithm.h"
#include "itkZeroFluxNeumannPadImageFilter.h"
#include "itkImageFileReader.h"
#include "itkImageSliceIteratorWithIndex.h"

#include "lddmm_common.h"
#include "lddmm_data.h"

struct StackParameters
{
  bool reuse;
  std::string output_dir;
  StackParameters()
    : reuse(false) {}
};


/** Parameters for splatting */
struct SplatParameters
{
  // Reference volume
  std::string reference;

  // When reference volume is not specified, a z-range must be given
  double z_first, z_last, z_step;

  // Which set of results to splat
  enum SplatSource {
    RAW, RECON, VOL_MATCH, VOL_ITER
  };

  SplatSource source_stage;

  // If relevant (VOL_ITER), the iteration
  unsigned int source_iter;

  // Splatting modes
  enum SplatMode {
    EXACT, NEAREST, LINEAR, PARZEN
  };

  SplatMode mode;

  // Tolerance on z-matching for EXACT
  double z_exact_tol;

  // Splatting sigma (in units of z)
  double sigma;

  // Output volume
  std::string fn_output;

  // Alternative manifest
  std::string fn_manifest;

  // Should alternative headers be ignored
  bool ignore_alt_headers;

  // Background values
  std::vector<double> background;

  // In-plane sigma
  double sigma_inplane;

  SplatParameters()
    : z_first(0.0), z_last(0.0), z_step(0.0),
      source_stage(RAW), source_iter(0),
      mode(EXACT), z_exact_tol(1e-6), sigma(0.0),
      ignore_alt_headers(false), background(1, 0.0),
      sigma_inplane(0.0) {}
};



/**
 * This class represents a reference to an image that may exist on disk, or may be
 * stored in memory. There is a limit on the amount of memory that can be used by
 * all the image refs, and images are rotated in and out of memory based on when
 * they were last accessed
 */
class ImageCache
{
public:

  ImageCache(unsigned long max_memory = 0l, unsigned int max_images = 0)
    : m_MaxMemory(max_memory), m_MaxImages(max_images), m_Counter(0l) {}

  template <typename TImage> typename TImage::Pointer GetImage(const std::string &filename)
  {
    // Check the cache for the image
    auto it = m_Cache.find(filename);
    if(it != m_Cache.end())
      {
      TImage *image = dynamic_cast<TImage *>(std::get<2>(it->second).GetPointer());
      if(!image)
        throw GreedyException("Type mismatch in image cache");
      typename TImage::Pointer image_ptr = image;
      return image_ptr;
      }

    // Image does not exist in cache, load it
    typedef itk::ImageFileReader<TImage> ReaderType;
    typename ReaderType::Pointer reader = ReaderType::New();
    reader->SetFileName(filename.c_str());
    reader->Update();
    typename TImage::Pointer image_ptr = reader->GetOutput();

    // Get the size of the image in bytes
    unsigned long img_size = image_ptr->GetPixelContainer()->Size()
                             * sizeof (typename TImage::PixelContainer::Element);

    // If the size of the image is too large, we need to reduce the size of the cache
    this->ShrinkCache(img_size, 1);

    // Add the new image
    m_Cache[filename] = std::make_tuple(m_Counter++, img_size, image_ptr);
    m_UsedMemory += img_size;

    // Return the image
    return image_ptr;
  }

  void ShrinkCache(unsigned long new_bytes, unsigned int new_images)
  {
    // Remove something from the cache until it's not empty and the constraints of the
    // cache are satisfied
    while(IsCacheFull(new_bytes, new_images) && m_Cache.size() > 0)
      {
      // Find the oldest entry in the cache
      std::map<unsigned long, std::string> sort_map;
      for(auto it : m_Cache)
        sort_map[std::get<0>(it.second)] = it.first;

      // Remove the first (oldest) entry
      auto it_erase = m_Cache.find(sort_map.begin()->second);
      m_UsedMemory -= std::get<1>(it_erase->second);
      m_Cache.erase(it_erase);
      }
  }

  bool IsCacheFull(unsigned long new_bytes, unsigned int new_images)
  {
    if(m_MaxMemory > 0 && m_UsedMemory + new_bytes > m_MaxMemory)
      return true;

    if(m_MaxImages > 0 && m_Cache.size() + new_images > m_MaxImages)
      return true;

    return false;
  }

  void PurgeCache()
  {
    m_Cache.clear();
    m_UsedMemory = 0;
  }

protected:

  // Cache entry (age, size, pointer)
  typedef std::tuple<unsigned long, unsigned long, itk::Object::Pointer> CacheEntry;
  typedef std::map<std::string, CacheEntry> CacheType;

  // Cache for images
  CacheType m_Cache;

  unsigned long m_MaxMemory, m_UsedMemory;
  unsigned int m_MaxImages;
  unsigned long m_Counter;
};


// How to specify how many neighbors a slice will be registered to?
// - minimum of one neighbor
// - maximum of <user_specified> neighbors
// - maximum offset

int usage(const std::string &stage = std::string())
{
  std::map<std::string, std::string> utext;
  utext[std::string()] = {
    #include "stackg_usage_main.h"
    0x00
  };

  utext["init"] = {
    #include "stackg_usage_init.h"
    0x00
  };

  utext["recon"] = {
    #include "stackg_usage_recon.h"
    0x00
  };

  utext["volmatch"] =  {
    #include "stackg_usage_volmatch.h"
    0x00
  };

  utext["voliter"] =  {
    #include "stackg_usage_voliter.h"
    0x00
  };

  utext["splat"] =  {
    #include "stackg_usage_splat.h"
    0x00
  };

  if(utext.find(stage) == utext.end())
    std::cout << utext[std::string()] << std::endl;
  else
    std::cout << utext[stage] << std::endl;
  return -1;
}


/**
 * A representation of the project
 */
class StackGreedyProject
{
public:

  // API typedefs
  typedef LDDMMData<double, 2> LDDMMType;
  typedef LDDMMData<double, 3> LDDMMType3D;
  typedef GreedyApproach<2, double> GreedyAPI;

  // Image typedefs
  typedef LDDMMType::CompositeImageType SlideImageType;
  typedef LDDMMType::CompositeImagePointer SlideImagePointer;
  typedef LDDMMType3D::CompositeImageType VolumeImage;
  typedef LDDMMType3D::CompositeImagePointer VolumePointer;


  /** Set of enums used to refer to files in the project directory */
  enum FileIntent {
    MANIFEST_FILE = 0, CONFIG_ENTRY, AFFINE_MATRIX, METRIC_VALUE, ACCUM_MATRIX, ACCUM_RESLICE,
    VOL_INIT_MATRIX, VOL_SLIDE, VOL_MASK_SLIDE, VOL_BEST_INIT_MATRIX,
    VOL_ITER_MATRIX, VOL_ITER_WARP, ITER_METRIC_DUMP
  };

  /** Constructor */
  StackGreedyProject(std::string project_dir, const StackParameters &param)
  {
    m_ProjectDir = project_dir;
    m_GlobalParam = param;
  }

  /** Initialize the project */
  void InitializeProject(std::string fn_manifest, std::string default_ext = "nii.gz")
  {
    // Read the manifest and write a copy to the project dir
    this->ReadManifest(fn_manifest);
    this->WriteManifest(GetFilenameForGlobal(MANIFEST_FILE));

    // Read the default extension and save it
    this->m_DefaultImageExt = default_ext;
    this->SaveConfigKey("DefaultImageExt", m_DefaultImageExt);

    // Report what has been done
    printf("stack_greedy: Project initialized in %s\n", m_ProjectDir.c_str());
  }

  /** Restore the initialized project */
  void RestoreProject()
  {
    this->ReadManifest(GetFilenameForGlobal(MANIFEST_FILE));
    m_DefaultImageExt = this->LoadConfigKey("DefaultImageExt", std::string(".nii.gz"));
  }


  template <typename T> void SaveConfigKey(const std::string &key, const T &value)
  {
    std::string fn = GetFilenameForGlobal(CONFIG_ENTRY, key.c_str());
    std::ofstream fout(fn);
    fout << value;
  }

  template <typename T> T LoadConfigKey(const std::string &key, const T &def_value)
  {
    std::string fn = GetFilenameForGlobal(CONFIG_ENTRY, key.c_str());
    std::ifstream fin(fn);
    T value;
    if(fin.good())
      fin >> value;
    else
      value = def_value;
    return value;
  }

  void ReadManifest(const std::string &fn_manifest)
  {
    // Reset the slices
    m_Slices.clear();
    m_SortedSlices.clear();

    // Read the manifest file
    std::ifstream fin(fn_manifest);
    std::string f_line;
    while(std::getline(fin, f_line))
      {
      // Read the values from the manifest
      std::istringstream iss(f_line);
      SliceData slice;
      if(!(iss >> slice.unique_id >> slice.z_pos >> slice.is_leader >> slice.raw_filename))
        throw GreedyException("Error reading manifest file, line '%s'", f_line.c_str());

      // Check that the manifest points to a real file
      if(!itksys::SystemTools::FileExists(slice.raw_filename.c_str(), true))
        throw GreedyException("File %s referenced in the manifest does not exist", slice.raw_filename.c_str());

      // Get an absolute filename
      slice.raw_filename = itksys::SystemTools::CollapseFullPath(slice.raw_filename.c_str());

      // Add to sorted list
      m_SortedSlices.insert(std::make_pair(slice.z_pos, m_Slices.size()));
      m_Slices.push_back(slice);
      }
  }

  bool CanSkipFile(const std::string &fn)
  {
    if(fn.length() == 0)
      return true;

    return m_GlobalParam.reuse && itksys::SystemTools::FileExists(fn.c_str(), true);
  }

  void WriteManifest(const std::string &fn_manifest)
  {
    std::ofstream fout(fn_manifest);
    for(auto slice : m_Slices)
      fout << slice.unique_id << " " << slice.z_pos << " " << slice.is_leader << slice.raw_filename << std::endl;
  }

  void ReconstructStack(double z_range, double z_epsilon, const GreedyParameters &gparam)
  {
    // Configure the threads
    GreedyAPI::ConfigThreads(gparam);

    // Store the z-parameters (although we probably do not need them)
    this->SaveConfigKey("Z_Range", z_range);
    this->SaveConfigKey("Z_Epsilon", z_epsilon);

    // This array represents the slice graph structure. Each slice uses one or more
    // neighbors as references for registration. However, some of the slices are
    // leaders and some are followers. The follower slices use leader slices as 
    // references, but leader slices ignore the follower slices. In the graph,
    // we represent this by having edges from reference images to moving images,
    // thus edges L->F, L->L but not F->F or F->L

    // We keep for each slice the list of z-sorted neigbors
    std::vector<slice_ref_set> slice_nbr(m_Slices.size());
    unsigned int n_edges = 0;

    // Forward pass
    for(auto it = m_SortedSlices.begin(); it != m_SortedSlices.end(); ++it)
      {
      // If this slide is a follower, don't add any neighbors
      if(m_Slices[it->second].is_leader)
        {
        // Add at least the following slice
        auto it_next = it; ++it_next;
        unsigned int n_added = 0;

        // Now add all the slices in the range
        while(it_next != m_SortedSlices.end()
          && (n_added < 1 || fabs(it->first - it_next->first) < z_range))
          {
          slice_nbr[it->second].insert(*it_next);
          n_added++; n_edges++; ++it_next;
          }
        }
      }

    // Backward pass
    for(auto it = m_SortedSlices.rbegin(); it != m_SortedSlices.rend(); ++it)
      {
      // If this slide is a follower, don't add any neighbors
      if(m_Slices[it->second].is_leader)
        {
        // Add at least the following slice
        auto it_next = it; ++it_next;
        unsigned int n_added = 0;

        // Now add all the slices in the range
        while(it_next != m_SortedSlices.rend()
              && (n_added < 1 || fabs(it->first - it_next->first) < z_range))
          {
          slice_nbr[it->second].insert(*it_next);
          n_added++; n_edges++; ++it_next;
          }
        }
      }

    // Set up a cache for loaded images. These images can be cycled in and out of memory
    // depending on need. TODO: let user configure cache sizes
    ImageCache slice_cache(0, 20);

    // At this point we can create a rigid adjacency structure for the graph-theoretic algorithm,
    vnl_vector<unsigned int> G_adjidx(m_SortedSlices.size()+1, 0u);
    vnl_vector<unsigned int> G_adj(n_edges, 0u);
    vnl_vector<double> G_edge_weight(n_edges, DijkstraShortestPath<double>::INFINITE_WEIGHT);

    for(unsigned int k = 0, p = 0; k < m_Slices.size(); k++)
      {
      G_adjidx[k+1] = G_adjidx[k] + (unsigned int) slice_nbr[k].size();
      for(auto it : slice_nbr[k])
        G_adj[p++] = it.second;
      }

    // Set up the graph of all registrations. Each slice is a node and edges are between each
    // slice and its closest slices, as well as between each slice and slices in the z-range
    typedef std::tuple<unsigned int, unsigned int, double> GraphEdge;
    std::set<GraphEdge> slice_graph;

    // Perform rigid registration between pairs of images. We should do this in a way that
    // the number of images loaded and unloaded is kept to a minimum, without filling memory.
    // The best way to do so would be to progress in z order and release images that are too
    // far behind in z to be included for the current 'reference' image
    for(auto it : m_SortedSlices)
      {
      // Skip the slide if it is a follower (followers are not used as reference slices)
      if(!m_Slices[it.second].is_leader)
        continue;

      const auto &nbr = slice_nbr[it.second];

      // Read the reference slide from the cache
      SlideImagePointer i_ref =
          slice_cache.GetImage<SlideImageType>(m_Slices[it.second].raw_filename);

      // Iterate over the neighbor slices
      unsigned int n_pos = 0;
      for(auto it_n : nbr)
        {
        // Load or retrieve the corresponding image
        SlideImagePointer i_mov =
            slice_cache.GetImage<SlideImageType>(m_Slices[it_n.second].raw_filename);

        // Get the filenames that will be generated by registration
        std::string fn_matrix = GetFilenameForSlicePair(m_Slices[it.second], m_Slices[it_n.second], AFFINE_MATRIX);
        std::string fn_metric = GetFilenameForSlicePair(m_Slices[it.second], m_Slices[it_n.second], METRIC_VALUE);
        double pair_metric = 1e100;

        // Perform registration or reuse existing registration results
        if(CanSkipFile(fn_matrix) && CanSkipFile(fn_metric))
          {
          std::ifstream fin(fn_metric);
          fin >> pair_metric;
          }
        else
          {
          // Perform the registration between i_ref and i_mov
          GreedyAPI greedy_api;

          // Make a copy of the template parameters
          GreedyParameters my_param = gparam;

          // Set up the image pair for registration
          ImagePairSpec img_pair(m_Slices[it.second].raw_filename, m_Slices[it_n.second].raw_filename);
          greedy_api.AddCachedInputObject(m_Slices[it.second].raw_filename, i_ref.GetPointer());
          greedy_api.AddCachedInputObject(m_Slices[it_n.second].raw_filename, i_mov.GetPointer());
          my_param.inputs.push_back(img_pair);

          // Set other parameters
          my_param.affine_dof = GreedyParameters::DOF_RIGID;
          my_param.affine_init_mode = IMG_CENTERS;

          // Set up the output of the affine
          my_param.output = fn_matrix;

          // Perform affine/rigid
          printf("#############################\n");
          printf("### Fixed :%s   Moving %s ###\n", m_Slices[it.second].unique_id.c_str(),m_Slices[it_n.second].unique_id.c_str());
          printf("#############################\n");
          greedy_api.RunAffine(my_param);

          // Get the metric for the affine registration
          pair_metric = greedy_api.GetLastMetricReport().TotalMetric;
          std::cout << "Last metric value: " << pair_metric << std::endl;

          // Normalize the metric to give the actual mean NCC
          pair_metric /= -10000.0 * i_ref->GetNumberOfComponentsPerPixel();
          std::ofstream f_metric(fn_metric);
          f_metric << pair_metric << std::endl;
          }

        // Map the metric value into a weight
        double weight = (1.0 - pair_metric) * pow(1 + z_epsilon, fabs(it_n.first - it.first));

        // Regardless of whether we did registration or not, record the edge in the graph
        G_edge_weight[G_adjidx[it.second] + n_pos++] = weight;
        }
      }

    // Run the shortest path computations
    DijkstraShortestPath<double> dijkstra((unsigned int) m_Slices.size(),
                                          G_adjidx.data_block(), G_adj.data_block(), G_edge_weight.data_block());

    // Compute the shortest paths from every slice to the rest and record the total distance. This will
    // help generate the root of the tree
    int i_root = -1;
    double best_root_dist = 0.0;
    for(unsigned int i = 0; i < m_Slices.size(); i++)
      {
      if(m_Slices[i].is_leader) 
        {
        dijkstra.ComputePathsFromSource(i);
        double root_dist = 0.0;
        for(unsigned int j = 0; j < m_Slices.size(); j++)
          root_dist += dijkstra.GetDistanceArray()[j];
        std::cout << "Root distance " << i << " : " << root_dist << std::endl;
        if(i_root < 0 || best_root_dist > root_dist)
          {
          i_root = i;
          best_root_dist = root_dist;
          }
        }
      }

    // No root? We have a problem
    if(i_root < 0)
      throw GreedyException("No root found for the registration graph. Did you specify any leader slices in the manifest?");

    // Store the root for reference
    SaveConfigKey("RootSlide", m_Slices[i_root].unique_id);

    // Compute the composed transformations between the root and each of the inputs
    dijkstra.ComputePathsFromSource(i_root);

    // Load the root image into memory
    LDDMMType::ImagePointer img_root;
    LDDMMType::img_read(m_Slices[i_root].raw_filename.c_str(), img_root);

    // Apply some padding to the root image.
    typedef itk::ZeroFluxNeumannPadImageFilter<LDDMMType::ImageType, LDDMMType::ImageType> PadFilter;
    PadFilter::Pointer fltPad = PadFilter::New();
    fltPad->SetInput(img_root);

    // Determine the amount of padding to add
    unsigned int max_dim =
        std::max(img_root->GetBufferedRegion().GetSize()[0],img_root->GetBufferedRegion().GetSize()[1]);
    itk::Size<2> pad_size;
    pad_size.Fill(max_dim / 4);
    fltPad->SetPadBound(pad_size);
    fltPad->Update();

    // Store the result
    LDDMMType::ImagePointer img_root_padded = fltPad->GetOutput();

    // The padded image has a non-zero index, which causes problems downstream for GreedyAPI.
    // To account for this, we save and load the image
    // TODO: handle this internally using a filter!
    LDDMMType::img_write(img_root_padded, "/tmp/padded.nii.gz");
    img_root_padded = LDDMMType::img_read("/tmp/padded.nii.gz");

    // Compute transformation for each slice
    for(unsigned int i = 0; i < m_Slices.size(); i++)
      {
      // Initialize the total transform matrix
      vnl_matrix<double> t_accum(3, 3, 0.0);
      t_accum.set_identity();

      // Traverse the path
      unsigned int i_curr = i, i_prev = dijkstra.GetPredecessorArray()[i];
      std::cout << "Chain for " << i << " : ";
      while(i_prev != DijkstraShortestPath<double>::NO_PATH && (i_prev != i_curr))
        {
        // Load the matrix
        std::string fn_matrix =
            GetFilenameForSlicePair(m_Slices[i_prev], m_Slices[i_curr], AFFINE_MATRIX);
        vnl_matrix<double> t_step = GreedyAPI::ReadAffineMatrix(TransformSpec(fn_matrix));

        // Accumulate the total transformation
        t_accum = t_accum * t_step;

        std::cout << i_prev << " ";

        // Go to the next edge
        i_curr = i_prev;
        i_prev = dijkstra.GetPredecessorArray()[i_curr];
        }

      std::cout << std::endl;

      // Store the accumulated transform
      std::string fn_accum_matrix = GetFilenameForSlice(m_Slices[i], ACCUM_MATRIX);
      GreedyAPI::WriteAffineMatrix(fn_accum_matrix, t_accum);

      // Write a resliced image
      std::string fn_accum_reslice = GetFilenameForSlice(m_Slices[i], ACCUM_RESLICE);

      // Hold the resliced image in memory
      LDDMMType::CompositeImagePointer img_reslice = LDDMMType::CompositeImageType::New();

      // Only do reslice if necessary
      if(!CanSkipFile(fn_accum_reslice))
        {
        // Perform the registration between i_ref and i_mov
        GreedyAPI greedy_api;

        // Make a copy of the template parameters
        GreedyParameters my_param = gparam;

        // Set up the image pair for registration
        my_param.reslice_param.ref_image = "root_slice_padded";
        my_param.reslice_param.images.push_back(ResliceSpec(m_Slices[i].raw_filename, fn_accum_reslice));
        my_param.reslice_param.transforms.push_back(TransformSpec(fn_accum_matrix));
        greedy_api.AddCachedInputObject("root_slice_padded", img_root_padded.GetPointer());
        greedy_api.AddCachedInputObject(m_Slices[i].raw_filename,
                                        slice_cache.GetImage<SlideImageType>(m_Slices[i].raw_filename));
        greedy_api.AddCachedOutputObject(fn_accum_reslice, img_reslice.GetPointer(), true);
        greedy_api.RunReslice(my_param);
        }
      else
        {
        // Just read the resliced image
        img_reslice = LDDMMType::cimg_read(fn_accum_reslice.c_str());
        }
      }
  }

  static SlideImagePointer ExtractSliceFromVolume(VolumePointer vol, double z_pos)
  {
    VolumePointer vol_slice = LDDMMType3D::CompositeImageType::New();
    typename LDDMMType3D::RegionType reg_slice = vol->GetBufferedRegion();
    reg_slice.GetModifiableSize()[2] = 1;
    vol_slice->CopyInformation(vol);
    vol_slice->SetRegions(reg_slice);
    vol_slice->Allocate();

    // Adjust the origin of the slice
    auto origin_slice = vol_slice->GetOrigin();
    origin_slice[2] = z_pos;
    vol_slice->SetOrigin(origin_slice);

    // Generate a blank deformation field
    LDDMMType3D::VectorImagePointer zero_warp = LDDMMType3D::new_vimg(vol_slice);

    // Sample the slice from the volume
    LDDMMType3D::interp_cimg(vol, zero_warp, vol_slice, false, true, 0.0);

    // Now drop the dimension of the slice to 2D
    LDDMMType::RegionType reg_slice_2d;
    LDDMMType::CompositeImageType::PointType origin_2d;
    LDDMMType::CompositeImageType::SpacingType spacing_2d;
    LDDMMType::CompositeImageType::DirectionType dir_2d;

    for(unsigned int a = 0; a < 2; a++)
      {
      reg_slice_2d.SetIndex(a, reg_slice.GetIndex(a));
      reg_slice_2d.SetSize(a, reg_slice.GetSize(a));
      origin_2d[a] = vol_slice->GetOrigin()[a];
      spacing_2d[a] = vol_slice->GetSpacing()[a];
      dir_2d(a,0) = vol_slice->GetDirection()(a,0);
      dir_2d(a,1) = vol_slice->GetDirection()(a,1);
      }

    SlideImagePointer vol_slice_2d = SlideImageType::New();
    vol_slice_2d->SetRegions(reg_slice_2d);
    vol_slice_2d->SetOrigin(origin_2d);
    vol_slice_2d->SetDirection(dir_2d);
    vol_slice_2d->SetSpacing(spacing_2d);
    vol_slice_2d->SetNumberOfComponentsPerPixel(vol_slice->GetNumberOfComponentsPerPixel());
    vol_slice_2d->Allocate();

    // Copy data between the pixel containers
    itk::ImageAlgorithm::Copy(vol_slice.GetPointer(), vol_slice_2d.GetPointer(),
                              vol_slice->GetBufferedRegion(), vol_slice_2d->GetBufferedRegion());

    return vol_slice_2d;
  }

  void InitialMatchToVolume(const std::string &fn_volume, const std::string &fn_mask,
                            const GreedyParameters &gparam)
  {
    // Configure the threads
    GreedyAPI::ConfigThreads(gparam);

    // Read the 3D volume into memory
    LDDMMType3D::CompositeImagePointer vol = LDDMMType3D::cimg_read(fn_volume.c_str());

    // Read the mask into memory
    LDDMMType3D::CompositeImagePointer mask;
    if(fn_mask.size())
      mask = LDDMMType3D::cimg_read(fn_mask.c_str());

    // Extract target slices from the 3D volume
    for(unsigned int i = 0; i < m_Slices.size(); i++)
      {
      // Filename for the volume slice corresponding to current slide
      std::string fn_vol_slide = GetFilenameForSlice(m_Slices[i], VOL_SLIDE);

      // Output matrix for this registration
      std::string fn_vol_init_matrix = GetFilenameForSlice(m_Slices[i], VOL_INIT_MATRIX);

      // Is there a mask?
      std::string fn_vol_mask_slide = (fn_mask.length())
                                      ? GetFilenameForSlice(m_Slices[i], VOL_MASK_SLIDE)
                                      : std::string();

      if(!CanSkipFile(fn_vol_slide) || !CanSkipFile(fn_vol_init_matrix) || !CanSkipFile(fn_vol_mask_slide))
        {
        // Extract the slice from the 3D image
        SlideImagePointer vol_slice_2d = ExtractSliceFromVolume(vol, m_Slices[i].z_pos);

        // Write the 2d slice
        LDDMMType::cimg_write(vol_slice_2d, fn_vol_slide.c_str());

        // Only consider leader slides for matching to volume
        if(m_Slices[i].is_leader)
          {
          // Try registration between resliced slide and corresponding volume slice with
          // a brute force search. This will be used to create a median transformation
          // between slide space and volume space. Since the volume may come with a mask,
          // we use volume slice as fixed, and the slide image as moving
          GreedyAPI greedy_api;
          GreedyParameters my_param = gparam;

          // Set up the image pair for registration
          std::string fn_accum_reslice = GetFilenameForSlice(m_Slices[i], ACCUM_RESLICE);

          ImagePairSpec img_pair("vol_slice", fn_accum_reslice, 1.0);
          greedy_api.AddCachedInputObject("vol_slice", vol_slice_2d.GetPointer());
          my_param.inputs.push_back(img_pair);

          // Handle the mask
          if(fn_vol_mask_slide.length())
            {
            // TODO: we are not caching because of different image types
            SlideImagePointer mask_slice_2d = ExtractSliceFromVolume(mask, m_Slices[i].z_pos);
            LDDMMType::cimg_write(mask_slice_2d, fn_vol_mask_slide.c_str());
            my_param.gradient_mask = fn_vol_mask_slide;
            }

          // Set other parameters
          my_param.affine_dof = GreedyParameters::DOF_AFFINE;
          my_param.affine_init_mode = IMG_CENTERS;

          // Set up the output of the affine
          my_param.output = fn_vol_init_matrix;

          // Run the affine registration
          greedy_api.RunAffine(my_param);
          }
        }
      }

    // Now we have a large set of per-slice matrices. We next try each matrix on each pair of
    // slices and store the metric, with the goal of finding a matrix that will provide the 
    // best possible match.
    std::vector<double> accum_metric(m_Slices.size(), 0.0);
    for(unsigned int i = 0; i < m_Slices.size(); i++)
      {
      if(!m_Slices[i].is_leader)
        continue;

      // Load the images to avoid N^2 IO operations
      LDDMMType::CompositeImagePointer vol_slice = 
        LDDMMType::cimg_read(GetFilenameForSlice(m_Slices[i], VOL_SLIDE).c_str());

      LDDMMType::CompositeImagePointer acc_slice = 
        LDDMMType::cimg_read(GetFilenameForSlice(m_Slices[i], ACCUM_RESLICE).c_str());

      // Loop over matrices
      for(unsigned int k = 0; k < m_Slices.size(); k++)
        {
        if(!m_Slices[k].is_leader)
          continue;

        // This API is for metric computation
        GreedyAPI greedy_api;
        GreedyParameters my_param = gparam;

        // Same image pairs as before
        ImagePairSpec img_pair("vol_slice", "acc_slice", 1.0);
        greedy_api.AddCachedInputObject("vol_slice", vol_slice.GetPointer());
        greedy_api.AddCachedInputObject("acc_slice", acc_slice.GetPointer());
        my_param.inputs.push_back(img_pair);

        // TODO: this is really bad, can't cache mask images
        if(fn_mask.length())
          my_param.gradient_mask = GetFilenameForSlice(m_Slices[i], VOL_MASK_SLIDE);

        // Set other parameters
        my_param.affine_init_mode = RAS_FILENAME;
        my_param.affine_init_transform = GetFilenameForSlice(m_Slices[k], VOL_INIT_MATRIX);

        // Run affine to get the metric value
        MultiComponentMetricReport metric_report;
        greedy_api.ComputeMetric(my_param, metric_report);
        printf("Slide %03d matrix %03d metric %8.4f\n", i, k, metric_report.TotalMetric);
        accum_metric[k] += metric_report.TotalMetric;
        }
      }

    // Now find the matrix with the best overall metric
    int k_best = -1; double m_best = 0.0;
    for(unsigned int k = 0; k < m_Slices.size(); k++)
      {
      if(!m_Slices[k].is_leader) 
        continue;

      printf("Across-slice metric for matrix %04d: %8.4f\n", k, accum_metric[k]);
      if(k_best < 0 || accum_metric[k] > m_best) 
        {
        k_best = k;
        m_best = accum_metric[k];
        }
      }

    // The median affine
    vnl_matrix<double> M_best =
      GreedyAPI::ReadAffineMatrix(GetFilenameForSlice(m_Slices[k_best], VOL_INIT_MATRIX));

    // Write the median affine to a file
    GreedyAPI::WriteAffineMatrix(GetFilenameForGlobal(VOL_BEST_INIT_MATRIX), M_best);

    // Now write the complete initial to-volume transform for each slide
    for(unsigned int i = 0; i < m_Slices.size(); i++)
      {
      vnl_matrix<double> M_root =
          GreedyAPI::ReadAffineMatrix(GetFilenameForSlice(m_Slices[i], ACCUM_MATRIX));

      vnl_matrix<double> M_vol = M_root * M_best;

      GreedyAPI::WriteAffineMatrix(
            GetFilenameForSlice(m_Slices[i], VOL_ITER_MATRIX, 0), M_vol);
      }
  }

  // Now that we have the affine initialization from the histology space to the volume space, we can
  // perform iterative optimization, where each slice is matched to its neighbors and to the
  // corresponding MRI slice. The only issue here is how do we want to use the graph in this
  // process: we don't want the bad neighbors to pull the registration away from the good
  // solution. On the other hand, we can expect the bad slices to eventually auto-correct. It seems
  // that the proper approach would be to down-weigh certain slices by their metric, but then
  // again, do we want to do this based on initial metric or current metric. For now, we can start
  // by just using same weights.
  void IterativeMatchToVolume(unsigned int n_affine, unsigned int n_deform, unsigned int i_first,
                              unsigned int i_last, double w_volume, const GreedyParameters &gparam)
  {
    // Set up a cache for loaded images. These images can be cycled in and out of memory
    // depending on need. TODO: let user configure cache sizes
    ImageCache slice_cache(0, 20);

    // What iteration?
    if(i_first > i_last || i_first == 0 || i_last > n_affine + n_deform)
      throw GreedyException("Iteration range (%d, %d) is out of range [1, %d]",
                            i_first, i_last, n_affine + n_deform);

    // Set the iteration specs
    SaveConfigKey("AffineIterations", n_affine);
    SaveConfigKey("DeformableIterations", n_deform);

    // Iterate
    for(unsigned int iter = i_first; iter <= i_last; ++iter)
      {
      // Randomly shuffle the order in which slices are considered
      std::vector<unsigned int> ordering(m_Slices.size());
      std::iota(ordering.begin(), ordering.end(), 0);
      std::random_shuffle(ordering.begin(), ordering.end());

      // Keep track of the total neighbor metric and total volume metric
      double total_to_nbr_metric = 0.0;
      double total_to_vol_metric = 0.0;

      // Iterate over the ordering
      for(unsigned int k : ordering)
        {
        // The output filename for this affine registration
        std::string fn_result =
            iter <= n_affine
            ? GetFilenameForSlice(m_Slices[k], VOL_ITER_MATRIX, iter)
            : GetFilenameForSlice(m_Slices[k], VOL_ITER_WARP, iter);

        // Has this already been done? Then on to the next!
        if(CanSkipFile(fn_result))
          continue;

        // Get the pointer to the current slide (used as moving image)
        SlideImagePointer img_slide = slice_cache.GetImage<SlideImageType>(m_Slices[k].raw_filename);

        // Get the corresponding slice from the 3D volume (it's already saved in the project)
        SlideImagePointer vol_slice_2d = slice_cache.GetImage<SlideImageType>(
                                           GetFilenameForSlice(m_Slices[k], VOL_SLIDE));

        // Load the mask as well
        GreedyAPI::ImagePointer mask_slice_2d;
        std::string fn_mask_slice = GetFilenameForSlice(m_Slices[k], VOL_MASK_SLIDE);
        if(itksys::SystemTools::FileExists(fn_mask_slice.c_str(), true))
          mask_slice_2d = slice_cache.GetImage<GreedyAPI::ImageType>(fn_mask_slice);

        // Set up the registration. We are registering to the volume and to the transformed
        // adjacent slices. We should do everything in the space of the MRI volume because
        // (a) it should be large enough to cover the histology and (b) there might be a mask
        // in this space, while we cannot expect there to be a mask in the other space.

        // Find the adjacent slices. TODO: there is all kinds of stuff that could be done here,
        // like allowing a z-range for adjacent slices registration, modulating weight by the
        // distance, and detecting and down-weighting 'bad' slices. For now just pick the slices
        // immediately below and above the current slice
        //
        // Correction (Sep 2019): we now search for the adjacent leader slices. Non-leader slices
        // are not used as they are assumed to be less reliable
        slice_ref_set k_nbr;

        // Find slice before k that is a leader slice
        for(auto itr = m_SortedSlices.rbegin(); itr != m_SortedSlices.rend(); itr++)
          {
          if(itr->first < m_Slices[k].z_pos && m_Slices[itr->second].is_leader)
            {
            k_nbr.insert(*itr);
            break;
            }
          }

        // Find slice after k that is a leader slice
        for(auto itf = m_SortedSlices.begin(); itf != m_SortedSlices.end(); itf++)
          {
          if(itf->first > m_Slices[k].z_pos && m_Slices[itf->second].is_leader)
            {
            k_nbr.insert(*itf);
            break;
            }
          }

        // Create the greedy API for the main registration task
        GreedyAPI api_reg;
        api_reg.AddCachedInputObject("moving", img_slide);
        api_reg.AddCachedInputObject("volume_slice", vol_slice_2d);

        // We need to hold on to the resliced image pointers, because otherwise they will be deallocated
        std::vector<SlideImagePointer> resliced_neighbors(m_Slices.size());

        // Set up the main registration pair
        GreedyParameters param_reg = gparam;
        param_reg.inputs.push_back(ImagePairSpec("volume_slice", "moving", w_volume));    

        // Handle mask
        if(mask_slice_2d)
          {
          api_reg.AddCachedInputObject("mask_slice", mask_slice_2d);
          param_reg.gradient_mask = "mask_slice";
          }

        // Handle each of the neighbors
        for(auto nbr : k_nbr)
          {
          unsigned int j = nbr.second;

          // Create an image pointer for the reslicing output
          resliced_neighbors[j] = SlideImageType::New();

          // Each of the neighbor slices needs to be resliced using last iteration's transform. We
          // could cache these images, but then again, it does not take so much to do this on the
          // fly. For now we will do this on the fly.
          GreedyAPI api_reslice;
          api_reslice.AddCachedInputObject("vol_slice", vol_slice_2d);
          api_reslice.AddCachedOutputObject("output", resliced_neighbors[j], false);

          GreedyParameters param_reslice = gparam;
          param_reslice.reslice_param.ref_image = "vol_slice";
          param_reslice.reslice_param.images.push_back(ResliceSpec(m_Slices[j].raw_filename, "output"));

          // Was the previous iteration a deformable iteration? If so, apply the warp
          if(iter - 1 <= n_affine)
            {
            param_reslice.reslice_param.transforms.push_back(
                  TransformSpec(GetFilenameForSlice(m_Slices[j], VOL_ITER_MATRIX, iter-1)));
            }
          else
            {
            param_reslice.reslice_param.transforms.push_back(
                  TransformSpec(GetFilenameForSlice(m_Slices[j], VOL_ITER_WARP, iter-1)));
            param_reslice.reslice_param.transforms.push_back(
                  TransformSpec(GetFilenameForSlice(m_Slices[j], VOL_ITER_MATRIX, iter-1)));
            }

          // Perform the reslicing
          api_reslice.RunReslice(param_reslice);

          // Add the image pair to the registration
          char fixed_fn[64];
          sprintf(fixed_fn, "neighbor_%03d", j);
          api_reg.AddCachedInputObject(fixed_fn, resliced_neighbors[j]);

          param_reg.inputs.push_back(ImagePairSpec(fixed_fn, "moving", 1.0));
          }

        printf("#############################\n");
        printf("### Iter :%d   Slide %s ###\n", iter, m_Slices[k].unique_id.c_str());
        printf("#############################\n");

        // What kind of registration are we doing at this iteration?
        if(iter <= n_affine)
          {
          // Specify the DOF, etc
          param_reg.affine_dof = GreedyParameters::DOF_AFFINE;
          param_reg.affine_init_mode = RAS_FILENAME;
          param_reg.affine_init_transform =
              TransformSpec(GetFilenameForSlice(m_Slices[k], VOL_ITER_MATRIX, iter-1));
          param_reg.rigid_search = RigidSearchSpec();

          // Specify the output
          param_reg.output = fn_result;

          // Run this registration!
          api_reg.RunAffine(param_reg);
          }
        else
          {
          // Apply the last affine transformation
          param_reg.affine_init_mode = VOX_IDENTITY;
          param_reg.moving_pre_transforms.push_back(
                TransformSpec(GetFilenameForSlice(m_Slices[k], VOL_ITER_MATRIX, iter-1)));

          // Specify the output
          param_reg.output = fn_result;

          // Make sure the warps are not truncated, since histology is full-resolution
          param_reg.warp_precision = 0.0;

          // Run the registration
          api_reg.RunDeformable(param_reg);
          }

        MultiComponentMetricReport last_metric_report = api_reg.GetLastMetricReport();
        total_to_vol_metric += last_metric_report.ComponentMetrics[0];
        for(unsigned int a = 1; a < last_metric_report.ComponentMetrics.size(); a++)
          total_to_nbr_metric += last_metric_report.ComponentMetrics[a];

        // Write the metric for this slide to file
        std::string fn_metric = GetFilenameForSlice(m_Slices[k], ITER_METRIC_DUMP, iter);
        std::ofstream fout(fn_metric);
        fout << api_reg.PrintIter(-1, -1, last_metric_report) << std::endl;

        // For deformable iterations, propagate the matrix from last iteration
        if(iter > n_affine)
          {
          GreedyAPI::WriteAffineMatrix(
                GetFilenameForSlice(m_Slices[k], VOL_ITER_MATRIX, iter),
                GreedyAPI::ReadAffineMatrix(
                  GetFilenameForSlice(m_Slices[k], VOL_ITER_MATRIX, iter-1)));
          }

        }

      printf("ITER %3d  TOTAL_VOL_METRIC = %8.4f  TOTAL_NBR_METRIC = %8.4f\n",
             iter, total_to_vol_metric, total_to_nbr_metric);
      }
  }

  void Splat(const SplatParameters &sparam, const GreedyParameters &gparam)
  {
    // The target volume into which we will be doing the splatting. It must either
    // be read from file or generated based on the 2D slices in the project
    LDDMMType3D::CompositeImagePointer target;

    // Use an image cache
    ImageCache icache(0, 20);

    // Before allocating the target, we need to know how many components to use. For
    // this we need to load the reference (root) slide
    unsigned int i_root = GetRootSlide();
    LDDMMType::CompositeImagePointer root_slide =
        icache.GetImage<LDDMMType::CompositeImageType>(m_Slices[i_root].raw_filename);

    // Set the number of components to that from the root slide. However, when using
    // an alternative manifest, we will get this from one of the input images instead
    unsigned int n_comp_out = root_slide->GetNumberOfComponentsPerPixel();

    // Determine which images and which ids to use for splatting.
    std::map<std::string, std::string> alt_source;
    if(sparam.fn_manifest.length())
      {
      std::ifstream fin(sparam.fn_manifest);
      std::string f_line;
      while(std::getline(fin, f_line))
        {
        std::istringstream iss(f_line);

        // Read an id from the manifest
        std::string id;
        if(!(iss >> id))
          throw GreedyException("Unable to read id from manifest file, line '%s'", f_line.c_str());

        // Find that id
        int index = FindSlideById(id);
        if(index < 0)
          throw GreedyException("Slide id '%s' is not in the project", id.c_str());

        // The user may request an alternative image to load
        std::string fn_alt_slice = m_Slices[index].raw_filename;
        if(iss >> fn_alt_slice)
          {
          if(!itksys::SystemTools::FileExists(fn_alt_slice.c_str(), true))
            throw GreedyException("File '%s' in manifest does not exist", fn_alt_slice.c_str());
          }

        // Insert the alternative slice source
        alt_source[id] = fn_alt_slice;
        }
      }

    // Check if the number of components should be updated
    if(alt_source.size())
      {
      std::string fn_alt = alt_source.begin()->second;
      n_comp_out = icache.GetImage<LDDMMType::CompositeImageType>(fn_alt)
                   ->GetNumberOfComponentsPerPixel();
      }

    // Was a referene volume specified?
    if(sparam.reference.length())
      {
      // If the reference volume is specified, use it
      LDDMMType3D::CompositeImagePointer ref = LDDMMType3D::cimg_read(sparam.reference.c_str());

      // Create a new reference volume
      if(ref->GetNumberOfComponentsPerPixel() == n_comp_out)
        target = ref;
      else
        target = LDDMMType3D::new_cimg(ref, (int) n_comp_out);
      }
    else
      {
      // When the reference volume is not specified, we will use 2D
      // slice metadata in the project as a reference for the x-y aspects of
      // the image, and have the user specify the spacing and origin in z.
      LDDMMType::ImageBaseType::Pointer ref_slide;

      // Are we reconstructing into slide space or volume space?
      if(sparam.source_stage == SplatParameters::RAW ||
         sparam.source_stage == SplatParameters::RECON)
        {
        // Read the reference slide
        ref_slide = root_slide;
        }
      else
        {
        // Read the reference volume slice
        std::string fn_slide_ref = GetFilenameForSlice(m_Slices[i_root], VOL_SLIDE);
        ref_slide = LDDMMType::img_read(fn_slide_ref.c_str());
        }

      // Create the 3D volume
      target = LDDMMType3D::CompositeImageType::New();

      // Set up the properties of the 3D volume
      LDDMMType3D::RegionType region_3d;
      LDDMMType3D::ImageType::PointType origin_3d;
      LDDMMType3D::ImageType::SpacingType spacing_3d;
      LDDMMType3D::ImageType::DirectionType dir_3d;

      dir_3d.SetIdentity();
      for(unsigned int a = 0; a < 2; a++)
        {
        region_3d.SetSize(a, ref_slide->GetBufferedRegion().GetSize(a));
        origin_3d[a] = ref_slide->GetOrigin()[a];
        spacing_3d[a] = ref_slide->GetSpacing()[a];
        for(unsigned int b = 0; b < 2; b++)
          dir_3d(a,b) = ref_slide->GetDirection()(a,b);
        }

      region_3d.SetSize(2, (unsigned long) ceil((sparam.z_last - sparam.z_first) / sparam.z_step));
      origin_3d[2] = sparam.z_first;
      spacing_3d[2] = sparam.z_step;

      target->SetRegions(region_3d);
      target->SetOrigin(origin_3d);
      target->SetSpacing(spacing_3d);
      target->SetDirection(dir_3d);
      target->SetNumberOfComponentsPerPixel(n_comp_out);
      target->Allocate();

      LDDMMType3D::CompositeImageType::PixelType cpix;
      cpix.SetSize(n_comp_out);
      cpix.Fill(gparam.current_interp.outside_value);
      target->FillBuffer(cpix);
      }


    // In addition to the target, we need a 2D reference image, which we will use
    // as the target for 2D reslice operations
    LDDMMType::CompositeImagePointer ref_2d = ExtractSliceFromVolume(target, target->GetOrigin()[2]);

    // In exact mode, we are not splatting, but rather just sampling along the z-axis.
    if(sparam.mode == SplatParameters::EXACT)
      {
      // We will iterate over the slices in the 3D volume
      typedef itk::ImageSliceIteratorWithIndex<LDDMMType3D::CompositeImageType> SliceIter;
      SliceIter it_slice(target, target->GetBufferedRegion());
      it_slice.SetFirstDirection(0);
      it_slice.SetSecondDirection(1);
      for(; !it_slice.IsAtEnd(); it_slice.NextSlice())
        {
        // Get the z position of this slide
        double z_pos = target->GetOrigin()[2] + target->GetSpacing()[2] * it_slice.GetIndex()[2];

        // Find the slice with tolerance
        int i_slice = FindSlideByZ(z_pos, sparam.z_exact_tol, alt_source);
        if(i_slice < 0)
          continue;

        // Get the filename to read
        std::string fn_source = alt_source.size() == 0
                                ? m_Slices[i_slice].raw_filename
                                : alt_source[m_Slices[i_slice].unique_id];

        // Read the source image
        LDDMMType::CompositeImagePointer img_source =
            icache.GetImage<LDDMMType::CompositeImageType>(fn_source);

        // Smooth the image if needed
        if(sparam.sigma_inplane > 0.0)
          {
          LDDMMType::Vec sigma_phys;
          for(unsigned int a = 0; a < 2; a++)
            sigma_phys[a] = sparam.sigma_inplane * img_source->GetSpacing()[a];
          LDDMMType::cimg_smooth(img_source, img_source, sigma_phys);
          }

        // Should we override the image header?
        if(fn_source != m_Slices[i_slice].raw_filename && sparam.ignore_alt_headers)
          {
          // Read the raw image
          LDDMMType::CompositeImagePointer img_raw =
              icache.GetImage<LDDMMType::CompositeImageType>(m_Slices[i_slice].raw_filename);

          // Adjust the header so that everything matches up
          LDDMMType::CompositeImageType::SpacingType spacing_adj;
          LDDMMType::CompositeImageType::PointType origin_adj;
          for(unsigned int a = 0; a < 2; a++)
            {
            double s1 = img_raw->GetSpacing()[a];
            double o1 = img_raw->GetOrigin()[a];
            unsigned long d1 = img_raw->GetBufferedRegion().GetSize()[a];
            unsigned long d2 = img_source->GetBufferedRegion().GetSize()[a];
            spacing_adj[a] = (d1 * s1) / d2;
            origin_adj[a] = o1 + 0.5 * (spacing_adj[a] - s1);
            }

          img_source->SetSpacing(spacing_adj);
          img_source->SetOrigin(origin_adj);
          img_source->SetDirection(img_raw->GetDirection());
          }

        // Check the number of components
        if(img_source->GetNumberOfComponentsPerPixel() != n_comp_out)
          throw GreedyException("Number of components in '%s' does not match %d",
                                fn_source.c_str(), n_comp_out);

        // Reslice into the target space
        GreedyAPI reslice_api;
        GreedyParameters my_param = gparam;

        // Print some progress
        printf("Splatting at z = %8.4f: %s\n", z_pos, fn_source.c_str());

        // Set the reslice specs
        my_param.reslice_param.ref_image = "ref";
        my_param.reslice_param.images.push_back(ResliceSpec("src", "out", my_param.current_interp));

        // Create an image to hold the output
        LDDMMType::CompositeImagePointer resliced = LDDMMType::CompositeImageType::New();

        // Add the cached images
        reslice_api.AddCachedInputObject("src", img_source.GetPointer());
        reslice_api.AddCachedInputObject("ref", ref_2d.GetPointer());
        reslice_api.AddCachedOutputObject("out", resliced.GetPointer());

        // Set up the transform chain
        if(sparam.source_stage == SplatParameters::RECON)
          {
          my_param.reslice_param.transforms.push_back(
                TransformSpec(GetFilenameForSlice(m_Slices[i_slice], ACCUM_MATRIX)));
          }
        else if(sparam.source_stage == SplatParameters::VOL_MATCH)
          {
          my_param.reslice_param.transforms.push_back(
                TransformSpec(GetFilenameForSlice(m_Slices[i_slice], VOL_ITER_MATRIX, 0)));
          }
        else if(sparam.source_stage == SplatParameters::VOL_ITER)
          {
          unsigned int n_affine = LoadConfigKey("AffineIterations", 0u);
          unsigned int n_deform = LoadConfigKey("DeformableIterations", 0u);

          if(sparam.source_iter < 1 || sparam.source_iter > n_affine + n_deform)
            throw GreedyException("Iteration parameter %d is out of range [1,%d]",
                                  sparam.source_iter, n_affine + n_deform);

          // Add the warp
          if(sparam.source_iter > n_affine)
            my_param.reslice_param.transforms.push_back(
                  TransformSpec(
                    GetFilenameForSlice(m_Slices[i_slice], VOL_ITER_WARP, sparam.source_iter)));

          // Add the affine matrix
          my_param.reslice_param.transforms.push_back(
                TransformSpec(
                  GetFilenameForSlice(m_Slices[i_slice], VOL_ITER_MATRIX, sparam.source_iter)));
          }

        // Run the reslice operation
        reslice_api.RunReslice(my_param);

        // Get an iterator into the result
        itk::ImageRegionConstIterator<LDDMMType::CompositeImageType> it(resliced, resliced->GetBufferedRegion());

        // Get the number of elements to copy
        unsigned long n_elts = resliced->GetPixelContainer()->Size();

        // Copy the pixels to destination. Some brute force pointer calculations here
        LDDMMType::CompositeImageType::InternalPixelType *p_src = resliced->GetBufferPointer();
        LDDMMType::CompositeImageType::InternalPixelType *p_trg =
            target->GetBufferPointer() + it_slice.GetIndex()[2] * n_elts;

        for(unsigned long i = 0; i < n_elts; i++)
          p_trg[i] = p_src[i];
        }
      }
    else throw GreedyException("Only exact mode is implemented.");

    // Write the image
    LDDMMType3D::cimg_write(target, sparam.fn_output.c_str());
  }

  int FindSlideById(const std::string &id)
  {
    for (int i = 0; i < m_Slices.size(); i++)
      {
      if(m_Slices[i].unique_id == id)
        return i;
      }

    return -1;
  }

  int FindSlideByZ(double z, double z_tol, const std::map<std::string, std::string> &alt_source)
  {
    int i_best = -1;
    double dz_best = 1e100;

    for (int i = 0; i < m_Slices.size(); i++)
      {
      // Check the filter
      if(alt_source.size() && alt_source.find(m_Slices[i].unique_id) == alt_source.end())
        continue;

      // Check the range
      if(m_Slices[i].z_pos >= z - z_tol && m_Slices[i].z_pos <= z + z_tol)
        {
        if(i_best < 0 || fabs(m_Slices[i].z_pos - z) < dz_best)
          {
          i_best = i;
          dz_best = fabs(m_Slices[i].z_pos - z);
          }
        }
      }

    return i_best;
  }

  unsigned int GetRootSlide()
  {
    // Determine the root slide
    std::string root_id = LoadConfigKey("RootSlide", std::string());
    if(root_id.length() == 0)
      throw GreedyException("Root slide has not been computed. Please run 'recon' first.");

    // Find the slide with that id
    int i_root = FindSlideById(root_id);
    if(i_root < 0)
      throw GreedyException("Root slide %s not found in manifest. Rerun 'recon'.", root_id.c_str());

    return (unsigned int) i_root;
  }

private:

  // Data associated with each slice, from the manifest
  struct SliceData
  {
    std::string raw_filename;
    std::string unique_id;
    double z_pos;
    bool is_leader;
  };

  // Path to the project
  std::string m_ProjectDir;

  // Default image file extension
  std::string m_DefaultImageExt;

  // Global parameters (parameters for the current run)
  StackParameters m_GlobalParam;

  // A flat list of slices (in manifest order)
  std::vector<SliceData> m_Slices;

  // A list of slices sorted by the z-position
  typedef std::pair<double, unsigned int> slice_ref;
  typedef std::set<slice_ref> slice_ref_set;
  slice_ref_set m_SortedSlices;

  std::string GetFilenameForSlicePair(
      const SliceData &ref, const SliceData &mov, FileIntent intent)
  {
    char filename[1024];
    const char *dir = m_ProjectDir.c_str(), *ext = m_DefaultImageExt.c_str();
    const char *rid = ref.unique_id.c_str(), *mid = mov.unique_id.c_str();

    switch(intent)
      {
      case AFFINE_MATRIX:
        sprintf(filename, "%s/recon/nbr/affine_ref_%s_mov_%s.mat", dir, rid, mid);
        break;
      case METRIC_VALUE:
        sprintf(filename, "%s/recon/nbr/affine_ref_%s_mov_%s_metric.txt", dir, rid, mid);
        break;
      default:
        throw GreedyException("Wrong intent in GetFilenameForSlicePair");
      }

    // Make sure the directory containing this exists
    itksys::SystemTools::MakeDirectory(itksys::SystemTools::GetFilenamePath(filename));

    return filename;
  }

  std::string GetFilenameForSlice(const SliceData &slice, int intent, ...)
  {
    char filename[1024];
    const char *dir = m_ProjectDir.c_str(), *ext = m_DefaultImageExt.c_str();
    const char *sid = slice.unique_id.c_str();

    va_list args;
    va_start(args, intent);

    int iter =
        (intent == VOL_ITER_MATRIX || intent == VOL_ITER_WARP || intent == ITER_METRIC_DUMP)
        ? va_arg(args, int) : 0;

    switch(intent)
      {
      case ACCUM_MATRIX:
        sprintf(filename, "%s/recon/accum/accum_affine_%s.mat", dir, sid);
        break;
      case ACCUM_RESLICE:
        sprintf(filename, "%s/recon/accum/accum_affine_%s_reslice.%s", dir, sid, ext);
        break;
      case VOL_INIT_MATRIX:
        sprintf(filename, "%s/vol/match/affine_refvol_mov_%s.mat", dir, sid);
        break;
      case VOL_SLIDE:
        sprintf(filename, "%s/vol/slides/vol_slide_%s.%s", dir, sid, ext);
        break;
      case VOL_MASK_SLIDE:
        sprintf(filename, "%s/vol/slides/vol_mask_slide_%s.%s", dir, sid, ext);
        break;
      case VOL_ITER_MATRIX:
        sprintf(filename, "%s/vol/iter%02d/affine_refvol_mov_%s_iter%02d.mat", dir, iter, sid, iter);
        break;
      case VOL_ITER_WARP:
        sprintf(filename, "%s/vol/iter%02d/warp_refvol_mov_%s_iter%02d.%s", dir, iter, sid, iter, ext);
        break;
      case ITER_METRIC_DUMP:
        sprintf(filename, "%s/vol/iter%02d/metric_refvol_mov_%s_iter%02d.txt", dir, iter, sid, iter);
        break;
      default:
        throw GreedyException("Wrong intent in GetFilenameForSlice");
      }

    va_end(args);

    // Make sure the directory containing this exists
    itksys::SystemTools::MakeDirectory(itksys::SystemTools::GetFilenamePath(filename));

    return filename;
  }

  std::string GetFilenameForGlobal(int intent, ...)
  {
    char filename[1024];
    const char *dir = m_ProjectDir.c_str(), *ext = m_DefaultImageExt.c_str();

    va_list args;
    va_start(args, intent);

    switch(intent)
      {
      case VOL_BEST_INIT_MATRIX:
        sprintf(filename, "%s/vol/match/affine_refvol_best.mat", dir);
        break;
      case MANIFEST_FILE:
        sprintf(filename, "%s/config/manifest.txt", dir);
        break;
      case CONFIG_ENTRY:
        sprintf(filename, "%s/config/dict/%s", dir, va_arg(args, char *));
        break;
      default:
        throw GreedyException("Wrong intent in GetFilenameForGlobal");
      }

    va_end(args);

    // Make sure the directory containing this exists
    itksys::SystemTools::MakeDirectory(itksys::SystemTools::GetFilenamePath(filename));

    return filename;
  }
};


/**
 * Initialize the project
 */
void init(StackParameters &param, CommandLineHelper &cl)
{
  // Parse the parameters
  std::string arg;
  std::string fn_manifest;
  std::string default_ext = "nii.gz";
  while(cl.read_command(arg))
    {
    if(arg == "-M")
      {
      fn_manifest = cl.read_existing_filename();
      }
    else if(arg == "-ext")
      {
      default_ext = cl.read_string();
      }
    else
      throw GreedyException("Unknown parameter to 'init': %s", arg.c_str());
    }

  // Check required parameters
  if(fn_manifest.size() == 0)
    throw GreedyException("Missing manifest file (-M) in 'init'");

  // Create the project
  StackGreedyProject sgp(param.output_dir, param);
  sgp.InitializeProject(fn_manifest, default_ext);
}

/**
 * Run the reconstruction module
 */
void recon(StackParameters &param, CommandLineHelper &cl)
{
  // List of greedy commands that are recognized by this mode
  std::set<std::string> greedy_cmd {
    "-m", "-n", "-threads", "-gm-trim", "-search"
  };

  // Greedy parameters for this mode
  GreedyParameters gparam;

  // Parse the parameters
  double z_range = 0.0;
  double z_epsilon = 0.1;
  std::string arg;
  while(cl.read_command(arg))
    {
    if(arg == "-z")
      {
      z_range = cl.read_double();
      z_epsilon = cl.read_double();
      }
    else if(greedy_cmd.find(arg) != greedy_cmd.end())
      {
      gparam.ParseCommandLine(arg, cl);
      }
    else
      throw GreedyException("Unknown parameter to 'init': %s", arg.c_str());
    }

  // Configure the threads
  GreedyApproach<2,double>::ConfigThreads(gparam);

  // Create the project
  StackGreedyProject sgp(param.output_dir, param);
  sgp.RestoreProject();
  sgp.ReconstructStack(z_range, z_epsilon, gparam);
}


/**
 * Run the volume matching module
 */
void volmatch(StackParameters &param, CommandLineHelper &cl)
{
  // List of greedy commands that are recognized by this mode
  std::set<std::string> greedy_cmd {
    "-m", "-n", "-threads", "-gm-trim", "-search"
  };

  // Greedy parameters for this mode
  GreedyParameters gparam;

  // Parse the parameters
  std::string fn_volume, fn_mask;
  std::string arg;
  while(cl.read_command(arg))
    {
    if(arg == "-i")
      {
      fn_volume = cl.read_existing_filename();
      }
    else if(arg == "-gm")
      {
      fn_mask = cl.read_existing_filename();
      }
    else if(greedy_cmd.find(arg) != greedy_cmd.end())
      {
      gparam.ParseCommandLine(arg, cl);
      }
    else
      throw GreedyException("Unknown parameter to 'volmatch': %s", arg.c_str());
    }

  // Check required parameters
  if(fn_volume.size() == 0)
    throw GreedyException("Missing volume file (-i) in 'volmatch'");

  // Configure the threads
  GreedyApproach<2,double>::ConfigThreads(gparam);

  // Create the project
  StackGreedyProject sgp(param.output_dir, param);
  sgp.RestoreProject();
  sgp.InitialMatchToVolume(fn_volume, fn_mask, gparam);
}


/**
 * Run the iterative module
 */
void voliter(StackParameters &param, CommandLineHelper &cl)
{
  // List of greedy commands that are recognized by this mode
  std::set<std::string> greedy_cmd {
    "-m", "-n", "-threads", "-gm-trim", "-s", "-e", "-sv", "-exp"
  };

  // Greedy parameters for this mode
  GreedyParameters gparam;

  // Parse the parameters
  unsigned int n_affine = 5, n_deform = 5;
  unsigned int i_first = 0, i_last = 0;
  double w_volume = 4.0;

  std::string arg;
  while(cl.read_command(arg))
    {
    if(arg == "-R")
      {
      i_first = cl.read_integer();
      i_last = cl.read_integer();
      }
    else if(arg == "-na")
      {
      n_affine = cl.read_integer();
      }
    else if(arg == "-nd")
      {
      n_deform = cl.read_integer();
      }
    else if(arg == "-w")
      {
      w_volume = cl.read_double();
      }
    else if(greedy_cmd.find(arg) != greedy_cmd.end())
      {
      gparam.ParseCommandLine(arg, cl);
      }
    else
      throw GreedyException("Unknown parameter to 'voliter': %s", arg.c_str());
    }

  // Default is to run all iterations
  if(i_first == 0 && i_last == 0)
    {
    i_first = 1;
    i_last = n_affine + n_deform;
    }

  // Configure the threads
  GreedyApproach<2,double>::ConfigThreads(gparam);

  // Create the project
  StackGreedyProject sgp(param.output_dir, param);
  sgp.RestoreProject();
  sgp.IterativeMatchToVolume(n_affine, n_deform, i_first, i_last, w_volume, gparam);
}


/**
 * Run the splatting module
 */
void splat(StackParameters &param, CommandLineHelper &cl)
{
  // List of greedy commands that are recognized by this mode
  std::set<std::string> greedy_cmd {
    "-threads", "-rb", "-ri"
  };

  // Greedy parameters for this mode
  GreedyParameters gparam;

  // Splatting parameters
  SplatParameters sparam;

  // Parse the parameters
  std::string arg;
  while(cl.read_command(arg))
    {
    if(arg == "-o")
      {
      sparam.fn_output = cl.read_output_filename();
      }
    else if(arg == "-i")
      {
      std::string mode = cl.read_string();
      if(mode == "raw")
        sparam.source_stage = SplatParameters::RAW;
      else if(mode == "recon")
        sparam.source_stage = SplatParameters::RECON;
      else if(mode == "volmatch")
        sparam.source_stage = SplatParameters::VOL_MATCH;
      else if(mode == "voliter")
        {
        sparam.source_stage = SplatParameters::VOL_ITER;
        if(cl.command_arg_count() > 0)
          sparam.source_iter = (unsigned int) cl.read_integer();
        }
      else throw GreedyException("Unknown stage specification %s", mode.c_str());
      }
    else if(arg == "-rf")
      {
      sparam.reference = cl.read_existing_filename();
      }
    else if(arg == "-z")
      {
      sparam.z_first = cl.read_double();
      sparam.z_step = cl.read_double();
      sparam.z_last = cl.read_double();
      }
    else if(arg == "-ztol")
      {
      sparam.z_exact_tol = cl.read_double();
      }
    else if(arg == "-S")
      {
      std::string mode = cl.read_string();
      if(mode == "exact")
        sparam.mode = SplatParameters::EXACT;
      else if(mode == "nearest")
        sparam.mode = SplatParameters::NEAREST;
      else if(mode == "linear")
        sparam.mode = SplatParameters::LINEAR;
      else if(mode == "parzen")
        sparam.mode = SplatParameters::PARZEN;
      }
    else if(arg == "-M")
      {
      sparam.fn_manifest = cl.read_existing_filename();
      }
    else if(arg == "-H")
      {
      sparam.ignore_alt_headers = true;
      }
    else if(arg == "-si")
      {
      sparam.sigma_inplane = cl.read_double();
      }
    else if(greedy_cmd.find(arg) != greedy_cmd.end())
      {
      gparam.ParseCommandLine(arg, cl);
      }
    else
      throw GreedyException("Unknown parameter to 'splat': %s", arg.c_str());
    }

  // Configure the threads
  GreedyApproach<2,double>::ConfigThreads(gparam);

  // Create the project
  StackGreedyProject sgp(param.output_dir, param);
  sgp.RestoreProject();
  sgp.Splat(sparam, gparam);
}

int main(int argc, char *argv[])
{
  // Parameters specifically for this application
  StackParameters param;

  // Parameters for running Greedy in general
  GreedyParameters gparam;
  GreedyParameters::SetToDefaults(gparam);
  if(argc < 2)
    return usage();

  // Read the first command
  try
  {
  CommandLineHelper cl(argc, argv);

  // Read the global commands
  while(!cl.is_at_end() && cl.peek_arg()[0] == '-')
    {
    std::string arg = cl.read_command();
    if(arg == "-N")
      {
      param.reuse = true;
      }
    else
      {
      std::cerr << "Unknown global option " << arg << std::endl;
      return -1;
      }
    }

  // Read the main command
  if(cl.is_at_end())
    {
    std::cerr << "Missing command. Run this program without parameters to see usage." << std::endl;
    return -1;
    }

  std::string cmd = cl.read_string();

  // Is the string known
  if(cmd == "help")
    {
    return usage(cl.is_at_end() ? std::string() : cl.read_string());
    }

  // All commands other than 'help' end with the project directory. So we should get that
  // as the last argument from the command-line
  CommandLineHelper cl_end = cl.take_end(1);
  param.output_dir = cl_end.read_output_filename();

  if(cmd == "init")
    {
    init(param, cl);
    }
  else if(cmd == "recon")
    {
    recon(param, cl);
    }
  else if(cmd == "volmatch")
    {
    volmatch(param, cl);
    }
  else if(cmd == "voliter")
    {
    voliter(param, cl);
    }
  else if(cmd == "splat")
    {
    splat(param, cl);
    }
  else
    {
    std::cerr << "Unknown command " << cmd << std::endl;
    return -1;
    }
  }
  catch(std::exception &exc)
  {
    std::cerr << "ERROR: exception thrown in the code:" << std::endl;
    std::cerr << exc.what() << std::endl;
    return -1;
  }

  return 0;
}


