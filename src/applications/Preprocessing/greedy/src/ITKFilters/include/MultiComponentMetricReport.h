#ifndef MultiComponentMetricReport_H
#define MultiComponentMetricReport_H

/**
 * A structure used to communicate metric information. Meant to be an atomic structure
 * that can be expanded with additional bits of information in the future
 */
struct MultiComponentMetricReport
{
  double TotalMetric;
  vnl_vector<double> ComponentMetrics;

  void Scale(double scale_factor)
  {
    TotalMetric *= scale_factor;
    ComponentMetrics *= scale_factor;
  }
};


#endif // MultiComponentMetricReport_H
