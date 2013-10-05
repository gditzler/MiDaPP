#!/usr/bin/env python

import numpy as np
from optparse import OptionParser
import cogent.cluster.metric_scaling as ms
from biom.parse import parse_biom_table
from qiime.parse import parse_mapping_file_to_dict
import cogent.maths.distance_transform as distance_transform

site_column = "HMPbodysubsite"

def hmp_pcoa(biom_path, map_path, distance="hellinger"):
  """
  @biom_path
  @map_path
  @distance
  """
  data,labn,labs,classes = load_data(biom_path, map_path)

  dist_mtrx_fcn = getattr(distance_transform, 'dist_'+distance)
  dist_mtrx = dist_mtrx_fcn(data)
  coords, eigvals=ms.principal_coordinates_analysis(dist_mtrx)
  
  returgvals = np.abs(eigvals)
  pcnts = (np.abs(eigvals)/sum(np.abs(eigvals)))*100
  idxs_descending = pcnts.argsort()[::-1]
  coords = coords[idxs_descending]
  eigvals = eigvals[idxs_descending]
  pcnts = pcnts[idxs_descending]
  return None 

def load_data(biom_path, map_path):
  """
  @biom_path
  @map_path
  @data
  @labn
  @labs
  @classes
  """
  obj,comm = parse_mapping_file_to_dict(open(map_path, "U"))
  labn,labs,classes = extract_labels(obj)
  data = extract_data(biom_path)
  return data, labn, labs, classes
 
def extract_data(biom_path):
  """
  @biom_path
  @data
  """
  biom_table = parse_biom_table(open(biom_path, "U"))
  data = []
  for x in biom_table.iterObservationData():
    data.append(x)
  return np.array(data).transpose()

def extract_labels(obj):
  """
  @obj
  @labels_numeric
  @labels_string
  @classes
  """
  classes = {}
  n = 0
  labels_numeric = []
  labels_string = [] 
  for id_name in obj.keys(): 
    labels_string.append(obj[id_name]["HMPbodysubsite"])
    if not classes.has_key(obj[id_name]["HMPbodysubsite"]):
      classes[obj[id_name]["HMPbodysubsite"]] = n
      n += 1 
    labels_numeric.append(classes[obj[id_name]["HMPbodysubsite"]])
  return labels_numeric, labels_string, classes

if __name__ == "__main__": 
  """
  """
  parser = OptionParser()
  parser.add_option("-a", "--biom",
    dest = "biom",
    help = "path to the biom file"
  )
  parser.add_option("-b", "--map",
    dest = "map",
    help = "path to the map file"
  )
  (opts, args) = parser.parse_args()
  hmp_pcoa(opts.biom, opts.map)

