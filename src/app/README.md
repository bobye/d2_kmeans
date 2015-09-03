## Usage

### Arguments
Input options
 - `--ifile <input_filename>, -i <input_filename>` : the main name of input file (required).
 - `--ofile <output_filename>, -o <output_filename>` : the main name of output file (optional), when it is default, the output_filename will be generated based on input filename and a random number.
 - `--metafile <meta_filename>, -D <meta_filename>` : you can optionally specify the meta filename (the .hist or .vacab file), when there is only one phase of interest. 
 - `-n <integer>` : number of instances read in per processor (required)
 - `--phase <integer>, -p <integer>` : the number of phases per instance in `<input_filename>` (default: 1).
 - `--phase_only <integer>, -t <integer>` : the phase to cluster upon (default: use all phases)
 - `--strides <integer array>, -s <integer array>` : the numbers of support points of computed centroids in each phase (required), integer array with comma delimiter and no spaces. 
 - `-d <integer array>` : the dimensions in each phase (required), integer array with comma delimiter and no spaces.
 - `--types <integer>, -E <integer>` : the type of D2 data (default: 0, see `include/d2_param.h` for details).
 
Algorithm options
 - `--clusters <integer>, -k <integer>` : number of clusters intended (default: 3, mostly required). If it is set to 1, the centroid of data is computed instead, which takes more ADMM steps (2000 steps) than that of clustering setting (100 steps). 
 - `--max_iters <integer>, -m <integer>` : the maximal number of iterations (default: 100).
 - `--non_triangle, -T` : disable the triangle inequality based acceleration (default: enabled).
 - `--eval <centroids_filename>, -e <centroids_filename>` : no clustering, but assigning instances to the pre-computed centroids (default: disabled).
 - `--load <centroids_filename>, -L <centroids_filename>` : load pre-computed centroids as initial start of D2 clustering. (optional, excluding the `--eval` option)
 - `--pre_process, -Q` : preprocess the input format data (more explanations TBA) (default: disabled).
 
Parallel computing options
 - `--prepare_batches <integer>, -P <integer>` : the number of batches that is equal to the number of processors in data pre-processing stage. It reads in `<input_filename> = data.d2` and generates files in say `data.d2.0, data.d2.1, data.d2.2, data.d2.3` containing randomly splitted parts of `data.d2`.
 
### Modes
1. The default mode is the clustering algorithm, which outputs results in two main files. For example, taking in `data.d2`, program outputs the centroids computed (`data.d2_123456_c.d2` which is again in D2 format) and the memberships of each instances (`data.d2_123456.label` in sequential run or `data.d2_123456.label_o` in parallel run). 
2. To preprocess a data file into multiple batches and later feed them into a parallel computing environment, one has to call `--prepare_batches`.
3. Given pre-computed centroids from a training set, one can assign cluster memberships to another testing set using `--eval`. 
