/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2017,2018 The plumed team
   (see the PEOPLE file at the root of the distribution for a list of names)

   See http://www.plumed.org for more information.

   This file is part of plumed, version 2.

   plumed is free software: you can redistribute it and/or modify
   it under the terms of the GNU Lesser General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

   plumed is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU Lesser General Public License for more details.

   You should have received a copy of the GNU Lesser General Public License
   along with plumed.  If not, see <http://www.gnu.org/licenses/>.
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
#include "colvar/Colvar.h"
#include "colvar/ActionRegister.h"
#include "core/PlumedMain.h"
#include "tools/Matrix.h"
#include "core/SetupMolInfo.h"
#include "core/ActionSet.h"
#include "tools/File.h"

#include <string>
#include <cmath>
#include <map>
#include <numeric>
#include <ctime>
#include "tools/Random.h"

using namespace std;

namespace PLMD {
namespace isdb {

//+PLUMEDOC ISDB_COLVAR EMMIVOX
/*
Calculate the fit of a structure or ensemble of structures with a cryo-EM density map.

This action implements the multi-scale Bayesian approach to cryo-EM data fitting introduced in  Ref. \cite Hanot113951 .
This method allows efficient and accurate structural modeling of cryo-electron microscopy density maps at multiple scales, from coarse-grained to atomistic resolution, by addressing the presence of random and systematic errors in the data, sample heterogeneity, data correlation, and noise correlation.

The experimental density map is fit by a Gaussian Mixture Model (GMM), which is provided as an external file specified by the keyword
GMM_FILE. We are currently working on a web server to perform
this operation. In the meantime, the user can request a stand-alone version of the GMM code at massimiliano.bonomi_AT_gmail.com.

When run in single-replica mode, this action allows atomistic, flexible refinement of an individual structure into a density map.
Combined with a multi-replica framework (such as the -multi option in GROMACS), the user can model an esemble of structures using
the Metainference approach \cite Bonomi:2016ip .

\warning
To use \ref EMMIVOX, the user should always add a \ref MOLINFO line and specify a pdb file of the system.

\note
To enhance sampling in single-structure refinement, one can use a Replica Exchange Method, such as Parallel Tempering.
In this case, the user should add the NO_AVER flag to the input line.

\note
\ref EMMIVOX can be used in combination with periodic and non-periodic systems. In the latter case, one should
add the NOPBC flag to the input line

\par Examples

In this example, we perform a single-structure refinement based on an experimental cryo-EM map. The map is fit with a GMM, whose
parameters are listed in the file GMM_fit.dat. This file contains one line per GMM component in the following format:

\plumedfile
#! FIELDS Id Weight Mean_0 Mean_1 Mean_2 Cov_00 Cov_01 Cov_02 Cov_11 Cov_12 Cov_22 Beta
     0  2.9993805e+01   6.54628 10.37820 -0.92988  2.078920e-02 1.216254e-03 5.990827e-04 2.556246e-02 8.411835e-03 2.486254e-02  1
     1  2.3468312e+01   6.56095 10.34790 -0.87808  1.879859e-02 6.636049e-03 3.682865e-04 3.194490e-02 1.750524e-03 3.017100e-02  1
     ...
\endplumedfile

To accelerate the computation of the Bayesian score, one can:
- use neighbor lists, specified by the keywords NL_CUTOFF and NL_STRIDE;
- calculate the restraint every other step (or more).

All the heavy atoms of the system are used to calculate the density map. This list can conveniently be provided
using a GROMACS index file.

The input file looks as follows:

\plumedfile
# include pdb info
MOLINFO STRUCTURE=prot.pdb

#  all heavy atoms
protein-h: GROUP NDX_FILE=index.ndx NDX_GROUP=Protein-H

# create EMMIVOX score
gmm: EMMIVOX NOPBC SIGMA_MIN=0.01 TEMP=300.0 NL_STRIDE=100 NL_CUTOFF=0.01 GMM_FILE=GMM_fit.dat ATOMS=protein-h

# translate into bias - apply every 2 steps
emr: BIASVALUE ARG=gmm.scoreb STRIDE=2

PRINT ARG=emr.* FILE=COLVAR STRIDE=500 FMT=%20.10f
\endplumedfile


*/
//+ENDPLUMEDOC

class EMMIVOX : public Colvar {

private:

// temperature in kbt
  double kbt_;
// model GMM - atom types
  vector<unsigned> GMM_m_type_;
// model GMM - list of atom sigmas - one per atom type
  vector<double> GMM_m_s_;
// model GMM - list of atom weights - one per atom type
  vector<double> GMM_m_w_;
// model GMM - map between residue number and list of atoms
  map< unsigned, vector<unsigned> > GMM_m_resmap_;
// model GMM - list of residue ids
  vector<unsigned> GMM_m_res_;
// model GMM - list of neighboring voxels per atom
  vector< vector<unsigned> > GMM_m_nb_;
// model GMM - map between res id and bfactor
  map<unsigned,double> GMM_m_b_;
// model overlap
  vector<double> ovmd_;

// data GMM - means, sigma2 + beta option
  vector<Vector> GMM_d_m_;
  vector<Vector> GMM_d_s_;
  vector<int>    GMM_d_beta_;
// data GMM - groups bookeeping
  vector < vector<int> > GMM_d_grps_;
// model GMM - list of neighboring atoms per voxel
  vector< vector<unsigned> > GMM_d_nb_;
// data GMM - overlap
  vector<double> ovdd_;

// derivatives
  vector<Vector> ovmd_der_;
  vector<Vector> atom_der_;
  vector<double> GMMid_der_;
// constants
  vector<double> cfact_;
  double inv_sqrt2_, sqrt2_pi_, inv_pi2_;
// metainference
  unsigned nrep_;
  unsigned replica_;
  vector<double> sigma_;
  vector<double> sigma_min_;
  vector<double> sigma_max_;
  vector<double> dsigma_;
// neighbor list
  double   nl_cutoff_;
  unsigned nl_stride_;
  bool first_time_;
  vector<unsigned> nl_;
// averaging
  bool no_aver_;
// parallel stuff
  unsigned size_;
  unsigned rank_;
// pbc
  bool pbc_;
// Monte Carlo stuff
  int      MCstride_;
  double   MCaccept_;
  double   MCtrials_;
  Random   random_;
// Bfact Monte Carlo
  int      MCBstride_;
  double   MCBaccept_;
  double   MCBtrials_;
  double   dbfact_;
  double   bfactmin_;
  double   bfactmax_;
  // status stuff
  unsigned int statusstride_;
  string       statusfilename_;
  OFile        statusfile_;
  bool         first_status_;
  // regression
  unsigned nregres_;
  double scale_;
  double scale_min_;
  double scale_max_;
  double dscale_;
  // tabulated exponential
  double dpcutoff_;
  double dexp_;
  unsigned nexp_;
  vector<double> tab_exp_;
  // simulated annealing
  unsigned nanneal_;
  double   kanneal_;
  double   anneal_;
  // prior exponent
  double prior_;
  // noise type
  unsigned noise_;
  // total score and virial;
  double ene_;
  Tensor virial_;
  // model overlap file
  unsigned int ovstride_;
  string       ovfilename_;

// write file with model overlap
  void write_model_overlap(long int step);
// get median of vector
  double get_median(vector<double> &v);
// annealing
  double get_annealing(long int step);
// do regression
  double scaleEnergy(double s);
  double doRegression();
// read and write status
  void read_status();
  void print_status(long int step);
// accept or reject
  bool doAccept(double oldE, double newE, double kbt);
// do MonteCarlo
  void doMonteCarlo();
// do MonteCarlo for Bfactor
  void doMonteCarloBfact();
// read error file
  vector<double> read_exp_errors(string errfile);
// calculate model GMM parameters
  vector<double> get_GMM_m(vector<AtomNumber> &atoms);
// read data file
  void get_exp_data(string datafile);
// auxiliary method
  void calculate_useful_stuff(double reso);
// calculate overlap between two Gaussians
  double get_overlap(const Vector &d_m, const Vector &m_m, const Vector &d_s,
                     double pref, double m_b, Vector &ov_der);
// update the neighbor list
  void update_neighbor_list();
// calculate overlap
  void calculate_overlap();
// Gaussian noise
  void calculate_Gauss();
// Outliers noise
  void calculate_Outliers();
// Marginal noise
  void calculate_Marginal();

public:
  static void registerKeywords( Keywords& keys );
  explicit EMMIVOX(const ActionOptions&);
// active methods:
  void prepare();
  virtual void calculate();
};

PLUMED_REGISTER_ACTION(EMMIVOX,"EMMIVOX")

void EMMIVOX::registerKeywords( Keywords& keys ) {
  Colvar::registerKeywords( keys );
  keys.add("atoms","ATOMS","atoms for which we calculate the density map, typically all heavy atoms");
  keys.add("compulsory","DATA_FILE","file with the experimental data");
  keys.add("compulsory","NL_CUTOFF","The cutoff in overlap for the neighbor list");
  keys.add("compulsory","NL_STRIDE","The frequency with which we are updating the neighbor list");
  keys.add("compulsory","SIGMA_MIN","minimum uncertainty");
  keys.add("compulsory","RESOLUTION", "Cryo-EM map resolution");
  keys.add("compulsory","NOISETYPE","functional form of the noise (GAUSS, OUTLIERS, MARGINAL)");
  keys.add("compulsory","NORM_DENSITY","integral of the experimental density");
  keys.add("compulsory","WRITE_STRIDE","write the status to a file every N steps, this can be used for restart");
  keys.add("optional","SIGMA0","initial value of the uncertainty");
  keys.add("optional","DSIGMA","MC step for uncertainties");
  keys.add("optional","MC_STRIDE", "Monte Carlo stride");
  keys.add("optional","DBFACT","MC step for bfactor");
  keys.add("optional","BFACT_MAX","Maximum value of bfactor");
  keys.add("optional","MCBFACT_STRIDE", "Bfactor Monte Carlo stride");
  keys.add("optional","ERR_FILE","file with experimental errors");
  keys.add("optional","STATUS_FILE","write a file with all the data useful for restart");
  keys.add("optional","REGRESSION","regression stride");
  keys.add("optional","REG_SCALE_MIN","regression minimum scale");
  keys.add("optional","REG_SCALE_MAX","regression maximum scale");
  keys.add("optional","REG_DSCALE","regression maximum scale MC move");
  keys.add("optional","SCALE","scale factor");
  keys.add("optional","ANNEAL", "Length of annealing cycle");
  keys.add("optional","ANNEAL_FACT", "Annealing temperature factor");
  keys.add("optional","TEMP","temperature");
  keys.add("optional","PRIOR", "exponent of uncertainty prior");
  keys.add("optional","WRITE_OV_STRIDE","write model overlaps every N steps");
  keys.add("optional","WRITE_OV","write a file with model overlaps");
  keys.addFlag("NO_AVER",false,"don't do ensemble averaging in multi-replica mode");
  componentsAreNotOptional(keys);
  keys.addOutputComponent("scoreb","default","Bayesian score");
  keys.addOutputComponent("acc",   "NOISETYPE","MC acceptance for uncertainty");
  keys.addOutputComponent("accB",  "default", "Bfactor MC acceptance");
  keys.addOutputComponent("scale", "REGRESSION","scale factor");
  keys.addOutputComponent("accscale", "REGRESSION","MC acceptance for scale regression");
  keys.addOutputComponent("enescale", "REGRESSION","MC energy for scale regression");
  keys.addOutputComponent("anneal","ANNEAL","annealing factor");
}

EMMIVOX::EMMIVOX(const ActionOptions&ao):
  PLUMED_COLVAR_INIT(ao),
  inv_sqrt2_(0.707106781186548),
  sqrt2_pi_(0.797884560802865),
  inv_pi2_(0.050660591821169),
  first_time_(true), no_aver_(false), pbc_(true),
  MCstride_(1), MCaccept_(0.), MCtrials_(0.),
  MCBstride_(1), MCBaccept_(0.), MCBtrials_(0.),
  dbfact_(0.0), bfactmax_(4.0),
  statusstride_(0), first_status_(true),
  nregres_(0), scale_(1.),
  dpcutoff_(15.0), nexp_(1000000), nanneal_(0),
  kanneal_(0.), anneal_(1.), prior_(1.), ovstride_(0)
{
  // periodic boundary conditions
  bool nopbc=!pbc_;
  parseFlag("NOPBC",nopbc);
  pbc_=!nopbc;

  // list of atoms
  vector<AtomNumber> atoms;
  parseAtomList("ATOMS", atoms);

  // file with experimental data
  string datafile;
  parse("DATA_FILE", datafile);

  // type of data noise
  string noise;
  parse("NOISETYPE",noise);
  if      (noise=="GAUSS")   noise_ = 0;
  else if(noise=="OUTLIERS") noise_ = 1;
  else if(noise=="MARGINAL") noise_ = 2;
  else error("Unknown noise type!");

  // minimum value for error
  double sigma_min;
  parse("SIGMA_MIN", sigma_min);
  if(sigma_min<0) error("SIGMA_MIN should be greater or equal to zero");

  // Monte Carlo in B-factors
  parse("DBFACT", dbfact_);
  parse("BFACT_MAX", bfactmax_);
  parse("MCBFACT_STRIDE", MCBstride_);
  if(dbfact_<0) error("DBFACT should be greater or equal to zero");
  if(dbfact_>0 && MCBstride_<=0) error("you must specify a positive MCBFACT_STRIDE");
  if(dbfact_>0 && bfactmax_<=0)  error("you must specify a positive BFACT_MAX");


  // the following parameters must be specified with noise type 0 and 1
  double sigma_ini, dsigma;
  if(noise_!=2) {
    // initial value of the uncertainty
    parse("SIGMA0", sigma_ini);
    if(sigma_ini<=0) error("you must specify a positive SIGMA0");
    // MC parameters
    parse("DSIGMA", dsigma);
    if(dsigma<0) error("you must specify a positive DSIGMA");
    parse("MC_STRIDE", MCstride_);
    if(dsigma>0 && MCstride_<=0) error("you must specify a positive MC_STRIDE");
  }

  // status file parameters
  parse("WRITE_STRIDE", statusstride_);
  if(statusstride_<=0) error("you must specify a positive WRITE_STRIDE");
  parse("STATUS_FILE",  statusfilename_);
  if(statusfilename_=="") statusfilename_ = "MISTATUS"+getLabel();
  else                    statusfilename_ = statusfilename_+getLabel();

  // error file
  string errfile;
  parse("ERR_FILE", errfile);

  // integral of the experimetal density
  double norm_d;
  parse("NORM_DENSITY", norm_d);

  // temperature
  double temp=0.0;
  parse("TEMP",temp);
  // convert temp to kbt
  if(temp>0.0) kbt_=plumed.getAtoms().getKBoltzmann()*temp;
  else kbt_=plumed.getAtoms().getKbT();

  // exponent of uncertainty prior
  parse("PRIOR",prior_);

  // simulated annealing stuff
  parse("ANNEAL", nanneal_);
  parse("ANNEAL_FACT", kanneal_);
  if(nanneal_>0 && kanneal_<=1.0) error("with ANNEAL, ANNEAL_FACT must be greater than 1");

  // regression stride
  parse("REGRESSION",nregres_);
  // other regression parameters
  if(nregres_>0) {
    parse("REG_SCALE_MIN",scale_min_);
    parse("REG_SCALE_MAX",scale_max_);
    parse("REG_DSCALE",dscale_);
    // checks
    if(scale_max_<=scale_min_) error("with REGRESSION, REG_SCALE_MAX must be greater than REG_SCALE_MIN");
    if(dscale_<=0.) error("with REGRESSION, REG_DSCALE must be positive");
  }

  // scale factor
  parse("SCALE", scale_);

  // read map resolution
  double reso;
  parse("RESOLUTION", reso);
  if(reso<=0.) error("RESOLUTION should be strictly positive");

  // neighbor list stuff
  parse("NL_CUTOFF",nl_cutoff_);
  if(nl_cutoff_<=0.0) error("NL_CUTOFF should be explicitly specified and positive");
  parse("NL_STRIDE",nl_stride_);
  if(nl_stride_<=0) error("NL_STRIDE should be explicitly specified and positive");

  // averaging or not
  parseFlag("NO_AVER",no_aver_);

  // write overlap file
  parse("WRITE_OV_STRIDE", ovstride_);
  parse("WRITE_OV", ovfilename_);
  if(ovstride_>0 && ovfilename_=="") error("With WRITE_OV_STRIDE you must specify WRITE_OV");

  checkRead();

  // set parallel stuff
  size_=comm.Get_size();
  rank_=comm.Get_rank();

  // get number of replicas
  if(rank_==0) {
    if(no_aver_) {
      nrep_ = 1;
    } else {
      nrep_ = multi_sim_comm.Get_size();
    }
    replica_ = multi_sim_comm.Get_rank();
  } else {
    nrep_ = 0;
    replica_ = 0;
  }
  comm.Sum(&nrep_,1);
  comm.Sum(&replica_,1);

  if(nrep_>1 && dbfact_>0) error("Bfactor sampling not supported with ensemble averaging");

  log.printf("  atoms involved : ");
  for(unsigned i=0; i<atoms.size(); ++i) log.printf("%d ",atoms[i].serial());
  log.printf("\n");
  log.printf("  experimental data file : %s\n", datafile.c_str());
  if(no_aver_) log.printf("  without ensemble averaging\n");
  log.printf("  type of data noise : %s\n", noise.c_str());
  log.printf("  neighbor list cutoff : %lf\n", nl_cutoff_);
  log.printf("  neighbor list stride : %u\n",  nl_stride_);
  log.printf("  minimum uncertainty : %f\n",sigma_min);
  log.printf("  scale factor : %lf\n",scale_);
  log.printf("  reading/writing to status file : %s\n",statusfilename_.c_str());
  log.printf("  with stride : %u\n",statusstride_);
  if(nregres_>0) {
    log.printf("  regression stride : %u\n", nregres_);
    log.printf("  regression minimum scale : %lf\n", scale_min_);
    log.printf("  regression maximum scale : %lf\n", scale_max_);
    log.printf("  regression maximum scale MC move : %lf\n", dscale_);
  }
  if(noise_!=2) {
    log.printf("  initial value of the uncertainty : %f\n",sigma_ini);
    log.printf("  max MC move in uncertainty : %f\n",dsigma);
    log.printf("  MC stride : %u\n", MCstride_);
  }
  if(dbfact_>0) {
    log.printf("  max MC move in bfactor : %f\n",dbfact_);
    log.printf("  Bfactor MC stride : %u\n", MCBstride_);
  }
  if(errfile.size()>0) log.printf("  reading experimental errors from file : %s\n", errfile.c_str());
  log.printf("  temperature of the system in energy unit : %f\n",kbt_);
  log.printf("  prior exponent : %f\n",prior_);
  log.printf("  number of replicas for averaging: %u\n",nrep_);
  log.printf("  id of the replica : %u\n",replica_);
  if(nanneal_>0) {
    log.printf("  length of annealing cycle : %u\n",nanneal_);
    log.printf("  annealing factor : %f\n",kanneal_);
  }
  if(ovstride_>0) {
    log.printf("  stride for writing model overlaps : %u\n",ovstride_);
    log.printf("  file for writing model overlaps : %s\n", ovfilename_.c_str());
  }

  // calculate model GMM constant parameters
  vector<double> GMM_m_w = get_GMM_m(atoms);

  // read data file
  get_exp_data(datafile);
  log.printf("  number of kernel functions : %u\n", static_cast<unsigned>(GMM_d_m_.size()));

  // normalize atom weight map
  double norm_m = accumulate(GMM_m_w.begin(),  GMM_m_w.end(),  0.0);
  // renormalization
  for(unsigned i=0; i<GMM_m_w_.size(); ++i) GMM_m_w_[i] *= norm_d / norm_m;
  // constant cfact
  for(unsigned i=0; i<GMM_m_w_.size(); ++i) cfact_.push_back(GMM_m_w_[i]/pow( 2.0*pi, 1.5 ));

  // read experimental errors
  vector<double> exp_err;
  if(errfile.size()>0) exp_err = read_exp_errors(errfile);

  log.printf("  number of GMM groups : %u\n", static_cast<unsigned>(GMM_d_grps_.size()));
  // cycle on GMM groups
  for(unsigned Gid=0; Gid<GMM_d_grps_.size(); ++Gid) {
    log.printf("    group %d\n", Gid);
    // calculate median overlap and experimental error
    vector<double> ovdd;
    vector<double> err;
    // cycle on the group members
    for(unsigned i=0; i<GMM_d_grps_[Gid].size(); ++i) {
      // GMM id
      int GMMid = GMM_d_grps_[Gid][i];
      // add to experimental error
      if(errfile.size()>0) err.push_back(exp_err[GMMid]);
      else                 err.push_back(0.);
      // add to GMM overlap
      ovdd.push_back(ovdd_[GMMid]);
    }
    // calculate median quantities
    double ovdd_m = get_median(ovdd);
    double err_m  = get_median(err);
    // print out statistics
    log.printf("     # of members : %u\n", GMM_d_grps_[Gid].size());
    log.printf("     median overlap : %lf\n", ovdd_m);
    log.printf("     median error : %lf\n", err_m);
    // add minimum value of sigma for this group of GMMs
    sigma_min_.push_back(sqrt(err_m*err_m+sigma_min*ovdd_m*sigma_min*ovdd_m));
    // these are only needed with Gaussian and Outliers noise models
    if(noise_!=2) {
      // set dsigma
      dsigma_.push_back(dsigma * ovdd_m);
      // set sigma max
      sigma_max_.push_back(10.0*ovdd_m + sigma_min_[Gid] + dsigma_[Gid]);
      // initialize sigma
      sigma_.push_back(std::max(sigma_min_[Gid],std::min(sigma_ini*ovdd_m,sigma_max_[Gid])));
    }
  }

  // calculate auxiliary stuff
  calculate_useful_stuff(reso);

  // read status file if restarting
  if(getRestart()) read_status();

  // prepare data and derivative vectors
  ovmd_.resize(ovdd_.size());
  atom_der_.resize(GMM_m_type_.size());
  GMMid_der_.resize(ovdd_.size());

  // add components
  addComponentWithDerivatives("scoreb"); componentIsNotPeriodic("scoreb");
  if(dbfact_>0) {addComponent("accB"); componentIsNotPeriodic("accB");}
  if(noise_!=2) {addComponent("acc"); componentIsNotPeriodic("acc");}
  if(nregres_>0) {
    addComponent("scale");     componentIsNotPeriodic("scale");
    addComponent("accscale");  componentIsNotPeriodic("accscale");
    addComponent("enescale");  componentIsNotPeriodic("enescale");
  }
  if(nanneal_>0) {addComponent("anneal"); componentIsNotPeriodic("anneal");}

  // initialize random seed
  unsigned iseed;
  if(rank_==0) iseed = time(NULL)+replica_;
  else iseed = 0;
  comm.Sum(&iseed, 1);
  random_.setSeed(-iseed);

  // request the atoms
  requestAtoms(atoms);

  // print bibliography
  log<<"  Bibliography "<<plumed.cite("Bonomi, Camilloni, Bioinformatics, 33, 3999 (2017)");
  log<<plumed.cite("Hanot, Bonomi, Greenberg, Sali, Nilges, Vendruscolo, Pellarin, bioRxiv doi: 10.1101/113951 (2017)");
  log<<plumed.cite("Bonomi, Pellarin, Vendruscolo, bioRxiv doi: 10.1101/219972 (2017)");
  if(!no_aver_ && nrep_>1)log<<plumed.cite("Bonomi, Camilloni, Cavalli, Vendruscolo, Sci. Adv. 2, e150117 (2016)");
  log<<"\n";
}

void EMMIVOX::write_model_overlap(long int step)
{
  OFile ovfile;
  ovfile.link(*this);
  std::string num; Tools::convert(step,num);
  string name = ovfilename_+"-"+num;
  ovfile.open(name);
  ovfile.setHeavyFlush();
  ovfile.fmtField("%10.7e ");
// write overlaps
  for(int i=0; i<ovmd_.size(); ++i) {
    ovfile.printField("Model", ovmd_[i]);
    ovfile.printField("ModelScaled", scale_ * ovmd_[i]);
    ovfile.printField("Data", ovdd_[i]);
    ovfile.printField();
  }
  ovfile.close();
}

double EMMIVOX::get_median(vector<double> &v)
{
// dimension of vector
  unsigned size = v.size();
// in case of only one entry
  if (size==1) {
    return v[0];
  } else {
    // reorder vector
    sort(v.begin(), v.end());
    // odd or even?
    if (size%2==0) {
      return (v[size/2-1]+v[size/2])/2.0;
    } else {
      return v[size/2];
    }
  }
}

void EMMIVOX::read_status()
{
  double MDtime;
// open file
  IFile *ifile = new IFile();
  ifile->link(*this);
  if(ifile->FileExist(statusfilename_)) {
    ifile->open(statusfilename_);
    while(ifile->scanField("MD_time", MDtime)) {
      // read sigma only if not marginal noise
      if(noise_!=2) {
        for(unsigned i=0; i<sigma_.size(); ++i) {
          // convert i to string
          std::string num; Tools::convert(i,num);
          // read entries
          ifile->scanField("s"+num, sigma_[i]);
        }
      }
      // always read bfactors
      for(unsigned i=0; i<GMM_m_type_.size(); ++i) {
        // convert i to string
        std::string num; Tools::convert(i,num);
        // read entries
        ifile->scanField("bfact"+num, GMM_m_b_[GMM_m_res_[i]]);
      }
      // new line
      ifile->scanField();
    }
    ifile->close();
  } else {
    error("Cannot find status file "+statusfilename_+"\n");
  }
  delete ifile;
}

void EMMIVOX::print_status(long int step)
{
// if first time open the file
  if(first_status_) {
    first_status_ = false;
    statusfile_.link(*this);
    statusfile_.open(statusfilename_);
    statusfile_.setHeavyFlush();
    statusfile_.fmtField("%6.3e ");
  }
// write fields
  double MDtime = static_cast<double>(step)*getTimeStep();
  statusfile_.printField("MD_time", MDtime);
  // write sigma only if not marginal noise
  if(noise_!=2) {
    for(unsigned i=0; i<sigma_.size(); ++i) {
      // convert i to string
      std::string num; Tools::convert(i,num);
      // print entry
      statusfile_.printField("s"+num, sigma_[i]);
    }
  }
  // always write bfactors
  for(unsigned i=0; i<GMM_m_type_.size(); ++i) {
    // convert i to string
    std::string num; Tools::convert(i,num);
    // print entry
    statusfile_.printField("bfact"+num, GMM_m_b_[GMM_m_res_[i]]);
  }
  statusfile_.printField();
}

bool EMMIVOX::doAccept(double oldE, double newE, double kbt) {
  bool accept = false;
  // calculate delta energy
  double delta = ( newE - oldE ) / kbt;
  // if delta is negative always accept move
  if( delta < 0.0 ) {
    accept = true;
  } else {
    // otherwise extract random number
    double s = random_.RandU01();
    if( s < exp(-delta) ) { accept = true; }
  }
  return accept;
}

void EMMIVOX::doMonteCarlo()
{
  // extract random GMM group
  unsigned nGMM = static_cast<unsigned>(floor(random_.RandU01()*static_cast<double>(GMM_d_grps_.size())));
  if(nGMM==GMM_d_grps_.size()) nGMM -= 1;

  // generate random move
  double shift = dsigma_[nGMM] * ( 2.0 * random_.RandU01() - 1.0 );
  // new sigma
  double new_s = sigma_[nGMM] + shift;
  // check boundaries
  if(new_s > sigma_max_[nGMM]) {new_s = 2.0 * sigma_max_[nGMM] - new_s;}
  if(new_s < sigma_min_[nGMM]) {new_s = 2.0 * sigma_min_[nGMM] - new_s;}
  // old s2
  double old_inv_s2 = 1.0 / sigma_[nGMM] / sigma_[nGMM];
  // new s2
  double new_inv_s2 = 1.0 / new_s / new_s;

  // cycle on GMM group and calculate old and new energy
  double old_ene = 0.0;
  double new_ene = 0.0;
  double ng = static_cast<double>(GMM_d_grps_[nGMM].size());

  // in case of Gaussian noise
  if(noise_==0) {
    double chi2 = 0.0;
    for(unsigned i=0; i<GMM_d_grps_[nGMM].size(); ++i) {
      // id GMM component
      int GMMid = GMM_d_grps_[nGMM][i];
      // deviation
      double dev = ( scale_*ovmd_[GMMid]-ovdd_[GMMid] );
      // add to chi2
      chi2 += dev * dev;
    }
    // final energy calculation: add normalization and prior
    old_ene = 0.5 * kbt_ * ( chi2 * old_inv_s2 - (ng+prior_) * std::log(old_inv_s2) );
    new_ene = 0.5 * kbt_ * ( chi2 * new_inv_s2 - (ng+prior_) * std::log(new_inv_s2) );
  }

  // in case of Outliers noise
  if(noise_==1) {
    for(unsigned i=0; i<GMM_d_grps_[nGMM].size(); ++i) {
      // id GMM component
      int GMMid = GMM_d_grps_[nGMM][i];
      // calculate deviation
      double dev = ( scale_*ovmd_[GMMid]-ovdd_[GMMid] );
      // add to energies
      old_ene += std::log( 1.0 + 0.5 * dev * dev * old_inv_s2);
      new_ene += std::log( 1.0 + 0.5 * dev * dev * new_inv_s2);
    }
    // final energy calculation: add normalization and prior
    old_ene = kbt_ * ( old_ene + (ng+prior_) * std::log(sigma_[nGMM]) );
    new_ene = kbt_ * ( new_ene + (ng+prior_) * std::log(new_s) );
  }

  // increment number of trials
  MCtrials_ += 1.0;

  // accept or reject
  bool accept = doAccept(old_ene/anneal_, new_ene/anneal_, kbt_);
  if(accept) {
    sigma_[nGMM] = new_s;
    MCaccept_ += 1.0;
  }
  // local communication
  if(rank_!=0) {
    for(unsigned i=0; i<sigma_.size(); ++i) sigma_[i] = 0.0;
    MCaccept_ = 0.0;
  }
  if(size_>1) {
    comm.Sum(&sigma_[0], sigma_.size());
    comm.Sum(&MCaccept_, 1);
  }
}

void EMMIVOX::doMonteCarloBfact()
{

// iterator over bfactor map
  map<unsigned,double>::iterator it;

// cycle over bfactor map
  for(it=GMM_m_b_.begin(); it!=GMM_m_b_.end(); ++it) {

    // residue id
    unsigned ires = it->first;
    // old bfactor
    double bfactold = it->second;

    // propose move in bfactor
    double bfactnew = bfactold + dbfact_ * ( 2.0 * random_.RandU01() - 1.0 );
    // check boundaries
    if(bfactnew > bfactmax_) {bfactnew = 2.0*bfactmax_ - bfactnew;}
    if(bfactnew < bfactmin_) {bfactnew = 2.0*bfactmin_ - bfactnew;}

    // useful quantities
    map<unsigned, double> deltaov;
    Vector pos, posn, der;
    double dist;
    set<unsigned> ngbs;

    // cycle over all the atoms belonging to residue ires
    for(unsigned ia=0; ia<GMM_m_resmap_[ires].size(); ++ia) {

      // get atom id
      unsigned im = GMM_m_resmap_[ires][ia];
      // get atom type, bs, weight and position
      unsigned atype = GMM_m_type_[im];
      double bold = GMM_m_s_[atype]+bfactold/4.0;
      double bnew = GMM_m_s_[atype]+bfactnew/4.0;
      double pref = cfact_[atype];
      pos = getPosition(im);

      // cycle on all the components affected
      for(unsigned i=0; i<GMM_m_nb_[im].size(); ++i) {
        // voxel id
        unsigned id = GMM_m_nb_[im][i];
        // get contribution before change
        double dold=get_overlap(GMM_d_m_[id], pos, GMM_d_s_[id], pref, bold, der);
        // get contribution after change
        double dnew=get_overlap(GMM_d_m_[id], pos, GMM_d_s_[id], pref, bnew, der);
        // update delta overlap
        deltaov[id] += dnew-dold;
        // look for neighbors
        for(unsigned j=0; j<GMM_d_nb_[id].size(); ++j) {
          // atom index of potential neighbor
          unsigned in = GMM_d_nb_[id][j];
          // residue index of potential neighbor
          unsigned iresn = GMM_m_res_[in];
          // check if same residue
          if(ires==iresn) continue;
          // distance
          posn = getPosition(in);
          if(pbc_) dist = pbcDistance(pos,posn).modulo();
          else     dist = delta(pos,posn).modulo();
          // if closer than 0.5 nm, add residue to lists
          if(dist>0 && dist<0.5) ngbs.insert(iresn);
        }
      }
    }

    // now calculate new and old score
    map<unsigned, double>::iterator itov;
    double old_ene = 0.0;
    double new_ene = 0.0;

    if(noise_==0) {
      for(itov=deltaov.begin(); itov!=deltaov.end(); ++itov) {
        // id of the component
        unsigned id = itov->first;
        // new value
        double ovmdnew = ovmd_[id]+itov->second;
        // scores
        double devold = ( scale_*ovmd_[id]-ovdd_[id] ) / sigma_[GMM_d_beta_[id]];
        double devnew = ( scale_*ovmdnew  -ovdd_[id] ) / sigma_[GMM_d_beta_[id]];
        old_ene += 0.5 * kbt_ * devold * devold;
        new_ene += 0.5 * kbt_ * devnew * devnew;
      }
    }
    if(noise_==1) {
      for(itov=deltaov.begin(); itov!=deltaov.end(); ++itov) {
        // id of the component
        unsigned id = itov->first;
        // new value
        double ovmdnew = ovmd_[id]+itov->second;
        // scores
        double devold = ( scale_*ovmd_[id]-ovdd_[id] ) / sigma_[GMM_d_beta_[id]];
        double devnew = ( scale_*ovmdnew  -ovdd_[id] ) / sigma_[GMM_d_beta_[id]];
        old_ene += kbt_ * std::log( 1.0 + 0.5 * devold * devold );
        new_ene += kbt_ * std::log( 1.0 + 0.5 * devnew * devnew );
      }
    }
    if(noise_==2) {
      for(itov=deltaov.begin(); itov!=deltaov.end(); ++itov) {
        // id of the component
        unsigned id = itov->first;
        // new value
        double ovmdnew = ovmd_[id]+itov->second;
        // scores
        double devold = scale_*ovmd_[id]-ovdd_[id];
        double devnew = scale_*ovmdnew  -ovdd_[id];
        old_ene += -kbt_ * std::log( 0.5 / devold * erf ( devold * inv_sqrt2_ / sigma_min_[GMM_d_beta_[id]] ));
        new_ene += -kbt_ * std::log( 0.5 / devnew * erf ( devnew * inv_sqrt2_ / sigma_min_[GMM_d_beta_[id]] ));
      }
    }

    // add restraint to keep Bfactor of close atoms close
    for(set<unsigned>::iterator is=ngbs.begin(); is!=ngbs.end(); ++is) {
      double gold = (bfactold-GMM_m_b_[*is])/sqrt(bfactold+GMM_m_b_[*is])/0.058;
      double gnew = (bfactnew-GMM_m_b_[*is])/sqrt(bfactnew+GMM_m_b_[*is])/0.058;
      old_ene += 0.5 * kbt_ * gold * gold;
      new_ene += 0.5 * kbt_ * gnew * gnew;
    }

    // increment number of trials
    MCBtrials_ += 1.0;

    // accept or reject
    bool accept = doAccept(old_ene, new_ene, kbt_);

    // in case of acceptance
    if(accept==1) {
      // update acceptance rate
      MCBaccept_ += 1.0;
      // update bfactor
      it->second = bfactnew;
      // change all the ovmd_ affected
      for(itov=deltaov.begin(); itov!=deltaov.end(); ++itov) ovmd_[itov->first] += itov->second;
    }

  } // end cycle on bfactors

// now communicate results
  vector<unsigned> ires;
  vector<double>   newb;
// make it right dimension
  ires.resize(GMM_m_b_.size()); newb.resize(GMM_m_b_.size());
// master rank fills them in
  if(rank_==0) {
    unsigned i=0;
    for(it=GMM_m_b_.begin(); it!=GMM_m_b_.end(); ++it) {
      ires[i] = it->first;
      newb[i] = it->second;
      ++i;
    }
  } else {
    MCBaccept_ = 0.0;
  }
// communicate
  if(size_>1) {
    comm.Sum(&MCBaccept_, 1);
    comm.Sum(&ires[0], ires.size());
    comm.Sum(&newb[0], newb.size());
  }
// put things back in map
  for(unsigned i=0; i<ires.size(); ++i) GMM_m_b_[ires[i]] = newb[i];

}

vector<double> EMMIVOX::read_exp_errors(string errfile)
{
  int nexp, idcomp;
  double err;
  vector<double> exp_err;
// open file
  IFile *ifile = new IFile();
  if(ifile->FileExist(errfile)) {
    ifile->open(errfile);
    // scan for number of experimental errors
    ifile->scanField("Nexp", nexp);
    // cycle on GMM components
    while(ifile->scanField("Id",idcomp)) {
      // total experimental error
      double err_tot = 0.0;
      // cycle on number of experimental overlaps
      for(unsigned i=0; i<nexp; ++i) {
        string ss; Tools::convert(i,ss);
        ifile->scanField("Err"+ss, err);
        // add to total error
        err_tot += err*err;
      }
      // new line
      ifile->scanField();
      // calculate RMSE
      err_tot = sqrt(err_tot/static_cast<double>(nexp));
      // add to global
      exp_err.push_back(err_tot);
    }
    ifile->close();
  } else {
    error("Cannot find ERR_FILE "+errfile+"\n");
  }
  return exp_err;
}

vector<double> EMMIVOX::get_GMM_m(vector<AtomNumber> &atoms)
{
  // list of weights - one per atom
  vector<double> GMM_m_w;

  vector<SetupMolInfo*> moldat=plumed.getActionSet().select<SetupMolInfo*>();
  // map of atom types to A and B coefficients of scattering factor
  // f(s) = A * exp(-B*s**2)
  // B is in Angstrom squared
  // map between an atom type and an index
  map<string, unsigned> type_map;
  type_map["C"]=0;
  type_map["O"]=1;
  type_map["N"]=2;
  type_map["S"]=3;
  // fill in sigma vector
  GMM_m_s_.push_back(0.01*15.146);  // type 0
  GMM_m_s_.push_back(0.01*8.59722); // type 1
  GMM_m_s_.push_back(0.01*11.1116); // type 2
  GMM_m_s_.push_back(0.01*15.8952); // type 3
  // fill in weight vector
  GMM_m_w_.push_back(2.49982); // type 0
  GMM_m_w_.push_back(1.97692); // type 1
  GMM_m_w_.push_back(2.20402); // type 2
  GMM_m_w_.push_back(5.14099); // type 3

  // check if MOLINFO line is present
  if( moldat.size()==1 ) {
    log<<"  MOLINFO DATA found, using proper atom names\n";
    for(unsigned i=0; i<atoms.size(); ++i) {
      // get atom name
      string name = moldat[0]->getAtomName(atoms[i]);
      char type;
      // get atom type
      char first = name.at(0);
      // GOLDEN RULE: type is first letter, if not a number
      if (!isdigit(first)) {
        type = first;
        // otherwise is the second
      } else {
        type = name.at(1);
      }
      // check if key in map
      std::string type_s = std::string(1,type);
      if(type_map.find(type_s) != type_map.end()) {
        // save atom type
        GMM_m_type_.push_back(type_map[type_s]);
        // this will be normalized in the final density
        GMM_m_w.push_back(GMM_m_w_[type_map[type_s]]);
        // get residue id
        unsigned ires = moldat[0]->getResidueNumber(atoms[i]);
        // add to map and list
        GMM_m_resmap_[ires].push_back(i);
        GMM_m_res_.push_back(ires);
        // initialize Bfactor map
        GMM_m_b_[ires] = 0.0;
      } else {
        error("Wrong atom type "+type_s+" from atom name "+name+"\n");
      }
    }
  } else {
    error("MOLINFO DATA not found\n");
  }
  return GMM_m_w;
}

// read experimental data file in PLUMED format:
void EMMIVOX::get_exp_data(string datafile)
{
  Vector pos, sigma;
  double ov, beta;
  int idcomp;

// open file
  IFile *ifile = new IFile();
  if(ifile->FileExist(datafile)) {
    ifile->open(datafile);
    while(ifile->scanField("Id",idcomp)) {
      ifile->scanField("Mean_0",pos[0]);
      ifile->scanField("Mean_1",pos[1]);
      ifile->scanField("Mean_2",pos[2]);
      ifile->scanField("Sigma_0",sigma[0]);
      ifile->scanField("Sigma_1",sigma[1]);
      ifile->scanField("Sigma_2",sigma[2]);
      ifile->scanField("Beta",beta);
      ifile->scanField("Overlap",ov);
      // check beta
      if(beta<0) error("Beta must be positive!");
      // center of the Gaussian
      GMM_d_m_.push_back(pos);
      // vector of sigma
      GMM_d_s_.push_back(sigma);
      // beta
      GMM_d_beta_.push_back(beta);
      // experimental overlap
      ovdd_.push_back(ov);
      // new line
      ifile->scanField();
    }
    ifile->close();
  } else {
    error("Cannot find DATA_FILE "+datafile+"\n");
  }
  delete ifile;
  // now create a set from beta (unique set of values)
  set<int> bu(GMM_d_beta_.begin(), GMM_d_beta_.end());
  // now prepare the group vector
  GMM_d_grps_.resize(bu.size());
  // and fill it in
  for(unsigned i=0; i<GMM_d_beta_.size(); ++i) {
    if(GMM_d_beta_[i]>=GMM_d_grps_.size()) error("Check Beta values");
    GMM_d_grps_[GMM_d_beta_[i]].push_back(i);
  }
}

void EMMIVOX::calculate_useful_stuff(double reso)
{
  // We use the following definition for resolution:
  // the Fourier transform of the density distribution in real space
  // f(s) falls to 1/e of its maximum value at wavenumber 1/resolution
  // i.e. from f(s) = A * exp(-B*s**2) -> Res = sqrt(B).
  // average value of B
  double Bave = 0.0;
  for(unsigned i=0; i<GMM_m_type_.size(); ++i) {
    Bave += GMM_m_s_[GMM_m_type_[i]];
  }
  Bave /= static_cast<double>(GMM_m_type_.size());
  // calculate blur factor
  bfactmin_ = 1.0e-5;
  if(reso*reso>Bave) bfactmin_ = 4.0 * ( reso*reso-Bave );
  // initialize B factor to reasonable value based on Gaussian width at half maximum height equal the resolution
  double bfactini = 4.0 * ( 2.0 * pow(0.425*pi*reso,2) - Bave );
  // check for min and max
  bfactini = min(bfactmax_, max(bfactmin_, bfactini));
  for(map<unsigned,double>::iterator it=GMM_m_b_.begin(); it!=GMM_m_b_.end(); ++it) {
    it->second = bfactini;
  }
  log.printf("  experimental map resolution : %3.2f\n", reso);
  log.printf("  minimum Bfactor value       : %3.2f\n", bfactmin_);
  log.printf("  maximum Bfactor value       : %3.2f\n", bfactmax_);
  log.printf("  initial Bfactor value       : %3.2f\n", bfactini);
  // tabulate exponential
  dexp_ = dpcutoff_ / static_cast<double> (nexp_-1);
  for(unsigned i=0; i<nexp_; ++i) {
    tab_exp_.push_back(exp(-static_cast<double>(i) * dexp_));
  }
}

// get overlap and derivatives
double EMMIVOX::get_overlap(const Vector &d_m, const Vector &m_m, const Vector &d_s,
                            double pref, double m_b, Vector &ov_der)
{
  Vector md, invs2;
  // calculate vector difference with/without pbc
  if(pbc_) md = pbcDistance(m_m, d_m);
  else     md = delta(m_m, d_m);
  // calculate invs2
  for(unsigned i=0; i<3; ++i) invs2[i] = 1.0/(d_s[i]+inv_pi2_*m_b);
  // calculate exponent
  double ov = md[0]*md[0]*invs2[0]+md[1]*md[1]*invs2[1]+md[2]*md[2]*invs2[2];
  // final calculation
  ov = pref * sqrt(invs2[0]*invs2[1]*invs2[2]) * exp(-0.5*ov);
  // derivatives
  ov_der = ov * Vector(md[0]*invs2[0],md[1]*invs2[1],md[2]*invs2[2]);
  return ov;
}

void EMMIVOX::update_neighbor_list()
{
  // dimension of atom vectors
  unsigned GMM_m_size = GMM_m_type_.size();
  // local neighbor list
  vector < unsigned > nl_l;
  // other useful stuff
  vector<double> ov_l;
  map<double, unsigned> ov_m;
  map<double, unsigned>::iterator it;
  Vector d_m, d_s, md, invs2;
  double ov_tot, expov, ov;
  double ov_cut, res;
  unsigned atype, itab;
  int tot_size, rank_size;

  // cycle on GMM components - in parallel
  for(unsigned id=rank_; id<ovdd_.size(); id+=size_) {
    // clear overlap lists and map
    ov_l.clear(); ov_m.clear();
    // Kernel functions
    d_m = GMM_d_m_[id];
    d_s = GMM_d_s_[id];
    // total overlap with id
    ov_tot = 0.0;
    // cycle on all atoms
    for(unsigned im=0; im<GMM_m_size; ++im) {
      // calculate vector difference m_m-d_m with/without pbc
      if(pbc_) md = pbcDistance(getPosition(im), d_m);
      else     md = delta(getPosition(im), d_m);
      // get atom type
      atype = GMM_m_type_[im];
      // total value of b
      double b = GMM_m_s_[atype]+GMM_m_b_[GMM_m_res_[im]]/4.0;
      // calculate invs2
      for(unsigned i=0; i<3; ++i) invs2[i] = 1.0/(d_s[i]+inv_pi2_*b);
      // calculate exponent
      expov = md[0]*md[0]*invs2[0]+md[1]*md[1]*invs2[1]+md[2]*md[2]*invs2[2];
      // get index of expov in tabulated exponential
      itab = static_cast<unsigned> (round( 0.5*expov/dexp_ ));
      // check boundaries and skip atom in case
      if(itab >= tab_exp_.size()) continue;
      // in case calculate overlap
      ov = cfact_[atype]*sqrt(invs2[0]*invs2[1]*invs2[2])*tab_exp_[itab];
      // add to list
      ov_l.push_back(ov);
      // and map to retrieve atom index
      ov_m[ov] = im;
      // increase ov_tot
      ov_tot += ov;
    }
    // check if zero size -> ov_tot = 0
    if(ov_l.size()==0) continue;
    // define cutoff
    ov_cut = ov_tot * nl_cutoff_;
    // sort ov_l in ascending order
    std::sort(ov_l.begin(), ov_l.end());
    // integrate ov_l
    res = 0.0;
    for(unsigned i=0; i<ov_l.size(); ++i) {
      res += ov_l[i];
      // if exceeding the cutoff for overlap, stop
      if(res >= ov_cut) break;
      else ov_m.erase(ov_l[i]);
    }
    // now add atoms to neighborlist
    for(it=ov_m.begin(); it!=ov_m.end(); ++it)
      nl_l.push_back(id*GMM_m_size+it->second);
    // end cycle on GMM components in parallel
  }
  // find total dimension of neighborlist
  vector <int> recvcounts(size_, 0);
  recvcounts[rank_] = nl_l.size();
  comm.Sum(&recvcounts[0], size_);
  tot_size = accumulate(recvcounts.begin(), recvcounts.end(), 0);
  // clear and resize neighbor stuff
  nl_.clear(); nl_.resize(tot_size);
  // calculate vector of displacement
  vector<int> disp(size_);
  disp[0] = 0;
  rank_size = 0;
  for(unsigned i=0; i<size_-1; ++i) {
    rank_size += recvcounts[i];
    disp[i+1] = rank_size;
  }
  // Allgather neighbor list
  comm.Allgatherv(&nl_l[0], recvcounts[rank_], &nl_[0], &recvcounts[0], &disp[0]);
  // now resize derivatives
  ovmd_der_.resize(tot_size);
  // now cycle over the neighbor list to creat a list of voxels per atom
  GMM_m_nb_.clear(); GMM_m_nb_.resize(GMM_m_size);
  GMM_d_nb_.clear(); GMM_d_nb_.resize(ovdd_.size());
  for(unsigned i=0; i<tot_size; ++i) {
    unsigned id = nl_[i] / GMM_m_size;
    unsigned im = nl_[i] % GMM_m_size;
    GMM_m_nb_[im].push_back(id);
    GMM_d_nb_[id].push_back(im);
  }
}

void EMMIVOX::prepare()
{
  if(getExchangeStep()) first_time_=true;
}

// overlap calculator
void EMMIVOX::calculate_overlap() {

  if(first_time_ || getExchangeStep() || getStep()%nl_stride_==0) {
    update_neighbor_list();
    first_time_=false;
  }

  // clean temporary vectors
  for(unsigned i=0; i<ovmd_.size(); ++i)     ovmd_[i] = 0.0;
  for(unsigned i=0; i<ovmd_der_.size(); ++i) ovmd_der_[i] = Vector(0,0,0);

  // we have to cycle over all model and data GMM components in the neighbor list
  unsigned id, im, atype;
  unsigned GMM_m_size = GMM_m_type_.size();
  for(unsigned i=rank_; i<nl_.size(); i=i+size_) {
    // get data (id) and atom (im) indexes
    id = nl_[i] / GMM_m_size;
    im = nl_[i] % GMM_m_size;
    // get atom type
    atype = GMM_m_type_[im];
    // total value of b
    double b = GMM_m_s_[atype]+GMM_m_b_[GMM_m_res_[im]]/4.0;
    // add overlap with im component of model GMM
    ovmd_[id] += get_overlap(GMM_d_m_[id], getPosition(im), GMM_d_s_[id],
                             cfact_[atype], b, ovmd_der_[i]);
  }
  // communicate stuff
  if(size_>1) {
    comm.Sum(&ovmd_[0], ovmd_.size());
    comm.Sum(&ovmd_der_[0][0], 3*ovmd_der_.size());
  }
}

double EMMIVOX::scaleEnergy(double s)
{
  double ene = 0.0;
  for(unsigned i=0; i<ovdd_.size(); ++i) {
    ene += std::log( abs ( s * ovmd_[i] - ovdd_[i] ) );
  }
  return ene;
}

double EMMIVOX::doRegression()
{
// standard MC parameters
  unsigned MCsteps = 100000;
  double kbtmin = 1.0;
  double kbtmax = 10.0;
  unsigned ncold = 5000;
  unsigned nhot = 2000;
  double MCacc = 0.0;
  double kbt, ebest, scale_best;

// initial value of scale factor and energy
  double scale = random_.RandU01() * ( scale_max_ - scale_min_ ) + scale_min_;
  double ene = scaleEnergy(scale);
// set best energy
  ebest = ene;

// MC loop
  for(unsigned istep=0; istep<MCsteps; ++istep) {
    // get temperature
    if(istep%(ncold+nhot)<ncold) kbt = kbtmin;
    else kbt = kbtmax;
    // propose move in scale
    double ds = dscale_ * ( 2.0 * random_.RandU01() - 1.0 );
    double new_scale = scale + ds;
    // check boundaries
    if(new_scale > scale_max_) {new_scale = 2.0 * scale_max_ - new_scale;}
    if(new_scale < scale_min_) {new_scale = 2.0 * scale_min_ - new_scale;}
    // new energy
    double new_ene = scaleEnergy(new_scale);
    // accept or reject
    bool accept = doAccept(ene, new_ene, kbt);
    // in case of acceptance
    if(accept) {
      scale = new_scale;
      ene = new_ene;
      MCacc += 1.0;
    }
    // save best
    if(ene<ebest) {
      ebest = ene;
      scale_best = scale;
    }
  }
// calculate acceptance
  double accscale = MCacc / static_cast<double>(MCsteps);
// global communication
  if(!no_aver_ && nrep_>1) {
    if(replica_!=0) {
      scale_best = 0.0;
      ebest = 0.0;
      accscale = 0.0;
    }
    if(rank_==0) {
      multi_sim_comm.Sum(&scale_best, 1);
      multi_sim_comm.Sum(&ebest, 1);
      multi_sim_comm.Sum(&accscale, 1);
    }
  }
  // local communication
  if(rank_!=0) {
    scale_best = 0.0;
    ebest = 0.0;
    accscale = 0.0;
  }
  if(size_>1) {
    comm.Sum(&scale_best, 1);
    comm.Sum(&ebest, 1);
    comm.Sum(&accscale, 1);
  }
// set scale parameters
  getPntrToComponent("accscale")->set(accscale);
  getPntrToComponent("enescale")->set(ebest);
// return scale value
  return scale_best;
}

double EMMIVOX::get_annealing(long int step)
{
// default no annealing
  double fact = 1.0;
// position in annealing cycle
  unsigned nc = step%(4*nanneal_);
// useful doubles
  double ncd = static_cast<double>(nc);
  double nn  = static_cast<double>(nanneal_);
// set fact
  if(nc>=nanneal_   && nc<2*nanneal_) fact = (kanneal_-1.0) / nn * ( ncd - nn ) + 1.0;
  if(nc>=2*nanneal_ && nc<3*nanneal_) fact = kanneal_;
  if(nc>=3*nanneal_)                  fact = (1.0-kanneal_) / nn * ( ncd - 3.0*nn) + kanneal_;
  return fact;
}

void EMMIVOX::calculate()
{

// calculate CV
  calculate_overlap();

  // rescale factor for ensemble average
  double escale = 1.0 / static_cast<double>(nrep_);

  // in case of ensemble averaging, calculate average overlap
  if(!no_aver_ && nrep_>1) {
    // if master node, calculate average across replicas
    if(rank_==0) {
      multi_sim_comm.Sum(&ovmd_[0], ovmd_.size());
      for(unsigned i=0; i<ovmd_.size(); ++i) ovmd_[i] *= escale;
    } else {
      for(unsigned i=0; i<ovmd_.size(); ++i) ovmd_[i] = 0.0;
    }
    // local communication
    if(size_>1) comm.Sum(&ovmd_[0], ovmd_.size());
  }

  // get time step
  long int step = getStep();

  // do regression
  if(nregres_>0) {
    if(step%nregres_==0 && !getExchangeStep()) scale_ = doRegression();
    // set scale component
    getPntrToComponent("scale")->set(scale_);
  }

  // write model overlap to file
  if(ovstride_>0 && step%ovstride_==0) write_model_overlap(step);

  // clear energy and virial
  ene_ = 0.0;
  virial_.zero();

  // Gaussian noise
  if(noise_==0) calculate_Gauss();

  // Outliers noise
  if(noise_==1) calculate_Outliers();

  // Marginal noise
  if(noise_==2) calculate_Marginal();

  // get annealing rescale factor
  if(nanneal_>0) {
    anneal_ = get_annealing(step);
    getPntrToComponent("anneal")->set(anneal_);
  }

  // annealing rescale
  ene_ /= anneal_;

  // in case of ensemble averaging
  if(!no_aver_ && nrep_>1) {
    // if master node, sum der_GMMid derivatives and ene
    if(rank_==0) {
      multi_sim_comm.Sum(&GMMid_der_[0], GMMid_der_.size());
      multi_sim_comm.Sum(&ene_, 1);
    } else {
      // set der_GMMid derivatives and energy to zero
      for(unsigned i=0; i<GMMid_der_.size(); ++i) GMMid_der_[i]=0.0;
      ene_ = 0.0;
    }
    // local communication
    if(size_>1) {
      comm.Sum(&GMMid_der_[0], GMMid_der_.size());
      comm.Sum(&ene_, 1);
    }
  }

  // clean temporary vector
  for(unsigned i=0; i<atom_der_.size(); ++i) atom_der_[i] = Vector(0,0,0);

  // get derivatives of bias with respect to atoms
  unsigned id, im;
  Vector tot_der, pos;
  for(unsigned i=rank_; i<nl_.size(); i=i+size_) {
    // get indexes of data and model component
    id = nl_[i] / GMM_m_type_.size();
    im = nl_[i] % GMM_m_type_.size();
    // chain rule + replica normalization
    tot_der = GMMid_der_[id] * ovmd_der_[i] * escale * scale_ / anneal_;
    if(pbc_) pos = pbcDistance(GMM_d_m_[id], getPosition(im)) + GMM_d_m_[id];
    else     pos = getPosition(im);
    // increment derivatives and virial
    atom_der_[im] += tot_der;
    virial_ += Tensor(pos, -tot_der);
  }

  // communicate local derivatives and virial
  if(size_>1) {
    comm.Sum(&atom_der_[0][0], 3*atom_der_.size());
    comm.Sum(virial_);
  }

  // set derivatives, virial, and score
  for(unsigned i=0; i<atom_der_.size(); ++i) setAtomsDerivatives(getPntrToComponent("scoreb"), i, atom_der_[i]);
  setBoxDerivatives(getPntrToComponent("scoreb"), virial_);
  getPntrToComponent("scoreb")->set(ene_);

  // This part is needed only for Gaussian and Outliers noise models
  if(noise_!=2) {
    // do Monte Carlo
    if(dsigma_[0]>0 && step%MCstride_==0 && !getExchangeStep()) doMonteCarlo();
    // calculate acceptance ratio
    double acc = MCaccept_ / MCtrials_;
    // set value
    getPntrToComponent("acc")->set(acc);
  }

  // Monte Carlo on b factors
  if(dbfact_>0) {
    // do Monte Carlo
    if(step%MCBstride_==0 && !getExchangeStep()) doMonteCarloBfact();
    // calculate acceptance ratio
    double acc = MCBaccept_ / MCBtrials_;
    // set value
    getPntrToComponent("accB")->set(acc);
  }

  // print status
  if(step%statusstride_==0) print_status(step);
}

void EMMIVOX::calculate_Gauss()
{
  double eneg, dev;
  int GMMid;
  // cycle on all the GMM groups
  for(unsigned i=0; i<GMM_d_grps_.size(); ++i) {
    eneg = 0.0;
    // cycle on all the members of the group
    for(unsigned j=0; j<GMM_d_grps_[i].size(); ++j) {
      // id of the GMM component
      GMMid = GMM_d_grps_[i][j];
      // calculate deviation
      dev = ( scale_*ovmd_[GMMid]-ovdd_[GMMid] ) / sigma_[i];
      // add to group energy
      eneg += 0.5 * dev * dev;
      // store derivative for later
      GMMid_der_[GMMid] = kbt_ * dev / sigma_[i];
    }
    // add to total energy along with normalizations and prior
    ene_ += kbt_ * ( eneg + (static_cast<double>(GMM_d_grps_[i].size())+prior_) * std::log(sigma_[i]) );
  }
}

void EMMIVOX::calculate_Outliers()
{
  double eneg, dev;
  int GMMid;
  // cycle on all the GMM groups
  for(unsigned i=0; i<GMM_d_grps_.size(); ++i) {
    // cycle on all the members of the group
    eneg = 0.0;
    for(unsigned j=0; j<GMM_d_grps_[i].size(); ++j) {
      // id of the GMM component
      GMMid = GMM_d_grps_[i][j];
      // calculate deviation
      dev = ( scale_*ovmd_[GMMid]-ovdd_[GMMid] ) / sigma_[i];
      // add to group energy
      eneg += std::log( 1.0 + 0.5 * dev * dev );
      // store derivative for later
      GMMid_der_[GMMid] = kbt_ / ( 1.0 + 0.5 * dev * dev ) * dev / sigma_[i];
    }
    // add to total energy along with normalizations and prior
    ene_ += kbt_ * ( eneg + (static_cast<double>(GMM_d_grps_[i].size())+prior_) * std::log(sigma_[i]) );
  }
}

void EMMIVOX::calculate_Marginal()
{
  double dev, errf;
  int GMMid;
  // cycle on all the GMM groups
  for(unsigned i=0; i<GMM_d_grps_.size(); ++i) {
    // cycle on all the members of the group
    for(unsigned j=0; j<GMM_d_grps_[i].size(); ++j) {
      // id of the GMM component
      GMMid = GMM_d_grps_[i][j];
      // calculate deviation
      dev = ( scale_*ovmd_[GMMid]-ovdd_[GMMid] );
      // calculate errf
      errf = erf ( dev * inv_sqrt2_ / sigma_min_[i] );
      // add to group energy
      ene_ += -kbt_ * std::log ( 0.5 / dev * errf ) ;
      // store derivative for later
      GMMid_der_[GMMid] = - kbt_/errf*sqrt2_pi_*exp(-0.5*dev*dev/sigma_min_[i]/sigma_min_[i])/sigma_min_[i]+kbt_/dev;
    }
  }
}

}
}