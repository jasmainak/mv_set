/************************************************
*
* ddt_core.h
*
* Reorganized and Re-Wrapped for MATLAB by:
* Jason Laska, Mark Davenport
*
*************************************************/

/************************************************
*
* --Header Contents--
*
* Classes:
*	DATAPOINT
*	PATTERN
*	tree_bin
*	dict_node
*
* Global Routines:
*	cmpdouble
*	FindInsertCavl
*	FindCavl
*	initialize_dictionary
*	Box_Data
*	Scale_Data_using_Quantiles
*	initialize_logtable
*	calc_optimal_tree
*
**************************************************/





//Global Variables
int DIMENSION=0;        // dimension of data
int NB_CLASSES=0;       // number of classes
double *RAD_TABLE=NULL;  // precomputed table of rademacher averages
// double *LOGTABLE=NULL;  // precomputed table of nat. logs of integers
int SIZE_LOGTABLE=-1;   
int EXPERMODE=0;        // indicator variable; how much to display 
// int CRITERION=0;        // loss function for classification problem
double CST_PEN=0;         // weight for constant penalty

// Added by Clay
int PENALTY=0;          // penalty indicator variable: 
    // 1->constant, 2->Rad upper bound, 3-> Exact Rademacher, 4-> Rad w/ volume term (MV/DL only) 
double EMPTY_PEN=0;     // penalty assigned to empty cell
double *PEN_WT=NULL;    // weights for Rademacher penalty
int *NB_TRAIN=NULL;     // [0]=total number, [i]=number in class i, 1<=i<=NB_CLASSES
int MODE=0;             // PE,CS,RL -> 1, NP,MM -> 2, and MV,DL -> 3

/***************************************
* Class
* 
* DATAPOINT
***************************************/
class DATAPOINT 
{
  public:

  //Public Variables
  double *feature;
  int class_point;

  //Public Functions
  DATAPOINT();
  DATAPOINT(const double *newfeat,const int newclass);
  ~DATAPOINT();
  
};



/***************************************
* Class
* 
* PATTERN
***************************************/
class PATTERN 
{
private:
  
  //Private Variables
  unsigned short *the_pattern;
  unsigned short int *pat_resol;

public:

  //Public Functions
  PATTERN();
  ~PATTERN();
  void pat_from_double(const double *datasrc, const unsigned short *resolutions);
  int diff_sign(const PATTERN &patcmp);
  void operator=(const PATTERN &patsrc);
  PATTERN *create_parent(const int dimension);
  PATTERN *create_sibling(const int dimension);
  bool is_left_bin(const int dimension);
  short int resolution(const int dimension);

};



/***************************************
* Class
* 
* tree_bin
***************************************/
class tree_bin 
{
private:
  
  //Private Variables
  int npoint; // pointer counting for memory management


public:
  
  //Public Variables
  tree_bin *left;
  tree_bin *right;
  int *count;
  double *cost; // total cost of all points, not just those misclassified
  double ctr;
  int cut_dir;
  short int class_bin;
  int nleaves;
  double vol;

  //Public Functions
  tree_bin(const int *new_count, const double *new_cost, double cur_vol);
  ~tree_bin();
  void inc_npoint();
  void dec_npoint();
  void compute_ctr();
  void compute_class();
  short classify_datapoint(DATAPOINT *data, unsigned short *resol);
  int compute_depth();
  void compute_stats(double *errors, double *pens, double *totalvol);
  
};




/***************************************
* Class
* 
* dict_node
***************************************/
class dict_node 
{
public:

  //Public Variables
  dict_node *left;
  dict_node *right;
  tree_bin *the_bin;
  PATTERN the_pattern;
  short int imbal;

  //Public Functions
  dict_node(const PATTERN *newpattern);
  ~dict_node();
  dict_node **big_fat_dump(const int dico_size); // ugly dump of pointers to dico nodes in a list

private:

  //Private Functions
  void private_dump(dict_node **fat_list, int *counter);
    
};



/***************************************
* Global Routine :: cmpdouble
***************************************/
int cmpdouble(const void *d1,const void *d2);

/***************************************
* Global Routine :: FindInsertCavl
***************************************/
int FindInsertCavl(const dict_node **pt_to_dico, const dict_node **target_node, const PATTERN *target_pattern);

/***************************************
* Global Routine :: FindCavl
***************************************/
dict_node *FindCavl(const dict_node *dico, const PATTERN *target_pattern);

/***************************************
* Global Routine :: initialize_dictionary
***************************************/
dict_node *initialize_dictionary(const DATAPOINT *data_list, 
				 const unsigned short *resolutions, double *cost_data, int *final_size);
				
/***************************************
* Global Routine :: initialize_logtable
***************************************/
// void initialize_logtable();

/***************************************
* Global Routine :: init_rad_table
***************************************/
void init_rad_table();

/***************************************
* Global Routine :: calc_optimal_tree
***************************************/
tree_bin *calc_optimal_tree(DATAPOINT *train_data,  
                unsigned short *max_dims, double *cost_data);
