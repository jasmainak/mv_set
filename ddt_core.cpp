/************************************************
*
* ddt_core.cpp
*
* type "help ddt_core" for information
*
****************************************************/


//Includes
#include <iostream>
#include <fstream>

#include <string.h>
#include "mex.h"

#include <stdlib.h>
#include <math.h>
#include <stdio.h>

using namespace std;

# include "ddt_core.h"


//Constants
#define TR_DATA prhs[0]
#define TR_LABEL prhs[1]
#define TE_DATA prhs[2]
#define MA_CUTS prhs[3]
#define MODE_VAR prhs[4]

#define TR_PREDICT plhs[0]
#define TE_PREDICT plhs[1]
#define ERRVAL plhs[2]
#define PENVAL plhs[3]
#define OPT_SIZE plhs[4]
#define OPT_DEPTH plhs[5]

//For Mex-File, since we dont have a main
//ease compiler's worries here:
extern void _main();


/***************************************
*
* Global Routines ---------------------------------------------------------------------
*
***************************************/

/***************************************
* Global Routine :: FindInsertCavl
***************************************/
int FindInsertCavl(dict_node **curnod, dict_node **target_node, PATTERN *target_pattern)
{

  int SHB;
  int flag;
  dict_node *prov;


  if (*curnod == NULL)
    {
      *curnod = new dict_node(target_pattern);
      *target_node = *curnod;
      return(1);
    }


  SHB = target_pattern->diff_sign( (*curnod)->the_pattern );

  if ( SHB == 0)
    {
      *target_node = *curnod;
      return(0);
    }
 
  if (SHB > 0)
    {
      flag = FindInsertCavl(&(*curnod)->right, target_node, target_pattern);
      if (flag == 0)              /* height has changed? */
	return(0);
      else
	if ((*curnod)->imbal == 1)   /* must make a rotation */
	  if ((*curnod)->right->imbal == 1) /* then simple rotation */
	    {
	      (*curnod)->imbal = 0;
	      prov = (*curnod)->right;
	      (*curnod)->right = prov->left;
	      prov->left = *curnod;
	      *curnod = prov;
	      (*curnod)->imbal = 0;
	      return(0);
	    }
	  else                              /* then more complicated rotation */
	    {
	      prov=(*curnod)->right->left;
	      if (prov->imbal == 1)
		(*curnod)->imbal = -1;
	      else
		(*curnod)->imbal = 0;

	      if (prov->imbal == -1)
		(*curnod)->right->imbal = 1;
	      else
		(*curnod)->right->imbal = 0;

	      (*curnod)->right->left = prov->right;
	      prov->right = (*curnod)->right;
	      (*curnod)->right = prov->left;
	      prov->left = (*curnod);
	      (*curnod) = prov;
	      (*curnod)->imbal=0;
	      return(0);
	    }
	else                  /* no rotation here but height may have increased */
	  {
	    (*curnod)->imbal++;
	    return( (*curnod)->imbal);
	  }
    }
  else                       /* symmetric of the previous block */
    {
      flag = FindInsertCavl(&(*curnod)->left, target_node, target_pattern);
      if (flag == 0)              /* height has changed? */
	return(0);
      else
	if ((*curnod)->imbal == -1)   /* must make a rotation */
	  if ((*curnod)->left->imbal == -1) /* then simple rotation */
	    {
	      (*curnod)->imbal = 0;
	      prov = (*curnod)->left;
	      (*curnod)->left = prov->right;
	      prov->right = *curnod;
	      *curnod = prov;
	      (*curnod)->imbal = 0;
	      return(0);
	    }
	  else                              /* then more complicated rotation */
	    {
	      prov=(*curnod)->left->right;
	      if (prov->imbal == -1)
		(*curnod)->imbal = 1;
	      else
		(*curnod)->imbal = 0;

	      if (prov->imbal == 1)
		(*curnod)->left->imbal = -1;
	      else
		(*curnod)->left->imbal = 0;

	      (*curnod)->left->right = prov->left;
	      prov->left = (*curnod)->left;
	      (*curnod)->left = prov->right;
	      prov->right = (*curnod);
	      (*curnod) = prov;
	      (*curnod)->imbal=0;
	      return(0);
	    }
	else                  /* no rotation here but height may have increased */
	  {
	    (*curnod)->imbal--;
	    return( -(*curnod)->imbal);
	  }
    }
 
  cerr << "Error, this part of the code must not be attained!\n";
  exit(1);

}


/***************************************
* Global Routine :: FindCavl
***************************************/
dict_node *FindCavl( dict_node *dico, PATTERN *target_pattern) {

  int SHB;

  if (dico == NULL)
    return(NULL);


  SHB = target_pattern->diff_sign( dico->the_pattern );

  if ( SHB == 0)
    return(dico);
 
  if (SHB > 0)
    return( FindCavl( dico->right, target_pattern) );
  else
    return( FindCavl( dico->left, target_pattern) );
}



/***************************************
* Global Routine :: cmpdouble
***************************************/
int cmpdouble(const void *d1,const void *d2) {

  if ( *((double *) d1) == *((double*) d2))
    return(0);
  else if ( *((double *) d1) > *((double *) d2))
    return(1);
  else
    return(-1);
}



/***************************************
* Global Routine :: initialize_dictionary
* Compute dictionary elements at deepest level
***************************************/
dict_node *initialize_dictionary(const DATAPOINT *data_list, 
				 const unsigned short *resolutions, double *cost_data, int *final_size) 
{

  PATTERN *thepatterns(NULL);
  dict_node *the_dico(NULL),*cur_node(NULL);
  tree_bin *cur_bin(NULL);

  int *empty_count;
  double *empty_cost;

  empty_count = new int[NB_CLASSES+1];
  empty_cost = new double[NB_CLASSES+1];

  for ( int i = 0; i <= NB_CLASSES; i++) 
  {
    empty_count[i]=0;
    empty_cost[i]=0;
  }

  thepatterns = new PATTERN[NB_TRAIN[0]];

  int size_dico = 0;
  
  double min_vol=1;    // volume of smallest cells
  for (int i=0; i<DIMENSION; i++)
      for (int j=0; j < resolutions[i]; j++)
          min_vol = min_vol/2.0;
  
  // printf("max_dep: %d\n",max_dep);

  for (int i = 0; i<NB_TRAIN[0]; i++)
    {
      thepatterns[i].pat_from_double(data_list[i].feature,resolutions);
      FindInsertCavl(&the_dico,&cur_node,thepatterns+i);

      cur_bin = cur_node->the_bin;

      if ( cur_bin == NULL ) {   // pattern not yet registered in dico
        cur_bin = new tree_bin(empty_count,empty_cost,min_vol);
        cur_bin->inc_npoint();
        cur_node->the_bin = cur_bin;
        size_dico++;
      }

      cur_bin->count[0]++;
      cur_bin->count[data_list[i].class_point]++;
      cur_bin->cost[0]+=cost_data[i];      // total cost of all points, not just those misclassified
      cur_bin->cost[data_list[i].class_point]+=cost_data[i];
    }

  *final_size = size_dico;

  dict_node **dico_list = the_dico->big_fat_dump(size_dico);
  
  for (int i = 0; i<size_dico; i++) {   // update contrast and pred. class
    dico_list[i]->the_bin->compute_class(); 
    dico_list[i]->the_bin->compute_ctr();
  }
  
  
  delete[] dico_list;
  delete[] empty_count;
  delete[] empty_cost;
  delete[] thepatterns;

  return(the_dico);
  
}


/***************************************
* Global Routine :: init_rad_table
***************************************/
void init_rad_table()
{
    int largest_num = NB_TRAIN[0];  
    double *normed_binomial_table=NULL;

    RAD_TABLE = new double[largest_num+1];
    for (int i=0; i<=largest_num; i++)
        RAD_TABLE[i]=0;

    // build temporary table of binomial probabilities
    int index, index_low;
    if(largest_num%2==0)
    {
     //largest_num ist gerade
    index=((largest_num+2)/2)*(largest_num+1);
    }
    else
    {  
    index=(largest_num+2)*((largest_num+1)/2);
    }
    normed_binomial_table=new double[index];
    normed_binomial_table[0]=1;

    // ex update-Formula: 1/(2^n) (n choose k) = \frac{n}{2k}* [1/(2^(n-1)) (n-1 choose k-1)]

    for(int upper=1;upper<=largest_num;upper++)
    {
        //upper =oberer Eintrag im Binomialkoeff
        //(upper ueber 0) =1/2 * (upper-1 ueber 0) setzen

        //Runden!!!

        if(upper%2==0)
        {
         //upper ist gerade
            index=(upper/2)*(upper+1);
            index_low=(upper/2)*(upper-1);
        }
        else
        {
            index=upper*((upper+1)/2);
            index_low=upper*((upper-1)/2);
        }
        //index=round( upper*(upper+1)/2 );
        //index_low=round( (upper-1)*upper/2 );

        normed_binomial_table[index]=normed_binomial_table[index_low]/2.0;
        RAD_TABLE[upper]+= (double) normed_binomial_table[index]*upper;

        for(int lower=1; lower <= upper; lower++)
        {
            normed_binomial_table[index+lower]=0.5*upper/lower*normed_binomial_table[index_low+lower-1];
            RAD_TABLE[upper]+= (double) normed_binomial_table[index+lower]*abs(upper - 2*lower);
        }
        //mexPrintf("%d %1.8f\n",upper,RAD_TABLE[upper]);
    }
    delete[] normed_binomial_table;
  
}



/***************************************
* Global Routine :: initialize_logtable
***************************************/
/*void initialize_logtable()
{
  LOGTABLE = new double[SIZE_LOGTABLE+1];

  for (int i=1; i<=SIZE_LOGTABLE ;i++)
    LOGTABLE[i] = log((long double)i);
  
  LOGTABLE[0]=0;
}*/



/***************************************
* Global Routine :: calc_optimal_tree
***************************************/
tree_bin *calc_optimal_tree(DATAPOINT *train_data, 
                unsigned short *max_dims, double *cost_data)
{
  
  dict_node *dico_new, *dico_old, *dico_node_cur, *dico_sib, *dico_par, **dico_list;
  tree_bin *parent_bin, *sibling_bin, *cur_bin;
  PATTERN *parent_pattern, *sibling_pattern;
  int old_dico_size = 0, new_dico_size = 0;

  dico_old = initialize_dictionary(train_data, max_dims, cost_data, &old_dico_size);

  int totaliter=0;
  for (int i = 0; i < DIMENSION; i++ )   
    totaliter += max_dims[i];

  double cur_vol = 2;
  for (int i = 0; i < totaliter; i++)
      cur_vol = cur_vol/2.0;
  
  int *zerocount;                   // initializer for empty bins when necessary
  double *zerocost;                   // initializer for empty bins when necessary

  zerocount = new int[NB_CLASSES+1];
  zerocost = new double[NB_CLASSES+1];
  for (int i = 0; i <= NB_CLASSES; i++ )
  {
    zerocount[i] = 0;
    zerocost[i] = 0;
  }

  // loop over depths starting with max depth
  for (int iter = 0; iter<totaliter; iter++)  {    // Main loop    

    if (MODE == 3 && PENALTY == 3)
        EMPTY_PEN = CST_PEN + sqrt(cur_vol / (double) NB_TRAIN[0]);
    
    if (EXPERMODE!=1)
	{
      //cout << "      resolution:" << iter << " dictionary size: " << old_dico_size << "\n";
	  mexPrintf("      depth: %2d dictionary size: %d\n",totaliter-iter-1,old_dico_size);
	}

    dico_new = NULL;
    new_dico_size = 0;

    dico_list = dico_old->big_fat_dump(old_dico_size);

    // Loop on current dictionary entries
    for ( int dico_index = 0; dico_index < old_dico_size ; dico_index++) {  
    
      dico_node_cur = dico_list[dico_index];
      PATTERN cur_pat;

      cur_pat = dico_node_cur->the_pattern;
      cur_bin = dico_node_cur->the_bin;
  
      for (int dimiter = 0; dimiter<DIMENSION; dimiter++) {  // Loop on dimension


	if ( dico_node_cur->the_pattern.resolution(dimiter)>0) {

	  parent_pattern = cur_pat.create_parent(dimiter);	  
	  sibling_pattern = cur_pat.create_sibling(dimiter);
      
      // Create new entry if necessary with null tree_bin
	  FindInsertCavl(&dico_new,&dico_par,parent_pattern);  

	  dico_sib = FindCavl(dico_old,sibling_pattern);       // Find Only!

	  delete parent_pattern;
	  delete sibling_pattern;

	  if (dico_par->the_bin == NULL) {  // parent bin does not yet exist in new dico

	    new_dico_size++;

	    if (dico_sib == NULL) {  // Parent exists but not sibling

	      parent_bin = new tree_bin(cur_bin->count,cur_bin->cost,cur_vol*2.0);    
                     // create parent bin with counts
	                // ( nleaves=1 and cut_dir=-1)

	      parent_bin->compute_ctr();          // this is important;

	      if (cur_bin->ctr + EMPTY_PEN < parent_bin->ctr) { 
                // sibling bin is empty so its contrast = 0 + EMPTY_PEN
                // in this case, better so split


            sibling_bin = new tree_bin(zerocount,zerocost,cur_vol);        
                // effectively create empty sibling bin in this case

            sibling_bin->ctr = EMPTY_PEN;

            parent_bin->cut_dir = dimiter;
            parent_bin->nleaves = cur_bin->nleaves + sibling_bin->nleaves;
            parent_bin->ctr = cur_bin->ctr+sibling_bin->ctr;


            if ( dico_node_cur->the_pattern.is_left_bin(dimiter) ) {
              parent_bin->left = cur_bin;
              parent_bin->right = sibling_bin;
            } else {
              parent_bin->right = cur_bin;
              parent_bin->left = sibling_bin;
            }

            sibling_bin->inc_npoint();          // update pointer counting
            cur_bin->inc_npoint();

              }


            }	

            else {      // Parent and sibling already exist
                // Here we have to compare merits with parent bin without cut

              sibling_bin = dico_sib->the_bin;

              int *addcount;       // we construct the parent bin (without split)
              double *addcost;

              addcount = new int[NB_CLASSES+1];
              addcost = new double[NB_CLASSES+1];

              for (int i=0;i<=NB_CLASSES;i++) {
                  addcount[i] = cur_bin->count[i] + sibling_bin->count[i];
                  addcost[i] = cur_bin->cost[i] + sibling_bin->cost[i];
              }

              parent_bin = new tree_bin(addcount,addcost,cur_vol);   
                    // by default, cut_dir = -1 and nleaves = 1

              delete[] addcount;                    // this is not very pretty, maybe I'll fix this later
              delete[] addcost;

              parent_bin->compute_ctr();          // compute new contrast

              if (cur_bin->ctr+sibling_bin->ctr < parent_bin->ctr) {  
                  // in this case finally it's better to split

                parent_bin->cut_dir = dimiter;
                parent_bin->nleaves = cur_bin->nleaves + sibling_bin->nleaves;
                parent_bin->ctr = cur_bin->ctr+sibling_bin->ctr;

                if ( dico_node_cur->the_pattern.is_left_bin(dimiter) ) {
                  parent_bin->left = cur_bin;
                  parent_bin->right = sibling_bin;
                }
                else {
                  parent_bin->right = cur_bin;
                  parent_bin->left = sibling_bin;
                }

                sibling_bin->inc_npoint();          // update pointer counting
                cur_bin->inc_npoint();
              }

            }

            dico_par->the_bin = parent_bin;     // store the new parent tree bin in dico tree
            parent_bin->inc_npoint();              // and update pointer counting


          }
          else {       // parent bin already exists in new dico

            parent_bin = dico_par->the_bin;

            if (dico_sib == NULL) {

              if ( cur_bin->count[0] != parent_bin->count[0] ) // For debugging
                cout << "Gruh2\n";                             // This should never happen!

              if (cur_bin->ctr + EMPTY_PEN < parent_bin->ctr) { 
                      // sibling bin is empty so its contrast = 0 + EMPTY_PEN
                      // in this case, better to split


                sibling_bin = new tree_bin(zerocount,zerocost,cur_vol);        
                        // effectively create empty sibling bin in this case
                sibling_bin->ctr = EMPTY_PEN;

                parent_bin->cut_dir = dimiter;
                parent_bin->nleaves = cur_bin->nleaves + sibling_bin->nleaves;
                parent_bin->ctr = cur_bin->ctr+sibling_bin->ctr;

                if ( parent_bin->left != NULL)     // bookkeeping for pointer counters
                  parent_bin->left->dec_npoint();

                if ( parent_bin->right != NULL)
                  parent_bin->right->dec_npoint();

                if ( dico_node_cur->the_pattern.is_left_bin(dimiter) ) {
                  parent_bin->left = cur_bin;
                  parent_bin->right = sibling_bin;
                }
                else {
                  parent_bin->right = cur_bin;
                  parent_bin->left = sibling_bin;
                }

                sibling_bin->inc_npoint();          // update pointer counting
                cur_bin->inc_npoint();

              }
            }
            else {              

              sibling_bin = dico_sib->the_bin;

              if ( cur_bin->count[0] + sibling_bin->count[0] != parent_bin->count[0] )                                            
                cout << "Gruh\n";          // For debugging-this shoud never happen!



              if (cur_bin->ctr+sibling_bin->ctr < parent_bin->ctr) {  
                    // in this case the new split is better
                parent_bin->cut_dir = dimiter;
                parent_bin->nleaves = cur_bin->nleaves + sibling_bin->nleaves;
                parent_bin->ctr = cur_bin->ctr+sibling_bin->ctr;

                if ( parent_bin->left != NULL)   // bookkeeping for pointer counters
                  parent_bin->left->dec_npoint();

                if ( parent_bin->right != NULL)
                  parent_bin->right->dec_npoint();


                if ( dico_node_cur->the_pattern.is_left_bin(dimiter) ) {
                  parent_bin->left = cur_bin;
                  parent_bin->right = sibling_bin;
                }
                else {
                  parent_bin->right = cur_bin;
                  parent_bin->left = sibling_bin;
                }

                sibling_bin->inc_npoint();          // update pointer counting
                cur_bin->inc_npoint();

              } 
            }
          }
        }  
      }
    }

    cur_vol = cur_vol*2.0;     // update cell volume for next level   

    delete dico_old; 
    delete[] dico_list;

    dico_old = dico_new;
    old_dico_size = new_dico_size;
  }

  tree_bin *ret_value = dico_old->the_bin;

  dico_old->the_bin->inc_npoint();  // to avoid that the optimal tree gets deleted!

  delete dico_old;

  delete[] zerocount;
  delete[] zerocost;

  return(ret_value);
  
}  












/***************************************
*
* Class Routines ---------------------------------------------------------------------
*
***************************************/


/***************************************
* DATAPOINT :: public :: Default Contructor
****************************************/
DATAPOINT::DATAPOINT() : class_point(-1) 
{
	feature =  new double[DIMENSION];
}

/***************************************
* DATAPOINT :: public :: Overloaded Contructor (with params)
****************************************/
DATAPOINT::DATAPOINT(const double *newfeat,const int newclass) : class_point(newclass) 
{

	if (newclass < 1) 
	{
      cerr << "Class must be >0\n";
      exit(1);
    }

    feature =  new double[DIMENSION];
    
    for (int i=0; i < DIMENSION; i++) 
      feature[i] = newfeat[i];
}

/***************************************
* DATAPOINT :: public :: Destructor
****************************************/
DATAPOINT::~DATAPOINT() 
{
    delete[] feature;
}


/***************************************
* PATTERN :: public :: Default Contructor
****************************************/
PATTERN::PATTERN() 
{
    the_pattern = new unsigned short[DIMENSION];
    pat_resol = new unsigned short[DIMENSION];
}
  
/***************************************
* PATTERN :: public :: Destructor
****************************************/
PATTERN::~PATTERN() {
    delete[] the_pattern;
    delete[] pat_resol;
}
  

  
/***************************************
* PATTERN :: public :: pat_from_double
****************************************/
void PATTERN::pat_from_double(const double *datasrc, const unsigned short *resolutions) 
{
	for (int i=0; i<DIMENSION; i++) 
	{
	  the_pattern[i] = (unsigned short) (datasrc[i]*( (double) (1 << resolutions[i]) ) );
      pat_resol[i] = resolutions[i];
    }
}

/***************************************
* PATTERN :: public :: diff_sign
****************************************/
int PATTERN::diff_sign(const PATTERN &patcmp) 
{

    for (int i=0; i<DIMENSION; i++) 
	{
      if (pat_resol[i] > patcmp.pat_resol[i])
	return( 1 );
      else if (pat_resol[i] < patcmp.pat_resol[i])
	return( -1);
      else if (the_pattern[i] > patcmp.the_pattern[i])
	return( 1 );
      else if (the_pattern[i] < patcmp.the_pattern[i])
	return( -1 );
    }    
    return( 0 );
}


/***************************************
* PATTERN :: public :: operator=
****************************************/
void PATTERN::operator=(const PATTERN &patsrc) 
{
    for (int i=0; i<DIMENSION; i++) 
	{
      pat_resol[i] = patsrc.pat_resol[i];
      the_pattern[i] = patsrc.the_pattern[i];
    }
}
  
/***************************************
* PATTERN :: public :: create_parent
****************************************/
PATTERN *PATTERN::create_parent(const int dimension) 
{
    
    if (pat_resol[dimension] == 0 ) 
	{
      cerr << "Asked for parent of resolution 0 cell\n";
      exit(1);
    }

    PATTERN *parent_pattern;

    parent_pattern = new PATTERN;
    *parent_pattern = *this;
    ( parent_pattern->pat_resol[dimension] )-- ;
    parent_pattern->the_pattern[dimension] >>= 1;

    return( parent_pattern );
}

/***************************************
* PATTERN :: public :: create_sibling
****************************************/
PATTERN *PATTERN::create_sibling(const int dimension)  
{

    PATTERN *sibling_pattern;

    sibling_pattern = new PATTERN;

    *sibling_pattern = *this;

    sibling_pattern->the_pattern[dimension] ^= 1;

    return( sibling_pattern );
}


/***************************************
* PATTERN :: public :: is_left_bin
****************************************/
bool PATTERN::is_left_bin(const int dimension) 
{
    return( (the_pattern[dimension] & 1) >0 );
}

/***************************************
* PATTERN :: public :: resolution
****************************************/
short int PATTERN::resolution(const int dimension) 
{
    return(pat_resol[dimension]);
}



/***************************************
* tree_bin :: public :: Constructor
****************************************/
tree_bin::tree_bin(const int *new_count, const double *new_cost, double cur_vol) :
    npoint(0),left(NULL),right(NULL),ctr(-.0001),cut_dir(-1),class_bin(-1),nleaves(1)
{
	count = new int[NB_CLASSES+1];
	cost = new double[NB_CLASSES+1];
	int totalcount = 0;
    double totalcost = 0;
	double max = 0;
	short int argmax = 1;

	for (short int i=1; i<=NB_CLASSES; i++) 
	{
		count[i] = new_count[i];
        cost[i] = new_cost[i];
		totalcount += new_count[i];
		totalcost += new_cost[i];

		if (new_cost[i] > max) 
		{
			max = new_cost[i];
			argmax = i;
		}
	}  

	count[0] = totalcount;
    cost[0] = totalcost;

    vol=cur_vol;
//    printf("vol = %1.4f\n",vol);

    if (MODE != 3)
    {
        if (totalcost >0 )         // leave class_bin at -1 if empty bin
            class_bin = argmax;
    }
    else    // min vol or level set
    {
        if (vol > cost[1])
            class_bin = 2;
        else
            class_bin = 1;
    }

    
}


/***************************************
* tree_bin :: public :: Destructor
****************************************/
tree_bin::~tree_bin() 
{    
    delete[] count;
    delete[] cost;
}
  
/***************************************
* tree_bin :: public :: inc_npoint
****************************************/
void tree_bin::inc_npoint() 
{
    npoint++;
}
  
  
/***************************************
* tree_bin :: public :: dec_npoint
****************************************/
void tree_bin::dec_npoint()
{
      npoint--;

      if (npoint<0) 
	  {
		cerr << "Npoint negative, this should not happen\n";
		exit(1);
      }

      if (npoint == 0 ) 
	  {
		if (left != NULL) 
			left->dec_npoint();
		if (right != NULL)
			right->dec_npoint();

		delete this;      // self - destruct if nobody points to this object
      }
}
	
/***************************************
* tree_bin :: public :: compute_ctr
****************************************/
void tree_bin::compute_ctr()
{
	if (nleaves>1)
		cerr << "Warning: should not call compute_ctr() on a composite bin\n";

	ctr = 0.0;

/*    if (CRITERION == 0 && MODE == 1)
    {    
		for (int i=1;i<=NB_CLASSES;i++) 
		{
			if (count[i]>SIZE_LOGTABLE) 
			{
				cerr << "Log Table not long enough\n";
				exit(1);
			}
			ctr -= count[i]*LOGTABLE[count[i]];
		}
	
		if (count[0]>SIZE_LOGTABLE) 
		{
			cerr << "Log Table not long enough\n";
			exit(1);
		}
	
		ctr += count[0]*LOGTABLE[count[0]];
      	ctr += CST_PEN;   // Here we have to assume it's a noncomposite bin

    } */
    
    // empirical error term
    if (MODE != 3)
    {
        ctr+= cost[0] - cost[class_bin];
    }
    else
    {
        if (class_bin ==1)
            ctr+=vol - cost[1];
    }

    // penalty term
/*
    double tmp1, tmp2;
    tmp1 = PEN_WT[0]*sqrt(count[1]/(double) (NB_TRAIN[1]*NB_TRAIN[1]));
    tmp2 = PEN_WT[0]*RAD_TABLE[count[1]]/NB_TRAIN[1];
    mexPrintf("%d %d %1.6f %1.6f \n",MODE,PENALTY,tmp1,tmp2);
*/
    
    ctr += CST_PEN;     // constant term
       
    if (PENALTY > 1)
    {

        switch (MODE)
        {
            case 1:
                switch (PENALTY)
                {
                    case 2:
                        ctr+= PEN_WT[0]*sqrt(count[0]/(double) (NB_TRAIN[0]*NB_TRAIN[0]));
                        break;
                    case 3:
                        ctr+= PEN_WT[0]*RAD_TABLE[count[0]]/NB_TRAIN[0];
                        break;
                }
                break;
            case 2:
                for (int i=1; i<=NB_CLASSES; i++)
                {
                    switch (PENALTY)
                    {
                        case 2:
                            ctr+= PEN_WT[i-1]*sqrt(count[i]/(double) (NB_TRAIN[i]*NB_TRAIN[i]));
                            break;
                        case 3:
                            ctr+= PEN_WT[i-1]*RAD_TABLE[count[i]]/NB_TRAIN[i];
                            break;
                    }
                }
                break;
            case 3:
                switch (PENALTY)
                {
                    case 2:     // rademacher + jensen upper bound
                        ctr+= PEN_WT[0]*sqrt(count[1]/(double) (NB_TRAIN[1]*NB_TRAIN[1]));
                        break;
                    case 3:     // exact rademacher
                        ctr+= PEN_WT[0]*RAD_TABLE[count[1]]/NB_TRAIN[1];
                        break;
                    case 4:     // rad upper bound with volume term 
                        ctr+= PEN_WT[0]*sqrt((((double) count[1]/ (double) NB_TRAIN[1]) + vol)/ (double)NB_TRAIN[1]);
                        break;
                }
        }
	}
}


/***************************************
* tree_bin :: public :: compute_class
****************************************/
void tree_bin::compute_class() 
{
    double max = cost[1];     // previously: = cost[1];
    short int argmax = 1;
    
    if (MODE != 3)
    {
        for (short int i=1; i<=NB_CLASSES; i++) 
        {
            if (cost[i] > max) 
            {
                max = cost[i];
                argmax = i;
            }
        }  
        class_bin = argmax;
    }
    else    // min vol or level set
    {
        if (vol > cost[1])
            class_bin = 2;
        else
            class_bin = 1;
    }
}

/***************************************
* tree_bin :: public :: classify_datapoint
****************************************/
short tree_bin::classify_datapoint(DATAPOINT *data, unsigned short *resol) 
{
    if (resol==NULL) 
	{
      resol = new unsigned short[DIMENSION];
      for (int i = 0; i < DIMENSION; i++) 
		resol[i] = 0;
    }


    if (cut_dir == -1) 
	{
      delete[] resol;
      return(class_bin);
    }
    else 
	{
      resol[cut_dir]++;
      short rep;
      PATTERN check_pat;
      check_pat.pat_from_double(data->feature, resol);

      if (check_pat.is_left_bin(cut_dir)) 
		rep =  left->classify_datapoint(data, resol);
      else
		rep = right->classify_datapoint(data, resol);

      if ( rep >=0 )          // Little hack: if empty bin, return class of parent bin
		return(rep);
      else 
		return(class_bin);
    }
}


/***************************************
* tree_bin :: public :: compute_depth
****************************************/
int tree_bin::compute_depth() 
{
    if (cut_dir == -1) 
      return(0);
    else 
	{
      int l = left->compute_depth()+1;
      int r = right->compute_depth()+1;
      return( (l>r) ? l :r);
    }
}

/***************************************
* tree_bin :: public :: compute_stats
*
* Compute empirical error(s); used for final tree
* Compute unweighted penalty; used to determine penalt(y/ies) of final tree
* When first called, errors and pens should be vectors of zeros of length 
* NB_CLASSES + 1
****************************************/
void tree_bin::compute_stats(double *errors, double *pens, double *totalvol) 
{
    double numer, denom;

    if (cut_dir == -1)  // terminal node
    {
        // update errors
        if (MODE !=3)
            errors[0]+=cost[0]-cost[class_bin];   
        for (int i=1; i<=NB_CLASSES; i++)
        {
            if (i != class_bin)               
            {
//            printf("i = %d\n",i);
  //          printf("errors[i] = %1.4f\n",errors[i]);
    //        printf("cost[i] = %1.4f\n",cost[i]);
              errors[i]+=cost[i];
            }
        }
                    
        // update pens
        pens[0]+=CST_PEN;

        if (PENALTY == 2 || PENALTY == 4) 
        {
            for (int i=0; i<=NB_CLASSES; i++)
            {                    
                numer = (double) count[i];
                if (PENALTY == 4)
                    numer += vol*NB_TRAIN[i];
                denom = (double) (NB_TRAIN[i]*NB_TRAIN[i]);
                pens[i] += sqrt(numer/denom);
            }
        } 
        else if (PENALTY == 3)
        {
            for (int i=0; i<=NB_CLASSES; i++)
            {                    
                pens[i] += RAD_TABLE[count[i]]/NB_TRAIN[i];
            }            
        }
        else if (PENALTY !=1 )
        {
            mexErrMsgTxt("PENALTY invalid in compute_stats");
            return;                                            
        }

            
       
        
        // udpate volume if cell is included
        if (MODE == 3)
        {
            if (class_bin == 1)    
               totalvol[0] += vol;
        }
    }
    else // non-terminal
    {
      left->compute_stats(errors,pens,totalvol);
      right->compute_stats(errors,pens,totalvol);
    }
}



/***************************************
* dict_node :: public :: Constructor
****************************************/
dict_node::dict_node(const PATTERN *newpattern) :
    left(NULL), right(NULL), the_bin(NULL), imbal(0) 
{
    the_pattern = *newpattern;
}
    
/***************************************
* dict_node :: public :: Destructor
****************************************/
dict_node::~dict_node() 
{
    if (left != NULL)
      delete left;
    if (right != NULL)
      delete right;
    the_bin->dec_npoint();    // update pointer counting
}

/***************************************
* dict_node :: public :: big_fat_dump
****************************************/
dict_node **dict_node::big_fat_dump(const int dico_size) // ugly dump of pointers to dico nodes in a list
{  
    int count_size=0;
    dict_node **fat_list;
    fat_list = new dict_node*[dico_size];
    private_dump(fat_list,&count_size);
    
    if (dico_size != count_size) 
	{
      cerr << "Uh-Oh\n";
      exit(1);
    }

    return( fat_list );
}

/***************************************
* dict_node :: private :: private_dump
****************************************/
void dict_node::private_dump(dict_node **fat_list, int *counter) {
    
    fat_list[*counter] = this;
    
    (*counter)++;

    if ( left !=NULL)
      left->private_dump(fat_list, counter);

    if ( right !=NULL)
      right->private_dump(fat_list, counter);
}















/***************************************
*
* MAIN --------------------------------------------------------------------------
*
***************************************/


/************************************************
*
* mexFunction - portal to MATLAB
*
* This function is in place of 'main()'
*
************************************************/
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{

  //Declarations
  int nb_test_data;
  DATAPOINT *train_data, *test_data;
  unsigned short *max_dims;
  double *cost_data;
  int cost_flag=0;


//---Jaska Code Begin

  //Mex Declarations
  char *SwitchName;
  int SwitchLength;
  double *tr_data_ptr, *te_data_ptr, *tr_label_ptr, *Dim_Point; *cost_data; // inputs
  double *tr_predict_ptr, *te_predict_ptr, *errors, *pens, *opt_sz, *opt_dp;  // outputs

  /*******
   * New Input Parsing - Retrieving
   ******/


  //Set Defaults
  PENALTY = 0;          // error produced if pen_wt set but not penalty
//  CRITERION = 1;        // not log loss
  EXPERMODE = 1;        // quiet
  EMPTY_PEN = 0;
  CST_PEN = 0;

  //Check Minimum Inputs (this is not the most thorough thing in the world)
  if(nrhs < 5)
  {
	mexErrMsgTxt("At least 5 input arguments required");
	return;
  }
  else if(nrhs%2 != 1)
  {
	mexErrMsgTxt("Must have odd number of inputs");
	return;
  }
  
  /*******
   * New Data Acquisition
   *
   * Set up dimension parameters, get data pointers, etc...
   *
   ******/
	
	//Get Data Pointers
	tr_data_ptr = mxGetPr(TR_DATA);
	te_data_ptr = mxGetPr(TE_DATA);
    tr_label_ptr = mxGetPr(TR_LABEL);
	

    //Get DIMENSION,NB_CLASSES
	DIMENSION = mxGetM(TR_DATA);
    int test_dim =  mxGetM(TE_DATA);
    if (DIMENSION != test_dim)
    {
        mexErrMsgTxt("Train and test dimensions disagree\n");
        return;
    }
    
    NB_CLASSES=1;
    int tmp_label;
    for (int i=0; i<mxGetN(TR_DATA); i++)
    { 
        tmp_label = (int) tr_label_ptr[i] + 1;
        if (tmp_label > NB_CLASSES)
            NB_CLASSES = tmp_label;
    }

    // count number of training examples per class
    NB_TRAIN = new int[NB_CLASSES+1];
	NB_TRAIN[0] = mxGetN(TR_DATA);
    for (int i=1; i<=NB_CLASSES; i++)
        NB_TRAIN[i]=0;
    for (int i=0; i<NB_TRAIN[0]; i++)
    { 
        tmp_label = (int) tr_label_ptr[i] + 1;
        NB_TRAIN[tmp_label]+=1;
        // printf("number = %d, label = %d, NB_TRAIN[%d]=%d\n",i,tmp_label,tmp_label,NB_TRAIN[tmp_label]);
    }
    
    // get total number of test examples
	nb_test_data = mxGetN(TE_DATA);
            
	//Set up LOGTABLE for data
//	SIZE_LOGTABLE = ( NB_TRAIN[0] > nb_test_data )  ? NB_TRAIN[0] : nb_test_data ;
//	initialize_logtable();

	//Initialize "DataPoints"
	train_data = new DATAPOINT[NB_TRAIN[0]];
	test_data = new DATAPOINT[nb_test_data];
	
	//Populate "DataPoints"
	int MATindex = 0;
	for (int i=0; i<NB_TRAIN[0]; i++) 
    {
		
		(train_data+i)->class_point = (int)tr_label_ptr[i];
		(train_data+i)->class_point++;
		if ( (train_data+i)->class_point < 1 || (train_data+i)->class_point > NB_CLASSES)
		{
            printf("i = %d\n",i);
            printf("class = %d\n",(train_data+i)->class_point);
			mexErrMsgTxt("Train: Class out of bounds\n");
			return;
		}
		
		for(int j=0;j<DIMENSION;j++)
		{
			(train_data+i)->feature[j] = tr_data_ptr[MATindex+j];
		}

		MATindex+=DIMENSION;
	}
	
	MATindex = 0;
	for (int i=0; i<nb_test_data; i++) 
    {
		
/*		(test_data+i)->class_point = (int)te_label_ptr[i];
		(test_data+i)->class_point++;
		if ( (test_data+i)->class_point < 1 || (test_data+i)->class_point > NB_CLASSES)
		{
            printf("i = %d\n",i);
            printf("class = %d\n",(test_data+i)->class_point);
			mexErrMsgTxt("Test: Class out of bounds\n");
			return;
		} */
		
		for(int j=0;j<DIMENSION;j++)
		{
			(test_data+i)->feature[j] = te_data_ptr[MATindex+j];
		}

		MATindex+=DIMENSION;
	}
	
	//Set up maximum dimensions; user can input scalar or vector
	max_dims = new unsigned short[DIMENSION];
	int total=0;
	
    if (mxGetM(MA_CUTS) > 1)
    {
		mexErrMsgTxt("max_cuts must be a row vector");
		return;
    }
        
    if (mxGetN(MA_CUTS)==1)
    {
        unsigned short tmp_dim = (unsigned short) mxGetScalar(MA_CUTS);
        for(int i=0; i<DIMENSION; i++)
        {
            max_dims[i] = tmp_dim;
            total += max_dims[i];
        }
    }
    else if (mxGetN(MA_CUTS)==DIMENSION)
    {
      	Dim_Point = mxGetPr(MA_CUTS);
        for(int i=0; i<DIMENSION; i++)
        {
            max_dims[i] = (unsigned short) Dim_Point[i];
            total += max_dims[i];
        }        
    }
    else
    {
		mexErrMsgTxt("max_cuts invalid");
		return;
    }
        
    
  // Get problem 
  MODE = (int) mxGetScalar(MODE_VAR);

  //Get Switches and Set Values
  for(int i=5; i<nrhs; i+=2)
  {
	if (mxIsChar(prhs[i]) != 1)
	{
		mexErrMsgTxt("Even inputs from the 6th onward must be strings");
		return;
	}
	else
	{
		//Get Switch Name
		SwitchLength = mxGetN(prhs[i])+1;
		SwitchName = (char *) mxCalloc(SwitchLength, sizeof(char));
		mxGetString(prhs[i],SwitchName,SwitchLength);
		
		//Compare Switches and Set
//		if(!strcmp(SwitchName,"criterion"))
//		{
//			CRITERION = (int) mxGetScalar(prhs[i+1]);
//		}
		if(!strcmp(SwitchName,"cst_pen"))
		{
			CST_PEN = (double) mxGetScalar(prhs[i+1]);
		}
		else if(!strcmp(SwitchName,"pen_wt"))
		{
            if (PENALTY==0)
            {
                 mexErrMsgTxt("Please define penalty before pen_wt");
                 return;                
            }
            else if (PENALTY==1)    // constant
            {
                double *tmp_ptr = mxGetPr(prhs[i+1]);
                if (mxGetN(prhs[i+1]) > 1)
                    cout << "Warning: only one weight needed for constant penalty\n";
                CST_PEN = tmp_ptr[0];
                EMPTY_PEN = CST_PEN;
            }
            else if (PENALTY==2 || PENALTY == 3 || PENALTY == 4)    // Rademacher and variants
            {
                int num_wt = mxGetN(prhs[i+1]);
                if (num_wt == NB_CLASSES)
                    PEN_WT = mxGetPr(prhs[i+1]);
                else if (num_wt == 1)
                {
                    double tmp_wt = mxGetScalar(prhs[i+1]);
                    PEN_WT = new double[NB_CLASSES];
                    for (int k=0; k<NB_CLASSES; k++)
                        PEN_WT[k]=tmp_wt;
                }
                else 
                {
                    mexErrMsgTxt("Number of penalty weights must be 1 or number of classes");
                    return;                
                }
            }
            else
            {
                mexErrMsgTxt("PENALTY variable invalid");
                return;                                
            }
		}
		else if(!strcmp(SwitchName,"expermode"))
		{
			EXPERMODE = (int)mxGetScalar(prhs[i+1]);
		}
		else if(!strcmp(SwitchName,"penalty"))
		{
			PENALTY = (int)mxGetScalar(prhs[i+1]);
            if (PENALTY == 4 && MODE != 3)
            {
                mexErrMsgTxt("PENALTY = 4 only valid when MODE == 3");
                return;                                
            }
		}
  		else if(!strcmp(SwitchName,"costs"))
		{
            int num_costs = mxGetN(prhs[i+1]);
            if (num_costs != NB_TRAIN[0])
            {
                mexErrMsgTxt("Number of costs must equal number of training examples");
                return;
            }
			cost_data = mxGetPr(prhs[i+1]);
            cost_flag=1;
		}

		else
		{
			mexErrMsgTxt("Unknown Parameter");
			return;
		}
	}
  }
    
  // define cost_data if not input by user
  if (cost_flag==0)
  {
    cost_data = new double[NB_TRAIN[0]];
    for (int i=0; i<NB_TRAIN[0]; i++)
        cost_data[i]=1.0/NB_TRAIN[0];
  }
    
  if (PENALTY == 3)
      init_rad_table();
  
  if (EXPERMODE == 2)
    {
        printf("\n");//debug
        printf("DIM: %i\n",DIMENSION);//debug
        printf("NB_C: %i\n",NB_CLASSES);//debug
        printf("EMPTY_PEN: %2.3f\n",EMPTY_PEN);//debug
        printf("nb_train: %d\n",NB_TRAIN[0]);//debug
        printf("nb_test: %d\n",nb_test_data);//debug*/
    }

  if (EXPERMODE == 2)
  {
      //Print the switch list
      printf("Paramter List:\n");
      printf("penalty: %d\n",(int)PENALTY);
//      printf("criterion: %d\n",(int)CRITERION);
      if (PENALTY==1 && MODE==1)
          printf("penalty weight: %2.3f\n",CST_PEN);
      if (PENALTY==2 && MODE==2)
      {
          for (int k=0; k<NB_CLASSES; k++)
          printf("penalty weight: %2.3f\n",PEN_WT[k]);
      }
      printf("exper_mode: %d\n",(int)EXPERMODE);
  }

//---Jaska Code End

  if (total > 25 && EXPERMODE ==2 )
  {
    //cerr << "WARNING: total nb of cut dims > 25 !\n";
	printf("WARNING: total nb of cut dims > 25 !\n");
  }

  // Call main subroutine
  tree_bin *Opt_tree = calc_optimal_tree(train_data, max_dims, cost_data);
  
  if (EXPERMODE == 2) {
    printf("1st cut direction: %d\n", Opt_tree->cut_dir);
    printf("Nb of leaves: %d\n" ,Opt_tree->nleaves);
    printf("Total nb of exmp: %d\n",Opt_tree->count[0]);
    printf("Contrast: %f\n", Opt_tree->ctr);
    printf("Depth: %d\n",Opt_tree->compute_depth() );
    for (int i=0; i<=NB_CLASSES; i++)
        printf("NB_TRAIN[%d]=%d\n",i,NB_TRAIN[i]);

  }

    
  // Set up outputs
  TR_PREDICT = mxCreateDoubleMatrix(1,NB_TRAIN[0],mxREAL);
  tr_predict_ptr = mxGetPr(TR_PREDICT);
  TE_PREDICT = mxCreateDoubleMatrix(1,nb_test_data,mxREAL);
  te_predict_ptr = mxGetPr(TE_PREDICT);
  ERRVAL = mxCreateDoubleMatrix(1,NB_CLASSES+1,mxREAL);
  errors = mxGetPr(ERRVAL);
  PENVAL = mxCreateDoubleMatrix(1,NB_CLASSES+1,mxREAL);
  pens = mxGetPr(PENVAL);
  OPT_SIZE = mxCreateDoubleMatrix(1,1,mxREAL);
  opt_sz = mxGetPr(OPT_SIZE);
  OPT_DEPTH = mxCreateDoubleMatrix(1,1,mxREAL);
  opt_dp = mxGetPr(OPT_DEPTH);
  
  double train_err = 0;

  for (int i=0; i < NB_TRAIN[0]; i++) {
    tr_predict_ptr[i]= Opt_tree->classify_datapoint(train_data+i,NULL)-1;
  }

  // Predict labels of test data
  for (int i=0; i < nb_test_data; i++) {
    te_predict_ptr[i]= Opt_tree->classify_datapoint(test_data+i,NULL)-1;
  }

  double *totalvol;
  totalvol = new double[1];
  totalvol[0]=0;
  for (int i=0; i<=NB_CLASSES; i++)
  {
      errors[i]=0;
      pens[i]=0;
  } 
  Opt_tree->compute_stats(errors,pens,totalvol);
  if (MODE==3)
      errors[0]=totalvol[0];

  *opt_sz=(double) Opt_tree->nleaves;
  *opt_dp=(double) Opt_tree->compute_depth();
  
//  printf("errors[1] = %2.5f\n",errors[1]);
//  printf("pens[1] = %2.9f\n",pens[1]);
  
  
  // printf("%f\n", (double) test_err/(double)nb_test_data);

  // final clean-up, really useful only when using Valgrind 
  // (to be sure that everything is OK on the memory side)

  delete[] max_dims;
  delete[] train_data;
  delete[] test_data;
  
  Opt_tree->dec_npoint();

//  delete[] LOGTABLE;

}

  
