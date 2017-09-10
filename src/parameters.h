/*
 * parameters.h
 *
 *  Created on: Feb 3, 2013
 *      Author: olga
 */

#ifndef PARAMETERS_H_
#define PARAMETERS_H_

#define kspring_cov 20.0 //kcal/molxA
#define R_limit  2.0  //A
#define Sigma	 3.8 // A
#define El		 1.0 //kcal/mol
#define covalent_cutoff 10.0 //A, source:soptop.c default value
#define pairs_cutoff	20.0
#define rLimitBond 8.0
#define Pi	3.14159265

#define BUF_SIZE	256

#define MAX_BONDS	6
#define MAX_NATIVE	128
#define MAX_PAIRS	2048
#define MAX_TENSORS 10
/*
#define DEFAULT_MODE			0
#define CRYSTAL_MODE			1
#define CAPSID_MODE				2
#define CRYSTAL_MODE_STRING		"crystal"
#define CAPSID_MODE_STRING		"capsid"

#define STRESS_CAUCHY_ON		"stressCauchy"
#define STRESS_PIOLAI_ON		"stressPiolaI"
#define STRESS_PIOLAII_ON		"stressPiolaII"
#define STRESS_KIRCHHOFF_ON		"stressKirchhoff"
#define STRESS_BIOT_ON			"stressBiot"

#define STRETCH_ON				"stretch"
#define STRAIN_GREEN_ON			"strainGreen"
#define STRAIN_ALMANSI_ON		"strainAlmansi"
#define STRAIN_BIOT_ON			"strainBiot"

#define OUTPUT_STRING			"output"
#define OUTPUT_PDB_STRING		"outputPDB"

#define PRINTING_ON				"printPDB"

#define MAX_FILES				9

#define INITIAL_FRAME			"initialFrame"
#define FINAL_FRAME				"finalFrame"

#define TENSOR_CALC_STRING		"tensorCalculation"
#define TENSOR_PRINT_STRING		"tensorPrint"

#define AVERAGE_MODE_CUTOFF		1
#define AVERAGE_MODE_CHAIN		2
#define AVERAGE_MODE_CUTOFF2	3
#define AVERAGE_MODE_DEFAULT	0

#define AVERAGE_MODE			"average"

#define VECTOR_STRING			"vector"
*/


#endif /* PARAMETERS_H_ */
