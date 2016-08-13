/* main include file for read_SAR_tape */
/* R. Mellors August 1997 */
/* defines structures used in reading SAR tapes */
/*----------------------------*/
#include "sarleader_dss.h"
#include "sarleader_fdr.h"
#include "sarleader_platform.h"
#include "write_fdr_fixseg.h"
#include "write_fdr_varseg.h"
#include "write_dss.h"
#include "write_platform.h"
#include "sardata.h"

struct SAR_info {
	struct sarleader_dss 		*dss; 
	struct sarleader_dpaf_dss	*dpaf_dss; 
	struct sarleader_fdr_fixseg 	*fixseg;
	struct sarleader_fdr_varseg 	*varseg;
	struct platform 		*platform;
	struct position_vector 		*position;
	struct sardata_header 		*dataheader;
	struct sardata_rec 		*datarec;
	};
