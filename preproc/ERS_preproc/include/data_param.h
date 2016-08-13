/* include file for CEOS data files parameters */
#define MIN_PRF	1640.0
#define MAX_PRF 1720.0
#define SEC_PER_PRI_COUNT 210.94e-09
#define LINESIZE 12060
struct lineparam {
int	ifc;
unsigned short	swst_dn;
unsigned short	pri_dn;
double 		pri;
double		swst;
};
