/************************************************************************
* Creator: David T. Sandwell    Scripps Institution of Oceanography    *
* Date   : 09/12/93             Copyright, David T. Sandwell           *
************************************************************************/

void swap32( in, out, n )	/* Swaps 4 bytes within each 32-bit word of
				array in. */
char	*in;		/* Input array */
char	*out;		/* Output array */
int	n;		/* # of short integers to swap */

{
register char *ip, *op;	/* Local register variables */

	if( n > 0 )
	{
		ip = in+4;
		op = out;
		while( n-- )
		{
			*op++ = *--ip;
			*op++ = *--ip;
			*op++ = *--ip;
			*op++ = *--ip;
			ip += 8;
		}
	}
}
