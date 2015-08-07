/************************************************************************
* Creator: David T. Sandwell    Scripps Institution of Oceanography    *
* Date   : 09/12/93             Copyright, David T. Sandwell           *
************************************************************************/

void swap16( in, out, n )	/* Swaps 2 bytes within each 16-bit word of
				array in. */
char	*in;		/* Input array */
char	*out;		/* Output array */
int	n;		/* # of short integers to swap */

{
register char *ip, *op;	/* Local register variables */

	if( n > 0 )	/* Make sure n is positive */
	{
		ip = in+2;	/* Load the pointers into temporary registers */
		op = out;
		while( n-- )
		{
			*op++ = *--ip;	/* Do the swap */
			*op++ = *--ip;
			ip += 4;
		}
	}
}
