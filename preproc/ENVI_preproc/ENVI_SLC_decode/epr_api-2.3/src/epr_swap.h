/*
 * $Id: epr_swap.h,v 1.1.1.1 2004-10-28 19:22:23 norman Exp $
 *
 * Copyright (C) 2002 by Brockmann Consult (info@brockmann-consult.de)
 *
 * This program is free software; you can redistribute it and/or modify it
 * under the terms of the GNU General Public License as published by the
 * Free Software Foundation. This program is distributed in the hope it will
 * be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
 * of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 * See the GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
 */

#ifndef EPR_SWAP_H_INCL
#define EPR_SWAP_H_INCL

#include <stdio.h> /* just to get the ANSI-C type FILE */

#ifdef __cplusplus
extern "C" 
{
#endif

void byte_swap_short(short *buffer, uint number_of_swaps);
void byte_swap_ushort(ushort* buffer, uint number_of_swaps);
void byte_swap_long(int *buffer, uint number_of_swaps);
void byte_swap_uint(uint* buffer, uint number_of_swaps);
void byte_swap_float(float* buffer, uint number_of_swaps);
void epr_swap_endian_order(const EPR_SField* field);
int epr_is_big_endian_order();
int epr_is_little_endian_order();


#ifdef __cplusplus
}
#endif


#endif /* EPR_SWAP_H_INCL */
