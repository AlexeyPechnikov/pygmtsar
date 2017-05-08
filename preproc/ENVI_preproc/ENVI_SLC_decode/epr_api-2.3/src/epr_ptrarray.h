/*
 * $Id: epr_ptrarray.h,v 1.1.1.1 2004-10-28 19:22:23 norman Exp $
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

#ifndef EPR_PTRARRAY_H_INCL
#define EPR_PTRARRAY_H_INCL

#ifdef __cplusplus
extern "C" 
{
#endif

#include <stdlib.h>

/**
 * The <code>EPR_PtrArray</code> structure represents a dynamic
 * array of pointers.
 */
struct EPR_PtrArray
{
	/** The current capacity */
    unsigned int capacity;
    /** The current length */
    unsigned int length;
    /** The pointer elements */
    void** elems;
};

typedef struct EPR_PtrArray EPR_SPtrArray;

/**
 * Creates a new dynamic pointer array instance.
 *
 * @param capacity the initial capacity
 * @return a new dynamic pointer array instance with the given capacity or
 *        <code>NULL</code> if memory could not be allocated
 */
EPR_SPtrArray* epr_create_ptr_array(unsigned int capacity);


/**
 * Frees the memory allocated through the given dynamic pointer array.
 *
 * <p> After calling this function the give record pointer array gets
 * invalid and should not be used anymore.
 * 
 * @param ptr_array the pointer array to be released, if <code>NULL</code>
 *        the function immediately returns
 */
void epr_free_ptr_array(EPR_SPtrArray* ptr_array);


/**
 * Special application of the <code>epr_free_ptr_array</code> for
 * arrays that contain dynamically allocated strings (type <code>char*</code>).
 *
 * <p>For each element in the given array the <code>epr_free_string</code> 
 * function is called.
 *
 * @param char_ptr_array an array containing strings
 */
void epr_free_char_ptr_array(EPR_SPtrArray* char_ptr_array);

/**
 * Adds a new pointer to the given pointer array.
 * The function automatically grows the array if necessary.
 *
 * @param ptr_array the pointer array to which to add the new element,
 *        must not be NULL.
 * @param elem the element to be added
 * @return zero for success, an error code otherwise
 */
int epr_add_ptr_array_elem(EPR_SPtrArray* ptr_array, void* elem);


/**
 * Grows the given pointer array so that is has the given capacity.
 * The length of the array is not touched by this function.
 *
 * @param ptr_array the pointer array to which to add the new element,
 *        must not be NULL.
 * @param capacity the new capacity
 * @return zero for success, an error code otherwise
 */
int epr_grow_ptr_array(EPR_SPtrArray* ptr_array, unsigned int capacity);


/**
 * Returns the length of the given pointer array.
 * @param ptr_array the pointer array, must not be NULL.
 * @return the length of the given pointer array
 */
unsigned int epr_get_ptr_array_length(const EPR_SPtrArray* ptr_array);


/**
 * Returns the capacity of the given pointer array.
 * @param ptr_array the pointer array, must not be NULL.
 * @return the capacity of the given pointer array
 */
unsigned int epr_get_ptr_array_capacity(const EPR_SPtrArray* ptr_array);


/**
 * Gets the element with the specified index of the given pointer array.
 * @param ptr_array the pointer array, must not be NULL.
 * @param index the zero-based index, must be less than the array's length
 * @return the element at the given index
 */
void* epr_get_ptr_array_elem_at(const EPR_SPtrArray* ptr_array, unsigned int index);


#ifdef __cplusplus
} /* extern "C" */
#endif
#endif /* #ifndef EPR_PTRARRAY_H_INCL */
