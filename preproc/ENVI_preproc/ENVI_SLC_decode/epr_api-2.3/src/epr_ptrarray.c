/*
 * $Id: epr_ptrarray.c,v 1.1.1.1 2004-10-28 19:22:23 norman Exp $
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

#include <assert.h>
#include <errno.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "epr_api.h"
#include "epr_ptrarray.h"
#include "epr_string.h"


EPR_SPtrArray* epr_create_ptr_array(unsigned int capacity)
{
    EPR_SPtrArray* ptr_array = NULL;
    void** elems = NULL;
    
    elems = (void**) calloc(capacity, sizeof (void*));
    if (elems == NULL) {
        return NULL;
    }
    
    ptr_array = (EPR_SPtrArray*) calloc(1, sizeof (EPR_SPtrArray));
    if (ptr_array == NULL) {
        free(elems);
        return NULL;
    }

    ptr_array->capacity = capacity;
    ptr_array->length = 0;
    ptr_array->elems = elems;
    return ptr_array;
}

void epr_free_ptr_array(EPR_SPtrArray* ptr_array)
{
    if (ptr_array == NULL)
        return;
    if (ptr_array->elems != NULL)
        free(ptr_array->elems);
    ptr_array->capacity = 0;
    ptr_array->length = 0;
    ptr_array->elems = NULL;
    free(ptr_array);
}



void epr_free_char_ptr_array(EPR_SPtrArray* char_ptr_array)
{
    uint i;
    for (i = 0; i < char_ptr_array->length; i++) {
        epr_free_string((char*) char_ptr_array->elems[i]);
    }
    epr_free_ptr_array(char_ptr_array);
}


int epr_add_ptr_array_elem(EPR_SPtrArray* ptr_array, void* elem)
{
    assert(ptr_array != NULL);
    assert(ptr_array->elems != NULL);

    if (ptr_array->length + 1 > ptr_array->capacity) {
        int status = epr_grow_ptr_array(ptr_array, 2 * ptr_array->capacity);
        if (status != 0)
            return status;
    }

    ptr_array->elems[ptr_array->length++] = elem;
    return e_err_none;
}


int epr_grow_ptr_array(EPR_SPtrArray* ptr_array, unsigned int capacity)
{
    void* elems = NULL;
    
    assert(ptr_array != NULL);
    assert(capacity >= ptr_array->capacity);

    if (capacity == ptr_array->capacity)
        return e_err_none;

    elems = (void**) realloc(ptr_array->elems, capacity * sizeof (void*));
    if (elems == NULL)
        return e_err_out_of_memory;

    memset(((char*)elems) + ptr_array->length * sizeof (void*), 
           0, 
           (capacity - ptr_array->length) * sizeof (void*));

    ptr_array->capacity = capacity;
    ptr_array->elems = elems;
    return e_err_none;
}


unsigned int epr_get_ptr_array_length(const EPR_SPtrArray* ptr_array)
{
    assert(ptr_array != NULL);
    return ptr_array->length;
}


void* epr_get_ptr_array_elem_at(const EPR_SPtrArray* ptr_array, unsigned int index)
{
    assert(ptr_array != NULL);
    assert(ptr_array->elems != NULL);
    assert(index < ptr_array->length);
    return ptr_array->elems[index];
}


