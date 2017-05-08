/*
 * $Id: epr_string.h,v 1.1.1.1 2004-10-28 19:22:23 norman Exp $
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

#ifndef EPR_STRING_H_INCL
#define EPR_STRING_H_INCL

#ifdef __cplusplus
extern "C"
{
#endif

#include <stdlib.h>
#include <string.h>
#include <assert.h>

char* epr_assign_string(char** str_clone, const char* str);
char* epr_create_string(unsigned int length);
char* epr_clone_string(const char* str);
void epr_cut_string(char** sub_str, const char* str, int start, int length);
char* epr_sub_string(const char* str, int start, int length);
epr_boolean epr_equal_names(const char* name1, const char* name2);
void epr_free_string(char* str);
char* epr_str_tok(const char* str, const char* seps, int* pos);
char* epr_str_tok_tok(const char* str, const char* seps, const char* exceptions, uint* pos);
int epr_find_first_not_white(const char* str);
int epr_find_last_not_white(const char* str);
char* epr_trim_string(char* str);
char* epr_strip_string_r(char* str);
int epr_if_no_letters(const char* str);
int epr_numeral_suspicion(const char* str);
int epr_get_positive_int(const char* str);
void epr_free_and_null_string(char** str);
int epr_if_no_scaling(const char* str);

/**
 * Non-ANSI string compare (ignores case).
 */
int stricmp(const char* s1, const char* s2);

#ifdef __cplusplus
} /* extern "C" */
#endif
#endif /* #ifndef EPR_STRING_H_INCL */
