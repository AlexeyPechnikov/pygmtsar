/*
 * $Id: epr_string.c,v 1.1.1.1 2004-10-28 19:22:23 norman Exp $
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
#include <ctype.h>
#include <string.h>

#include "epr_api.h"
#include "epr_core.h"
#include "epr_string.h"

char* epr_assign_string(char** str_clone, const char* str)
{
    assert(str_clone != NULL);
    epr_free_string(*str_clone);
    *str_clone = epr_clone_string(str);
    return *str_clone;
}

char* epr_create_string(unsigned int length)
{
    return (char*) calloc(length + 1, sizeof (char));
}


char* epr_clone_string(const char* str)
{
    char* str_clone = NULL;
    if (str != NULL) {
        str_clone = epr_create_string(strlen(str));
        strcpy(str_clone, str);
    }
    return str_clone;
}


void epr_cut_string(char** sub_str, const char* str, int start, int length)
{
    *sub_str = epr_create_string(length);
    strncpy(*sub_str, str + start, length);
    *sub_str[length] = '\0';
}


char* epr_sub_string(const char* str, int start, int length)
{
    char* str_sub_string = NULL;
    if (str != NULL) {
        str_sub_string = epr_create_string(length);
        strncpy(str_sub_string, str + start, length);
        str_sub_string[length] = '\0';
    }
    return str_sub_string;
}


/*
   Function:    epr_get_last_err_code
   Access:      private API implementation helper
   Changelog:   2002/01/05  nf initial version
 */
/**
 * Compares the two given names and returns <code>TRUE</code> if
 * they are equal ignoring the case of each letter.
 *
 * <p> This function is used to compare names throughout the
 * ENVISAT product reader API.
 *
 * @param name1 the first name, must not be NULL
 * @param name2 the second name, must not be NULL
 * @return  <code>TRUE</code> if the names are equal,
 *          <code>FALSE</code> otherwise
 */
epr_boolean epr_equal_names(const char* name1, const char* name2)
{
    assert(name1 != NULL);
    assert(name2 != NULL);

    return stricmp(name1, name2) == 0;
}


void epr_free_string(char* str)
{
    if (str == NULL)
        return;
    free(str);
}

/*
   Function:    epr_str_tok
   Access:      private API implementation helper
   Changelog:   2002/01/10  mp initial version
 */
/**
 * Findes substrings between separators.
 *
 * @param str the string to search
 * @param seps the separator simbols string
 * @param pos position with a search begin
 *
 * @return the next substring (or own string) of the found name or
 *         <code>(uint)NULL</code> if an error occured.
 */
char* epr_str_tok(const char* str, const char* seps, int* pos)
{
    char* token = NULL;
    int i, old_pos;
    int token_len = 0;

    assert(str != NULL);

    if (*pos >= (int)strlen(str)) return NULL;

    old_pos = *pos;

    for (i = *pos; str[i] != '\0'; i++) {
        if (strchr(seps, str[i]) != NULL)
        {
            token_len = i - *pos;
            token = epr_create_string(token_len);
            strncpy(token, str + *pos, token_len);
            token[token_len] = '\0';
            *pos = i + 1;
            return token;
        }
    }
    if (strlen(str) > 0) {
        if (old_pos == 0) {
            *pos = i + 1;
            token = epr_clone_string(str);
            return token;
        } else if (old_pos > 0) {
            token_len = i - old_pos;
            token = epr_create_string(token_len);
            strncpy(token, str + *pos, token_len);
            token[token_len] = '\0';
            *pos = (int)strlen(str);
            return token;
        }
    }
    return NULL;
}

/*
   Function:    epr_str_tok_tok
   Access:      private API implementation helper
   Changelog:   2002/01/10  mp initial version
 */
/**
 * Findes substrings between double separators.
 *
 * @param str the string to search
 * @param seps the separator simbols string
 * @param pos position with a search begin
 *
 * @return the next substring (or own string) of the found name or
 *         <code>(uint)NULL</code> if an error occured.
 */
char* epr_str_tok_tok(const char* str, const char* seps, const char* exceptions, uint* pos)
{
    char* token = NULL;
    uint i, old_pos;
    uint token_len = 0;

    assert(str != NULL);

    if (*pos >= (uint)strlen(str)) return NULL;

    old_pos = *pos;

    for (i = *pos; str[i] != '\0'; i++) {
        if (((strchr(seps, str[i]) != NULL)
                && (i == 0))
            || ((strchr(seps, str[i]) != NULL)
                && (i > 0)
                && (strchr(exceptions, str[i - 1]) == NULL))) {
            token_len = i - *pos;
            token = epr_create_string(token_len);
            strncpy(token, str + *pos, token_len);
            token[token_len] = '\0';
            *pos = i + 1;
            return token;
        }
    }
    if (strlen(str) > 0) {
        if (old_pos == 0) {
            *pos = i + 1;
            token = epr_clone_string(str);
            return token;
        } else if (old_pos > 0) {
            token_len = i - old_pos;
            token = epr_create_string(token_len);
            strncpy(token, str + *pos, token_len);
            token[token_len] = '\0';
            *pos = (int)strlen(str);
            return token;
        }
    }
    return NULL;
}


int epr_find_first_not_white(const char* str)
{
    int i;
    char white[] = {" "};

    for (i = 0; i < (int)strlen(str); i ++) {
        if (strchr(white, str[i]) == NULL) return i;
    }
    return (int)strlen(str);
}

int epr_find_last_not_white(const char* str)
{
    int i;
    char white[] = {" "};

    if ((int)strlen(str) == 0)
        return -1;

    for (i = (int)strlen(str) - 1; i >=  0; i --) {
        if (strchr(white, str[i]) == NULL) return i;
    }
    return -1;
}


char* epr_trim_string(char* str)
{
    int i, i1, i2, n;

    assert(str != NULL);

    n = strlen(str);
    if (n == 0)
        return str;

    i1 = -1;
    for (i = 0; str[i] != '\0'; i++) {
        if (!isspace(str[i])) {
            i1 = i;
            break;
        }
    }

    if (i1 == -1) {
        str[0] = '\0';
        return str;
    }

    i2 = -1;
    for (i = n - 1; i >= 0; i--) {
        if (!isspace(str[i])) {
            i2 = i;
            break;
        }
    }

    assert(i1 > -1 && i2 > -1);

    if (i1 > 0) {
        int j = 0;
        for (i = i1; i <= i2; i++) {
            str[j] = str[i];
            j++;
        }
    }

    str[i2 - i1 + 1] = '\0';

    return str;
}


char* epr_strip_string_r(char* str)
{
    int i, i1, n;

    assert(str != NULL);

    n = strlen(str);
    if (n == 0) {
        return str;
    }

    i1 = -1;
    for (i = n - 1; i >= 0; i--) {
        if (33 <= str[i] && str[i] <= 126) {
            i1 = i;
            break;
        }
    }

    if (i1 == -1) {
        str[0] = '\0';
    } else {
        str[i1+1] = '\0';
    }

    return str;
}

int epr_if_no_letters(const char* str)
{
    int i;
    char ciffer[] = {"0123456789+- "};

    for (i = 0; i < (int)strlen(str); i ++) {
        if (strchr(ciffer, str[i]) == NULL) return 0;
    }
    return 1;
}


int epr_get_positive_int(const char* str)
{
    int i;
    char ciffer[] = {"0123456789 "};

    for (i = 0; i < (int)strlen(str); i ++) {
        if (strchr(ciffer, str[i]) == NULL) return -1;
    }
    return atoi(str);
}


/*'Verdacht': if float-or-double*/
int epr_numeral_suspicion(const char* str)
{
    int i;
    char ciffer[] = {"0123456789+- .eE"};

    for (i = 0; i < (int)strlen(str); i ++) {
        if (strchr(ciffer, str[i]) == NULL) return 0;
    }
    /*if YES*/
    return 1;
}


int epr_if_no_scaling(const char* str)
{
    int result;

    result = strrchr(str, EPR_IRRELEVANCE_SYMBOL) - str + 1;
    if (result == 1) return 1;
    return 0;
}


void epr_free_and_null_string(char** str)
{
    if (*str == NULL) return;
    epr_free_string(*str);
    *str = NULL;
}


int stricmp(const char* s1, const char* s2)
{
    int d;
    assert(s1 != NULL);
    assert(s2 != NULL);
    while (*s1 != '\0' && *s2 != '\0') {
        d = tolower(*s1) - tolower(*s2);
        if (d != 0) {
            break;
        }
        s1++;
        s2++;
    }
    d = tolower(*s1) - tolower(*s2);
    return d;
}
