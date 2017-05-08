/*
 * $Id: epr_bitmask.c,v 1.2 2004-10-28 20:10:57 norman Exp $
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

#include <ctype.h>

#include "epr_api.h"
#include "epr_core.h"
#include "epr_string.h"
#include "epr_ptrarray.h"
#include "epr_swap.h"
#include "epr_field.h"
#include "epr_record.h"
#include "epr_param.h"
#include "epr_dsd.h"
#include "epr_msph.h"
#include "epr_band.h"
#include "epr_bitmask.h"

#include "epr_dddb.h"

void epr_resolve_bm_ref(EPR_SBmEvalContext* context, EPR_SBmTerm* term);



EPR_SBmEvalContext* epr_create_bm_eval_context(EPR_SProductId* product_id,
                                               int offset_x,
                                               int offset_y,
                                               EPR_SRaster* bitmask_raster)
{
    EPR_SBmEvalContext* context;

    context = (EPR_SBmEvalContext*) calloc(1, sizeof (EPR_SBmEvalContext));
    if (context == NULL) {
        epr_set_err(e_err_out_of_memory,
                    "epr_create_bm_eval_context: out of memory");
        return NULL;
    }
    context->product_id = product_id;
    context->offset_x = offset_x;
    context->offset_y = offset_y;
    context->bitmask_raster = bitmask_raster;
    context->flag_band_ids = epr_create_ptr_array(4);
    context->flag_rasters = epr_create_ptr_array(4);
    return context;
}


void epr_free_bm_eval_context(EPR_SBmEvalContext* context)
{
    EPR_SRaster* flag_raster;
    uint flag_index;

    if (context == NULL) {
        return;
    }

    if (context->flag_band_ids != NULL) {
        for (flag_index = 0; flag_index < context->flag_band_ids->length; flag_index++) {
            /*
             * Note that the release of band ID's is handled by the epr_close_product() function,
             * The rasters need to be released, because they are internally used by this
             * context.
             */
            flag_raster = (EPR_SRaster*)context->flag_rasters->elems[flag_index];
            epr_free_raster(flag_raster);
        }
        epr_free_ptr_array(context->flag_band_ids);
        epr_free_ptr_array(context->flag_rasters);
    }
    free(context);
}

/**
 * Reads bit-mask pixels of the given product for the given bit-mask expression
 * for the given region offset and raster.
 *
 * @param product_id the product ID
 * @param bm_expr the bit-mask expression
 * @param offset_x X-coordinate in pixel co-ordinates (zero-based) of the upper right raster corner to be searched for
 * @param offset_y Y-coordinate in pixel co-ordinates (zero-based) of the upper right raster corner to be searched for
 * @param raster the instance to the buffer information was used
 *
 * @return zero for success, an error code otherwise
 */
int epr_read_bitmask_raster(EPR_SProductId* product_id,
                            const char* bm_expr,
                            int offset_x,
                            int offset_y,
                            EPR_SRaster* bm_raster)
{
    EPR_SBmEvalContext* context;
    EPR_SBmTerm* term;
    uint x, y;
    uint pos;
    uchar* bm_buffer = NULL;
    EPR_EErrCode errcode;

    epr_clear_err();

    if (bm_raster->data_type != e_tid_uchar && bm_raster->data_type != e_tid_char) {
        epr_set_err(e_err_illegal_data_type,
            "epr_read_bitmask_raster: illegal raster datatype; must be 'char' or 'uchar'");
        return e_err_illegal_data_type;
    }

    bm_buffer = (uchar*) bm_raster->buffer;
    if (bm_buffer == NULL) {
         epr_set_err(e_err_out_of_memory,
             "epr_read_bitmask_raster: false memory allocation for a raster buffer");
        return e_err_out_of_memory;
    }

    context = epr_create_bm_eval_context(product_id, offset_x, offset_y, bm_raster);
    if (context == NULL) {
         epr_set_err(e_err_illegal_arg,
             "epr_read_bitmask_raster: the context cannot be created");
        return e_err_illegal_arg;
    }

    term = epr_parse_bm_expr_str(bm_expr);

    if (term == NULL) {
         epr_set_err(e_err_illegal_arg,
             "epr_read_bitmask_raster: the term was not build");
        return e_err_illegal_arg;
    }

    pos = 0;

    epr_clear_err();

    errcode = epr_get_last_err_code();
    for (y = 0; y < bm_raster->raster_height; y++) {
        for (x = 0; x < bm_raster->raster_width; x++) {
            bm_buffer[pos] = (uchar) epr_eval_bm_term(context, term, x, y);
            pos++;
            errcode = epr_get_last_err_code();
            if (errcode != 0) {
                break;
            }
        }
        if (errcode != 0) {
            break;
        }
    }

    epr_free_bm_term(term);
    epr_free_bm_eval_context(context);

    return errcode;
}

/**
 * Evaluates the given bitmask expression.
 *
 * @param term the bitmask term
 * @param x the pixel's x co-ordinate
 * @param y the pixel's y co-ordinate
 */
epr_boolean epr_eval_bm_term(EPR_SBmEvalContext* context, EPR_SBmTerm* term, int x, int y) {

	uint	temp;
	uint	eval;
    if (term == NULL) {
        return FALSE;
    }

    switch (term->op_code) {
    case BMT_REF:
		{
			EPR_SRaster* flag_raster = term->op.ref.flag_raster;
            uint flag_mask = term->op.ref.flag_mask;

            if (flag_raster == NULL) {
                epr_resolve_bm_ref(context, term);
                flag_raster = term->op.ref.flag_raster;
                flag_mask = term->op.ref.flag_mask;
                if (flag_raster == NULL) {
                    return FALSE;
                }
            }

            assert(flag_raster != NULL);
            assert(flag_mask != FLAG_MASK_NOT_COMPUTED);
			temp = epr_get_pixel_as_uint(flag_raster, x, y);
			eval = temp & flag_mask;
	        return (eval) == flag_mask;
        }
    case BMT_AND:
		{
			if (!epr_eval_bm_term(context, term->op.binary.arg1, x, y)) {
                return FALSE;
			}
            return epr_eval_bm_term(context, term->op.binary.arg2, x, y);
		}
    case BMT_OR:
		{
            if (epr_eval_bm_term(context, term->op.binary.arg1, x, y)) {
                return TRUE;
			}
            return epr_eval_bm_term(context, term->op.binary.arg2, x, y);
        }
    case BMT_NOT:
		{
            return !epr_eval_bm_term(context, term->op.unary.arg, x, y);
        }
    default:
        assert(0);
        return FALSE;
    }
}


void epr_resolve_bm_ref(EPR_SBmEvalContext* context, EPR_SBmTerm* term) {
    const char* band_name = term->op.ref.band_name;
    const char* flag_name = term->op.ref.flag_name;
    EPR_SBandId* flag_band_id = NULL;
    uint flag_band_index = (uint) -1;
    uint band_index = 0;
    uint num_bands = context->flag_band_ids->length;
    EPR_SRaster* flag_raster = NULL;
    uint flag_mask = 0;
	uint flag_computed = 0;

    /* Find the corresponding flag_band_id for band_name */
    for (band_index = 0; band_index < num_bands; band_index++) {
        flag_band_id = (EPR_SBandId*) context->flag_band_ids->elems[band_index];
        if (epr_equal_names(band_name, flag_band_id->band_name)) {
            flag_band_index = band_index;
            break;
        }
    }
    /* flag_band_id found? */
    if (flag_band_index != (uint) -1) {
        /* Yes, found: get flag_band_id  and the corresponding raster */
        flag_band_id = (EPR_SBandId*)(context->flag_band_ids->elems[flag_band_index]);
        flag_raster = (EPR_SRaster*)(context->flag_rasters->elems[flag_band_index]);
    } else {
        /* Not found: get flag_band_id from product and load the corresponding raster */
        flag_band_id = epr_get_band_id(context->product_id, band_name);
        if (flag_band_id != NULL) {
            flag_raster = epr_create_compatible_raster(flag_band_id,
                                            context->bitmask_raster->source_width,
                                            context->bitmask_raster->source_height,
                                            context->bitmask_raster->source_step_x,
                                            context->bitmask_raster->source_step_y);
            epr_read_band_raster(flag_band_id,
                                 context->offset_x,
                                 context->offset_y,
                                 flag_raster);
            /* register flag_band_id and flag_raster for later use */
            epr_add_ptr_array_elem(context->flag_band_ids, flag_band_id);
            epr_add_ptr_array_elem(context->flag_rasters, flag_raster);
        } else {
            epr_set_err(e_err_flag_not_found, "flags band not found");
            return;
        }
    }

    /* Now, compute flag_mask */

    /* Does the band have a flag coding? */
    if (flag_band_id->flag_coding != NULL) {
        /* Yes, now find flag definition for flag_name */
        EPR_SFlagDef* flag_def = NULL;
        uint flag_def_index;
        for (flag_def_index = 0; flag_def_index < flag_band_id->flag_coding->length; flag_def_index++) {
            flag_def = (EPR_SFlagDef*) flag_band_id->flag_coding->elems[flag_def_index];
            if (epr_equal_names(flag_def->name, flag_name)) {
				flag_computed = 1;
                flag_mask |= flag_def->bit_mask;
                /* TODO!!! */
                break;
            }
        }
    }
    if (flag_computed == 0) {
		flag_mask = FLAG_MASK_NOT_COMPUTED;
        epr_set_err(e_err_flag_not_found, "flag not found");
    }
    term->op.ref.flag_mask = flag_mask;
    term->op.ref.flag_raster = flag_raster;
}

/**
 * Parses the bit-mask expression given as character string.
 *
 * @param bm_expr the bit-mask expression given as character string
 * @return the bit-mask term tree representing the bit-mask expression
 * @throws BitmaskExpressionParseException if the given code could not be epr_parse'd
 * @throws IOException if an I/O error occurs
 */
EPR_SBmTerm* epr_parse_bm_expr_str(const char* bm_expr) {

    EPR_SBmTerm* term;
    EPR_SParseInfo parse_info;

    parse_info.bm_expr = bm_expr;
    parse_info.bm_expr_pos = 0;
    parse_info.pushed_back = FALSE;
    parse_info.token = NULL;
    parse_info.err_message = NULL;

    term = epr_parse_bm_expr(&parse_info, FALSE);

    epr_free_string(parse_info.token);
    parse_info.token = NULL;

    if (epr_is_bm_expr_error(&parse_info)) {
        char tmp[256] = {"bitmap-expression error: "};
        strcat(tmp, parse_info.err_message);
        epr_set_err(e_err_invalid_value, tmp);
    }

    return term;
}


EPR_SPtrArray* epr_create_flag_coding(EPR_SProductId* product_id, const char* flag_coding_name)
{
	int num_descr;
    int i, j;
    const struct FlagDescriptorTable* fc_tables;
    int fct_index;
    EPR_SPtrArray* flag_coding = NULL;

	if (product_id == NULL) {
        epr_set_err(e_err_null_pointer,
                    "epr_create_flag_coding: product_id must not be NULL");
        return NULL;
    }

    /* @DDDB */

    fc_tables = dddb_flag_coding_tables;
    fct_index = -1;
    for (i = 0; i < EPR_NUM_FLAG_CODING_TABLES; i++) {
        const char* id = fc_tables[i].name;
        if (epr_equal_names(id, flag_coding_name)) {
            fct_index = i;
            break;
        }
    }
    if (fct_index == -1) {
        epr_set_err(e_err_null_pointer,
                    "epr_create_flag_coding: unknown flag coding");
        return NULL;
    }

    flag_coding = epr_create_ptr_array(16);
    num_descr = fc_tables[fct_index].num_descriptors;
	for (i = 0; i < num_descr; i++) {
		EPR_SFlagDef* flag_def = (EPR_SFlagDef*) calloc(1, sizeof (EPR_SFlagDef));
		if (flag_def == NULL) {
			epr_set_err(e_err_out_of_memory,
						"epr_create_flag_coding: out of memory");
			return NULL;
		}
        /* 1: flag_name */
        epr_assign_string(&flag_def->name, fc_tables[fct_index].descriptors[i].id);
        if (flag_def->name == NULL) {
            epr_set_err(e_err_out_of_memory, "epr_get_flag_coding: out of memory");
            epr_free_flag_def(flag_def);
            return NULL;
        }
        /* 2: dataset_name */
        /* flag_def->bit_index = (uint)fc_tables[fct_index].descriptors[i].bit_index; */
		flag_def->bit_mask = 0;
		for (j = 0; j < fc_tables[fct_index].descriptors[i].num_indices; j++) {
			flag_def->bit_mask |= (1 << fc_tables[fct_index].descriptors[i].bit_indices[j]);
		}
         /* 3: sample_offset */
        epr_assign_string(&flag_def->description, fc_tables[fct_index].descriptors[i].description);

		epr_add_ptr_array_elem(flag_coding, flag_def);
	}
    return flag_coding;
}


void epr_free_flag_coding(EPR_SPtrArray* flag_coding)
{
    EPR_SFlagDef* flag_def = NULL;
    uint flag_index;

    if (flag_coding == NULL) {
        return;
    }

    for (flag_index = 0; flag_index < flag_coding->length; flag_index++)
    {
        flag_def = (EPR_SFlagDef*) flag_coding->elems[flag_index];
        epr_free_flag_def(flag_def);
        flag_def = NULL;
    }
    epr_free_ptr_array(flag_coding);
}


EPR_SFlagDef* epr_create_flag_def()
{
    EPR_SFlagDef* flag_def = NULL;

    flag_def = (EPR_SFlagDef*) calloc(1, sizeof (EPR_SFlagDef));
    if (flag_def == NULL) {
        epr_set_err(e_err_out_of_memory,
                    "epr_create_flag_def: out of memory");
        return NULL;
    }
    flag_def->magic = EPR_MAGIC_FLAG_DEF;

    return flag_def;
}


void epr_free_flag_def(EPR_SFlagDef* flag_def)
{
    if (flag_def == NULL)
        return;

    epr_free_string(flag_def->name);
    flag_def->name = NULL;
    flag_def->bit_mask = 0;
    epr_free_string(flag_def->description);
    flag_def->description = NULL;

    free(flag_def);
}


EPR_SBmTerm* epr_parse_bm_expr(EPR_SParseInfo* parse_info, epr_boolean term_required) {
    return epr_parse_bm_OR_expr(parse_info, term_required);
}



EPR_SBmTerm* epr_parse_bm_OR_expr(EPR_SParseInfo* parse_info, epr_boolean term_required) {
    EPR_SBmTerm* term1 = epr_parse_bm_AND_expr(parse_info, term_required);
    if (term1 == NULL) {
        return NULL;
    }

    while (!epr_is_bm_expr_error(parse_info)) {
        epr_next_bm_expr_token(parse_info);
        if (epr_is_bm_OR_keyword(parse_info) || epr_is_bm_OR_operator(parse_info)) {
            EPR_SBmTerm* term2 = epr_parse_bm_OR_expr(parse_info, TRUE);
            term1 = epr_create_bm_OR_term(term1, term2);
        } else {
            epr_push_back_bm_expr_token(parse_info);
            break;
        }
    }

    return term1;
}

EPR_SBmTerm* epr_parse_bm_AND_expr(EPR_SParseInfo* parse_info, epr_boolean term_required) {
    EPR_SBmTerm* term1 = epr_parse_bm_unary_expr(parse_info, term_required);
    if (term1 == NULL) {
        return NULL;
    }

    while (!epr_is_bm_expr_error(parse_info)) {
        epr_next_bm_expr_token(parse_info);
        if (epr_is_bm_AND_keyword(parse_info) || epr_is_bm_AND_operator(parse_info)) {
            EPR_SBmTerm* term2 = epr_parse_bm_AND_expr(parse_info, TRUE);
            term1 = epr_create_bm_AND_term(term1, term2);
        } else {
            epr_push_back_bm_expr_token(parse_info);
            break;
        }
    }

    return term1;
}


EPR_SBmTerm* epr_parse_bm_unary_expr(EPR_SParseInfo* parse_info, epr_boolean term_required)  {

    EPR_SBmTerm* term = NULL;

    epr_next_bm_expr_token(parse_info);
    if (epr_is_bm_NOT_keyword(parse_info) || epr_is_bm_NOT_operator(parse_info)) {
        term = epr_parse_bm_unary_expr(parse_info, TRUE);
        term = epr_create_bm_NOT_term(term);
    } else {
        epr_push_back_bm_expr_token(parse_info);
        term = epr_parse_bm_primary_expr(parse_info, term_required);
    }

    return term;
}

EPR_SBmTerm* epr_parse_bm_primary_expr(EPR_SParseInfo* parse_info, epr_boolean term_required) {

    EPR_SBmTerm* term = NULL;

    epr_next_bm_expr_token(parse_info);
    if (epr_get_token_char(parse_info) == '(') {
        term = epr_parse_bm_expr(parse_info, TRUE);
        epr_next_bm_expr_token(parse_info);
        if (epr_get_token_char(parse_info) != ')') {
            epr_set_bm_expr_error(parse_info, "')' expected");
        }
    } else if (epr_is_bm_name_token(parse_info)) {
        char* ds_name = epr_consume_token(parse_info);
        epr_next_bm_expr_token(parse_info);
        if (epr_get_token_char(parse_info) == '.') {
            epr_next_bm_expr_token(parse_info);
            if (epr_is_bm_name_token(parse_info)) {
                char* flag_name = epr_consume_token(parse_info);
                term = epr_create_bm_REF_term(ds_name, flag_name);
            } else {
                epr_set_bm_expr_error(parse_info, "flag name expected");
            }
        } else {
            epr_set_bm_expr_error(parse_info, "'.' expected");
        }
    } else if (epr_is_bm_EOS_token(parse_info)) {
        if (term_required) {
            epr_set_bm_expr_error(parse_info, "operator or flag name expected");
        }
    } else {
        epr_set_bm_expr_error(parse_info, "operator or flag name expected");
    }

    return term;
}

epr_boolean epr_is_bm_OR_keyword(EPR_SParseInfo* parse_info) {
    return epr_is_bm_name_token(parse_info) && stricmp("or", parse_info->token) == 0;
}

epr_boolean epr_is_bm_AND_keyword(EPR_SParseInfo* parse_info) {
    return epr_is_bm_name_token(parse_info) && stricmp("and", parse_info->token) == 0;
}

epr_boolean epr_is_bm_NOT_keyword(EPR_SParseInfo* parse_info) {
    return epr_is_bm_name_token(parse_info) && stricmp("not", parse_info->token) == 0;
}

epr_boolean epr_is_bm_AND_operator(EPR_SParseInfo* parse_info) {
    return epr_get_token_char(parse_info) == '&';
}

epr_boolean epr_is_bm_OR_operator(EPR_SParseInfo* parse_info) {
    return epr_get_token_char(parse_info) == '|';
}

epr_boolean epr_is_bm_NOT_operator(EPR_SParseInfo* parse_info) {
    return epr_get_token_char(parse_info) == '!';
}

epr_boolean epr_is_bm_name_token(EPR_SParseInfo* parse_info) {
    return parse_info->token_type == BME_NAME && parse_info->token != NULL;
}

epr_boolean epr_is_bm_EOS_token(EPR_SParseInfo* parse_info) {
    return parse_info->token_type == BME_EOS;
}

epr_boolean epr_is_bm_expr_error(EPR_SParseInfo* parse_info) {
    return parse_info->err_message != NULL;
}

int epr_get_token_char(EPR_SParseInfo* parse_info) {
    if (parse_info->token_type == BME_SPECIAL && parse_info->token != NULL) {
        return parse_info->token[0];
    }
    return '\0';
}

char* epr_consume_token(EPR_SParseInfo* parse_info) {
    char* token = parse_info->token;
    /* Prevent from being released by epr_free_string() */
    parse_info->token = NULL;
    parse_info->token_type = BME_UNKNOWN;
    parse_info->pushed_back = FALSE;
    return token;

}

void epr_next_bm_expr_token(EPR_SParseInfo* parse_info) {
    if (parse_info->pushed_back) {
        parse_info->pushed_back = FALSE;
        return;
    }
    epr_free_string(parse_info->token);
    parse_info->token_type = epr_tokenize_bm_expr(parse_info->bm_expr,
                                                  &(parse_info->bm_expr_pos),
                                                  &(parse_info->token));
}

void epr_push_back_bm_expr_token(EPR_SParseInfo* parse_info) {
    parse_info->pushed_back = TRUE;
}

void epr_set_bm_expr_error(EPR_SParseInfo* parse_info, const char* message) {
    static char msg_buf[2048];

    epr_push_back_bm_expr_token(parse_info);

    if (message != NULL) {
        if (!epr_is_bm_EOS_token(parse_info)) {
            sprintf(msg_buf, "%s, but found token '%s'", message, parse_info->token);
        } else {
            sprintf(msg_buf, "%s, but found 'end-of-string'", message);
        }
    } else {
        if (!epr_is_bm_EOS_token(parse_info)) {
            sprintf(msg_buf, "unexpected token '%s' found", parse_info->token);
        } else {
            sprintf(msg_buf, "unexpected 'end-of-string' found");
        }
    }

    parse_info->err_message = epr_clone_string(msg_buf);
}


int epr_tokenize_bm_expr(const char* bm_expr, int* bm_expr_pos, char** token)
{
    int pos = *bm_expr_pos;

    while (isspace(bm_expr[pos])) {
        pos++;
    }

    if (bm_expr[pos] == '\0') {
        *bm_expr_pos = pos;
        *token = NULL;
        return BME_EOS;
    }

    if (isalpha(bm_expr[pos]) || bm_expr[pos] == '_') {
        int pos0 = pos;
        size_t len;
        char* tok;

        pos++;
        while (isalnum(bm_expr[pos]) || bm_expr[pos] == '_') {
            pos++;
        }

        len = pos - pos0;
        tok = epr_create_string(len + 1);
        strncpy(tok, bm_expr + pos0, len);
        tok[len] = '\0';

        *token = tok;
        *bm_expr_pos = pos;
        return BME_NAME;
    }

    if (bm_expr[pos] == '(' ||
        bm_expr[pos] == ')' ||
        bm_expr[pos] == '.' ||
        bm_expr[pos] == '&' ||
        bm_expr[pos] == '|' ||
        bm_expr[pos] == '!') {
        char* tok;

        tok = epr_create_string(2);
        tok[0] = bm_expr[pos];
        tok[1] = '\0';

        pos++;

        *token = tok;
        *bm_expr_pos = pos;
        return BME_SPECIAL;
    }

    *token = NULL;
    *bm_expr_pos = pos;
    return BME_UNKNOWN;
}


EPR_SBmTerm* epr_create_bm_term(EPR_EBmOpCode op_code)
{
    EPR_SBmTerm* term = (EPR_SBmTerm*) malloc(sizeof (EPR_SBmTerm));
    term->op_code = op_code;
    return term;
}

EPR_SBmTerm* epr_create_bm_REF_term(char* band_name, char* flag_name)
{
    EPR_SBmTerm* term = epr_create_bm_term(BMT_REF);
    term->op.ref.band_name = band_name;
    term->op.ref.flag_name = flag_name;
    term->op.ref.flag_mask = FLAG_MASK_NOT_COMPUTED; /* not computed */
    term->op.ref.flag_raster = NULL;
    return term;
}

EPR_SBmTerm* epr_create_bm_NOT_term(EPR_SBmTerm* arg)
{
    EPR_SBmTerm* term = epr_create_bm_term(BMT_NOT);
    term->op.unary.arg = arg;
    return term;
}

EPR_SBmTerm* epr_create_bm_OR_term(EPR_SBmTerm* arg1, EPR_SBmTerm* arg2)
{
    EPR_SBmTerm* term = epr_create_bm_term(BMT_OR);
    term->op.binary.arg1 = arg1;
    term->op.binary.arg2 = arg2;
    return term;
}

EPR_SBmTerm* epr_create_bm_AND_term(EPR_SBmTerm* arg1, EPR_SBmTerm* arg2)
{
    EPR_SBmTerm* term = epr_create_bm_term(BMT_AND);
    term->op.binary.arg1 = arg1;
    term->op.binary.arg2 = arg2;
    return term;
}



uint epr_get_pixel_as_uint(const EPR_SRaster* raster, int x, int y)
{
    epr_clear_err();

    switch (raster->data_type) {
    case (e_tid_uchar) :
          return (uint) ((uchar*)raster->buffer)[y * raster->raster_width + x];
    case (e_tid_char)  :
          return (uint) ((char*)raster->buffer)[y * raster->raster_width + x];
    case (e_tid_ushort):
          return (uint) ((ushort*)raster->buffer)[y * raster->raster_width + x];
    case (e_tid_short) :
          return (uint) ((short*)raster->buffer)[y * raster->raster_width + x];
    case (e_tid_uint):
          return (uint) ((uint*)raster->buffer)[y * raster->raster_width + x];
    case (e_tid_int) :
          return (uint) ((int*)raster->buffer)[y * raster->raster_width + x];
    case (e_tid_float) :
          return (uint) ((float*)raster->buffer)[y * raster->raster_width + x];
    case (e_tid_double) :
          return (uint) ((double*)raster->buffer)[y * raster->raster_width + x];
    default:
          return 0;
    }
}

int epr_get_pixel_as_int(const EPR_SRaster* raster, int x, int y)
{
    epr_clear_err();

    switch (raster->data_type) {
    case (e_tid_uchar) :
          return (int) ((uchar*)raster->buffer)[y * raster->raster_width + x];
    case (e_tid_char)  :
          return (int) ((char*)raster->buffer)[y * raster->raster_width + x];
    case (e_tid_ushort):
          return (int) ((ushort*)raster->buffer)[y * raster->raster_width + x];
    case (e_tid_short) :
          return (int) ((short*)raster->buffer)[y * raster->raster_width + x];
    case (e_tid_uint):
          return (int) ((uint*)raster->buffer)[y * raster->raster_width + x];
    case (e_tid_int) :
          return (int) ((int*)raster->buffer)[y * raster->raster_width + x];
    case (e_tid_float) :
          return (int) ((float*)raster->buffer)[y * raster->raster_width + x];
    case (e_tid_double) :
          return (int) ((double*)raster->buffer)[y * raster->raster_width + x];
    default:
          return 0;
    }
}


float epr_get_pixel_as_float(const EPR_SRaster* raster, int x, int y)
{
    epr_clear_err();

    switch (raster->data_type) {
    case (e_tid_uchar) :
          return (float) ((uchar*)raster->buffer)[y * raster->raster_width + x];
    case (e_tid_char)  :
          return (float) ((char*)raster->buffer)[y * raster->raster_width + x];
    case (e_tid_ushort):
          return (float) ((ushort*)raster->buffer)[y * raster->raster_width + x];
    case (e_tid_short) :
          return (float) ((short*)raster->buffer)[y * raster->raster_width + x];
    case (e_tid_uint):
          return (float) ((uint*)raster->buffer)[y * raster->raster_width + x];
    case (e_tid_int) :
          return (float) ((int*)raster->buffer)[y * raster->raster_width + x];
    case (e_tid_float) :
          return (float) ((float*)raster->buffer)[y * raster->raster_width + x];
    case (e_tid_double) :
          return (float) ((double*)raster->buffer)[y * raster->raster_width + x];
    default:
          return 0;
    }
}

double epr_get_pixel_as_double(const EPR_SRaster* raster, int x, int y)
{
    epr_clear_err();

    switch (raster->data_type) {
    case (e_tid_uchar) :
          return (double) ((uchar*)raster->buffer)[y * raster->raster_width + x];
    case (e_tid_char)  :
          return (double) ((char*)raster->buffer)[y * raster->raster_width + x];
    case (e_tid_ushort):
          return (double) ((ushort*)raster->buffer)[y * raster->raster_width + x];
    case (e_tid_short) :
          return (double) ((short*)raster->buffer)[y * raster->raster_width + x];
    case (e_tid_uint):
          return (double) ((uint*)raster->buffer)[y * raster->raster_width + x];
    case (e_tid_int) :
          return (double) ((int*)raster->buffer)[y * raster->raster_width + x];
    case (e_tid_float) :
          return (double) ((float*)raster->buffer)[y * raster->raster_width + x];
    case (e_tid_double) :
          return (double) ((double*)raster->buffer)[y * raster->raster_width + x];
    default:
          return 0;
    }
}

/**
 * Releases a new bitmask term instance.
 */
void epr_free_bm_term(EPR_SBmTerm* term) {
    if (term == NULL)
        return;

    switch (term->op_code) {
    case BMT_REF:
        epr_free_string(term->op.ref.band_name);
        epr_free_string(term->op.ref.flag_name);
        term->op.ref.band_name = NULL;
        term->op.ref.flag_name = NULL;
        break;
    case BMT_AND: /* OR */ case BMT_OR:
        epr_free_bm_term(term->op.binary.arg1);
        epr_free_bm_term(term->op.binary.arg2);
        term->op.binary.arg1 = NULL;
        term->op.binary.arg2 = NULL;
        break;
    case BMT_NOT:
        epr_free_bm_term(term->op.unary.arg);
        term->op.unary.arg = NULL;
        break;
    default:
        assert(0);
    }

    free(term);
}


/**
 * Creates a new bitmask expression from the given bitmask term.
 * <p>The expression returned is a valid in the sense that the epr_parse_bm_expr()
 * applied to the returned string would return an equivalent term.
 *
 * @param term the term to be converted
 */
char* epr_create_bm_expr(EPR_SBmTerm* term)
{
    if (term == NULL)
        return NULL;

    switch (term->op_code) {
    case BMT_REF: {
            char* s0 = epr_create_string(strlen(term->op.ref.band_name) + strlen(term->op.ref.flag_name) + 16);
            sprintf(s0, "%s.%s", term->op.ref.band_name, term->op.ref.flag_name);
            return s0;
        }
    case BMT_AND: {
            char* s1 = epr_create_bm_expr(term->op.binary.arg1);
            char* s2 = epr_create_bm_expr(term->op.binary.arg2);
            char* s0 = epr_create_string(strlen(s1) + strlen(s2) + 16);
            sprintf(s0, "(%s) AND (%s)", s1, s2);
            epr_free_string(s1);
            epr_free_string(s2);
            return s0;
        }
    case BMT_OR: {
            char* s1 = epr_create_bm_expr(term->op.binary.arg1);
            char* s2 = epr_create_bm_expr(term->op.binary.arg2);
            char* s0 = epr_create_string(strlen(s1) + strlen(s2) + 16);
            sprintf(s0, "(%s) OR (%s)", s1, s2);
            epr_free_string(s1);
            epr_free_string(s2);
            return s0;
        }
    case BMT_NOT: {
            char* s1 = epr_create_bm_expr(term->op.unary.arg);
            char* s0 = epr_create_string(strlen(s1) + 16);
            sprintf(s0, "NOT (%s)", s1);
            epr_free_string(s1);
            return s0;
        }
    default:
        assert(FALSE);
        return NULL;
    }
}

/**
 * Prints the given term as an expression to the console.
 */
void epr_print_bm_term(EPR_SBmTerm* term)
{
    epr_write_bm_term(term, stdout);
}

/**
 * Writes the given term as an expression to the given output stream.
 */
void epr_write_bm_term(EPR_SBmTerm* term, FILE* ostream)
{
    char* bm_expr = epr_create_bm_expr(term);
    fprintf(ostream, "%s", bm_expr);
    epr_free_string(bm_expr);
}
