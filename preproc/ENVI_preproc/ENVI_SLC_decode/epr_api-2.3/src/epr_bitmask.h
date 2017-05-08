/*
 * $Id: epr_bitmask.h,v 1.1.1.1 2004-10-28 19:22:22 norman Exp $
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

#ifndef EPR_BITMASK_H_INCL
#define EPR_BITMASK_H_INCL

#ifdef __cplusplus
extern "C"
{
#endif


#define FLAG_MASK_NOT_COMPUTED ((uint) -1)

typedef struct EPR_BmTerm               EPR_SBmTerm;
typedef struct EPR_BmEvalContext        EPR_SBmEvalContext;
typedef struct EPR_BmFlagDataset        EPR_SBmFlagDataset;
typedef enum   EPR_BmOpCode             EPR_EBmOpCode;


/* private implementations */


enum EPR_BmOpCode {
    BMT_UNKNOWN = 0,
    BMT_REF,
    BMT_AND,
    BMT_OR,
    BMT_NOT
};


enum EPR_Tok {
    BME_UNKNOWN = 0,
    BME_EOS,
    BME_SPECIAL,
    BME_NAME
};

/**
 * The <code>EPR_BmTerm</code> structure is the union of structures:
 * each of them can contain either the subject (operand) for the logic operators
 * or this operators itself with referd operand(s). Thus they are recursive.
 * The example of term: <code>flags.WATER or flags.LAND</code>
 * here: <i>'flags'</i> is a band_name; <i>'WATER'</i> and <i>'LAND'</i> - flag_name's; <i>'or'</i> - logical operator.
 *
 * @see EPR_SRaster
 */
struct EPR_BmTerm {
    EPR_EBmOpCode op_code;
    union {
        struct /*BMT_REF*/ {
            char* band_name;
            char* flag_name;
            uint flag_mask;
            EPR_SRaster* flag_raster;
        } ref;
        struct /*BMT_NOT*/ {
            EPR_SBmTerm* arg;
        } unary;
        struct /*BMT_AND and BMT_OR*/ {
            EPR_SBmTerm* arg1;
            EPR_SBmTerm* arg2;
        } binary;
    } op;
};


/**
 * The <code>EPR_BmEvalContext</code> structure represents an evaluation context for bitmask
 * expressions. It is used internally only.
 * <p>
 * An instance of this structure holds the product ID, references to all flag datasets
 * required to evaluate the expression as well as information about the bitmask raster beeing
 * created.
 *
 * @see EPR_SProductId
 * @see EPR_SRaster
 * @see EPR_SPtrArray
 */
struct EPR_BmEvalContext
{
    /**
     * The ID of the product to which this band belongs to.
     */
    EPR_SProductId* product_id;

    /**
     * X-coordinate in pixel co-ordinates (zero-based) of the upper right corner raster.
     */
    int offset_x;

     /**
     * Y-coordinate in pixel co-ordinates (zero-based) of the upper right corner raster.
     */
    int offset_y;

    /**
     * The result 0/1 bit mask raster for all flag bands
     */
    EPR_SRaster* bitmask_raster;

    /**
     * The band_id of flags (can be 0, 1 or 4 - for AATSR)
     */
    EPR_SPtrArray* flag_band_ids;

    /**
     * The corresponding 0/1 bit mask raster for each flag bands
     */
    EPR_SPtrArray* flag_rasters;
};


/**
 * Represents a flag-field within a flag-record.
 *
 */
struct EPR_BmFlagDataset {
    /*The name of bitmask dataset*/
    char* name;
    /*The value of bitmask dataset (the number of the relevant bit in bitmask)*/
    uint bit_index;
    /*The description of bitmask dataset*/
    char* description;
};


/**
 * Creates bit-mask evaluation context
 * for the given raster and offsets of start corner
 *
 * @param product_id the product ID
 * @param offset_x X-coordinate in pixel co-ordinates (zero-based) of the upper right corner raster to search
 * @param offset_y Y-coordinate in pixel co-ordinates (zero-based) of the upper right corner raster to search
 * @param raster the bitmask_raster
 *
 * @return bit-mask evaluated context for success, and error code otherwise
 */
EPR_SBmEvalContext* epr_create_bm_eval_context(EPR_SProductId* product_id,
                                               int offset_x,
                                               int offset_y,
                                               EPR_SRaster* raster);


/**
 * Release the memory allocated through a EPR_SBmEvalContext.
 *
 * @param context the bit mask context, if <code>NULL</code> the function
 *        immediately returns zero.
 * @return zero for success, an error code otherwise
 */
void epr_free_bm_eval_context(EPR_SBmEvalContext* context);


/**
 * Reads bit-mask pixels of the given product for the given bit-mask expression
 * for the given region and and with the given sub-sampling.
 *
 * <blockquote>
 * <p><i>bit-mask-expression :=</i><br>
 * &nbsp;&nbsp;&nbsp;&nbsp;<i>or-expression</i>
 *
 * <p><i>or-expression :=</i><br>
 * &nbsp;&nbsp;&nbsp;&nbsp;<i>and-expression</i><br>
 * &nbsp;&nbsp;&nbsp;&nbsp;<i>or-expression</i> <b><code>or</code></b> <i>and-expression</i>
 *
 * <p><i>and-expression :=</i><br>
 * &nbsp;&nbsp;&nbsp;&nbsp;<i>not-expression</i><br>
 * &nbsp;&nbsp;&nbsp;&nbsp;<i>and-expression</i> <b><code>and</code></b> <i>not-expression</i>
 *
 * <p><i>not-expression :=</i><br>
 * &nbsp;&nbsp;&nbsp;&nbsp;<i>primary-expression</i><br>
 * &nbsp;&nbsp;&nbsp;&nbsp;<b><code>not</code></b> <i>not-expression</i>
 *
 * <p><i>primary-expression :=</i><br>
 * &nbsp;&nbsp;&nbsp;&nbsp;<i>flag-reference</i><br>
 * &nbsp;&nbsp;&nbsp;&nbsp;<b><code>(</code></b> <i>bit-mask-expression</i> <b><code>)</code></b>
 *
 * <p><i>flag-reference :=</i><br>
 * &nbsp;&nbsp;&nbsp;&nbsp;<i>dataset-name</i><b><code>.</code></b><i>flag-name</i></b>
 * </blockquote>
 *
 * <p>Where <i>dataset-name</i> and <i>flag-name</i> are names specific for a particular data product.
 * Names are in general resolved case-insenitively. The parser also accepts an alternate notation for
 * the boolean operators:
 * <blockquote>
 * The <b><code>|</code></b> character for the <b><code>or</code></b> operator,<br>
 * the <b><code>&</code></b> character for the <b><code>and</code></b> operator and finally<br>
 * the <b><code>!</code></b> character for the <b><code>not</code></b> operator.
 * </blockquote>
 *
 * @param product_id the product ID
 * @param bm_expr the bit-mask expression
 * @param xo X-coordinate in pixel co-ordinates (zero-based) of the upper right corner raster to search
 * @param yo Y-coordinate in pixel co-ordinates (zero-based) of the upper right corner raster to search
 * @param raster_width  the width in pixel co-ordinates of the raster to search
 * @param raster_height the height in pixel co-ordinates of raster to search
 * @param s_x X-step in pixel co-ordinates to get the next raster to search
 * @param s_y Y-step in pixel co-ordinates to get the next raster to search
 * @param raster_buffer [BYTE] the memory buffer to save information was read
 *
 * @return zero for success, and error code otherwise
 */
int epr_read_bitmask_data(const EPR_SProductId* product_id,
                          const char* bm_expr,
                          int xo,
                          int yo,
                          int raster_width,
                          int raster_height,
                          int s_x,
                          int s_y,
                          void* raster_buffer);



/**
 * Evaluates the given bitmask expression.
 *
 * @param term the bitmask term
 * @param x the x co-ordinate in pixels
 * @param y the y co-ordinate in pixels
 */
epr_boolean epr_eval_bm_term(EPR_SBmEvalContext* context,
                         EPR_SBmTerm* term,
                         int x,
                         int y);



/**
 * Parses a bitmask expression string.
 *
 * <p>The bit-mask expressions recognized by this parser must have the following syntax:
 *
 * <blockquote>
 * <p><i>bit-mask-expression :=</i><br>
 * &nbsp;&nbsp;&nbsp;&nbsp;<i>or-expression</i>
 *
 * <p><i>or-expression :=</i><br>
 * &nbsp;&nbsp;&nbsp;&nbsp;<i>and-expression</i><br>
 * &nbsp;&nbsp;&nbsp;&nbsp;<i>or-expression</i> <b><code>or</code></b> <i>and-expression</i>
 *
 * <p><i>and-expression :=</i><br>
 * &nbsp;&nbsp;&nbsp;&nbsp;<i>not-expression</i><br>
 * &nbsp;&nbsp;&nbsp;&nbsp;<i>and-expression</i> <b><code>and</code></b> <i>not-expression</i>
 *
 * <p><i>not-expression :=</i><br>
 * &nbsp;&nbsp;&nbsp;&nbsp;<i>primary-expression</i><br>
 * &nbsp;&nbsp;&nbsp;&nbsp;<b><code>not</code></b> <i>not-expression</i>
 *
 * <p><i>primary-expression :=</i><br>
 * &nbsp;&nbsp;&nbsp;&nbsp;<i>flag-reference</i><br>
 * &nbsp;&nbsp;&nbsp;&nbsp;<b><code>(</code></b> <i>bit-mask-expression</i> <b><code>)</code></b>
 *
 * <p><i>flag-reference :=</i><br>
 * &nbsp;&nbsp;&nbsp;&nbsp;<i>dataset-name</i><b><code>.</code></b><i>flag-name</i></b>
 * </blockquote>
 *
 * <p>Where <i>dataset-name</i> and <i>flag-name</i> are names specific for a particular data product.
 * Names are in general resolved case-insenitively. The parser also accepts an alternate notation for
 * the boolean operators:
 * <blockquote>
 * The <b><code>|</code></b> character for the <b><code>or</code></b> operator,<br>
 * the <b><code>&</code></b> character for the <b><code>and</code></b> operator and finally<br>
 * the <b><code>!</code></b> character for the <b><code>not</code></b> operator.
 * </blockquote>
 *
 * <p>For example, the following parseBitmaskExpression request will perform without errors:
 * <pre>
 *     BitmaskTerm term = BitmaskExpressionParser.parse("flags.LAND and not flags.DDV");
 * </pre>
 * <p>Another example for a valid expression in alternatate notation is:
 * <pre>
 *     BitmaskTerm term = BitmaskExpressionParser.parse("flags.LAND | (flags.COASTLINE & !flags.CLOUD)");
 * </pre>
 *
 * <p>The terms created in the examples above could successfully be evaluated in an evaluation context
 * provided by an ENVISAT MERIS Level 2 data product.
 *
 * @param bm_expr the bitmask expression
 * @return the bitmask term representing the given expression
 */
EPR_SBmTerm* epr_parse_bm_expr_str(const char* bm_expr);


struct EPR_ParseInfo {
    const char* bm_expr;
    int bm_expr_pos;
    epr_boolean pushed_back;
    int token_type;
    char* token;
    char* err_message;
};


typedef struct EPR_ParseInfo EPR_SParseInfo;

/**
 * This group of functions is for parsing the expression.
 *
 * @param parse_info parse_info structure
 * @param term_required the boolean value expression.
 *
 * @return the bit mask term (see EPR_BmTerm).
 */
/*@{*/
EPR_SBmTerm* epr_parse_bm_expr(EPR_SParseInfo* parse_info, epr_boolean term_required);
EPR_SBmTerm* epr_parse_bm_OR_expr(EPR_SParseInfo* parse_info, epr_boolean term_required);
EPR_SBmTerm* epr_parse_bm_AND_expr(EPR_SParseInfo* parse_info, epr_boolean term_required);
EPR_SBmTerm* epr_parse_bm_unary_expr(EPR_SParseInfo* parse_info, epr_boolean term_required);
EPR_SBmTerm* epr_parse_bm_primary_expr(EPR_SParseInfo* parse_info, epr_boolean term_required);
/*@}*/

/**
 * This group of functions is for recognizing the keyword.
 *
 * @param parse_info parse_info structure
 *
 * @return TRUE or FALSE.
 */
/*@{*/
epr_boolean epr_is_bm_OR_keyword(EPR_SParseInfo* parse_info);
epr_boolean epr_is_bm_AND_keyword(EPR_SParseInfo* parse_info);
epr_boolean epr_is_bm_NOT_keyword(EPR_SParseInfo* parse_info);
/*@}*/

/**
 * This group of functions is for recognizing the operator.
 *
 * @param parse_info parse_info structure
 *
 * @return TRUE or FALSE.
 */
/*@{*/
epr_boolean epr_is_bm_AND_operator(EPR_SParseInfo* parse_info);
epr_boolean epr_is_bm_OR_operator(EPR_SParseInfo* parse_info);
epr_boolean epr_is_bm_NOT_operator(EPR_SParseInfo* parse_info);
/*@}*/

/**
 * Tests the given expression for operand name only (not operator).
 *
 * @param parse_info parse_info structure
 *
 * @return TRUE if the term is not NULL and an operand, or FALSE otherwise
 */
epr_boolean epr_is_bm_name_token(EPR_SParseInfo* parse_info);

/**
 * Tests the given expression for the end of string.
 *
 * @param parse_info parse_info structure
 *
 * @return TRUE if the EOS occurs, or FALSE otherwise
 */
epr_boolean epr_is_bm_EOS_token(EPR_SParseInfo* parse_info);

/**
 * Tests the given expression for errors.
 *
 * @param parse_info parse_info structure
 *
 * @return TRUE if no error occurs, or FALSE otherwise
 */
epr_boolean epr_is_bm_expr_error(EPR_SParseInfo* parse_info);

/**
 * Gets the first character of token for the given expression or EOS.
 *
 * @param parse_info parse_info structure
 *
 * @return '(' , ')', '.' , '&' , '|' ,'!', or '\0' otherwise
 */
int epr_get_token_char(EPR_SParseInfo* parse_info);

/**
 * Releases the actual token given expression.
 *
 * @param parse_info parse_info structure
 *
 * @return token
 */
char* epr_consume_token(EPR_SParseInfo* parse_info);

/**
 * Selectss the next token given expression.
 *
 * @param parse_info parse_info structure
 */
void epr_next_bm_expr_token(EPR_SParseInfo* parse_info);

void epr_push_back_bm_expr_token(EPR_SParseInfo* parse_info);

void epr_set_bm_expr_error(EPR_SParseInfo* parse_info, const char* message);

int epr_tokenize_bm_expr(const char* bm_expr, int* bm_expr_pos, char** token);

/**
 * Creates a new bitmask term instance.
 */
EPR_SBmTerm* epr_create_bm_term(EPR_EBmOpCode op_code);

/**
 * Creates a new bitmask reference term instance.
 */
EPR_SBmTerm* epr_create_bm_REF_term(char* ds_name, char* flag_name);

/**
 * Creates a new bitmask NOT term instance.
 */
EPR_SBmTerm* epr_create_bm_NOT_term(EPR_SBmTerm* arg);

/**
 * Creates a new bitmask OR term instance.
 */
EPR_SBmTerm* epr_create_bm_OR_term(EPR_SBmTerm* arg1, EPR_SBmTerm* arg2);

/**
 * Creates a new bitmask reference AND term instance.
 */
EPR_SBmTerm* epr_create_bm_AND_term(EPR_SBmTerm* arg1, EPR_SBmTerm* arg2);

/**
 * Releases a new bitmask term instance.
 */
void epr_free_bm_term(EPR_SBmTerm* term);


/**
 * Creates a new bitmask expression from the given bitmask term.
 * <p>The expression returned is a valid in the sense that the epr_parse_bm_expr()
 * applied to the returned string would return an equivalent term.
 *
 * @param term the term to be converted
 */
char* epr_create_bm_expr(EPR_SBmTerm* term);
/**
 * Prints the given term as an expression to the console.
 */
void epr_print_bm_term(EPR_SBmTerm* term);

/**
 * Writes the given term as an expression to the given output stream.
 */
void epr_write_bm_term(EPR_SBmTerm* term, FILE* ostream);

/**
 * Creates the coding flag info
 *
 * @param str the local path to dddb
 *
 * @return the the pointer at the coding flag information.
 */
EPR_SPtrArray* epr_create_flag_coding(EPR_SProductId* product_id, const char* str);
/**
 * Creates the coding flag definition
 *
 * @return the the pointer at the coding flag definition information.
 */
EPR_SFlagDef* epr_create_flag_def();
/**
 * Releases the coding flag definition
 */
void epr_free_flag_def(EPR_SFlagDef* flag_def);

/**
 * Releases the coding flag info
 */
void epr_free_flag_coding(EPR_SPtrArray* flag_coding);

#ifdef __cplusplus
} /* extern "C" */
#endif

#endif
/* #ifndef EPR_BITMASK_H_INCL */
