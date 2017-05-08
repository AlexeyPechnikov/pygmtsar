#ifndef EPR_DDDB_H_INCL
#define EPR_DDDB_H_INCL

#include "epr_api.h"

#ifdef __cplusplus
extern "C"
{
#endif


struct RecordDescriptor {
    const char* id;
    const EPR_EDataTypeId type;
    const char* unit;
    const int elem_size;
    const char* num_elem;
    const char* description;
};

struct DatasetDescriptor {
    const char* id;
    const char* ds_name;
    const struct RecordDescriptor* rec_descriptor;
    const char* description;
};

struct BandDescriptor {
    const char* id;
    const char* rec_name;
    const EPR_ESampleModel sample_offset;
    const EPR_EDataTypeId type;
    const int spectral_index;
    const EPR_EScalingMethod scale_method;
    const char* scale_offset;
    const char* scale_factor;
    const char* bitmask_expr;
    const char* flag_coding_name;
    const char* unit;
    const char* description;
};

struct FlagDescriptor {
    const char* id;
    const int num_indices;
    const int bit_indices[2];
    const char* description;
};

struct DatasetDescriptorTable {
    const char* name;
    const char* description;
    int num_descriptors;
    const struct DatasetDescriptor* descriptors;
};

struct BandDescriptorTable {
    const char* name;
    const char* description;
    int num_descriptors;
    const struct BandDescriptor* descriptors;
};

struct FlagDescriptorTable {
    const char* name;
    const char* description;
    int num_descriptors;
    const struct FlagDescriptor* descriptors;
};

struct RecordDescriptorTable {
    const char* name;
    const char* description;
    int num_descriptors;
    const struct RecordDescriptor* descriptors;
};

extern const struct DatasetDescriptorTable dddb_product_tables[46];
extern const struct BandDescriptorTable dddb_band_tables[37];
extern const struct FlagDescriptorTable dddb_flag_coding_tables[6];
extern const struct RecordDescriptorTable dddb_meris_rec_tables[23];
extern const struct RecordDescriptorTable dddb_aatsr_rec_tables[20];
extern const struct RecordDescriptorTable dddb_asar_rec_tables[20];

#define EPR_NUM_PRODUCT_TABLES         46
#define EPR_NUM_BAND_TABLES            37
#define EPR_NUM_FLAG_CODING_TABLES     6
#define EPR_NUM_MERIS_REC_TABLES       23
#define EPR_NUM_AATSR_REC_TABLES       20
#define EPR_NUM_ASAR_REC_TABLES        20



#ifdef __cplusplus
} /* extern "C" */
#endif
#endif /* #ifndef EPR_DDDB_H_INCL */
