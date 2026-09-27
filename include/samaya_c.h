/* Stable C API for the samaya solver. All functions are safe to call from C and never throw. */
#ifndef SAMAYA_C_H
#define SAMAYA_C_H

#include <stddef.h>

#ifdef __cplusplus
extern "C" {
#endif

typedef struct samaya_model samaya_model;

/* Status codes, identical in value to samaya::Status. */
enum {
  SAMAYA_NOT_SOLVED = 0,
  SAMAYA_OPTIMAL = 1,
  SAMAYA_INFEASIBLE = 2,
  SAMAYA_UNBOUNDED = 3,
  SAMAYA_INFEASIBLE_OR_UNBOUNDED = 4,
  SAMAYA_TIME_LIMIT = 5,
  SAMAYA_ITERATION_LIMIT = 6,
  SAMAYA_NODE_LIMIT = 7,
  SAMAYA_NUMERICAL_ERROR = 8,
  SAMAYA_INVALID_MODEL = 9,
  SAMAYA_NOT_IMPLEMENTED = 10,
  SAMAYA_NOT_CONVEX = 11
};

/* Returns 0 on success. On failure *out is NULL and, if errbuf is non-NULL, a NUL-terminated
   message of at most errlen bytes is written to it. */
int samaya_read_mps(const char* path, samaya_model** out, char* errbuf, size_t errlen);
void samaya_model_free(samaya_model* model);

int samaya_model_num_rows(const samaya_model* model);
int samaya_model_num_cols(const samaya_model* model);
long long samaya_model_num_nonzeros(const samaya_model* model);

/* Solves with default parameters. Returns a status code; *objective receives the objective value
   (NaN if no solution) when objective is non-NULL. */
int samaya_solve(const samaya_model* model, double* objective);

const char* samaya_status_string(int status);
const char* samaya_version(void);

#ifdef __cplusplus
}
#endif

#endif /* SAMAYA_C_H */
