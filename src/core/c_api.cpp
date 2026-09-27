#include "samaya_c.h"

#include <cstring>
#include <exception>
#include <limits>

#include "samaya/io.hpp"
#include "samaya/solver.hpp"

struct samaya_model {
  samaya::Model model;
};

static_assert(SAMAYA_NOT_IMPLEMENTED == static_cast<int>(samaya::Status::kNotImplemented),
              "C status codes must mirror samaya::Status");
static_assert(SAMAYA_NOT_CONVEX == static_cast<int>(samaya::Status::kNotConvex),
              "C status codes must mirror samaya::Status");

namespace {

void write_error(char* errbuf, size_t errlen, const char* message) {
  if (errbuf == nullptr || errlen == 0) return;
  std::strncpy(errbuf, message, errlen - 1);
  errbuf[errlen - 1] = '\0';
}

}  // namespace

extern "C" {

int samaya_read_mps(const char* path, samaya_model** out, char* errbuf, size_t errlen) {
  if (out == nullptr) return 1;
  *out = nullptr;
  if (path == nullptr) {
    write_error(errbuf, errlen, "path is NULL");
    return 1;
  }
  try {
    *out = new samaya_model{samaya::read_mps(path)};
    return 0;
  } catch (const std::exception& e) {
    write_error(errbuf, errlen, e.what());
  } catch (...) {
    write_error(errbuf, errlen, "unknown error");
  }
  return 1;
}

void samaya_model_free(samaya_model* model) { delete model; }

int samaya_model_num_rows(const samaya_model* model) { return model ? model->model.num_rows() : 0; }

int samaya_model_num_cols(const samaya_model* model) { return model ? model->model.num_cols() : 0; }

long long samaya_model_num_nonzeros(const samaya_model* model) {
  return model ? static_cast<long long>(model->model.A.nnz()) : 0;
}

int samaya_solve(const samaya_model* model, double* objective) {
  if (objective != nullptr) *objective = std::numeric_limits<double>::quiet_NaN();
  if (model == nullptr) return SAMAYA_INVALID_MODEL;
  try {
    samaya::Params params;
    params.log_level = 0;
    const samaya::Result result = samaya::Solver(params).solve(model->model);
    if (objective != nullptr) *objective = result.objective;
    return static_cast<int>(result.status);
  } catch (...) {
    return SAMAYA_NUMERICAL_ERROR;
  }
}

const char* samaya_status_string(int status) {
  if (status < 0 || status > SAMAYA_NOT_CONVEX) return "unknown";
  return samaya::to_string(static_cast<samaya::Status>(status));
}

const char* samaya_version(void) { return samaya::version(); }

}  // extern "C"
