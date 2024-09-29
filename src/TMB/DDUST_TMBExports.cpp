// Generated by TMBtools: do not edit by hand

#define TMB_LIB_INIT R_init_DDUST_TMBExports
#include <TMB.hpp>
#include "DDUST.hpp"

template<class Type>
Type objective_function<Type>::operator() () {
  DATA_STRING(model);
  if(model == "DDUST") {
    return DDUST(this);
  } else {
    Rf_error("Unknown model.");
  }
  return 0;
}
