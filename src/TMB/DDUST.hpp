// Copyright 2024 Fisheries Queensland

// This file is part of DDUST.
// DDUST is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
// DDUST is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
// You should have received a copy of the GNU General Public License along with DDUST. If not, see <https://www.gnu.org/licenses/>.

#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR obj

template<class Type>
Type BevertonHolt(Type h, Type R0, Type S0, Type S) {
  Type alpha = (S0*(1-h))/(4*h*R0);
  Type beta = (5*h-1)/(4*h*R0);
  Type recruitment = S/(alpha+beta*S);
  return recruitment;
}

template<class Type>
struct vonMises {
  Type kappa;
  Type mu;
  vonMises (Type kappa_, Type mu_)
  : kappa (kappa_), mu (mu_) {}
  Type operator()(Type x){
    return exp(kappa*cos(x-mu))/(2*M_PI*besselI(kappa,Type(0)));
  }
};

template<class Type>
struct economic_list {
  vector<Type> cL;
  vector<Type> cM;
  vector<Type> cK;
  vector<Type> cF;
  vector<Type> cO;
  vector<Type> v;
  vector<Type> W;
  vector<Type> K;
  vector<Type> prop;
  vector<Type> B;
  vector<Type> dmean;
  vector<Type> o;
  vector<Type> d;
  economic_list(SEXP x){ // Constructor
    cL = asVector<Type>(getListElement(x,"cL"));
    cM = asVector<Type>(getListElement(x,"cM"));
    cK = asVector<Type>(getListElement(x,"cK"));
    cF = asVector<Type>(getListElement(x,"cF"));
    cO = asVector<Type>(getListElement(x,"cO"));
    v = asVector<Type>(getListElement(x,"v"));
    W = asVector<Type>(getListElement(x,"W"));
    K = asVector<Type>(getListElement(x,"K"));
    prop = asVector<Type>(getListElement(x,"prop"));
    B = asVector<Type>(getListElement(x,"B"));
    dmean = asVector<Type>(getListElement(x,"dmean"));
    o = asVector<Type>(getListElement(x,"o"));
    d = asVector<Type>(getListElement(x,"d"));
  }
};

template<class Type>
struct PopulationVectors {
  vector<Type> AnnualBiomass;
  vector<Type> AnnualF;
  vector<Type> AnnualRecruits;
  vector<Type> AnnualSpawners;
  vector<Type> AnnualSpawnersRatio;
  vector<Type> AnnualYield;
  vector<Type> Biomass;
  vector<Type> FishingMortality;
  vector<Type> MidBiomass;
  vector<Type> Numbers;
  vector<Type> Survivors;
  vector<Type> Spawners;
  vector<Type> Recruits;
  vector<Type> BevertonHoltRecruits;
  vector<Type> RecruitmentDeviation;
  vector<Type> Yield;
  vector<Type> Effort;
  vector<Type> HarvestRate;
};

template<class Type>
PopulationVectors<Type> PopulationDynamics(int num_years, vector<Type> B_initial, vector<Type> N_initial, vector<Type> R_initial, Type Sp_initial, Type sp_total,
                                           Type rho, vector<Type> weight_at_recruitment, Type mid_surv, vector<Type> sp_pattern, int dt, Type h, vector<Type> rec_pattern,
                                           Type M, Type F_mort, Type R0, array<Type> cpue, vector<Type> ctch, bool use_catch_data, bool equilibrium, vector<Type> eta,
                                           int recruitment_deviation_type, vector<Type> log_R_star, int first_year_rec_devs, int first_year_catch,
                                           vector<Type> F_pattern, vector<Type> q_pattern) {

  // Derived quantities
  int T = num_years*dt;
  Type surv = exp(-M/dt);

  // Vectors for every timestep
  vector<Type> BaseRecruitment(T+2);
  vector<Type> Biomass(T+2);
  vector<Type> FishingMortality(T+2);
  vector<Type> mid_biomass(T);
  vector<Type> Numbers(T+2);
  vector<Type> Survivors(T+2);
  vector<Type> Recruitment(T+2);
  vector<Type> RecruitmentDeviation(T/dt);
  vector<Type> Spawners(T);
  vector<Type> SS(dt);
  vector<Type> U(T+2);
  vector<Type> Yield(T+2);
  vector<Type> Effort(T+2);

  // Annual vectors
  vector<Type> AnnualRecruitment(T/dt);
  vector<Type> AnnualSpawners(T/dt+2);
  vector<Type> AnnualYield(num_years);
  vector<Type> AnnualBiomass(num_years);
  vector<Type> AnnualFishingMortality(num_years);
  vector<Type> AnnualSpawnersRatio(num_years);

  // Fill vectors to prevent accidental uninitialised memory access
  AnnualBiomass.fill(0.0);
  AnnualFishingMortality.fill(0.0);
  AnnualYield.fill(0.0);
  AnnualSpawnersRatio.fill(0.0);
  FishingMortality.fill(0.0);
  Yield.fill(0.0);
  Effort.fill(0.0);
  Biomass.fill(B_initial(0));
  Numbers.fill(N_initial(0));
  Survivors.fill(surv);
  Recruitment.fill(R_initial(0));
  mid_biomass.fill(B_initial(0));
  Spawners.fill(sp_total);
  SS.fill(Sp_initial/dt);
  U.fill(0);
  AnnualSpawners.fill(sp_total);
  AnnualRecruitment.fill(R_initial(0));
  RecruitmentDeviation.fill(0);
  BaseRecruitment.fill(R_initial(0));

  // T-1 and T-2
  Biomass(0) = B_initial(10);
  Biomass(1) = B_initial(11);
  Numbers(0) = N_initial(10);
  Numbers(1) = N_initial(11);
  Survivors(0) = surv;
  Survivors(1) = surv;
  Recruitment(0) = R_initial(10);
  Recruitment(1) = R_initial(11);

  int year = 0;
  for (int t=2;t<(T+2);t++){

    if ((t-2) % dt == 0){ // Calculate recruitment in the first month of the year
      AnnualSpawners(year) = SS.sum();

      if (equilibrium){ // Equilibrium uses R0 every year
        AnnualRecruitment(year) = R0;
      } else {
        AnnualRecruitment(year) = BevertonHolt(h, R0, sp_total, AnnualSpawners(year));
      }

      for (int month=0;month<dt;month++){
        BaseRecruitment(t+month) = AnnualRecruitment(year)*rec_pattern(month);
        if (recruitment_deviation_type == 0 || year < first_year_rec_devs - first_year_catch){
          // No recruitment deviation
          Recruitment(t+month) = BaseRecruitment(t+month);
          RecruitmentDeviation(year) = 0;
        } else if (recruitment_deviation_type == 1) {
          // Barycentric recruitment deviations
          Recruitment(t+month) = BaseRecruitment(t+month)*eta(year);
          RecruitmentDeviation(year) = log(eta(year));
        } else if (recruitment_deviation_type == 2) {
          // Redundant recruitment deviations
          Recruitment(t+month) = exp(log_R_star(year - (first_year_rec_devs - first_year_catch)))*rec_pattern(month);
          RecruitmentDeviation(year) = log_R_star(year - (first_year_rec_devs - first_year_catch)) - log(BevertonHolt(h, R0, sp_total, AnnualSpawners(year)));
        }
      }
      year += 1;
    }

    Biomass(t) = (1+rho)*Survivors(t-1)*Biomass(t-1) - rho*Survivors(t-1)*Survivors(t-2)*Biomass(t-2) - rho*Survivors(t-1)*weight_at_recruitment(0)*Recruitment(t-1) + weight_at_recruitment(1)*Recruitment(t);
    Numbers(t) = Survivors(t-1)*Numbers(t-1) + Recruitment(t);
    mid_biomass(t-2) = Biomass(t) * mid_surv;
    int month_index = (t-2) % dt;

    // Harvest and fishing mortality
    if (use_catch_data){
      Type catch_proportion = ctch(t-2)/mid_biomass(t-2);
      Type max_catch = 0.99;
      U(t) = std::min(catch_proportion,max_catch); // Minimum function might store the first evaluation and not change when optimisation happens for example
      FishingMortality(t) = -log(1-U(t));
      Yield(t) = U(t)*mid_biomass(t-2);
      if (cpue(0,t-2)>0) {
        Effort(t) = cpue(0,t-2)/ctch(t-2);
      }
    } else {
      U(t) = 1-exp(-F_mort*F_pattern(month_index));
      FishingMortality(t) = F_mort*F_pattern(month_index);
      Yield(t) = U(t)*mid_biomass(t-2);
      Type q = q_pattern(month_index);
      Effort(t) = U(t)/q;
    }

    // Spawners
    Survivors(t) = surv*(1-U(t));
    SS(month_index) = sp_pattern(month_index)*((1-Survivors(t))/(-log(Survivors(t))))*Numbers(t)*0.5;
    Spawners(t-2) = SS(month_index);

    // Annual summaries
    AnnualBiomass(year-1) += Biomass(t)/dt;
    AnnualFishingMortality(year-1) += FishingMortality(t);
    AnnualYield(year-1) += Yield(t);
    AnnualSpawnersRatio(year-1) = AnnualSpawners(year-1)/sp_total;
  }

  return PopulationVectors<Type> {
    .AnnualBiomass = AnnualBiomass,
    .AnnualF = AnnualFishingMortality,
    .AnnualRecruits = AnnualRecruitment,
    .AnnualSpawners = AnnualSpawners,
    .AnnualSpawnersRatio = AnnualSpawnersRatio,
    .AnnualYield = AnnualYield,
    .Biomass = Biomass,
    .FishingMortality = FishingMortality,
    .MidBiomass = mid_biomass,
    .Numbers = Numbers,
    .Survivors = Survivors,
    .Spawners = Spawners,
    .Recruits = Recruitment,
    .BevertonHoltRecruits = BaseRecruitment,
    .RecruitmentDeviation = RecruitmentDeviation,
    .Yield = Yield,
    .Effort = Effort,
    .HarvestRate = U,
  };
}

template<class Type>
struct Functor {
  vector<Type> F_pattern;
  Type mid_surv;
  int projection_years;
  int dt;
  vector<Type> R_initial;
  vector<Type> B_initial;
  vector<Type> N_initial;
  Type surv;
  Type sp_total;
  Type R0;
  vector<Type> rec_pattern;
  Type rho;
  Type h;
  vector<Type> weight_at_recruitment;
  vector<Type> sp_pattern;
  vector<Type> q_pattern;
  Type target_relative_biomass;
  int report_type;
  vector<Type> cL;
  vector<Type> cM;
  vector<Type> cK;
  vector<Type> cF;
  vector<Type> cO;
  vector<Type> v;
  vector<Type> W;
  vector<Type> K;
  vector<Type> prop;
  vector<Type> B;
  vector<Type> dmean;
  vector<Type> o;
  vector<Type> d;
  Functor (const vector<Type> &F_pattern,
           const Type &mid_surv,
           const int &projection_years,
           const int &dt,
           const vector<Type> &R_initial,
           const vector<Type> &B_initial,
           const vector<Type> &N_initial,
           const Type &surv,
           const Type &sp_total,
           const Type &R0,
           const vector<Type> &rec_pattern,
           const Type &rho,
           const Type &h,
           const vector<Type> &weight_at_recruitment,
           const vector<Type> &sp_pattern,
           const vector<Type> &q_pattern,
           const Type &target_relative_biomass,
           const int &report_type,
           const vector<Type> &cL,
           const vector<Type> &cM,
           const vector<Type> &cK,
           const vector<Type> &cF,
           const vector<Type> &cO,
           const vector<Type> &v,
           const vector<Type> &W,
           const vector<Type> &K,
           const vector<Type> &prop,
           const vector<Type> &B,
           const vector<Type> &dmean,
           const vector<Type> &o,
           const vector<Type> &d) :
    F_pattern(F_pattern),
    mid_surv(mid_surv),
    projection_years(projection_years),
    dt(dt),
    R_initial(R_initial),
    B_initial(B_initial),
    N_initial(N_initial),
    surv(surv),
    sp_total(sp_total),
    R0(R0),
    rec_pattern(rec_pattern),
    rho(rho),
    h(h),
    weight_at_recruitment(weight_at_recruitment),
    sp_pattern(sp_pattern),
    q_pattern(q_pattern),
    target_relative_biomass(target_relative_biomass),
    report_type(report_type),
    cL(cL),
    cM(cM),
    cK(cK),
    cF(cF),
    cO(cO),
    v(v),
    W(W),
    K(K),
    prop(prop),
    B(B),
    dmean(dmean),
    o(o),
    d(d) {}
  Type operator()(const vector<Type> &F_mort) {

    // Distribute fishing mortality and harvest rate throughout the year according to the fishing mortality pattern
    vector<Type> U_dt(dt);
    vector<Type> F_dt(dt);
    for (int i=0;i<dt;i++){
      F_dt(i) = F_mort(0)*F_pattern(i);
      U_dt(i) = 1-exp(-F_dt(i));
    }

    vector<Type> dummyvec(projection_years);
    dummyvec.fill(1.0);
    array<Type> dummyarray(projection_years);
    dummyarray.fill(1.0);
    int none = 0;
    Type M = -log(surv)*dt;

    PopulationVectors<Type> Population = PopulationDynamics(projection_years, B_initial, N_initial, R_initial, sp_total, sp_total, rho, weight_at_recruitment, mid_surv,
                                                       sp_pattern, dt, h, rec_pattern, M, F_mort(0), R0, dummyarray, dummyvec, false, false, dummyvec, none,
                                                       dummyvec, none, none, F_pattern, q_pattern);

    Type final_yield = Population.AnnualYield.sum(); //vector has length = projection_years, so -1 to access last index
    Type final_biomass_ratio = Population.AnnualBiomass(projection_years-1)/(sum(B_initial)/12);

    int T = projection_years*dt;
    vector<Type> Yield(T+2);
    vector<Type> Effort(T+2);
    Effort = Population.Effort;
    Yield = Population.Yield;

    vector<Type> CatchVariableCosts(T+2);
    vector<Type> EffortVariableCosts(T+2);
    CatchVariableCosts = (cL(0) * v(0) + cM(0))*Yield;
    EffortVariableCosts = (cK(0) + cF(0) + cO(0))*Effort;
    Type AnnualFixedCosts = W(0) + (o(0) + d(0)) * K(0);

    vector<Type> AnnualProfit(projection_years);
    AnnualProfit.fill(0.0);
    Type NPV = 0;

    for (int year=0;year<projection_years;year++){
      for (int m=0;m<dt;m++){
        AnnualProfit(year) += v(m)*Yield(year*dt+m+2) - CatchVariableCosts(year*dt+m+2) - EffortVariableCosts(year*dt+m+2) +
                              B(0)*Effort(year*dt+m+2) - AnnualFixedCosts/dt*Effort(year*dt+m+2)/dmean(0)*prop(0);
      }
      if (year < (projection_years-1)) {
        NPV += AnnualProfit(year)/pow(1+o(0),year);
      }
    }

    NPV += (AnnualProfit(projection_years-1)/o(0))/pow(1+o(0),projection_years-1);

    if (report_type == 0){ // return final biomass relative to initial biomass
      return final_biomass_ratio;
    } else if (report_type == 1){ // return negative yield for newton minimisation
      return -final_yield;
    } else if (report_type == 2){ // return discounted economic value
      return -NPV;
    } else if (report_type == 3){
      return pow(final_biomass_ratio - target_relative_biomass, 2); // return square difference for newton minimisation
    } else {
      return final_yield;
    }
  }
};

template <class Type>
Type DDUST(objective_function<Type>* obj) {

  using newton::newton_config_t;
  using newton::Newton;

  // Time-step Decisions
  DATA_INTEGER(Number_months_per_timestep); //12: annual model, 1: monthly model

  // Data
  DATA_VECTOR(proportion_spawning);
  DATA_VECTOR(weight_at_recruitment);
  DATA_SCALAR(weight_inf);
  DATA_INTEGER(last_year_catch);
  DATA_INTEGER(first_year_catch);
  DATA_INTEGER(first_year_rec_devs);
  DATA_VECTOR(ctch);
  DATA_ARRAY(cpue);
  DATA_ARRAY(cpue_sd);
  DATA_INTEGER(use_cpue_sd);
  DATA_VECTOR(absolute_biomass);
  DATA_VECTOR(absolute_biomass_sd);
  DATA_VECTOR(month_sequence);
  DATA_MATRIX(CoordBasis);
  DATA_INTEGER(num_years_warmup); // Number of years for equilibrium
  DATA_INTEGER(calculate_rho); // 1: Yes, 0: No
  DATA_SCALAR(rho_input);
  DATA_INTEGER(use_recruit_penalty);
  DATA_SCALAR(prior_mean_mu);
  DATA_SCALAR(prior_mean_k);
  DATA_SCALAR(prior_mean_xi);
  DATA_SCALAR(prior_sd_mu);
  DATA_SCALAR(prior_sd_k);
  DATA_SCALAR(prior_sd_xi);
  DATA_SCALAR(minimum_annual_harvest);
  DATA_SCALAR(recruit_penalty_strength);
  DATA_SCALAR(recruit_penalty_exponent);
  DATA_STRING(rec_dev_type); // "redundant", "off" or "barycentric"

  DATA_INTEGER(do_projections);
  DATA_SCALAR(target_relative_biomass);
  DATA_INTEGER(projection_years);
  DATA_VECTOR(projection_harvest);
  DATA_SCALAR(F_initial);
  DATA_STRUCT(cfg, newton_config_t);
  DATA_STRUCT(economic_data, economic_list);

  // Parameters
  PARAMETER(Rinit);
  PARAMETER(xi);
  PARAMETER(M);
  PARAMETER(k);
  PARAMETER(mu);
  PARAMETER(q1);
  PARAMETER(q2);
  PARAMETER(lsigmaR_sq);
  PARAMETER(lsigmaI_sq);
  PARAMETER_VECTOR(zeta);
  PARAMETER_VECTOR(log_R_star);

  // Derived quantities
  Type sigmaR_sq = exp(lsigmaR_sq); // recruitment variance
  Type R0 = exp(Rinit); // initial recruitment
  int dt = 12/Number_months_per_timestep; // number of timesteps per year
  Type surv = exp(-M/dt); // convert annual natural mortality to survivability in each timestep
  Type mid_surv = exp(-M/(2*dt));
  Type B0 = 0.5*R0; // initial guess for biomass
  Type N0 = 2*R0; // initial guess for numbers
  Type h = (1+exp(xi))/(5+exp(xi)); // steepness
  Type kappa = sqrt(pow(k,2)); // diffferentiable |k|

  // Flags for specifying rec dev type
  int none = 0;
  int barycentric = 1;
  int redundant = 2;

  // Growth parameter from weight-at-recruitment or user input
  Type rho = 0;
  if (calculate_rho == 1){
    rho = 1-(weight_at_recruitment(1)-weight_at_recruitment(0))/(weight_inf - weight_at_recruitment(0));
  } else {
    rho = rho_input;
  }

  // Define timescales
  int num_years = last_year_catch - first_year_catch + 1; // number of years
  int T = num_years*dt; // number of timesteps
  int T_warmup = num_years_warmup*dt; // number of warm-up timesteps

  // Set up barycentric recruitment deviations
  vector<Type> zeta_matrix = CoordBasis * zeta; // coordinate basis matrix to scale the distance of residuals (vertices of the simplex) from zero
  vector<Type> eta_temp = exp(zeta_matrix.array()); // exponentiate for multiplicative deviations
  vector<Type> eta(num_years);
  eta.fill(1.0);

  if (rec_dev_type == "barycentric"){
    // No recruitment deviations before first_year_rec_devs
    for (int t=0;t<(first_year_rec_devs-first_year_catch);t++){
      eta(t) = 1;
    }
    // Recruitment deviations begin
    for (int t=(first_year_rec_devs-first_year_catch);t<num_years;t++){
      eta(t) = eta_temp(t-(first_year_rec_devs-first_year_catch));
    }
  }

  // Seasonal recruitment pattern (by month)
  vonMises<Type> f(kappa, mu);
  vector<Type> phi_t(12);
  phi_t.fill(1/12);
  for (int t=0;t<12;t++){
    Type b1 = 2*M_PI/12*t - M_PI;
    Type b2 = 2*M_PI/12*(t+1) - M_PI;
    phi_t(t) = romberg::integrate(f, b1, b2);
  }

  // Aggregate recruitment pattern according to Number_months_per_timestep
  vector<Type> rec_pattern(T_warmup);
  rec_pattern.fill(0.0);
  for (int t=0;t<T_warmup;t++){
    rec_pattern(t) = 0;
    for (int m=t*Number_months_per_timestep;m<Number_months_per_timestep*(t+1);m++){
      int month_mod_12 = m % 12;
      rec_pattern(t) += phi_t(month_mod_12);
    }
  }
  vector<Type> rec = rec_pattern*R0; // Recruitment used in warmup

  // Aggregate spawning pattern according to Number_months_per_timestep
  vector<Type> sp_pattern(dt);
  for (int t=0;t<dt;t++){
    sp_pattern(t) = 0;
    for (int m=t*Number_months_per_timestep;m<Number_months_per_timestep*(t+1);m++){
      int month_mod_12 = m % 12;
      sp_pattern(t) += proportion_spawning(month_mod_12);
    }
  }

  // Equilibrium warmup
  Type zero = 0;
  vector<Type> B_initial(12);
  vector<Type> N_initial(12);
  vector<Type> R_initial(12);
  vector<Type> F_pattern_eq(dt);
  B_initial.fill(B0);
  N_initial.fill(N0);
  R_initial.fill(R0);
  F_pattern_eq.fill(1.0);

  PopulationVectors<Type> equilibrium = PopulationDynamics(num_years_warmup, B_initial, N_initial, R_initial, R0, R0, rho, weight_at_recruitment, mid_surv,
                                                     sp_pattern, dt, h, rec_pattern, M, zero, R0, cpue, ctch, false, true, eta, none,
                                                     log_R_star, first_year_rec_devs, first_year_catch, F_pattern_eq, F_pattern_eq); // sp_pattern is a dummy for F_pattern

  // // Outputs from equilibrium inform recruitment dynamics
  vector<Type> sp_eq(12);
  vector<Type> B_eq(12);
  vector<Type> N_eq(12);
  vector<Type> R_eq(12);
  for (int i=0;i<12;i++){
    sp_eq(i) = equilibrium.Spawners(T_warmup-12+i);
    B_eq(i) = equilibrium.Biomass(T_warmup-12+i);
    N_eq(i) = equilibrium.Numbers(T_warmup-12+i);
    R_eq(i) = equilibrium.Recruits(T_warmup-12+i);
  }

  Type sp_total = 0;
  for (int i=0;i<dt;i++){
    sp_total = sp_total + sp_eq(i);
  }

  // Set up for historical fishing
  vector<Type> R(T+2);
  vector<Type> B(T+2);
  vector<Type> SpawnB(T);
  vector<Type> N(T+2);
  vector<Type> S(T+2);
  vector<Type> SpYr(num_years+1);
  vector<Type> B_annual(num_years);
  vector<Type> SpawnB_annual_ratio(num_years);
  vector<Type> AnnualYield(num_years);
  vector<Type> F_annual(num_years);
  vector<Type> F_timestep(T+2);
  vector<Type> HarvestRate(T+2);
  vector<Type> MidBiomass(T);
  vector<Type> Yield(T+2);
  vector<Type> Effort(T+2);

  // Set up for recruitment and survival
  vector<Type> RecDev(num_years);
  vector<Type> Rec_BH(num_years);

  RecDev.fill(0.0);
  Rec_BH.fill(R0);

  if (rec_dev_type == "redundant"){
    PopulationVectors<Type> HistoricalPeriod = PopulationDynamics(num_years, B_eq, N_eq, R_eq, sp_total, sp_total, rho, weight_at_recruitment, mid_surv,
                                                            sp_pattern, dt, h, rec_pattern, M, F_initial, R0, cpue, ctch, true, false, eta, redundant,
                                                            log_R_star, first_year_rec_devs, first_year_catch, F_pattern_eq, F_pattern_eq); // F_pattern_eq does not need to be used
    B_annual = HistoricalPeriod.AnnualBiomass;
    F_annual = HistoricalPeriod.AnnualF;
    SpawnB_annual_ratio = HistoricalPeriod.AnnualSpawnersRatio;
    AnnualYield - HistoricalPeriod.AnnualYield;
    F_timestep = HistoricalPeriod.FishingMortality;

    B = HistoricalPeriod.Biomass;
    MidBiomass = HistoricalPeriod.MidBiomass;
    N = HistoricalPeriod.Numbers;
    SpawnB = HistoricalPeriod.Spawners;
    S = HistoricalPeriod.Survivors;
    R = HistoricalPeriod.BevertonHoltRecruits;
    RecDev = HistoricalPeriod.RecruitmentDeviation;
    SpYr = HistoricalPeriod.AnnualSpawners;
    Rec_BH = HistoricalPeriod.AnnualRecruits;
    HarvestRate = HistoricalPeriod.HarvestRate;
    Yield = HistoricalPeriod.Yield;
    Effort = HistoricalPeriod.Effort;

  } else if (rec_dev_type == "barycentric") {
    PopulationVectors<Type> HistoricalPeriod = PopulationDynamics(num_years, B_eq, N_eq, R_eq, sp_total, sp_total, rho, weight_at_recruitment, mid_surv,
                                                            sp_pattern, dt, h, rec_pattern, M, F_initial, R0, cpue, ctch, true, false, eta, barycentric,
                                                            log_R_star, first_year_rec_devs, first_year_catch, F_pattern_eq, F_pattern_eq); // F_pattern_eq does not need to be used
    B_annual = HistoricalPeriod.AnnualBiomass;
    F_annual = HistoricalPeriod.AnnualF;
    SpawnB_annual_ratio = HistoricalPeriod.AnnualSpawnersRatio;
    AnnualYield - HistoricalPeriod.AnnualYield;
    F_timestep = HistoricalPeriod.FishingMortality;

    B = HistoricalPeriod.Biomass;
    MidBiomass = HistoricalPeriod.MidBiomass;
    N = HistoricalPeriod.Numbers;
    SpawnB = HistoricalPeriod.Spawners;
    S = HistoricalPeriod.Survivors;
    R = HistoricalPeriod.Recruits;
    RecDev = HistoricalPeriod.RecruitmentDeviation;
    SpYr = HistoricalPeriod.AnnualSpawners;
    Rec_BH = HistoricalPeriod.AnnualRecruits;
    HarvestRate = HistoricalPeriod.HarvestRate;
    Yield = HistoricalPeriod.Yield;
    Effort = HistoricalPeriod.Effort;

  } else {
    PopulationVectors<Type> HistoricalPeriod = PopulationDynamics(num_years, B_eq, N_eq, R_eq, sp_total, sp_total, rho, weight_at_recruitment, mid_surv,
                                                            sp_pattern, dt, h, rec_pattern, M, F_initial, R0, cpue, ctch, true, false, eta, none,
                                                            log_R_star, first_year_rec_devs, first_year_catch, F_pattern_eq, F_pattern_eq); // F_pattern_eq does not need to be used
    B_annual = HistoricalPeriod.AnnualBiomass;
    F_annual = HistoricalPeriod.AnnualF;
    SpawnB_annual_ratio = HistoricalPeriod.AnnualSpawnersRatio;
    AnnualYield - HistoricalPeriod.AnnualYield;
    F_timestep = HistoricalPeriod.FishingMortality;

    B = HistoricalPeriod.Biomass;
    MidBiomass = HistoricalPeriod.MidBiomass;
    N = HistoricalPeriod.Numbers;
    SpawnB = HistoricalPeriod.Spawners;
    S = HistoricalPeriod.Survivors;
    R = HistoricalPeriod.Recruits;
    RecDev = HistoricalPeriod.RecruitmentDeviation;
    SpYr = HistoricalPeriod.AnnualSpawners;
    Rec_BH = HistoricalPeriod.AnnualRecruits;
    HarvestRate = HistoricalPeriod.HarvestRate;
    Yield = HistoricalPeriod.Yield;
    Effort = HistoricalPeriod.Effort;

  }

  // Calculate annual biomass ratio
  vector<Type> B_annual_ratio(num_years);
  B_annual_ratio = B_annual/(sum(B_eq)/12);

  // Count number of timesteps with non-zero catch rates
  vector<Type> count(cpue.rows());
  for (int f=0;f<cpue.rows();f++){
    count(f) = 0;
    for (int i=0;i<cpue.cols();i++){
      if (cpue(f,i) > 0){
        count(f) += 1;
      }
    }
  }

  // Calculate catchability for each catch rate timeseries
  vector<Type> qratio(cpue.rows());
  for (int f=0;f<cpue.rows();f++){
    qratio(f) = 0;
    for (int i=0;i<cpue.cols();i++){
      if (cpue(f,i) > 0){
        qratio(f) += log(cpue(f,i)/(B(2+i)*((1-S(2+i))/(-log(S(2+i))))));
      }
    }
  }

  // Mean catchability across years
  vector<Type> log_q(cpue.rows());
  log_q = qratio/count;

  // Set up predicted catch rate variables
  vector<Type> t_seq = 2*M_PI*month_sequence/12;
  matrix<Type> pred_cpue(cpue.rows(),cpue.cols());
  vector<Type> q(dt);
  vector<Type> cpueLL(cpue.rows());
  Type diff=0;

  for (int f=0;f<cpue.rows();f++){

    // Catchability pattern
    vector<Type> q_temp = exp(log_q(f)+q1*cos(t_seq)+q2*sin(t_seq)); // Courtney 2014 p135
    for (int t=0;t<dt;t++){
      q(t) = 0;
      for (int m=t*Number_months_per_timestep;m<Number_months_per_timestep*(t+1);m++){
        int month_mod_12 = m % 12;
        q(t) += q_temp(month_mod_12)/Number_months_per_timestep;
      }
    }

    int j = 0;
    for (int i=0;i<cpue.cols();i++){
      j = i % dt;
      pred_cpue(f,i) = q(j)*B(2+i)*((1-S(2+i))/(-log(S(2+i))));
    }

    // Likelihood components:
    // Catch rate log-likelihood
    if (use_cpue_sd){
      diff=0;
      for (int i=0;i<cpue.cols();i++){
        if (cpue(f,i) > 0){
          diff += log(cpue_sd(f,i))+0.5*pow(log(cpue(f,i))-log(pred_cpue(f,i)),2)/pow(cpue_sd(f,i),2);
        }}
      cpueLL(f) = diff;
    } else {
      diff=0;
      for (int i=0;i<cpue.cols();i++){
        if (cpue(f,i) > 0){
          diff += pow(log(cpue(f,i))-log(pred_cpue(f,i)),2);
        }}
      cpueLL(f) = count(f)*0.5*lsigmaI_sq+diff/(2*exp(lsigmaI_sq));
    }
  }

  Type sigmaI_sq = exp(lsigmaI_sq);
  Type sigmaI = sqrt(sigmaI_sq);

  // Recruitment deviations log-likelihood
  Type RecDevLL = 0;
  Type sigmaR = sqrt(sigmaR_sq);
  if (rec_dev_type == "barycentric"){
    for (int i=0;i<zeta.size();i++){
      RecDevLL -= dnorm(zeta(i),Type(0),sigmaR,true);
    }
  }
  // State-space process model
  if (rec_dev_type == "redundant"){
    for (int i=0;i<log_R_star.size();i++){
      RecDevLL -= dnorm(log_R_star(i),log(Rec_BH(i+first_year_rec_devs-first_year_catch)),sigmaR,true); // dnorm(RecDev(i),Type(0),sigmaR,true);
    }
  }

  // Absolute biomass log-likelihood
  Type biomassLL = 0;
  for (int i=0;i<absolute_biomass.size();i++){
    if (absolute_biomass(i) > Type(1)) {
      biomassLL -= dnorm(absolute_biomass(i),B(i+2),absolute_biomass_sd(i),true);
    }
  }

  // Penalty 1: prevent catch from exceeding exploitable biomass
  Type std1 = 0.1;
  Type a = 0;
  for (int i=1;i<T;i++){
    if (ctch(i) > MidBiomass(i)) {
      a += pow((log(ctch(i)/1000) - log(MidBiomass(i)/1000))/std1,2);
    }
  }
  Type penLL1 = 0.5*a;

  // Penalty 2: prevent catch from exceeding recruits
  Type lambda2 = 0;
  vector<Type> annual_catch(num_years);
  vector<Type> annual_recruitment(num_years);
  for (int i=0;i<num_years;i++){ // each year
    annual_catch(i) = 0;
    annual_recruitment(i) = 0;
    for (int m=0; m<dt;m++){ // each month
      annual_catch(i) += ctch(i*dt+m);
      annual_recruitment(i) += R(i*dt+m+2);
    }
    if (minimum_annual_harvest > annual_catch(i)/(annual_recruitment(i)*weight_at_recruitment(1))){
      lambda2 += pow(minimum_annual_harvest-annual_catch(i)/(annual_recruitment(i)*weight_at_recruitment(1)),recruit_penalty_exponent)/recruit_penalty_strength;
    }
  }
  Type penLL2 = use_recruit_penalty*lambda2;

  // Priors
  Type prior_xi = 0.5*(pow(xi-prior_mean_xi,2)/pow(prior_sd_xi,2));
  Type prior_mu = 0.5*(pow(mu-prior_mean_mu,2)/pow(prior_sd_mu,2));
  Type prior_k = 0.5*(pow(k-prior_mean_k,2)/pow(prior_sd_k,2));

  // End likelihood components

  // Start projections
  if (do_projections){


    // Extract economic data from list
    // Type econ_cL = economic_data.cL(0);
    // Type econ_cM = economic_data.cM(0);
    // Type econ_cK = economic_data.cK(0);
    // Type econ_cF = economic_data.cF(0);
    // Type econ_cO = economic_data.cO(0);
    // Type econ_v = economic_data.v(0);

    // Type econ_W = economic_data.W(0);
    // Type econ_K = economic_data.K(0);
    // Type econ_rho = economic_data.prop(0);

    // Type econ_B = economic_data.B(0);
    // Type econ_dmean = economic_data.dmean(0);
    // Type econ_o = economic_data.o(0);
    // Type econ_d = economic_data.d(0);

    // Five-year average seasonal F pattern
    vector<Type> F_pattern(dt);
    F_pattern.fill(0.0);
    for (int n=1;n<6;n++){
      for (int i=0;i<dt;i++){
        F_pattern(i) += F_timestep(T - n*dt + i);
      }
    }
    F_pattern = F_pattern/sum(F_pattern);

    //Set up projection
    int T_msy = projection_years*dt;
    int report_biomass = 0;
    int report_yield = 1;
    int report_economic_value = 2;
    int report_target = 3;
    vector<Type> Fvec(1);
    Fvec(0) = F_initial;

    //Define functor for yield optimisation (report_yield)
    Functor<TMBad::ad_aug> F(F_pattern,mid_surv,projection_years,dt,R_eq,B_eq,N_eq,surv,sp_total,R0,rec_pattern,rho,h,weight_at_recruitment,sp_pattern,q,target_relative_biomass,report_yield,economic_data.cL,economic_data.cM,economic_data.cK,economic_data.cF,economic_data.cO,economic_data.v,economic_data.W,economic_data.K,economic_data.prop,economic_data.B,economic_data.dmean,economic_data.o,economic_data.d);

    //Define functor for biomass optimisation (report_target)
    Functor<TMBad::ad_aug> G(F_pattern,mid_surv,projection_years,dt,R_eq,B_eq,N_eq,surv,sp_total,R0,rec_pattern,rho,h,weight_at_recruitment,sp_pattern,q,target_relative_biomass,report_target,economic_data.cL,economic_data.cM,economic_data.cK,economic_data.cF,economic_data.cO,economic_data.v,economic_data.W,economic_data.K,economic_data.prop,economic_data.B,economic_data.dmean,economic_data.o,economic_data.d);

    //Define functor for economic value optimisation (report_economic_value)
    Functor<TMBad::ad_aug> H(F_pattern,mid_surv,projection_years,dt,R_eq,B_eq,N_eq,surv,sp_total,R0,rec_pattern,rho,h,weight_at_recruitment,sp_pattern,q,target_relative_biomass,report_economic_value,economic_data.cL,economic_data.cM,economic_data.cK,economic_data.cF,economic_data.cO,economic_data.v,economic_data.W,economic_data.K,economic_data.prop,economic_data.B,economic_data.dmean,economic_data.o,economic_data.d);

    //Newton minimisation for Fmsy and Ftarg
    vector<Type> Fmsy = Newton(F, Fvec, cfg);
    vector<Type> Ftarg = Newton(G, Fvec, cfg);
    vector<Type> Fmey = Newton(H, Fvec, cfg);

    //Define functor for yield and biomass to return double
    Functor<Type> yield(F_pattern,mid_surv,projection_years,dt,R_eq,B_eq,N_eq,surv,sp_total,R0,rec_pattern,rho,h,weight_at_recruitment,sp_pattern,q,target_relative_biomass,report_yield,economic_data.cL,economic_data.cM,economic_data.cK,economic_data.cF,economic_data.cO,economic_data.v,economic_data.W,economic_data.K,economic_data.prop,economic_data.B,economic_data.dmean,economic_data.o,economic_data.d);
    Functor<Type> biomass(F_pattern,mid_surv,projection_years,dt,R_eq,B_eq,N_eq,surv,sp_total,R0,rec_pattern,rho,h,weight_at_recruitment,sp_pattern,q,target_relative_biomass,report_biomass,economic_data.cL,economic_data.cM,economic_data.cK,economic_data.cF,economic_data.cO,economic_data.v,economic_data.W,economic_data.K,economic_data.prop,economic_data.B,economic_data.dmean,economic_data.o,economic_data.d);
    Functor<Type> value(F_pattern,mid_surv,projection_years,dt,R_eq,B_eq,N_eq,surv,sp_total,R0,rec_pattern,rho,h,weight_at_recruitment,sp_pattern,q,target_relative_biomass,report_economic_value,economic_data.cL,economic_data.cM,economic_data.cK,economic_data.cF,economic_data.cO,economic_data.v,economic_data.W,economic_data.K,economic_data.prop,economic_data.B,economic_data.dmean,economic_data.o,economic_data.d);

    //Return yield and biomass @ Fmsy
    Type msy = -yield(Fmsy);
    Type Bmsy = biomass(Fmsy);

    //Return yield and biomass @ Fmey
    Type mey = -yield(Fmey);
    Type Vmey = -value(Fmey);
    Type Bmey = biomass(Fmey);

    //Return yield and biomass @ Ftarg
    Type Ytarg = -yield(Ftarg);
    Type Btarg = biomass(Ftarg);

    // Make a yield curve (MSY)
    vector<Type> F_values(100);
    vector<Type> yield_values(100);
    vector<Type> B_values(100);
    for (int i=0;i<100;i++){
      F_values(i) = 3*i*Fmsy(0)/100;
      vector<Type> Fmort(1); // functor wants a vector
      Fmort = F_values(i);
      yield_values(i) = -yield(Fmort);
      B_values(i) = biomass(Fmort);
    }

    // Make a value curve (MEY)
    vector<Type> MEY_F_curve(100);
    vector<Type> MEY_yield_curve(100);
    vector<Type> MEY_biomass_curve(100);
    vector<Type> MEY_value_curve(100);
    for (int i=0;i<100;i++){
      MEY_F_curve(i) = 3*i*Fmey(0)/100;
      vector<Type> Fmort(1); // functor wants a vector
      Fmort = F_values(i);
      MEY_yield_curve(i) = -yield(Fmort);
      MEY_biomass_curve(i) = biomass(Fmort);
      MEY_value_curve(i) = -value(Fmort);
    }

    vector<Type> B_start_projection(12);
    vector<Type> R_start_projection(12);
    vector<Type> N_start_projection(12);
    B_start_projection.fill(B(T));
    R_start_projection.fill(R(T));
    N_start_projection.fill(N(T));
    for (int t=0;t<12;t++){
      B_start_projection(t) = B(T+1-11+t);
      N_start_projection(t) = N(T+1-11+t);
      R_start_projection(t) = R(T+1-11+t);
    }

    ADREPORT(B_start_projection);
    ADREPORT(N_start_projection);
    ADREPORT(R_start_projection);
    ADREPORT(B_eq);
    ADREPORT(N_eq);
    ADREPORT(R_eq);

    PopulationVectors<Type> TargetProjection = PopulationDynamics(projection_years, B_start_projection, N_start_projection, R_start_projection,
                                                      SpYr(num_years-1), sp_total, rho, weight_at_recruitment, mid_surv, // Start projection from current biomass instead of equilibrium
                                                      sp_pattern, dt, h, rec_pattern, M, Ftarg(0), R0, cpue, ctch, false, false, eta, none,
                                                      log_R_star, first_year_rec_devs, first_year_catch, F_pattern, q); // B_initial etc need to be vectors

    PopulationVectors<Type> MSYProjection = PopulationDynamics(projection_years, B_start_projection, N_start_projection, R_start_projection,
                                                      SpYr(num_years-1), sp_total, rho, weight_at_recruitment, mid_surv, // Start projection from current biomass instead of equilibrium
                                                      sp_pattern, dt, h, rec_pattern, M, Fmsy(0), R0, cpue, ctch, false, false, eta, none,
                                                      log_R_star, first_year_rec_devs, first_year_catch, F_pattern, q); // B_initial etc need to be vectors

    PopulationVectors<Type> MEYProjection = PopulationDynamics(projection_years, B_start_projection, N_start_projection, R_start_projection,
                                                      SpYr(num_years-1), sp_total, rho, weight_at_recruitment, mid_surv, // Start projection from current biomass instead of equilibrium
                                                      sp_pattern, dt, h, rec_pattern, M, Fmey(0), R0, cpue, ctch, false, false, eta, none,
                                                      log_R_star, first_year_rec_devs, first_year_catch, F_pattern, q); // B_initial etc need to be vectors

    PopulationVectors<Type> HarvestProjection = PopulationDynamics(projection_years, B_start_projection, N_start_projection, R_start_projection,
                                                      SpYr(num_years-1), sp_total, rho, weight_at_recruitment, mid_surv, // Start projection from current biomass instead of equilibrium
                                                      sp_pattern, dt, h, rec_pattern, M, Fmsy(0), R0, cpue, projection_harvest, true, false, eta, none,
                                                      log_R_star, first_year_rec_devs, first_year_catch, F_pattern, q); // B_initial etc need to be vectors

    // Report projection results
    REPORT(msy);
    REPORT(Fmsy);
    REPORT(Bmsy);
    REPORT(Ftarg);
    REPORT(Btarg);
    REPORT(Ytarg);
    REPORT(mey);
    REPORT(Fmey);
    REPORT(Vmey);
    REPORT(Bmey);
    REPORT(F_pattern);
    REPORT(F_values);
    REPORT(yield_values);
    REPORT(B_values);
    REPORT(MEY_biomass_curve);
    REPORT(MEY_F_curve);
    REPORT(MEY_value_curve);
    REPORT(MEY_yield_curve);
    REPORT(TargetProjection.AnnualYield);
    REPORT(TargetProjection.AnnualBiomass);
    REPORT(TargetProjection.Biomass);
    REPORT(TargetProjection.Yield);
    REPORT(MSYProjection.AnnualYield);
    REPORT(MSYProjection.AnnualBiomass);
    REPORT(MSYProjection.Biomass);
    REPORT(MSYProjection.Effort);
    REPORT(MSYProjection.MidBiomass);
    REPORT(MSYProjection.Yield);
    REPORT(MEYProjection.AnnualYield);
    REPORT(MEYProjection.AnnualBiomass);
    REPORT(MEYProjection.Biomass);
    REPORT(MEYProjection.Effort);
    REPORT(MEYProjection.MidBiomass);
    REPORT(MEYProjection.Yield);
    REPORT(HarvestProjection.AnnualBiomass);
    REPORT(HarvestProjection.AnnualF);
    REPORT(HarvestProjection.AnnualYield);
    REPORT(HarvestProjection.Biomass);
    REPORT(HarvestProjection.Effort);
    REPORT(HarvestProjection.MidBiomass);
    REPORT(HarvestProjection.FishingMortality);
    REPORT(HarvestProjection.Yield);

    ADREPORT(msy);
    ADREPORT(Fmsy);
    ADREPORT(Bmsy);
    ADREPORT(Ftarg);
    ADREPORT(Btarg);
    ADREPORT(Ytarg);
    ADREPORT(mey);
    ADREPORT(Fmey);
    ADREPORT(Vmey);
    ADREPORT(Bmey);
    ADREPORT(F_pattern);
    ADREPORT(F_values);
    ADREPORT(yield_values);
    ADREPORT(B_values);
    ADREPORT(MEY_biomass_curve);
    ADREPORT(MEY_F_curve);
    ADREPORT(MEY_value_curve);
    ADREPORT(MEY_yield_curve);
    ADREPORT(TargetProjection.AnnualYield);
    ADREPORT(TargetProjection.AnnualBiomass);
    ADREPORT(TargetProjection.Biomass);
    ADREPORT(TargetProjection.Yield);
    ADREPORT(MSYProjection.AnnualYield);
    ADREPORT(MSYProjection.AnnualBiomass);
    ADREPORT(MSYProjection.Biomass);
    ADREPORT(MSYProjection.Effort);
    ADREPORT(MSYProjection.MidBiomass);
    ADREPORT(MSYProjection.Yield);
    ADREPORT(MEYProjection.AnnualYield);
    ADREPORT(MEYProjection.AnnualBiomass);
    ADREPORT(MEYProjection.Biomass);
    ADREPORT(MEYProjection.Effort);
    ADREPORT(MEYProjection.MidBiomass);
    ADREPORT(MEYProjection.Yield);
    ADREPORT(HarvestProjection.AnnualBiomass);
    ADREPORT(HarvestProjection.AnnualYield);
    ADREPORT(HarvestProjection.AnnualF);
    ADREPORT(HarvestProjection.Biomass);
    ADREPORT(HarvestProjection.Effort);
    ADREPORT(HarvestProjection.MidBiomass);
    ADREPORT(HarvestProjection.FishingMortality);
    ADREPORT(HarvestProjection.Yield);

  }

  // Objective function:
  Type LL = sum(cpueLL) + RecDevLL + biomassLL + prior_xi + prior_mu + prior_k + penLL1 + penLL2; // for use with minimizing algorithm

  // Reported values in alphabetical order
  ADREPORT(annual_recruitment);
  ADREPORT(B);
  ADREPORT(B_annual);
  ADREPORT(B_annual_ratio);
  ADREPORT(biomassLL);
  ADREPORT(cpueLL);
  ADREPORT(Effort);
  ADREPORT(F_annual);
  ADREPORT(F_timestep);
  ADREPORT(h);
  ADREPORT(HarvestRate);
  ADREPORT(k);
  ADREPORT(LL);
  ADREPORT(log_q);
  ADREPORT(log_R_star);
  ADREPORT(lsigmaI_sq);
  ADREPORT(lsigmaR_sq);
  ADREPORT(M);
  ADREPORT(MidBiomass);
  ADREPORT(mu);
  ADREPORT(N);
  ADREPORT(penLL1);
  ADREPORT(penLL2);
  ADREPORT(pred_cpue);
  ADREPORT(prior_k);
  ADREPORT(prior_mu);
  ADREPORT(prior_xi);
  ADREPORT(qratio);
  ADREPORT(q1);
  ADREPORT(q2);
  ADREPORT(R);
  ADREPORT(Rec_BH);
  ADREPORT(R0);
  ADREPORT(Rinit);
  ADREPORT(rec_pattern);
  ADREPORT(RecDev);
  ADREPORT(RecDevLL);
  ADREPORT(rho);
  ADREPORT(S);
  ADREPORT(sigmaI);
  ADREPORT(sigmaR);
  ADREPORT(sigmaI_sq);
  ADREPORT(sigmaR_sq);
  ADREPORT(sp_total);
  ADREPORT(SpYr);
  ADREPORT(SpawnB);
  ADREPORT(SpawnB_annual_ratio);
  ADREPORT(Yield);
  ADREPORT(xi);
  ADREPORT(zeta);
  ADREPORT(zeta_matrix);

  REPORT(annual_catch);
  REPORT(annual_recruitment);
  REPORT(B);
  REPORT(B_annual);
  REPORT(B_annual_ratio);
  REPORT(biomassLL);
  REPORT(cpueLL);
  REPORT(dt);
  REPORT(Effort);
  REPORT(eta);
  REPORT(F_annual);
  REPORT(F_timestep);
  REPORT(HarvestRate);
  REPORT(LL);
  REPORT(log_q);
  REPORT(M);
  REPORT(MidBiomass);
  REPORT(N);
  REPORT(penLL1);
  REPORT(penLL2);
  REPORT(phi_t);
  REPORT(pred_cpue);
  REPORT(qratio);
  REPORT(R);
  REPORT(Rec_BH);
  REPORT(R0);
  REPORT(rec_pattern);
  REPORT(RecDev);
  REPORT(RecDevLL);
  REPORT(rho);
  REPORT(S);
  REPORT(sigmaI);
  REPORT(sigmaR);
  REPORT(sigmaI_sq);
  REPORT(sigmaR_sq);
  REPORT(sp_pattern);
  REPORT(sp_total);
  REPORT(SpawnB);
  REPORT(SpawnB_annual_ratio);
  REPORT(SpYr);
  REPORT(Yield);
  REPORT(weight_at_recruitment);

  // Equilibrium
  REPORT(equilibrium.AnnualRecruits);
  REPORT(equilibrium.AnnualSpawners);
  REPORT(equilibrium.Biomass);
  REPORT(equilibrium.Numbers);
  REPORT(equilibrium.Survivors);
  REPORT(equilibrium.Spawners);
  REPORT(equilibrium.Recruits);
  REPORT(equilibrium.RecruitmentDeviation);

  ADREPORT(equilibrium.AnnualRecruits);
  ADREPORT(equilibrium.AnnualSpawners);
  ADREPORT(equilibrium.Biomass);
  ADREPORT(equilibrium.Numbers);
  ADREPORT(equilibrium.Survivors);
  ADREPORT(equilibrium.Spawners);
  ADREPORT(equilibrium.Recruits);
  ADREPORT(equilibrium.RecruitmentDeviation);

  return LL;
}

#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR this
