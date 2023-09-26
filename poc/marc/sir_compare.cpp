// [[odin.dust::compare_data(cases_inc = real_type)]]
// [[odin.dust::compare_function]]
template <typename T>
typename T::real_type
compare(const typename T::real_type * state,
        const typename T::data_type& data,
        const typename T::internal_type internal,
        std::shared_ptr<const typename T::shared_type> shared,
        typename T::rng_state_type& rng_state) {
  typedef typename T::real_type real_type;
  const real_type incidence_modelled = odin(cases_inc); // state[4]
  const real_type incidence_observed = data.cases_inc;
  return dust::density::poisson(incidence_observed, incidence_modelled, true);
}
