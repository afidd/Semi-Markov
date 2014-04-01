==========================
Observations of the State
==========================

After each transition, we observe the state of the
system using the GSPN object and the GSPNState object.

.. cpp:class:: afidd::smv::GSPNState\<Marking,UserSate\>
   
   This class is the state of the system.

.. cpp:member:: Marking\<Place,Tokens\> GSPNState::marking

   The marking is saved as a `marking` member. The member function
   `state.marking.modified()` will list which places changed
   during the last transition.

.. cpp:function:: double GSPNState::current_time() const

   The absolute time of the system.

.. cpp:member:: UserState GSPNState::user

   This member is an instance of the class `UserState.`


For example, an instance of the following class observes
when the number of infected individuals within a population
has passed a given threshold::

	template<typename GSPN>
	class SIROutputFunction
	{
	public:
	  typedef void result_type;

	private:
	  size_t _pop_cnt;
	  size_t _ind_cnt;
	  size_t _threshold;
	  std::vector<size_t> _infected_per_pop;
	  std::vector<double> _first_passage_time;
	  std::vector<bool> _passed;

	public:
	  SIROutputFunction(size_t population_cnt,
	    size_t individuals_per_metapopulation, size_t threshold_for_passage)
	  : _pop_cnt(population_cnt), _ind_cnt(individuals_per_metapopulation),
	    _threshold(threshold_for_passage),
	    _first_passage_time(population_cnt, 0), _passed(population_cnt, false),
	    _infected_per_pop(population_cnt, 0)
	  {

	  }



	  result_type operator()(const GSPN& gspn, const SIRState& state)
	  {
	    auto& modified=state.marking.modified();
	    for (auto place_idx : modified)
	    {
	      // Translate the internal place index into the PlaceKey
	      // used to define the GSPN when building the graph.
	      SIRPlace p=gspn.vertex_place(place_idx);
	      if (p.disease==1)
	      {
	        bool filled=(length<0>(state.marking, place_idx)>0);
	        if (filled)
	        {
	          _infected_per_pop[p.metapop]+=1;
	          if (_infected_per_pop[p.metapop]==_threshold)
	          {
	            _first_passage_time[p.metapop]=state.current_time();
	            _passed[p.metapop]=true;
	            BOOST_LOG_TRIVIAL(info)<<"Population "<<p.metapop
	              <<" at time "<<state.current_time();
	          }
	        }
	      }
	    }
	  }
	};

We use such an object by calling it within the main loop::

	  size_t transition_cnt=0;
	  auto nothing=[](SIRState&)->void {};
	  for ( ;
	    std::get<1>(next)<std::numeric_limits<double>::max();
	    next=propagate_competing_processes(system, nothing, rng))
	  {
	    BOOST_LOG_TRIVIAL(debug) << "trans " << std::get<0>(next) << " time " <<
	        std::get<1>(next);
	    ++transition_cnt;
	    output(gspn, state);
	  }
	  