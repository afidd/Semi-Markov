===================
Marking Reference
===================

Local Marking
================
The local marking is the set of ways a transition can examine
or modify the marking of its input and output tokens.
It is also responsible for keeping track of changes made to the
marking, so the system can make update affected transitions.
For this, it collaborates with the complete Marking.


Helper Classes
------------------

The marking can represent colored tokens in two ways.
If a place may contain only one type of token, then it is possible
to create a place in the marking which accepts only one C++ class
type. We call each of these classes a *token layer*.
The more common definition of color is that tokens at
a marking may be partitioned by a finite set of colors. This is
called a *colored token* container.

Two classes define storage for tokens at a place in a marking.

.. cpp:class:: afidd::smv::Uncolored\<TokenType\>

   This is a storage container for tokens of the same type and
   without color.
  
.. cpp:class:: afidd::smv::color_type\<TokenType\>

   Given an arbitrary token, this decares the type of its color.

LocalMarking Class
------------------

.. cpp:class:: afidd::smv::LocalMarking\<typename... TokenContainers\>

   This class is responsible for reading and changing the subset
   of the marking associated with the input places and output
   places of a transition. The Marking for the whole system
   puts values into this class and checks it to see what changes
   have been made.

   The template argument is a list of token containers. For example::

      typedef LocalMarking<Uncolored<CowToken>,Uncolored<BatToken>> MyLocal;

   With
   this definition, the `CowToken` is layer 0 and the `BatToken`
   is layer 1.

   All methods in this class refer to token containers with an
   ordinal `size_t`. This is exactly the order used when defining
   the transition in `BuildGraph`.

.. cpp:function:: afidd::smv::LocalMarking::LocalMarking()

   There is a default constructor.

.. cpp:function:: afidd::smv::LocalMarking::add<I>(size_t place_idx, const TokenType<I>& token)
   
   Add a token to a place. The template argument is a size_t which
   specifies the token layer.

.. cpp:function:: afidd::smv::LocalMarking::remove<I,RNG>(size_t place_idx, size_t count, RNG& random_generator)

   Remove `count` tokens from a place. This version chooses the tokens
   randomly using a random number generator of the kind defined in
   std::random and Boost.

.. cpp:function:: size_t afidd::smv::LocalMarking::length<I>(size_t place_idx) const

   Returns the number of tokens.

.. cpp:function:: int afidd::smv::LocalMarking::stochiometric_coefficient(size_t place_idx) const

   Returns the stochiometric coefficient associated with this place, for
   this transition.   

.. cpp:function:: afidd::smv::LocalMarking::layer(size_t place_idx) const

   Returns the token layer for this place.

.. cpp:function:: std::tuple\<F::result_type,bool\> afidd::smv::LocalMarking::get\<I,F\>(size_t place_idx, const F& functor) const

   This doesn't get a token. It operates on the container with a functor of
   type `F` and returns the result of that functor. For instance, you might
   want to know the average of the `birthday` property of a set of tokens
   on the second input arc of a transition::

      auto average=local_marking.template get<0>(1,
        [](const Uncolored<Children>& kids)->double {
            double total=0.0;
            for (const auto& k : kids) {
                total+=k.birthday;
            }
            return total/kids.size();
        });

   The `get` function returns a `std::tuple\<double,bool\>`, where the first
   type is the return type of your function and the second is whether
   there were any tokens at the place.

.. cpp:function:: afidd::smv::LocalMarking::get_token\<I,F\>(size_t place_idx, const F& functor) const

   This applied `functor` to the token at layer I, and returns an
   `std::tuple\<F::result_type,bool\>` where the first member of the tuple
   is the result type of the functor, `F`, and the second is whether there
   were any tokens at the place. It runs the functor on the *first* token
   in the container.

.. cpp:function:: afidd::smv::LocalMarking::move<I,J>(size_t from, size_t to, size_t count)

   Moves a token from one place to another. The two layers, I and J, are
   usually the same, but it is possible to move a token from one layer
   to another if the token types are the same for both layers. For instance,
   one may be colored an one uncolored.

.. cpp:function:: afidd::smv::LocalMarking::move<I,J,Modifier>(size_t from, size_t to, size_t count, const Modifier& functor)

   This moves `count` number of tokens from a place in the Ith token layer
   to a place in the Jth token layer (usually the same layer), and applies
   the function `functor` to each token. For example::

	   local_marking.template move<0,0>(1, 3, 1,
	     [](CowToken& bessie)->void {
	       bessie.parity+=1; // The number of times the cow gave birth.
	     });

.. cpp:function:: bool afidd::smv::LocalMarking::input_tokens_sufficient<I>() const

   Asks whether the number of tokens at each input place exceeds the
   requirements of the stochiometric coefficients to that place.
   This determines whether a transition is enabled.


.. cpp:function:: bool afidd::smv::LocalMarking::outputs_tokens_empty<I>() const

   Asks whether all output places are empty. A transition can have a policy
   that it will not enable unless output places have no tokens.

.. cpp:function:: void afidd::smv::LocalMarking::transfer_by_stochiometric_coefficient<I,RNG>(RNG& random_generator)


   This high level function looks at the stochiometric coefficient of each
   input and output arc for all tokens at layer I. It then moves tokens
   from inputs to outputs. If
   inputs have more than one token, they are chosen randomly, and input
   tokens are randomly assigned to outputs. Any extra inputs are removed
   and any extra outputs are created.

.. cpp:function:: void afidd::smv::LocalMarking::transfer_by_stochiometric_coefficient<I,RNG,AndModify>(RNG& random_generator, const AndModify& mod)

   As above, this moves tokens according to the stochiometric coefficients, but
   it also executes the AndModify functor on each token.


Marking Class
===============
The Marking class is much like the LocalMarking, but it specifies
places using a PlaceKey.

.. cpp:class:: afidd::smv::Marking<PlaceKey,typepname... TokenContainers>

   This class holds the marking for the whole system. Each place
   can contain one type of token container, specified by the
   list of TokenContainers.


.. cpp:function:: afidd::smv::Marking::Marking()

   The constructor takes no arguments.

.. cpp:function:: afidd::smv::Marking::modified()

   This returns the set of PlaceKeys of places whose tokens
   were modified in any way.

**Free Functions**

.. cpp:function:: void afidd::smv::add<I,Marking>(Marking& m, PlaceKey place, const TokenType<I>& token)

   Add a token to the container at a place. The size_t constant `I` specifies
   the layer. This is a free function.

.. cpp:function:: void afidd::smv::remove<I,Marking,RNG>(Marking& m, PlaceKey place, size_t cnt, RNG& rng)

   Remove `cnt` number of tokens from the marking at `place` in layer `I.`
   If there are more than `cnt` tokens at the place, then this uses the random
   number generator, `rng,` to select tokens.


.. cpp:function:: size_t afidd::smv::length<I,Marking>(Marking& m, PlaceKey place)

   Returns the number of tokens at a place.

.. cpp:function:: std::tuple\<F::result_type,bool\> afidd::smv::get<I,Marking,F>(const Marking& m, PlaceKey place, const F& functor)

   Apply a functor to the first token at `place.` Return a tuple with
   a) whether there was any token at the place and b) the result of the functor
   if there was.


.. cpp:function:: void afidd::smv::move<I,J,Marking>(const Marking& m, PlaceKey place_from, PlaceKey place_to, size_t count)

   Move `count` number of tokens from `place_from` in layer
   `I` to `place_to` in layer `J`.

.. cpp:function:: void afidd::smv::move<I,J,Marking,Modifier>(const Marking& m, PlaceKey place_from, PlaceKey place_to, size_t count, const Modifier& modify)

   Move `count` number of tokens from `place_from` in layer
   `I` to `place_to` in layer `J`. Additionally apply the functor `modify`
   to each moved token.

