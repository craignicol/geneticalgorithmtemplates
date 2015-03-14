# Introduction #
## variables required ##

  * population size
  * mutation rate
  * number of populations
  * population crossover
  * selection method - stochastic sampling, rank, map, ...
  * crossover rate
  * chromosome type
  * generation type - replace current, add to current, ...

## functions required ##

`ctor(set_all_parameters...)`


`run(int ngenerations)`
> run the GA for next `n` generations,
> or until the population is stable,
> whichever comes first.

`vector list_best(int nbest [, int pop])`
> returns a `vector` containing the `n` fittest chromosomes in the population
> numbered `pop`, or the fittest `n` overall if undefined (or `pop` is negative)

`run_once()`
> run the GA for one generation