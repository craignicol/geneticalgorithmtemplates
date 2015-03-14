# Introduction #
A chromosome needs the following operations:

`ctor(fn())`
> Generates a random value within the range defined by `fn()`
> default is across the whole range of the templated type

`operator+()`
> Crossover between two chromosomes

`operator~()`
> Mutation of a single chromosome

`f()`
> fitness of the current chromosome


---


## specialisations ##

numerical values:
> or should we? - is there any standard?

bitset:
> the classic GA

permutations:
> probably using containers - seperate template name
> container of any possible type

containers:
> optimized for vectors - random access required.


---


Template needs easy ability to change random function and fitness function
in template declaration. I need the latest gcc with STL to make this
compile.
optimize for vector... random access array