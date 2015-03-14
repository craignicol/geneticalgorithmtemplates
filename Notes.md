# Introduction #

The SUS (http://en.wikipedia.org/wiki/Stochastic_universal_sampling) selection method
is slow as treacle still, but I think I've got a decent implementation
of the same idea in USER0 which isn't much slower than completely
random selection.

Still needs a couple of updates, such as parameterising over the
randomising function, and I think it could do with some refactoring.
I'd like to use a design pattern to improve the population constructor
for a start, as it's a bit unwieldy with so many parameters. Don't you
have a book that talks about some design pattern that uses a class for
initialiser so you can do things like:

`population(popInit().popsize(100).mrate(0.2));`

The only other major update I can think of is to add the alternative
selection and mutation methods, as there's still a couple missing. It
should compile fine under windows using MinGW, as there's no platform
specific stuff, I think (not sure where timing functions fall on that
score).

See the other pages on [chromosome](chromosome.md), [population](population.md) and [sus\_select](sus_select.md) for technical notes