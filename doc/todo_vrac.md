# TODO (vrac)

* Refactor proxies (ie, ss interface) into toolbox lib.
    * Maybe make `gather` and `get` template functions instead of virtual methods?
    It would allow template params for everything (where virtual methods can't be template) and be more consistent overall.
    * Refactor `Proxy` so it supports various cases regarding `get` and `gather` parameters.
    E.g. allow proxies with `T get(int)` + `void gather()` and allow proxies with `T get(int)` + `void gather(int)`.
* Move everything that is not legacy-bayescode-dependent out of the repo
    * Goal is that e.g. simple non-phylogenetics MCMC projects don't need to depend on bayescode
    * Move debugging in a "basic toolbox" repo which would be a dependency of other libs repos (toolbox, minimpl, tagged_tuple)
    * Move gommette out of the repo i.e. extract all functionnality related to recursive traversal for tracing and mpi
* Add all dependencies directly (instead of recursively)?
* Doc
    * __write explanation for tagged tuples with pro/cons and justification in project__
    * write project structure overview with libs
    * __write toolbox explanations up to the point where a basic model can be made__
* Review of dead code and "to be killed" code in sources
    * Random.hpp/cpp should disappear at some point and everything random-related should be handled by the toolbox repo
* Create new repo from experimental branch in order to keep branch number under control? Maybe just clean branches in bayescode?