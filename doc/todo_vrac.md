# TODO (vrac)

* Refactor proxies into toolbox lib.
    * Maybe make `gather` and `get` template functions instead of virtual methods?
    It would allow template params for everything (where virtual methods can't be template) and be more consistent overall.
    * Refactor `Proxy` so it supports various cases regarding `get` and `gather` parameters.
    E.g. allow proxies with `T get(int)` + `void gather()` and allow proxies with `T get(int)` + `void gather(int)`.