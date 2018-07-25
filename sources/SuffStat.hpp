
#ifndef SUFFSTAT_H
#define SUFFSTAT_H

/**
 * \brief A dummy superclass for all sufficient statistics. Does not implement
 * anything in itself and does not provide any specific interface.
 */

class SuffStat {
  public:
    SuffStat() {}
    virtual ~SuffStat() {}
};

#endif
