#include "PathSuffStat.hpp"

ostream& operator<<(ostream& os, const PathSuffStat& suffstat)    {
    suffstat.ToStream(os);
    return os;
}

istream& operator>>(istream& is, PathSuffStat& suffstat)    {
    suffstat.FromStream(is);
    return is;
}

ostream& operator<<(ostream& os, const PathSuffStatArray& suffstat)    {
    suffstat.ToStream(os);
    return os;
}

istream& operator>>(istream& is, PathSuffStatArray& suffstat)    {
    suffstat.FromStream(is);
    return is;
}

ostream& operator<<(ostream& os, const PathSuffStatNodeArray& suffstat)    {
    suffstat.ToStream(os);
    return os;
}

istream& operator>>(istream& is, PathSuffStatNodeArray& suffstat)    {
    suffstat.FromStream(is);
    return is;
}

ostream& operator<<(ostream& os, const PathSuffStatBidimArray& suffstat)    {
    suffstat.ToStream(os);
    return os;
}

istream& operator>>(istream& is, PathSuffStatBidimArray& suffstat)    {
    suffstat.FromStream(is);
    return is;
}
