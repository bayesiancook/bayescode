#pragma once

class AcceptanceStats {
    double ntot{0};
    double nacc{0};

  public:
    void accept() {
        ntot++;
        nacc++;
    }
    void reject() { ntot++; }
    double ratio() const { return nacc / ntot; }
};