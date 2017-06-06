

class AbstractTransitionMatrix {
  public:
    virtual ~AbstractTransitionMatrix() = default;

    virtual void BackwardPropagate(const double *down, double *up, double length) const = 0;
    virtual void ForwardPropagate(const double *up, double *down, double length) const = 0;
    virtual const double *GetStationary() const = 0;
    virtual double Stationary(int i) const = 0;

    virtual int GetNstate() const = 0;
    virtual void CorruptMatrix() = 0;
    virtual double operator()(int, int) const = 0;
    virtual const double *GetRow(int i) const = 0;

    virtual bool check() const { return true; }
};
