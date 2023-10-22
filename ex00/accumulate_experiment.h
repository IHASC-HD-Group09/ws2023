#ifndef ACCUMULATE_EXPERIMENT_H
#define ACCUMULATE_EXPERIMENT_H

#include "experimet.h"

int accumulated_result;

class AccumulateExperiment : Experiment {
   private:
    int n;

   public:
    virtual ~AccumulateExperiment() = default;
    AccumulateExperiment(int n = 1024) : n(n) {}  // NOLINT(readability-magic-numbers)

    void run() const override {
        for (int i = 0; i < n; i++) {
            accumulated_result += i;
        }
    }
};

#endif