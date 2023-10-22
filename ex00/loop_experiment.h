#ifndef LOOP_EXPERIMENT_H
#define LOOP_EXPERIMENT_H

#include "experimet.h"

class LoopExperiment : Experiment {
   private:
    int n;

   public:
    virtual ~LoopExperiment() = default;
    LoopExperiment(int n = 1024) : n(n) {}  // NOLINT(readability-magic-numbers)

    void run() const override{
        for (int i = 0; i < n; i++) {
            // pass
        }
    }
};

#endif