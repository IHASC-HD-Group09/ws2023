#ifndef FILL_VECTOR_EXPERIMENT_H
#define FILL_VECTOR_EXPERIMENT_H

#include <vector>

#include "experimet.h"

std::vector<int> veccy;

class FillVectorExperiment : Experiment {
   private:
    int n;

   public:
    virtual ~FillVectorExperiment() = default;
    FillVectorExperiment(int n = 1024) : n(n) {  // NOLINT(readability-magic-numbers)
        veccy.resize(n);
    }

    void run() const override {
        for (int i = 0; i < n; i++) {
            veccy[i] = i;
        }
    }
};

#endif