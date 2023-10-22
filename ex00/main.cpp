#include "time_experiment.hh"

#include "accumulate_experiment.h"
#include "experimet.h"
#include "fill_vector_experiment.h"
#include "loop_experiment.h"

#include <iostream>

auto main(int argc, char* argv[]) -> int {  // NOLINT(misc-unused-parameters)
    const int N = 1024;
    LoopExperiment const loop_experiment(N);
    AccumulateExperiment const accumulate_experiment(N);
    FillVectorExperiment const fill_vector_experiment(N);


    std::cout << "LoopExperiment:       " << time_experiment(loop_experiment).second << "ms\n";
    std::cout << "AccumulateExperiment: " << time_experiment(accumulate_experiment).second << "ms\n";
    std::cout << "FillVectorExperiment: " << time_experiment(fill_vector_experiment).second << "ms\n";

    return 0;
}