#pragma once
#include <memory>
#include <string>

class IFitness;
class ISelection;
class IMutation;
class ICrossover;

class GA {
    std::unique_ptr<IFitness> fit;
    std::unique_ptr<ISelection> sel;
    std::unique_ptr<IMutation> mut;
    std::unique_ptr<ICrossover> cross;
public:
    GA(const std::string& fitness_name,
       const std::string& selection_name,
       const std::string& mutation_name,
       const std::string& crossover_name);
    ~GA();
};