#include <crsGA/GeneticAlgorithm.hpp>

#include <random>

#include "rf62Xtypes.h"

using namespace SDK::SCANNERS::RF62X;

std::shared_ptr<calib_table> m_table;

float computeProbability()
{
    std::random_device rd;
    std::uniform_real_distribution<float> distribution(-g_global * g_globalK / 100.0, g_global * g_globalK / 100.0);
    std::mt19937 engine(rd()); // Mersenne twister MT19937

    return distribution(engine);
}

/**
 * @brief Represents the smallest piece of information.
 * 
 * A bit of the bitstring of the onemax algorithm
 * 
 */
class Gen : public crsGA::IGen
{
  public:
    // Gen data (1 bit)
    float coefficient = 0;

  public:
    Gen() = default;
    virtual ~Gen() = default;
    // How to mutate this gen 
    virtual void mutate(const crsGA::UserData *, int index) override
    {
        if (index >= 0 && index <= 2)
            coefficient += 1.0 * computeProbability();
        else if (index >= 3 && index <= 5)
            coefficient += 3.0 * computeProbability();
        else if (index >= 6 && index <= 9)
            coefficient += 9.0 * computeProbability();
        else if (index >= 10 && index <= 14)
            coefficient += 27.0 * computeProbability();
        else if (index >= 15)
            coefficient += 27.0 * 3 * computeProbability();
    }
    // Random creation, used for initialization
    virtual void random(const crsGA::UserData *) override
    {
        coefficient = computeProbability();
    }
    // Just some friendly way of serializing the information
    friend std::ostream &operator<<(std::ostream &os, const Gen &gen)
    {
        os << static_cast<float>(gen.coefficient) << " ";
        return os;
    }
};

/**
 * @brief Defines how to compute the fitness of a chromosome
 * 
 * We return -sum(bitstring), since the algorithm tries to 
 * minimize the fittnes value
 */
class ComputeFitness : public crsGA::ComputeFitnessPolicy<Gen>
{
  public:
    virtual float computeFitness(const std::vector<Gen> &genes,
                                 const crsGA::UserData *) const override
    {
        float fitness = 0.0f;
        std::vector<int16_t> Xd = m_table->getX();
        //std::vector<int16_t> Zd = m_table->getZ();

        for (int z = 0; z < 2*m_table->getHeight(); z++)
        {
            for (int x = 0; x < m_table->getWidth(); x++)
            {
                if (Xd[z*m_table->getWidth() + x] == -32768)
                    continue;

                double xx = 0.000001 * (x*x);
                double xxx = 0.000001 * (xx * x);
                double xxxx = 0.000001 * (xx * xx);
                double zz = 0.000001 * (z*z);
                double zzz = 0.000001 * (zz * z);
                double zzzz = 0.000001 * (zz * zz);
                double xz = 0.000001 * (x*z);
                double xxz = 0.000001 * (xx * z);
                double xxxz = 0.000001 * (xxx * z);
                double xzz = 0.000001 * (x * zz);
                double xxzz = 0.000001 * (xx * zz);
                double xzzz = 0.000001 * (xz * zz);
                double zzzzz = 0.000001 * (zzzz * z);
                double xxxxz = 0.000001 * (xxxx * z);
                double xxxzz = 0.000001 * (xxx * zz);
                double xxxxx = 0.000001 * (xxxx * x);
                double xxzzz = 0.000001 * (xx * zzz);
                double xzzzz = 0.000001 * (xzz * zz);

                double reg_a00 = genes[0].coefficient - 10000;
                double reg_a10 = genes[1].coefficient * x;
                double reg_a01 = (genes[2].coefficient * z);
                double reg_a20 = (genes[3].coefficient * xx);
                double reg_a02 = (genes[4].coefficient * zz);
                double reg_a11 = (genes[5].coefficient * xz);
                double reg_a12 = (genes[6].coefficient * xzz);
                double reg_a30 = (genes[7].coefficient * xxx);
                double reg_a21 = (genes[8].coefficient * xxz);
                double reg_a03 = (genes[9].coefficient * zzz);
                double reg_a40 = (genes[10].coefficient * xxxx);
                double reg_a31 = (genes[11].coefficient * xxxz);
                double reg_a22 = (genes[12].coefficient * xxzz);
                double reg_a13 = (genes[13].coefficient * xzzz);
                double reg_a04 = (genes[14].coefficient * zzzz);
                double reg_a50 = (genes[15].coefficient * xxxxx);
                double reg_a41 = (genes[16].coefficient * xxxxz);
                double reg_a32 = (genes[17].coefficient * xxxzz);
                double reg_a23 = (genes[18].coefficient * xxzzz);
                double reg_a14 = (genes[19].coefficient * xzzzz);
                double reg_a05 = (genes[20].coefficient * zzzzz);

//                int64_t xx = (x*x) >> 17;
//                int64_t xxx = (xx * x) >> 17;
//                int64_t xxxx = (xx * xx) >> 17;
//                int64_t xxxxx = (xxxx * x) >> 17;
//                int64_t zz = (z*z) >> 17;
//                int64_t zzz = (zz * z) >> 17;
//                int64_t zzzz = (zz * zz) >> 17;
//                int64_t zzzzz = (zzzz * z) >> 17;
//                int64_t xz = (x*z) >> 17;
//                int64_t xxz = (xx * z) >> 17;
//                int64_t xxxz = (xxx * z) >> 17;
//                int64_t xxxxz = (xxxx * z) >> 17;
//                int64_t xzz = (x * zz) >> 17;
//                int64_t xxzz = (xx * zz) >> 17;
//                int64_t xxxzz = (xxx * zz) >> 17;
//                int64_t xzzz = (xz * zz) >> 17;
//                int64_t xxzzz = (xx * zzz) >> 17;
//                int64_t xzzzz = (xzz * zz) >> 17;

//                int64_t reg_a00 = genes[0].coefficient;
//                int64_t reg_a10 = int64_t(genes[1].coefficient * x) >> 18;
//                int64_t reg_a20 = int64_t(genes[2].coefficient * xx) >> 18;
//                int64_t reg_a30 = int64_t(genes[3].coefficient * xxx) >> 18;
//                int64_t reg_a40 = int64_t(genes[4].coefficient * xxxx) >> 18;
//                int64_t reg_a50 = int64_t(genes[5].coefficient * xxxxx) >> 18;
//                int64_t reg_a01 = int64_t(genes[6].coefficient * z) >> 18;
//                int64_t reg_a11 = int64_t(genes[7].coefficient * xz) >> 18;
//                int64_t reg_a21 = int64_t(genes[8].coefficient * xxz) >> 18;
//                int64_t reg_a31 = int64_t(genes[9].coefficient * xxxz) >> 18;
//                int64_t reg_a41 = int64_t(genes[10].coefficient * xxxxz) >> 18;
//                int64_t reg_a02 = int64_t(genes[11].coefficient * zz) >> 18;
//                int64_t reg_a12 = int64_t(genes[12].coefficient * xzz) >> 18;
//                int64_t reg_a22 = int64_t(genes[13].coefficient * xxzz) >> 18;
//                int64_t reg_a32 = int64_t(genes[14].coefficient * xxxzz) >> 18;
//                int64_t reg_a03 = int64_t(genes[15].coefficient * zzz) >> 18;
//                int64_t reg_a13 = int64_t(genes[16].coefficient * xzzz) >> 18;
//                int64_t reg_a23 = int64_t(genes[17].coefficient * xxzzz) >> 18;
//                int64_t reg_a04 = int64_t(genes[18].coefficient * zzzz) >> 18;
//                int64_t reg_a14 = int64_t(genes[19].coefficient * xzzzz) >> 18;
//                int64_t reg_a05 = int64_t(genes[20].coefficient * zzzzz) >> 18;

                double result = (reg_a00 + reg_a10 + reg_a20 + reg_a30 + reg_a40 + reg_a50
                        + reg_a01 + reg_a11 + reg_a21 + reg_a31 + reg_a41
                        + reg_a02 + reg_a12 + reg_a22 + reg_a32
                        + reg_a03 + reg_a13 + reg_a23
                        + reg_a04 + reg_a14
                        + reg_a05);


                double comp = abs(Xd[z*m_table->getWidth() + x] - result);
                if (fitness < comp)
                    fitness = comp;
            }
        }

        return fitness;
    }
};

/**
 * @brief Defines the initialization policy for the Population
 * 
 * @tparam ChromosomeT Type of Chromosome
 */
template <typename ChromosomeT>
class PopulationInitializationPolicy
{
  public:
    void initialize(std::vector<ChromosomeT> &chromosomes, const crsGA::UserData *data) const
    {
        for (auto &c : chromosomes)
        {
            for(auto &g: c.getGenes())
            {
                g.random(data);
            }
        }
    }
};

// Chromosome definition, internally a std::vector<Gen>
using Chromosome = crsGA::Chromosome<Gen, ComputeFitness>;
// Population definition
using Population = crsGA::Population<Chromosome, PopulationInitializationPolicy<Chromosome>>;
// Genetic algorithm implementation
using OneMaxGA = crsGA::GeneticAlgorithm<Gen, Chromosome, Population>;

int main(int, char **)
{
    m_table = calib_table::read_from_file("190068_XZ_2021.07.16_14.20.48.tbl");

    // Number of coefficients of the onemax bitstring
    auto numGenes = 21u;
    auto populationSize = 10u;
    auto fitnessGoal = static_cast<float>(1);
    // mutation rate
    float mutationFactor = 0.5f;
    OneMaxGA ga(populationSize, numGenes, fitnessGoal);
    ga.setMutationFactor(mutationFactor);
    ga.reset();
    // run up to 100 seconds
    ga.run(100000.0);
    return 0;
}
