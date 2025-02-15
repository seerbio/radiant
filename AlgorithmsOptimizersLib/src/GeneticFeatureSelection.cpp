//
// Created by andrewnichols on 2/11/25.
//

#include "GeneticFeatureSelection.h"

#include "PythiaDIAFFWorkflowAlgos/MsCalibratomaticSettertron.h"

#include <iostream>
#include <vector>
#include <algorithm>
#include <cmath>
#include <ctime>
#include <cstdlib>


namespace {

        // Generate an initial random population of feature masks
    QVector<QVector<bool>> initializePopulation(
        int populationSize,
        int numFeatures
        ) {

        QVector<QVector<bool>> population(
            populationSize,
            QVector<bool>(numFeatures, false)
            );

        for (QVector<bool> &individual : population) {
            for (int j = 0; j < numFeatures; ++j) {
                individual[j] = rand() % 2; // Randomly select or deselect features
            }
        }

        return population;
    }


    Err fitnessFunction(
        const QVector<QString> &msDataFilePaths,
        const QVector<bool> &featureMask,
        const QVector<Features> &features,
        PythiaParameters *pythiaParameters,
        TargetDecoyCandidatePairManager *targetDecoyCandidatePairManager,
        double *fitness
        ) {

        ERR_INIT

        *fitness = 0.0;

        const int selectedFeatures = std::accumulate(
            featureMask.begin(),
            featureMask.end(), 0,
            [](int a, bool b) { return a + (b ? 1 : 0); }
            );

        if (selectedFeatures == 0) {
            ERR_RETURN
        }

        for (const QString &msDataFilePath : msDataFilePaths) {

            MsReaderPointerAcc msReaderPointerAcc;
            e = msReaderPointerAcc.openFile(msDataFilePath); ree;

            TargetDecoyCandidatePairScoretron2 targetDecoyCandidatePairScoretron2;
            e = targetDecoyCandidatePairScoretron2.init(
                *pythiaParameters,
                &msReaderPointerAcc
                ); ree;

            QVector<Features> featureMaskVector;
            for (int i = 0; i < featureMask.size(); ++i) {
                if (!featureMask.at(i)) {
                    continue;
                }
                featureMaskVector.push_back(features.at(i));
            }

            MsCalibratomaticSettertron msCalibratomaticSettertron;
            e = msCalibratomaticSettertron.init(
                featureMaskVector,
                pythiaParameters,
                &msReaderPointerAcc,
                targetDecoyCandidatePairManager,
                &targetDecoyCandidatePairScoretron2,
                false
                ); ree;

            MsCalibratomatic msCalibratomatic;
            e = msCalibratomaticSettertron.buildCalibration(&msCalibratomatic);
            *fitness += msCalibratomaticSettertron.fdrWeightedMean();
        }

        ERR_RETURN
    }


    QVector<bool> crossover(
        const QVector<bool>& parent1,
        const QVector<bool>& parent2
        ) {

        int size = parent1.size();
        int crossoverPoint = rand() % size;
        QVector<bool> child(size);

        for (int i = 0; i < size; ++i) {
            if (i < crossoverPoint) {
                child[i] = parent1[i];
            } else {
                child[i] = parent2[i];
            }
        }
        return child;
    }


    void mutate(
        double mutationRate,
        QVector<bool> *individual
        ) {
        for (int i = 0; i < individual->size(); ++i) {
            const double randValue = (double)rand() / RAND_MAX;
            if (randValue < mutationRate) {
                (*individual)[i] = !(*individual)[i];
            }
        }
    }

    // Select an individual based on fitness (roulette wheel selection)
    QVector<bool> selectIndividual(
        const QVector<QVector<bool>> &population,
        const QVector<double>& fitnessValues
        ) {

        double totalFitness = 0.0;
        for (double fitness : fitnessValues) {
            totalFitness += fitness;
        }

        double randValue = ((double)rand() / RAND_MAX) * totalFitness;
        double cumulativeFitness = 0.0;

        for (int i = 0; i < population.size(); ++i) {
            cumulativeFitness += fitnessValues[i];
            if (cumulativeFitness >= randValue) {
                return population[i];
            }
        }

        return population.back(); // Fallback (should not reach)
    }

    QVector<int> qvectorBoolToQVectorInt(const QVector<bool>& boolVec) {
        QVector<int> intVec;
        for(bool b : boolVec) {
            intVec.push_back(static_cast<int>(b));
        }
        return intVec;
    }

    // Main Genetic Algorithm for Feature Selection
    QVector<bool> geneticAlgorithm(
        const QVector<QString> &msDataFilePaths,
        const QVector<Features> &features,
        int populationSize,
        int numGenerations,
        double mutationRate,
        PythiaParameters *pythiaParameters,
        TargetDecoyCandidatePairManager *targetDecoyCandidatePairManager
        ) {

        const int numFeatures = features.size();
        QVector<QVector<bool>> population = initializePopulation(
            populationSize,
            numFeatures
            );

        QVector<int> featureFrequencies(features.size(), 0);
        double bestFitnessOverall = 0.0;

        for (int generation = 0; generation < numGenerations; ++generation) {

            QVector<double> fitnessValues(populationSize);
            for (int i = 0; i < populationSize; ++i) {
                    fitnessFunction(
                        msDataFilePaths,
                        population[i],
                        features,
                        pythiaParameters,
                        targetDecoyCandidatePairManager,
                        &fitnessValues[i]
                        );

                qDebug() << "CURRENT POPULATION:" << qvectorBoolToQVectorInt(population[i]);
                qDebug() << "CURRENT FITNESS: " << fitnessValues;
            }

            QVector<QVector<bool>> newPopulation;
            for (int i = 0; i < populationSize; ++i) {

                const QVector<bool> parent1 = selectIndividual(population, fitnessValues);
                const QVector<bool> parent2 = selectIndividual(population, fitnessValues);

                QVector<bool> child = crossover(parent1, parent2);
                mutate(mutationRate, &child);
                newPopulation.push_back(child);
            }

            population = newPopulation;

            // Print best fitness in the current generation
            const double bestFitness = *std::max_element(fitnessValues.begin(), fitnessValues.end());
            bestFitnessOverall = std::max(bestFitnessOverall, bestFitness);
            std::cout << "Generation " << generation + 1 << ": Best Fitness = " << bestFitness << std::endl;

            for (int i = 0; i < fitnessValues.size(); ++i) {
                if (fitnessValues[i] >= bestFitnessOverall) {
                    for (int j = 0; j < population[i].size(); ++j) {

                        if (!population[i][j]) {
                            continue;
                        }

                        featureFrequencies[j]++;
                    }
                }
            }

            const int featureFrequenceMax = *std::max_element(featureFrequencies.begin(), featureFrequencies.end());

            QVector<double> featureFrequenciesNormalized;
            std::transform(
                featureFrequencies.begin(),
                featureFrequencies.end(),
                std::back_inserter(featureFrequenciesNormalized),
                [featureFrequenceMax](int i) { return MathUtils::pRound((double)i / featureFrequenceMax, 3); }
                );
            qDebug() << "Feature Frequencies: " << featureFrequencies;
            qDebug() << "Feature Frequencies Normalized: " << featureFrequenciesNormalized;
        }

        // Return the best individual from the final generation
        QVector<double> finalFitnessValues(populationSize);
        for (int i = 0; i < populationSize; ++i) {
            finalFitnessValues[i] = fitnessFunction(
                msDataFilePaths,
                population[i],
                features,
                pythiaParameters,
                targetDecoyCandidatePairManager,
                &finalFitnessValues[i]
                );
            qDebug() << "CURRENT POPULATION:" << qvectorBoolToQVectorInt(population[i]);
            qDebug() << "CURRENT FITNESS: " << finalFitnessValues;
        }
        const int bestIndex = std::max_element(finalFitnessValues.begin(), finalFitnessValues.end()) - finalFitnessValues.begin();
        return population[bestIndex];
    }

}//namespace



Err GeneticFeatureSelection::run(
    const PythiaParameters &pythiaParameters,
    const QVector<QString> &msDataFilePaths,
    const QVector<Features> &features,
    int maxGenerationIterations,
    int populationSize,
    double mutationRate,
    TargetDecoyCandidatePairManager *targetDecoyCandidatePairManager
    ) {

    ERR_INIT

    e = ErrorUtils::isNotEmpty(features); ree;
    e = ErrorUtils::isNotEmpty(msDataFilePaths); ree;

    srand(S_GLOBAL_SETTINGS.NUMBER_OF_THE_BEAST);

    PythiaParameters pythiaParametersCopy = pythiaParameters;

    const QVector<bool> bestFeatures = geneticAlgorithm(
        msDataFilePaths,
        features,
        populationSize,
        maxGenerationIterations,
        mutationRate,
        &pythiaParametersCopy,
        targetDecoyCandidatePairManager
        );

    ERR_RETURN
}
