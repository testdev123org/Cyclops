/*
 * BootstrapSelector.cpp
 *
 *  Created on: Nov 17, 2010
 *      Author: msuchard
 */

#include <cstdlib>
#include <iostream>
#include <sstream>

#include "BootstrapSelector.h"

namespace bsccs {

BootstrapSelector::BootstrapSelector(
		int replicates,
		std::vector<int>* inIds,
		SelectorType inType,
		long inSeed,
	    loggers::ProgressLoggerPtr _logger,
		loggers::ErrorHandlerPtr _error,		
		std::vector<real>* wtsExclude) : AbstractSelector(inIds, inType, inSeed, _logger, _error) {

    std::ostringstream stream;
	stream << "Performing bootstrap estimation with " << replicates
		<< " replicates [seed = " << seed << "]";
	logger->writeLine(stream);

	if(wtsExclude){
		for(size_t i = 0; i < wtsExclude->size(); i++){
			if(wtsExclude->at(i) == 0){
				indicesIncluded.push_back(i);
			}
		}
	}
	else{
		for(size_t i = 0; i < N; i++){
			indicesIncluded.push_back(i);
		}
	}

	permute();
}

BootstrapSelector::~BootstrapSelector() {
	// Nothing to do
}

void BootstrapSelector::permute() {
	selectedSet.clear();

	// Get non-excluded indices
	int N_new = indicesIncluded.size();
	if (type == SUBJECT) {
		for (int i = 0; i < N_new; i++) {
			int ind = rand() / (RAND_MAX / N_new + 1);
			int draw = indicesIncluded[ind];
			selectedSet.insert(draw);
		}
	} else {
        std::ostringstream stream;
        stream << "BootstrapSelector::permute is not yet implemented.";
        error->throwError(stream);	
	}

//	int total = 0;
//	for (int i = 0; i < N; i++) {
//		int count = selectedSet.count(i);
//		std::cout << i << " : " << count << std::endl;
//		total += count;
//	}
//	std::cout << "Total = " << total << std::endl;
//	exit(0);
}

void BootstrapSelector::getWeights(int batch, std::vector<real>& weights) {
	if (weights.size() != K) {
		weights.resize(K);
	}

	std::fill(weights.begin(), weights.end(), 0.0);
	if (batch == -1) {
		return;
	}

	if (type == SUBJECT) {
		for (size_t k = 0; k < K; k++) {
			int count = selectedSet.count(ids->at(k));
			weights[k] = static_cast<real>(count);
		}
	} else {
        std::ostringstream stream;
        stream << "BootstrapSelector::getWeights is not yet implemented.";
        error->throwError(stream);     
	}
}

void BootstrapSelector::getComplement(std::vector<real>& weights) {
    std::ostringstream stream;
    stream << "BootstrapSelector::getComplement is not yet implemented.";
    error->throwError(stream);
}

} // namespace

