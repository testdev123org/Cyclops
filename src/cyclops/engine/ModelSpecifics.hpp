/*
 * ModelSpecifics.hpp
 *
 *  Created on: Jul 13, 2012
 *      Author: msuchard
 */

#ifndef MODELSPECIFICS_HPP_
#define MODELSPECIFICS_HPP_

#include <cmath>
#include <cstdlib>
#include <algorithm>
#include <numeric>

#include "ModelSpecifics.h"
#include "Iterators.h"

#include "Recursions.hpp"

//#define USE_BIGNUM
#define USE_LONG_DOUBLE

namespace bsccs {

#ifdef USE_BIGNUM
	typedef bigNum DDouble;
#else
	#ifdef USE_LONG_DOUBLE
		typedef long double DDouble;
	#else
		typedef double DDouble;
	#endif
#endif
    
#ifdef DEBUG_COX
    using std::cerr;
    using std::endl;
#endif

//template <class BaseModel,typename WeightType>
//ModelSpecifics<BaseModel,WeightType>::ModelSpecifics(
//		const std::vector<real>& y,
//		const std::vector<real>& z) : AbstractModelSpecifics(y, z), BaseModel() {
//	// TODO Memory allocation here
//}

template <class BaseModel,typename WeightType>
ModelSpecifics<BaseModel,WeightType>::ModelSpecifics(const ModelData& input)
	: AbstractModelSpecifics(input), BaseModel() {
	// TODO Memory allocation here
}

template <class BaseModel,typename WeightType>
ModelSpecifics<BaseModel,WeightType>::~ModelSpecifics() {
	// TODO Memory release here
}

template <class BaseModel,typename WeightType>
bool ModelSpecifics<BaseModel,WeightType>::allocateXjY(void) { return BaseModel::precomputeGradient; }

template <class BaseModel,typename WeightType>
bool ModelSpecifics<BaseModel,WeightType>::allocateXjX(void) { return BaseModel::precomputeHessian; }

template <class BaseModel,typename WeightType>
bool ModelSpecifics<BaseModel,WeightType>::sortPid(void) { return BaseModel::sortPid; }

template <class BaseModel,typename WeightType>
bool ModelSpecifics<BaseModel,WeightType>::initializeAccumulationVectors(void) { return BaseModel::cumulativeGradientAndHessian; }

template <class BaseModel,typename WeightType>
bool ModelSpecifics<BaseModel,WeightType>::allocateNtoKIndices(void) { return BaseModel::hasNtoKIndices; }

template <class BaseModel,typename WeightType>
bool ModelSpecifics<BaseModel,WeightType>::hasResetableAccumulators(void) { return BaseModel::hasResetableAccumulators; }
   
template <class BaseModel,typename WeightType>
void ModelSpecifics<BaseModel,WeightType>::setWeights(real* inWeights, bool useCrossValidation) {
	// Set K weights
	if (hKWeight.size() != K) {
		hKWeight.resize(K);
	}
	if (useCrossValidation) {
		for (size_t k = 0; k < K; ++k) {
			hKWeight[k] = inWeights[k];
		}
	} else {
		std::fill(hKWeight.begin(), hKWeight.end(), static_cast<WeightType>(1));
	}
	// Set N weights (these are the same for independent data models
	if (hNWeight.size() != N) {
		hNWeight.resize(N);
	}
	std::fill(hNWeight.begin(), hNWeight.end(), static_cast<WeightType>(0));
	for (size_t k = 0; k < K; ++k) {
		WeightType event = BaseModel::observationCount(hY[k])*hKWeight[k];
		incrementByGroup(hNWeight.data(), hPid, k, event);
	}
}

template<class BaseModel, typename WeightType>
void ModelSpecifics<BaseModel, WeightType>::computeXjY(bool useCrossValidation) {
	for (size_t j = 0; j < J; ++j) {
		hXjY[j] = 0;
				
		GenericIterator it(*hXI, j);

		if (useCrossValidation) {
			for (; it; ++it) {
				const int k = it.index();				
				if (BaseModel::exactTies && hNWeight[BaseModel::getGroup(hPid, k)] > 1) {
					// Do not precompute
				} else {
					hXjY[j] += it.value() * hY[k] * hKWeight[k];
				}
			}
		} else {
			for (; it; ++it) {
				const int k = it.index();
				if (BaseModel::exactTies && hNWeight[BaseModel::getGroup(hPid, k)] > 1) {
					// Do not precompute					
				} else {							
					hXjY[j] += it.value() * hY[k];
				}
			}
		}
#ifdef DEBUG_COX
		cerr << "j: " << j << " = " << hXjY[j]<< endl;
#endif
	}
}

template<class BaseModel, typename WeightType>
void ModelSpecifics<BaseModel, WeightType>::computeXjX(bool useCrossValidation) {
	for (size_t j = 0; j < J; ++j) {
		hXjX[j] = 0;
		GenericIterator it(*hXI, j);

		if (useCrossValidation) {
			for (; it; ++it) {
				const int k = it.index();
				if (BaseModel::exactTies && hNWeight[BaseModel::getGroup(hPid, k)] > 1) {
					// Do not precompute
				} else {				
					hXjX[j] += it.value() * it.value() * hKWeight[k];
				}
			}
		} else {
			for (; it; ++it) {
				const int k = it.index();
				if (BaseModel::exactTies && hNWeight[BaseModel::getGroup(hPid, k)] > 1) {
					// Do not precompute
				} else {								
					hXjX[j] += it.value() * it.value();
				}
			}
		}
	}
}

template<class BaseModel, typename WeightType>
void ModelSpecifics<BaseModel, WeightType>::computeNtoKIndices(bool useCrossValidation) {

	hNtoK.resize(N+1);
	int n = 0;
	for (int k = 0; k < K;) {
		hNtoK[n] = k;
		int currentPid = hPid[k];
		do {
			++k;
		} while (k < K && currentPid == hPid[k]);
		++n;
	}
	hNtoK[n] = K;
}

template <class BaseModel,typename WeightType>
void ModelSpecifics<BaseModel,WeightType>::computeFixedTermsInLogLikelihood(bool useCrossValidation) {
	if(BaseModel::likelihoodHasFixedTerms) {
		logLikelihoodFixedTerm = 0.0;
		if(useCrossValidation) {
			for(size_t i = 0; i < K; i++) {
				logLikelihoodFixedTerm += BaseModel::logLikeFixedTermsContrib(hY[i], hOffs[i], hOffs[i]) * hKWeight[i];
			}
		} else {
			for(size_t i = 0; i < K; i++) {
				logLikelihoodFixedTerm += BaseModel::logLikeFixedTermsContrib(hY[i], hOffs[i], hOffs[i]);
			}
		}
	}
}

template <class BaseModel,typename WeightType>
void ModelSpecifics<BaseModel,WeightType>::computeFixedTermsInGradientAndHessian(bool useCrossValidation) {
	if (sortPid()) {
		doSortPid(useCrossValidation);
	}
	if (allocateXjY()) {
		computeXjY(useCrossValidation);
	}
	if (allocateXjX()) {
		computeXjX(useCrossValidation);
	}
	if (allocateNtoKIndices()) {
		computeNtoKIndices(useCrossValidation);
	}
}

template <class BaseModel,typename WeightType>
double ModelSpecifics<BaseModel,WeightType>::getLogLikelihood(bool useCrossValidation) {

	real logLikelihood = static_cast<real>(0.0);
	if (useCrossValidation) {
		for (size_t i = 0; i < K; i++) {
			logLikelihood += BaseModel::logLikeNumeratorContrib(hY[i], hXBeta[i]) * hKWeight[i];
		}
	} else {
		for (size_t i = 0; i < K; i++) {
			logLikelihood += BaseModel::logLikeNumeratorContrib(hY[i], hXBeta[i]);
		}
	}

	if (BaseModel::likelihoodHasDenominator) { // Compile-time switch
		if(BaseModel::cumulativeGradientAndHessian) {
			for (size_t i = 0; i < N; i++) {
				// Weights modified in computeNEvents()
				logLikelihood -= BaseModel::logLikeDenominatorContrib(hNWeight[i], accDenomPid[i]);
			}
		} else {  // TODO Unnecessary code duplication
			for (size_t i = 0; i < N; i++) {
				// Weights modified in computeNEvents()
				logLikelihood -= BaseModel::logLikeDenominatorContrib(hNWeight[i], denomPid[i]);
			}
		}
	}

	if (BaseModel::likelihoodHasFixedTerms) {
		logLikelihood += logLikelihoodFixedTerm;
	}

	return static_cast<double>(logLikelihood);
}

template <class BaseModel,typename WeightType>
double ModelSpecifics<BaseModel,WeightType>::getPredictiveLogLikelihood(real* weights) {
	real logLikelihood = static_cast<real>(0.0);

	if(BaseModel::cumulativeGradientAndHessian)	{
		for (size_t i = 0; i < N; ++i) {
			logLikelihood += BaseModel::logPredLikeContrib(hY[i], weights[i], hXBeta[i], &accDenomPid[0], hPid, i); // TODO Going to crash with ties
		}
	} else { // TODO Unnecessary code duplication
		for (size_t k = 0; k < K; ++k) { // TODO Is index of K correct?
			logLikelihood += BaseModel::logPredLikeContrib(hY[k], weights[k], hXBeta[k], denomPid, hPid, k);
		}
	}
	return static_cast<double>(logLikelihood);
}

template <class BaseModel,typename WeightType>
void ModelSpecifics<BaseModel,WeightType>::getPredictiveEstimates(real* y, real* weights){

	// TODO Check with SM: the following code appears to recompute hXBeta at large expense
//	std::vector<real> xBeta(K,0.0);
//	for(int j = 0; j < J; j++){
//		GenericIterator it(*hXI, j);
//		for(; it; ++it){
//			const int k = it.index();
//			xBeta[k] += it.value() * hBeta[j] * weights[k];
//		}
//	}
	if (weights) {
		for (size_t k = 0; k < K; ++k) {
			if (weights[k]) {
				BaseModel::predictEstimate(y[k], hXBeta[k]);
			}
		}
	} else {
		for (size_t k = 0; k < K; ++k) {
			BaseModel::predictEstimate(y[k], hXBeta[k]);
		}
	}
	// TODO How to remove code duplication above?
}

// TODO The following function is an example of a double-dispatch, rewrite without need for virtual function
template <class BaseModel,typename WeightType>
void ModelSpecifics<BaseModel,WeightType>::computeGradientAndHessian(int index, double *ogradient,
		double *ohessian, bool useWeights) {
	// Run-time dispatch, so virtual call should not effect speed
	if (useWeights) {
		switch (hXI->getFormatType(index)) {
			case INDICATOR :
				computeGradientAndHessianImpl<IndicatorIterator>(index, ogradient, ohessian, weighted);
				break;
			case SPARSE :
				computeGradientAndHessianImpl<SparseIterator>(index, ogradient, ohessian, weighted);
				break;
			case DENSE :
				computeGradientAndHessianImpl<DenseIterator>(index, ogradient, ohessian, weighted);
				break;
			case INTERCEPT :
				computeGradientAndHessianImpl<InterceptIterator>(index, ogradient, ohessian, weighted);
				break;
		}
	} else {
		switch (hXI->getFormatType(index)) {
			case INDICATOR :
				computeGradientAndHessianImpl<IndicatorIterator>(index, ogradient, ohessian, unweighted);
				break;
			case SPARSE :
				computeGradientAndHessianImpl<SparseIterator>(index, ogradient, ohessian, unweighted);
				break;
			case DENSE :
				computeGradientAndHessianImpl<DenseIterator>(index, ogradient, ohessian, unweighted);
				break;
			case INTERCEPT :
				computeGradientAndHessianImpl<InterceptIterator>(index, ogradient, ohessian, unweighted);
				break;
		}
	}
}

template <class BaseModel,typename WeightType> template <class IteratorType, class Weights>
void ModelSpecifics<BaseModel,WeightType>::computeGradientAndHessianImpl(int index, double *ogradient,
		double *ohessian, Weights w) {
	real gradient = static_cast<real>(0);
	real hessian = static_cast<real>(0);

	IteratorType it(*(sparseIndices)[index], N); // TODO How to create with different constructor signatures?

//std::cout << "YOYOYO" << std::endl;

	if (BaseModel::cumulativeGradientAndHessian) { // Compile-time switch
		
		real accNumerPid  = static_cast<real>(0);
		real accNumerPid2 = static_cast<real>(0);

		// This is an optimization point compared to iterating over a completely dense view:  
		// a) the view below starts at the first non-zero entry
		// b) we only access numerPid and numerPid2 for non-zero entries 
		// This may save time; should document speed-up in massive Cox manuscript
		
        // find start relavent accumulator reset point
        auto reset = begin(accReset);
        while( *reset < it.index() ) {
            ++reset;
        }
        
//         std::cout << "Will reset at " << *reset << std::endl;
				
		for (; it; ) {
			int i = it.index();
// 			std::cout << "i = " << i << std::endl;
			
// TODO CHECK		
			if (*reset <= i) {
			    accNumerPid  = static_cast<real>(0.0);
			    accNumerPid2 = static_cast<real>(0.0);
			    ++reset;
			} 			    
						
			if(w.isWeighted){ //if useCrossValidation
				accNumerPid  += numerPid[i] * hNWeight[i]; // TODO Only works when X-rows are sorted as well
				accNumerPid2 += numerPid2[i] * hNWeight[i];
			} else { // TODO Unnecessary code duplication
				accNumerPid  += numerPid[i]; // TODO Only works when X-rows are sorted as well
				accNumerPid2 += numerPid2[i];
			}
#ifdef DEBUG_COX
			cerr << "w: " << i << " " << hNWeight[i] << " " << numerPid[i] << ":" <<
					accNumerPid << ":" << accNumerPid2 << ":" << accDenomPid[i];
#endif			
			// Compile-time delegation
			BaseModel::incrementGradientAndHessian(it,
					w, // Signature-only, for iterator-type specialization
					&gradient, &hessian, accNumerPid, accNumerPid2,
					accDenomPid[i], hNWeight[i], it.value(), hXBeta[i], hY[i]); // When function is in-lined, compiler will only use necessary arguments
#ifdef DEBUG_COX		
			cerr << " -> g:" << gradient << " h:" << hessian << endl;	
#endif
			++it;
			
			if (IteratorType::isSparse) {
				const int next = it ? it.index() : N;
				for (++i; i < next; ++i) {
#ifdef DEBUG_COX
			cerr << "q: " << i << " " << hNWeight[i] << " " << 0 << ":" <<
					accNumerPid << ":" << accNumerPid2 << ":" << accDenomPid[i];
#endif	
// TODO CHECK                   
                    if (*reset <= i) {
			            accNumerPid  = static_cast<real>(0.0);
        			    accNumerPid2 = static_cast<real>(0.0);
		        	    ++reset;                   
                   } 		
					
					BaseModel::incrementGradientAndHessian(it,
							w, // Signature-only, for iterator-type specialization
							&gradient, &hessian, accNumerPid, accNumerPid2,
							accDenomPid[i], hNWeight[i], static_cast<real>(0), hXBeta[i], hY[i]); // When function is in-lined, compiler will only use necessary arguments
#ifdef DEBUG_COX		
			cerr << " -> g:" << gradient << " h:" << hessian << endl;	
#endif
					
				}						
			}
		}
	} else {
		for (; it; ++it) {
			const int i = it.index();
			
			if (BaseModel::exactTies && hNWeight[i] > 1) {
				int numSubjects = hNtoK[i+1] - hNtoK[i];
				int numCases = hNWeight[i];
				DenseView<IteratorType> x(IteratorType(*hXI, index), hNtoK[i], hNtoK[i+1]);
				
//				std::cerr << "Here " << hNtoK.size() << std::endl;			

				std::vector<DDouble> value = computeHowardRecursion<DDouble>(offsExpXBeta.begin() + hNtoK[i], x, numSubjects, numCases, hY + hNtoK[i]);

// 				std::cerr << i << std::endl;
// 				std::for_each(begin(value), end(value), [](DDouble d) {
// 					std::cerr << d << std::endl;
// 				});
// //				exit(-1);
// 				std::cerr << std::endl;
// 				
// 				if (index == 1) {
// 					for (; x; ++x) {
// 						std::cerr << " " << x.value();
// 					}
// 					std::cerr << std::endl;
// 					auto it = offsExpXBeta.begin() + hNtoK[i];
// 					for (int ii = 0; ii < numCases; ++ii) {					
// 						std::cerr << " " << *(it + ii);
// 					}
// 					std::cerr << std::endl;
// 					//exit(-1);
// 				}


				//MM
				//real banana = (real)((value[1]*value[1])/(value[0]*value[0]) - value[4]*value[4]);
				//real banana = (real)(pow(bigNum::div(value[1],value[0]).toDouble(), 2) - pow(value[4].toDouble(),2));
		
				//normal
				using namespace sugar;
								
// 				gradient += (real)(value[3] - value[1]/value[0]);				
// 				hessian += (real)((value[1]/value[0]) * (value[1]/value[0]) - value[2]/value[0]);				
				gradient += (real)(value[1]/value[0] - value[3]);	
				hessian += (real)(value[2]/value[0] - (value[1]/value[0]) * (value[1]/value[0]));	
			} else {			
				// Compile-time delegation
				BaseModel::incrementGradientAndHessian(it,
						w, // Signature-only, for iterator-type specialization
						&gradient, &hessian, numerPid[i], numerPid2[i],
						denomPid[i], hNWeight[i], it.value(), hXBeta[i], hY[i]); // When function is in-lined, compiler will only use necessary arguments		
			}
		}
	}
  // TODO Figure out how to handle these ...  do NOT pre-compute for ties????
	if (BaseModel::precomputeGradient) { // Compile-time switch
		gradient -= hXjY[index];
	}

	if (BaseModel::precomputeHessian) { // Compile-time switch
		hessian += static_cast<real>(2.0) * hXjX[index];
	}

	*ogradient = static_cast<double>(gradient);
	*ohessian = static_cast<double>(hessian);
}

template <class BaseModel,typename WeightType>
void ModelSpecifics<BaseModel,WeightType>::computeFisherInformation(int indexOne, int indexTwo,
		double *oinfo, bool useWeights) {

	if (useWeights) {
// 		std::cerr << "Weights are not yet implemented in Fisher Information calculations" << std::endl;
// 		exit(-1);
		throw new std::logic_error("Weights are not yet implemented in Fisher Information calculations");
	} else { // no weights
		switch (hXI->getFormatType(indexOne)) {
			case INDICATOR :
				dispatchFisherInformation<IndicatorIterator>(indexOne, indexTwo, oinfo, weighted);
				break;
			case SPARSE :
				dispatchFisherInformation<SparseIterator>(indexOne, indexTwo, oinfo, weighted);
				break;
			case DENSE :
				dispatchFisherInformation<DenseIterator>(indexOne, indexTwo, oinfo, weighted);
				break;
			case INTERCEPT :
				dispatchFisherInformation<InterceptIterator>(indexOne, indexTwo, oinfo, weighted);
				break;
		}
	}
}

template <class BaseModel, typename WeightType> template <typename IteratorTypeOne, class Weights>
void ModelSpecifics<BaseModel,WeightType>::dispatchFisherInformation(int indexOne, int indexTwo, double *oinfo, Weights w) {
	switch (hXI->getFormatType(indexTwo)) {
		case INDICATOR :
			computeFisherInformationImpl<IteratorTypeOne,IndicatorIterator>(indexOne, indexTwo, oinfo, w);
			break;
		case SPARSE :
			computeFisherInformationImpl<IteratorTypeOne,SparseIterator>(indexOne, indexTwo, oinfo, w);
			break;
		case DENSE :
			computeFisherInformationImpl<IteratorTypeOne,DenseIterator>(indexOne, indexTwo, oinfo, w);
			break;
		case INTERCEPT :
			computeFisherInformationImpl<IteratorTypeOne,InterceptIterator>(indexOne, indexTwo, oinfo, w);
			break;
	}
//	std::cerr << "End of dispatch" << std::endl;
}


template<class BaseModel, typename WeightType> template<class IteratorType>
SparseIterator ModelSpecifics<BaseModel, WeightType>::getSubjectSpecificHessianIterator(int index) {

	if (hessianSparseCrossTerms.find(index) == hessianSparseCrossTerms.end()) {
		// Make new
//		std::vector<int>* indices = new std::vector<int>();
        auto indices = make_shared<std::vector<int> >();
//		std::vector<real>* values = new std::vector<real>();
        auto values = make_shared<std::vector<real> >();      
//		CompressedDataColumn* column = new CompressedDataColumn(indices, values,
//				SPARSE);
    	CDCPtr column = bsccs::make_shared<CompressedDataColumn>(indices, values, SPARSE);
		hessianSparseCrossTerms.insert(std::make_pair(index, 
// 		    CompressedDataColumn(indices, values, SPARSE)));
		    column));

		IteratorType itCross(*hXI, index);
		for (; itCross;) {
			real value = 0.0;
			int currentPid = hPid[itCross.index()];  // TODO Need to fix for stratified Cox
			do {
				const int k = itCross.index();
				value += BaseModel::gradientNumeratorContrib(itCross.value(),
						offsExpXBeta[k], hXBeta[k], hY[k]);
				++itCross;
			} while (itCross && currentPid == hPid[itCross.index()]); // TODO Need to fix for stratified Cox
			indices->push_back(currentPid);
			values->push_back(value);
		}
	}
	return SparseIterator(*hessianSparseCrossTerms[index]);

}


template <class BaseModel, typename WeightType> template <class IteratorTypeOne, class IteratorTypeTwo, class Weights>
void ModelSpecifics<BaseModel,WeightType>::computeFisherInformationImpl(int indexOne, int indexTwo, double *oinfo, Weights w) {

	IteratorTypeOne itOne(*hXI, indexOne);
	IteratorTypeTwo itTwo(*hXI, indexTwo);
	PairProductIterator<IteratorTypeOne,IteratorTypeTwo> it(itOne, itTwo);

	real information = static_cast<real>(0);
	for (; it.valid(); ++it) {
		const int k = it.index();
		// Compile-time delegation

		BaseModel::incrementFisherInformation(it,
				w, // Signature-only, for iterator-type specialization
				&information,
				offsExpXBeta[k],
				0.0, 0.0, // numerPid[k], numerPid2[k], // remove
				denomPid[BaseModel::getGroup(hPid, k)],
				hKWeight[k], it.value(), hXBeta[k], hY[k]); // When function is in-lined, compiler will only use necessary arguments
	}

	if (BaseModel::hasStrataCrossTerms) {

		// Check if index is pre-computed
//#define USE_DENSE
#ifdef USE_DENSE
		if (hessianCrossTerms.find(indexOne) == hessianCrossTerms.end()) {
			// Make new
			std::vector<real> crossOneTerms(N);
			IteratorTypeOne crossOne(*hXI, indexOne);
			for (; crossOne; ++crossOne) {
				const int k = crossOne.index();
				incrementByGroup(crossOneTerms.data(), hPid, k,
						BaseModel::gradientNumeratorContrib(crossOne.value(), offsExpXBeta[k], hXBeta[k], hY[k]));
			}
			hessianCrossTerms[indexOne];
//			std::cerr << std::accumulate(crossOneTerms.begin(), crossOneTerms.end(), 0.0) << std::endl;
			hessianCrossTerms[indexOne].swap(crossOneTerms);
		}
		std::vector<real>& crossOneTerms = hessianCrossTerms[indexOne];

		// TODO Remove code duplication
		if (hessianCrossTerms.find(indexTwo) == hessianCrossTerms.end()) {
			std::vector<real> crossTwoTerms(N);
			IteratorTypeTwo crossTwo(*hXI, indexTwo);
			for (; crossTwo; ++crossTwo) {
				const int k = crossTwo.index();
				incrementByGroup(crossTwoTerms.data(), hPid, k,
						BaseModel::gradientNumeratorContrib(crossTwo.value(), offsExpXBeta[k], hXBeta[k], hY[k]));
			}
			hessianCrossTerms[indexTwo];
//			std::cerr << std::accumulate(crossTwoTerms.begin(), crossTwoTerms.end(), 0.0) << std::endl;
			hessianCrossTerms[indexTwo].swap(crossTwoTerms);
		}
		std::vector<real>& crossTwoTerms = hessianCrossTerms[indexTwo];

		// TODO Sparse loop
		real cross = 0.0;
		for (int n = 0; n < N; ++n) {
			cross += crossOneTerms[n] * crossTwoTerms[n] / (denomPid[n] * denomPid[n]);
		}
//		std::cerr << cross << std::endl;
		information -= cross;
#else
		SparseIterator sparseCrossOneTerms = getSubjectSpecificHessianIterator<IteratorTypeOne>(indexOne);
		SparseIterator sparseCrossTwoTerms = getSubjectSpecificHessianIterator<IteratorTypeTwo>(indexTwo);
		PairProductIterator<SparseIterator,SparseIterator> itSparseCross(sparseCrossOneTerms, sparseCrossTwoTerms);

		real sparseCross = 0.0;
		for (; itSparseCross.valid(); ++itSparseCross) {
			const int n = itSparseCross.index();
			sparseCross += itSparseCross.value() / (denomPid[n] * denomPid[n]);
		}
		information -= sparseCross;
#endif
	}

	*oinfo = static_cast<double>(information);
}



template <class BaseModel,typename WeightType>
void ModelSpecifics<BaseModel,WeightType>::computeNumeratorForGradient(int index) {
	// Run-time delegation
	switch (hXI->getFormatType(index)) {
		case INDICATOR : {
			IndicatorIterator it(*(sparseIndices)[index]);
			for (; it; ++it) { // Only affected entries
				numerPid[it.index()] = static_cast<real>(0.0);
			}
			incrementNumeratorForGradientImpl<IndicatorIterator>(index);
			}
			break;
		case SPARSE : {
			SparseIterator it(*(sparseIndices)[index]);
			for (; it; ++it) { // Only affected entries
				numerPid[it.index()] = static_cast<real>(0.0);
				if (BaseModel::hasTwoNumeratorTerms) { // Compile-time switch
					numerPid2[it.index()] = static_cast<real>(0.0); // TODO Does this invalid the cache line too much?
				}
			}
			incrementNumeratorForGradientImpl<SparseIterator>(index); }
			break;
		case DENSE :
			zeroVector(numerPid, N);
			if (BaseModel::hasTwoNumeratorTerms) { // Compile-time switch
				zeroVector(numerPid2, N);
			}
			incrementNumeratorForGradientImpl<DenseIterator>(index);
			break;
		case INTERCEPT :
			zeroVector(numerPid, N);
			if (BaseModel::hasTwoNumeratorTerms) { // Compile-time switch
				zeroVector(numerPid2, N);
			}
			incrementNumeratorForGradientImpl<InterceptIterator>(index);
			break;
		default : break;
			// throw error
			//exit(-1);
	}	
}

template <class BaseModel,typename WeightType> template <class IteratorType>
void ModelSpecifics<BaseModel,WeightType>::incrementNumeratorForGradientImpl(int index) {
	IteratorType it(*hXI, index);
	for (; it; ++it) {
		const int k = it.index();
		incrementByGroup(numerPid, hPid, k,
				BaseModel::gradientNumeratorContrib(it.value(), offsExpXBeta[k], hXBeta[k], hY[k]));
		if (!IteratorType::isIndicator && BaseModel::hasTwoNumeratorTerms) {
			incrementByGroup(numerPid2, hPid, k,
					BaseModel::gradientNumerator2Contrib(it.value(), offsExpXBeta[k]));
		}
		
#ifdef DEBUG_COX			
//			if (numerPid[BaseModel::getGroup(hPid, k)] > 0 && numerPid[BaseModel::getGroup(hPid, k)] < 1e-40) {
				cerr << "Increment" << endl;
				cerr << "hPid = " << hPid << ", k = " << k << ", index = " << BaseModel::getGroup(hPid, k) << endl;
				cerr << BaseModel::gradientNumeratorContrib(it.value(), offsExpXBeta[k], hXBeta[k], hY[k]) <<  " "
				<< it.value() << " " << offsExpXBeta[k] << " " << hXBeta[k] << " " << hY[k] << endl;
//				exit(-1);
//			}
#endif		
		
		
		
	}
}

template <class BaseModel,typename WeightType>
void ModelSpecifics<BaseModel,WeightType>::updateXBeta(real realDelta, int index, bool useWeights) {
	// Run-time dispatch to implementation depending on covariate FormatType
	switch(hXI->getFormatType(index)) {
		case INDICATOR :
			updateXBetaImpl<IndicatorIterator>(realDelta, index, useWeights);
			break;
		case SPARSE :
			updateXBetaImpl<SparseIterator>(realDelta, index, useWeights);
			break;
		case DENSE :
			updateXBetaImpl<DenseIterator>(realDelta, index, useWeights);
			break;
		case INTERCEPT :
			updateXBetaImpl<InterceptIterator>(realDelta, index, useWeights);
			break;
		default : break;
			// throw error
			//exit(-1);
	}	
}

template <class BaseModel,typename WeightType> template <class IteratorType>
inline void ModelSpecifics<BaseModel,WeightType>::updateXBetaImpl(real realDelta, int index, bool useWeights) {
	IteratorType it(*hXI, index);
	for (; it; ++it) {
		const int k = it.index();
		hXBeta[k] += realDelta * it.value(); // TODO Check optimization with indicator and intercept
		// Update denominators as well
		if (BaseModel::likelihoodHasDenominator) { // Compile-time switch
			real oldEntry = offsExpXBeta[k];
			real newEntry = offsExpXBeta[k] = BaseModel::getOffsExpXBeta(hOffs, hXBeta[k], hY[k], k);
			incrementByGroup(denomPid, hPid, k, (newEntry - oldEntry));
		}
	}
	computeAccumlatedNumerDenom(useWeights);
}

template <class BaseModel,typename WeightType>
void ModelSpecifics<BaseModel,WeightType>::computeRemainingStatistics(bool useWeights) {
	if (BaseModel::likelihoodHasDenominator) {
		fillVector(denomPid, N, BaseModel::getDenomNullValue());
		for (size_t k = 0; k < K; ++k) {
			offsExpXBeta[k] = BaseModel::getOffsExpXBeta(hOffs, hXBeta[k], hY[k], k);
			incrementByGroup(denomPid, hPid, k, offsExpXBeta[k]);
		}
		computeAccumlatedNumerDenom(useWeights);
	} 
#ifdef DEBUG_COX
	cerr << "Done with initial denominators" << endl;

	for (int i = 0; i < N; ++i) {
		cerr << denomPid[i] << " " << accDenomPid[i] << " " << numerPid[i] << endl;
	}
#endif
}

template <class BaseModel,typename WeightType>
void ModelSpecifics<BaseModel,WeightType>::computeAccumlatedNumerDenom(bool useWeights) {

	if (BaseModel::likelihoodHasDenominator && //The two switches should ideally be separated
		BaseModel::cumulativeGradientAndHessian) { // Compile-time switch
			if (accDenomPid.size() != N) {
				accDenomPid.resize(N, static_cast<real>(0));
			}
			if (accNumerPid.size() != N) {
				accNumerPid.resize(N, static_cast<real>(0));
			}
			if (accNumerPid2.size() != N) {
				accNumerPid2.resize(N, static_cast<real>(0));
			}

			// prefix-scan
			if(useWeights) { 
				//accumulating separately over train and validation sets
				real totalDenomTrain = static_cast<real>(0);
				real totalNumerTrain = static_cast<real>(0);
				real totalNumer2Train = static_cast<real>(0);
				real totalDenomValid = static_cast<real>(0);
				real totalNumerValid = static_cast<real>(0);
				real totalNumer2Valid = static_cast<real>(0);
				
// TODO CHECK   
                auto reset = begin(accReset);
				
				for (size_t i = 0; i < N; ++i) {
// TODO CHECK				
                    if (static_cast<unsigned int>(*reset) == i) { // TODO Check with sparse
			            totalDenomTrain = static_cast<real>(0);
				        totalNumerTrain = static_cast<real>(0);
				        totalNumer2Train = static_cast<real>(0);
				        totalDenomValid = static_cast<real>(0);
				        totalNumerValid = static_cast<real>(0);
				        totalNumer2Valid = static_cast<real>(0);				        
				        ++reset;
 				    }
 				    
					if(hKWeight[i] == 1.0){
						totalDenomTrain += denomPid[i];
						totalNumerTrain += numerPid[i];
						totalNumer2Train += numerPid2[i];
						accDenomPid[i] = totalDenomTrain;
						accNumerPid[i] = totalNumerTrain;
						accNumerPid2[i] = totalNumer2Train;
					} else {
						totalDenomValid += denomPid[i];
						totalNumerValid += numerPid[i];
						totalNumer2Valid += numerPid2[i];
						accDenomPid[i] = totalDenomValid;
						accNumerPid[i] = totalNumerValid;
						accNumerPid2[i] = totalNumer2Valid;
					}
				}
			} else {
				real totalDenom = static_cast<real>(0);
				real totalNumer = static_cast<real>(0);
				real totalNumer2 = static_cast<real>(0);
				
				auto reset = begin(accReset);
				
				for (size_t i = 0; i < N; ++i) {
// TODO CHECK				
                    if (static_cast<unsigned int>(*reset) == i) { // TODO Check with SPARSE
				        totalDenom = static_cast<real>(0);
				        totalNumer = static_cast<real>(0);
				        totalNumer2 = static_cast<real>(0);				    
                        ++reset;				    
				    }
				    
					totalDenom += denomPid[i];
					totalNumer += numerPid[i];
					totalNumer2 += numerPid2[i];
					accDenomPid[i] = totalDenom;
					accNumerPid[i] = totalNumer;
					accNumerPid2[i] = totalNumer2;
#ifdef DEBUG_COX
					cerr << denomPid[i] << " " << accDenomPid[i] << " (beta)" << endl;
#endif
				}

			}
	}
}

template <class BaseModel,typename WeightType>
void ModelSpecifics<BaseModel,WeightType>::doSortPid(bool useCrossValidation) {
/* For Cox model:
 *
 * We currently assume that hZ[k] are sorted in decreasing order by k.
 *
 */

//	cerr << "Copying Y" << endl;
//	// Copy y; only necessary if non-unique values in oY
//	nY.reserve(oY.size());
//	std::copy(oY.begin(),oY.end(),back_inserter(nY));
//	hY = const_cast<real*>(nY.data());

//	cerr << "Sorting PIDs" << endl;
//
//	std::vector<int> inverse_ranks;
//	inverse_ranks.reserve(K);
//	for (int i = 0; i < K; ++i) {
//		inverse_ranks.push_back(i);
//	}
//
//	std::sort(inverse_ranks.begin(), inverse_ranks.end(),
//			CompareSurvivalTuples<WeightType>(useCrossValidation, hKWeight, oZ));
//
//	nPid.resize(K, 0);
//	for (int i = 0; i < K; ++i) {
//		nPid[inverse_ranks[i]] = i;
//	}
//	hPid = const_cast<int*>(nPid.data());

//	for (int i = 0; i < K; ++i) {
//		cerr << oZ[inverse_ranks[i]] << endl;
//	}
//
//	cerr << endl;
//
//	for (int i = 0; i < K; ++i) {
//		cerr << oZ[i] << "\t" << hPid[i] << endl;
//	}
//
//	cerr << endl;
//
//	for (int i = 0; i < K; ++i) {
//		cerr << i << " -> " << hPid[i] << endl;
//	}
//
}

} // namespace

#endif /* MODELSPECIFICS_HPP_ */
