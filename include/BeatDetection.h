#ifndef BEAT_DETECTOR_H
#define BEAT_DETECTOR_H

#include <vector>
#include <math.h>

namespace beat_detection {

namespace BeatTypes {
	enum BeatType {
		LOW = 0,
		MID,
		HIGH,
		NUM_TYPES
	};
}

enum {
	DEFAULT_SPECTRUM_SIZE  = 1024,
	DEFAULT_BAND_SIZE      = 64,
	DEFAULT_HISTORY_SIZE   = 40,
	DEFAULT_DECIBEL_CUTOFF = 125,

	DEFAULT_LOW_CUTOFF = 4,
	DEFAULT_MID_CUTOFF = 16,
	DEFAULT_HIGH_CUTOFF = 32,

	DEFAULT_LOW_THRESHOLD = 150,
	DEFAULT_MID_THRESHOLD = 130,
	DEFAULT_HIGH_THRESHOLD = 80,
};

template<class T>
class SampleQueue
{
public:
	SampleQueue(size_t maxSize)
		: mSamples(maxSize, (T)0), mTotalSamples(0), mAverage((T)0), mTotal((T)0), mFull(false) { }
	SampleQueue(const SampleQueue& other)
		: mSamples(other.mSamples), mTotalSamples(other.mTotalSamples), mAverage(other.mAverage), mTotal(other.mTotal), mFull(other.mFull) { }

	void addSample(T sample)
	{
		size_t index = mTotalSamples++ % mSamples.size();
		if (!mFull && mTotalSamples >= mSamples.size())
			mFull = true;
		mTotal         += sample - mSamples[index];
		mSamples[index] = sample;
		mAverage        = mTotal / (T)mSamples.size();
		mVariance       = (T)0;
		for (size_t i = 0; i < getSampleCount(); ++i)
		{
			T delta = mSamples[i] - mAverage;
			//threshold of 145 works best with variance, and 157 with standard deviation
			//xmVariance += delta * delta;
			mVariance += std::abs(delta);			//standard deviation
		}
		mVariance /= getSampleCount();
	}

	inline T getAverage() const { return mAverage; }
	inline T getVariance() const { return mVariance; }
	inline size_t getSampleCount() const { return mFull ? mSamples.size() : mTotalSamples; }

protected:
	std::vector<T> mSamples;
	size_t         mTotalSamples;
	T              mTotal;
	T              mAverage;
	T              mVariance;
	bool           mFull;
};

template<class T>
struct BeatDetectionData
{
	typedef SampleQueue<T> Samples;

	BeatDetectionData(size_t sampleCountIn, size_t bandCountIn, size_t historyCountIn, T decibelCutoffIn = (T)125) 
		: sampleCount(sampleCountIn), bandCount(bandCountIn), historyCount(historyCountIn), sampleStorage(sampleCount, (T)0), bandStorage(bandCountIn, (T)0), dbStorage(sampleCount/2, (T)0), decibelCutoff(decibelCutoffIn) {
		
		for (size_t i = 0; i < bandCountIn; ++i) {
			history.push_back( Samples(historyCountIn) );
			}

		for (size_t i = 0; i < BeatTypes::NUM_TYPES; ++i) {
			energy[i]  = (T)0;
			counter[i] = 0;
		}

		cutoff[BeatTypes::LOW]  = (T)DEFAULT_LOW_CUTOFF;
		cutoff[BeatTypes::MID]  = (T)DEFAULT_MID_CUTOFF;
		cutoff[BeatTypes::HIGH] = (T)DEFAULT_HIGH_CUTOFF;

		threshold[BeatTypes::LOW]  = (T)DEFAULT_LOW_THRESHOLD;
		threshold[BeatTypes::MID]  = (T)DEFAULT_MID_THRESHOLD;
		threshold[BeatTypes::HIGH] = (T)DEFAULT_HIGH_THRESHOLD;
	}

	size_t         sampleCount;
	size_t         bandCount;
	size_t         historyCount;
	std::vector< Samples > history;
	std::vector<T> sampleStorage;
	std::vector<T> bandStorage;
	std::vector<T> dbStorage;

	T       energy[BeatTypes::NUM_TYPES];
	T       decibelCutoff;

	size_t  counter[BeatTypes::NUM_TYPES];
	size_t  cutoff[BeatTypes::NUM_TYPES];
	size_t  threshold[BeatTypes::NUM_TYPES];
};

template<class T>
class FFTTransform {
public:
	virtual void operator()(const T* in, T* out, size_t length) = 0;
};

template<class T>
class BeatCallback {
public:
	virtual void operator()(BeatTypes::BeatType beatType, T beatEnergy) = 0;
};

template<class T>
class SimpleBeatCallback : public BeatCallback<T> {
public:
	SimpleBeatCallback() {
		for (size_t i = 0; i < BeatTypes::NUM_TYPES; ++i)
		{
			beat[i]   = false;
			energy[i] = (T)0;
		}
	}
	virtual void operator()(BeatTypes::BeatType beatType, T beatEnergy) {
		beat[beatType]   = true;
		energy[beatType] = beatEnergy;
	}
	bool beat[BeatTypes::NUM_TYPES];
	T    energy[BeatTypes::NUM_TYPES];
};

template<class T>
class BeatDetection {
public:
	BeatDetection(FFTTransform<T>& fft,
	              size_t spectrumSize = DEFAULT_SPECTRUM_SIZE,
	              size_t bandSize     = DEFAULT_BAND_SIZE,
	              size_t historySize  = DEFAULT_HISTORY_SIZE,
	              T     decibelCutoff = (T)DEFAULT_DECIBEL_CUTOFF)
		: mFFT(fft), mData(spectrumSize, bandSize, historySize, decibelCutoff)
	{

	}

	void process(T* samples, BeatCallback<T>& callback)
	{
		/*
		if (*std::max_element(samples, samples + mData.sampleCount) > 0.001f)
		{
			for (size_t i = 0; i < mData.sampleCount/2; ++i)
			{
				T db = 10.0f  * static_cast<T>(std::log10(samples[i])) * 2.0f;
				if (db < -mData.decibelCutoff)
					db = -mData.decibelCutoff;

				db /= -static_cast<T>(mData.decibelCutoff);
				mData.dbStorage[i] = (T)1 - db;
			}
		}
		*/

		for (size_t i = 0; i < BeatTypes::NUM_TYPES; ++i)
			if (mData.counter[i] > 0)
				--mData.counter[i];

		mFFT(samples, &mData.sampleStorage[0], mData.sampleStorage.size());

		//collect the energy
		size_t j = 0;
		size_t samplesPerBand = mData.sampleCount / mData.bandCount;
		for (size_t i = 0; i < mData.sampleCount; ++i)
		{
			j = i / (samplesPerBand);
			mData.bandStorage[j] += mData.sampleStorage[i] * 10.0f;//normalized
		}

		//add energy to history
		for (size_t i = 0; i < mData.bandCount; ++i)
		{
			mData.bandStorage[i] /= (samplesPerBand);
			mData.history[i].addSample(mData.bandStorage[i]);
		}

		// check for beats
		size_t counts[BeatTypes::NUM_TYPES];
		T      peaks[BeatTypes::NUM_TYPES];
		T      means[BeatTypes::NUM_TYPES];

		for (size_t i = 0; i < BeatTypes::NUM_TYPES; ++i)
		{
			counts[i] = 0;
			peaks[i] = (T)0;
			means[i] = (T)0;
		}

		for (size_t i = 0 ; i < mData.bandCount; i++)
		{
			for (size_t j = 0; j < BeatTypes::NUM_TYPES; ++j)
			{
				means[j] += mData.history[i].getAverage();
				peaks[j] += mData.bandStorage[i];
				if (mData.bandStorage[i] > 
					(mData.history[i].getVariance()/mData.history[i].getAverage())
					+ (mData.history[i].getAverage()*(mData.threshold[j]/(T)100)))
				{
					++counts[j];
				}
			}
		}

		for (size_t j = 0; j < BeatTypes::NUM_TYPES; ++j)
		{
			if (0 == mData.counter[j])
			{
				if (counts[j] > (mData.cutoff[j]/2))
				{
					mData.counter[j] = 1;
					mData.energy[j] = peaks[j];
					callback(static_cast<BeatTypes::BeatType>(j), mData.energy[j]);
				}
			}
		}
	}


protected:
	FFTTransform<T>&     mFFT;
	BeatDetectionData<T> mData;
};

}

#endif
