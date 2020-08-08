#pragma once
#include <vector>

class SignalProcessing
{
public:

#pragma region Reordering helpers

	std::vector<double> ReverseVector( std::vector<double> input )
	{
		std::vector<double> output( input.size() );

		for (int i = 0; i < input.size(); ++i)
		{
			int revIndex = input.size() - i - 1;
			output[i] = input[revIndex];
		}

		return output;
	}

	std::vector<double> OffsetByHalf( std::vector<double> vec, int resolution )
	{
		std::vector<double> ret( resolution );
		int counter = resolution / 2;
		for (int i = 0; i < resolution; i++)
		{
			if (counter == resolution)
				counter = 0;

			ret[i] = vec[counter];

			++counter;
		}
		return ret;
	}

	std::vector<double> RollBackwards( std::vector<double> vec, int resolution )
	{
		std::vector<double> ret( resolution );
		double last = vec[0];
		for (int i = 0; i < resolution - 1; i++)
		{
			ret[i] = vec[i + 1];
		}
		ret[resolution - 1] = last;
		return ret;
	}

	std::vector<double> RollForwards( std::vector<double> vec, int resolution )
	{
		std::vector<double> ret( resolution );
		double last = vec[resolution - 1];
		for (int i = 0; i < resolution - 1; i++)
		{
			ret[i + 1] = vec[i];
		}
		ret[0] = last;
		return ret;
	}

#pragma endregion

#pragma region Signal generator

	double WaveFunction( double x )
	{
		double w0 = 10.0 * 3.14159265359;
		return 5.0 * sin( w0 * x ) + 3.0 * cos( 4.0 * w0 * x );
	}

	std::vector<double> WaveFunctionGenerator( double start, double end, int resolution )
	{
		double B = 2.0 * 3.14159265359 / ( end - start );

		std::vector<double> ret( resolution );

		int counter = 0;
		double x = start;

		while (counter < resolution)
		{
			ret[counter] = WaveFunction( x );

			x += ( end - start ) / resolution;

			++counter;
		}
		return ret;
	}

	std::vector<double> SquareWaveGeneratorSimple( double start, double end, double max, double min, double shift, int size )
	{
		double B = 2.0 * 3.14159265359 / ( end - start );

		std::vector<double> ret( size );

		double x = start;

		for (int i = 0; i < size; ++i)
		{
			double sinValue = sin( x + 2.0f * shift );
			if (sinValue >= 0)
				ret[i] = max;
			else
				ret[i] = min;

			x += ( end - start ) / size;
		}
		return ret;
	}

#pragma endregion

#pragma region Signal processing

	double Correlation( std::vector<double> f, std::vector<double> g, int resolution )
	{
		double signal_1_energy = 0.0;
		double signal_2_energy = 0.0;
		double scaling_factor = 0.0;

		double correlationReturn = 0.0;

		//scaling values
		for (int i = 0; i < resolution; ++i)
		{
			signal_1_energy += f[i] * f[i];
			signal_2_energy += g[i] * g[i];
		}

		scaling_factor = sqrt( signal_1_energy * signal_2_energy );

		for (int k = 0; k < resolution; ++k)
		{
			correlationReturn += f[k] * g[k] / scaling_factor;
		}
		return correlationReturn;
	}

	std::vector<double> CrossCorrelation( std::vector<double> f, std::vector<double> g, int resolution )
	{
		double signal_1_energy = 0.0;
		double signal_2_energy = 0.0;
		double scaling_factor = 0.0;

		std::vector<double> correlationReturn( resolution );

		std::vector<double> signal_2_lag_compensated = g;

		//scaling values
		for (int i = 0; i < resolution; ++i)
		{
			signal_1_energy += f[i] * f[i]; // energy signal 1
			signal_2_energy += g[i] * g[i]; // energy signal 2
		}
		scaling_factor = sqrt( signal_1_energy * signal_2_energy ); // calculate total scaling factor
		for (int j = 0; j < resolution; ++j)
		{
			double correlation = 0.0;
			for (int k = 0; k < resolution; ++k)
			{
				correlation += f[k] * signal_2_lag_compensated[k] / scaling_factor;
			}
			correlationReturn[j] = correlation;
			signal_2_lag_compensated = RollForwards( signal_2_lag_compensated, resolution );
		}
		return correlationReturn;
	}

	std::vector<double> FastFourierTransform( std::vector<double> fi, int resolution )
	{
		std::vector<double> f = ReverseVector( fi );
		std::vector<double> FFTReturn( resolution );

		double pi = 3.14159265359;

		for (int m = 0; m < resolution; ++m)
		{
			double FFT = 0.0;
			double FFT_img = 0.0;
			double FFT_real = 0.0;
			for (int k = 0; k < resolution; ++k)
			{
				FFT_real += f[k] * cos( 2.0 * pi * (double)m * (double)k / (double)resolution );
				FFT_img += f[k] * sin( -2.0 * pi * (double)m * (double)k / (double)resolution );
			}
			FFT = sqrt( ( FFT_real * FFT_real ) + ( FFT_img * FFT_img ) );
			FFTReturn[m] = FFT;
		}

		return FFTReturn;
	}

	std::vector<double> InverseFastFourierTransform( std::vector<double> f, int resolution )
	{
		std::vector<double> FFTReturn( resolution );

		double pi = 3.14159265359;

		for (int m = 0; m < resolution; ++m)
		{
			double FFT = 0.0;
			double FFT_img = 0.0;
			double FFT_real = 0.0;
			for (int k = 0; k < resolution; ++k)
			{
				FFT_real += f[k] * cos( 2.0 * pi * (double)m * (double)k / (double)resolution );
				FFT_img += f[k] * sin( 2.0 * pi * (double)m * (double)k / (double)resolution );
			}
			FFT = FFT_real - FFT_img;
			FFTReturn[m] = FFT / (double)resolution;
		}

		return FFTReturn;
	}

	std::vector<double> ConvolutionByFFT( std::vector<double> f, std::vector<double> g, int resolution )
	{
		std::vector<double> temp( resolution );
		std::vector<double> ret( resolution );
		auto f_fft = FastFourierTransform( f, resolution );
		auto g_fft = FastFourierTransform( g, resolution );

		for (int i = 0; i < resolution; ++i)
		{
			temp[i] = f_fft[i] * g_fft[i];
		}

		double pi = 3.14159265359;

		for (int m = 0; m < resolution; ++m)
		{
			double FFT = 0.0;
			double FFT_img = 0.0;
			double FFT_real = 0.0;
			for (int k = 0; k < resolution; ++k)
			{
				FFT_real += temp[k] * std::cos( 2.0 * pi * (double)m * (double)k / (double)resolution );
				FFT_img += temp[k] * std::sin( -2.0 * pi * (double)m * (double)k / (double)resolution );
			}
			FFT = sqrt( ( FFT_real * FFT_real ) + ( FFT_img * FFT_img ) );
			ret[m] = abs( FFT ) / (double)resolution;
		}

		return ret;
	}

#pragma endregion

};