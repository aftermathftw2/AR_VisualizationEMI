using DSPLib;
using FFTWSharp;
using System;
using System.Collections;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Numerics;
using System.Runtime.InteropServices;
using UnityEngine;

public class PostProcessing
{
    private static double _samplingFreq = 1000000;

    public Tuple<List<int>, List<double>> correction;
    public Tuple<List<int>, List<double>> standard;
    BandwidthsData bandwidths;

    public PostProcessing()
    {
        bandwidths = new BandwidthsData(_samplingFreq);
        correction = ReadFile(Path.Combine(Directory.GetCurrentDirectory(), "Correction557.csv"));
        standard = ReadFile(Path.Combine(Directory.GetCurrentDirectory(), "StandardNRE01LimitLandApp.csv"));
    }

    public Tuple<double, double> CalculateFFT(double[] dataInput)
    {
        //should always get 200,000 datapoints
        //do all three bands
        //band 1 only does 1 window, 2 and 3 do multiple with overlap

        //initialize values to band 1
        int currentBandAlreadyDoneIndex = 0;
        double[] band = bandwidths.band1.band;
        double Bif = bandwidths.band1.bif;

        Tuple<double[], int> window = bandwidths.band1.window;
        int zeros = (int)Math.Pow(2, 18) - window.Item1.Length;


        //TODO find better solution than this list of value, frequency. 
        //tuples are read only... custom class?
        List<double> minOfAllWindows = new List<double>() { 500.0, 0.0 };

        for (int bandIndex = 0; bandIndex < 3; bandIndex++)
        {
            //get information of current band
            band = bandwidths.GetBand(bandIndex).band;
            Bif = bandwidths.GetBand(bandIndex).bif;
            window = bandwidths.GetBand(bandIndex).window;
            zeros = bandwidths.GetBand(bandIndex).zeros;


            //if the window is the same size as dataInput it will do only one loop
            //if it is smaller it will loop until there is not enough data for a full window size
            for (int stfftIndex = 0; stfftIndex < dataInput.Length; stfftIndex += window.Item2)
            {
                double[] correctedResult = new double[band.Length];
                double[] currentWindowDataInput = new double[window.Item1.Length];

                if (dataInput.Length - stfftIndex >= window.Item1.Length)
                {
                    //copy current window with overlap
                    Array.Copy(dataInput, stfftIndex, currentWindowDataInput, 0, window.Item1.Length);
                    //calculate the FFT of a single window and map results to closest band frequencies
                    double[] mappedMagLog = SingleWindowCalculationNewLibrary(currentWindowDataInput, window.Item1, zeros, band);
                    //correct values

                    double minValue = 10000;
                    int indexOf = 0;

                    for (int i = 0; i < mappedMagLog.Length; i++)
                    {
                        correctedResult[i] = standard.Item2[i + currentBandAlreadyDoneIndex] - (mappedMagLog[i] + 120 + correction.Item2[i + currentBandAlreadyDoneIndex]);
                        //correctedResult[i] = (mappedMagLog[i] + 120 + correction.Item2[i + currentBandAlreadyDoneIndex]);

                        //moved the min calculation over here to reduce for loop amounts
                        if (correctedResult[i] < minValue)
                        {
                            minValue = correctedResult[i];
                            indexOf = i;
                        }

                    }
                    //corrected result has the same frequency steps as the current band of the calculation
                    //find max of that window, check if it is higher than previous
                    //save max and its frequency
                    //save that curve too? since we wipe it next loop?


                    //should search for the minimum. since the least distance to the limitline is the highest value
                    //the actual value of the measurementPoints is the limitline at that frequency - corrected value, which is negative

                    //moved min calculation into already existing for loop above
                    //Tuple<double, int> trialMinValueAndIndexOf = FindMin(correctedResult);

                    Tuple<double, int> trialMinValueAndIndexOf = Tuple.Create(minValue, indexOf);

                    if (trialMinValueAndIndexOf.Item1 < minOfAllWindows[0])
                    {
                        minOfAllWindows[0] = trialMinValueAndIndexOf.Item1;
                        minOfAllWindows[1] = trialMinValueAndIndexOf.Item2;
                    }
                }
                else
                {
                    //not enough data for a full window. stfft is done
                }
            }
            currentBandAlreadyDoneIndex += band.Length;
        }

        //now we have all data after correction we can find the max value
        //The max of all values of each band each window
        return Tuple.Create(standard.Item2[(int)minOfAllWindows[1] + currentBandAlreadyDoneIndex - band.Length] - minOfAllWindows[0], band[(int)minOfAllWindows[1]]);
    }

    private static double[] SingleWindowCalculation(double[] dataInput, double[] window, int zeros, double[] band)
    {
        //calculate the FFT
        Complex[] cpxResult = SegmentFFT(dataInput, window, zeros);
        //double[] magLog = ToMagnitudeDBV(cpxResult); //moved this into the for loop of mapmaglog

        //calculate the FrequencySpan
        double[] fSpan = FrequencySpan(_samplingFreq, cpxResult.Length);

        //Map frequencySpan to frequency Band, save the index
        //new array, mappedMagLog, with magLog at those indices.    
        //double[] mappedMagLog = MapMagLog(magLog, band, fSpan); //made alternative method which includes TomagnitudeDBV
        double[] mappedMagLog = MapMagLog(cpxResult, band, fSpan);
        //now we have an array with the same frequency steps as our correction file

        return mappedMagLog;
    }

    private static double[] SingleWindowCalculationNewLibrary(double[] dataInput, double[] window, int zeros, double[] band)
    {
        double[] cpxResult = FftW(Multiply(dataInput, window), true); //execute fft
        double[] mag = MultiplyHalfAndLog(Magnitude(cpxResult), Math.Sqrt(2) / dataInput.Length, dataInput.Length / 2); //to magnitude, scale, half spectrum and log
        double[] fSpan = FrequencySpan(_samplingFreq, dataInput.Length / 2);
        double[] mappedMagLog = MapMagLog(mag, band, fSpan);

        return mappedMagLog;
    }

    private static Tuple<double, int> FindMax(double[] input)
    {
        double maxValue = -10000;
        int indexOf = 0;

        for (int i = 0; i < input.Length; i++)
        {
            if (input[i] > maxValue)
            {
                maxValue = input[i];
                indexOf = i;
            }
        }
        return Tuple.Create(maxValue, indexOf);
    }

    private static Tuple<double, int> FindMin(double[] input)
    {
        double minValue = 10000;
        int indexOf = 0;

        for (int i = 0; i < input.Length; i++)
        {
            if (input[i] < minValue)
            {
                minValue = input[i];
                indexOf = i;
            }
        }
        return Tuple.Create(minValue, indexOf);
    }

    /// <summary>
    /// Maps frequency span fSpan to the frequencies in bandwidth.
    /// </summary>
    /// <param name="magLog"></param>
    /// <param name="bandwidth"></param>
    /// <param name="fSpan"></param>
    /// <returns>magLog at these frequencies in the bandwidth</returns>
    private static double[] MapMagLog(double[] magLog, double[] bandwidth, double[] fSpan)
    {
        double[] mappedMagLog = new double[bandwidth.Length];
        for (int i = 0; i < bandwidth.Length; i++)
        {
            double lowestDiff = 10000;
            int indexOfLowestDiff = 0;

            for (int ii = 0; ii < fSpan.Length; ii++)
            {
                if (Math.Abs(bandwidth[i] - fSpan[ii]) < lowestDiff)
                {
                    lowestDiff = bandwidth[i] - fSpan[ii];
                    indexOfLowestDiff = ii;
                }
                if (fSpan[ii] > bandwidth[bandwidth.Length - 1])
                {
                    //To limit the for loop to only frequencies within the bands (under 100k)
                    //if the span goes above 100,000 it will stop the loop
                    ii = fSpan.Length - 1;
                }
            }
            mappedMagLog[i] = magLog[indexOfLowestDiff];
        }
        return mappedMagLog;
    }
    /// <summary>
    /// Maps frequency span fSpan to the frequencies in bandwidth. Input cpxResult to limit for loops
    /// </summary>
    /// <param name="cpxResult"></param>
    /// <param name="bandwidth"></param>
    /// <param name="fSpan"></param>
    /// <returns>magLog at these frequencies in the bandwidth</returns>
    private static double[] MapMagLog(Complex[] cpxResult, double[] bandwidth, double[] fSpan)
    {
        double[] mappedMagLog = new double[bandwidth.Length];
        for (int i = 0; i < bandwidth.Length; i++)
        {
            double lowestDiff = 10000;
            int indexOfLowestDiff = 0;

            for (int ii = 0; ii < fSpan.Length; ii++)
            {
                if (Math.Abs(bandwidth[i] - fSpan[ii]) < lowestDiff)
                {
                    lowestDiff = bandwidth[i] - fSpan[ii];
                    indexOfLowestDiff = ii;
                }
                if (fSpan[ii] > bandwidth[bandwidth.Length - 1])
                {
                    //To limit the for loop to only frequencies within the bands (under 100k)
                    //if the span goes above 100,000 it will stop the loop
                    ii = fSpan.Length - 1;
                }
            }
            mappedMagLog[i] = ToMagnitudeDBV(cpxResult[indexOfLowestDiff]);
        }
        return mappedMagLog;
    }
    /// <summary>
    /// Calculate the FFT of a single Segment. segmentData and window should be the same size
    /// </summary>
    private static Complex[] SegmentFFT(double[] segmentData, double[] window, int zeros)
    {
        System.Diagnostics.Stopwatch stpwatch = new System.Diagnostics.Stopwatch();
        stpwatch.Start();
        // Instantiate & initialize the fft class
        DSPLib.FFT fft = new DSPLib.FFT();
        fft.Initialize((uint)window.Length, (uint)zeros);

        // Perform a FFT
        Complex[] cpxResult = fft.Execute(Multiply(segmentData, window));
        stpwatch.Stop();
        Debug.Log(String.Format("Actual fft time: {0}", (float) (stpwatch.ElapsedMilliseconds)));
        return cpxResult;
    }

    static double[] Magnitude(double[] x)
    {
        int n = x.Length / 2;
        var y = new double[n];
        for (int i = 0; i < n; i++)
        {
            y[i] = Math.Sqrt(x[2 * i] * x[2 * i] + x[2 * i + 1] * x[2 * i + 1]);
        }
        return y;
    }

    /// <summary>
    /// Return the Frequency Array for the currently defined FFT.
    /// Takes into account the total number of points and zero padding points that were defined.
    /// </summary>
    /// <param name="samplingFrequencyHz"></param>
    /// <param name="mLengthHalf"></param>
    /// <returns></returns>
    public static double[] FrequencySpanNew(double samplingFrequencyHz, UInt32 mLengthHalf)
    {
        UInt32 points = (UInt32)mLengthHalf;
        double[] result = new double[points];
        double stopValue = samplingFrequencyHz / 2.0;
        double increment = stopValue / ((double)points - 1.0);

        for (Int32 i = 0; i < points; i++)
            result[i] += increment * i;

        return result;
    }

    /// <summary>
    /// Computes the fast Fourier transform of a 1-D array of real or complex numbers.
    /// </summary>
    /// <param name="data">Input data.</param>
    /// <param name="real">Real or complex input flag.</param>
    /// <returns>Returns the FFT.</returns>
    private static double[] FftW(double[] data, bool real)
    {
        // If the input is real, make it complex
        if (real)
        {
            data = ToComplex(data);
        }
        // Get the length of the array
        int n = data.Length;
        /* Allocate an unmanaged memory block for the input and output data.
         * (The input and output are of the same length in this case, so we can use just one memory block.) */
        IntPtr ptr = fftw.malloc(n * sizeof(double));
        // Pass the managed input data to the unmanaged memory block
        Marshal.Copy(data, 0, ptr, n);
        // Plan the FFT and execute it (n/2 because complex numbers are stored as pairs of doubles)
        IntPtr plan = fftw.dft_1d(n / 2, ptr, ptr, fftw_direction.Forward, fftw_flags.Estimate);
        fftw.execute(plan);
        // Create an array to store the output values
        var fft = new double[n];
        // Pass the unmanaged output data to the managed array
        Marshal.Copy(ptr, fft, 0, n);
        // Do some cleaning
        fftw.destroy_plan(plan);
        fftw.free(ptr);
        fftw.cleanup();
        // Return the FFT output
        return fft;
    }
    /// <summary>
    /// Interlaces an array with zeros to match the FFTW convention of representing complex numbers.
    /// </summary>
    /// <param name="real">An array of real numbers.</param>
    /// <returns>Returns an array of complex numbers.</returns>
    private static double[] ToComplex(double[] real)
    {
        int n = real.Length;
        var comp = new double[n * 2];
        for (int i = 0; i < n; i++)
            comp[2 * i] = real[i];
        return comp;
    }

    /* OBSOLETE METHOD
    private static double[] convertToPeak(Complex[] cpxResult)
    {
        //cpxResult = cpxResult.Select(c => c / cpxResult.Length).ToArray(); // matlab does this, result is worse here,dsplib already scales          (These 3 lines are from a matlab script I was provided. I do not understand them.....)
        //double[] mag = ToMagnitude(cpxResult);   //take absolute value
        //mag = mag.Select(d => 2.0 * d * d).ToArray();                   //all values *2 because one sided spectrum                  (These 3 lines are from a matlab script I was provided. I do not understand them.....)
        //mag[0] = mag[0] / 2.0;                                      //dc term back to original                                  (These 3 lines are from a matlab script I was provided. I do not understand them.....)
        //double[] magLog = ToMagnitudeDBV(mag);

        //changed everything to 1 for loop.

        double[] magLog = ToMagnitudeDBV(cpxResult);
        return magLog;
    }
    */

    /// <summary>
    /// Convert Complex DFT/FFT Result to: Magnitude Vrms.
    /// </summary>
    /// Taken from DSPlib!!!
    /// <param name="rawFFT"></param>
    /// <returns>double[] Magnitude Format (Vrms)</returns>
    public static double[] ToMagnitude(Complex[] rawFFT)
    {
        UInt32 np = (UInt32)rawFFT.Length;
        double[] mag = new double[np];
        for (UInt32 i = 0; i < np; i++)
        {
            mag[i] = rawFFT[i].Magnitude;
        }

        return mag;
    }

    /// <summary>
    /// Convert Magnitude FT Result to: Magnitude dBVolts
    /// </summary>
    /// TAKEN FROM DSPlib!!!
    /// <param name="magnitude"></param>
    /// <returns>double[] array</returns>
    public static double[] ToMagnitudeDBV(double[] magnitude)
    {
        UInt32 np = (UInt32)magnitude.Length;
        double[] magDBV = new double[np];
        for (UInt32 i = 0; i < np; i++)
        {
            double magVal = magnitude[i];
            if (magVal <= 0.0)
                magVal = double.Epsilon; //to account for 0 which becomes infinity after log

            magDBV[i] = 20 * System.Math.Log10(magVal) - 59; //59 difference between perfect and PicoScope. Constant difference....
        }

        return magDBV;
    }

    /// <summary>
    /// Convert Complex DFT/FFT Result to: Log Magnitude dBV
    /// </summary>
    /// <param name="rawFFT"> Complex[] input array"></param>
    /// <returns>double[] Magnitude Format (dBV)</returns>
    /// TAKEN FROM DSPLIB
    public static double[] ToMagnitudeDBV(Complex[] rawFFT)
    {
        UInt32 np = (UInt32)rawFFT.Length;
        double[] mag = new double[np];
        for (UInt32 i = 0; i < np; i++)
        {
            double magVal = rawFFT[i].Magnitude;

            if (magVal <= 0.0)
                magVal = double.Epsilon;

            mag[i] = 20 * System.Math.Log10(magVal);
        }

        return mag;
    }

    public static double ToMagnitudeDBV(Complex cpxSingle)
    {
        double magVal = cpxSingle.Magnitude;
        if (magVal <= 0.0)
        {
            magVal = double.Epsilon;
        }
        return (20 * System.Math.Log10(magVal));
    }


    /// <summary>
    /// result[] = a[] * b[]
    /// </summary>
    public static Double[] Multiply(Double[] a, Double[] b)
    {
        Debug.Assert(a.Length == b.Length, "Length of arrays a[] and b[] must match.");

        double[] result = new double[a.Length];
        for (UInt32 i = 0; i < a.Length; i++)
            result[i] = a[i] * b[i];

        return result;

    }

    /// <summary>
    /// result[] = a[] * b
    /// </summary>
    public static Double[] MultiplyHalfAndLog(Double[] a, Double b, int endIndex)
    {
        double[] result = new double[endIndex];
        for (UInt32 i = 0; i < endIndex; i++)
        {
            result[i] = a[i] * b;
            if (result[i] <= 0.0)
            {
                result[i] = double.Epsilon;
            }
            result[i] = 20 * Math.Log10(result[i]) - 56.4; //offset...
        }
        

        return result;
    }

    /// <summary>
    /// Return the Frequency Array for a certain sampling frequency and number of points.
    /// </summary>
    /// <param name="samplingFrequencyHz"></param>
    /// <param name="nOfPoints"></param>
    /// TAKEN FROM DSPLib!!!
    /// <returns></returns>
    public static double[] FrequencySpan(double samplingFrequencyHz, int nOfPoints)
    {
        UInt32 points = (UInt32)(nOfPoints);
        double[] result = new double[points];
        double stopValue = samplingFrequencyHz / 2.0;
        double increment = stopValue / ((double)points - 1.0);

        for (Int32 i = 0; i < points; i++)
            result[i] += increment * i;

        return result;
    }


    public static Tuple<double[], int> GenerateNormalizedWindow(double[] bandwidth, double sampling_freq, double bif)
    {
        double df = bandwidth[1] - bandwidth[0]; //frequency step size
        int N = (int)(sampling_freq / df - 1); //length of segment
        double[] n = GetSteppedSequence(-(double)N / 2.0, (double)N / 2.0, N).ToArray(); //array from -N/2 to N/2 with stepsize 1
        double Ts = (double)1.0 / sampling_freq; //period

        double[] window = new double[n.Length];
        for (int i = 0; i < n.Length; i++)
        {
            window[i] = (1.0 / ((N * Ts) * N)) * Math.Exp(-0.5 * Math.Pow(Math.PI / Math.Sqrt(2 * Math.Log(2)) * bif * Ts * n[i], 2));
        }
        double coherentGain = (1.0 / N) * window.Sum(); //calculate coherentGain to normalize window
        double[] normWindow = window.Select(d => d / coherentGain).ToArray(); //Normalize window

        //calculation of dt...??
        double Foverlap = 1 - (df / bif); //percentage overlap
        int Noverlap = (int)Math.Ceiling(Foverlap * N); //number of overlapped samples

        return Tuple.Create(normWindow, Noverlap);
    }

    public static IEnumerable<double> GetSteppedSequence(double from, double to, int numberOfSteps)
    {
        if (numberOfSteps < 1)
        {
            throw new ArgumentOutOfRangeException(nameof(numberOfSteps), "Number of steps must be greater than zero.");
        }

        var stepSize = (to - from) / numberOfSteps;
        return Enumerable.Range(0, numberOfSteps + 1).Select(stepIndex => from + stepIndex * stepSize);
    }

    public static Tuple<List<int>, List<double>> ReadFile(string filePath)
    {
        List<int> item1 = new List<int>();
        List<double> item2 = new List<double>();

        using (var reader = new StreamReader(filePath))
        {
            while (!reader.EndOfStream)
            {
                var splitLine = reader.ReadLine().Split(';');

                int currentItem1;
                double currentItem2;

                if (Int32.TryParse(splitLine[0], out currentItem1))
                {
                    item1.Add(currentItem1);
                }
                else
                {
                    Console.WriteLine("not an integer in csv file");
                    //TODO: add custom exception
                    item1.Add(-100000000);
                }

                if (Double.TryParse(splitLine[1], System.Globalization.NumberStyles.Any, System.Globalization.CultureInfo.InvariantCulture, out currentItem2))
                {
                    item2.Add(currentItem2);
                }
                else
                {
                    Console.WriteLine("not a double in csv file: {0}", splitLine[1]);
                    //TODO: add custom exception
                    item2.Add(-10000000);
                }
            }
        }
        return Tuple.Create(item1, item2);
    }
}

public class BandwidthsData
{
    public Band band1;
    public Band band2;
    public Band band3;
    public double[] totalband;

    public BandwidthsData(double samplingFreq)
    {
        band1 = new Band(30, 1000, 5, 10, samplingFreq);
        band2 = new Band(1000, 10000, 50, 100, samplingFreq);
        band3 = new Band(10000, 100000, 500, 1000, samplingFreq);

        totalband = new double[band1.band.Length + band2.band.Length + band3.band.Length];
        band1.band.CopyTo(totalband, 0);
        band2.band.CopyTo(totalband, band1.band.Length);
        band3.band.CopyTo(totalband, band2.band.Length);
    }

    public Band GetBand(int index)
    {
        Band outBand = band1;
        if (index == 0)
        {
            outBand = band1;
        }
        else if (index == 1)
        {
            outBand = band2;
        }
        else if (index == 2)
        {
            outBand = band3;
        }
        else
        {
            Console.WriteLine("Uncompatible band index, defaulting to band1. index 0 = band1, 1 = band2, 2 = band 3");
        }

        return outBand;
    }
}

public class Band
{
    public double[] band;
    public int bif;
    public int df;
    public Tuple<double[], int> window;
    public int zeros;

    public Band(int from, int to, int df, int Bif, double samplingFreq)
    {
        this.bif = Bif;
        this.df = df;
        band = PostProcessing.GetSteppedSequence(from, to, (to - from) / df).ToArray();
        window = PostProcessing.GenerateNormalizedWindow(band, samplingFreq, bif);

        //could make function to detect which 2 power is closest to the window length.
        //Since we know the bands we want to analyse i decided to just use an ugly else if.
        if (window.Item1.Length == 200000)
        {
            //window is 200,000 long. nearest 2 power is 2^18
            zeros = (int)Math.Pow(2, 18) - window.Item1.Length;
        }
        else if (window.Item1.Length == 20000)
        {
            //window is 20,000 long. nearest 2 power is 2^15
            zeros = (int)Math.Pow(2, 15) - window.Item1.Length;
        }
        else if (window.Item1.Length == 2000)
        {
            //window is 2,000 long. nearest 2 power is 2^11
            zeros = (int)Math.Pow(2, 11) - window.Item1.Length;
        }
    }
}
