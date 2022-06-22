using System.Collections;
using System.Collections.Generic;
using System.Text;
using System;
using UnityEngine;

using PS5000AImports;
using PicoPinnedArray;
using PicoStatus;
using System.Threading;

public class PicoScopeHandler
{
    public readonly short _handle;
    int _channelCount;
    private ChannelSettings[] _channelSettings;

    short[][] appBuffers;
    short[][] buffers;

    ushort[] inputRanges = { 10, 20, 50, 100, 200, 500, 1000, 2000, 5000, 10000, 20000, 50000 };

    bool _autoStop;
    bool _powerSupplyConnected;
    bool _ready = false;

    short _maxValue;
    short _trig = 0;

    int _noEnabledChannels;
    int _sampleCount;

    uint _trigAt = 0;
    uint _startIndex;

    Imports.DeviceResolution _resolution = Imports.DeviceResolution.PS5000A_DR_8BIT;

    Imports.Range _firstRange = Imports.Range.Range_10mV;
    Imports.Range _lastRange = Imports.Range.Range_20V;

    private bool running;
    short[][] current200kPackage;
    public bool doneCopying;
    public PostProcessing post;

    public bool usePingPong;

    /****************************************************************************
    * Callback
    * Used by ps5000a data streaming collection calls, on receipt of data.
    * Used to set global flags etc checked by user routines
    ****************************************************************************/
    void StreamingCallback(short handle,
                            int noOfSamples,
                            uint startIndex,
                            short overflow,
                            uint triggerAt,
                            short triggered,
                            short autoStop,
                            IntPtr pVoid)
    {
        // used for streaming
        _sampleCount = noOfSamples;
        _startIndex = startIndex;
        _autoStop = autoStop != 0;

        _ready = true;

        // flags to show if & where a trigger has occurred
        _trig = triggered;
        _trigAt = triggerAt;

        if (_sampleCount != 0)
        {
            for (int ch = 0; ch < _channelCount * 2; ch += 2)
            {
                if (_channelSettings[(int)(Imports.Channel.ChannelA + (ch / 2))].enabled)
                {

                    Array.Copy(buffers[ch], _startIndex, appBuffers[ch], _startIndex, _sampleCount); //max
                    Array.Copy(buffers[ch + 1], _startIndex, appBuffers[ch + 1], _startIndex, _sampleCount);//min

                }
            }
        }
    }


    /****************************************************************************
     * adc_to_mv
     *
     * Convert an 16-bit ADC count into millivolts
     ****************************************************************************/
    double adc_to_mv(int raw, int ch)
    {
        return (double)((double)raw * (double)inputRanges[ch]) / (double) _maxValue;
    }

    /****************************************************************************
     * mv_to_adc
     *
     * Convert a millivolt value into a 16-bit ADC count
     *
     *  (useful for setting trigger thresholds)
     ****************************************************************************/
    short mv_to_adc(short mv, short ch)
    {
        return (short)((mv * _maxValue) / inputRanges[ch]);
    }

    /****************************************************************************
   * Stream Data Handler
   * - Used by the two stream data examples - untriggered and triggered
   * Inputs:
   * - unit - the unit to sample on
   * - preTrigger - the number of samples in the pre-trigger phase 
   *					(0 if no trigger has been set)
   ***************************************************************************/
    short[][] StreamDataHandler(uint preTrigger, uint rawDataSamples)
    {
        int sampleCount = 2048 * 100; /*  *100 is to make sure buffer large enough */

        appBuffers = new short[_channelCount * 2][];
        buffers = new short[_channelCount * 2][];
        current200kPackage = new short[_channelCount * 2][];

        uint postTriggerSamples = rawDataSamples; // Number of raw data samples
        int totalSamples = 0;
        uint triggeredAt = 0;
        uint sampleInterval = 1;
        uint status;

        for (int ch = 0; ch < _channelCount * 2; ch += 2) // create data buffers
        {
            buffers[ch] = new short[sampleCount];
            buffers[ch + 1] = new short[sampleCount];

            appBuffers[ch] = new short[sampleCount];
            appBuffers[ch + 1] = new short[sampleCount];

            current200kPackage[ch] = new short[200000];
            current200kPackage[ch + 1] = new short[200000];

            status = Imports.SetDataBuffers(_handle, (Imports.Channel)(ch / 2), buffers[ch], buffers[ch + 1], sampleCount, 0, Imports.RatioMode.None);
        }

        _autoStop = false;
        status = Imports.RunStreaming(_handle, ref sampleInterval, Imports.ReportedTimeUnits.MicroSeconds, preTrigger, postTriggerSamples, 1, 1, Imports.RatioMode.None, (uint)sampleCount);
        string tHing1 = buffers[0][0].ToString();

        while (!_autoStop)
        {
            /* Poll until data is received. Until then, GetStreamingLatestValues wont call the callback */
            //wait until new data is in buffers. then copy it to appBuffers
            Thread.Sleep(0);
            _ready = false;
            status = Imports.GetStreamingLatestValues(_handle, StreamingCallback, IntPtr.Zero);


            Console.Write((status > StatusCodes.PICO_OK && status != StatusCodes.PICO_BUSY /*PICO_BUSY*/) ? "Status =  {0}\n" : "", status);

            if (_ready && _sampleCount > 0) /* can be ready and have no data, if autoStop has fired */
            {
                if (_trig > 0)
                {
                    triggeredAt = (uint)totalSamples + _trigAt;
                }

                totalSamples += _sampleCount;
                Console.Write("\nCollected {0,4} samples, index = {1,5}, Total = {2,5}", _sampleCount, _startIndex, totalSamples);

                if (totalSamples > 199999)
                {
                    //latest package is now in app buffers and ready to be used
                    doneCopying = false;
                    Array.Copy(appBuffers[0], 0, current200kPackage[0], 0, 200000);
                    // top line appBuffers[0] is max adc, is min important?:
                    //Array.Copy(appBuffers[1], _startIndex, currentPackage[1], _startIndex, _sampleCount);
                    doneCopying = true;

                    //can be made to include multiple channels, see copy statement in callback function. for loop over channels

                    totalSamples = 0;
                }



            }
        }

        Imports.Stop(_handle);

        if (!_autoStop)
        {
            Debug.Log(String.Format("\nData collection aborted"));
        }

        return appBuffers;
    }


    /****************************************************************************
    * Initialise unit' structure with Variant specific defaults
    ****************************************************************************/
    void GetDeviceInfo()
    {
        string[] description = {
                                       "Driver Version",
                                       "USB Version",
                                       "Hardware Version",
                                       "Variant Info",
                                       "Serial",
                                       "Calibration Date",
                                       "Kernel Version",
                                       "Digital Hardware",
                                       "Analogue Hardware",
                                       "Firmware 1",
                                       "Firmware 2"
                                    };

        System.Text.StringBuilder line = new System.Text.StringBuilder(80);

        if (_handle >= 0)
        {
            for (int i = 0; i < description.Length; i++)
            {
                short requiredSize;
                Imports.GetUnitInfo(_handle, line, 80, out requiredSize, (uint)i);
                Debug.Log(String.Format("{0}: {1}", description[i], line));

                if (_powerSupplyConnected)
                {
                    if (i == 3)
                    {
                        _channelCount = int.Parse(line[1].ToString());

                    }
                }
                else
                {
                    _channelCount = 2;
                }
                _channelSettings = new ChannelSettings[_channelCount];
                _noEnabledChannels = _channelCount;

            }

            Imports.Range voltageRange = Imports.Range.Range_200mV;

            for (int ch = 0; ch < _channelCount; ch++)
            {
                Imports.SetChannel(_handle, Imports.Channel.ChannelA + ch, 1, Imports.Coupling.PS5000A_DC, voltageRange, 0);
                _channelSettings[ch].enabled = true;
                _channelSettings[ch].range = voltageRange;
                Debug.Log(String.Format("{0}: VoltageRange: {1}", Imports.Channel.ChannelA + ch, voltageRange));
            }
        }

    }

    /****************************************************************************
    * CollectStreamingImmediate
    *  this function demonstrates how to collect a stream of data
    *  from the unit (start collecting immediately)
    ***************************************************************************/
    void CollectStreamingImmediate()
    {
        /* Trigger disabled	*/
        Imports.SetSimpleTrigger(_handle, 0, Imports.Channel.ChannelA, 0, Imports.ThresholdDirection.None, 0, 0);

        var previousTime = DateTime.Now;
        uint numberRawSamples = 1000000;
        StreamDataHandler(0, numberRawSamples);
        //Debug.Log(String.Format("number of samples: {0}   elapsed time: {1}", numberRawSamples, DateTime.Now - previousTime));
    }
    
    /// <summary>
    /// Returns the max of the fft data of channel A
    /// returns null if not done copying
    /// </summary>
    /// use measurement index to write to measurementPoint? so it can be done on different thread than main?
    /// how to write from here? events cant be used outside main thread?
    /// <returns></returns>
    public Tuple<double, double> GetData()
    {
        if (doneCopying && !usePingPong)
        {
            double[] dataCha_mv = new double[200000];
            //this used to be till i < current200kPackage.length
            //however it sometimes gives a nullpointer exception.
            //to not let it completely wipe dataCHa_mv with zeroes just explicitely
            //setting it to 200000 seemed better.
            for (int i = 0; i < 200000; i++)
            {
                try
                {
                    //for different channels see build file body of original picoscope example script:
                    //adc_to_mv(appBuffers[ch][i], (int)_channelSettings[(int)(Imports.Channel.ChannelA + (ch / 2))].range)
                    dataCha_mv[i] = (double)adc_to_mv(current200kPackage[0][i], (int)_channelSettings[0].range);
                }
                catch (NullReferenceException e)
                {
                    Debug.Log(String.Format("something was null but what....: {0}", e));
                    if (i - 1 < 0)
                    {
                        dataCha_mv[i] = 0.0;
                    }
                    else
                    {
                        dataCha_mv[i] = dataCha_mv[i - 1];
                    }
                }
            }

            //if there are 200,000 samples. PostProcessing can begin
            return post.CalculateFFT(dataCha_mv);
        }
        else if (usePingPong)
        {
            System.Random rnd = new System.Random();
            return Tuple.Create(150 * rnd.NextDouble(), 20000.0 * rnd.Next(1,5));
        }
        else
        {
            Debug.Log(String.Format("Not done copying from appBuffers to currentPackage"));
            return null;
        }
    }

    public void Run()
    {
        while (running)
        {
            CollectStreamingImmediate();
        }
        Debug.Log(String.Format("DataCollectionThread stopping"));
    }

    public void StartThread()
    {
        Thread dataCollectionThread = new Thread(Run)
        {
            Name = "dataCollectionThread"
        };
        running = true;
        dataCollectionThread.Start();
    }

    public PicoScopeHandler(short handle, bool powerSupplyConnected)
    {
        _handle = handle;
        _powerSupplyConnected = powerSupplyConnected;
        GetDeviceInfo();
        Imports.MaximumValue(_handle, out _maxValue);
        post = new PostProcessing();
        doneCopying = false; //not allowing grabbing of values until first package
    }

    public PicoScopeHandler(bool pingPong)
    {
        doneCopying = true; //Otherwise data post thread will sleep forever
        usePingPong = pingPong;
        post = new PostProcessing();
    }

    /// <summary>
    /// Opens the picoscope unit. Returns the handle and the bool powerSupplyConnected.
    /// 
    /// </summary>
    /// <returns></returns>
    public static Tuple<short, bool> PerformSetup()
    {
        Debug.Log(String.Format("Enumerating devices...\n"));

        short count = 0;
        short serialsLength = 40;
        StringBuilder serials = new StringBuilder(serialsLength);

        uint status = Imports.EnumerateUnits(out count, serials, ref serialsLength);

        if (status != StatusCodes.PICO_OK)
        {
            Debug.Log(String.Format("No devices found.\n"));
            Debug.Log(String.Format("Error code : {0}", status));
            throw new Exception(String.Format("Error code : {0}", status));
        }
        else
        {
            if (count == 1)
            {
                Debug.Log(String.Format("Found {0} device:", count));
            }
            else
            {
                Debug.Log(String.Format("Found {0} devices", count));
            }

            Debug.Log(String.Format("Serial(s) {0}", serials));

        }

        // Open unit and show splash screen
        Debug.Log(String.Format("\nOpening device..."));

        short handle;

        status = Imports.OpenUnit(out handle, null, Imports.DeviceResolution.PS5000A_DR_8BIT);

        bool powerSupplyConnected = true;

        Debug.Log(String.Format("Handle: {0}", handle));

        if (status == StatusCodes.PICO_POWER_SUPPLY_NOT_CONNECTED || status == StatusCodes.PICO_USB3_0_DEVICE_NON_USB3_0_PORT)
        {
            status = Imports.ChangePowerSource(handle, status);
            powerSupplyConnected = false;
        }
        else if (status != StatusCodes.PICO_OK)
        {
            Debug.Log(String.Format("Cannot open device error code: " + status.ToString()));
            throw new Exception(String.Format("Error code : {0}", status));
        }
        else
        {
            // Do nothing - power supply connected
        }

        if (status != StatusCodes.PICO_OK)
        {
            Debug.Log(String.Format("Unable to open device."));
            Debug.Log(String.Format("Error code : {0}", status));
        }
        else
        {
            Debug.Log(String.Format("Device opened successfully.\n"));
        }
        return Tuple.Create(handle, powerSupplyConnected);
    }

    /// <summary>
    /// Closes the connection with the picoscope with given handle. sets running to false to stop datacollection thread.
    /// </summary>
    /// <param name="handle"></param>
    public void StopPico(short handle)
    {
        Debug.Log("Stopping PicoScope");
        running = false;
        Imports.CloseUnit(handle); 
    }
}
    struct ChannelSettings
{
    public Imports.Range range;
    public bool enabled;
}
